import json, math, os
from pymatgen.core import Species
import traceback
import logging
_CN_ROMAN = {1: "I", 2: "II", 3: "III", 4: "IV", 5: "V", 6: "VI",
             7: "VII", 8: "VIII", 9: "IX", 10: "X", 11: "XI", 12: "XII"}
_QUALITY_RANK = {"shannon": 0, "highspin": 1, "fallback_A": 2, "fallback_B": 3,
                 "partial": 4, "failed": 5}

class checkTolerance():
    def __init__(self, NEIGH_PATH="NEIGHBOR_QUERY.json",
                OUTPUT_JSON = "TOLERANCE_QUERY.json", logger = None):
        """
        Computes perovskite tolerance factors (octahedral factor mu,
        Goldschmidt t, Bartel tau) from checkOxidation.py and
        checkNeighbor.py results. 
        Every radius/oxidation-state value is a lookup by element, not a
        measurement off the actual structure.
 
        NEIGH_PATH : str, filename of checkNeighbor.py's output
        """
        self.NEIGH_PATH = NEIGH_PATH
        self.logger  = logger if logger is not None else logging.getLogger(__name__)

    def _collect_sites(self, entry):
        """
        Split neighbor's cn dict into A and B site records by role.
    
        Returns lists of {element, oxi_state, CN} -- o
        ne per crystallographic
        site, multiplicity preserved.
        """
        out = {"A": [], "B": [], "X": []}
        for d in entry["cn"].values():
            role = d.get("role")
            if role in out:
                out[role].extend([d] * int(round(d.get("weight", 1))))
        return out

    def _snap_to_shannon(self, actual_cn):
        """
        Snaps a float/measured CN to the nearest available Shannon CN integer.
        """
        valid_cns = sorted(_CN_ROMAN.keys())
        return min(valid_cns, key=lambda x: abs(x - actual_cn))

    def _get_radius(self, element, oxi_state, cn):
            """
            Shannon ionic radius for one site.
    
            Returns (radius, quality) where quality is:
              "shannon"  -- tabulated at this exact (element, oxi, CN)
              "highspin" -- tabulated, but only via the High Spin variant
              "fallback_A" -- Species.ionic_radius with known CN (A-site) >= 8
              "fallback_B" -- Species.ionic_radius with known CN (A-site) < 8
              "missing"  -- no data; radius is None
            """
            if oxi_state is None:
                return None, "missing"
    
            sp = Species(element, oxi_state)
            cn_roman = _CN_ROMAN[cn]
    
            try:
                return sp.get_shannon_radius(cn_roman), "shannon"
            except KeyError:
                pass
    
            try:
                return sp.get_shannon_radius(cn_roman, spin="High Spin"), "highspin"
            except KeyError:
                pass
    
            if sp.ionic_radius is not None:
                if cn >= 8:
                    return float(sp.ionic_radius), "fallback_A"
                else:
                    return float(sp.ionic_radius), "fallback_B"
    
            return None, "missing"
 
    def _group_radius(self, sites, cn_override=None):
            """
            Mean radius over one site group (A, B or X), with provenance.

            cn_override: force a CN (6 for B and X). None means use each site's
            measured CN, snapped to the nearest Shannon value.

            Returns (mean_radius | None, quality, detail) where quality is the
            WORST tier seen: shannon < highspin < fallback < partial < failed.
            detail["tiers"] records one (element, oxi_state, cn, quality) tuple
            per resolved site, so a fallback can be attributed to the exact
            lookup that caused it.
            """
            radii, qualities, missing, tiers = [], [], [], []

            for d in sites:
                cn = cn_override if cn_override is not None else self._snap_to_shannon(d["CN"])
                r, q = self._get_radius(d["element"], d["oxi_state"], cn)
                if r is None:
                    missing.append(d["element"])
                else:
                    radii.append(r)
                    qualities.append(q)
                    tiers.append((d["element"], d["oxi_state"], cn, q))

            detail = {"missing": missing, "n_used": len(radii), "tiers": tiers}

            if not radii:
                return None, "failed", detail

            tier = max(qualities, key=lambda q: _QUALITY_RANK[q])
            if missing:
                tier = "partial"

            return sum(radii) / len(radii), tier, detail

    def _n_a(self, a_sites):
        """Mean A-site oxidation state (Bartel's n_A). None if any site untagged."""
        states = [d["oxi_state"] for d in a_sites]
        if not states or any(s is None for s in states):
            return None
        return sum(states) / len(states)
        
    def _bartel_tau(self, rA, rB, rX, nA):
        if nA is None or rA <= rB:
            return None
        ratio = rA/rB 
        return rX/rB-nA*(nA-ratio/math.log(ratio))

    def get_tolerance_factor(self, entry):
        mpid = entry.get("material_id")
        site = self._collect_sites(entry)

        blank = {
            "material_id": mpid,
            "a_elements": None, "b_elements": None, "x_elements": None,
            "r_A": None, "r_B": None, "r_X": None, "n_A": None,
            "octahedral_factor": None, "goldschmidt_t": None, "bartel_tau": None,
            "radius_quality": "failed",
            "radius_detail": {"A": None, "B": None, "X": None},
        }

        if not entry.get("is_perovskite"):
            self.logger.info(f"{mpid} is not perovskite; no tolerance factors")
            return blank

        site = self._collect_sites(entry) 
        a_sites, b_sites, x_sites = site["A"], site["B"], site["X"]

        r_A, qA, dA = self._group_radius(a_sites, cn_override=12)
        r_B, qB, dB = self._group_radius(b_sites, cn_override=6)
        r_X, qX, dX = self._group_radius(x_sites, cn_override=6)

        quality = max((qA, qB, qX), key=lambda q: _QUALITY_RANK[q])
        detail = {"A": dA, "B": dB, "X": dX}

        if None in (r_A, r_B, r_X):
            self.logger.warning(
                f"{mpid}: radius group failed (A={qA}, B={qB}, X={qX}); "
                f"missing {dA['missing']} / {dB['missing']} / {dX['missing']}"
            )
            return {**blank,
                    "a_elements": entry["a_true"],
                    "b_elements": entry["b_true"],
                    "x_elements": entry["x_elements"],
                    "radius_quality": quality,
                    "radius_detail": detail}

        n_A = self._n_a(a_sites)
        mu = r_B / r_X
        t = (r_A + r_X) / (math.sqrt(2) * (r_B + r_X))
        tau = self._bartel_tau(r_A, r_B, r_X, n_A)

        if tau is None:
            self.logger.warning(
                f"{mpid}: Bartel tau undefined (n_A={n_A}, r_A={r_A:.3f}, "
                f"r_B={r_B:.3f}); mu and t still valid"
            )

        return {
            "material_id": mpid,
            "a_elements": entry["a_true"],
            "b_elements": entry["b_true"],
            "x_elements": entry["x_elements"],
            "r_A": r_A, "r_B": r_B, "r_X": r_X, "n_A": n_A,
            "octahedral_factor": mu,
            "goldschmidt_t": t,
            "bartel_tau": tau,
            "radius_quality": quality,
            "radius_detail": detail,
        }