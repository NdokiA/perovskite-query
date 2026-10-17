from pymatgen.core import Structure
from collections import Counter
from pymatgen.analysis.local_env import CrystalNN
import os, json, logging
import numpy as np

class checkNeighbor():
    def __init__(self, JSON_PATH = "TEST_QUERY.json", OXI_CIF_DIR = "TEST_CIF",
                 OUTPUT_JSON = "NEIGHBOR_QUERY.json", logger = None):
        """
        Given a crystal structure (.cif), checks whether its cation sites
        form BX6 octahedra (6 anion neighbors) and whether those octahedra
        connect to each other by sharing corners — the defining feature of
        a perovskite-type structure.

        JSON_PATH : str, optional
            Filename of the query file listing materials
            to check (each entry needs at least "material_id", and either
            "cif_file" or a CIF named after its material_id).


        OXI_CIF_DIR : str, optional
            Folder containing the .cif files referenced by JSON_PATH.
        """

        self.JSON_PATH = JSON_PATH
        self.OXI_CIF_DIR = OXI_CIF_DIR
        self.OUTPUT_PATH = OUTPUT_JSON
        self.logger = logger if logger is not None else logging.getLogger(__name__)

    def load_json(self):
        with open(self.JSON_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    
    def save_json(self, metadata):
        """
        Save metadata to JSON
        """

        if os.path.exists(self.OUTPUT_PATH):
            with open(self.OUTPUT_PATH, "r") as f:
                try:
                    results = json.load(f)
                    if not isinstance(results, list):
                        results = [results]
                except json.JSONDecodeError:
                    results = []
        else:
            results = []
        
        results.append(metadata)
        with open(self.OUTPUT_PATH, "w") as f:
            json.dump(results, f, indent=2)
    
    def _extract_cif(self,doc):
        """
        Loads the structure from its .cif file
        """

        if doc.get("cif_file"):
            cif_path = os.path.join(self.OXI_CIF_DIR, os.path.basename(doc["cif_file"]))
        else:
            cif_path = os.path.join(self.OXI_CIF_DIR, str(doc["material_id"]) + ".cif")
        
        structure = Structure.from_file(cif_path)

        return structure

    def _classify_sites(self, structure):
        """
        Classifies each site as anion (X) or cation (B-candidate)

        """
        x_idx, b_idx = [], []
 
        for i, site in enumerate(structure):
            species_list = [site.specie] if site.is_ordered else list(site.species.keys())
 
            signs = set()
            for sp in species_list:
                oxi = getattr(sp, "oxi_state", None)
                if oxi is None:
                    continue
                signs.add(1 if oxi > 0 else (-1 if oxi < 0 else 0))
 
            if not signs:
                self.logger.warning(
                    f"Site {i} ({self._site_label(structure, i)}) has no oxidation "
                    f"state tagged; skipping from anion/cation classification"
                )
                continue
 
            if -1 in signs:
                x_idx.append(i)
            if 1 in signs:
                b_idx.append(i)
            if len(signs) > 1:
                self.logger.warning(
                    f"Site {i} ({self._site_label(structure, i)}) hosts species of "
                    f"mixed sign ({sorted(signs)}); classified as both anion and "
                    f"cation candidate"
                )
 
        return x_idx, b_idx
    
    def _get_site_cn(self, structure, site_idx, anion_indices, cnn=None):
        """
        Finds the anion neighbors of one cation site using CrystalNN, and
        returns them as a set of (site_index, image) pairs. The count of
        this set is the site's coordination number (CN).

        Using (site_index, image) instead of just site_index because
        two neighbors can be the SAME crystallographic site but different
        periodic copies
        """

        cnn = cnn or CrystalNN()
        nn = cnn.get_nn_info(structure, site_idx)
        #key by (site_index, image) so two different periodic copies of the 
        # same crystallographic site are not collapsed into one neighbor
        anion_nn = {
            (n['site_index'], n['image'])
            for n in nn if n['site_index'] in anion_indices
        }
        return anion_nn
    
    def _build_polyhedra_map(self, structure, cation_indices, anion_indices, cnn = None):
        """
        Runs _get_site_cn for every candidate cation site 
        output:     {cation_site_index: set of anion neighbors}.
        """

        cnn = cnn or CrystalNN()
        return {ci: self._get_site_cn(structure, ci, anion_indices, cnn)
                for ci in cation_indices}
    
    def _check_cn(self, structure, poly_map, cn_target=6, tol=0):
        """
        CN pass/fail check (default tol=0, strict for BX6)
        """
        return {ci: {
                    "element": self._site_label(structure, ci),
                    "oxi_state": self._site_oxi(structure, ci),
                    "CN": len(anions), 
                    "CN_ok": abs(len(anions) - cn_target)<=tol
                    }
                     for ci, anions in poly_map.items()
                     }

    def _check_cn_anion(self, structure, poly_map, x_idx, b_sites):
        '''
        WARNING: "CN" here is NOT a coordination number.
        Count of confirmed B sites this anion ligates
        checkTolerance ignores it and forces CN=6 for the Shannon
        radius lookup;
        '''


        bonds = Counter()
        for ci in b_sites:
            for (a, _img) in poly_map[ci]:
                bonds[a] += 1

        return {i: {"element": self._site_label(structure, i),
                    "oxi_state": self._site_oxi(structure, i),
                    "CN": bonds.get(i, 0),
                    "CN_ok": None,
                    "role": "X",
                    }
                for i in x_idx}
    
    def _classify_share_mode(self, n_shared):
        """

        Turns a shared-anion count into the usual crystallography name:
        1 shared anion = corner-sharing, 
        2 = edge-sharing, 
        3 = face-sharing,
        0 = not touching
        >3 = unusual — labeled "over-face", maybe geometry issue, double-check
        """

        return {0: 'none', 1: 'corner', 2: 'edge', 3: 'face'
                }.get(n_shared, 'over-face')
    
    def _site_label(self, structure, idx):
            """Species symbol for a site (majority species if disordered)."""
            site = structure[idx]
            sp = site.specie if site.is_ordered else max(site.species.items(),
                                                         key=lambda kv: kv[1])[0]
            return sp.symbol

    def _site_oxi(self, structure, idx):
        """Oxidation state of a site's (majority) species, None if untagged."""
        site = structure[idx]
        sp = site.specie if site.is_ordered else max(site.species.items(),
                                                     key=lambda kv: kv[1])[0]
        oxi = getattr(sp, "oxi_state", None)
        return float(oxi) if oxi is not None else None
        
    def _find_edges(self, poly_map, i, j):

        """
        Find all distinct periodic 'edges' (specific corner/edge/face-sharing
        contacts)  site i and site j.

        Returns dict: T (image offset) -> shared anion set for that specific
        instance pair.
        """

        by_anion_j = {}

        for (a, img) in poly_map[j]:
            by_anion_j.setdefault(a, []).append(img)

        edges = {}
        for (a, imgA) in poly_map[i]:
            for imgB in by_anion_j.get(a, []):
                T = tuple(round(x - y) for x, y in zip(imgA, imgB))
                if i == j:
                    if T == (0, 0, 0):
                        continue  # same instance, not a neighbor

                    if T > tuple(-x for x in T):
                        continue  # dedupe: T and -T are the same physical bond

                edges.setdefault(T, set()).add((a, imgA))
        return edges



    def _check_pairwise_connectivity(self, structure, poly_map):
        """
        Goes through every pair of cation sites (including a site paired
        with itself) and records each touching point found by
        _find_edges as one entry: which two sites, how many anions they
        share, what kind of sharing it is (corner/edge/face), and
        whether the two sites are different elements ("hetero").
        """

        pairs = []
        keys = list(poly_map.keys())

        for idx_a in range(len(keys)):
            for idx_b in range(idx_a, len(keys)):
                i, j = keys[idx_a], keys[idx_b]
                edges = self._find_edges(poly_map, i, j)

                for T, shared in edges.items():

                    n = len(shared)
                    if n == 0:
                        continue

                    pairs.append({
                        "i": i, "j": j, "image": T, "n_shared": n,
                        "mode": self._classify_share_mode(n),
                        "hetero": self._site_label(structure, i) != self._site_label(structure, j)
                    })
        return pairs

    def is_perovskite(self, cn_results, conn_results, corner_frac_threshold=0.9):
            """
            Computes BX6 connectivity structure by corner-sharing check > threshold
    
            Also computes the connectivity statistics over that same confirmed
            B-B subset, so callers don't re-filter conn_results.
    
            Returns
            -------
            has_net : bool
                The perovskite verdict. True if a CN-ok site exists AND the
                corner-sharing criterion is met 
            b_sites : list[int]
                Site indices that passed the CN check
            corner_frac : float
                Fraction of B-B contacts that are corner-sharing, in [0, 1].
                note this is indistinguishable from a genuine all-edge/face
                structure, so read it together with len(b_conn).
            hetero_corner_frac : float | None
                Fraction of the CORNER contacts (not all contacts) that link
                two different B species (the B-X-B' signature)
            b_conn : list[dict]
                The connectivity entries from conn_results restricted to pairs
                where BOTH sites are confirmed B sites. 
            """
            b_sites = [ci for ci, r in cn_results.items() if r["CN_ok"]]
            if not b_sites:
                return False, b_sites, 0.0, None, []
    
            b_set = set(b_sites)
            b_conn = [p for p in conn_results if p["i"] in b_set and p["j"] in b_set]
            if not b_conn:
                return False, b_sites, 0.0, None, []
    
            n_corner = sum(p["mode"] == "corner" for p in b_conn)
            n_hetero_corner = sum(p["mode"] == "corner" and p["hetero"] for p in b_conn)
    
            corner_frac = n_corner / len(b_conn)
            hetero_corner_frac = (n_hetero_corner/ n_corner
                if n_corner else None
            )
    
            has_net = corner_frac >= corner_frac_threshold
            return has_net, b_sites, corner_frac, hetero_corner_frac, b_conn

    def _network_rank(self, b_conn):
        """
        Rank of the lattice-translation span of the corner-sharing B-B
        network: 3 = 3D framework, 2 = layered (RP/DJ/Aurivillius),
        1 = chains, 0 = isolated or no corner contacts.

        Lower bound only. Contacts internal to one cell carry T=(0,0,0)
        and contribute no translation, so a genuine 3D framework described
        in a large supercell can report < 3. Report, don't gate.
        """
        T = [p["image"] for p in b_conn if p["mode"] == "corner"]
        if not T:
            return 0
        return int(np.linalg.matrix_rank(np.array(T, dtype=float)))
        
    def verify_bx6(self, doc, cn_target=6, cn_tol=0, corner_threshold=0.9):
            """
            Main entry point. Classifies sites into anions (X) and cation
            candidates, builds each cation's anion-neighbor set, checks which
            hit CN==6, then checks how those octahedra connect including through periodic 
            images.
            Sharing stats are restricted to confirmed CN==6 sites
            """
            structure = self._extract_cif(doc)
            x_idx, b_idx = self._classify_sites(structure)
            x_elements = sorted({self._site_label(structure, i) for i in x_idx})
            b_elements = sorted({self._site_label(structure, i) for i in b_idx})

            x_oxi = {self._site_label(structure, i): self._site_oxi(structure, i)
                      for i in x_idx}
        
            cnn = CrystalNN()
            poly_map = self._build_polyhedra_map(structure, b_idx, x_idx, cnn)
            cn_results = self._check_cn(structure, poly_map, cn_target, cn_tol)
            conn_results = self._check_pairwise_connectivity(structure, poly_map)
    
            perov, b_sites, corner_frac, hetero_corner_frac, b_conn = self.is_perovskite(
                cn_results, conn_results, corner_frac_threshold=corner_threshold)
    
            b_set = set(b_sites)

            x_results = self._check_cn_anion(structure, poly_map, x_idx, b_sites)

            for ci, r in cn_results.items():                       
                r["role"] = "B" if ci in b_set else "A"

            overlap = set(cn_results) & set(x_results)              # mixed-sign sites
            for i in overlap:
                self.logger.warning(
                    f"Site {i} ({self._site_label(structure, i)}) is both anion and "
                    f"cation candidate; keeping cation role, excluded from X group")
            cn_results.update({i: r for i, r in x_results.items() if i not in overlap})

            b_true = sorted({self._site_label(structure, ci) for ci in b_sites})
            a_true = sorted({self._site_label(structure, ci)
                             for ci in b_idx if ci not in b_set})

            network = self._network_rank(b_conn)
    
            return {
                "material_id": doc.get("material_id"),
                "formula_pretty": doc.get("formula_pretty"),
                "is_perovskite": perov,
                "x_elements": x_elements,
                "x_oxi": x_oxi,
                "b_candidates": b_elements,
                "b_true": b_true,
                "a_true": a_true,
                "n_b_sites": len(b_sites),
                "cn": cn_results,
                "corner_share_fraction": corner_frac,
                "hetero_corner_fraction": hetero_corner_frac,
                "network_rank": network,
                "n_b_contacts": len(b_conn),
                "share_modes": {m: sum(p["mode"] == m for p in b_conn)
                            for m in ("corner", "edge", "face", "over-face")},
            }