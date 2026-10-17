from itertools import product
import json, os, logging
from pymatgen.core import Species, Element
from pymatgen.core.structure import Structure
from pymatgen.analysis.bond_valence import BVAnalyzer


class checkOxidation():
    def __init__(self, tolerance=0.01, TEST_JSON_PATH="TEST_QUERY.json", OUTPUT_JSON_PATH="OXIDATION_QUERY.json",
                 CIF_DIR="QUERY/CIF", OXI_CIF_DIR="QUERY/OXI_CIF",
                 logger=None):
        """
        Checks whether a material's composition can satisfy charge
        neutrality.

        Primary method: per-SITE oxidation-state assignment via
        pymatgen's BVAnalyzer on actual structure

        Fallback method: if BVAnalyzer can't find a valid bond-valence
        solution for the structure, falls back to the old per-ELEMENT
        combinatorial search over Materials Project's `possible_species`
        field

        tolerance : float, optional
            How far off zero total charge (per formula unit) is still
            counted as "neutral". Default 0.01.

        TEST_JSON_PATH : str, optional
            Filename of the input query file to read from.

        OUTPUT_JSON_PATH : str, optional
            Filename to write oxidation-state results to.

        CIF_DIR : str, optional
            Folder containing the raw (undecorated) .cif files.

        OXI_CIF_DIR : str, optional
            Folder to write oxidation-decorated .cif files to. Only
            written when is_neutral is True
        """

        self.tolerance = tolerance
        self.json_path = TEST_JSON_PATH
        self.output_json = OUTPUT_JSON_PATH
        self.cif_dir = CIF_DIR
        self.oxi_cif_dir = OXI_CIF_DIR

        os.makedirs(self.oxi_cif_dir, exist_ok=True)

        self.logger = logger if logger is not None else logging.getLogger(__name__)

        self._TIER_SEVERITY = {"common": 0, "icsd": 1, "exotic": 2, "unlisted": 3}

    def load_json(self):
        with open(self.json_path, "r", encoding="utf-8") as f:
            return json.load(f)

    def save_json(self, metadata):
        """
        Save metadata to JSON
        """
        if os.path.exists(self.output_json):
            with open(self.output_json, "r") as f:
                try:
                    results = json.load(f)
                    if not isinstance(results, list):
                        results = [results]
                except json.JSONDecodeError:
                    results = []
        else:
            results = []

        results.append(metadata)
        with open(self.output_json, "w") as f:
            json.dump(results, f, indent=2)

    def _load_structure(self, material_id):
        cif_path = os.path.join(self.cif_dir, f"{material_id}.cif")
        return Structure.from_file(cif_path)

    def _write_oxi_cif(self, structure, material_id):
        oxi_cif_path = os.path.join(self.oxi_cif_dir, f"{material_id}.cif")
        structure.to(filename=oxi_cif_path)

    def _rank_states(self, element_symbol, states):
        """
        Order one element's oxidation states based on its chemical
        plausibility tier and tags with their label.
        common (1): own preferred order
        icsd (2): covered in icsd
        exotic (3): else
        """

        try:
            el = Element(element_symbol)
            common = [s for s in el.common_oxidation_states if s in states]

            icsd = [
                s for s in el.icsd_oxidation_states
                if s in states and s not in common
            ]
        except ValueError:
            common, icsd = [], []

        rest = sorted(
            s for s in states
            if s not in common and s not in icsd
        )
        return ([(s, "common") for s in common] +
                [(s, "icsd") for s in icsd] +
                [(s, "exotic") for s in rest])

    #per-site path (primary)

    def _compute_site_charge(self, structure):
        """
        Sums (oxidation_state * occupancy) across every site. Disordered
        sites contribute an occupancy-weighted sum over their species.
        """
        total = 0.0
        for site in structure:
            if site.is_ordered:
                total += site.specie.oxi_state
            else:
                for sp, occu in site.species.items():
                    #example: {Species("Fe", 2): 0.6} -> sp.element and sp.oxi_state
                    total += sp.oxi_state * occu
        return total

    def _tag_state_tier_per_site(self, structure, doc):
        """
        Cross-checks each site's BVAnalyzer-assigned oxidation state
        against that element's possible_species-derived candidates (from
        MP) as a confidence label
        """
        candidates = self._get_oxidation_states(doc)
        tiers = []

        for site in structure:
            species_list = [(site.specie, 1.0)] if site.is_ordered else list(site.species.items())
            for sp, _occu in species_list:
                element = str(sp.element)
                el_candidates = candidates.get(element)
                if el_candidates is None:
                    tiers.append("unlisted")
                    continue
                tier_map = dict(self._rank_states(element, el_candidates))
                tiers.append(tier_map.get(sp.oxi_state, "unlisted"))

        if not tiers:
            return None
        return max(tiers, key=lambda t: self._TIER_SEVERITY.get(t, 99))

    def _check_charge_per_site(self, structure, doc):
        """
        Primary path: BVAnalyzer per-site bond-valence assignment.
        Mutates `structure` in place (tags oxidation states by site).
        Lets any exception from BVAnalyzer/pymatgen propagate -- the
        caller catches it and falls back to the per-element method.
        """
        valences = BVAnalyzer().get_valences(structure)
        structure.add_oxidation_state_by_site(valences)

        total_charge = self._compute_site_charge(structure)
        _, z = structure.composition.get_reduced_composition_and_factor()
        charge_per_formula_unit = total_charge / z if z else total_charge

        is_neutral = abs(charge_per_formula_unit) <= self.tolerance
        state_tier = self._tag_state_tier_per_site(structure, doc)

        return is_neutral, round(charge_per_formula_unit, 4), state_tier

    # per-element path 

    def _apply_oxidation_element(self, structure, assignment):
        """
        Replaces every site's species with the per-element oxidation
        state chosen by _find_best_assigment, mutating `structure` in
        place.
        """
        for i, site in enumerate(structure):
            if site.is_ordered:
                el = site.specie.symbol
                new_species = Species(el, assignment[el])
            else:
                new_species = {
                    Species(sp.symbol, assignment[sp.symbol]): occu
                    for sp, occu in site.species.items()
                }

            structure.replace(
                i,
                new_species,
                coords=site.frac_coords,
                coords_are_cartesian=False,
                properties=site.properties
            )
        return structure

    def _get_composition(self, structure):
        """
        Extract reduced (per-formula-unit) composition. Charges are
        compared against self.tolerance on a per-formula-unit basis
        """
        reduced_comp, _z = structure.composition.get_reduced_composition_and_factor()
        return {
            el.symbol: float(amount)
            for el, amount in reduced_comp.items()
        }

    def _get_oxidation_states(self, doc):
        """
        Reads doc["possible_species"], groups the oxidation states by
        element. Non-integer oxidation states are dropped.
        """
        candidates = {}
        dropped = []
        for species in doc.get("possible_species", []):
            sp = Species.from_str(species)
            if not float(sp.oxi_state).is_integer():
                dropped.append(species)
                continue
            element = str(sp.element)
            candidates.setdefault(element, set()).add(sp.oxi_state)

        if dropped:
            self.logger.warning(
                f"{doc.get('material_id')}: dropped fractional oxidation "
                f"state candidate(s) {dropped}"
            )
        return candidates


    def _evaluate_combination(self, composition, elements, oxidation_states):
        """
        Adds up total charge for one specific choice of oxidation state
        per element: sum(oxidation_state * how_many_atoms) across all
        elements.
        """
        return round(sum(ox * composition[element]
                for element, ox in zip(elements, oxidation_states))
                , 4)

    def _find_best_assigment(self, composition, candidates):
        """
        Tries every combination of oxidation states (one per element, all
        combos from `candidates`, lookup from rank states) and looks for
        one that adds up to ~0 total charge.
        """
        elements = list(composition.keys())
        option_lists = [self._rank_states(e, candidates[e]) for e in elements]

        best_charge = None
        best_assignment = None
        best_tier = None

        for combo in product(*option_lists):
            states = tuple(state for state, _tier in combo)
            tiers = [tier for _state, tier in combo]
            overall_tier = max(tiers, key=lambda t: self._TIER_SEVERITY[t])

            charge = self._evaluate_combination(composition, elements, states)
            assignment = dict(zip(elements, states))

            if best_charge is None or abs(charge) < abs(best_charge):
                best_charge = charge
                best_assignment = assignment
                best_tier = overall_tier

            if abs(charge) <= self.tolerance:
                return True, charge, assignment, overall_tier

        return False, best_charge, best_assignment, best_tier

    def _check_charge_per_element(self, structure, doc):
        """
        Fallback path: one oxidation state per element, applied
        uniformly to every site of that element. Only reached when
        BVAnalyzer can't resolve a per-site assignment.
        """
        material_id = doc.get("material_id")

        if not doc.get("possible_species"):
            self.logger.error(f"Possible species of {material_id} is not found")
            return False, None, None, None

        composition = self._get_composition(structure)
        candidates = self._get_oxidation_states(doc)

        if not all(element in candidates for element in composition):
            self.logger.error(f"Inconsistency in {material_id} oxidation states and composition")
            return False, None, None, None

        is_neutral, charge, assignment, state_tier = self._find_best_assigment(composition, candidates)
        return is_neutral, charge, assignment, state_tier

    # ---------------- main entry point ----------------

    def check_charge(self, doc):
        """
        Main entry point: given one Materials Project doc, tries to
        assign oxidation states that balance to zero charge (per formula
        unit, within self.tolerance).

        Tries BVAnalyzer's per-site bond-valence method first. If fails, 
        falls back to the per-element combinatorial search over possible_species 

        Writes an oxidation-decorated CIF to OXI_CIF_DIR only when a
        neutral assignment is found.

        """
        material_id = doc.get("material_id")
        formula_pretty = doc.get("formula_pretty")

        result = {
            "material_id": material_id,
            "formula_pretty": formula_pretty,
            "is_neutral": False,
            "total_charge": None,
            "state_tier": None,
            "method": None,
        }

        try:
            structure = self._load_structure(material_id)
        except Exception as e:
            self.logger.error(f"Could not load structure for {material_id}: {e!r}")
            return result

        try:
            is_neutral, charge, state_tier = self._check_charge_per_site(structure, doc)
            result.update({
                "is_neutral": is_neutral,
                "total_charge": charge,
                "state_tier": state_tier,
                "method": "bvanalyzer",
            })
            if is_neutral:
                self._write_oxi_cif(structure, material_id)
            return result

        except Exception as e:
            self.logger.warning(
                f"BVAnalyzer could not assign per-site valences for {material_id} "
                f"({e!r}); falling back to per-element assignment"
            )

        # fallback
        try:
            structure = self._load_structure(material_id)
            is_neutral, charge, assignment, state_tier = self._check_charge_per_element(structure, doc)
            result.update({
                "is_neutral": is_neutral,
                "total_charge": charge,
                "state_tier": state_tier,
                "method": "element_fallback",
            })
            if is_neutral:
                self._apply_oxidation_element(structure, assignment)
                self._write_oxi_cif(structure, material_id)
            return result

        except (KeyError, TypeError, ValueError) as e:
            self.logger.error(f"Oxidation-state check failed ({material_id}): {e!r}")
            return result