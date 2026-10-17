import os, logging, json, time
from tqdm import tqdm
from config import mpAPI
from mp_api.client import MPRester
from .queryStructure import queryStructure as qS
from pathlib import Path

FIELDS  = [
    # Identifiers
    "material_id", "formula_pretty", "formula_anonymous",
    "chemsys", "nelements", "elements",
    "deprecated", "theoretical",
    # Structural
    "symmetry",
    "composition",
    "nsites",
    "structure",
    # Stability
    "formation_energy_per_atom",
    "energy_above_hull",
    "is_stable",
    # Electronic
    "band_gap", "is_gap_direct", "is_metal",
    "cbm", "vbm", "efermi",
    # Oxidation states
    "possible_species",
    "composition_reduced"
]

class queryPerovskite_Robocryst(qS):
    def __init__(self, fields=FIELDS, CIF_DIRS="QUERY/CIF", JSON_PATH="QUERY/QUERY.json",
                 ROBO_PATH = "QUERY/ROBOQUERY.json", LOG_PATH="QUERY/LOG_QUERY/queryRobo.log", BATCH_QUERY = 9999):
        """
        Query and structurally verify candidate perovskite structures from Materials Project.
        Candidates are identified via a robocrystallographer keyword search (see `obtain_ID`)

        Parameters
        ----------
        fields : list(str), optional
            Fields of the projected output for MP queries on each structure, default is FIELDS.
        CIF_DIRS : str, optional
            Directory path for CIF files of the queried structure.
        JSON_PATH : str, optional
            File path for the general query JSON output.
        LOG_PATH : str, optional
            File path for the run log (parent directory is created automatically if
            it doesn't exist).
        """

        CIF_DIRS = Path(CIF_DIRS)
        JSON_PATH = Path(JSON_PATH)
        LOG_PATH = Path(LOG_PATH)
        ROBO_PATH = Path(ROBO_PATH)

        for path in [JSON_PATH, LOG_PATH]:
            path.parent.mkdir(parents=True, exist_ok=True)
        os.makedirs(CIF_DIRS, exist_ok=True)
        
        self.log_path = LOG_PATH
        self.robo_path = ROBO_PATH
        self.BATCH_QUERY = BATCH_QUERY

        super().__init__(fields, CIF_DIRS, JSON_PATH, BATCH_QUERY)

    def obtain_info(self):
        """
        Query candidate perovskite entries from Materials Project Robocrystallographer.
    
        Returns
        -------
        list[dict]
            List of Robocrystallographer documents containing all returned fields
            except ``fields_not_requested``.
        """
    
        self._log("Starting robocrys keyword search for 'perovskite'")
    
        with MPRester(mpAPI) as mpr:
            robocrys_docs = list(
                mpr.materials.robocrys.search(
                    keywords=["perovskite", "corner-sharing octahedra", "ABX3 framework", "BO6 octahedra"],
                    chunk_size=1000,   # 1000 for production
                    num_chunks=None,    # None for production
                )
            )
        
        self._log(f"Robocrys search returned {len(robocrys_docs)} candidate documents")
        self._log(f"Total {len(robocrys_docs)} candidate documents successfully obtained")

        self._timed("Processing queried structures", self.processing_info, robocrys_docs)
    
        return robocrys_docs

    def obtain_ID(self, robocrys_docs):
        IDs = {d.material_id.string for d in robocrys_docs}
        return list(IDs)

    def processing_info(self, docs):
        for doc in tqdm(
            docs,
            desc="Processing Robocrys documents",
            unit="doc"
        ):
            data = {
                k: v
                for k, v in doc.model_dump().items()
                if k not in {"fields_not_requested", "last_updated"}
            }
    
            # Convert MPID object to string for filename/JSON
            data["material_id"] = doc.material_id.string
    
            self.save_json(data, self.robo_path)

if __name__ == "__main__":
    query = queryPerovskite_Robocryst()
    IDs = query.obtain_ID()


    
