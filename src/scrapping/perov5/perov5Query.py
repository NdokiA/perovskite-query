import os, logging, json, time
import pandas as pd
from tqdm import tqdm
from .queryStructure import queryStructure as qS
from pymatgen.core import Structure
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pathlib import Path
from config import mpAPI
from mp_api.client import MPRester


class queryPerov5(qS):
    def __init__(self, DIRS="PQUERY"):

        self.DIRS = Path(DIRS)
        self.CIF_DIRS = self.DIRS / "CIF"
        self.JSON_PATH = self.DIRS / "QUERY.json"
        self.LOG_PATH = self.DIRS / "LOG_QUERY" / "query.log"

        for path in [self.JSON_PATH, self.LOG_PATH]:
            path.parent.mkdir(parents=True, exist_ok=True)

        self.CIF_DIRS.mkdir(parents=True, exist_ok=True)

        self.log_path = self.LOG_PATH

        super().__init__(CIF_DIRS=self.CIF_DIRS, JSON_PATH=self.JSON_PATH)

    def load_json(path):
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)

    def _extract_cif(self,doc):
        """
        Loads the structure from its .cif file
        """

        cif_path = os.path.join(self.CIF_DIRS, doc["material_id"] + ".cif")
        structure = Structure.from_file(cif_path)
        return structure
    
    def save_cif(self, data):
        """
        Save structure to CIF files
        """
        os.makedirs(self.cif_dirs, exist_ok=True)

        structure = data["structure"]

        cif_path = os.path.join(
            self.cif_dirs,
            f"{data['material_id']}.cif"
        )

        with open(cif_path, "w") as f:
            f.write(structure)

    def _obtain_data(self):
        dataset_path = [
            path for path in self.DIRS.iterdir()
            if path.is_file() and path.suffix.lower() == ".csv"
        ]

        perov5_data = pd.concat(
            [pd.read_csv(path) for path in dataset_path],
            ignore_index=True
        )

        perov5_data.rename(columns={'cif': 'structure', 'formula': 'formula_pretty'}, inplace=True)

        perov5_dict = perov5_data.to_dict(orient="records")

        return perov5_dict
    
    def processing_perov5(self):
        docs = self._obtain_data()
        for doc in tqdm(
            docs,
            desc="Processing query results",
            unit="material"
        ):
    
            self.save_cif(doc)
            self.save_json(doc)
