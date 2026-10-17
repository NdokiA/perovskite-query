import os, json, argparse, ast, logging, time
from tqdm import tqdm
from pymatgen.core import Structure
from config import mpAPI
from mp_api.client import MPRester

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

class queryStructure():
    def __init__(self, fields=FIELDS, CIF_DIRS = "TEST_CIF", JSON_PATH = "TEST_QUERY.json",
                 BATCH_QUERY = 9999, logger = None):
        """
        Query structures from Materials Project based on its MPIDs
        structural_ID   :   str or list(str)
            MPIDs for a given structure, enter using argument parser
        fields          :   list(str), optional
            fields of the projected output for MP queries on each structure,
            default is listed
        CIF_DIRS        :   str, optional
            directory path for cif files of the queried structure
        JSON_PATH       :   str, optional
            file path for the json files. 
        """

        self.fields = fields
        self.cif_dirs = CIF_DIRS
        self.json_path = JSON_PATH
        self.logger  = logger if logger is not None else logging.getLogger(__name__)
        self._setup_logger()


    def _setup_logger(self):
        """
        Configure a logger that writes timestamped action logs to `self.log_path` and
        echoes them to console.

        Uses a unique logger name per instance so that multiple `queryPerovskite`
        instances pointed at different `LOG_PATH`s don't accumulate each other's
        file handlers.
        """
        log_dir = os.path.dirname(self.log_path)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)

        self.logger = logging.getLogger(f"{__name__}.queryPerovskite.{id(self)}")
        self.logger.setLevel(logging.INFO)
        self.logger.propagate = False

        if not self.logger.handlers:
            formatter = logging.Formatter("%(asctime)s\t\t%(message)s", datefmt="%Y-%m-%d %H:%M:%S")

            file_handler = logging.FileHandler(self.log_path)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(formatter)
            self.logger.addHandler(stream_handler)

    def _log(self, action):
        """
        Log a single action with a timestamp, formatted as `time<TAB><TAB>action`.

        Parameters
        ----------
        action : str
            Description of the action being logged.
        """
        self.logger.info(action)

    def _timed(self, name, func, *args, **kwargs):
        """
        Execute a function, log its execution time, and return the result.
        """
        start = time.perf_counter()
        try:
            result = func(*args, **kwargs)
            elapsed = time.perf_counter() - start
            self._log(f"{name} completed in {elapsed:.3f} s")
            return result
        
        except Exception as e:
            elapsed = time.perf_counter() - start
            self._log(f"ERROR in {name} after {elapsed:.3f} s: {e}")
            return None

    def save_cif(self, data):
        """
        Save structure to CIF files
        """
        os.makedirs(self.cif_dirs, exist_ok=True)
        structure = data.get("structure")
        if structure is None:
            raise ValueError("No structure available")
        if isinstance(structure, dict):
            structure = Structure.from_dict(structure)
        
        cif_path = os.path.join(self.cif_dirs,f"{data['material_id']}.cif")
        with open(cif_path, "w") as f:
            f.write(structure.to(fmt="cif"))
    
    def save_json(self, data, path = None):
        """
        Save metadata to JSON
        """
        if path is None:
            path = self.json_path
        metadata = data.copy()
        metadata.pop("structure", None)

        if os.path.exists(path):
            with open(path, "r") as f:
                try:
                    results = json.load(f)
                    if not isinstance(results, list):
                        results = [results]
                except json.JSONDecodeError:
                    results = []
        else:
            results = []
        
        results.append(metadata)
        with open(path, "w") as f:
            json.dump(results, f, indent=2)
            
    def query(self, ID):
        """
        Main function for querying MPIDs
        """
        IDs = ID if isinstance(ID, list) else [ID]
        with MPRester(mpAPI) as mpr:
            docs = mpr.materials.summary.search(
                material_ids=IDs,
                fields = self.fields,
            )
        if len(docs) == 0:
            self.logger.error(f"No material found for {IDs}")

        return docs

    def batch_ID(self, mpids):
        batches = [mpids[i:i+self.BATCH_QUERY] 
                   for i in range(0, len(mpids), self.BATCH_QUERY)
                    ]
        return batches
    
    def query_ID(self, mpids):
        """
        Run the full verification pipeline (oxidation-state check, CrystalNN BX6
        connectivity check, tolerance-factor calculation) on each candidate MPID and
        save each result to its respective JSON output.
        """
        total_start = time.perf_counter()
        if not mpids:
            self._log("No candidate MPIDs provided. Aborting pipeline")
            return
        
        batches = self._timed(f"Batching queried structures", self.batch_ID, mpids)
        self._log(f"Querying in {len(batches)} batch(es) of up to {self.BATCH_QUERY} MPIDs")

        docs = []
        self._log(f"Starting verification pipeline for {len(mpids)} MPIDs")
        for batch_num, batch in enumerate(batches, start = 1):
            batch_docs = self._timed(
                f"Querying MPIDs (batch {batch_num}/{len(batches)})",
                self.query,
                batch
            )

            if batch_docs is None:
                self._log(
                    f"Batch {batch_num}/{len(batches)} failed -- "
                    f"skipping its {len(batch)} MPIDs"
                )
                continue
            docs.extend(batch_docs)


        self._timed("Processing queried structures", self.processing_query, docs)
        self._log(f"Finished Querying. Total downloaded files: {len(docs)} of requested {len(mpids)}")

        total_elapsed = time.perf_counter() - total_start
        self._log(f"Total execution time: {total_elapsed:.3f} s")
        
    def processing_query(self, docs):
        for doc in tqdm(
            docs,
            desc="Processing query results",
            unit="material"
        ):
            data = doc.model_dump(exclude={"fields_not_requested"})
    
            # Convert MPID object to string
            data["material_id"] = doc.material_id.string
    
            self.save_cif(data)
            self.save_json(data)


