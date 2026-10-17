from config import mpAPI
from mp_api.client import MPRester
import os, logging, json, time

class queryOrigin():
    def __init__(self, JSON_PATH="QUERY/QUERY.json",
                 ORIGIN_PATH="QUERY/ORIGIN.json",
                 TASK_PATH="QUERY/ELECTRONIC_ORIGIN.json",
                 LOG_PATH="QUERY/LOG_QUERY/queryOrigin.log",
                 BATCH_QUERY = 9999):
        
        self.json_path = JSON_PATH
        self.origin_path = ORIGIN_PATH
        self.log_path = LOG_PATH
        self.task_path = TASK_PATH

        docs = self.load_json(self.json_path)
        self.docs = {doc["material_id"]: doc for doc in docs}

        self._setup_logger()
        self.BATCH_QUERY = BATCH_QUERY
    
    def load_json(self, path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

    def save_json(self, metadata, path, update=True):

        if update and os.path.exists(path):
            with open(path, "r") as f:
                try:
                    results = json.load(f)
                    if not isinstance(results, list):
                        results = [results]
                except json.JSONDecodeError:
                    results = []
        else:
            results = []

        if isinstance(metadata, list):
            results.extend(metadata)
        else:
            results.append(metadata)

        with open(path, "w") as f:
            json.dump(results, f, indent=2)

    def batch_ID(self, mpids):
        """Split mpids into chunks of at most self.BATCH_QUERY."""
        batches = [mpids[i:i + self.BATCH_QUERY]
                   for i in range(0, len(mpids), self.BATCH_QUERY)]
        return batches

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
        
    def _query_origin_batch(self, batch):
        with MPRester(mpAPI) as mpr:
            return mpr.materials.summary.search(
                material_ids=batch,
                fields=["material_id", "origins"]
            )
 
    def _query_task_batch(self, batch):
        with MPRester(mpAPI) as mpr:
            return mpr.materials.tasks.search(
                task_ids=batch,
                fields=["task_id", "run_type", "calc_type"]
            )


    def get_origin(self, mpids):
        mpids = mpids if isinstance(mpids, list) else [mpids]
        total = len(mpids)
        metadata_origin = []
        docs = []
 
        start_total = time.time()
        self._log(f"Starting Origin Query for {total} MPIDs...")
 
        batches = self.batch_ID(mpids)
        for batch_num, batch in enumerate(batches, start=1):
            batch_docs = self._timed(
                f"Querying origins (batch {batch_num}/{len(batches)})",
                self._query_origin_batch,
                batch
            )
            if batch_docs is None:
                self._log(
                    f"Batch {batch_num}/{len(batches)} failed -- "
                    f"skipping its {len(batch)} MPIDs"
                )
                continue
            docs.extend(batch_docs)
 
        obtained_ids = set()
        for doc in docs:
            metadata = {origin.name: origin.task_id for origin in doc.origins}
            metadata["material_id"] = doc["material_id"]
            metadata_origin.append(metadata)
            obtained_ids.add(doc["material_id"])
 
        missing_ids = [m for m in mpids if m not in obtained_ids]
        for mpid in missing_ids:
            self._log(f"{mpid} origin not available")
 
        obtained = len(obtained_ids)
        self._log(f"Finished Origin Query for all {total} MPIDs")
        self._log(
            f"Origin query: obtained {obtained}/{total} "
            f"({obtained/total:.1%}), missing {total - obtained}"
        )
 
        self.save_json(metadata_origin, self.origin_path, False)
 
        total_elapsed = time.time() - start_total
        self._log(f"Total execution time: {total_elapsed:.3f} s")
        return metadata_origin
 
    def get_task(self, target="electronic_structure", data_origins=None):
        if data_origins is None:
            data_origins = self.load_json(self.origin_path)
 
        start_total = time.time()
        self._log("Starting Task Origin Query...")
 
        task_to_material = {
            str(origin[target]): origin["material_id"]
            for origin in data_origins
            if target in origin
        }
 
        no_target = [
            origin["material_id"] for origin in data_origins
            if target not in origin
        ]
        for mpid in no_target:
            self._log(f"{mpid} has no '{target}' task")
 
        task_ids = list(task_to_material.keys())
        total = len(task_ids)
        metadata_task = []
        docs = []
 
        batches = self.batch_ID(task_ids)
        for batch_num, batch in enumerate(batches, start=1):
            batch_docs = self._timed(
                f"Querying tasks (batch {batch_num}/{len(batches)})",
                self._query_task_batch,
                batch
            )
            if batch_docs is None:
                self._log(
                    f"Batch {batch_num}/{len(batches)} failed -- "
                    f"skipping its {len(batch)} task IDs"
                )
                continue
            docs.extend(batch_docs)
 
        obtained_ids = set()
        for doc in docs:
            tid = str(doc.task_id)
            obtained_ids.add(tid)
            metadata = {
                "material_id": task_to_material[tid],
                "task_id": tid,
                "run_type": doc.run_type.value,
                "calc_type": doc.calc_type.value,
            }
            metadata_task.append(metadata)
 
        missing_task_ids = [t for t in task_ids if t not in obtained_ids]
        for tid in missing_task_ids:
            self._log(f"{task_to_material[tid]} task not available")
 
        obtained = len(obtained_ids)
        self._log(f"Finished Task Query for all {total} task IDs")
        if total:
            self._log(
                f"Task query: obtained {obtained}/{total} "
                f"({obtained/total:.1%}), missing {total - obtained}"
            )
 
        self.save_json(metadata_task, self.task_path, False)
 
        total_elapsed = time.time() - start_total
        self._log(f"Total execution time: {total_elapsed:.3f} s")
        return metadata_task

if __name__ == "__main__":
    qO = queryOrigin()
    mpids = list(qO.docs.keys())
    metadata_origin = qO.get_origin(mpids)
    metadata_task = qO.get_task(data_origins=metadata_origin)
    

