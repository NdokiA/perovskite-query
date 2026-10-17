import os, logging, json, time
from config import mpAPI
from .queryStructure import queryStructure as qS
from mp_api.client import MPRester

def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)
    
class queryPerovskite_Formula(qS):
    def __init__(self, FORMULA = None,
                 CIF_DIRS="QUERY/CIF", JSON_PATH="QUERY/QUERY.json",
                 LOG_PATH="QUERY/LOG_QUERY/queryFormula.log",
                 QUERIED = None,
                 BATCH_QUERY = 9999):
        """
        Query and structurally verify candidate perovskite structures from Materials Project.

        Candidates are first identified based on formula_anonymous keyword and various others
        """
        
        self.log_path = LOG_PATH
        self.formula = FORMULA 
        self.queried = set(QUERIED)
        self.BATCH_QUERY = BATCH_QUERY

        super().__init__(BATCH_QUERY=BATCH_QUERY, CIF_DIRS=CIF_DIRS, JSON_PATH=JSON_PATH)
        

    def obtain_ID(self):
        """
        Query candidate perovskite MPIDs manually based on the formula_anonymous
        from previous run (robocryst Query) or initialized in class

        Returns
        -------
        list(str)
            Union of MPIDs
        """
        self._log(f"Starting formula keyword search for formula: {self.formula}")

        with MPRester(mpAPI) as mpr:
            formula_docs = mpr.materials.summary.search(
                formula  = list(self.formula),
                energy_above_hull = (0,0.4),
                fields = ["material_id"],
                chunk_size = 1000, #set to 1000 for prod, 1 for test
                num_chunks = None, #set to None for prod, 1 for test

            )

            candidate_mpids = {d.material_id.string for d in formula_docs} - self.queried
        
        self._log(f"Formula search returned {len(candidate_mpids)} candidate MPIDs")
 
        self._log(f"Total {len(candidate_mpids)} candidate MPIDs successfully obtained")
        return list(candidate_mpids)


if __name__ == "__main__":
    query = queryPerovskite_Formula()
    IDs = query.obtain_ID()
    if not IDs:
        raise ValueError("query() called with an empty list of material IDs. Nothing to query.")
    query.query_ID(IDs)
