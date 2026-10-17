import os, logging, json, time
from tqdm.auto import tqdm
from config import mpAPI
from .queryStructure import queryStructure as qS
from mpcontribs.client import Client
client = Client(apikey=mpAPI, project="Multinary_Oxides")

class queryPerovskite_Contribs(qS):
    def __init__(self, CIF_DIRS = "QUERY/CIF", JSON_PATH = "QUERY/QUERY.json",
                 LOG_PATH = "QUERY/LOG_QUERY/queryContrib.log"):
        
        self.log_path = LOG_PATH
        super().__init__(CIF_DIRS=CIF_DIRS, JSON_PATH=JSON_PATH)

    def query_data(self):
        query = client.download_contributions(
            query = {"data__perovskite__exact": "Yes"},
            include=["structures"]
        )
        return query

    def obtain_json(self, query):

        query_dict = {
            "material_id": query["id"],
            "formula_pretty": query["formula"],
            "tolerance_t": query["data"]["tolerance"]["t"]["value"],
            "tolerance_tau": query["data"]["tolerance"]["tau"]["value"],
            "formation_energy": query["data"]["dH"]["formation"]["value"],
            "decomposition_energy": query["data"]["dH"]["decomposition"]["value"],
            "is_perovskite":  query["data"]["perovskite"] == "Yes",
            "structure": query["structures"][0]
            }

        return query_dict

    def query_ID(self):
        start = time.perf_counter()
        self._log(f"Querying Contrib MP Perovskite")

        query = self._timed("Querying MP-Contribs Data", self.query_data)
        self._log("Saving Data...")
        for doc in tqdm(query, desc = "Processing query results", unit = "material"):
            doc_dict = self.obtain_json(doc)
            self.save_cif(doc_dict)
            self.save_json(doc_dict)

            elapsed = time.perf_counter() - start
        self._log(f"Query completed in {elapsed:.3f} s")

if __name__ == "__main__":
    qPC = queryPerovskite_Contribs(
        CIF_DIRS="CQUERY/CIF", 
        JSON_PATH="CQUERY/QUERY.json",
        LOG_PATH="CQUERY/LOG_QUERY/queryContrib.log"
    )
    qPC.query_ID()