"""Builds the community models."""

from os.path import isfile
import micom
from micom import Community
import pandas as pd
from tqdm import tqdm
from multiprocessing import Process

micom.logger.file_logger("micom.log")
max_procs = 6
processes = []

taxonomy = pd.read_csv("species.csv").dropna(subset=["agora_id"])
taxonomy["file"] = taxonomy.agora_id.apply(
    lambda ids: ["agora/" + i + ".xml" for i in ids.split("|")])
taxonomy.name = taxonomy.name.replace(r"[^A-Za-z0-9_\s]", "", regex=True)
taxonomy.name = taxonomy.name.replace(r"\s+", "_", regex=True)
assert not taxonomy.name.str.contains(" ").any()
taxonomy = taxonomy.rename(columns={
    "id": "samples",
    "name": "id",
    "reads": "abundance"
})


def build_and_save(tax, filename):
    com = Community(tax, id=s, progress=False)
    com.to_pickle(filename)


def consume():
    global processes
    for p in processes:
        p.join()
    processes = []


samples = tqdm(taxonomy.samples.unique(), unit="sample(s)")
for s in samples:
    filename = "models/" + s + ".pickle"
    if isfile(filename):
        continue
    tax = taxonomy[taxonomy.samples == s]
    p = Process(target=build_and_save, args=(tax, filename))
    p.start()
    processes.append(p)
    if len(processes) >= max_procs:
        consume()
consume()
