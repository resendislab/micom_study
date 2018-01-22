"""Builds the community models."""

from os.path import isfile
import micom
from micom import Community
import pandas as pd
from micom.workflows import workflow

micom.logger.file_logger("micom.log")
max_procs = 16

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


def build_and_save(args):
    s, tax = args
    filename = "models/" + s + ".pickle"
    if isfile(filename):
        return
    com = Community(tax, id=s, progress=False)
    com.to_pickle(filename)


samples = taxonomy.samples.unique()
args = [(s, taxonomy[taxonomy.samples == s]) for s in samples]
workflow(build_and_save, args, max_procs)
