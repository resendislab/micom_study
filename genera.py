import pandas as pd

keep = ["samples", "ncbi_genus_id", "kingdom", "phylum", "class", "order",
        "family", "genus", "oxygen_status", "metabolism", "gram", "type",
        "pubseed_id", "genome_size", "genes", "file", "name", "reads",
        "relative"]


def reduce_group(df):
    new = df.iloc[0,:]
    new["file"] = "|".join(df.id.apply(lambda id: f"{id}.xml"))
    return new


agora = pd.read_csv("agora.csv")
agora_genus = agora.groupby("ncbi_genus_id").apply(reduce_group)

genera = pd.read_csv("taxonomy.csv")
genera = genera[genera["rank"] == "genus"]

genus_models = pd.merge(genera, agora_genus, left_on="taxid",
                        right_on="ncbi_genus_id")
genus_models = genus_models.rename(columns={"id_x": "samples"})[keep]
genus_models.to_csv("genera.csv", index=False)
