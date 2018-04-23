import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from adjustText import adjust_text

sample_keep = ["run_accession", "subset", "status", "type"]


def export_rates_plot(data, subset, filter):
    subset = data[data.subset == subset]
    filtered = subset[subset.name.str.contains(filter)]
    log_fluxes = filtered.groupby(["sample", "type"]).flux.apply(
        lambda x: x.abs().sum() + 1e-6
    ).reset_index()
    pl = sns.stripplot(x="type", y="flux", data=log_fluxes,
                       jitter=True, order=["CTRL", "metformin-", "metformin+"])
    return pl


media = pd.read_csv("../results/minimal_media.csv", index_col=0).fillna(0.0)
media["sample"] = media.index
media = media.melt(id_vars="sample", var_name="reaction", value_name="flux")
metabolites = pd.read_csv("../results/metabolites.csv", index_col=0)
media["id"] = media.reaction.str.lstrip("EX_")
media = pd.merge(media, metabolites, on="id")
samples = pd.read_csv("../recent.csv")[sample_keep]
samples = samples.rename(columns={"run_accession": "sample"})
media = pd.merge(media, samples, on="sample")

fig = plt.figure(figsize=(16, 3))
plt.tight_layout()
combinations = [("MHD", "acetate"), ("MHD", "butyrate"),
                ("SWE", "acetate"), ("SWE", "butyrate")]
for i, comb in enumerate(combinations):
    ax = plt.subplot(1, 4, i+1)
    plt.yscale("log")
    export_rates_plot(media, comb[0], comb[1])
    plt.xlabel("")
    plt.ylabel("")
fig.savefig("exports.svg")
plt.close()

mat = media.pivot("id", "sample", "flux")
mat = mat.apply(lambda x: x / x.abs().max(), axis=1)
g = sns.clustermap(mat, cmap="seismic", figsize=(40, 42))
plt.savefig("media.png")
plt.close()

fluxes = pd.read_csv("../results/min_media_fluxes.csv")
fluxes = fluxes.melt(id_vars=["sample", "compartment"], var_name="reaction",
                     value_name="flux")
fluxes = fluxes[fluxes.reaction.str.startswith("EX_") &
                (fluxes.compartment != "medium")].dropna()
fluxes["taxa"] = fluxes.compartment + "_" + fluxes["sample"]
mat = fluxes.pivot("taxa", "reaction", "flux").fillna(0.0)
taxa = mat.index.str.split("_").str[0]
tsne = TSNE(n_components=2).fit_transform(mat)
tsne = pd.DataFrame(tsne, columns=["x", "y"], index=mat.index)
tsne["taxa"] = taxa
g = sns.FacetGrid(tsne, hue="taxa", size=10, aspect=1)
gm = g.map(plt.scatter, "x", "y", alpha=0.25)
means = tsne.groupby(taxa).agg("mean").reset_index()
texts = means.apply(lambda df: plt.text(df.x, df.y, df.taxa, alpha=0.5),
                    axis=1)
texts = adjust_text(texts, force_text=(0.02, 0.5),
                    arrowprops=dict(arrowstyle='-|>', alpha=0.5))
plt.savefig("individual_media.png", dpi=200)
plt.close()
