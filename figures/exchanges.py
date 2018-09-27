import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from adjustText import adjust_text

sample_keep = ["run_accession", "subset", "status", "type"]

SCFAs = {"butyrate": "EX_but_", "acetate": "EX_ac_", "propionate": "EX_ppa_"}


def box_jitter(x, y, **kwargs):
    sns.boxplot(x=x, y=y, **kwargs)
    sns.stripplot(x=x, y=y, **kwargs)


def export_rates_plot(fluxes, groups, samples):
    dfs = []
    for name, filt in groups.items():
        df = fluxes[fluxes.reaction.str.contains(filt)]
        res = samples.copy()
        df = df.groupby(["sample", "compartment"]).tot_flux.sum().reset_index()
        res["production"] = df[df.tot_flux < 0].groupby("compartment").tot_flux.sum()
        res["consumption"] = df[df.tot_flux > 0].groupby("compartment").tot_flux.sum()
        res["net"] = df[df.tot_flux > 0].groupby("compartment").tot_flux.sum()
        res = res.melt(id_vars=["sample", "status", "type", "subset"],
                       var_name="class", value_name="flux")
        dfs.append(res)
    fluxes = pd.concat(dfs)
    fluxes.status[fluxes.status == "ND"] = ""
    fluxes["name"] = fluxes.status + " " + fluxes.type.fillna("")
    print(fluxes)
    grid = sns.FacetGrid(fluxes, row="class", col="subset")
    g = grid.map(box_jitter, "name", "flux")
    return g


media = pd.read_csv("../results/minimal_media.csv", index_col=0).fillna(0.0)
media["sample"] = media.index
media = media.melt(id_vars="sample", var_name="reaction", value_name="flux")
metabolites = pd.read_csv("../results/metabolites.csv", index_col=0)
media["id"] = media.reaction.str.lstrip("EX_")
media = pd.merge(media, metabolites, on="id")
samples = pd.read_csv("../recent.csv")[sample_keep]
samples = samples.rename(columns={"run_accession": "sample"})
media = pd.merge(media, samples, on="sample")

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

samples = pd.read_csv("../recent.csv")[
    ["run_accession", "status", "subset", "type"]]
genera = pd.read_csv("../genera.csv")[["samples", "name", "reads"]]
totals = genera.groupby("samples").reads.sum().reset_index().reads
genera["relative"] = genera.reads / totals[genera.samples].values
fluxes = pd.merge(fluxes, genera, left_on=["sample", "compartment"],
                  right_on=["samples", "name"])
fluxes["tot_flux"] = fluxes.flux * fluxes.relative
samples = samples.rename(columns={"run_accession": "sample"})
samples.index = samples["sample"]
fig = plt.figure(figsize=(16, 12))
plt.tight_layout()
export_rates_plot(fluxes, SCFAs, samples)
fig.savefig("scfas.svg")
plt.close()


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
