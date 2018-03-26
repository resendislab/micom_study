import pandas as pd
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext

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
