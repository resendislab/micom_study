"""Analyze interactions between taxa."""

import pandas as pd
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

ko = pd.read_csv("../knockouts.csv", index_col=0)
ko["knocked"] = ko.index
ko = ko.melt(id_vars=["knocked", "sample"], var_name="genus",
             value_name="change").dropna()
maxs = ko.groupby("genus").change.apply(lambda x: x.abs().max())
normalized = ko.change.values / maxs[ko.genus]
normalized.index = ko.index
ko["norm_change"] = normalized


pos = ko[ko.change.abs() > 0.01]
pos["change"] = pos.change.abs()
graph = nx.from_pandas_dataframe(pos,
                                 "knocked", "genus", "norm_change")

nx.write_gexf(graph, "interactions.gexf")

mean_change = ko.groupby(["knocked", "genus"]).norm_change.mean().reset_index()
change_mat = mean_change.pivot("knocked", "genus", "norm_change").replace(np.nan, 0)
sns.clustermap(change_mat, cmap="seismic", figsize=(20, 20), vmin=-1, vmax=1)
plt.savefig("mean_change.png")
plt.close()
