"""Analyze elasticities for reference samples."""

import pandas as pd
import networkx as nx
import nxviz
import seaborn as sns
import matplotlib.pyplot as plt

samples = ["ERR260275", "ERR260214", "ERR260174"]


elast = []
for sa in samples:
    e = pd.read_csv("../results/" + sa + ".csv")
    e["id"] = sa
    maxs = e.groupby("reaction").elasticity.apply(lambda x: x.abs().max())
    e = e.assign(norm_elasticity=e.elasticity.values /
                 (maxs[e.reaction].values + 1e-6))
    elast.append(e)
elast = pd.concat(elast)
elast = elast[elast.direction != "zero"]
elast.effector = elast.effector.str.replace("EX_", "").str.replace("_m", "")

for sa in samples:
    sns.kdeplot(elast[elast.id == sa].elasticity, bw=4,
                shade=True, label=sa)
plt.legend()
plt.xlabel("elasticity [a.u.]")
plt.ylabel("density")
plt.savefig("elast_densities.svg")
plt.close()

for sa in samples:
    e = elast[elast.id == sa].copy()
    e.elasticity = e.elasticity.abs()
    e = e[(e.elasticity > 0.5) & (e.norm_elasticity > 0.5)]
    graph = nx.from_pandas_edgelist(e, source="effector", target="reaction",
                                    edge_attr="elasticity")
    for idx, _ in graph.nodes(data=True):
        if idx.startswith("EX_"):
            cl = "exchange flux"
        elif idx[0].isupper():
            cl = "abundance"
        else:
            cl = "diet"
        graph.node[idx]["class"] = cl

    circos = nxviz.CircosPlot(graph, node_labels=True, rotate_labels=True,
                              edge_color="elasticity", node_color="class",
                              node_grouping="class", node_order="class",
                              figsize=(20, 18))
    circos.draw()
    plt.tight_layout(rect=[0.1, 0.1, 0.8, 0.9])
    plt.savefig(sa + ".png")
    plt.close()

but_ac = elast[(elast.reaction == "EX_but_m") | (elast.reaction == "EX_ac_m")]
