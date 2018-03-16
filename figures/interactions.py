"""Analyze interactions between taxa."""

import pandas as pd
import networkx as nx
import nxviz
import matplotlib.pyplot as plt
import numpy as np

ko = pd.read_csv("../knockouts.csv", index_col=0)
ko["knocked"] = ko.index
ko = ko.melt(id_vars=["knocked", "sample"], var_name="species",
             value_name="change").dropna()

pos = ko[ko.change.abs() > 0.01]
pos["change"] = pos.change.abs()
graph = nx.from_pandas_dataframe(pos,
                                 "knocked", "species", "change")

nx.write_gexf(graph, "interactions.gexf")

ma = ko.change.abs().max()
circ = nxviz.CircosPlot(graph, edge_color="change", rotate=True,
                        node_labels=True,
                        nodeprops={"size": 10, "radius": 0.01})
circ.figure.set_size_inches(16, 12)
circ.figure.set_dpi(80)
circ.draw()
plt.tight_layout(rect=(0.1, 0.15, 0.75, 0.85))
plt.savefig("circos.png")
