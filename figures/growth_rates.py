import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext
import seaborn as sns
import joypy as jp
import pickle

community_mass = 200

with open("../tradeoff.pickle", "rb") as pick:
    solutions = pickle.load(pick)


def get_members(sol):
    return sol.apply(lambda s: s.members.drop("medium")).iloc[0]


rates = (solutions.groupby(["tradeoff", "sample"]).solution.
         apply(get_members).reset_index())
rates.growth_rate /= community_mass
rates["log_rates"] = np.log10(rates.growth_rate)
rates.loc[rates.growth_rate <= 0, "log_rates"] = -18

fig, axes = jp.joyplot(rates,
                       by="tradeoff", column="log_rates",
                       color="cornflowerblue")
lf = LogFormatterMathtext()
xax = axes[-1].xaxis
xax.set_ticklabels(lf(10.0**ex) for ex in xax.get_ticklocs())
plt.xlabel("growth rate [1/h]")
plt.ylabel("tradeoff")
plt.savefig("dists.svg")
plt.close()

non_zero = rates.groupby(["sample", "tradeoff"]).growth_rate. \
           apply(lambda x: sum(x > 1e-6) / len(x)).reset_index(name="non_zero")
g = sns.boxplot("tradeoff", "non_zero", data=non_zero, color="skyblue")
plt.xlabel("tradeoff")
plt.ylabel("fraction growing")
plt.savefig("non_zero.svg")
plt.close()

pos = rates.query("growth_rate > 1e-6 and tradeoff == 0.5")
o = pos.groupby("compartments").log_rates.mean().sort_values().index
ax = sns.stripplot("log_rates", "compartments", data=pos,
                   order=o, alpha=0.33)
xax = ax.xaxis
xax.set_ticklabels(lf(10.0**ex) for ex in xax.get_ticklocs())
ax.figure.set_figwidth(5)
ax.figure.set_figheight(10)
plt.xlabel("growth rate [1/h]")
plt.ylabel("")
plt.tight_layout()
plt.savefig("gcs.svg", width=10, height=20)
plt.close()

fig = plt.figure()
community_growth = rates[rates.tradeoff == 0.5].groupby("sample"). \
                   apply(lambda df: sum(df.abundance * df.growth_rate))
sns.distplot(community_growth)
plt.savefig("community_growth.svg")
plt.close()
