
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import matplotlib.patheffects as path_effects
from mpl_toolkits.axes_grid1 import make_axes_locatable

ab1data = pd.read_csv(snakemake.input[0])

dicto = {"RS" : "Blue", "NR": "Red"}
df2=ab1data.replace({"Group": dicto})

timepointlist = ["d0", "d2", "d4", "d6"]
timepointlist2 = ["Day 0", "Day 2", "Day 4", "Day 6"]

def get_PCs(ab1_data, timepoint):
    day0 = ab1_data[ab1_data["Timepoint"] == timepoint]
    X = day0["PC1"]
    Y = day0["PC2"]
    c = day0["Group"]
    return X, Y, c

max_lims = 15000

if 'ab1' in snakemake.input[0]:
    chartlabel = "AB1"
else:
    chartlabel = "Renca"

rows, cols = 2, 2
fig = plt.figure(figsize=(cols * 2, rows * 2.45))

for row in range(rows):
    for col in range(cols):
        index = row * cols + col + 1

        ax = plt.subplot(rows, cols, index, aspect=1, frameon=True)
        ax.tick_params(which="both", direction="in")
        ax.tick_params(which="both", right=True)
        ax.set_axisbelow(True)

        ax.set_xlim(-max_lims, max_lims)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax.set_xticklabels([])

        ax.set_ylim(-max_lims, max_lims)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1000))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax.set_yticklabels([])

        ax.grid(color=".9", linestyle="--")

        #####
        
        X, Y, c = get_PCs(df2, timepointlist[(index - 1)])
        ax.scatter(X, Y, alpha=0.5, color = c, edgecolor="None")
        
        #####

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="15%", pad=0)
        cax.get_xaxis().set_visible(False)
        cax.get_yaxis().set_visible(False)
        cax.set_facecolor("black")
        cax.text(
            0.05,
            0.45,
            timepointlist2[index-1],
            size=10,
            color="white",
            ha="left",
            va="center",
            weight="bold",
        )
        cax.text(
            0.95,
            0.45,
            chartlabel,
            size=10,
            color="white",
            ha="right",
            va="center",
            weight="bold",
        )

plt.tight_layout()
plt.savefig(snakemake.output[0], dpi = 600)
