#!/bin/python
import datatable as dt
import pandas as pd
import timeit
from nancorrmp.nancorrmp import NaNCorrMp
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Read in genotype data
print("Reading in data...")
geno = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv")
geno = geno.to_pandas()
geno = geno.set_index("ID")

# Calculate correlation between isolates
print("Calculating genotype pearson correlation coefficients...")
geno_corr = NaNCorrMp.calculate(geno, n_jobs=4, chunks=10000)

# Generate heatmap
print("Generating heatmap...")
fig = sns.clustermap(geno_corr, figsize=(12,12), cmap="RdBu_r", center=0, yticklabels=False, xticklabels=False, dendrogram_ratio=(.1, .1))
ax = fig.ax_heatmap
ax.set_xlabel("Isolates")
ax.set_ylabel("Isolates")
ax.tick_params(right=False, bottom=False)
plt.tight_layout(rect=[0, 0.7, 0, 0.7]) # left, bottom, right, top
plt.savefig("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/geno_isolates_corr.pdf", bbox_inches = "tight")
