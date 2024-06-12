#!/usr/bin/env python3
############################################################################
# Figure 2
############################################################################

import os
import pandas as pd
import matplotlib.pyplot as plt

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

# Single environment RF (after feature selection) model performances
res = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_FS.txt", sep="\t")
pc = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t")

# Get the difference in performance compared to PCs
snp_diff = res.loc[(res.Data=="SNP") & (res.Alg=="RF"), ["new_cond", "r2_test"]].set_index("new_cond")
snp_diff = snp_diff.r2_test - pc.set_index("new_cond").r2_test
pav_diff = res.loc[(res.Data=="PAV") & (res.Alg=="RF"), ["new_cond", "r2_test"]].set_index("new_cond")
pav_diff = pav_diff.r2_test - pc.set_index("new_cond").r2_test
cnv_diff = res.loc[(res.Data=="CNV") & (res.Alg=="RF"), ["new_cond", "r2_test"]].set_index("new_cond")
cnv_diff = cnv_diff.r2_test - pc.set_index("new_cond").r2_test

# Plot Figure 2A and 2B
fig, ax = plt.subplots(1, 2, sharey=True, figsize=(8.75, 7.5))
ax[0].barh(y=pc["new_cond"], width=pc["r2_test"], color="skyblue") # 2A
ax[0].set_xlabel("PCs Test R2")
ax[1].scatter(x=snp_diff, y=snp_diff.index, color="tomato", label="SNP") # 2B
ax[1].scatter(x=pav_diff, y=pav_diff.index, color="deepskyblue", label="PAV") # 2B
ax[1].scatter(x=cnv_diff, y=cnv_diff.index, color="mediumorchid", label="CNV") # 2B
ax[1].set_xlabel("Difference in test R2 compared to PCs")
ax[1].legend()
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_2/Figure_2.pdf")
plt.close()
