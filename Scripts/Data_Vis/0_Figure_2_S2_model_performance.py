#!/usr/bin/env python3
############################################################################
# Figure 2
############################################################################

import os
import pandas as pd
import matplotlib.pyplot as plt

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

# Single environment RF (after feature selection) model performances
snp = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_FS.txt", sep="\t")
pav = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_FS.txt", sep="\t")
cnv = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_FS.txt", sep="\t")
pc = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t")

snp = snp.set_index('new_cond').loc[pc['new_cond']]
pav = pav.set_index('new_cond').loc[pc['new_cond']]
cnv = cnv.set_index('new_cond').loc[pc['new_cond']]

# Get the difference in performance compared to PCs
snp_diff = snp.r2_test - pc.set_index("new_cond").r2_test
pav_diff = pav.r2_test - pc.set_index("new_cond").r2_test
cnv_diff = cnv.r2_test - pc.set_index("new_cond").r2_test

# Plot Figure 2
fig, ax = plt.subplots(1, 1, sharey=True, figsize=(5.5, 7.5))
ax.barh(y=pc["new_cond"], width=pc["r2_test"], color="lightgray")
ax.set_xlabel("PCs Test R2")
# ax.axvline(x=0, color="black", linestyle="--")
ax.scatter(x=snp.r2_test, y=snp.index, color="tomato", label="SNP")
ax.scatter(x=pav.r2_test, y=pav.index, color="deepskyblue", label="PAV")
ax.scatter(x=cnv.r2_test, y=cnv.index, color="mediumorchid", label="CNV")
# ax.scatter(x=snp_diff, y=snp_diff.index, color="tomato", label="SNP") # 2B
# ax.scatter(x=pav_diff, y=pav_diff.index, color="deepskyblue", label="PAV") # 2B
# ax.scatter(x=cnv_diff, y=cnv_diff.index, color="mediumorchid", label="CNV") # 2B
ax.set_xlabel("Difference in test R2 compared to PCs")
ax.legend()
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_2/Figure_2_v2.pdf")
plt.close()

# How many envs were SNPs, PAVs, and CNVs better than PCs?
sum(snp_diff > 0) # 23
sum(pav_diff > 0) # 21
sum(cnv_diff > 0) # 21

# How many envs did CNVs perform better than SNPs and PAVs?
sum((cnv.r2_test - snp.r2_test) > 0) # 18
sum((cnv.r2_test - pav.r2_test) > 0) # 18
cnv.loc[(cnv.r2_test > pav.r2_test) & (cnv.r2_test > snp.r2_test) & \
        (cnv.r2_test > pc.set_index('new_cond').r2_test),:].shape[0] # 13

# How many envs did SNPs perform better than PAVs and CNVs?
snp.loc[(snp.r2_test > pav.r2_test) & (snp.r2_test > cnv.r2_test) & \
        (snp.r2_test > pc.set_index('new_cond').r2_test),:].shape[0] # 12

# How many envs did PAVs perform better than SNPs and CNVs?
pav.loc[(pav.r2_test > snp.r2_test) & (pav.r2_test > cnv.r2_test) & \
        (pav.r2_test > pc.set_index('new_cond').r2_test),:].shape[0] # 7

# How many envs did PCs perform better than SNPs, PAVs, and CNVs?
pc.set_index('new_cond').loc[(pc.set_index('new_cond').r2_test > snp.r2_test) &\
        (pc.set_index('new_cond').r2_test > pav.r2_test) & \
        (pc.set_index('new_cond').r2_test > cnv.r2_test),:].shape[0] # 3

############################################################################
# Figure S2
############################################################################
h2 = pd.read_csv("Scripts/Data_Vis/Section_2/Heritability_h2_H2_sommer.csv")
h2 = pd.concat([h2.set_index("Conditions"), pc.set_index("Y")["new_cond"]],
               axis=1, ignore_index=False).reset_index()

h2.set_index("new_cond", inplace=True)
h2.sort_values("h2", ascending=False, inplace=True)
pc.set_index("new_cond", inplace=True)

fig, ax = plt.subplots(1, 1, sharey=True, figsize=(5.5, 7.5))
ax.barh(y=h2.index, width=h2["h2"], color="white", edgecolor="black")
ax.set_xlabel("Narrow-sense heritability (bars)")
ax.scatter(x=pc.r2_test, y=pc.index, color="lightgray", label="PC")
ax.scatter(x=snp.r2_test, y=snp.index, color="tomato", label="SNP")
ax.scatter(x=pav.r2_test, y=pav.index, color="deepskyblue", label="PAV")
ax.scatter(x=cnv.r2_test, y=cnv.index, color="mediumorchid", label="CNV")
ax.legend(title="Test R2")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_2/Figure_S2_h2_explained.pdf")
plt.close()
