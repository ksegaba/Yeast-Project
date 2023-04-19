"""
Script to create additional training and testing sets to train models with.
1st method: The original test set is 1/6 of the instances. Thus, I will make 
            5 more test sets with the remaining 5/6.
2nd method: Classify each instance into a subpopulation category based on 
            kinship
sets. 
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.model_selection import StratifiedKFold
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from scipy.spatial.distance import cdist

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018")

pheno = pd.read_csv("pheno.csv") # original diploid fitness data
test = pd.read_csv("Test.txt", header=None) # test set I've been using for all models
kin = pd.read_csv("kinship.csv", index_col=0) # kinship
eco = pd.read_excel("Peter_2018_Supplementary_Tables.xls", header=None) # ecological origins
eco = eco.iloc[3:,:] # drop extra rows
eco.columns = eco.iloc[0,:] # set column names
eco = eco.iloc[1:,:] # drop first row
eco.reset_index(drop=True, inplace=True)
eco.dtypes

### Split the training instances into 5 testing sets
# See original code for Test.txt: Data_Processing/10_holdout_test_stratified.py
pheno["Class"] = 1
pheno.set_index("ID", inplace=True)
Tr_te = StratifiedKFold(n_splits=6, random_state=42,shuffle=True)
trait="YPACETATE"
nn = 0
for train_index, test_index in Tr_te.split(pheno.loc[:,trait],pheno.Class):
    if nn == 0: # this is the original Test.txt file
        train_ID = pheno.iloc[train_index]
        test_ID = pheno.iloc[test_index]
    # create Test_split{#}.txt files for the remaining 5 splits
    if nn != 0:
        train_ID = pheno.iloc[train_index]
        test_ID = pheno.iloc[test_index]
        test = test_ID[test_ID.Class==1].index.to_list()
        out = open(f"Test_split{nn}.txt", "w")
        for t in test:
            out.write('%s\n'%t)
        out.close()
    nn += 1

### Split by taking kinship-based subpopulations into account
# Get subpopulations by clustering kinship matrix
distortions = []
inertias = []
mapping1 = {}
mapping2 = {}
for k in range(2,20):
    # Create a subplot with 1 row and 2 columns
    fig, ax1 = plt.subplots(1, 1)
    fig.set_size_inches(7, 18)

    # The 1st subplot is the silhouette plot
    ax1.set_xlim([-1, 1]) # silhouette score ranges from [-1,1]
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(kin) + (k + 1) * 10])

    # Initialize the clusterer
    clusterer = KMeans(n_clusters=k, random_state=10).fit(kin)
    cluster_labels = clusterer.labels_
    distortions.append(sum(
            np.min(
            cdist(kin, clusterer.cluster_centers_, "euclidean"), axis=1))/kin.shape[0])
    inertias.append(clusterer.inertia_)
    mapping1[k] = sum(np.min(
            cdist(kin, clusterer.cluster_centers_, "euclidean"), axis=1)) / kin.shape[0]
    mapping2[k] = clusterer.inertia_

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(kin, cluster_labels)
    print(f"For k = {k}. The average silhouette_score is : {silhouette_avg}")

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(kin, cluster_labels)

    y_lower = 10 # margin to separate cluster silhouettes
    for i in range(k):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels==i]
        ith_cluster_silhouette_values.sort()
        size_cluster_i = len(ith_cluster_silhouette_values)
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / k)
        ax1.fill_betweenx(np.arange(y_lower, y_upper), 0, 
                          ith_cluster_silhouette_values, 
                          facecolor=color, 
                          edgecolor=color, 
                          alpha=0.7,)

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax1.set_title("Silhouette plot")
    ax1.set_xlabel("silhouette coefficient values")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
    plt.suptitle(
        f"Silhouette analysis for KMeans clustering with k = {k}",
        fontsize=14,
        fontweight="bold",)
    plt.savefig(f"/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/kinship_kmeans_silhouette_k{k}.pdf")

# Elbow subplots
fig, (ax1,ax2) = plt.subplots(1, 2)
fig.set_size_inches(10, 5)
ax1.plot(range(2,20), distortions, 'bx-')
ax1.set_xlabel('Values of K')
ax1.set_ylabel('Distortion')
ax1.set_title('The Elbow Method using Distortion')
ax2.plot(range(2,20), inertias, 'bx-')
ax2.set_xlabel('Values of K')
ax2.set_ylabel('Distortion')
ax2.set_title('The Elbow Method using Inertia')
plt.savefig(f"/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/kinship_kmeans_elbow.pdf")

# Determine if clusters are based on ecological origin


# train-test split by cluster