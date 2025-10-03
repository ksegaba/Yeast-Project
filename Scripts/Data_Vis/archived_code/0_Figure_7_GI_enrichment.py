#!/usr/bin/env python3
'''
Enrichment of known genetic interactions in SHAP-based interactions
'''
# %%
import os, glob, re, swifter
import pandas as pd
import datatable as dt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import networkx as nx
from venn import venn
from tqdm import tqdm
from scipy import sparse

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

def get_unique_gp(df, Gene1="Gene1", Gene2="Gene2"):
	'''
	Get the unique gene pairs from the dataframe
	'''
	df_gp = df.apply(lambda x: set([x[Gene1], x[Gene2]]), axis=1).values # gene pairs
	df_gp = {frozenset(sorted(set)) for set in df_gp} # get unique interactions
	return df_gp


# %%
## BioGRID database genetic interactions
biogrid = dt.fread("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/BioGRID/yeast_gi_biogrid.txt").to_pandas()
biogrid = biogrid.iloc[:,[5,6,7,8,11,13,14,17,36]]
biogrid.columns = ["Systematic Name Interactor A", "Systematic Name Interactor B",
				   "Standard Name Interactor A", "Standard name Interactor B",
				   "Evidence", "Author", "PMID", "Throughput", "Organism"]
biogrid = biogrid.loc[biogrid.Organism=="Saccharomyces cerevisiae (S288c)",:]
biogrid = biogrid.loc[biogrid.Evidence.str.strip().isin(["Synthetic Growth Defect",\
			"Synthetic Lethality", "Synthetic Rescue", "Negative Genetic",\
			"Positive Genetic"]),:] # remove overexpression gene pairs
biogrid.to_csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/BioGRID/yeast_gi_biogrid_filtered.txt", sep="\t", index=False)

biogrid_gp = get_unique_gp(biogrid, "Systematic Name Interactor A", "Systematic Name Interactor B")
len(biogrid_gp) # 438546

## Costanzo 2021 benomyl genetic interactions
costanzo = pd.read_excel("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Costanzo_2021/2021_Costanzo_Data File S3_Raw interaction dataset.xlsx",\
	engine="openpyxl", sheet_name="Genome-scale_Benomyl")

# Filter the genetic interaction networks using stringent criteria
costanzo_cond = costanzo.loc[(costanzo.mean_condition_epsilon.abs() > 0.12) & \
	(costanzo.condition_p_value < 0.05),:] # benomyl interactions
costanzo_ctrl = costanzo.loc[(costanzo.mean_reference_epsilon.abs() > 0.12) & \
	(costanzo.reference_p_value < 0.05),:] # control condition interactions
costanzo_diff = costanzo_cond.loc[(costanzo.mean_differential_epsilon.abs() > 0.12) & \
	(costanzo.differential_p_value < 0.05),:] # differential interaction network

# Get the unique gene pairs from each network
costanzo_cond.insert(0, "array_gene", costanzo.apply(lambda x: x.array_orf.split("_")[0], axis=1))
costanzo_cond.insert(0, "query_gene", costanzo.apply(lambda x: x.query_orf.split("_")[0], axis=1))
costanzo_ctrl.insert(0, "array_gene", costanzo.apply(lambda x: x.array_orf.split("_")[0], axis=1))
costanzo_ctrl.insert(0, "query_gene", costanzo.apply(lambda x: x.query_orf.split("_")[0], axis=1))
costanzo_diff.insert(0, "array_gene", costanzo.apply(lambda x: x.array_orf.split("_")[0], axis=1))
costanzo_diff.insert(0, "query_gene", costanzo.apply(lambda x: x.query_orf.split("_")[0], axis=1))

costanzo_cond_gp = get_unique_gp(costanzo_cond, "query_gene", "array_gene")
costanzo_ctrl_gp = get_unique_gp(costanzo_ctrl, "query_gene", "array_gene")
costanzo_diff_gp = get_unique_gp(costanzo_diff, "query_gene", "array_gene")
len(costanzo_cond_gp) # 3472
len(costanzo_ctrl_gp) # 3417
len(costanzo_diff_gp) # 588
len(biogrid_gp.union(costanzo_cond_gp)) # 440463 total unique
len(biogrid_gp.union(costanzo_ctrl_gp)) # 440116 total unique
len(biogrid_gp.union(costanzo_diff_gp)) # 438842 total unique
len(biogrid_gp.union(costanzo_cond_gp).union(costanzo_ctrl_gp)) # 441520

costanzo_ctrl.to_csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Costanzo_2021/strict_costanzo_ctrl_interactions.tsv", sep="\t", index=False)
costanzo_cond.to_csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Costanzo_2021/strict_costanzo_benomyl_interactions.tsv", sep="\t", index=False)

## Venn diagram of the experimentally verified genetic interactions
sets = {"Costanzo Control": costanzo_ctrl_gp,
		"Costanzo Benomyl": costanzo_cond_gp,\
		# "Costanzo Differential": costanzo_diff_gp,\
		"BioGRID": biogrid_gp}
venn(sets) # total: 441520 unique GIs
plt.title("Experimentally verified genetic interactions")
plt.tight_layout()
plt.savefig("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Scripts/Data_Vis/Section_6/Experimentally_verified_GIs_venn_diagram.pdf")
plt.close()
del sets

# %%
####### Compare SHAP-based interactions to known GIs from the literature #######
# How many known GIs are identified by the SHAP-based interactions?
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction"
tosave = "Scripts/Data_Vis/Section_6/shap_interaction"

# SNP and ORF gene maps
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")

def count_GIs(model_dir):
	''' Count the total number of GIs identified by SHAP interaction scores,
	determine how many of these are experimentally verified according to BIOGRID
	or Costanzo et al., 2021, and plot SHAP interaction scores into a violin.'''
	
	# dictionary to collect GI counts
	gi_counts = {'costanzo_ctrl':{"SNP":{}, "PAV":{}, "CNV":{}},
				'costanzo_cond':{"SNP":{}, "PAV":{}, "CNV":{}},
				'costanzo_diff':{"SNP":{}, "PAV":{}, "CNV":{}},
				'biogrid':{"SNP":{}, "PAV":{}, "CNV":{}},
				'total_shap_gi':{"SNP":{}, "PAV":{}, "CNV":{}}} # number of known GIs identified by SHAP
	collect_gi = {'costanzo_ctrl':{"SNP":{}, "PAV":{}, "CNV":{}},
				'costanzo_cond':{"SNP":{}, "PAV":{}, "CNV":{}},
				'costanzo_diff':{"SNP":{}, "PAV":{}, "CNV":{}},
				'biogrid':{"SNP":{}, "PAV":{}, "CNV":{}},
				'total_shap_gi':{"SNP":{}, "PAV":{}, "CNV":{}}} # the known GIs identified by SHAP
	
	violin_data = {} # violin plot data
	points_data = {'biogrid': {}, 'costanzo_ctrl': {}, 'costanzo_cond': {}, 'costanzo_diff': {}} # stripplot data
	
	for env in target_envs:
		for data_type in ["SNP", "PAV", "CNV"]:
			shap_sum = pd.read_csv(
				glob.glob(f"{d}/{data_type}/{model_dir}/shap_interaction_scores_{data_type.lower()}_{env}_*_summed.txt")[0],
				sep="\t") # SHAP interaction scores summed across isolates
			
			shap_sum = shap_sum.loc[shap_sum.Interaction > 0,:] # drop rows with zero values
			print(env, data_type, shap_sum.shape)
			
			# add violin plot data to dictionary
			violin_data[f"{data_type}_{env}"] = shap_sum.Interaction.to_numpy()
			
			if len(shap_sum) != 0:
				# map the feature names to the gene names
				if data_type == "SNP":
					shap_sum["Gene1"] = shap_sum.Feature1.map(map_snps.set_index("snp").gene)
					shap_sum["Gene2"] = shap_sum.Feature2.map(map_snps.set_index("snp").gene)
					# drop interactions with at least one intergenic snp, we are only looking for gene-gene interactions
					shap_sum = shap_sum.loc[(shap_sum.Gene1 != "intergenic") & (shap_sum.Gene2 != "intergenic"),:]
					# drop snps with multiple gene matches
					shap_sum = shap_sum.loc[~(shap_sum.Gene1.str.contains(",")) & ~(shap_sum.Gene2.str.contains(",")),:]
				else:
					shap_sum.Feature1 = shap_sum.Feature1.apply(lambda x: re.sub("^X", "", x))
					shap_sum.Feature1 = shap_sum.Feature1.apply(lambda x: re.sub("\.", "-", x))
					shap_sum.Feature2 = shap_sum.Feature2.apply(lambda x: re.sub("^X", "", x))
					shap_sum.Feature2 = shap_sum.Feature2.apply(lambda x: re.sub("\.", "-", x))
					shap_sum["Gene1"] = shap_sum.Feature1.map(map_orfs.set_index("orf").gene)
					shap_sum["Gene2"] = shap_sum.Feature2.map(map_orfs.set_index("orf").gene)
					shap_sum["Gene1"] = shap_sum["Gene1"].fillna(shap_sum["Feature1"])
					shap_sum["Gene2"] = shap_sum["Gene2"].fillna(shap_sum["Feature2"])
				
				# exclude shap_sum rows where Gene1 == Gene2
				shap_sum = shap_sum.loc[shap_sum.Gene1 != shap_sum.Gene2]
				
				# count the number of known GIs identified by the SHAP-based interactions
				shap_gp = get_unique_gp(shap_sum, "Gene1", "Gene2")
				gi_counts['biogrid'][data_type][env]= len(shap_gp.intersection(biogrid_gp))
				gi_counts['costanzo_ctrl'][data_type][env]= len(shap_gp.intersection(costanzo_ctrl_gp))
				gi_counts['costanzo_cond'][data_type][env]= len(shap_gp.intersection(costanzo_cond_gp))
				gi_counts['costanzo_diff'][data_type][env]= len(shap_gp.intersection(costanzo_diff_gp))
				gi_counts['total_shap_gi'][data_type][env]= len(shap_gp)
				
				# collect the known GIs identified by SHAP-based interactions
				collect_gi['biogrid'][data_type][env] = shap_gp.intersection(biogrid_gp)
				collect_gi['costanzo_ctrl'][data_type][env] = shap_gp.intersection(costanzo_ctrl_gp)
				collect_gi['costanzo_cond'][data_type][env] = shap_gp.intersection(costanzo_cond_gp)
				collect_gi['costanzo_diff'][data_type][env] = shap_gp.intersection(costanzo_diff_gp)
				collect_gi['total_shap_gi'][data_type][env] = shap_gp
				
				# collect the interaction scores for the known GIs to plot onto the violinplot
				for key in ['biogrid', 'costanzo_ctrl', 'costanzo_cond', 'costanzo_diff']:
					points_data[key][f"{data_type}_{env}"] = []
					
					if len(collect_gi[key][data_type][env]) > 0:
						for g1, g2 in collect_gi[key][data_type][env]:
							points_data[key][f"{data_type}_{env}"].extend(
								shap_sum.loc[((shap_sum.Gene1 == g1) &\
											  (shap_sum.Gene2 == g2)) |\
											 ((shap_sum.Gene1 == g2) &\
											  (shap_sum.Gene2 == g1)), "Interaction"].values
							)
				
				del shap_gp
			else:
				gi_counts['biogrid'][data_type][env] = 0
				gi_counts['costanzo_ctrl'][data_type][env] = 0
				gi_counts['costanzo_cond'][data_type][env] = 0
				gi_counts['costanzo_diff'][data_type][env] = 0
				gi_counts['total_shap_gi'][data_type][env] = 0
				gi_counts['biogrid'][data_type][env] = {}
				gi_counts['costanzo_ctrl'][data_type][env] = {}
				gi_counts['costanzo_cond'][data_type][env] = {}
				gi_counts['costanzo_diff'][data_type][env] = {}
				gi_counts['total_shap_gi'][data_type][env] = {}
			
			del shap_sum#, mean_vals
	
	# plot violin of log transformed SHAP interaction scores
	colors = ['#3274A1', '#3274A1', '#3274A1', '#E1812C', '#E1812C', '#E1812C',
			  '#3B923B', '#3B923B', '#3B923B', '#C03D3E', '#C03D3E', '#C03D3E',
			  '#9472B2', '#9472B2', '#9472B2']
	
	points_data = pd.DataFrame.from_dict(points_data).stack().explode().dropna()
	points_data = points_data.reset_index()
	points_data.columns = ["model", "known_gi", "shap_interaction_score"]
	
	sns.violinplot(violin_data, inner="box", fill=False, log_scale=10)
	sns.stripplot(x='model', y='shap_interaction_score', data=points_data,
				  linewidth=.1, palette=colors, alpha=.4, log_scale=10, size=3,
				  hue='known_gi', dodge=True)
	
	plt.ylabel("SHAP Interaction Scores")
	plt.xticks(rotation=55, ha="right")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_6/shap_interaction/violin_log_summed_shap_interaction.pdf")
	plt.close()
	del violin_data
	
	return gi_counts, collect_gi


best_rf_fs_gi_counts, best_rf_fs_gi = count_GIs("best_rf_fs")
pd.DataFrame.from_dict({(i, j): best_rf_fs_gi_counts[i][j]
						for i in best_rf_fs_gi_counts.keys()
						for j in best_rf_fs_gi_counts[i].keys()}).to_csv(
	"Scripts/Data_Vis/Section_6/shap_interaction/shap_int_best_rf_fs_gi_counts.tsv", sep="\t", index=True, header=True)

# Check for overlap between GIs found by SNPs, PAVs, and CNVs
for env in target_envs:
	for key in best_rf_fs_gi.keys():
		print(env, key)
		print(f"SNP vs PAV: {len(best_rf_fs_gi[key]['SNP'][env].intersection(best_rf_fs_gi[key]['PAV'][env]))}")
		print(f"SNP vs CNV: {len(best_rf_fs_gi[key]['SNP'][env].intersection(best_rf_fs_gi[key]['CNV'][env]))}")
		print(f"PAV vs CNV: {len(best_rf_fs_gi[key]['PAV'][env].intersection(best_rf_fs_gi[key]['CNV'][env]))}")

for key in best_rf_fs_gi.keys():
	print(key)
	for env in target_envs:
		for env2 in target_envs:
			if env != env2:
				print(f"SNP {env} vs {env2}: {len(best_rf_fs_gi[key]['SNP'][env].intersection(best_rf_fs_gi[key]['SNP'][env2]))}")
				print(f"PAV {env} vs {env2}: {len(best_rf_fs_gi[key]['PAV'][env].intersection(best_rf_fs_gi[key]['PAV'][env2]))}")
				print(f"CNV {env} vs {env2}: {len(best_rf_fs_gi[key]['CNV'][env].intersection(best_rf_fs_gi[key]['CNV'][env2]))}")


# %%
'''According to the gi counts, only the following models identified known GIs:
- Costanzo_ctrl; SNP; all except YPDCAFEIN40
- Costanzo_cond; SNP; all except YPDCAFEIN40
- Costanzo_diff; SNP; all except YPDCAFEIN40 & YPDCUSO410MM
- Biogrid; SNP; all
- Biogrid; PAV; all except YPDCAFEIN50 & YPDBENOMYL500
- Biogrid; CNV; all except YPDCAFEN40 & 50

And, there is overlap (non-unique GIs) between genetic variant types and between
environments.
'''

## Identify which of the known GIs identified by SHAP involve at least one benchmark gene

# Draw the networks of the known GIs identified by SHAP-based interactions
def draw_networks(best_rf_fs_gi, key, d, model_dir, env, data_type):
	# Get the identified GIs
	gi = best_rf_fs_gi[key][data_type][env]
	
	# Load the SHAP-based interaction scores
	shap = pd.read_csv(
		glob.glob(f"{d}/{data_type}/{model_dir}/shap_interaction_scores_{data_type.lower()}_{env}_*_summed.txt")[0],
		sep="\t")
	
	# map the feature names to the gene names
	if data_type == "SNP":
		shap["Gene1"] = shap.Feature1.map(map_snps.set_index("snp").gene)
		shap["Gene2"] = shap.Feature2.map(map_snps.set_index("snp").gene)
		# drop interactions with at least one intergenic snp, we are only looking for gene-gene interactions
		shap = shap.loc[(shap.Gene1 != "intergenic") & (shap.Gene2 != "intergenic"),:]
		# drop snps with multiple gene matches
		shap = shap.loc[~(shap.Gene1.str.contains(",")) & ~(shap.Gene2.str.contains(",")),:]
	else:
		shap.Feature1 = shap.Feature1.apply(lambda x: re.sub("^X", "", x))
		shap.Feature1 = shap.Feature1.apply(lambda x: re.sub("\.", "-", x))
		shap.Feature2 = shap.Feature2.apply(lambda x: re.sub("^X", "", x))
		shap.Feature2 = shap.Feature2.apply(lambda x: re.sub("\.", "-", x))
		shap["Gene1"] = shap.Feature1.map(map_orfs.set_index("orf").gene)
		shap["Gene2"] = shap.Feature2.map(map_orfs.set_index("orf").gene)
		shap["Gene1"] = shap["Gene1"].fillna(shap["Feature1"])
		shap["Gene2"] = shap["Gene2"].fillna(shap["Feature2"])
	
	# get the known GIs identified by SHAP-based interactions
	idx = []
	for pair in gi:
		pair = list(pair)
		idx.append(shap.loc[(shap.Gene1 == pair[0]) & (shap.Gene2 == pair[1]) |\
			(shap.Gene1 == pair[1]) & (shap.Gene2 == pair[0])].index[0])
	
	identified_gi = shap.loc[idx,:]
	
	# Build the network
	G = nx.from_pandas_edgelist(identified_gi, "Gene1", "Gene2", "Interaction",
		create_using=nx.Graph())
	
	# Extract edge colors
	edge_weights = identified_gi["Interaction"].values
	norm = mcolors.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))  # Normalize color scale
	cmap = cm.Blues  # Choose colormap
	sm = cm.ScalarMappable(norm=norm, cmap=cmap)  # Create a ScalarMappable for the colorbar
	sm.set_array([])  # Required for colorbar to work
	
	# Draw the network
	fig = plt.figure(figsize=(10,10))
	gs = gridspec.GridSpec(1, 2, width_ratios=[10, 0.5]) # allocate space for colorbar
	ax = plt.subplot(gs[0]) # main plot area
	cax = plt.subplot(gs[1]) # colorbar axis
	pos = nx.spring_layout(G)  # Compute node positions
	nx.draw_networkx_edges(G, pos, edge_color=edge_weights, edge_cmap=cmap, edge_vmin=min(edge_weights), edge_vmax=max(edge_weights), ax=ax)
	nx.draw_networkx_nodes(G, pos, node_size=12, ax=ax)
	nx.draw_networkx_labels(G, pos, font_size=12, ax=ax)
	
	# color the nodes that correspond to benchmark genes in red
	if data_type == "SNP":
		benchmark_genes = map_snps.loc[map_snps[["Benomyl", "Caffeine", "CuSO4", "Sodium_meta-arsenite"]].sum(axis=1) > 0,"gene"].unique()
	if data_type in ["PAV", "CNV"]:
		benchmark_genes = map_ords.loc[map_orfs[["Benomyl", "Caffeine", "CuSO4", "Sodium_meta-arsenite"]].sum(axis=1) > 0,"gene"].unique()
	benchmark_genes = set(benchmark_genes).intersection(G.nodes()) # search for benchmark_genes values in G.nodes()
	node_colors = ["red" if node in benchmark_genes else "blue" for node in G.nodes()]
	nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=12, ax=ax)
	
	cbar = plt.colorbar(sm, cax=cax) # Add a colorbar legend
	cbar.set_label("SHAP Interaction Score")
	ax.set_title(f"{data_type} {env} {model_dir}")
	plt.tight_layout()
	
	tosave = "Scripts/Data_Vis/Section_6/shap_interaction/"
	plt.savefig(f"{tosave}/shap_interaction_network_{data_type.lower()}_{env}_{model_dir}_{key}.pdf")
	plt.close()
	
	return G


model_dir = "best_rf_fs"
for key in best_rf_fs_gi.keys():
	for env in target_envs:
		print(key, env)
		if (key == 'costanzo_ctrl' and env != 'YPDCAFEIN40') |\
			(key == 'costanzo_cond' and env != 'YPDCAFEIN40') |\
			(key == 'costanzo_diff' and env not in ['YPDCAFEIN40', 'YPDCUSO410MM']):
			draw_networks(best_rf_fs_gi, key, d, model_dir, env, 'SNP')
		# elif key == 'biogrid':
		# 	draw_networks(best_rf_fs_gi, key, d, model_dir, env, 'SNP')
		# elif (key == 'biogrid' and env not in ['YPDCAFEIN50', 'YPDBENOMYL500']):
		# 	draw_networks(best_rf_fs_gi, key, d, model_dir, env, 'PAV')
		# elif (key == 'biogrid' and env not in ['YPDCAFEIN40', 'YPDCAFEIN50']):
		# 	draw_networks(best_rf_fs_gi, key, d, model_dir, env, 'CNV')
		# biogrid and total_shap_gi take too long to plot

# %%
# Are experimentally verified GIs enriched in the lit models compared to the top non-lit models?

# %%
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction"
num_int = {}
model_dir = "best_rf_fs" # or should it be the only_... folders?
for env in target_envs:
	print(env)
	
	# Load the SHAP-based interaction scores
	snp = dt.fread(glob.glob(f"{d}/SNP/{model_dir}/shap_interaction_scores_snp_{env}_top_*_lit_genes_summed.txt")[0]).to_pandas()
	pav = dt.fread(glob.glob(f"{d}/PAV/{model_dir}/shap_interaction_scores_pav_{env}_top_*_lit_genes_summed.txt")[0]).to_pandas()
	cnv = dt.fread(glob.glob(f"{d}/CNV/{model_dir}/shap_interaction_scores_cnv_{env}_top_*_lit_genes_summed.txt")[0]).to_pandas()
	
	# Is there experimental evidence for these SHAP interactions?
	try:
		snp["Known_biogrid"] = snp.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in biogrid_gp else 0, axis=1)
		snp["Known_costanzo_2021"] = snp.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in costanzo_gp else 0, axis=1)
		print('snp biogrid', snp["Known_biogrid"].value_counts())
		print('snp costanzo', snp["Known_costanzo_2021"].value_counts())
		num_int[("biogrid", "snp", env)] = snp["Known_biogrid"].value_counts()
		num_int[("costanzo", "snp", env)] = snp["Known_costanzo_2021"].value_counts()
	except:
		pass
	pav["Known_biogrid"] = pav.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in biogrid_gp else 0, axis=1)
	pav["Known_costanzo_2021"] = pav.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in costanzo_gp else 0, axis=1)
	print('pav biogrid', pav["Known_biogrid"].value_counts())
	print('pav costanzo', pav["Known_costanzo_2021"].value_counts())
	num_int[("biogrid", "pav", env)] = pav["Known_biogrid"].value_counts()
	num_int[("costanzo", "pav", env)] = pav["Known_costanzo_2021"].value_counts()
	cnv["Known_biogrid"] = cnv.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in biogrid_gp else 0, axis=1)
	cnv["Known_costanzo_2021"] = cnv.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in costanzo_gp else 0, axis=1)
	print('cnv biogrid', cnv["Known_biogrid"].value_counts())
	print('cnv costanzo', cnv["Known_costanzo_2021"].value_counts())
	num_int[("biogrid", "cnv", env)] = cnv["Known_biogrid"].value_counts()
	num_int[("costanzo", "cnv", env)] = cnv["Known_costanzo_2021"].value_counts()
	#
	# make sure it's saved as 0s and 1s not a boolean
	try:
		snp[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
			"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]] = \
			snp[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
			"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]].astype(int)
		snp.to_csv(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_plus_comb_lit_genes_summed_lit_genes_known.txt", index=False, sep="\t")
	except:
		pass
	pav[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]] = \
		pav[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]].astype(int)
	cnv[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]] = \
		cnv[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]].astype(int)
	pav.to_csv(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_plus_comb_lit_genes_summed_lit_genes_known.txt", index=False, sep="\t")
	cnv.to_csv(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_plus_comb_lit_genes_summed_lit_genes_known.txt", index=False, sep="\t")
	# try:
	# 	del snp, pav, cnv
	# except:
	# 	del pav, cnv

num_int = pd.DataFrame(num_int).T
num_int.fillna(0, inplace=True)

biogrid_sub = num_int.loc[num_int.index.get_level_values(0)=="biogrid",:]
sns.barplot(y = biogrid_sub[1].values,\
	x = biogrid_sub.index.get_level_values(2).values,\
	hue = biogrid_sub.index.get_level_values(1).values,
	palette="viridis")
plt.xticks(rotation=45, size=7)
plt.ylabel("Number of known GIs identified")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_5/Num_known_GIs_biogrid_identified_barplot.pdf")
plt.close()

costanzo_sub = num_int.loc[num_int.index.get_level_values(0)=="costanzo",:]
sns.barplot(y = costanzo_sub[1].values,\
	x = costanzo_sub.index.get_level_values(2).values,\
	hue = costanzo_sub.index.get_level_values(1).values,
	palette="viridis")
plt.xticks(rotation=45, size=7)
plt.ylabel("Number of known GIs identified")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_5/Num_known_GIs_costanzo_identified_barplot.pdf")
plt.close()