# By Scott
# NOTE I am not as proficient in R as I am in Python, however I tried my best to provide
# ample notes in this script.
# Obviously not all genes can be mapped onto the GO, we never did anything with the genes that were lost,
# FUTURE maybe investigate more...
# NOTE we did not take GO hierarchy into account because since this dataset is already downstream
# of a lot of filtering, we did not want to penalize the model further. 

# NOTE each installation command must be done separately
#install.packages("BiocManager")
#BiocManager::install(lib='/home/scott/R/x86_64-pc-linux-gnu-library/4.1')
#BiocManager::install("topGO", lib='/home/scott/R/x86_64-pc-linux-gnu-library/4.1', force=TRUE)
#BiocManager::install("ALL")
#BiocManager::install("Rgraphviz", lib='/home/scott/R/x86_64-pc-linux-gnu-library/4.1')

suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(Rgraphviz))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tools))


############################################
# Check to see you have the args
args= commandArgs(trailingOnly=TRUE)

# NOTE get input files from input arguments
# MAGIC represents folder of modules in Arabidopsis gene format
files_w_tsv = Sys.glob(file.path(args[1], "*.tsv"))

# Define your mapping (gene universe), MAGIC input arg 2
geneID2GO = readMappings(file=args[2], sep='\t') # NB global
master_genes = names(geneID2GO)  # NB global

# output directory
output_dir = args[3]
doc_dir = args[4]
############################################



function_run_topgo = function(master_genes, geneID2GO, my_interesting_genes, module_name, ontology_group, output_dir){
	geneList = factor(as.integer(master_genes %in% my_interesting_genes))
	names(geneList) = master_genes # NB give it names

	GOdata = new('topGOdata',
		     ontology = ontology_group,
		     allGenes = geneList,
		     annot = annFUN.gene2GO,
		     gene2GO = geneID2GO)

	# NOTE get the number of genes in this interesting list
	# NOTE how many genes am I losing when I perform this step?
	# Do some genes not have the 'MF' ontology or something?
	print(GOdata)  # NB shows the number of genes lost

	sig_genes = sigGenes(GOdata)  # the significant genes
	num_sig_genes = numSigGenes(GOdata) # NB the no. of signifcant genes,
	# not necessarily the nodes in the graph
	num_nodes = length(usedGO(object=GOdata))  # NB all nodes in the graph,
	# not necessarily the significant nodes
	# FUTURE would this be useful information to store?

	#resultFisher_weight = runTest(GOdata, algorithm='weight01', statistic='fisher')
	resultFisher_classic = runTest(GOdata, algorithm='classic', statistic='fisher')


	# NB from Nolan:
		# algorithm="classic" WILL NOT take GO hierarchy into account
	    	# ^ The limitation of using "classic" is that all genes
       		# annotated to a GO terms will be automatically annotated
       		# to its parents as well, therefore a GO term might look
       		# enriched just because its children are enriched
	    	# algorithm="weight01" WILL take GO hierarchy into account
		# These p-values have not been corrected for multiple testing

	# NB from Pat:
       		# After discussion with Pat, since I already have a cutoff applied
       		# to my modules, and the GO hieracrhy is useful information, I don't want
       		# to penalize the stats further, so I am going with classic and won't do
       		# an FDR.

	# list the top significant GO terms
    	allRes = GenTable(GOdata, classicFisher=resultFisher_classic, 
			  orderBy='classicFisher',
			  ranksOf='classicFisher',
			  topNodes=num_nodes,
			  numChar=1000) # NB avoid truncation of string
			      
			      
	# MAGIC p-val threshold
       	pval_threshold = 0.05
	# NB apply the p-val threshold
	allRes_significant = allRes[which(allRes$classicFisher < pval_threshold),]


	# NB add a string to a new column that says 'yes' if we have
       	# more significant than expected
	allRes_significant$Overrepresented = ifelse(allRes_significant$Significant > allRes_significant$Expected, 'Y', 'N')

	# NB A GO term is synonymoous with node
	# NB node size defaults to 1, no pruning is performed
	# This is the minimum number of genes annotated to a go

	# MAGIC filename
	write.table(allRes_significant, file=file.path(output_dir, paste(module_name, '_', ontology_group, "_",pval_threshold,"_pval_sig_GO_Summary_Table.tsv", sep='')), sep='\t', quote=FALSE, row.names=FALSE)

	# Get the graphic of connected nodes colored by significance
	outpdf = file.path(output_dir, paste(module_name,"_",ontology_group,"_",pval_threshold,"_pval_sig_GO_tree.pdf", sep=''))
    	pdf(file=outpdf)
	# NB the nodes shown (rectanlges) are significant, and the
        # ovals are the connecting information? Unsure. 
	# NB arbitrary number of sig nodes shown for display/legibility reasons
    	showSigOfNodes(GOdata, score(resultFisher_classic), firstSigNodes = 15,
		       useInfo ='all')	 
    	dev.off()
}

for(input_module in files_w_tsv){	
	module_name = file_path_sans_ext(basename(input_module))

	# Import a list of predefined interesting genes
	my_interesting_genes = read.table(input_module, sep = '\t')

	# MAGIC convert to character vector and use magic column-name
	my_interesting_genes = as.character(my_interesting_genes$V1)

	# Run the code and do each ontology group	
	function_run_topgo(master_genes, geneID2GO, my_interesting_genes, module_name, 'MF', output_dir)
	function_run_topgo(master_genes, geneID2GO, my_interesting_genes, module_name, 'BP', output_dir)
	function_run_topgo(master_genes, geneID2GO, my_interesting_genes, module_name, 'CC', output_dir)
}

# NB this outputs the session information for easy package management.
# Perhaps better than outputting a conda env
# MAGIC get nice doc dir
sink(file=paste(doc_dir, "TopGO_sessionInfo.txt", sep='/'))
sessionInfo()
sink()