# The goal of this workflow is to find the genes in common among the analysis of all three 
# data sets: GSE198449, GSE215865, and GSE212041. The assumption is that genes found in common
# among these three analyses, each of which searched for genes that exhibited the most significant
# differences between/among groups that were unambiguously labeled as testing positive for Covid or
# not, represent a robust list of candidate genes - suitable for further investigation. 

# Load all required packages ----------------------------------------------
# Note, Dependency installation on Ubuntu : sudo apt-get install libzmq3-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev build-essential libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev

library(gprofiler2)

# Source all additional functions from Kevin's github repository ----------------------------------------------

source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_data.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_metadata.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/export_data.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/heatmap_dendrogram.r")

# Create a directory for working (if it doesn't already exist) and move to it --------
# Specify the directory path
dir_path <- "/Users/kosso/OneDrive/Documents/combine_covid/"
# Check if the directory exists
if (dir.exists(dir_path)) {
  stop("Error: The directory already exists.")
} else {
  # Create the directory if it does not exist
  dir.create(dir_path)
  print("Directory created successfully.")
}
# move to that directory 
setwd(dir_path)

# Import data calculated from the three studies ----------------------------------------------
# Note: The workflow below assumes that you have completed the analysis workflows for each of the three
# datasets. It also assumes that you saved your work 

# GSE198449
# Import the statistically subselected data with stats removed
GSE198449_stat <- import_data("~/GSE198449/GSE198449_stat_results_subselected_and_cleaned.txt")
# The rownames of this dataset contain the ENSG stable gene IDs along with a version tag. The tag is removed
# to allow for comparison of just the stable gene ID
rownames(GSE198449_stat) <- gsub(pattern="\\..*", replacement = "", x=rownames(GSE198449_stat))
ncol(GSE198449_stat) # i.e., 506 samples

# GSE215865
# Import the statistically subselected data with stats removed
GSE215865_stat <- import_data("~/GSE215865/GSE215865_stat_results_subselected_and_cleaned.txt")
# The rownames of this dataset contain the ENSG stable gene IDs along with a version tag. The tag is removed
# to allow for comparison of just the stable gene ID
rownames(GSE215865_stat) <- gsub(pattern="\\..*", replacement = "", x=rownames(GSE215865_stat))
ncol(GSE215865_stat) # i.e., 1391 samples

# GSE212041
# Import the statistically subselected data with stats removed
GSE212041_stat <- import_data("~/GSE212041/GSE212041_stat_results_subselected_and_cleaned.txt")
# The rownames of this dataset contain the ENSG stable gene ID WITHOUT a version tag. No modification
# necessary.
ncol(GSE212041_stat) # i.e., 781 samples

# Use an intersection to find the ENSG stable gene IDs that are in common 
# between the first two analyses
first_pair <- intersect( rownames(GSE198449_stat),rownames(GSE215865_stat) )
length(first_pair) # 134
# Now among all three analyses
all_three <- intersect( first_pair,rownames(GSE212041_stat) )
length(all_three) # 61

# ANNOTATE RESULTS --------------------------------------------------------
# In the first analysis (of GSE198449) the results were annotated using data downloaded
# from https://useast.ensembl.org/biomart/martview/ We use that annotation file again here
# You may have to change the path of the file below.
# Load the annotations downloaded in the analysis of GSE198449
annotations <- read.table(file="~/GSE198449/mart_export.txt",row.names=NULL,header=TRUE,sep="\t", #changed row.names from 1 to NULL
                          colClasses = "character", check.names=FALSE,
                          comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE)

# Prepare a matrix from the gene IDs found in common that will hold the IDs along with annotation information
# For the sake of simplicity we just include "Gene description" and "Gene name"
my_results <- matrix(ncol=ncol(annotations), nrow=length(all_three))
colnames(my_results) <- colnames(annotations)
# Use a simple loop to build a table that contains the "Gene stable ID", "Gene description" and "Gene name"
# for each gene in the list of IDs that were found in common.
for (i in 1:length(all_three)){
  # Find the annotations for our genes
  my_annotation <- annotations[annotations$`Gene stable ID` == all_three[i],]
  # Add them to the matrix for the results
  my_results[i,"Gene stable ID"] <- my_annotation[1,"Gene stable ID"]
  my_results[i,"Gene description"] <- my_annotation[1,"Gene description"]
  my_results[i,"Gene name"] <- my_annotation[1,"Gene name"]
}
# Clean up the results by removing NA rows - introduced by the initial construction of the matrix, 
# and an ID that was not annotated -- we'll look into the latter below
# Identify rows with "NA" in the selected column.
na_rows <- which(is.na(my_results[, "Gene stable ID"]))
# Remove identified "NA" rows
my_results <- my_results[-na_rows, ]
# Check to make sure that we still have 61 genes
dim(my_results) # 60 3 -- one wasn't annotated. Figure out which one it was.
setdiff( all_three, my_results[,"Gene stable ID"] ) # ENSG00000234290, not clear why it was not annotated
# This test indicate that there is no annotation for this ID in annotations we downloaded from Ensembl
annotations[annotations$`Gene stable ID` == "ENSG00000234290",] # <0 rows> (or 0-length row.names)
# A quick search of the internet reveals that this ID is annotated elsewhere
# The solution is here:
# https://useast.ensembl.org/Homo_sapiens/Gene/Idhistory?g=ENSG00000234290
# Status Retired 
# Add the unannotated ENSG back to the results 
my_results <- rbind(my_results, c("ENSG00000234290", NA, NA))
# Write the mostly annotated list of genes to file
write.table(file = "In_common.txt", x = my_results, sep="\t",  quote=FALSE, col.names=TRUE, row.names=FALSE)

# Heatmap-Dendrograms of the selected genes across all three datasets --------

# The data from each of the datasets is subselected to just the genes found in common
# Those data are saved and used to create heatmap-dendrograms
# Start with GSE198449
GSE198449_stat_all_three <- GSE198449_stat[all_three,]
# make sure that data colnames and metadata rownames are the same  
GSE198449_sorted_samples <- sort(colnames(GSE198449_stat_all_three))
GSE198449_stat_all_three <- GSE198449_stat_all_three[,GSE198449_sorted_samples]
GSE198449_metadata <- import_metadata("~/GSE198449/GSE198449.metadata.txt")
GSE198449_metadata <- GSE198449_metadata[GSE198449_sorted_samples,]
# Export sorted data
export_data(GSE198449_stat_all_three, "GSE198449_stat_all_three.txt")
# Export sorted metadata
export_data(GSE198449_metadata, "GSE198449.metadata.txt")
# Create hetamap dendrogram
heatmap_dendrogram(file_in = "GSE198449_stat_all_three.txt",
                   metadata_table = "GSE198449.metadata.txt",
                   metadata_column="PCR.test.for.SARS.Cov.2"
)
shell("start GSE198449_stat_all_three.txt.HD.png")

# Now GSE215865
GSE215865_stat_all_three <- GSE215865_stat[all_three,]
# make sure that data colnames and metadata rownames are the same  
GSE215865_sorted_samples <- sort(colnames(GSE215865_stat_all_three))
GSE215865_stat_all_three <- GSE215865_stat_all_three[,GSE215865_sorted_samples]
GSE215865_metadata <- import_metadata("~/GSE215865/GSE215865.metadata.txt")
GSE215865_metadata <- GSE215865_metadata[GSE215865_sorted_samples,]
# export sorted data
export_data(GSE215865_stat_all_three, "GSE215865_stat_all_three.txt")
# export sorted metadata
export_data(GSE215865_metadata, "GSE215865.metadata.txt")
# create hetamap dendrogram
heatmap_dendrogram(file_in = "GSE215865_stat_all_three.txt",
                   metadata_table = "GSE215865.metadata.txt",
                   metadata_column="covid.19_positive"
)
shell("start GSE215865_stat_all_three.txt.HD.png")

# Now GSE212041
GSE212041_stat_all_three <- GSE212041_stat[all_three,]
# make sure that data colnames and metadata rownames are the same  
GSE212041_sorted_samples <- sort(colnames(GSE212041_stat_all_three))
GSE212041_stat_all_three <- GSE212041_stat_all_three[,GSE212041_sorted_samples]
GSE212041_metadata <- import_metadata("~/GSE212041/GSE212041.metadata.txt")
GSE212041_metadata <- GSE212041_metadata[GSE212041_sorted_samples,]
# export sorted data
export_data(GSE212041_stat_all_three, "GSE212041_stat_all_three.txt")
# export sorted metadata
export_data(GSE212041_metadata, "GSE212041.metadata.txt")
# create hetamap dendrogram
heatmap_dendrogram(file_in = "GSE212041_stat_all_three.txt",
                   metadata_table = "GSE212041.metadata.txt",
                   metadata_column="patient_category"
)
shell ("start GSE212041_stat_all_three.txt.HD.png")

# PRELIMINARY PATHWAY ANALYSIS --------------------------------------------

# Use the list of gene ids from above
all_three

# Use the gprofiler2 package to perform preliminary pathway analysis
gostres <- gost(query = all_three, organism = 'hsapiens', significant = TRUE)
gostplot(gostres, capped = FALSE, interactive = TRUE)

# ADDITIONAL PATHWAY ANALYSIS --------------------------------------------

# Adapted from https://bioconductor.org/packages/devel/bioc/vignettes/rWikiPathways/inst/doc/Pathway-Analysis.html

# Install and load necessary packages
if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways", update = FALSE)
  BiocManager::install("clusterProfiler", update = FALSE)
  BiocManager::install("org.Hs.eg.db", update = FALSE)
  install.packages("ggplot2")
}
library(rWikiPathways)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Source a script from Kevin's githib repository
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_metadata.r")

# Get the gene ids from Kevin's study that were found to be significantly deferentially expressed
# among the three considered studies
# Use the list of genes from above
my_genes <- as.character(all_three)

# Convert Ensemble ids to Entrez
my_genes.entrez <- clusterProfiler::bitr(my_genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
cat("\n\nWhich column contains my new Entrez IDs?\n")
head(my_genes.entrez) # The second column

# Get the gene Ontology # This will take a minute
my_ont <- clusterProfiler::enrichGO(
  gene     = my_genes.entrez[[2]],
  #universe = NULL,
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = TRUE)

# plot the results
barplot(my_ont, showCategory = 20)
dotplot(my_ont, showCategory = 20)
goplot(my_ont)

# barplot with ggplot
ggplot(my_ont[1:20], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))

# Pathway analysis
pathway_analysis <- clusterProfiler::enrichWP(
  my_genes.entrez[[2]],
  #universe = bkgd.genes.entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff; relaxed for demo purposes
)

head(pathway_analysis)

# For some reason, enricher doesn’t automatically add gene symbols to the result object, 
# but there is a handy function in DOSE that does…
pathway_analysis <- DOSE::setReadable(pathway_analysis, org.Hs.eg.db, keyType = "ENTREZID")
head(pathway_analysis)

# Plot the results
barplot(pathway_analysis, showCategory = 20)
# with ggplot
ggplot(pathway_analysis[1:20], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))

dotplot(pathway_analysis, showCategory = 20)

# That's it for now. I hope you have enjoyed and/or learned something new from this
# basement project. Cheers, Kevin K 