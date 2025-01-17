# First, got to the directory that contains the 
setwd("~/covid_play_R/Pathway_analysis_play/")


# Additional pathway analysis
# Adapted from https://bioconductor.org/packages/devel/bioc/vignettes/rWikiPathways/inst/doc/Pathway-Analysis.html

# Install and load necessary packages
if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways", update = FALSE)
  BiocManager::install("clusterProfiler", update = FALSE)
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
library(rWikiPathways)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Source script from Kevin's githib repository
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_metadata.r")

# Get the gene ids from Kevin's study that were found to be significantly deferentially expressed
# among the three considered studies
# First import the results table
my_data <- import_metadata(group_table ="~/covid_play_R/manuscipt/Supplemental_Table_1.csv")
# Then extract just the gene ids (Ensemble gene ids)
my_genes <- as.character(rownames(my_data))

# Convert Ensemble ids to Entrez
my_genes.entrez <- clusterProfiler::bitr(my_genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
cat("\n\nWhich column contains my new Entrez IDs?\n")
head(up.genes.entrez)

# Get the gene Ontology
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

# with ggplot
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

barplot(pathway_analysis, showCategory = 20)
# with ggplot
ggplot(pathway_analysis[1:20], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))

dotplot(pathway_analysis, showCategory = 20)
