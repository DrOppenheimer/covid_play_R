setwd("~/covid_play_R/Pathway_analysis_play/")


# Adapted from https://bioconductor.org/packages/devel/bioc/vignettes/rWikiPathways/inst/doc/Pathway-Analysis.html

# Install main package
if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways", update = FALSE)
}
library(rWikiPathways)

# Install all dependent packages
load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
install.packages(load.libs)

#options(install.packages.check.source = "no")
#options(install.packages.compile.from.source = "never")
#if (!require("pacman")) install.packages("pacman"); library(pacman)
#p_load(load.libs, update = TRUE, character.only = TRUE)
#status <- sapply(load.libs,require,character.only = TRUE)
#if(all(status)){
#  print("SUCCESS: You have successfully installed and loaded all required libraries.")
#} else{
#  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
#  status
#}


# Load the data
lung.expr <- read.csv(system.file("extdata","data-lung-cancer.csv", package="rWikiPathways"),stringsAsFactors = FALSE)
nrow(lung.expr)
head(lung.expr)
# Get up and down regulated lists of genes
up.genes <- lung.expr[lung.expr$log2FC > 1 & lung.expr$adj.P.Value < 0.05, 1] 
dn.genes <- lung.expr[lung.expr$log2FC < -1 & lung.expr$adj.P.Value < 0.05, 1]
bkgd.genes <- lung.expr[,1]

# Convert Ensemble IDs to Entrez
BiocManager::install("clusterProfiler", update = FALSE)
library(clusterProfiler)
BiocManager::install("org.Hs.eg.db", update = FALSE)
library(org.Hs.eg.db)
up.genes.entrez <- clusterProfiler::bitr(up.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez <- bitr(dn.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
cat("\n\nWhich column contains my new Entrez IDs?\n")
head(up.genes.entrez)

# Here’s the complete list of identifiers that this particular tool can convert across. You have to spell these precisely and in all caps for the bitr function to work:
keytypes(org.Hs.eg.db)


# Gene Ontology
egobp <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

head(egobp,10)

# plot the results
barplot(egobp, showCategory = 20)
dotplot(egobp, showCategory = 20)
goplot(egobp)

# with ggplot
library(ggplot2)
ggplot(egobp[1:20], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))

# Skip enrichment map
# Skip Over-representation Analysis (ORA)

# Try Gene Set Enrichment Analysis (GSEA)

lung.expr$fcsign <- sign(lung.expr$log2FC)
lung.expr$logfdr <- -log10(lung.expr$P.Value)
lung.expr$sig <- lung.expr$logfdr/lung.expr$fcsign
sig.lung.expr.entrez<-merge(lung.expr, bkgd.genes.entrez, by.x = "GeneID", by.y = "ENSEMBL")
gsea.sig.lung.expr <- sig.lung.expr.entrez[,8]
names(gsea.sig.lung.expr) <- as.character(sig.lung.expr.entrez[,9])
gsea.sig.lung.expr <- sort(gsea.sig.lung.expr,decreasing = TRUE)

gwp.sig.lung.expr <- clusterProfiler::gseWP(
  gsea.sig.lung.expr,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff
  organism = "Homo sapiens"
)

gwp.sig.lung.expr.df = as.data.frame(gwp.sig.lung.expr)
gwp.sig.lung.expr.df[which(gwp.sig.lung.expr.df$NES > 1),] #pathways enriched for upregulated lung cancer genes
gwp.sig.lung.expr.df[which(gwp.sig.lung.expr.df$NES < -1),] #pathways enriched 




