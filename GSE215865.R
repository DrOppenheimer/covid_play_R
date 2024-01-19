# Data were found by a search in:
# https://www.ncbi.nlm.nih.gov/geo/browse/
# The keyword "covid" was used and studies were sorted by "Series Type"
# to find "Expression profiling by high throughput sequencing", and by
# "Samples" to find the studies with the largest number of samples.
# As of 1/15/2024 three studies were found with more than 500 samples 
# The analysis below is for one of those three datasets.
# All data in this analysis come from:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215865
# Note that the metadata for this dataset could not be downloaded programmatically.
# The procedure for downloading the metadata is outlined below.

# load all required packages ----------------------------------------------

library(readxl)
library(tidyverse)
library(gprofiler2)
library(plotly)

# source all additional functions from Kevin's github repository ----------------------------------------------

source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_data.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_metadata.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/preprocessing_tool.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/export_data.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/calculate_pco.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/render_calculated_pcoa.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/sigtest.R")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/remove_last_n_columns.R")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/heatmap_dendrogram.r")

# Create a directory for working (if it doesn't already exist) and move to it --------
# Specify the directory path
dir_path <- "~/Downloads/GSE215865/"
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

# Download data and metadata ----------------------------------------------

# First the data ...
# Specify the FTP URL of the file
data_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE215nnn/GSE215865/suppl/GSE215865_rnaseq_raw_count_matrix.csv.gz"
# Specify the local file path to save the downloaded content
local_file_path_data <- "GSE215865_rnaseq_raw_count_matrix.csv.gz"
# Download the file
download.file(data_url, destfile = local_file_path_data, mode = "wb")
# now we have to unzip the gz file, easiest to just do this with a system call
system( paste("gunzip", local_file_path_data) )

# Now the metadata ...
# Go to this page:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215865
# Click on the "SRA Run Selector" link close to the very bottom of the page
# On the page you are taken to in the "Download" column of the "Select" section click on "Metadata" 
# This will download "SraRunTable.txt", move it to this directory and you can proceed below

# Import and explore the data and the metadata -----------------------------------------

# import raw data and metadata
GSE215865_metadata <- read_csv("SraRunTable.txt")
GSE215865_data <- read_csv("GSE215865_rnaseq_raw_count_matrix.csv")

# see how many samples there are at the start
ncol(GSE215865_data) - 1 # 1392, first column has "Ensembl_Gene_ID" 

# look at sample ids in the data and metadata
colnames(GSE215865_data)[1:3] #  They look like this in the data: "Subj_16ae238fT0_Plate_2"
length(colnames(GSE215865_data)) # a column for the Ensembl Gene IDs followed by 1392 samples

# reconstruct sample ID in metadata that matches those found in colnames of the data
GSE215865_metadata$blood_sample_id[1] # "Subj_5d362febT1"
GSE215865_metadata$library_prep_plate[1] # "Plate_8"
paste( GSE215865_metadata$blood_sample_id[1], GSE215865_metadata$library_prep_plate[1], sep="_" ) # "Subj_5d362febT1_Plate_8"
# check to make sure this reconstructed sample name is actually in the data
"Subj_5d362febT1_Plate_8" %in% colnames(GSE215865_data) # TRUE
# Now add a column with the reconstructed sample ID to the metadata
length(colnames(GSE215865_metadata))
GSE215865_metadata <- GSE215865_metadata |>
  mutate( GSE215865_sample_id = paste( blood_sample_id, library_prep_plate, sep = "_") )
length(colnames(GSE215865_metadata))   # we added a 40th column to the metadata

# Try to get IDs to match between the data and metadata, cull each accordingly ------------------------
  
# Find samples that are in common between data and metadata
common_samples <- intersect( GSE215865_metadata$GSE215865_sample_id, colnames(GSE215865_data) )
common_samples <- sort(common_samples)
length(common_samples) # 1391
common_samples[1:3] # take a look at the first three sample IDs

# subselect data and metadata for the 1391 samples for which there are data and metadata
# and then sort/order
# First the metadata 
GSE215865_selected_metadata <- GSE215865_metadata |>
  filter( GSE215865_sample_id %in% common_samples ) |>
  arrange( GSE215865_sample_id )
dim(GSE215865_selected_metadata) # 1391 40
# Create rownames from the "GSE215865_sample_id" column
GSE215865_selected_metadata <- column_to_rownames(GSE215865_selected_metadata, var="GSE215865_sample_id")
# (rownames_to_column() and column_to_rownames() # Two useful Tibble functions for dealing with rownames
rownames(GSE215865_selected_metadata) # double check to make sure that the rownames are there

# Now the data
# subselect with indices
GSE215865_selected_data <- GSE215865_data[,c("Ensembl_Gene_ID",common_samples)]
dim(GSE215865_selected_data)
# Create rownames from the "Ensembl_Gene_ID" column
GSE215865_selected_data <- column_to_rownames(GSE215865_selected_data, var="Ensembl_Gene_ID")
dim(GSE215865_selected_data)

# SAVE THE DATA AND THE METADATA TO FILE (First attempt) ----------------------------------
# Now save the metadata to file
write.table(file = "GSE215865.metadata.txt", x = GSE215865_selected_metadata, sep="\t",  quote=FALSE, col.names=NA)  
# Write the data to file
write.table(file = "GSE215865.data.txt", x = GSE215865_selected_data, sep="\t", quote=FALSE, col.names=NA) 

# Check that subselecting "worked" ----------------------------------------

# import data and metadata from files
GSE215865_data <- import_data("GSE215865.data.txt")
GSE215865_metadata <- import_metadata("GSE215865.metadata.txt")
# check that samples are in common
length( intersect( colnames(GSE215865_data), rownames(GSE215865_metadata) ) ) # 1391 samples are in common
# check that samples are identically ordered
identical( colnames(GSE215865_data), rownames(GSE215865_metadata) ) #FALSE - but ordering did not work, try again

# try to reorder
# create an ordered index of the samples
ordered_sample_names <- sort( colnames(GSE215865_data) )
# order the data columns by this index
GSE215865_data <- GSE215865_data[,ordered_sample_names]
# order the metadata rows by this index
GSE215865_metadata <- GSE215865_metadata[ordered_sample_names,]
# check for identity once again
identical( colnames(GSE215865_data), rownames(GSE215865_metadata) ) # TRUE

# SAVE THE DATA AND THE METADATA TO FILE (Second attempt) ----------------------------------
write.table(file = "GSE215865.data.txt", x = GSE215865_data, sep="\t", quote=FALSE, col.names=NA) 
write.table(file = "GSE215865.metadata.txt", x = GSE215865_metadata, sep="\t",  quote=FALSE, col.names=NA)  

# OPTIONAL Clean House ----------------------------------------------------

# At this point you may want to clean house the dreaded "rm(list = ls())" command
# Note the if you do, you will have to run the bit of code from above that sources
# functions from github. Sorry they're not in a package as of yet.

# FUN STUFF ---------------------------------------------------------------

# The analysis really starts at the this point, everything up until now has been
# data wrangling to make the analyses below possible

# import the data and metadata  -------------------------------------------
GSE215865_data <- import_data("GSE215865.data.txt")
GSE215865_metadata <- import_metadata("GSE215865.metadata.txt")

# look at distribution of raw data ----------------------------------------

# explore distribution of data
all_GSE215865_data <- as.vector(GSE215865_data)
hist(all_GSE215865_data, breaks =100)
hist(log10(all_GSE215865_data), breaks =100) # not really log normal
# with ggplot
all_GSE215865_data <- data.frame(Values = as.vector(GSE215865_data))
ggplot(data=all_GSE215865_data,
       mapping=aes(x=Values))+
  geom_histogram()


# look at the distribution of summary values in the raw data
GSE215865_raw_mean <- apply(GSE215865_data, 2, mean) 
GSE215865_raw_median <- apply(GSE215865_data, 2, median) 
GSE215865_raw_min <- apply(GSE215865_data, 2, min) 
GSE215865_raw_max <- apply(GSE215865_data, 2, max) 
GSE215865_raw_sd <- apply(GSE215865_data, 2, sd)

plot(GSE215865_raw_mean, pch=20) # many outliers, roughly anything above 800 or below 200
plot(GSE215865_raw_median, pch=20) # surprising amount of structure in this data, appears to be tiered
plot(GSE215865_raw_min, pch=20) # no surprises
plot(GSE215865_raw_max, pch=20) # one extreme outlier, above 2.5e+07
plot(GSE215865_raw_sd, pch=20) # one extreme outlier, above 120,000

# make these data tidy for ggplot
GSE215865_raw_descriptive_stats <- as.matrix(names(GSE215865_raw_mean))
GSE215865_raw_descriptive_stats <- cbind(GSE215865_raw_descriptive_stats, as.matrix(GSE215865_raw_mean))
GSE215865_raw_descriptive_stats <- cbind(GSE215865_raw_descriptive_stats, as.matrix(GSE215865_raw_median))
GSE215865_raw_descriptive_stats <- cbind(GSE215865_raw_descriptive_stats, as.matrix(GSE215865_raw_min))
GSE215865_raw_descriptive_stats <- cbind(GSE215865_raw_descriptive_stats, as.matrix(GSE215865_raw_max))
GSE215865_raw_descriptive_stats <- cbind(GSE215865_raw_descriptive_stats, as.matrix(GSE215865_raw_sd))
colnames(GSE215865_raw_descriptive_stats) <- c("GSE215865_sample_ids","GSE215865_mean", "GSE215865_median", "GSE215865_min", "GSE215865_max", "GSE215865_sd")
rownames(GSE215865_raw_descriptive_stats) <- names(GSE215865_raw_mean)
GSE215865_raw_descriptive_stats <- as.data.frame(GSE215865_raw_descriptive_stats)
GSE215865_raw_descriptive_stats$GSE215865_mean <- as.numeric(GSE215865_raw_descriptive_stats$GSE215865_mean)
GSE215865_raw_descriptive_stats$GSE215865_median <- as.numeric(GSE215865_raw_descriptive_stats$GSE215865_median)
GSE215865_raw_descriptive_stats$GSE215865_min <- as.numeric(GSE215865_raw_descriptive_stats$GSE215865_min)
GSE215865_raw_descriptive_stats$GSE215865_max <- as.numeric(GSE215865_raw_descriptive_stats$GSE215865_max)
GSE215865_raw_descriptive_stats$GSE215865_sd <- as.numeric(GSE215865_raw_descriptive_stats$GSE215865_sd)

ggplot(
  data=GSE215865_raw_descriptive_stats,
  mapping=aes(x=1:nrow(GSE215865_raw_descriptive_stats))) +
  scale_y_log10() +
  geom_line(aes(y=GSE215865_mean, colour="GSE215865_mean")) +
  geom_line(aes(y=GSE215865_median, colour="GSE215865_median")) +
  geom_line(aes(y=GSE215865_min, colour="GSE215865_min")) +
  geom_line(aes(y=GSE215865_max, colour="GSE215865_max")) +
  geom_line(aes(y=GSE215865_sd, colour="GSE215865_sd")) +
  scale_color_manual(name = "Stat Values", values = c("GSE215865_mean" = "black", "GSE215865_median" = "red", "GSE215865_min"="blue", "GSE215865_max"="green", "GSE215865_sd"="purple"))
# The outliers don't look so bad in log space, for now will keep all samples


# Preprocess and look at distribution of normed data -----------------------------------------------------

# Preprocess the data -- attempt to normalize it
preprocessing_tool("GSE215865.data.txt", pseudo_count=0.1)

# look at distribution of normed data
GSE215865_data_n <- import_data("GSE215865.data.txt.quantile.PREPROCESSED.txt")
all_GSE215865_data_n <- as.vector(GSE215865_data_n)
hist(all_GSE215865_data_n, breaks =100)
hist(log10(all_GSE215865_data_n), breaks =100) # data have a peculiar distribution, definitely non-normal

# look at the distribution of summary values in the normed data
GSE215865_norm_mean <- apply(GSE215865_data_n, 2, mean)
GSE215865_norm_median <- apply(GSE215865_data_n, 2, median) 
GSE215865_norm_min <- apply(GSE215865_data_n, 2, min) 
GSE215865_norm_max <- apply(GSE215865_data_n, 2, max) 
GSE215865_norm_sd <- apply(GSE215865_data_n, 2, sd)

plot(GSE215865_norm_mean) # anything below 0.10 is pretty anomalous 
plot(GSE215865_norm_median) # data still has an odd tiered nature
plot(GSE215865_norm_min) # tiered structure is still apparent
plot(GSE215865_norm_max) # a group of somewhat aberrant points around 0.98, two extreme at around 0.1
plot(GSE215865_norm_sd) # anything below 0.15 seems to be fairly anomalous
# Not sure what to make of the tiered nature of the data; not understanding it I'm 
# hesitatnt to toss any of the samples. We may be able to see how they correlate with
# the metadata in the PCoAs below. My guess is non-biological structure (bias, batach effects)
# that can be identified and removed later

# PCoA --------------------------------------------------------------------

# calculate PCoA on the preprocessed data
calculate_pco(file_in="GSE215865.data.txt.quantile.PREPROCESSED.txt")

# Iterate through the metadata creating a static 3d PCoA for each metadata column
plot_static_colored_3d_pcoas(
  pcoa_filename = "GSE215865.data.txt.quantile.PREPROCESSED.txt.euclidean.PCoA",
  selected_eigen_vectors = c(1,2,3),
  metadata_filename = "GSE215865.metadata.txt",
  debug = TRUE
)
system("open *.PCoA.*.png")
# Some of the outlier samples give the PCoA an unusual appearance.
# Some correlation between the distribution of the points and "covid-19_positive"
# Also some correlation to "patient_classification_at_first_sample"
# and "sampling_time_point_label"

# Render an interactive 3d PCoA plot from the PCoA and the corresponding metadata
source("~/Documents/GitHub/workflow_play/render_calculated_pcoa.r")
plot_interactive_colored_3d_pcoa(
  pcoa_data_file = "GSE215865.data.txt.quantile.PREPROCESSED.txt.euclidean.PCoA",
  selected_eigen_vectors = c(1,2,3),
  pcoa_metadata_file = "GSE215865.metadata.txt",
  metadata_column = "covid-19_positive" 
)
# The tails in the PCoA may make it difficult to see any other obvious relationships
# Can try to identify the samples that have extreme PCoA coordinates and remove them

# try to find samples that have extreme PCoA coordinates
GSE215865_pcoa <- load_pcoa_data("GSE215865.data.txt.quantile.PREPROCESSED.txt.euclidean.PCoA")

# Take a look at the distribution of values in the first three vectors of the PCoA
boxplot( GSE215865_pcoa$eigen_vectors[,1:3] )

hist( GSE215865_pcoa$eigen_vectors[,1:3] )
# looks like values less than -10 are the outliers, so will remove them and look again

# create a data.frame that contains the vectors
GSE215865_vectors <- as.data.frame(GSE215865_pcoa$eigen_vectors)
# move the rownames to a column names "GSE215865_rownames" to keep the data tidy
GSE215865_vectors <- rownames_to_column(GSE215865_vectors,var="GSE215865_rownames")
# remove the quotes from "GSE215865_rownames" values
GSE215865_vectors$GSE215865_rownames <- str_replace_all(GSE215865_vectors$GSE215865_rownames, '"', '')
dim(GSE215865_vectors) # 1391 1392
# look at the first few
GSE215865_vectors[1:3,1:3]
# Remove rows with values less than -10 in the first three vectors
GSE215865_vectors_filtered <- GSE215865_vectors |>
  filter(`"PCO1"` > -10) |>
  filter(`"PCO2"` > -10) |>
  filter(`"PCO3"` > -10)
dim(GSE215865_vectors_filtered) # 1351 1392, 40 rows/samples were removed
# look at the distribution of the remaining values, first with boxplots
boxplot( GSE215865_vectors_filtered[,2:4] ) # still a large number of outliers
# take a look at the dsitribution of values in the individual vectors
hist( as.numeric(GSE215865_vectors_filtered[,2]) ) # less than -2 appear to be outliers
hist( as.numeric(GSE215865_vectors_filtered[,3]) ) # less than -5 outliers interesting that this vector exhibits a normal distribution when the others do not
hist( as.numeric(GSE215865_vectors_filtered[,4]) ) # less than -4 outlier

# create an interactive plot of the filtered vector values
plot_ly(
  x = GSE215865_vectors_filtered[,2], 
  y = GSE215865_vectors_filtered[,3], 
  z = GSE215865_vectors_filtered[,4], 
  type = "scatter3d", 
  mode = "markers"
)

# Get rid of quotes in the rownames
GSE215865_vectors_filtered$GSE215865_rownames <- str_replace_all(GSE215865_vectors_filtered$GSE215865_rownames, '"', '')
# see how many samples were culled
length(intersect( GSE215865_vectors_filtered$GSE215865_rownames, GSE215865_vectors$GSE215865_rownames ))
length(setdiff( GSE215865_vectors$GSE215865_rownames, GSE215865_vectors_filtered$GSE215865_rownames ))
# lost just 40 samples

# The data may have some horseshoe artifact, try another PCoA with Bray-Curtis distance 
# call function without args to see available distance options
calculate_pco()
# run again using "bray-curtis" as the "dist_method"
calculate_pco(file_in="GSE215865.data.txt.quantile.PREPROCESSED.txt", dist_method = "bray-curtis")

# Render an interactive 3d PCoA plot from the PCoA calculated with Bray-Curtis distance
plot_interactive_colored_3d_pcoa(
  pcoa_data_file = "GSE215865.data.txt.quantile.PREPROCESSED.txt.bray-curtis.PCoA",
  selected_eigen_vectors = c(1,2,3),
  pcoa_metadata_file = "GSE215865.metadata.txt",
  metadata_column = "covid-19_positive" 
)
# Horseshoe is actually worse with bc distance, try another 
calculate_pco(file_in="GSE215865.data.txt.quantile.PREPROCESSED.txt", dist_method = "minkowski")
# Render an interactive 3d PCoA plot from the PCoA calculated with Bray-Curtis distance
plot_interactive_colored_3d_pcoa(
  pcoa_data_file = "GSE215865.data.txt.quantile.PREPROCESSED.txt.minkowski.PCoA",
  selected_eigen_vectors = c(1,2,3),
  pcoa_metadata_file = "GSE215865.metadata.txt",
  metadata_column = "covid-19_positive" 
)
# Looks similar to the bc PCoA, no horseshoe but there are one or more tails 
# For now, use all of the data, but clearly some of the samples should be culled to make
# results more interpretable 

# perform a stat test on the data -----------------------------------------

# We want to identify the genes that exhibit the most significant difference between
# "Negative" and "Positive" samples listed under the "covid-19_positive" metadata.
# We need an unpaired test for two groups that are non-normal, Mann-Witney was selected as
# the most appropriate test

# Try stats with all of the samples, including the outliers
sigtest(data_file="GSE215865.data.txt.quantile.PREPROCESSED.txt", 
        metadata_file="GSE215865.metadata.txt",  
        metadata_column="covid-19_positive", 
        stat_test="Mann-Whitney-unpaired-Wilcoxon",
        p_adjust_method = "BH"
)

# Load the stat test results and subselect data based on the stat results
GSE215865_stat_results <- import_data("GSE215865.data.txt.quantile.PREPROCESSED.Mann-Whitney-unpaired-Wilcoxon.covid-19_positive.STAT_RESULTS.txt")
# make the results "tidy" for ggplot
GSE215865_stat_results <- as_tibble(GSE215865_stat_results,rownames=NA)

# take a look at the stat results
ggplot(
  data=GSE215865_stat_results,
  mapping=aes(x=1:length(bonferroni_p))) +
  geom_line(aes(y=p, colour="p")) +
  geom_line(aes(y=bonferroni_p, colour="bonferroni_p")) +
  geom_line(aes(y=BH_p, colour="BH_p")) +
  scale_color_manual(name = "Stat Values", values = c("p" = "black", "bonferroni_p" = "red", "BH_p"="blue"))

# Filter the data to retain just the 5% ~most significant   
dim(GSE215865_stat_results)
GSE215865_stat_results_subselected <- GSE215865_stat_results |> 
  rownames_to_column("GSE215865_rowname") |> # do this to retain the rownanes somewhere in the matrix
  filter( bonferroni_p < 0.1 ) # Bonferroni may be a little too strict for this dataset
dim(GSE215865_stat_results_subselected)
dim(GSE215865_stat_results_subselected)[1]/dim(GSE215865_stat_results)[1]*100 # find a value that retains ~5% most significant
# put the rownames back where they belong
GSE215865_stat_results_subselected <- column_to_rownames(GSE215865_stat_results_subselected, var="GSE215865_rowname")

# export results with stats
export_data(data_object = GSE215865_stat_results_subselected, file_name = "GSE215865_stat_results_subselected.txt")
# remove the columns that contain the stats(this is hacky)
GSE215865_stat_results_subselected_and_cleaned <- remove_last_n_columns(GSE215865_stat_results_subselected, n=7)
# export the data ready to create a HD
export_data(data_object = GSE215865_stat_results_subselected_and_cleaned, file_name = "GSE215865_stat_results_subselected_and_cleaned.txt")

# heatmap dendrograms -----------------------------------------------------

# Here I originally ran into a bit of trouble.
# Heatmap didn't work for a while -- seems like issue was unwanted characters
# in the column name that R did not directly complain about. Per note below
# removed all of the potentially offending characters and it worked
# from https://github.com/Teichlab/cellphonedb/issues/219
# "I replaced the symbols(/-_ et al) with dot(.) in cell/cluster names, and then the heatmap plot function works again. May it can help you."

GSE215865_metadata <- import_metadata("GSE215865.metadata.txt")
# We replace the existing column names with syntactically correct ones like this
colnames(GSE215865_metadata) <- make.names(colnames(GSE215865_metadata))
# Save the metadata back to file
write.table(file = "GSE215865.metadata.txt", x = GSE215865_metadata, sep="\t",  quote=FALSE, col.names=NA)


# First, look at all the data with HD
#tic()
#heatmap_dendrogram(file_in = "GSE215865.data.txt.quantile.PREPROCESSED.txt",
#                   metadata_table = "GSE215865.metadata.txt",
#                   metadata_column="covid.19_positive"
#)
#toc() # This may fail depending on your system memory. It failed for me,
# so I've commented this code out.
#system("open GSE215865_stat_results_subselected_and_cleaned.txt.HD.png")


# Now look at just the statistically subselected genes
heatmap_dendrogram(file_in = "GSE215865_stat_results_subselected_and_cleaned.txt",
                   metadata_table = "GSE215865.metadata.txt",
                   metadata_column="covid.19_positive"
)
# This worked fine on my local machine
system("open GSE215865_stat_results_subselected_and_cleaned.txt.HD.png")

# ANNOTATE RESULTS --------------------------------------------------------

# Got the current ENSG annotations from here
# https://useast.ensembl.org/biomart/martview/
# There is a simple to use user interface that allows you to choose the data to download
# Choose Database --> Ensembl Genes
# Choose Dataset --> Human Genes (GRCh38.p14)
# Click on "Attributes" and then "GENE"
# Check "Gene stable ID", "Gene description", "Gene name"
# Click on the "Results" button
# Export all results to "File" "TSV"
# Click on "Go"
# This will download a file called "mart_export.txt" that we use below
# I assume that you move "mart_export.txt" to this directory before proceeding
# The goal here is to add two columns to the stat results that contain the 
# gene name and gene descirption for each ENSG id 
annotations <- read.table(file="mart_export.txt",row.names=NULL,header=TRUE,sep="\t", # compared to import_data, changed row.names from 1 to NULL
                          colClasses = "character", check.names=FALSE,
                          comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE)

GSE215865_stat_results_subselected <- import_data("GSE215865_stat_results_subselected.txt")
# The ENSG annotations in the data set have the version tag.
# Remove the version tag and you can match to the Gene stable ID
GSE215865_genes <- rownames(GSE215865_stat_results_subselected)
GSE215865_genes[1]
GSE215865_genes <- gsub(pattern="\\..*", replacement = "", x=GSE215865_genes) # have a .num tag, part of the stable ID version
GSE215865_genes[1]

# add columns for the annotation information, NA to start
GSE215865_stat_results_subselected <- as.data.frame(GSE215865_stat_results_subselected)
GSE215865_stat_results_subselected$gene_descriptions <- NA 
GSE215865_stat_results_subselected$gene_names <- NA

# now use a simple loop to look up and add the annotation values to the stat results
for ( i in 1:length(GSE215865_genes)){ 
  GSE215865_gene_annotation <- annotations[annotations$`Gene stable ID` == GSE215865_genes[i],]
  GSE215865_stat_results_subselected[i,"gene_names"] <- GSE215865_gene_annotation[1,"Gene name"]
  GSE215865_stat_results_subselected[i,"gene_descriptions"] <- GSE215865_gene_annotation[1,"Gene description"]
}
write.table(file = "GSE215865.annotations.txt", x = GSE215865_stat_results_subselected, sep="\t", quote=FALSE, col.names=NA) 

# PRELIMINARY PATHWAY ANALYSIS --------------------------------------------

# Use the list of gene ids from above
GSE215865_genes

# Perform pathway analysis and generate an interactive visualization
gostres <- gost(query = GSE215865_genes, organism = 'hsapiens', significant = FALSE)
gostplot(gostres, capped = TRUE, interactive = TRUE)

# JUNK BELOW HERE ---------------------------------------------------------
# JUNK BELOW HERE ---------------------------------------------------------
# JUNK BELOW HERE ---------------------------------------------------------

# https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&search=covid&display=20&zsort=samples
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215865
# Submission date	Oct 16, 2022
# Last update date	Feb 08, 2023
# Contact name	Ryan Conrad Thompson
# E-mail(s)	rct@thompsonclan.org
# Organization name	Icahn School of Medicine at Mount Sinai
# Department	Charles Bronfman Institute for Personalized Medicine
# Lab	Beckmann Lab
# Street address	1 Gustave L. Levy Place
# City	New York
# State/province	NY
# ZIP/Postal code	10029-5674
# Country	USA