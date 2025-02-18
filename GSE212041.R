# Data were found by a search in:
# https://www.ncbi.nlm.nih.gov/geo/browse/
# The keyword "covid" was used and studies were sorted by "Series Type"
# to find "Expression profiling by high throughput sequencing", and by
# "Samples" to find the studies with the largest number of samples.
# As of 1/15/2024 three studies were found with more than 500 samples 
# The analysis below is for one of those three datasets.
# All data in this analysis come from:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212041
# Note that the metadata used for this analysis was not readily available. The metadata
# were obtained from the accession pages as outlined below.

# Load all required packages ----------------------------------------------
# Note, Dependency installation on Ubuntu : sudo apt-get install libzmq3-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev build-essential libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev

library(readxl)
library(tidyverse)
library(gprofiler2)
library(plotly)
library(RCurl)
library(R.utils)
library(preprocessCore) # To install this package, library(BiocManager), BiocManager::install("preprocessCore")
library(matlab) # To install this package, library(BiocManager), BiocManager::install("matlab")
library(ecodist) # To install this package, library(BiocManager), BiocManager::install("ecodist")

# Source all additional functions from Kevin's github repository ----------------------------------------------

source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_data.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_metadata.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/preprocessing_tool.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/export_data.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/calculate_pco.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/refs/heads/master/render_calculated_pcoa.2-17-2025.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/sigtest.R")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/remove_last_n_columns.R")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/heatmap_dendrogram.r")

# Create a directory for working (if it doesn't already exist) and move to it --------
# Specify the directory path
dir_path <- "~/GSE212041/"
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
data_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212041/suppl/GSE212041_Neutrophil_RNAseq_TPM_Matrix.txt.gz"
# Specify the local file path to save the downloaded content
local_file_path_data <- "GSE212041_Neutrophil_RNAseq_TPM_Matrix.txt.gz"
# Download the file
download.file(data_url, destfile = local_file_path_data, mode = "wb")
# Now we have to unzip the gz file, easiest to just do this with a system call
gunzip(local_file_path_data, remove = TRUE) # unzip the file and remove the orignal gz file

# Now the metadata ...
# Metadata for this dataset had to be collected from multiple accession pages.
# The following procedure was used.
# Go to this page:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212041
# Under "Samples (781)" click on the "More" link
# At the bottom of the accession list click on the following link
# "Accession list truncated, click here to browse through all related public accessions"
# This will bring you to "https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=212041"
# Click on the "Export" button at the top of the page
# Under "Amount" select "All search results", Under "Format" select "CSV"
# Now click on the "Export" button
# This will download a file called "sample.csv" that contains the accession numbers for 781 samples
# Move that file into this directory and proceed below.
current_path <- "C:/Users/[your_username]/Downloads/sample.csv"
# Define the new path of the file
new_path <- paste(dir_path, "sample_GSE212041.csv", sep = "")# rename the sample table with something less ambiguous 
# Move (rename) the file
file.rename(current_path, new_path)

# Load file that has sample accession information
GSE212041_metadata <- read_csv("sample_GSE212041.csv")

# add columns to the existing metadata to accept new values that will be collected from the accession sites
GSE212041_metadata <- GSE212041_metadata |>
  mutate( time=NA ) |>
  mutate( patient_category=NA ) |>
  mutate( covid_status=NA ) |>
  mutate( acuity_max=NA )

# Loop that uses RCurl to collect the values of interest from accession websites
# The details of parsing the pages is a little beyond the scope of this analysis,
# but briefly, individual accession pages were inspected to identify the desired metadata.
# The script below was written to collect the desired values from the rest of the
# accession pages. Note that this tool was specifically written for these meadata and is not
# intended as a general solution. You can write something far more elegant.
for (i in 1:nrow(GSE212041_metadata)){
  # Get the accession numbers
  GSE212041_accession <- GSE212041_metadata[i,"Accession"]
  # Use the accession numbers to query sites that contain the data
  GSE212041_url <- paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GSE212041_accession, sep="")
  print( paste("i", i, GSE212041_url) )
  # Collect data for a single accession
  GSE212041_url_content <- getURL(GSE212041_url) 
  # Parse the page by line
  GSE212041_line_sep_content <- strsplit(GSE212041_url_content, split="\n")
  # Parse further based on presence of unique text, "covid-19"
  for (j in 1:length(GSE212041_line_sep_content[[1]])){
    if( grepl(pattern="covid-19", x = GSE212041_line_sep_content[[1]][j]) ==TRUE ){
      print(paste("j", j, GSE212041_line_sep_content[[1]][j]))
      # Parse further based on line breaks
      GSE212041_br_split_content <- strsplit( GSE212041_line_sep_content[[1]][j], split="<br>"  ) 
      GSE212041_time <- strsplit( GSE212041_br_split_content[[1]][2], split=": " )[[1]][2]
      GSE212041_patient_category <- strsplit( GSE212041_br_split_content[[1]][3], split=": " )[[1]][2]
      GSE212041_covid_status <- strsplit( GSE212041_br_split_content[[1]][4], split=": " )[[1]][2] 
      GSE212041_acuity_max <- strsplit( GSE212041_br_split_content[[1]][5], split=": " )[[1]][2]
      # Load parsed values into the metadata tibble
      GSE212041_metadata[i, "time"] <- GSE212041_time
      GSE212041_metadata[i, "patient_category"] <- GSE212041_patient_category
      GSE212041_metadata[i, "covid_status"] <- GSE212041_covid_status
      GSE212041_metadata[i, "acuity_max"] <- GSE212041_acuity_max
    } 
  }
} # This loop takes a few minutes to run; there are 781 sites to be queried

# Use the "Title" column as the rownames for the tibble
GSE212041_metadata <- column_to_rownames(GSE212041_metadata, var="Title")
# Check to see how many unique values there are for "patient_category"
unique( GSE212041_metadata$patient_category ) # Three, "COVID+", "COVID- symptomatic", and "healthy"   

# Now load the data
GSE212041_data <- import_data("GSE212041_Neutrophil_RNAseq_TPM_Matrix.txt")

# See how many samples there are at the start
ncol(GSE212041_data) - 1 # 781, first column has "Symbol", not sample data 

# Check that sample ids are the same in data and metadata
length( intersect(rownames(GSE212041_metadata), colnames(GSE212041_data)) ) # 781, They all appear to be in common
identical(rownames(GSE212041_metadata), colnames(GSE212041_data)) # FALSE, But, there is a difference
dim(GSE212041_metadata) # 781 15
dim(GSE212041_data) # 60640 782 -- data for 781 samples plus an additional column
setdiff(colnames(GSE212041_data), rownames(GSE212041_metadata)) # "Symbol" is the column that does not exist in the metadata
# See what "Symbol" is --  
GSE212041_data[1:3,1:3] # looks like an additional ID for the ENSG number? Not entirely sure.
# I do not think I need it to proceed so I will remove it.
# remove "Symbol" column from GSE212041_data ...
# Specify the column name to be removed
column_to_remove <- "Symbol"
# Find the index of the column to be removed
col_index <- which(colnames(GSE212041_data) == column_to_remove)
# Remove the specified column
GSE212041_data <- GSE212041_data[, -col_index]
dim(GSE212041_data) # 60640   781
# check for identity again
length( intersect(rownames(GSE212041_metadata), colnames(GSE212041_data)) ) # 781, They are all the same
identical(rownames(GSE212041_metadata), colnames(GSE212041_data)) # FALSE,  But the order must be different
# Sort the data and metadata by the samples ID's taken from the colnames of the data
sorted_indices <- sort( colnames(GSE212041_data) )
# Use the sorted indices to sort the data by column
GSE212041_data <- GSE212041_data[,sorted_indices]
# And sort the metadata by row
GSE212041_metadata <- GSE212041_metadata[sorted_indices,]
# Check the ordering one more time
identical(rownames(GSE212041_metadata), colnames(GSE212041_data)) # Now it is TRUE, time to save both 

# SAVE THE DATA AND THE METADATA TO FILE ----------------------------------
# First the data
export_data(data_object = GSE212041_data, file_name = "GSE212041.data.txt")
# Now the metadata
export_data(data_object = GSE212041_metadata, file_name = "GSE212041.metadata.txt")

# OPTIONAL Clean House ----------------------------------------------------

# At this point you may want to clean house the dreaded "rm(list = ls())" command
# Note the if you do, you will have to run the bit of code from above that sources
# functions from github. Sorry they're not in a package as of yet.

# FUN STUFF ---------------------------------------------------------------

# The analysis really starts at the this point, everything up until now has been
# data wrangling to make the analyses below possible

# Import the data and metadata  -------------------------------------------
GSE212041_data <- import_data("GSE212041.data.txt")
GSE212041_metadata <- import_metadata("GSE212041.metadata.txt")

# Look at distribution of raw data ----------------------------------------

# Look at TPM data distribution
# Explore distribution of data
all_GSE212041_data <- as.vector(GSE212041_data)
hist(all_GSE212041_data, breaks =100)
hist(log10(all_GSE212041_data), breaks =100) # The data appear to be roughly log normal
# Look at with ggplot
all_GSE212041_data <- data.frame(Values = as.vector(GSE212041_data))
ggplot(data=all_GSE212041_data,mapping=aes(x=Values))+
  geom_histogram(bins = 100, fill = "purple", color = "black") +
  scale_x_log10() +
  labs(title = "Histogram of expression values",
  x = "Log-Transformed Values",
  y = "Frequency")

# Look at the distribution of summary values in the raw data
GSE212041_raw_mean <- apply(GSE212041_data, 2, mean) 
GSE212041_raw_median <- apply(GSE212041_data, 2, median) 
GSE212041_raw_min <- apply(GSE212041_data, 2, min) 
GSE212041_raw_max <- apply(GSE212041_data, 2, max) 
GSE212041_raw_sd <- apply(GSE212041_data, 2, sd)

plot(GSE212041_raw_mean, pch=20) # no extreme outliers
plot(GSE212041_raw_median, pch=20) # a few outliers, roughly those above 1
plot(GSE212041_raw_min, pch=20) # 0, no surprises
plot(GSE212041_raw_max, pch=20) # a few fairly extreme outliers, those above 3e+05
plot(GSE212041_raw_sd, pch=20) # values above 1500 are somewhat anomalous 

# Make these data tidy for ggplot
GSE212041_raw_descriptive_stats <- as.matrix(names(GSE212041_raw_mean))
GSE212041_raw_descriptive_stats <- cbind(GSE212041_raw_descriptive_stats, as.matrix(GSE212041_raw_mean))
GSE212041_raw_descriptive_stats <- cbind(GSE212041_raw_descriptive_stats, as.matrix(GSE212041_raw_median))
GSE212041_raw_descriptive_stats <- cbind(GSE212041_raw_descriptive_stats, as.matrix(GSE212041_raw_min))
GSE212041_raw_descriptive_stats <- cbind(GSE212041_raw_descriptive_stats, as.matrix(GSE212041_raw_max))
GSE212041_raw_descriptive_stats <- cbind(GSE212041_raw_descriptive_stats, as.matrix(GSE212041_raw_sd))
colnames(GSE212041_raw_descriptive_stats) <- c("GSE212041_sample_ids","GSE212041_mean", "GSE212041_median", "GSE212041_min", "GSE212041_max", "GSE212041_sd")
rownames(GSE212041_raw_descriptive_stats) <- names(GSE212041_raw_mean)
GSE212041_raw_descriptive_stats <- as.data.frame(GSE212041_raw_descriptive_stats)
GSE212041_raw_descriptive_stats$GSE212041_mean <- as.numeric(GSE212041_raw_descriptive_stats$GSE212041_mean)
GSE212041_raw_descriptive_stats$GSE212041_median <- as.numeric(GSE212041_raw_descriptive_stats$GSE212041_median)
GSE212041_raw_descriptive_stats$GSE212041_min <- as.numeric(GSE212041_raw_descriptive_stats$GSE212041_min)
GSE212041_raw_descriptive_stats$GSE212041_max <- as.numeric(GSE212041_raw_descriptive_stats$GSE212041_max)
GSE212041_raw_descriptive_stats$GSE212041_sd <- as.numeric(GSE212041_raw_descriptive_stats$GSE212041_sd)

ggplot(
  data=GSE212041_raw_descriptive_stats,
  mapping=aes(x=1:nrow(GSE212041_raw_descriptive_stats))) +
  scale_y_log10() +
  geom_line(aes(y=GSE212041_mean, colour="GSE212041_mean")) +
  geom_line(aes(y=GSE212041_median, colour="GSE212041_median")) +
  geom_line(aes(y=GSE212041_min, colour="GSE212041_min")) +
  geom_line(aes(y=GSE212041_max, colour="GSE212041_max")) +
  geom_line(aes(y=GSE212041_sd, colour="GSE212041_sd")) +
  scale_color_manual(name = "Stat Values", values = c("GSE212041_mean" = "black", "GSE212041_median" = "red", "GSE212041_min"="blue", "GSE212041_max"="green", "GSE212041_sd"="purple"))
# The outliers don't look so bad in log space, for now will keep all samples

# Preprocess and look at distribution of normed data -----------------------------------------------------

# Preprocess the data
preprocessing_tool("GSE212041.data.txt")

# Look at distribution of normed data
GSE212041_data_n <- import_data("GSE212041.data.txt.quantile.PREPROCESSED.txt")
all_GSE212041_data_n <- as.vector(GSE212041_data_n)
hist(log10(all_GSE212041_data_n), breaks =100) # no longer looks log normal

# look at the distribution of summary values in the normed data
GSE212041_norm_mean <- apply(GSE212041_data_n, 2, mean)
GSE212041_norm_median <- apply(GSE212041_data_n, 2, median) 
GSE212041_norm_min <- apply(GSE212041_data_n, 2, min) 
GSE212041_norm_max <- apply(GSE212041_data_n, 2, max) 
GSE212041_norm_sd <- apply(GSE212041_data_n, 2, sd)

plot(GSE212041_norm_mean) # outliers less than 0.10
plot(GSE212041_norm_median) # outliers above 0.005, some clear structure/tiering in the data
plot(GSE212041_norm_min) # outliers above and below structured data
plot(GSE212041_norm_max) # Given treatment, all have max of 1.
plot(GSE212041_norm_sd) # majority in the same range, a scattered few are below
# Overall assessment - leave the data alone for now. Not sure how to interpret the structure
# or outliers in the data. Let's see if outliers correlate with metadata at all - to see if they
# are the product of batch effects.

# PCoA --------------------------------------------------------------------

# calculate PCoA on the preprocessed data
calculate_pco(file_in="GSE212041.data.txt.quantile.PREPROCESSED.txt")

# Iterate through the metadata creating a static 3d PCoA for each metadata column
plot_static_colored_3d_pcoas(
  pcoa_filename = "GSE212041.data.txt.quantile.PREPROCESSED.txt.euclidean.PCoA",
  selected_eigen_vectors = c(1,2,3),
  metadata_filename = "GSE212041.metadata.txt",
  debug = TRUE
)
# Open all of the the generataed PCoA images.
# system("open *.PCoA.*.png") # prior to windows 11 this worked, 
# with windows 11
# something like this - please let me know if you find a different approach
my_pngs <- dir(pattern = ".*\\.png$")
for (i in my_pngs){
  my_command <- paste("start", i, sep=)
  print(my_command)
  shell(my_command)
}
# Clear horseshoe artifact. Correlations not obvious. Some pattern for "patient_category"
# and "covid_status", but the separation is merky at best.

# render interactive 3d PCoA plot from the PCoA and the corresponding metadata
plot_interactive_colored_3d_pcoa(
  pcoa_data_file = "GSE212041.data.txt.quantile.PREPROCESSED.txt.euclidean.PCoA",
  selected_eigen_vectors = c(1,2,3),
  pcoa_metadata_file = "GSE212041.metadata.txt",
  metadata_column = "patient_category" # "covid_status"
)
# There are several metadata/data correlations that can be observed with the colored PCoAs.
# There is clear correlation with "time".
# "covid_status" and "patient_category" exhibit an overlapping correlation.
# "acuity_max" also has a discernible pattern.
# Further examination of all these correlations is warranted; however, we're just interested in
# metadata that unambiguously express Covid related disease state.
# Therefore, I selected to use "patient_category" as the group designator for the statistical test below.
# It has three groups, but these are unambiguously labeled, whereas 
# "covid_status" is a boolean. Not completely sure what 0/1 represent, 
# although I can guess based on overlap with "patient_category".
# see unique(GSE212041_metadata[,"patient_category"]) vs unique(GSE212041_metadata[,"covid_status"])
# There is a fairly obvious horseshoe artifact in the data. We try to address this below by 
# using another distance metric.

# See if the horseshoe can be improved by changing the distance metric to Bray-Curtis
calculate_pco(file_in="GSE212041.data.txt.quantile.PREPROCESSED.txt", dist_method = "bray-curtis")

# Render the Bray-Curtis calculated PCoA
plot_interactive_colored_3d_pcoa(
  pcoa_data_file = "GSE212041.data.txt.quantile.PREPROCESSED.txt.bray-curtis.PCoA",
  selected_eigen_vectors = c(1,2,3),
  pcoa_metadata_file = "GSE212041.metadata.txt",
  metadata_column = "patient_category" # "covid_status"
)
# Does not help with the horseshoe, suggestions to address this are welcome.

# perform a stat test on the data -----------------------------------------

# Perform a stat test on the data/metadata (selected metadata column is used to create groups for chosen stat test)
# After normalization/preprocessing data exhibit a non-normal distribution. There are three groups. 
# Kruskal-Wallis is an appropriate test given these criteria.
sigtest(data_file="GSE212041.data.txt.quantile.PREPROCESSED.txt", 
        metadata_file="GSE212041.metadata.txt",  
        metadata_column="patient_category",
        stat_test="Kruskal-Wallis",
        p_adjust_method = "BH"
)

# Load the stat test results and subselect data based on the stat results
GSE212041_stat_results <- import_data("GSE212041.data.txt.quantile.PREPROCESSED.Kruskal-Wallis.patient_category.STAT_RESULTS.txt")
GSE212041_stat_results <- as_tibble(GSE212041_stat_results,rownames=NA) # convert to tibble BUT keep the rownames

# Take a look at the stat results
ggplot(
  data=GSE212041_stat_results,
  mapping=aes(x=1:length(bonferroni_p))) +
  geom_line(aes(y=p, colour="p")) +
  geom_line(aes(y=bonferroni_p, colour="bonferroni_p")) +
  geom_line(aes(y=BH_p, colour="BH_p")) +
  scale_color_manual(name = "Stat Values", values = c("p" = "black", "bonferroni_p" = "red", "BH_p"="blue"))
# In this study Bonferroni appears to be a very stringent filter. Our goal is to fine as robust a set
# of genes related to Covid-19 expression as possible, As such, I've decided to still use Bonferroni
# but reduce the stringency to allow for a reasonable number of genes to pass our filtering.

# Filter the data to retain just the top ~5% ~most significant 
dim(GSE212041_stat_results) # 50811   788
GSE212041_stat_results_subselected <- GSE212041_stat_results |> 
  rownames_to_column("GSE212041_rowname") |> # do this to retain the rownanes somewhere in the matrix
  filter( bonferroni_p < 0.1 ) # Bonferroni seems to be a little too stringent for this dataset
dim(GSE212041_stat_results_subselected) #  994 789
dim(GSE212041_stat_results_subselected)[1]/dim(GSE212041_stat_results)[1]*100 # find a value that retains 5% most significant
# The parameters used retain just about 2% of the genes; however, I don't want to use a Bf that is any
# less stringent. You can adjust this parameter or use the BH p or recalculate the stats with another
# p adjustment method ( see the p_adjust_method argument in sigtest() ).
# Put the rownames back where they belong
GSE212041_stat_results_subselected <- column_to_rownames(GSE212041_stat_results_subselected, var="GSE212041_rowname")

# Export results with stats
export_data(data_object = GSE212041_stat_results_subselected, file_name = "GSE212041_stat_results_subselected.txt")
# Remove the columns that contain the stats
GSE212041_stat_results_subselected_and_cleaned <- remove_last_n_columns(GSE212041_stat_results_subselected, n=7)
# Export the data ready to create a HD
export_data(data_object = GSE212041_stat_results_subselected_and_cleaned, file_name = "GSE212041_stat_results_subselected_and_cleaned.txt")

# Heatmap-dendrogram -----------------------------------------------------

# Produce a heatmap-dendrogram of the entire dataset
# heatmap_dendrogram(file_in = "GSE212041.data.txt.quantile.PREPROCESSED.txt",
#                    metadata_table = "GSE212041.metadata.txt",
#                    metadata_column="patient_category"
# )
# This requires computational resources and time(). Given your particular system, it may fail.
# It failed on my laptop, hence commenting the code out.
# system("open GSE212041_stat_results_subselected_and_cleaned.txt.HD.png")

# Produce a heatmap dendrogram of the statistically subselected data - this works with no issues
heatmap_dendrogram(file_in = "GSE212041_stat_results_subselected_and_cleaned.txt",
                   metadata_table = "GSE212041.metadata.txt",
                   metadata_column="patient_category"
)
# system("open GSE212041_stat_results_subselected_and_cleaned.txt.HD.png") # This worked prior to windows 11
shell("start GSE212041_stat_results_subselected_and_cleaned.txt.HD.png") # This works under windows 11

# Calculate and render a PCoA on the statistically subselected data
calculate_pco(file_in="GSE212041_stat_results_subselected_and_cleaned.txt")

# Render the statistically subselected PCoA
plot_interactive_colored_3d_pcoa(
  pcoa_data_file = "GSE212041_stat_results_subselected_and_cleaned.txt.euclidean.PCoA",
  selected_eigen_vectors = c(1,2,3),
  pcoa_metadata_file = "GSE212041.metadata.txt",
  metadata_column = "patient_category" # "covid_status"
)
# Separation of the diseased and healthy groups is much more pronounced

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
# Define the current path of the file
current_path <- "C:/Users/[your_username]/Downloads/"
# Define the new path of the file
new_path <- paste(dir_path, "mart_export.txt", sep = "")#
# Move (rename) the file
file.rename(current_path, new_path)

# The goal here is to add two columns to the stat results that contain the 
# gene name and gene description for each ENSG id 
annotations <- read.table(file="mart_export.txt",row.names=NULL,header=TRUE,sep="\t", # compared to import_data, changed row.names from 1 to NULL
                          colClasses = "character", check.names=FALSE,
                          comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE)

GSE212041_stat_results_subselected <- import_data("GSE212041_stat_results_subselected.txt")
# The ENSG annotations in the data set have the version tag.
# Remove the version tag and you can match to the Gene stable ID
GSE212041_genes <- rownames(GSE212041_stat_results_subselected)
GSE212041_genes[1] # no version tag, can skip this step -- commented out the line below
#GSE212041_genes <- gsub(pattern="\\..*", replacement = "", x=GSE212041_genes) # have a .num tag, part of the stable ID version

# add columns for the annotation information, NA to start
GSE212041_stat_results_subselected <- as.data.frame(GSE212041_stat_results_subselected)
GSE212041_stat_results_subselected$gene_descriptions <- NA 
GSE212041_stat_results_subselected$gene_names <- NA

# now use a simple loop to look up and add the annotation values to the stat results
for ( i in 1:length(GSE212041_genes)){ 
  GSE212041_gene_annotation <- annotations[annotations$`Gene stable ID` == GSE212041_genes[i],]
  GSE212041_stat_results_subselected[i,"gene_names"] <- GSE212041_gene_annotation[1,"Gene name"]
  GSE212041_stat_results_subselected[i,"gene_descriptions"] <- GSE212041_gene_annotation[1,"Gene description"]
}
write.table(file = "GSE212041.annotations.txt", x = GSE212041_stat_results_subselected, sep="\t", quote=FALSE, col.names=NA) 

# PRELIMINARY PATHWAY ANALYSIS --------------------------------------------

# Use the list of gene ids from above
GSE212041_genes

# Perform pathway analysis and generate an interactive visualization
gostres <- gost(query = GSE212041_genes, organism = 'hsapiens', significant = FALSE)
gostplot(gostres, capped = TRUE, interactive = TRUE)

# That's it for the analysis of GSE212041.
# Continue on to one of the other three datasets and combine_covid when you're done with all three

# JUNK BELOW HERE ---------------------------------------------------------
# JUNK BELOW HERE ---------------------------------------------------------
# JUNK BELOW HERE ---------------------------------------------------------

# Use RCurl to collect metadata from accession sites
# install.packages("RCurl")
# library(RCurl)
# ls(package:RCurl)
# html_content <- getURL(url)
# test <- strsplit(html_content, split="\n")