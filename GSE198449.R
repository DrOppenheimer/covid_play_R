# Data were found by a search in:
# https://www.ncbi.nlm.nih.gov/geo/browse/
# The keyword "covid" was used and studies were sorted by "Series Type"
# to find "Expression profiling by high throughput sequencing", and by
# "Samples" to find the studies with the largest number of samples.
# As of 1/15/2024 three studies were found with more than 500 samples 
# The analysis below is for one of those three datasets.
# All data in this analysis come from:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198449
# The script below goes through my entire analysis, from download to final figures
# The only step that cannot be completed by following this workflow is to download
# annotations from Ensembl. Directions to do this are below in the "Annotate Results" section

# Load all required packages ----------------------------------------------

library(readxl)
library(tidyverse)
library(gprofiler2)
library(plotly)
library(R.utils)

# Source all additional functions from Kevin's github repository ----------------------------------------------

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
dir_path <- "~/GSE198449/"
# Check if the directory exists
if (dir.exists(dir_path)) {
  stop("Error: The directory already exists.")
} else {
  # Create the directory if it does not exist
  dir.create(dir_path)
  print("Directory created successfully.")
}
# Move to that directory 
setwd(dir_path)

# Download data and metadata ----------------------------------------------

# First the data ...
# Specify the FTP URL of the file
data_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198449/suppl/GSE198449_featureCounts.txt.gz"
# Specify the local file path to save the downloaded content
local_file_path_data <- "GSE198449_featureCounts.txt.gz"
# Download the file
download.file(data_url, destfile = local_file_path_data, mode = "wb")
# now we have to unzip the gz file, easiest to just do this with a system call
gunzip(local_file_path_data, remove = TRUE) # unzip the file and remove the orignal gz file

# Now the metadata
metadata_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198449/suppl/GSE198449_PCR_and_Symptoms_Full_Info.xlsx"
local_file_path_metadata <- "GSE198449_PCR_and_Symptoms_Full_Info.xlsx"
download.file(metadata_url, destfile = local_file_path_metadata, mode = "wb")

# Import and explore the metadata -----------------------------------------

# Import metadata
GSE198449_metadata <- read_excel(path="GSE198449_PCR_and_Symptoms_Full_Info.xlsx", skip=21, na="NA")
# see how many samples there are to start
dim(GSE198449_metadata) # 3049 6
# edit the colnames to make them syntactically correct
colnames(GSE198449_metadata)
colnames(GSE198449_metadata) <- make.names(colnames(GSE198449_metadata))
colnames(GSE198449_metadata)
# NOTE: had to do this because later read.table will covert " ", "-", "(", or ")" to "." 
# when it reads in metadata that I have saved to file. Apparently read_excel does not 
# have this issue. There is an option to prevent this (check.names=FALSE in read.table), 
# but it doesn't seem to work - found some sparse documentation to support this.
# Also the HD function fails on column names in metadata that contain a "_"

# Only retain metadata for patients with a PCR test
unique(GSE198449_metadata[,"PCR.test.for.SARS.Cov.2"])
GSE198449_filtered_metadata <- GSE198449_metadata |>
  filter(!is.na(PCR.test.for.SARS.Cov.2))
dim(GSE198449_metadata) # 3049 6
dim(GSE198449_filtered_metadata) # 2830    6, lost 219 samples 

# Take a look at the PCR test results
ggplot(
  GSE198449_filtered_metadata,
  aes(x=PCR.test.for.SARS.Cov.2))+
  geom_bar() # vast majority are "Detected" or "Not"

# Remove samples from metadata that are not "Detected" or "Not"
GSE198449_filtered_metadata_culled <- GSE198449_filtered_metadata |>
  filter(PCR.test.for.SARS.Cov.2 == "Detected"|PCR.test.for.SARS.Cov.2 == "Not")
# check
dim(GSE198449_filtered_metadata_culled) # 2775    6, lost 55 more samples
ggplot(
  GSE198449_filtered_metadata_culled,
  aes(x=PCR.test.for.SARS.Cov.2))+
  geom_bar()

# Take a look at the metadata to see if there are any obvious questions we can address
colnames(GSE198449_filtered_metadata_culled)

# How do "enrollment.batch" and "PCR.test.for.SARS.Cov.2" relate
ggplot(
  GSE198449_filtered_metadata_culled,
  aes(x=enrollment.batch,fill=PCR.test.for.SARS.Cov.2))+
  geom_bar() # interesting - detected are all in later surveillance

# How do "sample.collection.time.point..days.since.T0." and "PCR.test.for.SARS.Cov.2" relate
ggplot(
  GSE198449_filtered_metadata_culled,
  aes(x=sample.collection.time.point..days.since.T0.,fill=PCR.test.for.SARS.Cov.2))+
  geom_bar() # complex relationship, beyond what we want to accomplish

# See the uniqueness of the symptoms
ggplot(
  GSE198449_filtered_metadata_culled,
  aes(x=`symptom`))+
  geom_bar()
dim( unique( GSE198449_filtered_metadata_culled[,"symptom"] ) ) # 260 unique combinations of symptoms, very complicated

# Sort metadata by sample id
GSE198449_filtered_metadata_culled_sorted <- GSE198449_filtered_metadata_culled |>
  arrange(biocollection.ID)
GSE198449_filtered_metadata_culled_sorted[1:3,]

# Import and explore the data ------------------------------

# Start by looking at the first line
GSE198449_data_file <- "GSE198449_featureCounts.txt"
readLines(GSE198449_data_file, n = 1) # tab separated with a rowname in the first column

# Import data using the first column ("Geneid") as the rowname
GSE198449_data <- import_data(GSE198449_data_file)
# see how many samples we're starting with
ncol(GSE198449_data) # 1858

# Look at the distribution of the raw data
hist( as.vector(GSE198449_data) )
hist( log10(as.vector(GSE198449_data)) ) # Roughly log normal

# Figure out how important the plate stat is ------------------------------
# i.e. I want to see if there is obvious bias introduced that corresponds to each plate
# The plate information is in each column name of the data
colnames(GSE198449_data)[1] # It's the "P1" part of this example

# Get the part of the date colname in the data that comes after the "-"
GSE198449_data_colnames <- colnames(GSE198449_data) |>
  strsplit("-")
dim(GSE198449_data_colnames) # it's a list
GSE198449_data_colnames <- as.data.frame(GSE198449_data_colnames) # cønvert to data frame
GSE198449_data_colnames <- t(GSE198449_data_colnames) # transpose the df
dim(GSE198449_data_colnames) # 1858 2
GSE198449_data_colnames[1:3,] # look at the first few

# Now get the "T" and "P" portions of the sample name separated
GSE198449_p <- GSE198449_data_colnames[,2] |> # this creates a list
  strsplit("_")
GSE198449_p <- as.data.frame(GSE198449_p) # convert to dataframe
GSE198449_p <- t(GSE198449_p) # transpose
colnames(GSE198449_p) <- c("T", "P") # add convenient column names

# See how many plates there are
unique(GSE198449_p[,"P"]) 
length(unique(GSE198449_p[,"P"])) # there are 27 plates

# and this is how they are distributed about the samples
GSE198449_p <- as.data.frame(GSE198449_p)
ggplot(
  GSE198449_p,
  aes(x=P)) +
  geom_bar() # seems pretty random with respect to "Detected" and "Not"

# Try to get IDs to match between the data and metadata, cull each accordingly ------------------------
dim(GSE198449_data) # 58929  1858
colnames(GSE198449_data)[1:3]
# "20_5063-T00_P1" "20_5134-T00_P1" "20_5218-T00_P1"
dim(GSE198449_filtered_metadata_culled_sorted) #2775    6
GSE198449_filtered_metadata_culled_sorted[1:3,]
# A tibble: 3 × 6
#enrollment.batch   biocollection.ID participant.ID sample.collection.time.point..days.since.T0. PCR.test.for.SARS.Cov.2 symptom
#<chr>              <chr>            <chr>          <chr>                                        <chr>                   <chr>  
# 1 Later Surveillance 20_5002-T00      5002           00                                           Not                     NA     
# 2 Later Surveillance 20_5002-T07      5002           07                                           Not                     NA     
# 3 Later Surveillance 20_5002-T14      5002           14                                           Not                     NA 
# The assumption I've made looking at the data and metadata is that the column name of the data corresponds
# to the "biocollection.ID" column of the metadata. However, "biocollection.ID" has an additional tag for
# plate that we will have to take care of

# First, make data colnames syntactically correct
colnames(GSE198449_data)[1:3]
colnames(GSE198449_data) <- make.names( colnames(GSE198449_data) )
colnames(GSE198449_data)[1:3]

# Now make the biocollection.ID values in the metadata syntactically correct
GSE198449_filtered_metadata_culled_sorted[1:3,"biocollection.ID"]
# GSE198449_filtered_metadata_culled_sorted[,"biocollection.ID"] <- make.names( GSE198449_filtered_metadata_culled_sorted[,"biocollection.ID"] )
# NOTE: The above line does not work like I would expect. It produced a huge string and applies it to each row
# instead do this, it works as expected
for (i in 1:nrow(GSE198449_filtered_metadata_culled_sorted[,"biocollection.ID"])){
  GSE198449_filtered_metadata_culled_sorted[i,"biocollection.ID"] <- make.names( GSE198449_filtered_metadata_culled_sorted[i,"biocollection.ID"] )
}
GSE198449_filtered_metadata_culled_sorted[1:3,"biocollection.ID"]

# See how many biocollection_IDs in the metadata are unqiue
dim(GSE198449_filtered_metadata_culled_sorted) # 2775    6
dim(unique( GSE198449_filtered_metadata_culled_sorted[,"biocollection.ID"] )) # all of them are unique, 2775

# How many column names in the data are unique
dim(GSE198449_data)
length(unique(colnames(GSE198449_data))) # all, but just 1858, so there is a discrepancy in the number of samples

# Now find samples that have data and metadata - matching based on column name for the data and
# biocollection.ID for the metadata
# Also, keep track of the plate variable to add it to the metadata
data_match_list <- list()
metadata_match_list <- list()
plate_list <- list()
for (i in 1:ncol(GSE198449_data)){
  GSE198449_sample <- gsub(pattern="_P.", replacement="", x=colnames(GSE198449_data)[i]) # remove the _P* tag
  GSE198449_split_string <- strsplit(x = colnames(GSE198449_data)[i], split = "_")
  GSE198449_plate <- GSE198449_split_string[[1]][3]
  is_present <- GSE198449_sample %in% GSE198449_filtered_metadata_culled_sorted$biocollection.ID # way to reference tibble column as a list
  if( is_present == TRUE ){
    print(paste(i,is_present))
    data_match_list <- c(data_match_list, colnames(GSE198449_data)[i])
    metadata_match_list <- c(metadata_match_list, GSE198449_sample)
    plate_list <- c(plate_list, GSE198449_plate)
    }
} 
length(data_match_list) # just 507 of them, this is surprising
length(metadata_match_list) # same
length(plate_list) # same

# Cull data and metadata to just those 507 that have data and metadata
# Also, add the plate var to the metadata

# First cull the metadata
GSE198449_metadata <- data.frame()
for (i in 1:length(metadata_match_list)){
  GSE198449_sample <- metadata_match_list[[i]]
  print(GSE198449_sample)
  selected_row <- GSE198449_filtered_metadata_culled_sorted[GSE198449_filtered_metadata_culled_sorted[,"biocollection.ID"] == GSE198449_sample,]
  GSE198449_metadata <- rbind(GSE198449_metadata, selected_row)
}
dim(GSE198449_metadata) # 507   6
# Add the plate var to the metadata
GSE198449_metadata$plate <- as.character(plate_list)

# Explore how plates are distributed
ggplot(
  data=GSE198449_metadata,
  aes(x=PCR.test.for.SARS.Cov.2, fill=plate))+
  geom_bar() # looks pretty random

# Name the rows of metadata with reconstructed full sample names
# i.e get this 
colnames(GSE198449_data)[1] # "X20_5068.T00_P3"
# from these
GSE198449_metadata[1,"biocollection.ID"] # "X20_5068.T00"
GSE198449_metadata[1,"plate"][1] # "P3"
# paste them together
metadata_rownames <- paste(GSE198449_metadata$biocollection.ID, "_", GSE198449_metadata$plate, sep="")
# Now use them and check them
rownames(GSE198449_metadata) <- metadata_rownames # Tibble will complain here

# Now cull the data
GSE198449_data <- GSE198449_data[,as.character(data_match_list)]

# Double check to make sure that the colnames of data and rownames of metadata(i.e. the sample names)
# are the same and...
common_sample_names <- intersect(
  colnames(GSE198449_data),
  rownames(GSE198449_metadata)
)
length(common_sample_names) # 507, They contain the same entries
# and are identically ordered
identical_sample_names <- identical(
  colnames(GSE198449_data),
  rownames(GSE198449_metadata)
)
identical_sample_names # TRUE, They are ordered the same, but may not have been

# So order/sort them to be sure
GSE198449_ordered_sample_names <- sort(colnames(GSE198449_data))

# Make sure that data are the correct sizes -- 
# 507 columns for data
# 507 rows for metadata
dim(GSE198449_data) # 58929   507
dim(GSE198449_metadata) # 507   7
# now order them both
GSE198449_data <- GSE198449_data[ ,GSE198449_ordered_sample_names ]
GSE198449_metadata <- GSE198449_metadata[GSE198449_ordered_sample_names, ] 
# NOTE: subsetting a tibble like this eliminates the rownames
# This is a feature of tibbles, not a bug. Add them back below
rownames(GSE198449_metadata) <- GSE198449_ordered_sample_names # Tibble will complain
dim(GSE198449_data) # 58929   507
dim(GSE198449_metadata) # 507   7
# try this again
identical_sample_names <- identical(
  colnames(GSE198449_data),
  rownames(GSE198449_metadata)
)
identical_sample_names # TRUE, Now it seems ok

# SAVE THE DATA AND THE METADATA TO FILE ----------------------------------

export_data(data_object = GSE198449_data, file_name = "GSE198449.data.txt" )
export_data(data_object = GSE198449_metadata, file_name ="GSE198449.metadata.txt"  )

# re-open data and metadata and make sure that each is sorted by the id
##source("~/Documents/GitHub/workflow_play/import_metadata.r")
GSE198449_data2 <- import_data("GSE198449.data.txt")
GSE198449_metadata2 <- import_metadata("GSE198449.metadata.txt")

identical_sample_names <- identical(
  colnames(GSE198449_data2),
  rownames(GSE198449_metadata2)
)
identical_sample_names # TRUE, Now we can proceed with the real meat of the analyses
rm(GSE198449_data2, GSE198449_metadata2)

# OPTIONAL Clean House ----------------------------------------------------

# At this point you may want to clean house with the dreaded "rm(list = ls())" command.
# Note the if you do, you will have to run the bit of code from above that sources
# functions from github again. Sorry they're not in a package as of yet.

# FUN STUFF ---------------------------------------------------------------

# The analysis really starts at the this point, everything up until now has been
# data wrangling to make the analyses below possible
# It is concerning that a majority (1352 from 1858) of samples were eliminated based on the 
# apparent absence of metadata. I think I must have missed something, Please let me know
# if you find it!

# import the data and metadata  -------------------------------------------
GSE198449_data <- import_data("GSE198449.data.txt")
GSE198449_metadata <- import_metadata("GSE198449.metadata.txt")

# look at distribution of raw data ----------------------------------------
all_GSE198449_data <- as.vector(GSE198449_data)
hist(all_GSE198449_data, breaks =100)
hist(log10(all_GSE198449_data), breaks =100) # is roughly log normal
# with ggplot
all_GSE198449_data <- data.frame(Values = as.vector(GSE198449_data))
ggplot(data=all_GSE198449_data,mapping=aes(x=Values)) +
  geom_histogram(bins = 100, fill = "purple", color = "black") +
  scale_x_log10() +
  labs(title = "Histogram of expression values",
       x = "Log-Transformed Values",
       y = "Frequency")

# look at the distribution of summary values in the raw data
GSE198449_raw_mean <- apply(GSE198449_data, 2, mean)
GSE198449_raw_median <- apply(GSE198449_data, 2, median) 
GSE198449_raw_min <- apply(GSE198449_data, 2, min) 
GSE198449_raw_max <- apply(GSE198449_data, 2, max) 
GSE198449_raw_sd <- apply(GSE198449_data, 2, sd)

plot(GSE198449_raw_mean, pch=20) # several outliers, roughly any above 400
plot(GSE198449_raw_median, pch=20) # no apparent outliers
plot(GSE198449_raw_min, pch=20) # no apparent outliers
plot(GSE198449_raw_max, pch=20) # one extreme outlier, few other minor ones
plot(GSE198449_raw_sd, pch=20) # one extreme outlier, few other minor ones

# make these data tidy for ggplot
GSE198449_raw_descriptive_stats <- as.matrix(names(GSE198449_raw_mean))
GSE198449_raw_descriptive_stats <- cbind(GSE198449_raw_descriptive_stats, as.matrix(GSE198449_raw_mean))
GSE198449_raw_descriptive_stats <- cbind(GSE198449_raw_descriptive_stats, as.matrix(GSE198449_raw_median))
GSE198449_raw_descriptive_stats <- cbind(GSE198449_raw_descriptive_stats, as.matrix(GSE198449_raw_min))
GSE198449_raw_descriptive_stats <- cbind(GSE198449_raw_descriptive_stats, as.matrix(GSE198449_raw_max))
GSE198449_raw_descriptive_stats <- cbind(GSE198449_raw_descriptive_stats, as.matrix(GSE198449_raw_sd))
colnames(GSE198449_raw_descriptive_stats) <- c("GSE198449_sample_ids","GSE198449_mean", "GSE198449_median", "GSE198449_min", "GSE198449_max", "GSE198449_sd")
rownames(GSE198449_raw_descriptive_stats) <- names(GSE198449_raw_mean)
GSE198449_raw_descriptive_stats <- as.data.frame(GSE198449_raw_descriptive_stats)
GSE198449_raw_descriptive_stats$GSE198449_mean <- as.numeric(GSE198449_raw_descriptive_stats$GSE198449_mean)
GSE198449_raw_descriptive_stats$GSE198449_median <- as.numeric(GSE198449_raw_descriptive_stats$GSE198449_median)
GSE198449_raw_descriptive_stats$GSE198449_min <- as.numeric(GSE198449_raw_descriptive_stats$GSE198449_min)
GSE198449_raw_descriptive_stats$GSE198449_max <- as.numeric(GSE198449_raw_descriptive_stats$GSE198449_max)
GSE198449_raw_descriptive_stats$GSE198449_sd <- as.numeric(GSE198449_raw_descriptive_stats$GSE198449_sd)

# Use ggplot to look at all of the descriptive stats at one time
# First try it geom_line
ggplot(
  data=GSE198449_raw_descriptive_stats,
  mapping=aes(x=1:nrow(GSE198449_raw_descriptive_stats))) +
  scale_y_log10() +
  geom_line(aes(y=GSE198449_mean, colour="GSE198449_mean")) +
  geom_line(aes(y=GSE198449_median, colour="GSE198449_median")) +
  geom_line(aes(y=GSE198449_min, colour="GSE198449_min")) +
  geom_line(aes(y=GSE198449_max, colour="GSE198449_max")) +
  geom_line(aes(y=GSE198449_sd, colour="GSE198449_sd")) +
  scale_color_manual(name = "Stat Values", values = c("GSE198449_mean" = "black", "GSE198449_median" = "red", "GSE198449_min"="blue", "GSE198449_max"="green", "GSE198449_sd"="purple"))

# Now try it with geom_point
ggplot(
  data=GSE198449_raw_descriptive_stats,
  mapping=aes(x=1:nrow(GSE198449_raw_descriptive_stats))) +
  scale_y_log10() +
  geom_point(aes(y=GSE198449_mean, colour="GSE198449_mean")) +
  geom_point(aes(y=GSE198449_median, colour="GSE198449_median")) +
  geom_point(aes(y=GSE198449_min, colour="GSE198449_min")) +
  geom_point(aes(y=GSE198449_max, colour="GSE198449_max")) +
  geom_point(aes(y=GSE198449_sd, colour="GSE198449_sd")) +
  scale_color_manual(name = "Stat Values", values = c("GSE198449_mean" = "black", "GSE198449_median" = "red", "GSE198449_min"="blue", "GSE198449_max"="green", "GSE198449_sd"="purple"))
# In both visualizations there appears to be at least one extreme outlier, around index 30

# Find row indices where Value is larger than 4e+06 - this is the extreme outlier
row_indices <- which(GSE198449_raw_descriptive_stats$GSE198449_max > 4e+06) # 37

GSE198449_raw_descriptive_stats[row_indices,]
#                  GSE198449_sample_ids  GSE198449_mean GSE198449_median GSE198449_min   GSE198449_max    GSE198449_sd
#X20_5217.T00_P5 X20_5217.T00_P5 576.3335         0      0 10855021 58296.76

# remove and see how things look
GSE198449_raw_descriptive_stats <- GSE198449_raw_descriptive_stats[-row_indices,]

ggplot(
  data=GSE198449_raw_descriptive_stats,
  mapping=aes(x=1:nrow(GSE198449_raw_descriptive_stats))) +
  scale_y_log10() +
  geom_point(aes(y=GSE198449_mean, colour="GSE198449_mean")) +
  geom_point(aes(y=GSE198449_median, colour="GSE198449_median")) +
  geom_point(aes(y=GSE198449_min, colour="GSE198449_min")) +
  geom_point(aes(y=GSE198449_max, colour="GSE198449_max")) +
  geom_point(aes(y=GSE198449_sd, colour="GSE198449_sd")) +
  scale_color_manual(name = "Stat Values", values = c("GSE198449_mean" = "black", "GSE198449_median" = "red", "GSE198449_min"="blue", "GSE198449_max"="green", "GSE198449_sd"="purple"))
# Things look pretty good once this one outlier sample is removed, so I will remove it 
# from the data and metadata before further analysis

# Remove the outlier
GSE198449_data <- GSE198449_data[, !colnames(GSE198449_data) %in% "X20_5217.T00_P5"]
GSE198449_metadata <- GSE198449_metadata[!rownames(GSE198449_metadata) %in% "X20_5217.T00_P5",]

# save the "edited" data
export_data(data_object = GSE198449_data, file_name = "GSE198449.data.edit.txt")
export_data(data_object = GSE198449_metadata, file_name = "GSE198449.metadata.edit.txt")

# Preprocess and look at distribution of normed data -------------------------------------

# Preprocess the data -- attempt to normalize it
preprocessing_tool("GSE198449.data.edit.txt", pseudo_count=0.1) 
# NOTE: all the viz methods break without pseudo counts

# Look at distribution of normed data
GSE198449_data_n <- import_data("GSE198449.data.edit.txt.quantile.PREPROCESSED.txt")
all_GSE198449_data_n <- as.vector(GSE198449_data_n)
hist(all_GSE198449_data_n, breaks =100)
hist(log10(all_GSE198449_data_n), breaks =100) # Data are approximately, but not quite, normal
# The non-normality of the data will affect our selection of a statistical test below.

# Look at the distribution of summary values in the normed data
GSE198449_norm_mean <- apply(GSE198449_data_n, 2, mean)
GSE198449_norm_median <- apply(GSE198449_data_n, 2, median) 
GSE198449_norm_min <- apply(GSE198449_data_n, 2, min) 
GSE198449_norm_max <- apply(GSE198449_data_n, 2, max) 
GSE198449_norm_sd <- apply(GSE198449_data_n, 2, sd)

plot(GSE198449_norm_mean) # several outliers, roughly any below 0.232
plot(GSE198449_norm_median) # several apparent outliers, those below 0.15
plot(GSE198449_norm_min) # no apparent outliers
plot(GSE198449_norm_max) # no apparent outliers
plot(GSE198449_norm_sd) # a few outliers, those above 0.256

# First blush, identified additional outlier samples but did not 
# eliminate/remove any more. Will retain them to see if they correspond 
# with any metadata in the PCoA and heatmap-dendrograms.

# PCoA --------------------------------------------------------------------

# calculate PCoA on the preprocessed data
calculate_pco(file_in="GSE198449.data.edit.txt.quantile.PREPROCESSED.txt")
# NOTE: This will take a couple minutes

# render interactive 3d PCoA plot from the PCoA and the corresponding metadata
plot_interactive_colored_3d_pcoa(
  pcoa_data_file = "GSE198449.data.edit.txt.quantile.PREPROCESSED.txt.euclidean.PCoA",
  selected_eigen_vectors = c(1,2,3),
  pcoa_metadata_file = "GSE198449.metadata.edit.txt",
  metadata_column = "PCR.test.for.SARS.Cov.2" # "plate" or any other value from colnames GSE198449_metadata will work
)
# There are clearly some outliers that could have been removed.
# More important, while not perfect, segregation between 
# Covid "Not" and "Detected" is strong and pretty clear  

# Iterate through all of the metadata, creating a colored PCoA per column
# Do this to see if there are any obvious correlations beween the data and
# metadata
plot_static_colored_3d_pcoas(
  pcoa_filename = "GSE198449.data.edit.txt.quantile.PREPROCESSED.txt.euclidean.PCoA",
  selected_eigen_vectors = c(1,2,3),
  metadata_filename = "GSE198449.metadata.edit.txt",
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
# There is some pattern with respect to "plate", interpret as batch effect
# and also to "sample.collection.time point.days.since.TO.", interpret the latter 
# as a biological signal. "biocollection.ID", "enrollment.batch", "participant.ID", and
# "sympton" all look pretty random; I think this is as we would expect.
# For now, we're just interested in "PCR.test.for.SARS.Cov.2". It contains the Covid state
# related information that we are interested in. We examine it further below.

# perform a stat test on the data -----------------------------------------

# perform a stat test on the data/metadata (selected metadata column is used to create groups for chosen stat test)
# We'll use Mann-Whitney as the data do not appear to be normally distributed, even after our preprocessing procedure.
# We're testing for significance of differential expression between "Not" and "Detected" under "PCR.test.for.SARS.Cov.2"
sigtest(data_file="GSE198449.data.edit.txt.quantile.PREPROCESSED.txt", 
        metadata_file="GSE198449.metadata.edit.txt",  
        metadata_column="PCR.test.for.SARS.Cov.2", 
        stat_test="Mann-Whitney-unpaired-Wilcoxon", # two groups of non-normally distributed data
        p_adjust_method = "BH"
)

# Load the stat test results and subselect data based on the stat results
GSE198449_stat_results <- import_data("GSE198449.data.edit.txt.quantile.PREPROCESSED.Mann-Whitney-unpaired-Wilcoxon.PCR.test.for.SARS.Cov.2.STAT_RESULTS.txt")
GSE198449_stat_results <- as_tibble(GSE198449_stat_results,rownames=NA) # Convert to a tibble but KEEP the rownames

# Take a look at the stat results
ggplot(
  data=GSE198449_stat_results,
  mapping=aes(x=1:length(bonferroni_p))) +
  geom_line(aes(y=p, colour="p")) +
  geom_line(aes(y=bonferroni_p, colour="bonferroni_p")) +
  geom_line(aes(y=BH_p, colour="BH_p")) +
  scale_color_manual(name = "Stat Values", values = c("p" = "black", "bonferroni_p" = "red", "BH_p"="blue"))
# In this case, Bonferroni (0.05) does not appear to be an overly aggressive filter
# More than 5,000 genes pass it -- filter them a bit below based on more arbitrary criterea

# Filter the data to retain just the top 5% "most" significant 
dim(GSE198449_stat_results) # 40237   513
GSE198449_stat_results_subselected <- GSE198449_stat_results |> 
  rownames_to_column(var="GSE198449_rowname") |> # do this to retain the rownames somewhere in the tibble/matrix
  filter( bonferroni_p < 0.00000000001 ) |> # Arbitrary, selected to retain ~5% most significant as calculated below
  column_to_rownames(var="GSE198449_rowname") # put the rownames back where they were; this is a feature not a bug
dim(GSE198449_stat_results_subselected) # 2007  513
dim(GSE198449_stat_results_subselected)[1]/dim(GSE198449_stat_results)[1]*100 # find a value that retains 5% most significant

# Export results with stats
export_data(data_object = GSE198449_stat_results_subselected, file_name = "GSE198449_stat_results_subselected.txt")
# remove the columns that contain the stats(this is hacky)
GSE198449_stat_results_subselected_and_cleaned <- remove_last_n_columns(GSE198449_stat_results_subselected, n=7)
# make sure that samples are sorted/ordered as they are in the metadata above 
length(GSE198449_ordered_sample_names) # 507

# Export the data ready to create a new PCoA or HD
export_data(data_object = GSE198449_stat_results_subselected_and_cleaned, file_name = "GSE198449_stat_results_subselected_and_cleaned.txt")

# Calculate and render a new PCoA just on the stat filtered values
calculate_pco(file_in="GSE198449_stat_results_subselected_and_cleaned.txt")

# render interactive 3d PCoA plot from the PCoA and the corresponding metadata
plot_interactive_colored_3d_pcoa(
  pcoa_data_file = "GSE198449_stat_results_subselected_and_cleaned.txt.euclidean.PCoA",
  selected_eigen_vectors = c(1,2,3),
  pcoa_metadata_file = "GSE198449.metadata.edit.txt",
  metadata_column = "PCR.test.for.SARS.Cov.2"
)
# Outliers still apparent but Covid dependent segregation is very clear

# heatmap dendrograms -----------------------------------------------------

# create heatmap dendrograms of original and stat subselected data
# First, the original data
heatmap_dendrogram(file_in = "GSE198449.data.edit.txt.quantile.PREPROCESSED.txt",
                   metadata_table = "GSE198449.metadata.edit.txt",
                   metadata_column="PCR.test.for.SARS.Cov.2"
) # This operation may fail, the image is complex and requires a lot of resources
# It also takes a considerable amount of time, almost an hour on my laptop.
shell("start GSE198449.data.edit.txt.quantile.PREPROCESSED.txt.HD.png")

# Now just the statistically subselected data - this is much quicker.
heatmap_dendrogram(file_in = "GSE198449_stat_results_subselected_and_cleaned.txt",
                   metadata_table = "GSE198449.metadata.edit.txt",
                   metadata_column="PCR.test.for.SARS.Cov.2"
)
shell("start GSE198449_stat_results_subselected_and_cleaned.txt.HD.png")

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
# gene name and gene description for each ENSG id 
annotations <- read.table(file="mart_export.txt",row.names=NULL,header=TRUE,sep="\t", # compared to import_data, changed row.names from 1 to NULL
                          colClasses = "character", check.names=FALSE,
                          comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE)

GSE198449_stat_results_subselected <- import_data("GSE198449_stat_results_subselected.txt")
# The ENSG annotations in the data set have the version tag.
# Remove the version tag and you can match to the Gene stable ID
GSE198449_genes <- rownames(GSE198449_stat_results_subselected)
GSE198449_genes[1]
GSE198449_genes <- gsub(pattern="\\..*", replacement = "", x=GSE198449_genes) # have a .num tag, part of the stable ID version
GSE198449_genes[1]

# add columns for the annotation information, NA to start
GSE198449_stat_results_subselected <- as.data.frame(GSE198449_stat_results_subselected)
GSE198449_stat_results_subselected$gene_descriptions <- NA 
GSE198449_stat_results_subselected$gene_names <- NA

# now use a simple loop to look up and add the annotation values to the stat results
for ( i in 1:length(GSE198449_genes)){ 
  GSE198449_gene_annotation <- annotations[annotations$`Gene stable ID` == GSE198449_genes[i],]
  GSE198449_stat_results_subselected[i,"gene_names"] <- GSE198449_gene_annotation[1,"Gene name"]
  GSE198449_stat_results_subselected[i,"gene_descriptions"] <- GSE198449_gene_annotation[1,"Gene description"]
}
export_data(data_object = GSE198449_stat_results_subselected, file_name = "GSE198449.annotations.txt")

# PRELIMINARY PATHWAY ANALYSIS --------------------------------------------

# Use the list of gene ids from above
GSE198449_genes

# Perform pathway analysis and generate an interactive visualization
gostres <- gost(query = GSE198449_genes, organism = 'hsapiens', significant = FALSE)
gostplot(gostres, capped = TRUE, interactive = TRUE)

# That's it for the analysis of GSE212041.
# Continue on to one of the other three datasets and combine_covid when you're done with all three

# JUNK BELOW HERE ---------------------------------------------------------
# JUNK BELOW HERE ---------------------------------------------------------
# JUNK BELOW HERE ---------------------------------------------------------

# Get the index of the column
column_index <- which(colnames(GSE198449_data) == column_name)

# remove column or row by index
GSE198449_data <- GSE198449_data[, !colnames(GSE198449_data) %in% "X20_5217.T00_P5"]
GSE198449_metadata <- GSE198449_metadata[!rownames(GSE198449_metadata) %in% "X20_5217.T00_P5",]

# rownames with Tibble
rownames_to_column(x,var="GSE198449_rowname") 
column_to_rownames(x,var="GSE198449_rowname")

# https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&search=covid&display=20&zsort=samples
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198449
# Submission date	Mar 11, 2022
# Last update date	Jul 14, 2023
# Contact name	Yongchao Ge
# E-mail(s)	yongchao.ge@gmail.com
# Organization name	Icahn School of Medicine at Mount Sinai
# Department	Neurology
# Street address	One Gustave L Levy Place #1137
# City	New York
# State/province	NY
# ZIP/Postal code	10029
# Country	USA