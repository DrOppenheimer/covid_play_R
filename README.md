# covid_play_R
A repository for R code related to my re-analysis of Covid-19 RNASeq data

This repository contains R scripts (simple workflow files) that are sufficient to reproduce my re-analysis of three Covid related RNASeq count data sets.

You will need to complete the individual analysis scripts first:

GSE198449.R, GSE215865.R, acd ~/DOwnd GSE212041.R 

before you can successfully complete the final analysis script that combines the results into a single, statistically robust list of genes:

combine_covid.R

I would recommend using RStudio to run through the analyses step by step.
A few notes:
- Not every step can be completed in R, those steps that cannot have explicit instructions, so please read the comments carefully.
- I completed this work on a Mac so commands like system() will not work on a PC
- All users will likely have to edit the paths in the scripts; I’ve tried to keep this as easy as possible with a “dir_path” argument at the top of the three GSE scripts as the only path that needs changing. The combine_covid.R script necessarily has multiple paths that are required to combine the data. I can’t think of a way to simplify this script any more – let me know if you think of something.

This analysis produces many data products and visualizations that were not included in the FORTHCOMING STUDY HERE. 
