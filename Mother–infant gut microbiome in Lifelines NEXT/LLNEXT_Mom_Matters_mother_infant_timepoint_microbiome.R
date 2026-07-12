#############################################################
# DYNAMICS OF MICROBIOME : Timepoint Microbiome Associations
#############################################################

# Author: Trishla Sinha 
# Date last update: 9th July, 2026
# Adding annotations to mention links to supplementary tables and source data


###############################
#Functions
##############################
source("/Users/trishlasinha/Desktop/LLNEXT/Analysis/scripts/new_scripts_2026/Functions_associations_phenotypes_LLNEXT_new.R")

##############################
# Loading libraries
##############################
library(gplots)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(MetBrewer)
library(RLRsim)
library(lmerTest)
library(ggforce)
library(pheatmap)

##############################
# MOTHERS 
##############################

metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
row.names(metadata)<-metadata$NG_ID
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_mothers<-metadata[metadata$Type=="mother", ]
metadata_mothers <- subset(metadata_mothers, Timepoint_categorical != "M1" & Timepoint_categorical != "M2")
metadata_mothers$Timepoint_categorical<-factor(metadata_mothers$Timepoint_categorical, levels = c("P12", "P28", "B", "M3"))

taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/NEXT_metaphlan_4_CLR_transformed_fil_30_percent_SGB_mothers_03_08_2023.txt")


time_on_microbiome_mother <- mixed_models_without_time_correction(metadata_mothers, 
                                                 "NG_ID", 
                                                taxa, 
                                                 c("Timepoint_categorical"))

write.table(time_on_microbiome_mother, "~/Desktop/LLNEXT/Analysis/results/timepoint_SGB_mother_infant/timepoint_mother_SGBs_2026.txt", sep = "\t", row.names = F) # Supplementary Table S17


#################################
# INFANTS
####################################
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
row.names(metadata)<-metadata$NG_ID
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))

taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")


time_on_microbiome_infant <- mixed_models_without_time_correction(metadata_infants, 
                                        "NG_ID", 
                                        taxa, 
                                        c("Timepoint_categorical"))

write.table(time_on_microbiome_infant, "~/Desktop/LLNEXT/Analysis/results/timepoint_SGB_mother_infant/timepoint_infant_SGB's.txt", sep = "\t", row.names = F) # Supplementary Table S19


