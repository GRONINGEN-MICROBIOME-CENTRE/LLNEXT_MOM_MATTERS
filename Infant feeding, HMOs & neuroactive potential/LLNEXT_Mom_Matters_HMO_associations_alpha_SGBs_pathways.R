################## Associations HMOs with bacterial taxa and metabolic pathways ####### 
# Author: Trishla Sinha 
# Last update: 7th October, 2025

# Load packages 
library(tidyverse)
library(wesanderson)

# Load all functions from script : 
source("Functions_associations_phenotypes_LLNEXT_new.R")

## Associations with HMOs in breast milk in exclusively breastfed infants  ##
hmo_all <-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/hmo/HMOs_matched_infant_MGS_breastfeeding_only_early_timepoints.txt")
row.names(hmo_all)<-hmo_all$NG_ID

hmo_colnames <- names(hmo_all)[grepl("^mother_milk_HMO", names(hmo_all))]
# Tested HMO's 
# "mother_milk_HMO_milk_group"     "mother_milk_HMO_Le"             "mother_milk_HMO_Se"            
# "mother_milk_HMO_2FL_ugml"       "mother_milk_HMO_3FL_ugml"       "mother_milk_HMO_LDFT_ugml"     
# "mother_milk_HMO_3GL_ugml"       "mother_milk_HMO_6GL_ugml"       "mother_milk_HMO_A_tetra_ugml"  
# "mother_milk_HMO_3SL_ugml"       "mother_milk_HMO_LNT_ugml"       "mother_milk_HMO_LNnT_ugml"     
# "mother_milk_HMO_6SL_ugml"       "mother_milk_HMO_3F3SL_ugml"     "mother_milk_HMO_LNFP_V_ugml"   
# "mother_milk_HMO_LNnFP_V_ugml"   "mother_milk_HMO_LNFP_I_ugml"    "mother_milk_HMO_LNFP_III_ugml" 
# "mother_milk_HMO_LNFP_II_ugml"   "mother_milk_HMO_LNnDFH_ugml"    "mother_milk_HMO_LSTb_ugml"     
# "mother_milk_HMO_LNDH_I_ugml"    "mother_milk_HMO_LSTc_ugml"      "mother_milk_HMO_LNH_ugml"      
# "mother_milk_HMO_MFLNH_III_ugml" "mother_milk_HMO_DSLNT_ugml"     "mother_milk_HMO_DFLNHa_ugml"   
# "mother_milk_HMO_Total_ugml"     "mother_milk_HMO_Fuc_ugml"       "mother_milk_HMO_Neut_ugml"     
# "mother_milk_HMO_Sia_ugml"    

# Inverse transformation 
hmo_all[hmo_colnames] <- run_invrank_dataFrame(hmo_all[hmo_colnames])
num_unique_infants <- length(unique(hmo_all$NEXT_ID))

# With infant alpha diversity 
alpha<-hmo_all %>% select(shannon)

# With Infant SGB's 
taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
taxa <- taxa[match(rownames(hmo_all), rownames(taxa)),]

# With pathways 
pathways<-read.delim("~/Desktop/LLNEXT/Analysis/pathways/metacyc_infants_ALR_zeroTreated_20250722.tsv")
pathways <- pathways[match(rownames(hmo_all), rownames(pathways)),]


# Cazyme analysis 
metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                      "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint")

result_HMO_alpha <- gam_function(alpha, hmo_all,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER"))
result_HMO_taxa <- gam_function(taxa, hmo_all,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER"))
result_HMO_pathway <- gam_function(pathways, hmo_all,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER"))


result_HMO_alpha$FDR<-p.adjust (result_HMO_alpha$p, method = "fdr")
result_HMO_taxa$FDR<-p.adjust (result_HMO_taxa$p, method = "fdr")
result_HMO_pathway$FDR<-p.adjust (result_HMO_pathway$p, method = "fdr")
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/hmo")

write.table(result_HMO_alpha, "result_HMO_alpha_oct_2025.txt", sep = "\t", row.names = F)
write.table(result_HMO_taxa, "result_HMO_taxa_SGB_oct_2025.txt", sep = "\t", row.names = F)
write.table(result_HMO_pathway, "result_HMO_pathway_oct_2025.txt", sep = "\t", row.names = F)


###############################
# Figure 1: Shannon diversity
###############################

plot_shannon <- ggplot(
  hmo_all,
  aes(
    mother_milk_HMO_milk_group,
    y = shannon,
    fill = mother_milk_HMO_milk_group,
    color = mother_milk_HMO_milk_group
  )
) +
  scale_fill_manual(values = wes_palette("BottleRocket2", n = 4)) +
  scale_color_manual(values = wes_palette("BottleRocket2", n = 4)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA) +
  geom_point(alpha = 0.6,
             position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0)) +
  theme_bw() +
  labs(
    x = "",
    y = "Shannon Diversity Index",
    fill = "Mother Milk HMO Milk Group",
    color = "Mother Milk HMO Milk Group"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 60, hjust = 1),
    legend.position = "none"
  )

ggsave(
  "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/hmo/shannon_diversity_milk_group.pdf",
  plot = plot_shannon,
  width = 4,
  height = 4,
  useDingbats = FALSE
)

source_data_shannon <- hmo_all %>%
  dplyr::select(
    mother_milk_HMO_milk_group,
    shannon
  )

write.csv(
  source_data_shannon,
  "/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_post_acceptance_2026/Source_data/Source_Data_HMO_group_vs_Shannon.csv",
  row.names = FALSE
)


###############################
# Figure 2: Clostridium perfringens
###############################

all <- merge(hmo_all, taxa, by = "row.names")

plot_clostridium <- ggplot(
  all,
  aes(
    mother_milk_HMO_milk_group,
    y = Clostridium_perfringens.t__SGB6191,
    fill = mother_milk_HMO_milk_group,
    color = mother_milk_HMO_milk_group
  )
) +
  scale_fill_manual(values = wes_palette("BottleRocket2", n = 4)) +
  scale_color_manual(values = wes_palette("BottleRocket2", n = 4)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA) +
  geom_point(alpha = 0.6,
             position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0)) +
  theme_bw() +
  labs(
    x = "",
    y = "Clostridium perfringens SGB6191",
    fill = "Mother Milk HMO Milk Group",
    color = "Mother Milk HMO Milk Group"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 60, hjust = 1),
    legend.position = "none"
  )

ggsave(
  "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/hmo/HMO_group_vs_Clostridium_perfringens.pdf",
  plot = plot_clostridium,
  width = 4,
  height = 4,
  useDingbats = FALSE
)

source_data_clostridium <- all %>%
  dplyr::select(
    mother_milk_HMO_milk_group,
    Clostridium_perfringens.t__SGB6191
  )

write.csv(
  source_data_clostridium,
  "/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_post_acceptance_2026/Source_data/Source_Data_HMO_group_vs_Clostridium_perfringens.csv",
  row.names = FALSE
)
