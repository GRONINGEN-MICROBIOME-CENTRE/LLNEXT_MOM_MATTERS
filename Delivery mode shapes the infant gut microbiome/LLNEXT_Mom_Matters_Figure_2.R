####################### Figure 2 ########################

# Author: Trishla Sinha
# This script gives the code for Main figure 2 and related supplementary tables 
# Last update: 2nd July, 2026, adding annotations to mention links to supplementary tables and source data

# Load libraries 
library (ggplot2)
library(tidyverse)


# Figure 2e: alpha diversity associations with cross-sectional and dynamic traits 

# Load files
alpha_dir <- "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/alpha_diversity"

alpha_dynamic <- read.delim(
  file.path(alpha_dir, "Alpha_diversity_dynamic_phenotypes_results_all_23_09_2025.txt")  # Supplementary table S24
)

alpha_cross <- read.delim(
  file.path(alpha_dir, "Alpha_diversity_cross_phenotypes_results_23_09.txt") # Supplementary table S23
)

# Dynamic phenotypes: calculate t-value and rename effect.level to levels
sig_alpha_dynamic <- alpha_dynamic %>%
  filter(!is.na(FDR_base), FDR_base < 0.05) %>%
  mutate(
    levels = effect.level,
    t.value_cor_base = Estimate_cor_base / SE_cor_base
  ) %>%
  select(
    outcome,
    trait,
    levels,
    Estimate_cor_base,
    SE_cor_base,
    t.value_cor_base,
    p_cor_base,
    FDR_base,
    p_cor_delivery_breastfed,
    trait_group,
    type_phenotype
  )

# Cross-sectional phenotypes
sig_alpha_cross <- alpha_cross %>%
  filter(!is.na(FDR_base), FDR_base < 0.05) %>%
  select(
    outcome,
    trait,
    levels,
    Estimate_cor_base,
    SE_cor_base,
    t.value_cor_base,
    p_cor_base,
    FDR_base,
    p_cor_delivery_breastfed,
    trait_group,
    type_phenotype
  )

# Merge dynamic and cross-sectional significant traits
merged_alpha_sig <- bind_rows(sig_alpha_dynamic, sig_alpha_cross) %>%
  mutate(
    variable = if_else(
      is.na(levels) | levels == "",
      trait,
      paste0(trait, "_", levels)
    )
  )

# Source data 2e 
merged_alpha_sig %>%
  select(type_phenotype, trait, levels, variable, t.value_cor_base, FDR_base)

write.table(merged_alpha_sig, "/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_post_acceptance_2026/Source_data/source_data_fig_2_e_merged_alpha_sig.txt", sep="\t", row.names = F)


# Plot 2e 
merged_alpha <- ggplot(
  merged_alpha_sig,
  aes(
    x = reorder(variable, t.value_cor_base),
    y = t.value_cor_base,
    fill = t.value_cor_base > 0
  )
) +
  geom_bar(stat = "identity", color = "black") +
  xlab("") +
  ylab("T-value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.position = "top",
    legend.justification = c(1.2, 0.1)
  ) +
  scale_fill_manual(
    values = c("red", "blue"),
    labels = c("Decreased", "Increased")
  ) +
  coord_flip() +
  guides(fill = guide_legend(title = ""))

merged_alpha

ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_2/shannon_diversity_infants_all.pdf", 
       width = 6, height = 5)

# Figure 2f
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint")
ResultsAdonis<-read.delim("combined_results_timepoint_adonis.txt") # Supplementary table S25
ReSultsTCAM<-read.delim("combined_results_TCAM_adonis.txt") # Supplementary table S26
ReSultsTCAM$Timepoint <-"Overall"
merged_adonis <- bind_rows(ResultsAdonis, ReSultsTCAM)
significant<-merged_adonis[merged_adonis$p.value_cor_delivery_breastfed<0.001,]
significant$Timepoint <- factor(significant$Timepoint, levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12", "Overall"))
significant <- significant %>%
  mutate(R2 = as.numeric(as.character(R2_cor_delivery_breastfed)),  
         R2 = round(R2, 3))  

significant2 <- significant %>% filter(
  !grepl("vitaminK|after_birth|mode_simple",Phenotype))

significant2$Timepoint <- factor(as.character(significant2$Timepoint),
                                 levels = levels(significant2$Timepoint),
                                 ordered = T)
segments <- significant2 %>% group_by(Phenotype) %>% 
  summarise(Timepoint2=min(Timepoint),
            xend=max(Timepoint)) %>% rename(Timepoint=Timepoint2)


segments <- rbind(segments,data.frame(Phenotype="birth_deliverybirthcard_mode_binary",
                                      Timepoint="M6",xend="M6"))
segments <- rbind(segments,data.frame(Phenotype="birth_deliverybirthcard_mode_binary",
                                      Timepoint="M9",xend="M9"))


adonis_infant_per_timepoint<-ggplot()+
  geom_segment(data=segments, aes(x=Timepoint, y=Phenotype, 
                                  xend=xend,yend=Phenotype),color="black",linewidth=0.2)+
  geom_point(data = significant2,
             aes(x = Timepoint, y=Phenotype,
                 size = R2, fill = Timepoint),shape=21,color="black") +
  # facet_wrap(~ Timepoint) + 
  # coord_flip() +  
  xlab("Timepoint") +  
  ylab("") +  
  theme_bw() +
  theme(axis.text.y=element_text(size=8),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())+
  scale_size_continuous(range = c(1,4))+
  scale_fill_manual(values = c( "#b5dd88","#41c0b4", "#4397bb", "#eca4c9","#a42097", "#390962", "black"))
adonis_infant_per_timepoint

ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_2_new/infant_adonis.pdf", 
       width = 8, height = 4)

# Source data 2f
write.table(significant2, "/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_post_acceptance_2026/Source_data/source_data_fig_2_f_beta_diversity.txt", sep="\t", row.names = F)

