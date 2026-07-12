############################### Extended Data 1 ############################
# Author: Trishla Sinha
# Code for Extended Data 1 and related source data
# Last update: 9th July, 2026

library(tidyverse)
library(pheatmap)
library(grid)
library(gridExtra)

base_dir <- "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/timepoint_SGB_mother_infant/"
out_dir  <- "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/supplementary/"

setwd(base_dir)

###############################
# MOTHERS
###############################

mother_data <- read.delim("timepoint_mother_SGBs_2026.txt")

mother_sig <- mother_data %>%
  filter(FDR < 0.05) %>%
  mutate(
    Timepoint = gsub("Timepoint_categorical", "", Feature),
    Source_data = "Extended Data Fig. 1a",
    Sample_type = "Mother"
  )

source_data_mother <- mother_sig %>%
  select(Source_data, Sample_type, Bug, Timepoint, t.value, FDR)

write.csv(
  source_data_mother,
  file.path(out_dir, "Source_Data_Extended_Data_Fig_1a_mothers.csv"),
  row.names = FALSE
)

mother_wide <- mother_sig %>%
  select(Bug, Timepoint, t.value) %>%
  pivot_wider(names_from = Timepoint, values_from = t.value) %>%
  select(Bug, P28, B, M3)

mother_matrix <- mother_wide %>%
  column_to_rownames("Bug") %>%
  as.matrix()

time_mother_SGB <- pheatmap(
  mother_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  display_numbers = FALSE,
  angle_col = 90,
  fontsize_col = 8,
  fontsize_row = 6,
  cellwidth = 5,
  silent = TRUE,
  main = "Mothers"
)

time_mother_SGB

###############################
# INFANTS
###############################

infant_data <- read.delim("timepoint_infant_SGB's.txt")

infant_sig <- infant_data %>%
  filter(FDR < 0.05) %>%
  mutate(
    Timepoint = gsub("Timepoint_categorical", "", Feature),
    Source_data = "Extended Data Fig. 1b",
    Sample_type = "Infant"
  )

source_data_infant <- infant_sig %>%
  select(Source_data, Sample_type, Bug, Timepoint, t.value, FDR)

write.csv(
  source_data_infant,
  file.path(out_dir, "Source_Data_Extended_Data_Fig_1b_infants.csv"),
  row.names = FALSE
)

infant_wide <- infant_sig %>%
  select(Bug, Timepoint, t.value) %>%
  pivot_wider(names_from = Timepoint, values_from = t.value) %>%
  select(Bug, M1, M2, M3, M6, M9, M12)

infant_matrix <- infant_wide %>%
  column_to_rownames("Bug") %>%
  as.matrix()

time_infant_SGB <- pheatmap(
  infant_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  display_numbers = FALSE,
  angle_col = 90,
  fontsize_col = 8,
  fontsize_row = 6,
  cellwidth = 6,
  silent = TRUE,
  main = "Infants"
)


###############################
# SAVE PDFs
###############################

# Mother heatmap
pdf(
  file.path(out_dir, "Extended_Data_Fig_1a_mothers.pdf"),
  width = 95 / 25.4,
  height = 210 / 25.4
)

grid::grid.newpage()
grid::grid.draw(time_mother_SGB$gtable)

dev.off()


# Infant heatmap
pdf(
  file.path(out_dir, "Extended_Data_Fig_1b_infants.pdf"),
  width = 100 / 25.4,
  height = 210 / 25.4
)

grid::grid.newpage()
grid::grid.draw(time_infant_SGB$gtable)

dev.off()
