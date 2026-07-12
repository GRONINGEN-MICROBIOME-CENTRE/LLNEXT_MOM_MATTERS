heatmap<-readRDS("/Users/trishlasinha/Downloads/Source_fig_clusters/Heatmap_source.rds")


# Metadata (first data frame)
metadata <- heatmap[[2]]
row.names(metadata)<-metadata$NG_ID
metadata$NG_ID=NULL

# Heatmap values (second data frame or matrix)
heatmap_matrix <- as.data.frame(heatmap[[1]])

source_data_ED_2a<-merge(metadata, heatmap_matrix,by="row.names")
source_data_ED_2a$Row.names=NULL

# Save source data 
write.csv(
  source_data_ED_2a,
  "/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_post_acceptance_2026/Source_data/source_data_ED_2a.csv",
  row.names = FALSE
)
