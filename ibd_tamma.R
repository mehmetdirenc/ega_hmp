library(readr)
library(dplyr)
# library(fitdistrplus)
library(performance)
library(correlation)
choose_correlation <- function(x, y)
{
  # Check normality of both variables
  normal_x <- check_normality(x)
  normal_y <- check_normality(y)
  # Use Pearson if both are normal, otherwise use Spearman
  if (normal_x > 0.05 & normal_y > 0.05)
  {
    method <- "pearson"
  }
   else
  {
    method <- "spearman"
  }
  return(method)
}

metadata_tmp <- readr::read_csv("/home/direnc/tools/ibd_tamma/all_metadata.csv")
bacterial_counts <- read.table("/home/direnc/tools/ibd_tamma/data/bacteria_species/counts/Clostridium_perfringens.tsv", sep = "\t", header = TRUE)
S100A8 <- read.table("/home/direnc/tools/ibd_tamma/data/human/counts/S100A8.tsv", sep = "\t", header = TRUE)
S100A9 <- read.table("/home/direnc/tools/ibd_tamma/data/human/counts/S100A9.tsv", sep = "\t", header = TRUE)
metadata <- metadata_tmp %>%
  dplyr::select(group, tissue, sample)

unique_tissues <- unique(metadata$tissue)
tissue_list <- as.list(unique_tissues)
unique_groups <- unique(metadata$group)
groups_list <- as.list(unique_groups)
colnames(bacterial_counts)[colnames(bacterial_counts) == "counts"] <- "bacterial_counts"
colnames(S100A8)[colnames(S100A8) == "counts"] <- "S100A8_counts"
colnames(S100A9)[colnames(S100A9) == "counts"] <- "S100A9_counts"

bacterial_counts <- bacterial_counts[, !colnames(bacterial_counts) %in% "Feature"]
S100A8 <- S100A8[, !colnames(S100A8) %in% "Gene"]
S100A9 <- S100A9[, !colnames(S100A9) %in% "Gene"]

combined_data <- merge(S100A8, S100A9, by = "sample")
combined_data <- merge(combined_data, bacterial_counts, by = "sample")
combined_data <- merge(combined_data, metadata, by = "sample")

s18 <- "S100A8_counts"
s19 <- "S100A9_counts"
bc  <- "bacterial_counts"
combined_data <- combined_data[!is.na(combined_data[[bc]]), ]

for (it in seq_along(tissue_list))
  {
    selected_tissue <- tissue_list[[it]]
    for (ig in seq_along(groups_list))
    {
      selected_group <- groups_list[[ig]]
      filtered_df <- combined_data %>% dplyr::filter(tissue == selected_tissue & group == selected_group)
      cor_method <- choose_correlation(filtered_df[[bc]], filtered_df[[s18]])
      if (cor_method == "spearman")
      {
        cor_val <- correlation::correlation(filtered_df[[bc]], filtered_df[[s18]], method = cor_method)
      }
      else
      {
        cor_val <- correlation::correlation(filtered_df[[bc]], filtered_df[[s18]], method = cor_method)
      }
      print(paste("Group:", selected_group, ", Tissue:", selected_tissue, ", Correlation method:", cor_method, ", Correlation:", cor_val))
    }
  }


























