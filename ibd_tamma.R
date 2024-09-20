library(readr)
library(dplyr)
library(ggplot2)
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

res_folder <- "/home/direnc/results/ega/external_datasets/"
metadata_tmp <- readr::read_csv("/home/direnc/tools/ibd_tamma/all_metadata.csv")
bacterial_counts <- read.table("/home/direnc/tools/ibd_tamma/data/bacteria_species/counts/Clostridium_perfringens.tsv", sep = "\t", header = TRUE)
S100A8 <- read.table("/home/direnc/tools/ibd_tamma/data/human/counts/S100A8.tsv", sep = "\t", header = TRUE)
S100A9 <- read.table("/home/direnc/tools/ibd_tamma/data/human/counts/S100A9.tsv", sep = "\t", header = TRUE)
TNF <- read.table("/home/direnc/tools/ibd_tamma/data/human/counts/TNF.tsv", sep = "\t", header = TRUE)
IL6 <- read.table("/home/direnc/tools/ibd_tamma/data/human/counts/IL6.tsv", sep = "\t", header = TRUE)
 # <- read.table("/home/direnc/tools/ibd_tamma/data/human/counts/TNF.tsv", sep = "\t", header = TRUE)
metadata <- metadata_tmp %>%
  dplyr::select(group, tissue, sample)

unique_tissues <- unique(metadata$tissue)
# tissue_list <- as.list(unique_tissues)
tissue_list <- c("Colon", "Ileum", "Rectum")
unique_groups <- unique(metadata$group)
groups_list <- as.list(unique_groups)
colnames(bacterial_counts)[colnames(bacterial_counts) == "counts"] <- "CP_counts"
colnames(S100A8)[colnames(S100A8) == "counts"] <- "S100A8"
colnames(S100A9)[colnames(S100A9) == "counts"] <- "S100A9"
colnames(TNF)[colnames(TNF) == "counts"] <- "TNF"
colnames(IL6)[colnames(IL6) == "counts"] <- "IL6"

bacterial_counts <- bacterial_counts[, !colnames(bacterial_counts) %in% "Feature"]
S100A8 <- S100A8[, !colnames(S100A8) %in% "Gene"]
S100A9 <- S100A9[, !colnames(S100A9) %in% "Gene"]
TNF <- TNF[, !colnames(TNF) %in% "Gene"]
IL6 <- IL6[, !colnames(IL6) %in% "Gene"]

combined_data <- merge(S100A8, S100A9, by = "sample")
combined_data <- merge(combined_data, bacterial_counts, by = "sample")
combined_data <- merge(combined_data, TNF, by = "sample")
combined_data <- merge(combined_data, IL6, by = "sample")
combined_data <- merge(combined_data, metadata, by = "sample")

genes <- c("S100A8", "S100A9", "TNF", "IL6")
bc  <- "CP_counts"
combined_data <- combined_data[!is.na(combined_data[[bc]]), ]

cor_results <- data.frame(
  Gene = character(),
  Group = character(),
  Tissue = character(),
  Correlation_Method = character(),
  Correlation_Value = numeric(),
  Correlation_p_value = numeric(),
  stringsAsFactors = FALSE
)

for (gene in genes)
{
  for (selected_tissue in tissue_list)
  {
    # selected_tissue <- tissue_list[[it]]
    for (selected_group in groups_list)
    {
      # selected_group <- groups_list[[ig]]
      filtered_df <- combined_data %>% dplyr::filter(tissue == selected_tissue & group == selected_group)
      cor_method <- tryCatch(
      {
        choose_correlation(filtered_df[[bc]], filtered_df[["TNF"]])
      },
      error = function(cond)
      {
        # message(paste("Group:", selected_group, ", Tissue:", selected_tissue))
        ""
      })
      if (cor_method != "")
      {
        cor_test <- cor.test(filtered_df[[bc]], filtered_df[[gene]], method = cor_method)
        print(paste("Group:", selected_group, "Gene:", gene, ", Tissue:", selected_tissue, ", Correlation method:", cor_method, ", Correlation:", cor_test$estimate))
      }
      else
      {
        next
        # cor_test <- cor.test(filtered_df[[bc]], filtered_df[[gene]], method = "spearman")
      }
      cor_val <- cor_test$estimate
      cor_p <- cor_test$p.value
      cor_results <- rbind(cor_results, data.frame(
        Gene = gene,
        Group = selected_group,
        Tissue = selected_tissue,
        Correlation_Method = ifelse(cor_method != "", cor_method, "pearson"),
        Correlation_Value = cor_val,
        Correlation_p_value = cor_p,
        stringsAsFactors = FALSE
      ))
    }
  }
}
# filtered_data <- combined_data %>%
#   filter(tissue == "Colon" & group == "CD" & !is.na(CP_counts) & !is.na(TNF))
#
# # Generate the dot plot for TNF vs CP_counts for the filtered data
# ggplot(filtered_data, aes(x = CP_counts, y = TNF)) +
#   geom_point(size = 3, color = "#0072B2") +  # Customize dot size and color
#   theme_minimal() +
#   labs(title = "Dot Plot of TNF vs CP_counts (Colon, CD)",
#        x = "CP_counts",
#        y = "TNF Expression") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Generate dot plot for each gene
for (gene in genes) {
  # Filter data for the current gene
  plot_data_tmp <- cor_results %>% dplyr::filter(Gene == gene)
  plot_data <- plot_data_tmp %>% dplyr::filter(Correlation_p_value <= 0.05)

  # Create the dot plot with custom y-axis ticks and colors
  p <- ggplot(plot_data, aes(x = Group, y = Correlation_Value, color = Tissue)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = paste("Correlation for", gene),
         x = "Group",
         y = "Correlation Value",
         color = "Tissue") +
    scale_y_continuous(breaks = c(-1, -0.7, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.7, 1)) +
    # scale_color_manual(values = custom_colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(res_folder, gene, ".png"), plot = p, bg = "white")
  # Display the plot
  # print(p)
}