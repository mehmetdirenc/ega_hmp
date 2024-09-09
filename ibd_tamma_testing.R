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
  Correlation_Value_spearman = numeric(),
  p_Value_spearman = numeric(),
  Correlation_Value_pearson = numeric(),
  p_Value_pearson = numeric(),
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
        cor_test_p <- cor.test(filtered_df[[bc]], filtered_df[[gene]], method = "pearson")
        cor_test_s <- cor.test(filtered_df[[bc]], filtered_df[[gene]], method = "spearman")
        print(paste("Group:", selected_group, "Gene:", gene, ", Tissue:", selected_tissue, ", Correlation method:", cor_method, ", Correlation:", cor_val))

      },
      error = function(cond)
      {
        # message(paste("Group:", selected_group, ", Tissue:", selected_tissue))
        ""
      })
      # if (cor_method != "")
      # {
      # cor_val_p <- cor(filtered_df[[bc]], filtered_df[[gene]], method = "pearson")
      # cor_val_s <- cor(filtered_df[[bc]], filtered_df[[gene]], method = "spearman")
      # cor_test_p <- cor.test(filtered_df[[bc]], filtered_df[[gene]], method = "pearson")
      # cor_test_s <- cor.test(filtered_df[[bc]], filtered_df[[gene]], method = "spearman")
      # print(paste("Group:", selected_group, "Gene:", gene, ", Tissue:", selected_tissue, ", Correlation method:", cor_method, ", Correlation:", cor_val))
      # # }
      # else
      # {
      #   cor_val <- cor(filtered_df[[bc]], filtered_df[[gene]], method = "pearson")
      #   # cor_test <- cor.test(filtered_df[[bc]], filtered_df[[gene]], method = "pearson")
      # }
      correlation_value_p <- cor_test_p$estimate
      p_value_p <- cor_test_p$p.value
      correlation_value_s <- cor_test_s$estimate
      p_value_s <- cor_test_s$p.value
      cor_results <- rbind(cor_results, data.frame(
        Gene = gene,
        Group = selected_group,
        Tissue = selected_tissue,
        # Correlation_Method = ifelse(cor_method != "", cor_method, "pearson"),
        Correlation_Value_spearman = correlation_value_s,
        p_Value_spearman = p_value_s,
        Correlation_Value_pearson = correlation_value_p,
        p_Value_pearson = p_value_p,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Filter the data for rows where Tissue is "Colon" and group is "CD"
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
  plot_data <- cor_results %>% dplyr::filter(Gene == gene)

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

























