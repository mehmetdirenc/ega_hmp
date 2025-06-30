library(ggplot2)
library(dplyr)
library(gt)
library(reshape2)
library(patchwork)
library(ggpubr)
library(tidyr)

summarize_bacteria_counts <- function(data) {
  # Rename for consistency
  colnames(data) <- c("Sample", "Study.Group", "CP", "HM")

  # List of bacteria and thresholds
  bacteria <- c("CP", "HM")
  study_groups <- c("UC", "CD", "Control")

  # Initialize result container
  result <- data.frame()

  for (bact in bacteria) {
    for (grp in study_groups) {
      subset <- data %>% filter(Study.Group == grp)
      total_count <- sum(subset[[bact]]  > 0, na.rm = TRUE)
      count_above_0.001 <- sum(subset[[bact]] > 0.001, na.rm = TRUE)
      count_above_0.01  <- sum(subset[[bact]] > 0.01, na.rm = TRUE)

      result <- bind_rows(result, data.frame(
        Bacterium = bact,
        `Study Group` = grp,
        `Total Count` = total_count,
        `Total Count above 0.001` = count_above_0.001,
        `Total Count above 0.01` = count_above_0.01
      ))
    }
  }

  return(result)
}

load_datasets <- function()
{
  # source("/mnt/lustre/home/mager/magmu818/tools/microbiome-metabolome-curated-data/scripts/data_organization/utils.R")
  setwd('/home/direnc/tools/microbiome-metabolome-curated-data/scripts/')
  # setwd("/mnt/lustre/home/mager/magmu818/tools/microbiome-metabolome-curated-data/scripts")
  source("data_organization/utils.R")
  source("data_analysis/hmdb_utils.R")
  options(scipen = 999)

  all.data <- load.all.datasets()
  # for(i in 1:length(all.data)) assign(names(all.data)[i], all.data[[i]])
  # rm(all.data)
  return(all.data)
}


generate_filtered_tables <- function (updated_dataset)
{
  species_counts <- updated_dataset$species.counts
  #bacteria of interest
  cp <- "perfringens"
  hm <- "hathewaya massiliensis"
  filtered_tables_list <- list()

  for (i in seq_along(species_counts))
  {
    df <- species_counts[[i]]
    ##get the keys of the list
    df_name <- names(species_counts)[i]
    # Filter the columns
    filtered_df <- df[, c(TRUE, grepl(cp, colnames(df)[-1], ignore.case = TRUE)
      | grepl(hm, colnames(df)[-1], ignore.case = TRUE))]
    new_colnames <- colnames(filtered_df)
    new_colnames[-1] <- sapply(strsplit(new_colnames[-1], ";"), function(x)
    {
      last_item <- tail(x, 1)
      sub("^s__", "", last_item)
    })
    colnames(filtered_df) <- new_colnames
    for (col_name in colnames(filtered_df)[-1])
    {
      percentage_col_name <- paste0(col_name, "_percentage")
      filtered_df[[percentage_col_name]] <- (filtered_df[[col_name]] / rowSums(df[,-1])) * 100
    }
    # study_group <- updated_dataset$metadata[[df_name]]$Study.Group
    matched_indices <- match(filtered_df$Sample, updated_dataset$metadata[[df_name]]$Sample)
    filtered_df$Study.Group <- updated_dataset$metadata[[df_name]]$Study.Group[matched_indices]
    # filtered_df$Study.Group <- updated_dataset$metadata[[df_name]]$Study.Group
    filtered_tables_list[[df_name]] <- filtered_df
  }
  filtered_tables_list$specific_metadata <- updated_dataset$specific_metadata
  return(filtered_tables_list)
}

generate_dataset_metadata <- function (original_dataset)
{
  ##MARS_IBS_2020 this is an IBS dataset

  ibd <- c("Yes", "Yes")
  datasets <- c("FRANZOSA_IBD_2019",
                "iHMP_IBDMDB_2019")
  sequencing <- c("WGS", "WGS")

  healthy_control_labels <- list("FRANZOSA_IBD_2019" = "Control",
                                 "iHMP_IBDMDB_2019" = "nonIBD")
  
  #   ibd <- c("No", "No", "No", "Yes", "No", "No", "No", "No", "Yes", "No", "No")
  # datasets <- c("ERAWIJANTARI_GASTRIC_CANCER_2020",
  #               "YACHIDA_CRC_2019",
  #               "KIM_ADENOMAS_2020",
  #               "FRANZOSA_IBD_2019",
  #               "MARS_IBS_2020",
  #               "KANG_AUTISM_2017",
  #               "JACOBS_IBD_FAMILIES_2016",
  #               "SINHA_CRC_2016",
  #               "iHMP_IBDMDB_2019",
  #               "WANG_ESRD_2020",
  #               "POYET_BIO_ML_2019")
  # sequencing <- c("16S", "16S", "16S", "WGS", "16S", "16S", "16S", "16S", "WGS", "16S", "16S")

    # healthy_control_labels <- list("ERAWIJANTARI_GASTRIC_CANCER_2020" = "Healthy",
    #                              "YACHIDA_CRC_2019" = "Healthy",
    #                              "KIM_ADENOMAS_2020" = "Control",
    #                              "FRANZOSA_IBD_2019" = "Control",
    #                              "MARS_IBS_2020" = "H",
    #                              "KANG_AUTISM_2017" = "Neurotypical",
    #                              "JACOBS_IBD_FAMILIES_2016" = "Normal",
    #                              "SINHA_CRC_2016" = 0,
    #                              "iHMP_IBDMDB_2019" = "nonIBD",
    #                              "WANG_ESRD_2020" = "Control",
    #                              "POYET_BIO_ML_2019" = "ALL_HEALTHY")
  experiments_metadata <- data.frame(
    study_names = datasets,
    ibd = ibd,
    sequencing = sequencing,
    labels = healthy_control_labels,
    stringsAsFactors = FALSE
  )
  original_dataset$specific_metadata <- experiments_metadata
  return(original_dataset)
}

# Function to create dot plots
create_dotplot <- function(data, study_group_label, table_name) {
  # Filter only the columns that contain "percentage" in their names
  percentage_columns <- grep("percentage", colnames(data), value = TRUE, ignore.case = TRUE)

  # Ensure that the necessary columns are kept: Sample, Study.Group, and percentage columns
  data <- data[, c("Sample", "Study.Group", percentage_columns)]

  ibd_datasets <- c("FRANZOSA_IBD_2019", "iHMP_IBDMDB_2019")
  if (table_name %in% ibd_datasets)
  {
    data <- data %>%
  mutate(Study.Group = ifelse(Study.Group == "nonIBD", "Control", Study.Group))
    print(table_name)
    if (startsWith(table_name, "F"))
    {
      wilcoxon_test(data, table_name, "Control", 0)
      wilcoxon_test(data, table_name, "Control", 0.001)
      fisher_test_result_cd <- fishers_exact_test(data, "CD", "Control")
      fisher_test_result_uc <- fishers_exact_test(data, "UC", "Control")
    }
    else
    {
      wilcoxon_test(data, table_name, "Control", 0)
      wilcoxon_test(data, table_name, "Control", 0.001)
      fisher_test_result_cd <- fishers_exact_test(data, "CD", "Control")
      fisher_test_result_uc <- fishers_exact_test(data, "UC", "Control")
    }
    print("cd")
    print(fisher_test_result_cd$p.value)
    print(fisher_test_result_cd$estimate)
    print(fisher_test_result_cd$conf.int)
    print("uc")
    print(fisher_test_result_uc$p.value)
    print(fisher_test_result_uc$estimate)
    print(fisher_test_result_cd$conf.int)
    # Apply the threshold, replacing values < 0.01 with NA
    data[, percentage_columns] <- lapply(data[, percentage_columns], function(x) ifelse(x < 0.001, NA, x))
    # fisher_test_result <- fishers_exact_test(data)
    # Melt the data for plotting
    melted_data_ibd <- melt(data, id.vars = c("Sample", "Study.Group"))
    melted_data_ibd$value <- log10(melted_data_ibd$value)
    # Calculate number of non-NA values for each Study.Group
    valid_counts <- melted_data_ibd %>%
        group_by(Study.Group) %>%
        summarize(valid_count = sum(!is.na(value)))

    # Calculate total number of patients in each Study.Group
    total_counts <- data %>%
        group_by(Study.Group) %>%
        summarize(total_count = n())

    # Merge the valid and total counts to prepare for custom labels
    counts <- merge(valid_counts, total_counts, by = "Study.Group")

    # Create custom x-axis labels with count and total on separate lines
    new_labels <- paste0(counts$Study.Group, "\nn=", counts$valid_count, "\n(total=", counts$total_count, ")")
    melted_data_ibd$variable <- recode(melted_data_ibd$variable,
                                   "Clostridium_P perfringens_percentage" = "CP",
                                   "Hathewaya massiliensis_percentage" = "HM")
    # Create the plot
    ibd_plot <- ggplot(melted_data_ibd, aes(x = Study.Group, y = value)) +
        geom_point(aes(color = variable), size = 3, alpha = 0.7) +
        labs(x = table_name, y = "Abundance Percentages (Log10)", title = "Both IBD patients and healthy feces samples") +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  # Make x-axis labels horizontal

    # Set custom x-axis labels
    scale_x_discrete(labels = new_labels)
    # Save the plot
    return(ibd_plot)
    # ggsave(paste0(plots_path, table_name, "_ibd_log.png"), plot = ibd_plot)
  }
  # Filter the data for the healthy label
  data <- data[data$Study.Group == study_group_label, ]
  data[ , percentage_columns] <- lapply(data[ , percentage_columns], function(x) ifelse(x < 0.001, NA, x))

  # Melt the data for plotting
  melted_data <- melt(data, id.vars = c("Sample", "Study.Group"))
 # Step 1: Calculate counts for > 0.001 (before log transformation) and total counts
  species_counts <- melted_data %>%
  group_by(variable) %>%
  summarize(
    above_0_001 = sum(10^value > 0.001, na.rm = TRUE),  # Count of samples with values > 0.001 (before log10)
    total_count = n()                                   # Total number of samples per species
  )
  new_labels <- species_counts %>%
  mutate(label = paste0(variable, "\n>0.001 n=", above_0_001, "\ntotal n=", total_count))
  melted_data$value <- log10(melted_data$value)
# Step 3: Generate the plot with customized x-axis labels
  p <- ggplot(melted_data, aes(x = variable, y = value)) +
  geom_point() +
  labs(title = paste("Percentages of healthy feces samples"),
       x = table_name,
       y = "Abundance Percentages (Log10)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_x_discrete(labels = new_labels$label)
  return(p)
}

# Function to map labels and generate plots
generate_plots <- function(filtered_tables, specific_metadata, plots_path) {
  # plot_list <- list()  # Initialize list to hold plots
  # ibd_datasets <- c("FRANZOSA_IBD_2019", "iHMP_IBDMDB_2019")
  for (i in seq_along(filtered_tables[1:6])) {
    table_name <- names(filtered_tables)[i]  # Get the name of the table
    table <- filtered_tables[[i]]

    # Find the corresponding healthy label in specific_metadata
    healthy_label_column <- paste0("labels.", table_name)

    if (healthy_label_column %in% colnames(specific_metadata)) {
      healthy_label <- specific_metadata[1, healthy_label_column]

      # If a matching label is found, generate the plot
      if (!is.na(healthy_label)) {
        plot <- create_dotplot(table, healthy_label, table_name)

        # if (!(table_name %in% ibd_datasets))
        # {
        #   plot_list[[i]] <- plot  # Add plot to list
        # }

        # Save or display the plot
        # print(plot)
        ggsave(paste0(plots_path, table_name, ".png"), plot = plot)
      }
    } else {
      message(paste("Column", healthy_label_column, "not found in metadata"))
    }
  }
  # combined_plot <- wrap_plots(plot_list)  # Combine all plots
  # ggsave(paste0(plots_path, "combined_plot", ".png"), plot = combined_plot)
}


fishers_exact_test <- function (data, condition, healthy_label)
  {
    data_fisher <- data %>%
      mutate(Perfringens_Present = ifelse(`Hathewaya massiliensis_percentage` > 0.001, 1, 0))

    # Summarize the counts of Perfringens presence and absence for each Study Group
    summary_table <- data_fisher %>%
      group_by(Study.Group) %>%
      summarize(Perfringens_Present = sum(Perfringens_Present),
                Perfringens_Absent = sum(Perfringens_Present == 0),
                Total = n())

    # View the summary table (optional)
    # print(summary_table)

    # Extract the counts for CD and Healthy groups
    cd_counts <- summary_table %>% filter(Study.Group == condition)
    healthy_counts <- summary_table %>% filter(Study.Group == healthy_label)

    # Create a contingency table for Fisher's exact test or Chi-square test
    contingency_table <- matrix(c(cd_counts$Perfringens_Present, cd_counts$Total - cd_counts$Perfringens_Present,
                              healthy_counts$Perfringens_Present, healthy_counts$Total - healthy_counts$Perfringens_Present),
                            nrow = 2, byrow = TRUE)

    # Set row and column names for better clarity
    rownames(contingency_table) <- c("CD", "Healthy")
    colnames(contingency_table) <- c("Perfringens Present", "Perfringens Absent")

    # Perform Fisher's Exact Test (recommended for small sample sizes)
    fisher_test_result <- fisher.test(contingency_table)
    # Print the result
    # print(fisher_test_result)
    return(fisher_test_result)
  }


log10_safe <- function(x) {
  ifelse(x == 0, 0, log10(x))
}

wilcoxon_test <- function(data, table_name, healthy_label, threshold) {
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  if (threshold != 0)
  {
    data <- data %>% filter(`Hathewaya massiliensis_percentage` >= threshold)
  }
  else
  {
    print("asd")
    data <- data %>% filter(`Hathewaya massiliensis_percentage` != 0)
    summary_table <- summarize_bacteria_counts(data)
    print(summary_table)
  }
  # Log-transform the data (keeping zeros as log1p(0) = 0)
  data <- data %>%
    mutate(Log_Abundance = log10_safe(`Hathewaya massiliensis_percentage`))

  # Set the x-axis order to CD, Control, UC
  data$Study.Group <- factor(data$Study.Group, levels = c(healthy_label, "CD", "UC"))

  # Calculate the maximum y-value for scaling
  max_y <- max(data$Log_Abundance, na.rm = TRUE) * 1.5

  # Count the number of points per group
  group_counts <- data %>%
    group_by(Study.Group) %>%
    summarise(n = n())

  # Update Study.Group labels to include counts
  data <- data %>%
    mutate(Study.Group = paste0(Study.Group, " (n=", group_counts$n[match(Study.Group, group_counts$Study.Group)], ")"))

  # Perform pairwise Wilcoxon tests
  cd_vs_control <- data %>% filter(grepl("CD", Study.Group) | grepl(healthy_label, Study.Group))
  uc_vs_control <- data %>% filter(grepl("UC", Study.Group) | grepl(healthy_label, Study.Group))

  wilcox_cd_control <- wilcox.test(`Hathewaya massiliensis_percentage` ~ Study.Group, data = cd_vs_control)
  wilcox_uc_control <- wilcox.test(`Hathewaya massiliensis_percentage` ~ Study.Group, data = uc_vs_control)

  # Collect and adjust p-values
  p_values <- c(wilcox_cd_control$p.value, wilcox_uc_control$p.value)
  adjusted_p_values <- p.adjust(p_values, method = "BH")

  # Define annotation labels
  annotated_p_values <- paste0("Adj. p = ", signif(adjusted_p_values, digits = 3))
  # annotated_p_values <- paste0("Adj. p = ", format(adjusted_p_values, scientific = TRUE, digits = 2))

  annotations <- data.frame(
    x = c(2, 3),  # CD vs Control (left) and Control vs UC (right)
    y = rep(max_y * 0.45, 2),
    label = annotated_p_values
  )

# Add formatted labels (Study.Group with n counts)
group_counts <- group_counts %>%
  mutate(label = paste0(Study.Group, " (n=", n, ")"))

# Dynamically identify the reference group (e.g., control/nonIBD)
reference_group <- group_counts$label[group_counts$Study.Group %in% c("nonIBD", "Control")]

# Define pairwise comparisons dynamically
pairwise_comparisons <- group_counts %>%
  filter(!Study.Group %in% c("nonIBD", "Control")) %>%  # Exclude the reference group
  pull(label) %>%                                       # Get labels of other groups
  lapply(function(x) c(x, reference_group))             # Pair with the reference group

group_means_raw <- data %>%
  group_by(Study.Group) %>%
  summarise(Mean_Raw_Abundance = mean(`Hathewaya massiliensis_percentage`, na.rm = TRUE))
group_means_raw <- group_means_raw %>%
  arrange(
    # Custom sorting: first "Control"-starting groups, then all others alphabetically
    factor(Study.Group, levels = c(
      sort(unique(Study.Group[grepl("^Control", Study.Group)])),  # Control-like groups first
      sort(unique(Study.Group[!grepl("^Control", Study.Group)]))  # The rest alphabetically
    ))
  )
  groups <- unique(data$Study.Group)
  sorted_groups <- c(
  sort(groups[grepl("^Control", groups)]),  # Control-like groups first
  sort(groups[!grepl("^Control", groups)])  # The rest alphabetically
)
# Add mean annotations for each group
  group_means_annotations <- data.frame(
  x = unique(sorted_groups),  # Use updated x-axis labels
  y = rep(max_y * 0.50, length(groups)),  # Position of means
  label = paste0("Mean = ", signif(group_means_raw$Mean_Raw_Abundance, digits = 3))
)
# group_means_annotations <- group_means_annotations %>%
#   arrange(x)

data$Study.Group <- factor(data$Study.Group,
                            levels = c(
                              # Control-related groups first
                              grep("^Control", unique(data$Study.Group), value = TRUE),

                              # Then the rest of the groups in alphabetical order
                              sort(setdiff(unique(data$Study.Group), grep("^Control", unique(data$Study.Group), value = TRUE)))
                            ))

# Now plot the data with the updated x-axis order
p <- ggplot(data, aes(x = Study.Group, y = Log_Abundance)) +
  geom_jitter(aes(color = Study.Group),
              width = 0.2, size = 2, alpha = 0.6) +
  geom_violin(aes(fill = Study.Group), alpha = 0.3, outlier.shape = NA) +
  theme_minimal() +
  labs(title = paste(table_name, "_hochberg, threshold: ", threshold, sep = ""),
       x = "Condition",
       y = "Log Scale HM") +
  scale_y_continuous(limits = c(min(data$Log_Abundance, na.rm = TRUE), max_y), expand = c(0, 0))
  # print(p)
  p <- p +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +  # Increase top margin for labels
  geom_text(data = annotations,
            aes(x = x, y = y, label = label),
            color = "red", size = 4, vjust = -0.5)
p <- p +
  # Add group means
  geom_text(data = group_means_annotations, aes(x = x, y = y, label = label),
            color = "darkgreen", size = 4, vjust = -0.5)
  # print(p)

  # Save the plot
  filepath = paste("/home/direnc/results/microbiome_metabolome_curated_data/HM_wilcoxon_", table_name, "_", threshold, ".pdf", sep = "")
  ggsave(filepath, p, width = 10, height = 8, device=cairo_pdf, dpi = 600)

  # Return the test results
  # return(list(
  #   raw_p_values = p_values,
  #   adjusted_p_values_bonferroni = adjusted_p_values
  # ))
}


# main <- function ()
# {
plots_path <- "/home/direnc/results/microbiome_metabolome_curated_data/"
original_dataset <- load_datasets()
updated_dataset <- generate_dataset_metadata(original_dataset)
filtered_tables <- generate_filtered_tables(updated_dataset)
generate_plots(filtered_tables, filtered_tables$specific_metadata, plots_path)
# }

# main()