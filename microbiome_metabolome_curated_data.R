library(ggplot2)
library(dplyr)
library(gt)
library(reshape2)

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

  ibd <- c("No", "No", "No", "Yes", "No", "No", "No", "No", "Yes", "No", "No")
  datasets <- c("ERAWIJANTARI_GASTRIC_CANCER_2020",
                "YACHIDA_CRC_2019",
                "KIM_ADENOMAS_2020",
                "FRANZOSA_IBD_2019",
                "MARS_IBS_2020",
                "KANG_AUTISM_2017",
                "JACOBS_IBD_FAMILIES_2016",
                "SINHA_CRC_2016",
                "iHMP_IBDMDB_2019",
                "WANG_ESRD_2020",
                "POYET_BIO_ML_2019")
  sequencing <- c("16S", "16S", "16S", "WGS", "16S", "16S", "16S", "16S", "WGS", "16S", "16S")

  healthy_control_labels <- list("ERAWIJANTARI_GASTRIC_CANCER_2020" = "Healthy",
                                 "YACHIDA_CRC_2019" = "Healthy",
                                 "KIM_ADENOMAS_2020" = "Control",
                                 "FRANZOSA_IBD_2019" = "Control",
                                 "MARS_IBS_2020" = "H",
                                 "KANG_AUTISM_2017" = "Neurotypical",
                                 "JACOBS_IBD_FAMILIES_2016" = "Normal",
                                 "SINHA_CRC_2016" = 0,
                                 "iHMP_IBDMDB_2019" = "nonIBD",
                                 "WANG_ESRD_2020" = "Control",
                                 "POYET_BIO_ML_2019" = "ALL_HEALTHY")
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
    data[ , percentage_columns] <- lapply(data[ , percentage_columns], function(x) ifelse(x < 0.01, NA, x))
    melted_data_ibd <- melt(data, id.vars = c("Sample", "Study.Group"))
    ibd_plot <- ggplot(melted_data_ibd, aes(x = Study.Group, y = value)) +
    geom_point(aes(color = variable), size = 3, alpha = 0.7) + # Convert to percentage format
    labs(x = "Study Group", y = "Percentage", title = "IBD Patients") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste0(plots_path, table_name, "_ibd.png"), plot = ibd_plot)
  }
  # Filter the data for the healthy label
  data <- data[data$Study.Group == study_group_label, ]
  data[ , percentage_columns] <- lapply(data[ , percentage_columns], function(x) ifelse(x < 0.01, NA, x))

  # Melt the data for plotting
  melted_data <- melt(data, id.vars = c("Sample", "Study.Group"))
  # melted_data[ , percentage_columns] <- lapply(melted_data[ , percentage_columns], function(x) ifelse(x < 0.01, NA, x))
  # Generate the plot
  p <- ggplot(melted_data, aes(x = variable, y = value)) +
    geom_point() +
    labs(title = paste("Percentages of healthy feces samples"),
         x = "Species",
         y = "Percentage") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}

# Function to map labels and generate plots
generate_plots <- function(filtered_tables, specific_metadata, plots_path) {
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

        # Save or display the plot
        # print(plot)
        ggsave(paste0(plots_path, table_name, ".png"), plot = plot)
      }
    } else {
      message(paste("Column", healthy_label_column, "not found in metadata"))
    }
  }
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