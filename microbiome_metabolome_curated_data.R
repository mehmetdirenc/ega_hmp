library(ggplot2)
library(dplyr)
library(gt)

load_datasets <- function()
{
  # source("/mnt/lustre/home/mager/magmu818/tools/microbiome-metabolome-curated-data/scripts/data_organization/utils.R")
  setwd('/home/direnc/tools/microbiome-metabolome-curated-data/scripts/')
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
    # Calculate the row sums of the original data frame (excluding the first column)
    row_sums <- rowSums(df[,-1])
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
    filtered_tables_list[[df_name]] <- filtered_df
  }
}

generate_dataset_metadata <- function (original_dataset)
{
  ibd_datasets <- c("FRANZOSA_IBD_2019", "iHMP_IBDMDB_2019")
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



  original_dataset$metadata$YACHIDA_CRC_2019 <- original_dataset$metadata$YACHIDA_CRC_2019 %>%
    filter(! Shared.w.ERAWIJANTARI_2020)
  updated.yachida.sample.list <- original_dataset$metadata$YACHIDA_CRC_2019$Sample
  original_dataset$mtb$YACHIDA_CRC_2019 <- original_dataset$mtb$YACHIDA_CRC_2019 %>% filter(Sample %in% updated.yachida.sample.list)
  original_dataset$genera$YACHIDA_CRC_2019 <- original_dataset$genera$YACHIDA_CRC_2019 %>% filter(Sample %in% updated.yachida.sample.list)

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
  # file.lm.results.raw <- "data_analysis/linear_models_genus_metabolite.tsv"
  # file.rem.results <- "data_analysis/rem_results.tsv"
  # file.hmdb.data <- "data_analysis/hmdb_info.tsv"
  # file.cytoscape.network <- "data_analysis/cytoscape_network.tsv"
  # file.cytoscape.nodes <- "data_analysis/cytoscape_node_attributes.tsv"
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

main <- function ()
{
  original_dataset <- load_datasets()
  updated_dataset <- generate_dataset_metadata(original_dataset)
  generate_filtered_tables(updated_dataset)
}

main()