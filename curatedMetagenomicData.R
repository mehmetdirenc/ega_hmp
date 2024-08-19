# BiocManager::install("curatedMetagenomicData")
# browseVignettes("curatedMetagenomicData")
# install.packages("table1")
# BiocManager::install("ExperimentHub")
#
# BiocManager::install("Biobase")
library(SummarizedExperiment)
library(ExperimentHub)
library(curatedMetagenomicData)
library(dplyr)
library(table1)
library(DT)
library(Biobase)
# all_stool_healthy_metadata <- sampleMetadata
# body_sites <- as.list(unique(sampleMetadata$body_site))
# study_names <- as.list(unique(sampleMetadata$study_name))
# curatedMetagenomicData("")
# curatedMetagenomicData("AsnicarF_20.+")
#

eh <- ExperimentHub()

# Query for curatedMetagenomicData
query_results <- query(eh, "curatedMetagenomicData")

# Display the available datasets
query_results

# Assuming we choose a specific dataset ID from the query results
dataset_id <- query_results$ah_id[1]

# Load the dataset
eset <- eh[[dataset_id]]

# View the sample metadata
sample_metadata <- pData(eset)

# Filter for healthy stool samples (assuming health status is indicated in the metadata)
healthy_samples <- sample_metadata[sample_metadata$health_status == "healthy", ]

# Extract the abundance data for healthy samples
abundance_data <- exprs(eset)

# Subset abundance data to only include healthy samples
healthy_abundance_data <- abundance_data[, rownames(healthy_samples)]

# Find the row corresponding to "Clostridium perfringens"
taxa <- rownames(healthy_abundance_data)
clostridium_index <- grep("perfringens", taxa)

# Check if "Clostridium perfringens" is present in healthy samples
if (length(clostridium_index) > 0) {
  clostridium_data <- healthy_abundance_data[clostridium_index, ]
  print(clostridium_data)
} else {
  print("Clostridium perfringens not found in the data.")
}



#
# stool_healthy_marker_presence <- filter(sampleMetadata,study_name == "AsnicarF_2017", disease == "healthy", body_site == "stool") %>%
#   returnSamples(dataType = "marker_presence",
#                 rownames = "short")
#
# table1::table1( ~ disease + age + gender + subject_id + study_name,
#                 data = colData(stool_healthy_marker_presence))
#
# prevalences <- rowSums(assay(stool_healthy_marker_presence) > 0) / ncol(stool_healthy_marker_presence)
# prevalences <- tibble(species = names(prevalences), prevalence = signif(prevalences, 2)) %>%
#   filter(prevalence > 0) %>%
#   arrange(-prevalence)
# DT::datatable(prevalences)
#
#
# stool_healthy <- filter(sampleMetadata,study_name == "AsnicarF_2017", disease == "healthy", body_site == "stool") %>%
#   returnSamples(dataType = "relative_abundance",
#                 rownames = "short")
#
# table1::table1( ~ disease + age + gender + subject_id + study_name,
#                 data = colData(stool_healthy))
#
# prevalences <- rowSums(assay(stool_healthy) > 0) / ncol(stool_healthy)
# prevalences <- tibble(species = names(prevalences), prevalence = signif(prevalences, 2)) %>%
#   filter(prevalence > 0) %>%
#   arrange(-prevalence)
# DT::datatable(prevalences)



# crc_subset <- filter(sampleMetadata, disease == "healthy", body_site == "stool") %>%
#   returnSamples(dataType = "relative_abundance",
#                 rownames = "short")
#
# table1::table1( ~ disease + age + gender + country + study_name,
#                 data = colData(crc_subset))
#
# prevalences <- rowSums(assay(crc_subset) > 0) / ncol(crc_subset)
# prevalences <- tibble(species = names(prevalences), prevalence = signif(prevalences, 2)) %>%
#   filter(prevalence > 0) %>%
#   arrange(-prevalence)
# DT::datatable(prevalences)

# stool_healthy <- filter(sampleMetadata, disease == "healthy", body_site == "stool") %>%
#   returnSamples(dataType = "relative_abundance",
#                 rownames = "short")
#
# table1::table1( ~ disease + age + gender + country + study_name + subject_id,
#                 data = colData(stool_healthy))
# x <- colData(stool_healthy)
#
# prevalences_stool_healthy <- rowSums(assay(stool_healthy) > 0) / ncol(stool_healthy)
# prevalences_stool_healthy <- tibble(species = names(prevalences_stool_healthy),
#                                     prevalence = signif(prevalences_stool_healthy, 2)) %>%
#   filter(prevalences_stool_healthy > 0) %>%
#   arrange(-prevalences_stool_healthy)
# DT::datatable(prevalences_stool_healthy)


#
#
# metadata_healthy_stool <- sampleMetadata %>%
#   filter(disease == "healthy", body_site == "stool", study_name = "AsnicarF_2017.relative_abundance")
# example_dataset <- curatedMetagenomicData("AsnicarF_20.+.relative_abundance", dryrun = FALSE, counts = TRUE, rownames = "short")
# crc_subset <- filter(sampleMetadata, study_condition == "CRC") %>%
#   returnSamples(dataType = "relative_abundance",
#                 rownames = "short")
# samples_of_interest <- returnSamples(metadata_healthy_stool, "relative_abundance", counts = FALSE, rownames = "long")
#
# table1::table1( ~ disease + disease_subtype + age + gender + country + study_name,
#                 data = colData(example_dataset))