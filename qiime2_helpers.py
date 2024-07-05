import os
import numpy as np
import pandas as pd
import itertools
from parse_metadata import *



def create_metadata(manifest_filepath, metadata_filepath, patients_metadata_path):
    patient_phenotype_dict, header = create_patient_phenotype_dict(patients_metadata_path)
    sample_to_patient_dict = create_sample_to_patient_dict(samples_tsv_path)
    with open(manifest_filepath, mode='r') as manifest_file, open(metadata_filepath, mode='w') as metadata_file:
        metadata_file.write(patient_phenotype_dict["ID_1000IBD"].replace("ID_1000IBD", "sampleid"))
        for line in manifest_file:
            if line.startswith('sampleid'):
                continue
            sample_id = line.split()[0]
            patient_id = sample_to_patient_dict[sample_id]
            if patient_id in patient_phenotype_dict:
                metadata_file.write(patient_phenotype_dict[patient_id].replace(patient_id, sample_id))
            else:
                print(sample_id)
    return



def parse_seq_to_taxa(tax_tsv):
    feature_id_taxa_dict = {}
    with open(tax_tsv, "r") as tax_feature:
        for line in tax_feature:
            if line.startswith("Feature") or line.startswith("#"):
                continue
            split_line = line.split("\t")
            feature_id = split_line[0]
            full_tax = split_line[1]
            feature_id_taxa_dict[feature_id] = full_tax


def parse_seq_to_taxa(tax_tsv):
    feature_id_taxa_dict = {}
    with open(tax_tsv, "r") as tax_feature:
        for line in tax_feature:
            if line.startswith("Feature") or line.startswith("#"):
                continue
            split_line = line.split("\t")
            feature_id = split_line[0]
            full_tax = split_line[1]
            feature_id_taxa_dict[feature_id] = full_tax
    return feature_id_taxa_dict

def sample_taxa_abundance(qzv_table_metadata_file):
    df = pd.read_csv(qzv_table_metadata_file, sep='\t', header=[0])
    df = df[~df.id.str.startswith('#')]
    return df



def abundance_boi(qzv_table_metadata_file, tax_tsv):
    boi_features = {}
    boi_features_only_final = {}
    abundance_df = sample_taxa_abundance(qzv_table_metadata_file)
    feature_id_taxa_dict = parse_seq_to_taxa(tax_tsv)
    for feature in feature_id_taxa_dict:
        if "g__Sarcina" in feature_id_taxa_dict[feature] or "g__Hathewaya" in feature_id_taxa_dict[feature] or  "g__Bacillus" in feature_id_taxa_dict[feature]:
            boi_features[feature] = feature_id_taxa_dict[feature]
            boi_features_only_final[feature] = feature_id_taxa_dict[feature].split(";")[-1]
    df_boi = abundance_df[itertools.chain(["id"], list(boi_features.keys()))]
    df_boi.rename(columns=boi_features_only_final, inplace=True)
    df_boi.to_csv('/home/direnc/results/ega/16s/qiime/all_results/boi.tsv', sep="\t")
    return




if __name__ == '__main__':
    patients_metadata_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003991/EGAF00002487099/EGA_Phenotypes_1000IBD_release_2.txt"
    manifest_filepath = "/mnt/lustre/projects/mager-1000ibd/input/ega/16s/EGAD00001008215/manifest.tsv"
    metadata_filepath = "/mnt/lustre/projects/mager-1000ibd/input/ega/16s/EGAD00001008215/metadata.tsv"
    samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008215/metadata/samples.tsv"
    # abundance_table = "/mnt/lustre/projects/mager-1000ibd/results/ega/16s/EGAD00001008215/qiime/all_results/exported_abundance_table/metadata.tsv"
    # tax_tsv = "/mnt/lustre/projects/mager-1000ibd/results/ega/16s/EGAD00001008215/qiime/all_results/exported_taxonomy/taxonomy.tsv"
    abundance_table = "/home/direnc/results/ega/16s/qiime/all_results/exported_abundance_table/metadata.tsv"
    tax_tsv = "/home/direnc/results/ega/16s/qiime/all_results/exported_taxonomy/taxonomy.tsv"
    # sample_taxa_abundance(abundance_table)
    abundance_boi(abundance_table, tax_tsv)
    # create_metadata(manifest_filepath, metadata_filepath, patients_metadata_path)
    # create_sample_to_patient_dict(samples_tsv_path)
    print("a")