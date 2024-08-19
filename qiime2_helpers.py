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
    df = pd.read_csv(qzv_table_metadata_file, sep='\t', header=[0], dtype='unicode')
    df = df[~df.id.str.startswith('#')]
    return df



def abundance_boi(qzv_table_metadata_file, tax_tsv, patients_metadata_path, samples_tsv_path, boi_path_sum):
    non = []
    patient_phenotype_dict, header = create_patient_phenotype_dict(patients_metadata_path)
    sample_to_patient_dict = create_sample_to_patient_dict(samples_tsv_path)
    sample_patient_phenotype_dict = {}
    for sample in sample_to_patient_dict:
        if sample_to_patient_dict[sample] in patient_phenotype_dict:
            sample_patient_phenotype_dict[sample] = patient_phenotype_dict[sample_to_patient_dict[sample]].replace(sample_to_patient_dict[sample], sample)
        else:
            non.append(sample)
    header = header.replace("patient_id", "id")
    metadata_df = pd.DataFrame([x.strip().split('\t') for x in sample_patient_phenotype_dict.values()])
    metadata_df.columns = header.split("\t")
    non_len = len(non)
    boi_features = {}
    boi_features_only_final = {}
    # boi_path = os.path.join(os.path.dirname(qzv_table_metadata_file), "boi.tsv")
    # boi_path_sum = os.path.join(os.path.dirname(qzv_table_metadata_file), "summed_boi.tsv")
    abundance_df = sample_taxa_abundance(qzv_table_metadata_file)
    feature_id_taxa_dict = parse_seq_to_taxa(tax_tsv)
    abundance_df = abundance_df.apply(pd.to_numeric, errors='coerce').fillna(abundance_df)
    abundance_df['patient_read_sum'] = abundance_df.sum(axis=1, numeric_only=True)
    feature_id_taxa_dict['patient_read_sum'] = 'patient_read_sum'
    for feature in feature_id_taxa_dict:
        if "patient_read_sum" in feature_id_taxa_dict[feature] or "g__Sarcina" in feature_id_taxa_dict[feature] or "g__Hathewaya" in feature_id_taxa_dict[feature] or  "g__Bacillus" in feature_id_taxa_dict[feature]:
            boi_features[feature] = feature_id_taxa_dict[feature]
            boi_features_only_final[feature] = feature_id_taxa_dict[feature].split(";")[-1]
    df_boi = abundance_df[itertools.chain(["id"], list(boi_features.keys()))]
    df_boi.rename(columns=boi_features_only_final, inplace=True)
    df_boi2 = df_boi.apply(pd.to_numeric, errors='coerce').fillna(df_boi)
    # for col in df_boi:
    #     for val in df_boi[col]:
    #         type_val = type(val)

    # df_boi.apply(pd.to_numeric, errors='ignore')
    first_column = df_boi2[['id']]
    dropped = df_boi2.drop(columns= "id")
    df2_filtered = dropped.loc[:, dropped.sum() >= 20]
    df2_filtered = pd.concat([first_column, df2_filtered], axis=1)
    # df2_filtered.to_csv(boi_path, sep="\t")
    df2_filtered["perfringens_sum"] = df2_filtered.filter(like="Sarcina").sum(axis=1)
    df2_filtered["perfringens_sum_percentage"] = (df2_filtered['perfringens_sum']/df2_filtered['patient_read_sum']).mul(100)
    df2_filtered["hathewaya_sum"] = df2_filtered.filter(like="Hathewaya").sum(axis=1)
    df2_filtered["hathewaya_sum_percentage"] = (df2_filtered['hathewaya_sum']/df2_filtered['patient_read_sum']).mul(100)
    df2_filtered.drop(list(df2_filtered.filter(regex='s__')), axis=1, inplace=True)
    merged = pd.merge(df2_filtered, metadata_df, on='id', how='left')
    merged_no_na = merged.dropna()
    merged.to_csv(boi_path_sum, sep="\t", na_rep="NA", index=False)
    return




if __name__ == '__main__':
    patients_metadata_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003991/EGAF00002487099/EGA_Phenotypes_1000IBD_release_2.txt"
    manifest_filepath = "/mnt/lustre/projects/mager-1000ibd/input/ega/16s/EGAD00001008215/manifest.tsv"
    metadata_filepath = "/mnt/lustre/projects/mager-1000ibd/input/ega/16s/EGAD00001008215/metadata.tsv"
    samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008215/metadata/samples.tsv"
    # abundance_table = "/mnt/lustre/projects/mager-1000ibd/results/ega/16s/EGAD00001008215/qiime/all_results/exported_abundance_table/metadata.tsv"
    # tax_tsv = "/mnt/lustre/projects/mager-1000ibd/results/ega/16s/EGAD00001008215/qiime/all_results/exported_taxonomy/taxonomy.tsv"
    # abundance_table = "/home/direnc/results/ega/16s/EGAD00001003936/qiime/all_results/abundance_table.tsv"
    # tax_tsv = "/home/direnc/results/ega/16s/EGAD00001003936/qiime/all_results/exported_taxonomyNoFilt_gtdb_barcodes_v3_v4/taxonomy.tsv"
    # abundance_boi(abundance_table, tax_tsv)
    abundance_table = "/mnt/lustre/projects/mager-1000ibd/results/ega/16s/EGAD00001008215/qiime/all_results/metadata_abundance/metadata.tsv"
    tax_tsv = "/mnt/lustre/projects/mager-1000ibd/results/ega/16s/EGAD00001008215/qiime/all_results/exported_taxonomyNoFilt_gtdb_barcodes_v3_v4/taxonomy.tsv"
    boi_path_sum = "/mnt/lustre/projects/mager-1000ibd/results/ega/summaries/16s/EGAD00001008215/summary_result_with_percentages.tsv"
    abundance_boi(abundance_table, tax_tsv, patients_metadata_path, samples_tsv_path, boi_path_sum)
    # create_metadata(manifest_filepath, metadata_filepath, patients_metadata_path)
    # create_sample_to_patient_dict(samples_tsv_path)
    print("a")