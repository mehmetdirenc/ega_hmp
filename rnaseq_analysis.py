import os, pathlib
from parse_metadata import *
import pandas as pd



def extract_gene_counts(rnaseq_results_folder):
#     x = "/mnt/lustre/projects/mager-1000ibd/results/ega/rnaseq"\
# "/EGAD00001008214/EGAN00003342087/star_salmon/salmon.merged.gene_counts.tsv"
    samples = os.listdir(rnaseq_results_folder)
    final_df = pd.DataFrame()
    for sample in samples:
        file_path = os.path.join(rnaseq_results_folder, sample, "star_salmon", "salmon.merged.gene_counts.tsv")
        df = pd.read_csv(file_path, sep="\t")
        df = df.drop(columns=['gene_name'])
        # sample  # Assuming the last column is the identifier (e.g., 'EGAN00003342087')
        df = df.rename(columns={sample: sample})
        if final_df.empty:
            final_df = df
        else:
            # Merge the current DataFrame with the final DataFrame on 'gene_id'
            final_df = pd.merge(final_df, df, on='gene_id', how='outer')
    return








if __name__ == '__main__':
    patients_metadata_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003991/EGAF00002487099/EGA_Phenotypes_1000IBD_release_2.txt"
    rna_muc_samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008214/metadata/samples.tsv"
    rna_muc_sample_to_patient_dict = create_sample_to_patient_dict(rna_muc_samples_tsv_path)
    muc_16s_samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008215/metadata/samples.tsv"
    muc_16s_sample_to_patient_dict = create_sample_to_patient_dict(muc_16s_samples_tsv_path)
    int_16s_samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008214/metadata/samples.tsv"
    int_16s_sample_to_patient_dict = create_sample_to_patient_dict(int_16s_samples_tsv_path)
    meta_feces_samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008214/metadata/samples.tsv"
    meta_feces_sample_to_patient_dict = create_sample_to_patient_dict(meta_feces_samples_tsv_path)
    patient_phenotype_dict, header = create_patient_phenotype_dict(patients_metadata_path)
    # reads_to_samples(dataset_folder_path)
    rnaseq_results_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/rnaseq/EGAD00001008214"
    extract_gene_counts(rnaseq_results_folder)
    pass