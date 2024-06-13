import os
import numpy as np
from parse_metadata import *



def create_metadata(manifest_filepath, metadata_filepath, patients_metadata_path):
    patient_phenotype_dict, header = create_patient_phenotype_dict(patients_metadata_path)
    sample_to_patient_dict = create_sample_to_patient_dict(samples_tsv_path)
    return






if __name__ == '__main__':
    patients_metadata_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003991/EGAF00002487099/EGA_Phenotypes_1000IBD_release_2.txt"
    manifest_filepath = "/mnt/lustre/projects/mager-1000ibd/input/ega/16s/EGAD00001008215/manifest.tsv"
    metadata_filepath = "/mnt/lustre/projects/mager-1000ibd/input/ega/16s/EGAD00001008215/metadata.tsv"
    samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008215/metadata/samples.tsv"
    create_metadata(manifest_filepath, metadata_filepath, patients_metadata_path)
    create_sample_to_patient_dict(samples_tsv_path)
    print("a")