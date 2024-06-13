import os

def reads_to_samples(dataset_folder_path):
    sample_file_tsv = os.path.join(dataset_folder_path, "metadata", "sample_file.tsv")
    sample_to_reads_dict = {}
    with open(sample_file_tsv, "r") as f:
        for line in f:
            if line.startswith("sample_accession"):
                continue
            split_line = line.strip().split("\t")
            sample_accession = split_line[0]
            folder_accession = split_line[-1]
            if sample_accession not in sample_to_reads_dict:
                read_1 = split_line[2]
                read_1_filepath = os.path.join(dataset_folder_path, folder_accession, read_1)
                sample_to_reads_dict[sample_accession] = [read_1_filepath]
            else:
                read_2 = split_line[2]
                read_2_filepath = os.path.join(dataset_folder_path, folder_accession, read_2)
                sample_to_reads_dict[sample_accession].append(read_2_filepath)
    return sample_to_reads_dict


def create_sample_to_patient_dict(samples_tsv):
    sample_to_patient_dict = {}
    with open(samples_tsv, "r") as f:
        for line in f:
            if line.startswith("accession"):
                continue
            split_line = line.strip().split("\t")
            sample_id = split_line[0]
            patient_id = split_line[5]
            sample_to_patient_dict[sample_id] = patient_id
    return sample_to_patient_dict

def create_patient_phenotype_dict(phenotypes_tsv):
    patient_phenotype_dict = {}
    with open(phenotypes_tsv, "r") as f:
        for line in f:
            if line.startswith("ID"):
                header = line.replace("ID_1000IBD", "patient_id")
            split_line = line.strip().split("\t")
            patient_id = split_line[0]
            patient_phenotype_dict[patient_id] = line
    return patient_phenotype_dict, header




if __name__ == '__main__':
    dataset_folder_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001004194/"
    patients_metadata_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003991/EGAF00002487099/EGA_Phenotypes_1000IBD_release_2.txt"
    samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001004194/metadata/samples.tsv"
    create_sample_to_patient_dict(samples_tsv_path)
    # reads_to_samples(dataset_folder_path)
