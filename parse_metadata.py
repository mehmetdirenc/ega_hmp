import os
from pyvenn import *
import pandas as pd
from matplotlib import pyplot as plt

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


def create_venn_patient_ids(datasets, all_phenotypes):
    res_filepath = os.path.join(os.path.dirname(all_phenotypes), "venn.png")
    ds_dict_set = {}
    genes = []
    # ds_frame = pd.read_csv(all_phenotypes, sep="\t")
    for dataset in datasets:
        ds_path = datasets[dataset]
        ds_frame = pd.read_csv(ds_path, sep="\t")
        ds_dict_set[dataset] = set(ds_frame["subject_id"])
    # exp_dict = dict(sorted(exp_dict.items()))
    # for gene in exp_dict:
    #     exp_dict_set[gene] = set(exp_dict[gene])
    #     genes.append(gene)
    labels1 = list(ds_dict_set.values())
    labels, all_collections = get_labels(labels1)
    fig, ax = venn4(labels, names=list(ds_dict_set.keys()))
    fig = fig.tight_layout()
    # plt.show()
    plt.savefig(res_filepath)
    plt.close()
    return all_collections


def create_nfcore_rnaseq_sheet(samples_tsv_path, dataset_folder_path, samplesheet_path):
    sample_to_reads_dict = reads_to_samples(dataset_folder_path)
    sample_to_patient_dict = create_sample_to_patient_dict(samples_tsv_path)
    new_sample_to_patient_dict = {}
    with open(samplesheet_path, "w") as f:
        f.write("sample,fastq_1,fastq_2,strandedness\n")
        for sample in sample_to_reads_dict:

            f.write(",".join([sample, sample_to_reads_dict[sample][0], sample_to_reads_dict[sample][1], "auto\n"]))
            # if sample_to_patient_dict[sample] in new_sample_to_patient_dict:
            #     print(sample_to_patient_dict[sample])
            # new_sample_to_patient_dict[sample_to_patient_dict[sample]] = sample_to_reads_dict[sample]
    return


if __name__ == '__main__':
    samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008214/metadata/samples.tsv"
    dataset_folder_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008214/"
    samplesheet_path = "/mnt/lustre/projects/mager-1000ibd/input/ega/rnaseq/EGAD00001008214/samplesheet.csv"
    create_nfcore_rnaseq_sheet(samples_tsv_path, dataset_folder_path, samplesheet_path)
    print("DONE")
    # samplesheet_path = ""
    # create_nfcore_rnaseq_sheet(out_file_path)
    # dataset_folder_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001004194/"
    # patients_metadata_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003991/EGAF00002487099/EGA_Phenotypes_1000IBD_release_2.txt"
    # samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001004194/metadata/samples.tsv"
    # create_sample_to_patient_dict(samples_tsv_path)
    # reads_to_samples(dataset_folder_path)
    # datasets =\
    # {
    #     "16s_int_bio" : "/home/direnc/inputs/ega/EGAD00001003936_metadata/samples.tsv",
    #     "metagenomics_feces" : "/home/direnc/inputs/ega/EGAD00001004194_metadata/samples.tsv",
    #     "RNA_muc_bio" : "/home/direnc/inputs/ega/EGAD00001008214_metadata/samples.tsv",
    #     "16s_muc_bio" : "/home/direnc/inputs/ega/EGAD00001008215_metadata/samples.tsv"
    # }
    # all_phenotypes = "/home/direnc/inputs/ega/EGA_Phenotypes_1000IBD_release_2.txt"
    # create_venn_patient_ids(datasets, all_phenotypes)