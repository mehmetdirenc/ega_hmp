import os, sys
from parse_metadata import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from parse_metadata import *
def parse_kraken_report(kreport_filepath):
    tax_to_attributes = {}
    tax_of_interest = {"9999999" : "barcode9_HM", "9999998" : "barcode8_CP", "9999997" : "barcode7_BP",
                       "2026187" : "Bacillus pacificus", "1964382" : "Hathewaya massiliensis",
                       "1502" : "Clostridium perfringens"}
    for tax in tax_of_interest:
        tax_to_attributes[tax] = {"assigned_reads": "0", "percentage" : "0"}

    with open(kreport_filepath) as kf:
        for line in kf:
            # lower_line = line.lower()
            split_line = line.split()
            taxid = split_line[4]
            if line.split()[4] == "2":
                bacterial_reads = split_line[1]
                tax_to_attributes["bacterial_reads"] = bacterial_reads
            if taxid in tax_of_interest:
                # if taxid not in tax_to_attributes:
                #     tax_to_attributes[taxid] = {"assigned_reads": 0, "percentage" : 0}
                total_fragments_assigned = split_line[1]
                tax_to_attributes[taxid]["assigned_reads"] = total_fragments_assigned
                tax_to_attributes[taxid]["percentage"] = str(float(total_fragments_assigned) / float(bacterial_reads) * 100)
    return tax_to_attributes

def write_summary(results_folder, samples_tsv_path, patients_metadata_path, summary_result_file):
    results = os.listdir(results_folder)
    sample_to_patient_dict = create_sample_to_patient_dict(samples_tsv_path)
    patient_phenotype_dict, header = create_patient_phenotype_dict(patients_metadata_path)
    with open(summary_result_file, "w") as summary_file:
        summary_file.write("#sample_accession_id\tbarcode7_BP read count\tbarcode7_BP percentage\tbarcode8_CP read count\tbarcode8_CP percentage"
                           "\tbarcode9_HM read count\tbarcode9_HM percentage\tBacillus pacificus read count\tBacillus pacificus percentage\tClostridium perfringens read count"
                           "\tClostridium perfringens percentage\tHathewaya massiliensis read count\tHathewaya massiliensis percentage\t%s"%header)
        for egan_result in results:
            kraken_report = os.path.join(results_folder, egan_result, "kraken", "kraken_out.kreport")
            tax_to_attributes = parse_kraken_report(kraken_report)
            summary_file.write(egan_result + "\t")
            summary_file.write("\t".join([tax_to_attributes["9999997"]["assigned_reads"], tax_to_attributes["9999997"]["percentage"],
                                          tax_to_attributes["9999998"]["assigned_reads"], tax_to_attributes["9999998"]["percentage"],
                                          tax_to_attributes["9999999"]["assigned_reads"], tax_to_attributes["9999999"]["percentage"],
                                          tax_to_attributes["2026187"]["assigned_reads"], tax_to_attributes["2026187"]["percentage"],
                                          tax_to_attributes["1502"]["assigned_reads"], tax_to_attributes["1502"]["percentage"],
                                          tax_to_attributes["1964382"]["assigned_reads"], tax_to_attributes["1964382"]["percentage"]]))
            summary_file.write("\t%s"%patient_phenotype_dict[sample_to_patient_dict[egan_result]])
            # summary_file.write("\n")
            # return


def analyze_summary(summary_result_file):
    df = pd.read_csv(summary_result_file, sep='\t', header=0)
    header = list(df)
    b8_cp_perc = df["barcode8_CP percentage"]
    b8_cp_count = df["barcode8_CP read count"]
    diagnosis = df["diagnosis_last_record"]
    print(diagnosis)
    sns.scatterplot(data=df, x=df.index, y=b8_cp_count, hue=diagnosis)
    plt.savefig('/mnt/lustre/projects/mager-1000ibd/results/ega/summaries/output.png')
    for bacterial_reads in header[1:13]:
    return

if __name__ == '__main__':
    results_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/metagenomics/EGAD00001004194"
    summary_result_file = "/mnt/lustre/projects/mager-1000ibd/results/ega/summaries/EGAD00001004194_summary.tsv"
    samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001004194/metadata/samples.tsv"
    patients_metadata_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003991/EGAF00002487099/EGA_Phenotypes_1000IBD_release_2.txt"
    # patient_phenotype_dict, header = create_patient_phenotype_dict(patients_metadata_path)
    # sample_to_patient_dict = create_sample_to_patient_dict(samples_tsv_path)
    # create_sample_to_patient_dict(samples_tsv_path)
    # write_summary(results_folder, samples_tsv_path, patients_metadata_path, summary_result_file)
    analyze_summary(summary_result_file)
    # print("a")