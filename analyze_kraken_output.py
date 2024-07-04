import os, sys, time
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


def analyze_summary(summary_result_file, min_read_threshold, min_percentage_threshold):
    df = pd.read_csv(summary_result_file, sep='\t', header=0)
    scatter_plot(df, min_read_threshold, min_percentage_threshold)
    return

def scatter_plot(df, min_read_threshold, min_percentage_threshold):
    # explanations = {"Montreal_A" : {"NA" : "NA",
    #     "Not applicable" : "Not applicable", "A1" : "<17 age",
    #                                 "A2" : "17-40 age",
    #                                 "A3" : ">40 age"
    #                  }}
    if min_read_threshold == 0:
        result_folder = os.path.join(os.path.dirname(summary_result_file), "normal_scale", "raw")
        result_folder_log = os.path.join(os.path.dirname(summary_result_file), "log_scale", "raw")
    else:
        result_folder = os.path.join(os.path.dirname(summary_result_file), "normal_scale", "min_%sr_%sp"%(str(min_read_threshold), str(min_percentage_threshold).split(".")[1]))
        result_folder_log = os.path.join(os.path.dirname(summary_result_file), "log_scale",
                                     "min_%sr_%sp" % (str(min_read_threshold), str(min_percentage_threshold).split(".")[1]))
    # if os.path.exists(result_folder) == False:
    #     os.mkdir(result_folder)
    os.makedirs(result_folder_log, exist_ok=True)
    os.makedirs(result_folder, exist_ok=True)
    header = list(df)
    diagnosis_list = df['diagnosis_last_record'].unique()
    for bacterial_reads in header[1:13]:
        if "_BP" in bacterial_reads or "Bacillus" in bacterial_reads or "Hathewa" in bacterial_reads:
            continue
        if "barcode" in bacterial_reads:
            fp = bacterial_reads.split()[0] + bacterial_reads.split()[1][0]
        else:
            fp = bacterial_reads.split()[0][0] + bacterial_reads.split()[1][0] + bacterial_reads.split()[2][0]
        for diagnosis in diagnosis_list:
            for condition in [header[24], header[25], header[26], header[27], header[28], header[21]]:
                plt_path = os.path.join(result_folder, fp + "_" + diagnosis + "_" + condition.replace(" ", "_")) + ".png"
                plt_path_log = os.path.join(result_folder_log, fp + "_" + diagnosis + "_" + condition.replace(" ", "_")) + ".png"
                # if "barcode8_CPp_CD_Montreal_B" not in plt_path:
                #     if "barcode8_CPr_CD_Montreal_B" not in plt_path:
                #         continue
                if min_read_threshold == 0:
                    if "read count" in bacterial_reads:
                        # & (df['#sample_accession_id'] != "EGAN00001762065")
                        # & (df['#sample_accession_id'] != "EGAN00001762065")
                        new_df = df.loc[df['diagnosis_last_record'].eq(diagnosis) & df[bacterial_reads].ge(min_read_threshold) & df[bacterial_reads].le(2000)]
                    else:
                        new_df = df.loc[df['diagnosis_last_record'].eq(diagnosis) & df[bacterial_reads].ge(min_percentage_threshold) & df[bacterial_reads].le(0.01)]
                else:
                    if "read count" in bacterial_reads:
                        new_df = df.loc[df['diagnosis_last_record'].eq(diagnosis) & df[bacterial_reads].ge(min_read_threshold)]
                    else:
                        new_df = df.loc[df['diagnosis_last_record'].eq(diagnosis) & df[bacterial_reads].ge(min_percentage_threshold)]
                x_u_list = list(map(str, new_df[condition].unique().tolist()))
                if "nan" in x_u_list:
                    x_u_list.remove("nan")
                if "Not applicable" in x_u_list:
                    x_u_list.remove("Not applicable")
                if len(x_u_list) == 0:
                    continue
                x_u_list.sort()
                fig, ax = plt.subplots()
                ax.set(yscale="log")
                sns.boxplot(data=new_df, x=condition, y=bacterial_reads, order=x_u_list, ax=ax).set_title(diagnosis)
                sns.stripplot(data=new_df, x=condition, y=bacterial_reads, order=x_u_list, ax=ax, color="red", s=5).set_title(diagnosis)
                plt.savefig(plt_path_log, bbox_inches='tight')
                plt.clf()
                plt.cla()
                plt.close()
                fig, ax = plt.subplots()
                sns.boxplot(data=new_df, x=condition, y=bacterial_reads, order=x_u_list).set_title(diagnosis)
                sns.stripplot(data=new_df, x=condition, y=bacterial_reads, order=x_u_list, color="red", s=5).set_title(diagnosis)
                plt.savefig(plt_path, bbox_inches='tight')
                plt.clf()
                plt.cla()
                # plt.show()
                plt.close()



if __name__ == '__main__':
    results_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/metagenomics/EGAD00001004194"
    summary_result_file = "/mnt/lustre/projects/mager-1000ibd/results/ega/summaries/EGAD00001004194_summary_sr_corrected.tsv"
    samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001004194/metadata/samples.tsv"
    patients_metadata_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003991/EGAF00002487099/EGA_Phenotypes_1000IBD_release_2.txt"
    # patient_phenotype_dict, header = create_patient_phenotype_dict(patients_metadata_path)
    # sample_to_patient_dict = create_sample_to_patient_dict(samples_tsv_path)
    # create_sample_to_patient_dict(samples_tsv_path)
    write_summary(results_folder, samples_tsv_path, patients_metadata_path, summary_result_file)
    min_read_threshold = 50
    min_percentage_threshold = 0.001
    analyze_summary(summary_result_file, min_read_threshold, min_percentage_threshold)
    # min_read_threshold = 0
    # min_percentage_threshold = 0
    # analyze_summary(summary_result_file, min_read_threshold, min_percentage_threshold)
    # print("a")