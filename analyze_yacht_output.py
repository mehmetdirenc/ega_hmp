import argparse
import os
import sys


def parse_yacht_result(yacht_result_filepath, org_to_presence_dict):
    specs = os.path.basename(yacht_result_filepath).split('result_')[1].strip(".csv")
    org_to_presence_dict[specs] = {}
    with open(yacht_result_filepath) as yrp:
        for line in yrp:
            if line.startswith("organism_name"):
                continue
            if line.startswith("---"):
                type = line.split(" -")[1].strip()
                if type not in org_to_presence_dict[specs]:
                    org_to_presence_dict[specs][type] = {}
                continue
            split_line = line.strip().split(',')
            org = split_line[0].split("|")[0].split("_contig")[0]
            presence = split_line[7]
            pvalue = split_line[8]
            if org in org_to_presence_dict[specs][type]:
                if org_to_presence_dict[specs][type][org][0] != "TRUE":
                    org_to_presence_dict[specs][type][org] = [presence, pvalue]
            else:
                org_to_presence_dict[specs][type][org] = [presence, pvalue]
    return org_to_presence_dict

def analyze_yacht_results(yacht_results_folder):
    yacht_results = os.listdir(yacht_results_folder)
    org_to_presence_dict = {}
    for yacht_folder in yacht_results:
        if "scale_" not in yacht_folder:
            continue
        specs = os.path.basename(yacht_folder).split("scale_")[1].strip(".csv")
        result_file = os.path.join(yacht_results_folder, yacht_folder, "result_%s.csv" % specs)
        org_to_presence_dict = parse_yacht_result(result_file, org_to_presence_dict)
    return org_to_presence_dict

def write_results(yacht_results_folder, results_summary_filepath):
    org_to_presence_dict = analyze_yacht_results(yacht_results_folder)
    with open(results_summary_filepath, mode='w') as results_summary:
        results_summary.write("specs\torganism_name\tpresence\tpvalue\tscale\tsig\tcoverage\tani\n")
        for specs in org_to_presence_dict:
            scale = specs.split("_")[0]
            sig = specs.split("_")[1]
            ani = specs.split("_")[3]
            for c in org_to_presence_dict[specs]:
                coverage = c.replace("min_coverage", "").replace("_result", "")
                for org in org_to_presence_dict[specs][c]:
                    presence, pvalue = org_to_presence_dict[specs][c][org]
                    results_summary.write("\t".join([specs, org, presence, pvalue, scale, sig, coverage, ani]) + "\n")
    #
    # with open(results_summary_filepath.replace("results_summary.tsv", "results_summary_raw.tsv"), mode='w') as results_summary:
    #     results_summary.write("specs\torganism_name\tpresence\tpvalue\tscale\tsig\tcoverage\tani\n")
    #     for specs in org_to_presence_dict:
    #         scale = specs.split("_")[0]
    #         sig = specs.split("_")[1]
    #         coverage = specs.split("_")[2]
    #         ani = specs.split("_")[3]


if __name__ == '__main__':
    yacht_results_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/metagenomics/yacht/only_barcodes/"
    # analyze_yacht_results(yacht_results_folder)
    write_results(yacht_results_folder, os.path.join(yacht_results_folder, "results_summary.tsv"))
    # parse_yacht_result("/mnt/lustre/projects/mager-1000ibd/results/ega/metagenomics/yacht/only_barcodes/
    # scale_100_50_001_80/result_100_50_001_80.csv")