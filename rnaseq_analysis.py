import os.path
import sys
# from multiprocessing.reduction import duplicate

import numpy as np
from parse_metadata import *
import pandas as pd
from scipy.stats import pearsonr, spearmanr, shapiro
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import re, pickle
pd.options.mode.copy_on_write = True

def extract_gene_counts(rnaseq_results_folder, rnaseq_final_res_path):
    if os.path.exists(rnaseq_final_res_path):
        return rnaseq_final_res_path
#     x = "/mnt/lustre/projects/mager-1000ibd/results/ega/rnaseq"\
# "/EGAD00001008214/EGAN00003342087/star_salmon/salmon.merged.gene_counts.tsv"
    samples = os.listdir(rnaseq_results_folder)
    final_df = pd.DataFrame()
    for sample in samples:
        file_path = os.path.join(rnaseq_results_folder, sample, "star_salmon", "salmon.merged.gene_tpm.tsv")
        df = pd.read_csv(file_path, sep="\t")
        df = df.drop(columns=['gene_name'])
        # sample  # Assuming the last column is the identifier (e.g., 'EGAN00003342087')
        df = df.rename(columns={sample: sample})
        if final_df.empty:
            final_df = df
        else:
            # Merge the current DataFrame with the final DataFrame on 'gene_id'
            final_df = pd.merge(final_df, df, on='gene_id', how='outer')
            final_df2 = final_df.set_index('gene_id').T
    final_df2.to_csv(rnaseq_final_res_path, sep='\t')


def choose_correlation_method(x, y):
    # Perform Shapiro-Wilk test for normality (p < 0.05 means non-normal distribution)
    _, p_value_x = shapiro(x)
    _, p_value_y = shapiro(y)

    # If both variables are normally distributed, use Pearson, otherwise use Spearman
    if p_value_x > 0.05 and p_value_y > 0.05:
        return "pearson"
    else:
        return "spearman"

def check_scatter_plot_df_final(x, y, gene, organism):
    # SCATTER PLOT #
    plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
    sns.scatterplot(x=x, y=y)

    # Set plot labels and title
    plt.xlabel(organism)
    plt.ylabel('%s Gene Expression'%gene)
    plt.title(f'Scatter Plot: {gene} vs {organism}')

    # Display the plot
    plt.show()
    plt.close()
    # return
    # SCATTER PLOT

def generate_rnaseq_df(rnaseq_final_res_path, rna_muc_sample_to_patient_dict):
    df_rna = pd.read_csv(rnaseq_final_res_path, sep="\t")
    df_rna = df_rna.rename(columns={'Unnamed: 0': 'sample'})
    df_rna['sample'] = df_rna['sample'].replace(rna_muc_sample_to_patient_dict)
    return df_rna


def generate_meta_df(meta_summary_results, meta_feces_sample_to_patient_dict):
    df = pd.read_csv(meta_summary_results, sep="\t")
    df = df.rename(columns={'#sample_accession_id': 'sample'})
    df['sample'] = df['sample'].replace(meta_feces_sample_to_patient_dict)
    return df

####EGAD...3936 has almost none of our bacteria of interest
# def generate_16s_int(int_16s_samples_tsv_path, int_16s_sample_to_patient_dict):
#     df = pd.read_csv(int_16s_samples_tsv_path, sep="\t")
#     df = df.rename(columns={'#sample_accession_id': 'sample'})
#     df['sample'] = df['sample'].replace(int_16s_sample_to_patient_dict)
#     return df

def generate_16s_muc(muc_16s_samples_tsv_path, muc_16s_sample_to_patient_dict):
    df = pd.read_csv(muc_16s_samples_tsv_path, sep="\t")
    df = df.rename(columns={'id': 'sample'})
    df = df.rename(columns={'perfringens_sum_percentage': 'Clostridium perfringens percentage'})
    df = df.rename(columns={'hathewaya_sum_percentage': 'Hathewaya massiliensis percentage'})
    df = df.rename(columns={'Duodenibacillus_sum': 'Duodenibacillus massiliensis percentage'})
    selected_orgs = {"vulgatus" : "Phocaeicola vulgatus percentage", "gnavus" : 'Ruminococcus gnavus percentage',
                     "fragilis" : 'Bacteroides fragilis percentage', "prausnitzii_I" : "Faecalibacterium prausnitzii_I percentage",
                     "uniformis" : "Bacteroides uniformis percentage", "bifidum" : "Bifidobacterium bifidum percentage",
                     "Prevotella copri" : "Prevotella copri percentage",
                 "Mediterraneibacter lactaris" : "Mediterraneibacter lactaris percentage",
                 "Fusicatenibacter saccharivorans" : "Fusicatenibacter saccharivorans percentage",
                 "Hungatella effluvii" : "Hungatella effluvii percentage"}
    for so in selected_orgs:
        df = df.rename(columns={f'{so}_sum_percentage': selected_orgs[so]})
    # df = df.rename(columns={'vulgatus_sum_percentage': 'Phocaeicola vulgatus percentage'})
    # df = df.rename(columns={'fragilis_sum_percentage': 'Bacteroides fragilis percentage'})
    # df = df.rename(columns={'prausnitzii_I_sum_percentage': 'Bacteroides fragilis percentage'})
    df['sample'] = df['sample'].replace(muc_16s_sample_to_patient_dict)
    return df

def sum_genes(df_rna, gene_dict):
    for main_gene in list(gene_dict.keys()):
        if main_gene not in gene_dict[main_gene]:
            gene_dict[main_gene].insert(0, main_gene)
    aggregated_data = {'sample': df_rna['sample']}
    for main_gene, sub_genes in gene_dict.items():
        existing_genes = [gene for gene in sub_genes if gene in df_rna.columns]

        # Sum the columns specified in existing_genes for each main_gene
        aggregated_data[main_gene] = df_rna[existing_genes].sum(axis=1)
    df_rna_new = pd.DataFrame(aggregated_data)
    return df_rna_new

def get_100_genes(scatter_folder):
    files = os.listdir(os.path.join(scatter_folder, "threshold_0.001"))
    genes = []
    for file in files:
        gene = file.split("_")[-2]
        if gene not in genes:
            genes.append(gene)
    return genes

def extract_ega_genes(ega_paper_correlation):
    df = pd.read_csv(ega_paper_correlation)
    df_filtered = df[df['Group'].isin([1, 2])]
    df_filtered_uc = df[df['Group'].isin([3])]
    # top_50 = df_filtered.nlargest(500, 'Beta-coefficient')
    # # Get the 50 rows with the smallest values in the 'correlation' column
    # bottom_50 = df_filtered.nsmallest(500, 'Beta-coefficient')
    # # Concatenate the two results to keep only the top 50 and bottom 50 rows
    # df_filtered = pd.concat([top_50, bottom_50]).drop_duplicates().reset_index(drop=True)
    df_filtered['Group'] = df_filtered['Group'].replace({1: 'CD', 2: 'CD'})
    # df = df.drop_duplicates().reset_index(drop=True)
    # top_50_uc = df_filtered_uc.nlargest(500, 'Beta-coefficient')
    # # Get the 50 rows with the smallest values in the 'correlation' column
    # bottom_50_uc = df_filtered_uc.nsmallest(500, 'Beta-coefficient')
    # # Concatenate the two results to keep only the top 50 and bottom 50 rows
    # df_filtered_uc = pd.concat([top_50_uc, bottom_50_uc]).drop_duplicates().reset_index(drop=True)
    df_filtered_uc['Group'] = df_filtered_uc['Group'].replace({3: 'UC'})
    # df = df.drop_duplicates().reset_index(drop=True)
    gene_has_both = df_filtered.groupby('Gene')['Beta-coefficient'].apply(
        lambda x: (x > 0).any() and (x < 0).any()
    )
    genes_with_both_signs = gene_has_both[gene_has_both].index.tolist()
    df_filtered_cd = df_filtered[~df_filtered['Gene'].isin(genes_with_both_signs)]
    gene_has_both = df_filtered_uc.groupby('Gene')['Beta-coefficient'].apply(
        lambda x: (x > 0).any() and (x < 0).any()
    )
    genes_with_both_signs = gene_has_both[gene_has_both].index.tolist()
    df_filtered_uc = df_filtered_uc[~df_filtered_uc['Gene'].isin(genes_with_both_signs)]
    # if len(genes_with_both_signs) > 0:
    #     print("incompatible genes, exiting")
    #     return 1
    return df_filtered_cd, df_filtered_uc

def get_gene_aliases(human_gene_info):
    alias_dict = {}
    with open(human_gene_info, 'r') as f:
        for line in f:
            split_line = line.split('\t')
            if split_line[1] == "-":
                continue
            genes = split_line[1].strip().split('|')
            for gene in genes:
                if gene == "-":
                    continue
                if gene not in alias_dict:
                    alias_dict[gene] = [split_line[0]]
                else:
                    alias_dict[gene].append(split_line[0])
    return alias_dict

def get_correlation(gene_aliases, dlr_type, df_final, gene, ega_paper_cd_100, ega_paper_uc_100, tpm_threshold, organism):
    # for dlr_type in dlrs:
    result_to_add = {"ega_coefficient" : [], "alias" : [], "trustworthy" : [], "gene" : [], "correlation_type" : [], 'dlr_type': [], 'correlation': [], 'p_value': []}
    ega_coefficient = "Na"
    trustworthy = "Yes"
    if dlr_type == "nan" or dlr_type == "Diagnosis reconsidered":
        return {"Na"}
    # Filter the DataFrame for the current dlr type
    df_subset = df_final[df_final['diagnosis_last_record'] == dlr_type]
    all_columns = df_subset.columns
    alias = gene
    if gene not in all_columns:
        if gene not in gene_aliases:
            print(gene)
            return {"Na"}
        else:
            aliases = gene_aliases[gene]

            for alias_t in aliases:
                if alias_t in all_columns:
                    df_subset = df_subset[df_subset[alias_t] >= tpm_threshold]
                    corr, p_value = spearmanr(df_subset[organism],
                                              df_subset[alias_t])
                    alias = alias_t
                    break
        # return results
    else:
        df_subset = df_subset[df_subset[gene] >= tpm_threshold]
        x = df_subset[organism]
        y = df_subset[gene]
        if len(y) < 10:
            return {"Na"}
        try:
            corr, p_value = spearmanr(x, y)
        except:
            print("problematic gene")
            sys.exit()
            return {"Na"}
    if p_value > 0.05 or p_value == "Nan" or p_value == "NaN":
        trustworthy = "No"
    if dlr_type == "UC":
        if gene in ega_paper_uc_100["Gene"].tolist():
            ega_coefficient = ega_paper_uc_100.loc[ega_paper_uc_100['Gene'] == gene, 'Beta-coefficient'].tolist()[0]
    elif dlr_type == "CD":
        if gene in ega_paper_cd_100["Gene"].tolist():
            ega_coefficient = ega_paper_cd_100.loc[ega_paper_cd_100['Gene'] == gene, 'Beta-coefficient'].tolist()[0]
    else:
        if gene in ega_paper_cd_100['Gene'].tolist():
            ega_coefficient = ega_paper_cd_100.loc[ega_paper_cd_100['Gene'] == gene, 'Beta-coefficient'].tolist()[0]
        elif gene in ega_paper_uc_100['Gene'].tolist():
            ega_coefficient = ega_paper_uc_100.loc[ega_paper_uc_100['Gene'] == gene, 'Beta-coefficient'].tolist()[0]
        else:
            ega_coefficient = "Na"
    # if dlr_type == "IBDU":
    #     ega_coefficient = "Na"
    result_to_add['ega_coefficient'].append(ega_coefficient)
    result_to_add['trustworthy'].append(trustworthy)
    result_to_add['correlation_type'].append("spearmann")
    result_to_add['dlr_type'].append(dlr_type)
    result_to_add['gene'].append(gene)
    result_to_add['alias'].append(alias)
    result_to_add['correlation'].append(corr)
    result_to_add['p_value'].append(p_value)
    return result_to_add

def combine_rnaseq_metagenomics(rnaseq_final_res_path, meta_summary_results,
                                correlation_folder, rna_muc_sample_to_patient_dict,
                                meta_feces_sample_to_patient_dict,
                                muc_16s_sample_to_patient_dict, muc_16s_summary_results,
                                threshold, genes_dict, scatter_folder, ega_paper_correlation,
                                human_gene_info, ega_based, tpm_threshold, organism, cp_chosen, abbr):
    chosen_pickle_path = os.path.join(correlation_folder, f"{ega_based}_chosen_genes.pickle")
    for_scatter_pickle_path = os.path.join(correlation_folder, f"{ega_based}for_scatter.pickle")
    # if os.path.exists(chosen_pickle_path) and os.path.exists(for_scatter_pickle_path):
    #     with open(chosen_pickle_path, "rb") as file:
    #         chosen_genes = pickle.load(file)
    #     with open(for_scatter_pickle_path, "rb") as file:
    #         for_scatter = pickle.load(file)
    #     create_scatter_plots(chosen_genes, for_scatter, scatter_folder, threshold, ega_based, tpm_threshold, organism, abbr)
    #     return
    if not os.path.exists(correlation_folder):
        os.mkdir(correlation_folder)
    if not os.path.exists(scatter_folder):
        os.mkdir(scatter_folder)
    gene_aliases = get_gene_aliases(human_gene_info)
    ega_paper_cd_100, ega_paper_uc_100 = extract_ega_genes(ega_paper_correlation)

    ega_cd_genes = ega_paper_cd_100["Gene"].tolist()
    ega_uc_genes = ega_paper_uc_100["Gene"].tolist()
    df_rna = generate_rnaseq_df(rnaseq_final_res_path, rna_muc_sample_to_patient_dict)
    # h_genes = get_100_genes(scatter_folder)
    # Generate means for multiple samples from one patient
    df_rna = df_rna.groupby('sample', as_index=False).mean()
    df_rna = sum_genes(df_rna, genes_dict)
    # df_meta = generate_meta_df(meta_summary_results, meta_feces_sample_to_patient_dict)
    # df_meta = df_meta[df_meta['Clostridium perfringens percentage'] >= threshold]
    df_16s_muc = generate_16s_muc(muc_16s_summary_results, muc_16s_sample_to_patient_dict)
    df_16s_muc = df_16s_muc[df_16s_muc[organism] >= threshold]
    df_results = {}
    dfs = {"muc_16s" : df_16s_muc}
    # genes = ["ATG16L1", "PTGER4", "TNF", "IL6", "IL17", "CXCL8",
    #          "IL23", "IL10", "S100A8", "S100A9", "NOD2", "STAT3",
    #          "REG1A", "REG1B", "DUOXA2", "ANXA10", "MUC5AC", "DUOX2",
    #          "REG1B", "MMP3", "AQP8", "CLDN8", "CDHR1", "SLC38A4", "FMO1"]
    print(len(genes_dict))
    for_scatter = {}
    all_results = {}
    # for cor_type in cor_types:
    for exp in dfs:
        if exp == "metagenomics_feces":
            continue
        df = dfs[exp]
        # all_results = {}
        all_results[exp] = []
        df_final = df_rna.merge(df[['sample', organism, "diagnosis_last_record"]], on='sample', how='inner')
        # df_final = df_rna.merge(df[['sample', 'Clostridium perfringens percentage', 'Hathewaya massiliensis percentage',
        #                             'Ruminococcus gnavus percentage', 'Phocaeicola vulgatus percentage', 'Bacteroides fragilis percentage',
        #                             'Duodenibacillus massiliensis percentage', "Bacteroides uniformis percentage",
        #                             "Bifidobacterium bifidum percentage" , "Faecalibacterium prausnitzii_I percentage",
        #                             "diagnosis_last_record"]], on='sample', how='inner')
        disease_counts_before = df_final['diagnosis_last_record'].value_counts()
        df_final = extra_disease_check(rna_muc_samples_tsv_path, df_final)
        disease_counts_after = df_final['diagnosis_last_record'].value_counts()
        for_scatter[exp] = df_final
        dlrs = df_final['diagnosis_last_record'].unique()
        # columns = df_final.columns
        df = None
        if ega_based == "ega_based":
            cp_cd_genes = ega_cd_genes
            cp_uc_genes = ega_uc_genes
        elif ega_based == "cp_based":
            cp_genes = pd.read_csv(cp_chosen, sep="\t")
            cp_cd_genes = cp_genes[cp_genes["dlr_type"] == "CD"]["gene"].tolist()
            cp_uc_genes = cp_genes[cp_genes["dlr_type"] == "UC"]["gene"].tolist()
        for dlr_type in dlrs:
            if dlr_type == "nan" or dlr_type == "Diagnosis reconsidered":
                continue
            elif dlr_type == "CD" or dlr_type == "IBDU":
                # for gene in ega_cd_genes: #### CHANGED IT TO CLOSTRIDIUM PERFRINGENS GENES #########
                for gene in cp_cd_genes:
                    gene_subtype = gene +"_" + dlr_type
                    result_to_add = get_correlation(gene_aliases, dlr_type, df_final, gene,
                                                    ega_paper_cd_100, ega_paper_uc_100, tpm_threshold, organism)
                    if result_to_add == {"Na"}:
                        continue
                    # try:
                    result_to_add["gene_subtype"] = gene_subtype
                    # except:
                    #     pass
                    all_results[exp].append(result_to_add)
                    # all_results[exp][gene_subtype] = result_to_add
                # print("a")
            elif dlr_type == "UC" or dlr_type == "IBDU":
                # for gene in ega_uc_genes #### CHANGED IT TO CLOSTRIDIUM PERFRINGENS GENES #########:
                for gene in cp_uc_genes:
                    # if gene == "ASB11":
                    #     print(gene)
                    gene_subtype = gene +"_" + dlr_type
                    result_to_add = get_correlation(gene_aliases, dlr_type, df_final, gene,
                                                    ega_paper_cd_100, ega_paper_uc_100, tpm_threshold, organism)
                    if result_to_add == {"Na"}:
                        continue
                    result_to_add["gene_subtype"] = gene_subtype
                    all_results[exp].append(result_to_add)
                    # all_results[exp][gene_subtype] = result_to_add
# CONVERT LISTS TO STRING
        for entry in all_results[exp]:
            for key, value in entry.items():
                entry[key] = flatten_value(value)
#         for entry in all_results[exp]:
#             for key, value in entry.items():
#                 entry[key] = ', '.join(map(str, value)) if isinstance(value, list) else str(value)
        df_results[exp] = pd.DataFrame(all_results[exp]).dropna()
        # try:
        df_results[exp]['adjusted_p_value'] = multipletests(df_results[exp]['p_value'], method='fdr_bh')[1]
        # df_results[exp]['adjusted_p_value'] = multipletests(df_results[exp]['p_value'], method='bonferroni')[1]
        # except:
        #     print("PROBLEMATIC GENE")
        # df_results[exp]['adjusted_p_value_bf'] = multipletests(df_results[exp]['p_value'], method='bonferroni')[1]
        # df_results[exp]['adjusted_p_value_holm'] = multipletests(df_results[exp]['p_value'], method='holm')[1]
        #### GOT RID OF THE 0.05 PVALUE THRESHOLD ###
        # if organism == "Clostridium perfringens percentage":
        if ega_based != "cp_based":
            df_results[exp] = df_results[exp][df_results[exp]['adjusted_p_value'] <= 0.05]
        all_sum_file = os.path.join(correlation_folder, f"{ega_based}_{exp}_{str(threshold)}_all_correlation_summary.tsv")
        df_results[exp].drop_duplicates().to_csv(all_sum_file, index=False, header=True, sep='\t')
        compare_coefficient_to_correlation(all_sum_file)
        # print(df_results[exp].head())
    chosen_genes = create_correlation_plots(df_results, correlation_folder, threshold, ega_based, tpm_threshold, organism)
    with open(chosen_pickle_path, "wb") as file:
        pickle.dump(chosen_genes, file)
    with open(for_scatter_pickle_path, "wb") as file:
        pickle.dump(for_scatter, file)
    create_scatter_plots(chosen_genes, for_scatter, scatter_folder, threshold, ega_based, tpm_threshold, organism, abbr)
    return

def flatten_value(value):
    # If the value is a list and has one element
    if isinstance(value, list) and len(value) == 1:
        element = value[0]
        # Return as int, float, or str based on the element type
        if isinstance(element, (int, float)):
            return element
        else:
            return str(element)
    # If not a list, return as is
    return value

def log10_with_zeros(x):
    return np.where(x == 0, 0, np.log10(x))

def create_scatter_plots(chosen_genes, for_scatter, scatter_folder, threshold, ega_based, tpm_threshold, organism, abbr):
    scatter_res_folder = os.path.join(scatter_folder, 'threshold_' + str(threshold) + "_" + ega_based)
    if not os.path.exists(scatter_folder):
        os.mkdir(scatter_folder)
    if not os.path.exists(scatter_res_folder):
        os.mkdir(scatter_res_folder)
    for exp in chosen_genes:
        cor_sum_file = os.path.join(scatter_res_folder, f"{ega_based}{exp}_{str(threshold)}_all_correlation_summary.tsv")
        df_x = chosen_genes[exp]
        exp = exp[0:-2]
        df_x.to_csv(cor_sum_file, index=False, header=True, sep='\t')
        try:
            df_y = for_scatter[exp]
        except:
            print(exp)
        for _, row in df_x.iterrows():
            gene_name = row['gene']
            dlr_type = row['dlr_type']
            alias = row['alias']
            # Filter df_y based on diagnosis_last_record matching dlr_type
            filtered_df_y = df_y[df_y['diagnosis_last_record'] == dlr_type]
            from scipy.stats import linregress
            if gene_name not in filtered_df_y:
                gene_name = alias
                # x = filtered_df_y[alias]
            # else:
            # x = filtered_df_y[gene_name]
            # e = filtered_df_y[gene_name]
            # e2 = [filtered_df_y[gene_name] >= tpm_threshold]
            filtered_df_y = filtered_df_y[filtered_df_y[gene_name] >= tpm_threshold]
            x = log10_with_zeros(filtered_df_y[gene_name])
            # except:
            #     print(gene_name)
            y = log10_with_zeros(filtered_df_y[organism])
            # y = log10_with_zeros(filtered_df_y["Clostridium perfringens percentage"])
            # Perform linear regression
            slope, intercept, r_value, p_value, std_err = linregress(x, y)

            # Create the regression line
            regression_line = slope * x + intercept
            # Check if the gene_name exists in df_y columns before plotting
            if gene_name not in filtered_df_y.columns:
                print(f"Gene '{gene_name}' not found in df_y columns. Skipping.")
                continue
            # if gene_name == "GPN1" and dlr_type == "CD" and exp == "muc_16s" and threshold == 0.001:
            #     print("asd")
            result_file = os.path.join(scatter_res_folder, f"{exp}_{dlr_type}_{gene_name}_{str(threshold)}.png")
            # Create scatter plot
            plt.figure(figsize=(8, 6))
            plt.scatter(x, y, alpha=0.7)
            # plt.plot(x, regression_line, color='red', label="Regression Line")
            plt.xlabel(gene_name)
            plt.ylabel(organism)
            num_points = str(len(filtered_df_y))
            plt.title(f"Scatter Plot for {gene_name} to {abbr} Abundance in {dlr_type} with {tpm_threshold}, n = {num_points} (x y log10)")
            # plt.grid(True)
            # plt.show()
            plt.savefig(result_file)
            plt.close()
    return


def create_correlation_plots(df_results, correlation_folder, threshold, ega_based, tpm_threshold, organism):
    organisms = {'Clostridium perfringens percentage' : 'Clostridium perfringens percentage',
                 'Phocaeicola vulgatus percentage': "Phocaeicola vulgatus percentage UC_pro_inf CD_anti_inf",
                 'Ruminococcus gnavus percentage': "Ruminococcus gnavus percentage pro_inf",
                 'Bacteroides fragilis percentage': "Bacteroides fragilis percentage pro_inf",
                 "Faecalibacterium prausnitzii_I percentage": "Faecalibacterium prausnitzii_I percentage anti_inf",
                 "Bacteroides uniformis percentage" : "Bacteroides uniformis percentage anti_inf",
                 "Bifidobacterium bifidum percentage" : "Bifidobacterium bifidum anti_inf",
                 "Prevotella copri percentage" : "Prevotella copri pro_inf",
                 "Mediterraneibacter lactaris percentage" : "Mediterraneibacter lactaris mixed",
                 "Fusicatenibacter saccharivorans percentage" : "Fusicatenibacter saccharivorans anti?_inf",
                 "Hungatella effluvii percentage" : "Hungatella effluvii pro_inf"}
    # Mediterraneibacter lactaris sometimes pro sometimes anti
    # Fusicatenibacter saccharivorans should be anti
    # Hungatella effluvii
    organism = organisms[organism]
    uc_genes_up = ["S100A9", "TNF", "IL6", "IL17", "CXCL8", "IL8", "S100A8", "IL1B", "S100A12", "LILRA5", "MMP9"]
    res_folder = os.path.join(correlation_folder, 'threshold_' + str(threshold))
    if not os.path.exists(res_folder):
        os.mkdir(res_folder)
    # sns.set(style="whitegrid")
    chosen_genes ={}
    for exp in df_results:
        df = df_results[exp].drop_duplicates().reset_index(drop=True)
        for dlr in df['dlr_type'].unique():
            if dlr == "IBDU":
                continue
            df_subset = df[df['dlr_type'] == dlr]
            # if dlr == "UC":
            #     extra_rows = df_subset[df_subset['gene'].isin(uc_genes_up)]
            # else:
            extra_rows = df_subset[df_subset['gene'].isin(uc_genes_up)]
            # if ega_based == "cp_based":
            top_50 = df_subset.nlargest(10, 'ega_coefficient')
            bottom_50 = df_subset.nsmallest(10, 'ega_coefficient')
            # elif ega_based == "ega_based":
            #     top_50 = df_subset.nlargest(10, 'ega_coefficient')
            #     bottom_50 = df_subset.nsmallest(10, 'ega_coefficient')
            df_subset = pd.concat([top_50, bottom_50, extra_rows]).drop_duplicates().reset_index(drop=True)
            df_subset = df_subset.sort_values('correlation', ascending=False)
            chosen_genes[exp+dlr] = df_subset
            # Filter the dataframe for each dlr_type
            df_subset['color'] = df_subset.apply(lambda row: 'p' if row["adjusted_p_value"] > 0.05 else 'b' if row['correlation'] * row['ega_coefficient'] > 0 else 'r', axis=1)

            size_scale = 200  # Adjust scale factor for dot size
            df_subset['dot_size'] = np.abs(df_subset['ega_coefficient']) * size_scale

            # if len(df_subset["correlation_type"].unique()) > 1:
            #     print("lNETSNDFOMKASDFAFG,LADFN")
            cor_type = df_subset["correlation_type"].unique()[0]
            # Create the dot plot using seaborn
            plt.figure(figsize=(13, 7))  # Adjust figure size
            min_size = df_subset['dot_size'].min()
            median_size = df_subset['dot_size'].median()
            max_size = df_subset['dot_size'].max()
            sns.scatterplot(data=df_subset, x='gene', y='correlation',
                            palette={'b': 'blue', 'r': 'red', "p" : "purple"},
                            hue='color',
                            legend=False,
                            size='dot_size',
                            sizes=(min_size, max_size),
                            color=df_subset['color'])  # Use manual color assignment
            # sns.scatterplot(data=df_subset, x='gene', y='correlation',
            #                 hue='color', palette={'g': 'green', 'r': 'red'}, legend=False)

            # Add labels and title and result file
            # for size in [min_size, median_size, max_size]:
            plt.scatter([], [], s=min_size, color='black', label=f"{min_size / size_scale:.2f}")
            plt.scatter([], [], s=1, color='none', edgecolor='none', label=" ")
            plt.scatter([], [], s=max_size, color='black', label=f"{max_size / size_scale:.2f}")
            plt.ylim(-1, 1)
            # Place the legend outside the plot
            plt.legend(
                title=f"Abundance threshold:{str(threshold)}\n"
                      f"TPM threshold: {str(tpm_threshold)}\nBeta-coefficient values \n(Original EGA dataset, t-test)",  # Legend title
                loc='center left',  # Position on the right side
                bbox_to_anchor=(1, 0.5),# Fine-tune the placement
                fontsize=10,
                title_fontsize=12,
                frameon=False  # Add a frame around the legend
            )
            if threshold == 0:
                result_file = os.path.join(res_folder, f"{ega_based}_{exp}_{dlr}_{cor_type}.png")
                plt.title(f'{ega_based} {dlr} {organism}', fontsize=16)
            else:
                result_file = os.path.join(res_folder, f"{ega_based}_{exp}_{dlr}_{cor_type}_{str(threshold)}.png")
                plt.title(f'{ega_based} {dlr} {organism}', fontsize=16)
            plt.xlabel('Gene', fontsize=12)
            plt.ylabel('Correlation Value', fontsize=12)

            # Rotate gene labels if necessary for readability
            plt.xticks(rotation=90)
            for label in plt.gca().get_xticklabels():
                if label.get_text() in uc_genes_up:
                    label.set_fontweight('bold')  # Make the label bold
            # Show the plot
            plt.tight_layout()
            # plt.show()
            plt.savefig(result_file)
            plt.close()
    return chosen_genes

def parse_attributes(attributes_str):
    # Regular expression to extract key-value pairs
    attribute_dict = {}
    # Find all key="value" pairs in the attributes string
    for match in re.finditer(r'(\S+) "([^"]+)"', attributes_str):
        key, value = match.groups()
        attribute_dict[key] = value
    return attribute_dict


def get_genes(gtf_file):
    gene_geneid_dict = {}
    gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

    # Rename columns for easier access
    gtf_df.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
    gtf_df['attributes_parsed'] = gtf_df['attributes'].apply(parse_attributes)
    genes = gtf_df[gtf_df['feature'] == 'gene']
    genes['gene_name'] = genes['attributes_parsed'].apply(lambda x: x.get('gene', None))
    genes['gene_id'] = genes['attributes_parsed'].apply(lambda x: x.get('gene_id', None))
    genes = genes[['gene_name', 'gene_id']]
    gene_dict = genes.groupby('gene_name')['gene_id'].apply(list).to_dict()
    return gene_dict
    # open_file = open(gtf_file)
    # parsed_gff = GFF.parse(open_file)
    # for rec in parsed_gff:
    #     print("asd")
    #     continue
    # with open(gtf_file) as gtf:
    #     for line in gtf:
    #         if line.startswith('#') or "gene_id" not in line or "gene" not in line:
    #             continue


def extra_disease_check(rna_muc_samples_tsv_path, df_final):
    df_final['diagnosis_last_record'].replace('nan', pd.NA, inplace=True)
    patient_to_disease = {}
    with open(rna_muc_samples_tsv_path, 'r') as rna_muc_samples_tsv:
        for line in rna_muc_samples_tsv:
            if line.startswith('accession_id'):
                continue
            patient_id = line.split('\t')[5]
            extra_attributes = re.findall(r"\[(.*?)\]", line)[0]
            # for tag in extra_attributes:

            tags = extra_attributes.split("},")
            for i in tags:
                if "disease" in i:
                    attributes = i.split(",")[-1].split("\"")
                    disease = attributes[-2]
                    patient_to_disease[patient_id] = disease
    # dlrs = df_final['diagnosis_last_record'].unique()
    df_final['diagnosis_last_record'] = df_final['diagnosis_last_record'].fillna(
        df_final['sample'].map(patient_to_disease)
    )


    return df_final


def compare_coefficient_to_correlation(correlation_summary):
    dir = os.path.dirname(correlation_summary)
    result_filepath = os.path.join(dir, "ega_coefficient_vs_correlation_percentage.tsv")
    df = pd.read_csv(correlation_summary, sep='\t')
    # uc_genes = cor_genes_df[cor_genes_df['dlr_type'] == "UC"]
    # cd_genes = cor_genes_df[cor_genes_df['dlr_type'] == "CD"]
    UC_accurate = []
    UC_inaccurate = []
    CD_accurate = []
    CD_inaccurate = []
    cd_genes = []
    uc_genes = []
    uc_total_genes = len(df[df["dlr_type"] == "UC"])
    cd_total_genes = len(df[df["dlr_type"] == "CD"])
    # all_gene_count = len(df)
    # Check for sign agreement and populate lists
    for _, row in df.iterrows():
        if row["dlr_type"] == "UC":
            if (row["ega_coefficient"] > 0 and row["correlation"] > 0) or (row["ega_coefficient"] < 0 and row["correlation"] < 0):
                UC_accurate.append(row["gene"])
            else:
                UC_inaccurate.append(row["gene"])
            uc_genes.append(row["gene"])
        elif row["dlr_type"] == "CD":
            if (row["ega_coefficient"] > 0 and row["correlation"] > 0) or (row["ega_coefficient"] < 0 and row["correlation"] < 0):
                CD_accurate.append(row["gene"])
            else:
                CD_inaccurate.append(row["gene"])
            cd_genes.append(row["gene"])

    results = {
        "Category": ["UC_accurate", "UC_inaccurate", "CD_accurate", "CD_inaccurate"],
        "Count": [
            len(UC_accurate),
            len(UC_inaccurate),
            len(CD_accurate),
            len(CD_inaccurate),
        ],
    }
    # results["Percentage"] = [
    #     count / uc_total_genes * 100 if "UC" in category else count / cd_total_genes * 100
    #     for count, category in zip(results["Count"], results["Category"])
    # ]
    results["Percentage"] = [
        count / uc_total_genes * 100 if "UC" in category and uc_total_genes != 0
        else count / cd_total_genes * 100 if "CD" in category and cd_total_genes != 0
        else "N/A"
        for count, category in zip(results["Count"], results["Category"])
    ]
    result_df = pd.DataFrame(results)
    result_df.to_csv(result_filepath, sep='\t', index=False)
    return






if __name__ == '__main__':
    ega_paper_correlation = "/mnt/lustre/projects/mager-1000ibd/results/ega/correlated_genes_tables_ega.csv"
    patients_metadata_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003991/EGAF00002487099/EGA_Phenotypes_1000IBD_release_2.txt"
    rna_muc_samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008214/metadata/samples.tsv"
    rna_muc_sample_to_patient_dict = create_sample_to_patient_dict(rna_muc_samples_tsv_path)
    muc_16s_samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008215/metadata/samples.tsv"
    muc_16s_sample_to_patient_dict = create_sample_to_patient_dict(muc_16s_samples_tsv_path)
    ####EGAD...3936 (int_16s_bio) has almost none of our bacteria of interest
    # int_16s_samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003936/metadata/samples.tsv"
    # int_16s_sample_to_patient_dict = create_sample_to_patient_dict(int_16s_samples_tsv_path)
    meta_feces_samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001004194/metadata/samples.tsv"
    meta_feces_sample_to_patient_dict = create_sample_to_patient_dict(meta_feces_samples_tsv_path)
    patient_phenotype_dict, header = create_patient_phenotype_dict(patients_metadata_path)
    # reads_to_samples(dataset_folder_path)
    rnaseq_results_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/rnaseq/EGAD00001008214_all_options"
    rnaseq_final_res_path = "/mnt/lustre/projects/mager-1000ibd/results/ega/rnaseq/summary_table.tsv"
    summary_results_meta = "/mnt/lustre/projects/mager-1000ibd/results/ega/summaries/metagenomics/EGAD00001004194/EGAD00001004194_summary_final.tsv"
    muc_16s_summary_results = "/mnt/lustre/projects/mager-1000ibd/results/ega/summaries/16s/EGAD00001008215/summary_result_with_percentages.tsv"
    # extract_gene_counts(rnaseq_results_folder, rnaseq_final_res_path)
    gtf_filepath = "/mnt/lustre/home/mager/magmu818/datasets/public_databases/human_genome/GCF_000001405.40_GRCh38.p14_genomic.gtf"
    genes_dict = get_genes(gtf_filepath)
    human_gene_info = "/mnt/lustre/projects/mager-1000ibd/results/ega/human_gene_info.tsv"
    # threshold = 0
    # combine_rnaseq_metagenomics(rnaseq_final_res_path, summary_results_meta,
    #                             correlation_folder, rna_muc_sample_to_patient_dict,
    #                             meta_feces_sample_to_patient_dict,
    #                             muc_16s_sample_to_patient_dict, muc_16s_summary_results,
    #                             threshold, genes_dict, scatter_folder, ega_paper_correlation ,human_gene_info)
    # tpm_thresholds = [0, 0.01, 0.1, 1]
    tpm_thresholds = [0.01, 0.1]
    # tpm_thresholds = [0.1]
    threshold = 0.001
    # correlation_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/correlation_plots_tpm_threshold_"
    # scatter_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/scatter_plots_tpm_threshold_"
    # organisms = {'Clostridium perfringens percentage' : "CP",
    #              'Ruminococcus gnavus percentage' : "RG",
    #              'Phocaeicola vulgatus percentage': "PV"}
    # organisms = {'Phocaeicola vulgatus percentage': "PV", 'Duodenibacillus massiliensis percentage': "DM"}
    organisms = {
                 'Clostridium perfringens percentage' : "CP"
                 # 'Phocaeicola vulgatus percentage': "PV",
                 # 'Ruminococcus gnavus percentage': "RG",
                 # 'Bacteroides fragilis percentage': "BF",
                 # "Faecalibacterium prausnitzii_I percentage": "FP",
                 # "Bacteroides uniformis percentage" : "BU",
                 # "Bifidobacterium bifidum percentage" : "BB",
                 # "Prevotella copri percentage" : "PC",
                 # "Mediterraneibacter lactaris percentage": "ML",
                 # "Fusicatenibacter saccharivorans percentage": "FS"
                 }
    # organisms = {'Bacteroides fragilis percentage': "BF",
    #              "Faecalibacterium prausnitzii_I percentage": "FP",
    #              "Bacteroides uniformis percentage" : "BU",
    #              "Bifidobacterium bifidum percentage" : "BB"}

    # "fragilis", "prausnitzii_I", "uniformis"
    # "Duodenibacillus massiliensis"
    for organism in organisms:
        print(organism)
        abbr = organisms[organism]
        correlation_folder = f"/mnt/lustre/projects/mager-1000ibd/results/ega/correlation_plots/{abbr}"
        scatter_folder = f"/mnt/lustre/projects/mager-1000ibd/results/ega/scatter_plots/{abbr}"
        if not os.path.exists(correlation_folder):
            os.mkdir(correlation_folder)
        if not os.path.exists(scatter_folder):
            os.mkdir(scatter_folder)
        for tpm_threshold in tpm_thresholds:
            cp_chosen = (f"/mnt/lustre/projects/mager-1000ibd/results/ega/CP/correlation_plots_tpm_threshold_{str(tpm_threshold)}/"
                         f"ega_based_muc_16s_{str(threshold)}_all_correlation_summary.tsv")
            correlation_folder = os.path.join(correlation_folder, "correlation_plots_tpm_threshold_" + str(tpm_threshold))
            scatter_folder = os.path.join(scatter_folder, "scatter_plots_tpm_threshold_" + str(tpm_threshold))
            combine_rnaseq_metagenomics(rnaseq_final_res_path, summary_results_meta,
                                        correlation_folder, rna_muc_sample_to_patient_dict,
                                        meta_feces_sample_to_patient_dict,
                                        muc_16s_sample_to_patient_dict, muc_16s_summary_results,
                                        threshold, genes_dict, scatter_folder, ega_paper_correlation ,human_gene_info,
                                        "ega_based", tpm_threshold, organism, cp_chosen, abbr)
            if organism != "Clostridium perfringens percentage":
                combine_rnaseq_metagenomics(rnaseq_final_res_path, summary_results_meta,
                                            correlation_folder, rna_muc_sample_to_patient_dict,
                                            meta_feces_sample_to_patient_dict,
                                            muc_16s_sample_to_patient_dict, muc_16s_summary_results,
                                            threshold, genes_dict, scatter_folder, ega_paper_correlation ,human_gene_info,
                                            "cp_based", tpm_threshold, organism, cp_chosen, abbr)
            # cp_chosen = f"/mnt/lustre/projects/mager-1000ibd/results/ega/CP/{exp}_{str(threshold)}_all_correlation_summary.tsv"
            correlation_folder = f"/mnt/lustre/projects/mager-1000ibd/results/ega/correlation_plots/{abbr}"
            scatter_folder = f"/mnt/lustre/projects/mager-1000ibd/results/ega/scatter_plots/{abbr}"
            # correlation_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/correlation_plots_tpm_threshold_"
            # scatter_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/scatter_plots_tpm_threshold_"
        # break
    # compare_coefficient_to_correlation("/mnt/lustre/projects/mager-1000ibd/results/ega/correlation_plots_tpm_threshold_01/ega_based_muc_16s_0.001_all_correlation_summary.tsv")
    pass



