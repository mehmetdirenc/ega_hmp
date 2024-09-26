from parse_metadata import *
import pandas as pd
from scipy.stats import pearsonr, spearmanr, shapiro
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

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

def check_scatter_plot_df_final(df_final):
    # SCATTER PLOT #
    plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
    sns.scatterplot(x=df_final['Clostridium perfringens percentage'], y=df_final['TNF'])

    # Set plot labels and title
    plt.xlabel('Clostridium perfringens percentage')
    plt.ylabel('TNF Gene Expression')
    plt.title('Scatter Plot: TNF vs Clostridium perfringens percentage')

    # Display the plot
    plt.show()
    return
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
    df['sample'] = df['sample'].replace(muc_16s_sample_to_patient_dict)
    return df

def combine_rnaseq_metagenomics(rnaseq_final_res_path, meta_summary_results,
                                correlation_folder, rna_muc_sample_to_patient_dict,
                                meta_feces_sample_to_patient_dict,
                                muc_16s_sample_to_patient_dict, muc_16s_summary_results, threshold):
    df_rna = generate_rnaseq_df(rnaseq_final_res_path, rna_muc_sample_to_patient_dict)
    # Generate means for multiple samples from one patient
    df_rna = df_rna.groupby('sample', as_index=False).mean()
    df_meta = generate_meta_df(meta_summary_results, meta_feces_sample_to_patient_dict)
    df_meta = df_meta[df_meta['Clostridium perfringens percentage'] >= threshold]
    df_16s_muc = generate_16s_muc(muc_16s_summary_results, muc_16s_sample_to_patient_dict)
    df_16s_muc = df_16s_muc[df_16s_muc['Clostridium perfringens percentage'] >= threshold]
    df_results = {}
    dfs = {"metagenomics_feces" : df_meta, "muc_16s" : df_16s_muc}
    genes = ["ATG16L1", "PTGER4", "TNF", "IL6", "IL17", "CXCL8",
             "IL23", "IL10", "S100A8", "S100A9", "NOD2", "STAT3",
             "REG1A", "REG1B", "DUOXA2", "ANXA10", "MUC5AC", "DUOX2",
             "REG1B", "MMP3", "AQP8", "CLDN8", "CDHR1", "SLC38A4", "FMO1"]
    print(len(df_16s_muc['sample']))
    print(len(set(df_16s_muc['sample'])))
    for exp in dfs:
        df = dfs[exp]
        results = {"gene" : [], "recommended" : [], 'dlr_type': [], 'correlation': [], 'p_value': []}
        df_final = df_rna.merge(df[['sample', 'Clostridium perfringens percentage', 'Hathewaya massiliensis percentage', "diagnosis_last_record"]], on='sample', how='inner')
        # check_scatter_plot_df_final(df_final)
        for gene in genes:
            if gene not in df_final.columns:
                print(gene)
            for dlr_type in df_final['diagnosis_last_record'].unique():
                # Filter the DataFrame for the current dlr type
                df_subset = df_final[df_final['diagnosis_last_record'] == dlr_type]
                all_columns = [col for col in df_subset.columns if col.startswith(gene)]
                # Perform Pearson/Spearman correlation
                for gene_s in all_columns:
                    pearson_corr = "N/A"
                    pearson_p_value = "N/A"
                    spearman_corr = "N/A"
                    spearman_p_value = "N/A"
                    x = df_subset['Clostridium perfringens percentage']
                    y = df_subset[gene_s]
                    correlation_type = choose_correlation_method(x, y)

                    if correlation_type == "pearson":
                        corr, p_value = pearsonr(df_subset['Clostridium perfringens percentage'],
                                                                 df_subset[gene_s])
                        if p_value > 0.05 or p_value == "Nan" or p_value == "NaN":
                            continue
                    if correlation_type == "spearman":
                        corr, p_value = spearmanr(df_subset['Clostridium perfringens percentage'],
                                                                    df_subset[gene_s])
                        if p_value > 0.05 or p_value == "Nan" or p_value == "NaN":
                            continue
                    results['recommended'].append(correlation_type)
                    results['dlr_type'].append(dlr_type)
                    results['gene'].append(gene_s)
                    results['correlation'].append(corr)
                    results['p_value'].append(p_value)
                    # results['pearson_corr'].append(pearson_corr)
                    # results['pearson_p_value'].append(pearson_p_value)
                    # results['spearman_corr'].append(spearman_corr)
                    # results['spearman_p_value'].append(spearman_p_value)

            # Convert the results dictionary to a DataFrame for better readability
        df_results[exp] = pd.DataFrame(results).dropna()
        df_results[exp]['adjusted_p_value'] = multipletests(df_results[exp]['p_value'], method='fdr_bh')[1]
        df_results[exp] = df_results[exp][df_results[exp]['adjusted_p_value'] <= 0.05]
        print(df_results[exp].head())
    create_correlation_plots(df_results, correlation_folder, threshold)
    return


def create_correlation_plots(df_results, correlation_folder, threshold):
    res_folder = os.path.join(correlation_folder, 'threshold_' + str(threshold))
    if not os.path.exists(res_folder):
        os.mkdir(res_folder)
    sns.set(style="whitegrid")
    for exp in df_results:
        df = df_results[exp]
        for dlr in df['dlr_type'].unique():
            # Filter the dataframe for each dlr_type
            df_subset = df[df['dlr_type'] == dlr]
            df_subset = df_subset.sort_values('correlation', ascending=False)
            # Create the dot plot using seaborn
            plt.figure(figsize=(10, 6))  # Adjust figure size
            sns.scatterplot(data=df_subset, x='gene', y='correlation', legend=False,
                            color='b')

            # Add labels and title and result file

            if threshold == 0:
                result_file = os.path.join(res_folder, f"{exp}_{dlr}.png")
                plt.title(f'{exp}  {dlr} with no threshold', fontsize=16)
            else:
                result_file = os.path.join(res_folder, f"{exp}_{dlr}_{str(threshold)}.png")
                plt.title(f'{exp}  {dlr} with threshold {str(threshold)}', fontsize=16)
            plt.xlabel('Gene', fontsize=12)
            plt.ylabel('Correlation Value', fontsize=12)

            # Rotate gene labels if necessary for readability
            plt.xticks(rotation=90)

            # Show the plot
            plt.tight_layout()
            # plt.show()
            plt.savefig(result_file)
            plt.close()




if __name__ == '__main__':
    correlation_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/correlation_plots"
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
    threshold = 0
    combine_rnaseq_metagenomics(rnaseq_final_res_path, summary_results_meta,
                                correlation_folder, rna_muc_sample_to_patient_dict,
                                meta_feces_sample_to_patient_dict,
                                muc_16s_sample_to_patient_dict, muc_16s_summary_results, threshold)

    threshold = 0.01
    combine_rnaseq_metagenomics(rnaseq_final_res_path, summary_results_meta,
                                correlation_folder, rna_muc_sample_to_patient_dict,
                                meta_feces_sample_to_patient_dict,
                                muc_16s_sample_to_patient_dict, muc_16s_summary_results, threshold)
    pass

