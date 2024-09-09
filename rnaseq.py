from parse_metadata import *
import subprocess




def create_nfcore_rnaseq_sheet(sample, sample_to_reads_dict):

    samplesheets_path = "/mnt/lustre/projects/mager-1000ibd/input/ega/rnaseq/EGAD00001008214/"
    samplesheet_folder = os.path.join(samplesheets_path, sample)
    if not os.path.isdir(samplesheet_folder):
        os.mkdir(samplesheet_folder)
    samplesheet_path = os.path.join(samplesheet_folder, "samplesheet.csv")
    with open(samplesheet_path, "w") as samplesheet:
        samplesheet.write("sample,fastq_1,fastq_2,strandedness\n")
        samplesheet.write(",".join([sample, sample_to_reads_dict[sample][0], sample_to_reads_dict[sample][1], "auto\n"]))
    return samplesheet_path

def nf_core_write_slurm(samplesheet_path, sample_id):
    results_dir = "/mnt/lustre/projects/mager-1000ibd/results/ega/rnaseq/EGAD00001008214"
    result_dir = os.path.join(results_dir, sample_id)
    if not os.path.isdir(os.path.join("/mnt/lustre/projects/mager-1000ibd/slurm/ega/rnaseq/EGAD00001008214", sample_id)):
        os.mkdir(os.path.join("/mnt/lustre/projects/mager-1000ibd/slurm/ega/rnaseq/EGAD00001008214", sample_id))
    slurm_path = os.path.join("/mnt/lustre/projects/mager-1000ibd/slurm/ega/rnaseq/EGAD00001008214", sample_id, "job.sh")
    with open(slurm_path, "w") as slurm_job_sh:
        job = \
        (
            "#!/bin/bash\n"
            "#SBATCH --ntasks=1\n"
            "#SBATCH --cpus-per-task=12\n"
            "#SBATCH --job-name=ega_rnaseq\n\n\n"
            "source /mnt/lustre/home/mager/magmu818/anaconda3/etc/profile.d/conda.sh\n"
            "conda activate /mnt/lustre/groups/mager/magmu818/.conda/envs/bacass\n\n"
            "cd /mnt/lustre/projects/mager-1000ibd/nextflow_base\n\n"
            "/mnt/lustre/projects/mager-1000ibd/nextflow_base/ega_rnaseq_nfcore.sh "
            "%s "
            "%s "
            "/mnt/lustre/home/mager/magmu818/datasets/public_databases/human_genome/GCF_000001405.40_GRCh38.p14_genomic.gtf "
            "/mnt/lustre/home/mager/magmu818/datasets/public_databases/human_genome/GCF_000001405.40_GRCh38.p14_genomic.fna"%(samplesheet_path, result_dir)
        )
        slurm_job_sh.write(job)
    return slurm_path

def run_slurm(dataset_folder_path):
    sample_to_reads_dict = reads_to_samples(dataset_folder_path)
    for sample in sample_to_reads_dict:
        samplesheet_path = create_nfcore_rnaseq_sheet(sample, sample_to_reads_dict)
        # slurm_main_folder = "/mnt/lustre/projects/mager-1000ibd/slurm/ega/rnaseq/EGAD00001008214"
        # slurm_accession_folder = os.path.join(slurm_main_folder, sample)
        # if not os.path.exists(slurm_accession_folder):
        #     os.mkdir(slurm_accession_folder)
        slurm_sh_path = nf_core_write_slurm(samplesheet_path, sample)
        current_slurm_dir = os.path.dirname(slurm_sh_path)
        subprocess.run("sbatch %s" % slurm_sh_path, shell=True, cwd=current_slurm_dir)
        break
        # sample_to_reads_dict = reads_to_samples(dataset_folder_path)
        # for sample_accession_id in sample_to_reads_dict:
        #     read1_path = sample_to_reads_dict[sample_accession_id][0]
        #     read2_path = sample_to_reads_dict[sample_accession_id][1]
        #     slurm_sh_path = write_slurm_job(main_slurm_folder, main_result_folder, read1_path, read2_path, sample_accession_id)
        #     current_slurm_dir = os.path.dirname(slurm_sh_path)
        #     # print(slurm_sh_path)
        #     subprocess.run("sbatch %s"%slurm_sh_path, shell=True, cwd=current_slurm_dir)
            # break




if __name__ == '__main__':
    samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008214/metadata/samples.tsv"
    dataset_folder_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008214/"
    samplesheets_path = "/mnt/lustre/projects/mager-1000ibd/input/ega/rnaseq/EGAD00001008214/"
    run_slurm(dataset_folder_path)
    # create_nfcore_rnaseq_sheet(samples_tsv_path, dataset_folder_path, samplesheets_path)


