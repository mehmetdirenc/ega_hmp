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
            "#SBATCH --cpus-per-task=16\n"
            "#SBATCH --partition=cpu3-long\n"
            "#SBATCH --time=1-00:55\n"
            "#SBATCH --array=0-19%%3\n"
            "#SBATCH --job-name=ega_rnaseq_%s\n\n\n"
            "source /mnt/lustre/home/mager/magmu818/anaconda3/etc/profile.d/conda.sh\n"
            "conda activate /mnt/lustre/groups/mager/magmu818/.conda/envs/bacass\n\n"
            "cd /mnt/lustre/projects/mager-1000ibd/nextflow_base\n\n"
            "sample=$(sed -n \"$((SLURM_ARRAY_TASK_ID + 1))p\" %s)\n\" %s "
            "/mnt/lustre/projects/mager-1000ibd/nextflow_base/ega_rnaseq_nfcore.sh "
            "%s "
            "%s "
            "/mnt/lustre/home/mager/magmu818/datasets/public_databases/human_genome/GCF_000001405.40_GRCh38.p14_genomic.gtf "
            "/mnt/lustre/home/mager/magmu818/datasets/public_databases/human_genome/GCF_000001405.40_GRCh38.p14_genomic.fna"%(sample_id, samples_list_path, samplesheet_path, result_dir)
        )
        slurm_job_sh.write(job)
    return slurm_path


def nf_core_write_slurm_all(sample_to_reads_dict):
    all_slurm_path = "/mnt/lustre/projects/mager-1000ibd/slurm/ega/rnaseq/ega_all_slurm_job.sh"

    samples_list_path = '/mnt/lustre/projects/mager-1000ibd/slurm/ega/rnaseq/ega_slurm_samples.txt'
    with open(samples_list_path, 'w') as f:
        c = 0
        for sample in sample_to_reads_dict:
            c += 1
            f.write(f"{sample}\n")
    print(str(c))
    with open(all_slurm_path, "w") as slurm_job_sh:
        job = \
        (
            "#!/bin/bash\n"
            "#SBATCH --ntasks=1\n"
            "#SBATCH --cpus-per-task=1\n"
            "#SBATCH --partition=cpu3-long\n"
            "#SBATCH --time=1-10:55\n"
            "#SBATCH --array=0-439%44\n"
            "#SBATCH --job-name=ega_rnaseq_all\n\n\n"
            "source /mnt/lustre/home/mager/magmu818/anaconda3/etc/profile.d/conda.sh\n"
            "conda activate /mnt/lustre/groups/mager/magmu818/.conda/envs/bacass\n\n"
            "echo \"SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID\"\n"
            "cd /mnt/lustre/projects/mager-1000ibd/nextflow_base\n\n"
            "sample=$(sed -n \"$((SLURM_ARRAY_TASK_ID + 1))p\" /mnt/lustre/projects/mager-1000ibd/slurm/ega/rnaseq/ega_slurm_samples.txt)\n"
            "echo \"Sample: $sample\"\n"
            "sample_sheet_path=/mnt/lustre/projects/mager-1000ibd/input/ega/rnaseq/EGAD00001008214/$sample/samplesheet.csv\n"
            "result_dir=/mnt/lustre/projects/mager-1000ibd/results/ega/rnaseq/EGAD00001008214_all_options/$sample\n"
            "/mnt/lustre/projects/mager-1000ibd/nextflow_base/ega_rnaseq_nfcore_all_options.sh "
            "$sample_sheet_path "
            "$result_dir "
            "/mnt/lustre/home/mager/magmu818/datasets/public_databases/human_genome/GCF_000001405.40_GRCh38.p14_genomic.gtf "
            "/mnt/lustre/home/mager/magmu818/datasets/public_databases/human_genome/GCF_000001405.40_GRCh38.p14_genomic.fna"
        )
        slurm_job_sh.write(job)
    return all_slurm_path



def run_slurm(dataset_folder_path):
    sample_to_reads_dict = reads_to_samples(dataset_folder_path)
    print(len(sample_to_reads_dict))
    all_slurm_path = nf_core_write_slurm_all(sample_to_reads_dict)
    subprocess.run("sbatch %s" % all_slurm_path, shell=True, cwd="/mnt/lustre/projects/mager-1000ibd/slurm/ega/rnaseq/slurm_logs_all_options")




if __name__ == '__main__':
    samples_tsv_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008214/metadata/samples.tsv"
    dataset_folder_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008214/"
    samplesheets_path = "/mnt/lustre/projects/mager-1000ibd/input/ega/rnaseq/EGAD00001008214/"
    run_slurm(dataset_folder_path)
