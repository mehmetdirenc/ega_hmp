import subprocess

from parse_metadata import *
def write_slurm_job(slurm_main_folder, result_folder, read1_path, read2_path, sample_accession_id):
    sample_accession_result_folder = os.path.join(result_folder, sample_accession_id)
    if not os.path.exists(sample_accession_result_folder):
        os.mkdir(sample_accession_result_folder)
    kneaddata_result_folder = os.path.join(result_folder, sample_accession_id, 'kneaddata')
    kneaddata_read_1_path = os.path.join(kneaddata_result_folder, os.path.basename(read1_path)[:-10] + "1_kneaddata_paired_1.fastq")
    kneaddata_read_2_path = os.path.join(kneaddata_result_folder, os.path.basename(read2_path)[:-10] + "1_kneaddata_paired_2.fastq")
    unfiltered_kneaddata_read_1_path = os.path.join(kneaddata_result_folder, os.path.basename(read1_path)[:-10] + "1_kneaddata.trimmed.1.fastq")
    unfiltered_kneaddata_read_2_path = os.path.join(kneaddata_result_folder, os.path.basename(read2_path)[:-10] + "1_kneaddata.trimmed.2.fastq")
    kraken_result_folder = os.path.join(result_folder, sample_accession_id, 'kraken')
    header = \
    (
        "#!/bin/bash\n"
        "#SBATCH --ntasks=1\n" 
        "#SBATCH --cpus-per-task=63\n"
        "#SBATCH --job-name=ega_map_abundance\n\n\n"
        "source /mnt/lustre/home/mager/magmu818/anaconda3/etc/profile.d/conda.sh\n"
    )
    body_kneaddata = \
    (

        "conda activate /mnt/lustre/groups/mager/magmu818/.conda/envs/kneaddata\n"
        "kneaddata "
        "--input1 %s "
        "--input2 %s "
        "-db /mnt/lustre/home/mager/magmu818/datasets/public_databases/kneaddata/g38/kndt_human_g38_p14 "
        "--trimmomatic /mnt/lustre/groups/mager/magmu818/.conda/envs/kneaddata/share/trimmomatic-0.39-2/ "
        "--output %s\n"
        "\n\n\n"%(read1_path, read2_path, kneaddata_result_folder)
    )
    body_kraken = \
    (
        "mkdir %s\n"
        "conda activate /mnt/lustre/groups/mager/magmu818/.conda/envs/kraken2\n"
        "kraken2 --db /mnt/lustre/home/mager/magmu818/datasets/public_databases/kraken/bacteria --use-names --paired "
        "--report %s/kraken_out.kreport --threads 63 "
        "--output %s/kraken_out "
        "%s "
        "%s"
        "\n"%(kraken_result_folder, kraken_result_folder, kraken_result_folder, kneaddata_read_1_path, kneaddata_read_2_path)
        +
        "kraken2 --db /mnt/lustre/home/mager/magmu818/datasets/public_databases/kraken/bacteria --confidence 0.5 --use-names --paired "
        "--report %s/unfiltered_kraken_out.kreport --threads 63 "
        "--output %s/unfiltered_kraken_out "
        "%s "
        "%s"
        "\n\n\n"%(kraken_result_folder, kraken_result_folder, unfiltered_kneaddata_read_1_path, unfiltered_kneaddata_read_2_path)
    )

    body_bracken = \
    (
        "conda activate /mnt/lustre/groups/mager/magmu818/.conda/envs/kraken2\n"
        "/mnt/lustre/home/mager/magmu818/bracken/bracken "
        "-d /mnt/lustre/home/mager/magmu818/datasets/public_databases/kraken/bacteria "
        "-i %s/kraken_out.kreport -t 50 "
        "-o %s/bracken_out"%(kraken_result_folder, kraken_result_folder)
    )


    slurm_accession_folder = os.path.join(slurm_main_folder, sample_accession_id)
    if not os.path.exists(slurm_accession_folder):
        os.mkdir(slurm_accession_folder)
    slurm_sh_path = os.path.join(slurm_accession_folder, "run.sh")
    with open(slurm_sh_path, "w") as w:
        w.write(header)
        # w.write(body_kneaddata)
        w.write(body_kraken)
        w.write(body_bracken)
    return slurm_sh_path

def run_slurm_jobs(dataset_folder_path, main_result_folder, main_slurm_folder):
    sample_to_reads_dict = reads_to_samples(dataset_folder_path)
    for sample_accession_id in sample_to_reads_dict:
        read1_path = sample_to_reads_dict[sample_accession_id][0]
        read2_path = sample_to_reads_dict[sample_accession_id][1]
        slurm_sh_path = write_slurm_job(main_slurm_folder, main_result_folder, read1_path, read2_path, sample_accession_id)
        current_slurm_dir = os.path.dirname(slurm_sh_path)
        # print(slurm_sh_path)
        subprocess.run("sbatch %s"%slurm_sh_path, shell=True, cwd=current_slurm_dir)
        # break




if __name__ == '__main__':
    # dataset_folder_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001004194/"
    # main_result_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/metagenomics/EGAD00001004194/"
    # main_slurm_folder = "/mnt/lustre/projects/mager-1000ibd/slurm/ega/metagenomics/EGAD00001004194"
    # run_slurm_jobs(dataset_folder_path, main_result_folder, main_slurm_folder)
    dataset_folder_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001004194/"
    main_result_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/metagenomics/EGAD00001004194/"
    main_slurm_folder = "/mnt/lustre/projects/mager-1000ibd/slurm/ega/metagenomics/EGAD00001004194"
    run_slurm_jobs(dataset_folder_path, main_result_folder, main_slurm_folder)



