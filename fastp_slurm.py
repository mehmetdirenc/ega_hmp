import subprocess

from parse_metadata import *
def write_slurm_job(slurm_main_folder, result_folder, read1_path, read2_path, sample_accession_id):
    sample_accession_result_folder = os.path.join(result_folder, sample_accession_id)
    if not os.path.exists(sample_accession_result_folder):
        os.mkdir(sample_accession_result_folder)
    fastp_result_folder = os.path.join(sample_accession_result_folder, 'fastp')
    if not os.path.exists(fastp_result_folder):
        os.mkdir(fastp_result_folder)
    fastp_read_1_path = os.path.join(fastp_result_folder, os.path.basename(read1_path))
    fastp_read_2_path = os.path.join(fastp_result_folder, os.path.basename(read2_path))
    header = \
    (
        "#!/bin/bash\n"
        "#SBATCH --ntasks=1\n" 
        "#SBATCH --cpus-per-task=8\n"
        "#SBATCH --job-name=ega_map_abundance\n\n\n"
        "source /mnt/lustre/home/mager/magmu818/anaconda3/etc/profile.d/conda.sh\n"
    )
    if read2_path == "":
        body_fastp = \
            (
                    "conda activate /mnt/lustre/groups/mager/magmu818/.conda/envs/fastp\n"
                    "fastp -w 8 "
                    "-i %s "
                    "-o %s\n"
                    "\n\n\n" % (read1_path, fastp_read_1_path)
            )
    else:
        body_fastp = \
        (

            "conda activate /mnt/lustre/groups/mager/magmu818/.conda/envs/fastp\n"
            "fastp -w 8 "
            "-i %s "
            "-I %s "
            "-o %s "
            "-O %s\n"
            "\n\n\n"%(read1_path, read2_path, fastp_read_1_path, fastp_read_2_path)
        )


    slurm_accession_folder = os.path.join(slurm_main_folder, sample_accession_id)
    if not os.path.exists(slurm_accession_folder):
        os.mkdir(slurm_accession_folder)
    slurm_sh_path = os.path.join(slurm_accession_folder, "run.sh")
    with open(slurm_sh_path, "w") as w:
        w.write(header)
        w.write(body_fastp)
    return slurm_sh_path

def run_slurm_jobs(dataset_folder_path, main_result_folder, main_slurm_folder, paired_or_single):
    sample_to_reads_dict = reads_to_samples(dataset_folder_path)
    for sample_accession_id in sample_to_reads_dict:
        if paired_or_single == "single":
            read1_path = sample_to_reads_dict[sample_accession_id][0]
            read2_path = ""
        else:
            read1_path = sample_to_reads_dict[sample_accession_id][0]
            read2_path = sample_to_reads_dict[sample_accession_id][1]
        slurm_sh_path = write_slurm_job(main_slurm_folder, main_result_folder, read1_path, read2_path, sample_accession_id)
        current_slurm_dir = os.path.dirname(slurm_sh_path)
        # print(slurm_sh_path)
        subprocess.run("sbatch %s"%slurm_sh_path, shell=True, cwd=current_slurm_dir)
        # break

if __name__ == '__main__':
    # dataset_folder_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001008215/"
    # main_result_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/16s/EGAD00001008215/"
    # main_slurm_folder = "/mnt/lustre/projects/mager-1000ibd/slurm/ega/16s/EGAD00001008215"
    # run_slurm_jobs(dataset_folder_path, main_result_folder, main_slurm_folder)
    dataset_folder_path = "/mnt/lustre/projects/mager-1000ibd/datasets/EGAD00001003936/"
    main_result_folder = "/mnt/lustre/projects/mager-1000ibd/results/ega/16s/EGAD00001003936/"
    main_slurm_folder = "/mnt/lustre/projects/mager-1000ibd/slurm/ega/16s/EGAD00001003936"
    paired_or_single = "single"
    run_slurm_jobs(dataset_folder_path, main_result_folder, main_slurm_folder, paired_or_single)



