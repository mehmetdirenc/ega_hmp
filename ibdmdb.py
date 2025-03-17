import os
from qiime2_helpers import abundance_boi


def create_manifest(manifest_filepath, data_folder, sra_runtable_filepath):
    srr_patient_dict = parse_runtable(sra_runtable_filepath)
    srr_list = os.listdir(data_folder)
    with open(manifest_filepath, 'w') as manifest_file:
        manifest_file.write('sampleid\tforward-absolute-filepath\treverse-absolute-filepath\n')
        for srr in srr_patient_dict:
            patient = srr_patient_dict[srr]
            if srr + "_1.fastq.gz" not in srr_list or srr + "_2.fastq.gz" not in srr_list:
                print(srr)
            else:
                fw = os.path.join(data_folder, srr + "_1.fastq.gz")
                rw = os.path.join(data_folder, srr + "_2.fastq.gz")
                manifest_file.write(patient + "\t" + fw + "\t" + rw + "\n")
    return





def parse_metadata(metadata_filepath):
    return

def parse_runtable(sra_runtable_filepath):
    srr_id_patient_dict = {}
    with open(sra_runtable_filepath) as sra_runtable:
        for line in sra_runtable:
            split_line = line.split(',')
            if line.startswith('Run') or split_line[1] != "AMPLICON":
                continue
            patient_id = split_line[29]
            srr_id = split_line[0]
            srr_id_patient_dict[srr_id] = patient_id
    return srr_id_patient_dict





if __name__ == '__main__':
    # manifest_filepath = "/mnt/lustre/home/mager/magmu818/inputs/ibdmdb/16s/manifest.tsv"
    # metadata_filepath = "/mnt/lustre/home/mager/magmu818/inputs/ibdmdb/16s/metadata.tsv"
    # sra_runtable_filepath = "/mnt/lustre/home/mager/magmu818/inputs/ibdmdb/16s/SraRunTable.txt"
    # data_folder = "/mnt/lustre/home/mager/magmu818/datasets/ibdmdb/raw_data/16s"
    # create_manifest(manifest_filepath, data_folder, sra_runtable_filepath)
    patients_metadata_path = "/home/direnc/inputs/ibdmdb/hmp2_metadata.csv"
    abundance_table = "/home/direnc/results/ibdmdb/metadata/metadata.tsv"
    tax_tsv = "/home/direnc/results/ibdmdb/exported_2/taxonomyNoFilt_gtdb_barcodes_v4/taxonomy.tsv"
    boi_path_sum = "/home/direnc/results/ibdmdb/exported_2/summary_result_with_percentages.tsv"
    abundance_boi(abundance_table, tax_tsv, patients_metadata_path, samples_tsv_path, boi_path_sum)
    pass