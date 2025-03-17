#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=63
#SBATCH --job-name=16s_qiime_v4
#SBATCH --time=3-00:05

manifest=$1
result_path_without_slash_at_the_end=$2
#metadata=$3

source /mnt/lustre/home/mager/magmu818/anaconda3/etc/profile.d/conda.sh
conda activate /mnt/lustre/groups/mager/magmu818/.conda/envs/qiime_updated

qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path $1 \
  --output-path $2/imported_seqs.qza

qiime demux summarize \
        --i-data $2/imported_seqs.qza \
        --o-visualization $2/imported_seqs.qzv


qiime dada2 denoise-paired \
        --i-demultiplexed-seqs $2/imported_seqs.qza \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-trunc-len-f 230 \
        --p-trunc-len-r 220 \
        --o-table $2/tableNoFilt.qza \
        --o-representative-sequences $2/repseqsNoFilt.qza \
        --o-denoising-stats $2/denoising-statsNoFilt.qza \
        --p-n-threads 63

qiime metadata tabulate \
        --m-input-file $2/repseqsNoFilt.qza \
        --o-visualization $2/repseqsNoFilt.qzv

qiime metadata tabulate \
        --m-input-file $2/denoising-statsNoFilt.qza \
        --o-visualization $2/denoising-statsNoFilt.qzv

qiime metadata tabulate \
        --m-input-file $2/tableNoFilt.qza \
        --o-visualization $2/tableNoFilt.qzv




qiime feature-classifier classify-sklearn \
        --i-reads $2/repseqsNoFilt.qza \
        --i-classifier /mnt/lustre/home/mager/magmu818/datasets/public_databases/qiime/gtdb/with_barcodes/gtdb_220_with_barcodes_classifier_v4.qza \
        --o-classification $2/taxonomyNoFilt_gtdb_barcodes_v4.qza --p-n-jobs 63


qiime metadata tabulate \
        --m-input-file $2/taxonomyNoFilt_gtdb_barcodes_v4.qza \
        --o-visualization $2/taxonomyNoFilt_gtdb_barcodes_v4.qzv

qiime taxa barplot \
        --i-table $2/tableNoFilt.qza \
        --i-taxonomy $2/taxonomyNoFilt_gtdb_barcodes_v4.qza \
        --o-visualization $2/taxa-bar-plots_gtdb_barcodes_v4.qzv


qiime tools export \
        --input-path $2/taxonomyNoFilt_gtdb_barcodes_v4.qza \
        --output-path $2/exported_2/taxonomyNoFilt_gtdb_barcodes_v4
q