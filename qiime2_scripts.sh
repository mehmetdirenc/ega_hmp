#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=63
#SBATCH --job-name=16s_qiime_v4
#SBATCH --time=3-00:05

#manifest=$1
#result_path_without_slash_at_the_end=$2
#metadata=$3

source /mnt/lustre/home/mager/magmu818/anaconda3/etc/profile.d/conda.sh
conda activate /mnt/lustre/groups/mager/magmu818/.conda/envs/qiime_updated


cd /mnt/lustre/home/mager/magmu818/slurm/mahana/16s_v4/test_moving


qiime tools import \
  --type 'EMPSingleEndSequences' \
  --input-path emp-single-end-sequences \
  --output-path emp-single-end-sequences.qza

qiime tools peek emp-single-end-sequences.qza

qiime demux emp-single \
  --i-seqs emp-single-end-sequences.qza \
  --m-barcodes-file sample-metadata.tsv \
  --m-barcodes-column barcode-sequence \
  --o-per-sample-sequences demux.qza \
  --o-error-correction-details demux-details.qza


qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 120 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza

qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv


qiime metadata tabulate \
  --m-input-file demux-filter-stats.qza \
  --o-visualization demux-filter-stats.qzv


qiime feature-table summarize \
  --i-table table.qza \
  --m-sample-metadata-file sample-metadata.tsv \
  --o-visualization table.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --output-dir phylogeny-align-to-tree-mafft-fasttree

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny phylogeny-align-to-tree-mafft-fasttree/rooted_tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1103 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir diversity-core-metrics-phylogenetic

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity-core-metrics-phylogenetic/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity-core-metrics-phylogenetic/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity-core-metrics-phylogenetic/shannon_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization shannon-group-significance.qzv


qiime diversity beta-group-significance \
  --i-distance-matrix diversity-core-metrics-phylogenetic/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column body-site \
  --p-pairwise \
  --o-visualization unweighted-unifrac-body-site-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix diversity-core-metrics-phylogenetic/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column subject \
  --p-pairwise \
  --o-visualization unweighted-unifrac-subject-group-significance.qzv

qiime emperor plot \
  --i-pcoa diversity-core-metrics-phylogenetic/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-custom-axes days-since-experiment-start \
  --o-visualization unweighted-unifrac-emperor-days-since-experiment-start.qzv

qiime emperor plot \
  --i-pcoa diversity-core-metrics-phylogenetic/bray_curtis_pcoa_results.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-custom-axes days-since-experiment-start \
  --o-visualization bray-curtis-emperor-days-since-experiment-start.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/lustre/home/mager/magmu818/datasets/public_databases/qiime/gtdb/with_barcodes/gtdb_220_with_barcodes_classifier_v4.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv


qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-where '[body-site]="gut"' \
  --o-filtered-table gut-table.qza


qiime composition ancombc \
  --i-table gut-table.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-formula subject \
  --o-differentials ancombc-subject.qza
qiime composition da-barplot \
  --i-data ancombc-subject.qza \
  --p-significance-threshold 0.001 \
  --o-visualization da-barplot-subject.qzv


qiime taxa collapse \
  --i-table gut-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table gut-table-l6.qza
qiime composition ancombc \
  --i-table gut-table-l6.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-formula subject \
  --o-differentials l6-ancombc-subject.qza
qiime composition da-barplot \
  --i-data l6-ancombc-subject.qza \
  --p-significance-threshold 0.001 \
  --o-visualization l6-da-barplot-subject.qzv






#
#qiime tools import \
#  --type "SampleData[PairedEndSequencesWithQuality]" \
#  --input-format PairedEndFastqManifestPhred33V2 \
#  --input-path $1 \
#  --output-path $2/imported_seqs.qza

#qiime demux summarize \
#        --i-data $2/imported_seqs.qza \
#        --o-visualization $2/imported_seqs.qzv
#
#
#qiime dada2 denoise-paired \
#        --i-demultiplexed-seqs $2/imported_seqs.qza \
#        --p-trim-left-f 0 \
#        --p-trim-left-r 0 \
#        --p-trunc-len-f 230 \
#        --p-trunc-len-r 220 \
#        --o-table $2/tableNoFilt.qza \
#        --o-representative-sequences $2/repseqsNoFilt.qza \
#        --o-denoising-stats $2/denoising-statsNoFilt.qza \
#        --p-n-threads 63
#
#qiime metadata tabulate \
#        --m-input-file $2/repseqsNoFilt.qza \
#        --o-visualization $2/repseqsNoFilt.qzv
#
#qiime metadata tabulate \
#        --m-input-file $2/denoising-statsNoFilt.qza \
#        --o-visualization $2/denoising-statsNoFilt.qzv
#
#qiime metadata tabulate \
#        --m-input-file $2/tableNoFilt.qza \
#        --o-visualization $2/tableNoFilt.qzv
#
#
#
#
#qiime feature-classifier classify-sklearn \
#        --i-reads $2/repseqsNoFilt.qza \
#        --i-classifier /mnt/lustre/home/mager/magmu818/datasets/public_databases/qiime/gtdb/with_barcodes/gtdb_220_with_barcodes_classifier_v4.qza \
#        --o-classification $2/taxonomyNoFilt_gtdb_barcodes_v4.qza --p-n-jobs 63
#
#
#qiime metadata tabulate \
#        --m-input-file $2/taxonomyNoFilt_gtdb_barcodes_v4.qza \
#        --o-visualization $2/taxonomyNoFilt_gtdb_barcodes_v4.qzv
#
#qiime taxa barplot \
#        --i-table $2/tableNoFilt.qza \
#        --i-taxonomy $2/taxonomyNoFilt_gtdb_barcodes_v4.qza \
#        --o-visualization $2/taxa-bar-plots_gtdb_barcodes_v4.qzv
#
#
#qiime tools export \
#        --input-path $2/taxonomyNoFilt_gtdb_barcodes_v4.qza \
#        --output-path $2/exported_2/taxonomyNoFilt_gtdb_barcodes_v4