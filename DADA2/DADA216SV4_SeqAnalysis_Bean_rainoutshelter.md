# Raw Sequences Analysis of V4 region of 16S rRNA  gene from common bean seed under rainoutshelter treatment (field experiment at MSU Agraonomy & UPRECH Farms)

## Analysis of 16S Miseq Data
Here we used the DADA2 denoising pipeline for pre-processing raw sequences data and UPARSE method for OTU clustering (Edgar 2013). Additional analyses were conducted using QIIME2 (qiime2-2021.4).

raw sequence data stored on HPCC:
/mnt/research/ShadeLab/Sequence/raw_sequence/Bean_rainoutshelter/

Moved/copy raw sequences to the working space:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/

Renaming the fastq files name to sample id, i.e.:
```
for f in *_01_*.gz; do mv -v "$f" "${f/AFB_20210429_01/101}"; done;

output:
‘AFB_20210429_01_S55_L001_R1_001.fastq.gz’ -> ‘101_S55_L001_R1_001.fastq.gz’
‘AFB_20210429_01_S55_L001_R2_001.fastq.gz’ -> ‘101_S55_L001_R2_001.fastq.gz’

# adding 'r1' or 'r2' as sequencing run replicate, i.e.:
for f in *.gz; do mv -v "$f" "${f/_S/_r2_S}"; done;

# sequences of r1 and r2 will be analyzed separately 
````

## All analysis results are stored on hpcc: 

Sequencing Run 1:

"/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/16S_Bean_rainoutshelter/20210601_16SV4_PE250_renamed/DADA2"

#########################################################################################################################################################

Sequencing Run 2:

"/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/16S_Bean_rainoutshelter/20210604_16SV4_PE250_renamed/DADA2"

# Part I: Quality control and denoising using DADA2 to make ASV

## 1) Making manifest file using Excel (tab delimited file, .tsv)

```
## example:

sample-id forward-absolute-filepath reverse-absolute-filepath
1001  /mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210601_16SV4_PE250_renamed/raw_reads/1001_r1_S5_L001_R1_001.fastq   /mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210601_16SV4_PE250_renamed/raw_reads/1001_r1_S5_L001_R2_001.fastq

## save the file as .tsv  file

Sequencing Run 1: 20210601_16SV4_manifest.tsv 
Sequencing Run 1: 20210604_16SV4_manifest.tsv

## transfer the manifest file to your working directory on msu hpcc

Sequencing Run 1:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210601_16SV4_PE250_renamed/DADA2
20210601_16SV4_manifest.tsv 

Sequencing Run 2:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210604_16SV4_PE250_renamed/DADA2
```

## 2) Load Qiime2
```
conda activate qiime2-2021.4
```

## 3) Convert manifest file to q2 artifacts
```
## It can take a while if you have many sample. Thus it is better to submit as a job on MSU HPCC. 


#### Sequencing Run 1 ####

## Path to the job file:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210601_16SV4_PE250_renamed/DADA2/importmanifest.sb

## command for importing manifest file
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path 20210601_16SV4_manifest.tsv \
--output-path q2/paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

## create a visualization file from 'paired-end-demux.qza' (optional)
qiime demux summarize \
--i-data paired-end-demux.qza \
--o-visualization q2/demux.qzv

#########################################################################################################################################################

#### Sequencing Run 2 ####

## Path to the job file:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210604_16SV4_PE250_renamed/DADA2/importmanifest.sb

## command for importing manifest file
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path 20210601_16SV4_manifest.tsv \
--output-path q2/paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

## create a visualization file from 'paired-end-demux.qza' (optional)

qiime demux summarize \
--i-data q2/paired-end-demux.qza \
--o-visualization q2/demux.qzv
```

## 4) Denoise reads with DADA2
```
## Prior to denoising  using DADA2, you need to investigate your Q-score distribution of your reads to determine the 'trunc-len-f' and 'trunc-len-r' values. The goal is to trim the reads/sequences by removing as much of the lower quality portions of the reads as possible and still leave enough overlap so that merging of the forward and reverse reads can be optimized.

##run the figaro activation and tutorial from John Quensen's website (http://john-quensen.com/tutorials/figaro/) and the reference can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7424690/

######## deactivate qiime2 before use figaro because figaro need miniconda3, while qiime2 environment need anaconda2 #############
conda deactivate

## activate figaro
cd
source ~/miniconda3/etc/profile.d/conda.sh
conda activate figaro

#### Sequencing Run 1 ####

## run command
mkdir figaro

python /mnt/home/bintarti/figaro/figaro/figaro.py -i /mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210601_16SV4_PE250_renamed/raw_reads/ -o /mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210601_16SV4_PE250_renamed/DADA2/figaro/ -f 1 -r 1 -a 253 -F illumina

## Based on the figaro result, the recommended forward truncation position is 135 and the recommended reverse truncation position is 140. After trimming and truncation, the expected number of errors in the forward and reverse read is 1. By using this parameters should result in merging 92.12% of the reads.

## deactivate figaro
conda deactivate

##  activate qiime2
conda activate /mnt/home/bintarti/anaconda2/envs/qiime2-2021.4

## It can take a while if you have many sample. Thus it is better to submit as a job on MSU HPCC. Path to the job file:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210601_16SV4_PE250_renamed/DADA2/denoise_dada2.sb

## command for denoising using DADA2
qiime dada2 denoise-paired \
--i-demultiplexed-seqs q2/paired-end-demux.qza \
--p-trunc-len-f 135 \
--p-trunc-len-r 140 \
--o-table q2/table.qza \
--o-representative-sequences q2/rep-seqs.qza \
--o-denoising-stats q2/denoising-stats.qza \
--p-n-threads 10

## create visualization files
qiime metadata tabulate --m-input-file q2/denoising-stats.qza --o-visualization q2/denoising-stats.qzv
qiime feature-table summarize --i-table q2/table.qza --o-visualization q2/table.qzv
qiime feature-table tabulate-seqs --i-data q2/rep-seqs.qza --o-visualization q2/rep-seqs.qzv

#########################################################################################################################################################

#### Sequencing Run 2 ####

## run command
mkdir figaro

python /mnt/home/bintarti/figaro/figaro/figaro.py -i /mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210604_16SV4_PE250_renamed/raw_reads/ -o /mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210604_16SV4_PE250_renamed/DADA2/figaro/ -f 1 -r 1 -a 253 -F illumina

## Based on the figaro result, the recommended forward truncation position is 139 and the recommended reverse truncation position is 136. After trimming and truncation, the expected number of errors in the forward and reverse read is 1. By using this parameters should result in merging 91.07% of the reads.

## deactivate figaro
conda deactivate

##  activate qiime2
conda activate /mnt/home/bintarti/anaconda2/envs/qiime2-2021.4

## It can take a while if you have many sample. Thus it is better to submit as a job on MSU HPCC. Path to the job file:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/20210604_16SV4_PE250_renamed/DADA2/denoise_dada2.sb

## command for denoising using DADA2
qiime dada2 denoise-paired \
--i-demultiplexed-seqs q2/paired-end-demux.qza \
--p-trunc-len-f 139 \
--p-trunc-len-r 136 \
--o-table q2/table.qza \
--o-representative-sequences q2/rep-seqs.qza \
--o-denoising-stats q2/denoising-stats.qza \
--p-n-threads 10

## create visualization files:

qiime metadata tabulate --m-input-file q2/denoising-stats.qza --o-visualization q2/denoising-stats.qzv

qiime feature-table summarize --i-table q2/table.qza --o-visualization q2/table.qzv

qiime feature-table tabulate-seqs --i-data q2/rep-seqs.qza --o-visualization q2/rep-seqs.qzv
```

## 5) Assign taxonomy using pre-trained (515F-806R) SILVA 138 reference database downloaded from the QIIME2 website: https://docs.qiime2.org/2021.4/data-resources/. Reference: Bokulich et al. 2018 (https://doi.org/10.1186/s40168-018-0470-z, http://doi.org/10.5281/zenodo.3891931). Classify-sklearn is a machine learning based classification method with a Naive Bayes classifier.
```
##  Download this file from the QIIME2 website:
Silva 138 99% OTUs from 515F/806R region of sequences (MD5: e05afad0fe87542704be96ff483824d4)
## You also can submit it as a job on the MSU HPCC because it may take a while to run


#### Sequencing Run 1 ####

qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/home/bintarti/qiime2/silva-138-99-515-806-nb-classifier.qza \
  --i-reads q2/rep-seqs.qza \
  --o-classification q2/taxonomy.qza \
  --p-confidence .8 \
  --p-n-jobs -10

## create visualization files

qiime metadata tabulate --m-input-file q2/taxonomy.qza --o-visualization q2/taxonomy.qzv

### output ###
taxonomy.qza
taxonomy.qzv

#########################################################################################################################################################

#### Sequencing Run 2 ####

qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/home/bintarti/qiime2/silva-138-99-515-806-nb-classifier.qza \
  --i-reads q2/rep-seqs.qza \
  --o-classification q2/taxonomy.qza \
  --p-confidence .8 \
  --p-n-jobs -10

## create visualization files

qiime metadata tabulate --m-input-file q2/taxonomy.qza --o-visualization q2/taxonomy.qzv

### output ###
taxonomy.qza
taxonomy.qzv
```

## 6) Filter non-bacteria/archcd from tables and sequences
```
#### Sequencing Run 1 ####

# from ASV table

qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table ASVtable-no-mitoch-no-chloro.qza

# from rep sequences

qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences rep-seqs-no-mitoch-no-chloro.qza

### output ###
ASVtable-no-mitoch-no-chloro.qza
rep-seqs-no-mitoch-no-chloro.qza

## create visualization files

qiime feature-table summarize --i-table q2/ASVtable-no-mitoch-no-chloro.qza --o-visualization q2/ASVtable-no-mitoch-no-chloro.qzv
qiime feature-table tabulate-seqs --i-data q2/rep-seqs-no-mitoch-no-chloro.qza --o-visualization q2/rep-seqs-no-mitoch-no-chloro.qzv

#########################################################################################################################################################

#### Sequencing Run 2 ####

# from ASV table

qiime taxa filter-table \
  --i-table q2/table.qza \
  --i-taxonomy q2/taxonomy.qza \
  --p-exclude mitochondria,chloroplast,unassigned \
  --o-filtered-table q2/ASVtable-no-mitoch-no-chloro.qza

# from rep sequences

qiime taxa filter-seqs \
  --i-sequences q2/rep-seqs.qza \
  --i-taxonomy q2/taxonomy.qza \
  --p-exclude mitochondria,chloroplast,unassigned \
  --o-filtered-sequences q2/rep-seqs-no-mitoch-no-chloro.qza

### output ###
ASVtable-no-mitoch-no-chloro.qza
rep-seqs-no-mitoch-no-chloro.qza

## create visualization files

qiime feature-table summarize --i-table q2/ASVtable-no-mitoch-no-chloro.qza --o-visualization q2/ASVtable-no-mitoch-no-chloro.qzv

qiime feature-table tabulate-seqs --i-data q2/rep-seqs-no-mitoch-no-chloro.qza --o-visualization q2/rep-seqs-no-mitoch-no-chloro.qzv
```

## 7) Align the representative sequences and construct a phylogenetic tree from the alignment using MAFFT (Multiple Alignment using Fast Fourier Transform)
```
#### Sequencing Run 1 ####

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences q2/rep-seqs-no-mitoch-no-chloro.qza \
--o-alignment q2/aligned-rep-seqs-no-mitoch-no-chloro.qza \
--o-masked-alignment q2/masked-aligned-rep-seqs-no-mitoch-no-chloro.qza \
--o-tree q2/unrooted-tree.qza \
--o-rooted-tree q2/rooted-tree.qza

### output ###
aligned-rep-seqs-no-mitoch-no-chloro.qza
masked-aligned-rep-seqs-no-mitoch-no-chloro.qza
unrooted-tree.qza
rooted-tree.qza

#########################################################################################################################################################

#### Sequencing Run 2 ####

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences q2/rep-seqs-no-mitoch-no-chloro.qza \
--o-alignment q2/aligned-rep-seqs-no-mitoch-no-chloro.qza \
--o-masked-alignment q2/masked-aligned-rep-seqs-no-mitoch-no-chloro.qza \
--o-tree q2/unrooted-tree.qza \
--o-rooted-tree q2/rooted-tree.qza

### output ###
aligned-rep-seqs-no-mitoch-no-chloro.qza
masked-aligned-rep-seqs-no-mitoch-no-chloro.qza
unrooted-tree.qza
rooted-tree.qza
```
## 8) Export table to biom
```
mkdir export_file

cp table.qza export_file/
cp ASVtable-no-mitoch-no-chloro.qza export_file/
cp taxonomy.qza export_file/ 

qiime tools export \
    --input-path export_file/table.qza \
    --output-path Biom/

mv Biom/feature-table.biom Biom/ASVtable_unfiltered.biom

qiime tools export \
    --input-path export_file/ASVtable-no-mitoch-no-chloro.qza \
    --output-path Biom/

mv Biom/feature-table.biom Biom/ASVtable-no-mitoch-no-chloro.biom

biom convert -i Biom/ASVtable_unfiltered.biom -o Biom/ASVtable_unfiltered.tsv --to-tsv

biom convert -i Biom/ASVtable-no-mitoch-no-chloro.biom -o Biom/ASVtable-no-mitoch-no-chloro.tsv --to-tsv

qiime tools export --input-path export_file/taxonomy.qza --output-path exported_tax

cp exported_tax/taxonomy.tsv Biom/biom-taxonomy.tsv
```

## 9) Add taxonomy to the ASV table
```
biom add-metadata -i Biom/ASVtable_unfiltered.biom -o Biom/ASVtable_unfiltered_tax.biom --observation-metadata-fp=Biom/biom-taxonomy.tsv --sc-separated=taxonomy --observation-header=OTUID,taxonomy,confidence

biom add-metadata -i Biom/ASVtable-no-mitoch-no-chloro.biom -o Biom/ASVtable-no-mitoch-no-chloro_tax.biom --observation-metadata-fp=Biom/biom-taxonomy.tsv --sc-separated=taxonomy --observation-header=OTUID,taxonomy,confidence

biom convert -i Biom/ASVtable_unfiltered_tax.biom -o TSV/ASVtable_unfiltered_tax.tsv --header-key taxonomy --to-tsv

biom convert -i Biom/ASVtable-no-mitoch-no-chloro_tax.biom -o TSV/ASVtable-no-mitoch-no-chloro_tax.tsv --header-key taxonomy --to-tsv
```









