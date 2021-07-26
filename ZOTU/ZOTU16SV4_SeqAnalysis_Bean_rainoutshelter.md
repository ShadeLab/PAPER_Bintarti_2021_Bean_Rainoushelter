# Raw Sequences Analysis of V4 region of 16S rRNA  gene from common bean seed under rainoutshelter treatment (field experiment at MSU Agraonomy & UPRECH Farms)

## Analysis of 16S Miseq Data
Here we used the USEARCH pipeline (v10.0.240_i86linux64) for pre-processing raw sequences data and UPARSE method for OTU clustering (Edgar 2013). Additional analyses were conducted using QIIME2

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
```

## All analysis results are stored on hpcc: 

Sequencing Run 1:

"/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/16S_Bean_rainoutshelter/20210601_16SV4_PE250_renamed/ZOTU"

Sequencing Run 2:

"/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/16S_Bean_rainoutshelter/20210604_16SV4_PE250_renamed/ZOTU"

# Part I: Denoising using UNoise 3 to make zOTU

usearch v10.0.240_i86linux64, 16.3Gb RAM, 4 cores
(C) Copyright 2013-17 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

## 1) Merge Paired End Reads
```
# decompress the reads
gunzip *.gz

# make directory called "mergedfastq" in 
mkdir mergedfastq

# merge paired end reads
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs raw_reads/*R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -fastq_merge_maxee 1.0 -tabbedout mergedfastq/merged.report.txt -alnout mergedfastq/merged_aln.txt -log mergedfastq/merge_log.txt


# -fastq_maxdiffs 10: Allow 10 max differences in overlap region

### Output ###

#### for r1: ####

100.0% 93.7% merged
Totals:
   8776805  Pairs (8.8M)
   8220965  Merged (8.2M, 93.67%)
   5563200  Alignments with zero diffs (63.39%)
    544653  Too many diffs (> 10) (6.21%)
     11187  No alignment found (0.13%)
         0  Alignment too short (< 16) (0.00%)
         0  Exp.errs. too high (max=1.0) (0.00%)
      4887  Staggered pairs (0.06%) merged & trimmed
    247.33  Mean alignment length
    252.64  Mean merged length
      0.39  Mean fwd expected errors
      0.52  Mean rev expected errors
      0.07  Mean merged expected errors

#### for r2: ####

04:55 702Mb   100.0% 91.5% merged
Totals:
   7240182  Pairs (7.2M)
   6625368  Merged (6.6M, 91.51%)
   4696419  Alignments with zero diffs (64.87%)
    600148  Too many diffs (> 10) (8.29%)
         0  Fwd tails Q <= 2 trimmed (0.00%)
         5  Rev tails Q <= 2 trimmed (0.00%)
     14666  No alignment found (0.20%)
         0  Alignment too short (< 16) (0.00%)
         0  Exp.errs. too high (max=1.0) (0.00%)
      4325  Staggered pairs (0.06%) merged & trimmed
    247.36  Mean alignment length
    252.61  Mean merged length
      0.32  Mean fwd expected errors
      0.52  Mean rev expected errors
      0.06  Mean merged expected errors

```
## 2) Check Sequence Quality of Merged Seqs
```
mkdir fastq_info

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 mergedfastq/merged.fq -output fastq_info/eestats.txt

### output ###

#### for r1: ####

8220965 reads, max len 475, avg 252.6

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    8181754( 99.5%)    8218924(100.0%)    8220841(100.0%)
   100    8133989( 98.9%)    8212626( 99.9%)    8220701(100.0%)
   150    8086760( 98.4%)    8203159( 99.8%)    8219884(100.0%)
   200    8031988( 97.7%)    8191242( 99.6%)    8219187(100.0%)
   250    7942471( 96.6%)    8158227( 99.2%)    8212880( 99.9%)
   300        317(  0.0%)        511(  0.0%)        643(  0.0%)
   350        257(  0.0%)        432(  0.0%)        585(  0.0%)
   400        209(  0.0%)        379(  0.0%)        550(  0.0%)
   450          0(  0.0%)          0(  0.0%)          3(  0.0%)

#### for r2: ####

6625368 reads, max len 480, avg 252.6

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    6599866( 99.6%)    6623869(100.0%)    6625173(100.0%)
   100    6568195( 99.1%)    6619644( 99.9%)    6624956(100.0%)
   150    6529400( 98.6%)    6611081( 99.8%)    6624095(100.0%)
   200    6499366( 98.1%)    6604114( 99.7%)    6623503(100.0%)
   250    6457336( 97.5%)    6587906( 99.4%)    6618743( 99.9%)
   300        220(  0.0%)        373(  0.0%)        541(  0.0%)
   350        165(  0.0%)        316(  0.0%)        456(  0.0%)
   400        149(  0.0%)        279(  0.0%)        399(  0.0%)
```
## 3) Filter and Truncate the Merged Seqs
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter mergedfastq/merged.fq -fastq_maxee 1 -fastq_trunclen 250 -fastqout filtered_merged.fq -log filter_log.txt

### output ###

#### for r1: ####

100.0% Filtering, 99.2% passed
   8220965  Reads (8.2M)                    
      4159  Discarded reads length < 250
     58579  Discarded reads with expected errs > 1.00
   8158227  Filtered reads (8.2M, 99.2%)


#### for r2: ####

00:41 627Mb   100.0% Filtering, 99.4% passed
   6625368  Reads (6.6M)                    
      3562  Discarded reads length < 250
     33900  Discarded reads with expected errs > 1.00
   6587906  Filtered reads (6.6M, 99.4%)
```
## 4) Dereplicate Sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques filtered_merged.fq -fastqout uniques_filtered_merged.fastq -sizeout -log uniques_log.txt

### output ###

#### for r1: ####

8158227 seqs, 302184 uniques, 174248 singletons (57.7%)
00:32 7.1Gb  Min size 1, median 1, max 3845324, avg 27.00

#### for r2: ####

6587906 seqs, 244200 uniques, 136716 singletons (56.0%)
00:26 5.8Gb  Min size 1, median 1, max 3081804, avg 26.98

```
## 5) Remove Singeltons 
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize uniques_filtered_merged.fastq -fastqout nosigs_uniques_filtered_merged.fastq -minsize 2

### output ###

#### for r1: ####

Sorting 127936 sequences

#### for r2: ####

Sorting 107484 sequences
```

## 6) Predicting biological sequences = ZOTUs
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -unoise3 nosigs_uniques_filtered_merged.fastq -zotus ZOTU/zotus.fa -tabbedout ZOTU/zotus_report.txt -log ZOTU/unoise_log.txt

### output ###

#### for r1: ####

682 amplicons, 2370441 bad (size >= 8)
00:01 117Mb   100.0% 122 good, 560 chimeras 

#### for r2: ####

636 amplicons, 1869152 bad (size >= 8)
00:01 104Mb   100.0% 115 good, 521 chimeras  
```

## 7) Capitalize All Zotu Names To ZOTU (Necessary Downstream)
```
sed -i 's/Zotu/ZOTU/g' ZOTU/zotus.fa
```

## 8) Mapping reads to ZOTU
```
# It take sometime. Thus, it is better to submit it as a job on the MSU HPCC

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -otutab mergedfastq/merged.fq -zotus ZOTU/zotus.fa -uc ZOTU/ZOTU_map.uc -otutabout ZOTU/unoise3_ZOTU_table.txt -biomout ZOTU/unoise3_ZOTU_table.biom -notmatchedfq ZOTU/ZOTU_unmapped.fq -log ZOTU/zotutab_log.txt

### output ###
zotus.fa
unoise3_ZOTU_table.txt
unoise3_ZOTU_table.biom
ZOTU_unmapped.fq
ZOTU_map.uc
zotutab_log.txt
```

# Part II: Switch to QIIME2 (qiime2-2021.4)

## 1) Load Qiime2
```
conda activate qiime2-2021.4
```
## 2) Convert rep-seqs to q2 artifacts
```
mkdir q2

qiime tools import \
--input-path zotus.fa \
--output-path q2/rep-seqs.qza \
--type "FeatureData[Sequence]"

### output ###
rep-seqs.qza
```

## 3) Convert unoise3_ZOTU_table.txt to unoise3_ZOTU_table_hdf5.biom (BIOM v21 format)
```
biom convert  --input-fp unoise3_ZOTU_table.txt  --output-fp unoise3_ZOTU_table_hdf5.biom  --to-hdf5  --table-type="OTU table"

### output ###
unoise3_ZOTU_table_hdf5.biom
```

## 4) Convert ZOTU table to q2 artifacts
```
qiime tools import \
  --input-path unoise3_ZOTU_table_hdf5.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path q2/ZOTU_table.qza

### output ###
ZOTU_table.qza
```

## 5) Assign taxonomy using pre-trained (515F-806R) SILVA 138 reference database downloaded from Qiime2 website: https://docs.qiime2.org/2021.4/data-resources/. Reference: Bokulich et al. 2018 (https://doi.org/10.1186/s40168-018-0470-z, http://doi.org/10.5281/zenodo.3891931). Classify-sklearn is a machine learning based classification method with a Naive Bayes classifier.
```
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/home/bintarti/qiime2/silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza \
  --p-confidence .8 \
  --p-n-jobs -10

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

### output ###
taxonomy.qza
taxonomy.qzv
```

## 6) Filter non-bacteria/archaea from tables and sequences
```
# from ZOTU table

qiime taxa filter-table \
  --i-table ZOTU_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table ZOTUtable-no-mitoch-no-chloro.qza

# from rep sequences

qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences rep-seqs-no-mitoch-no-chloro.qza

### output ###
ZOTUtable-no-mitoch-no-chloro.qza
rep-seqs-no-mitoch-no-chloro.qza

## create visualization files

qiime feature-table summarize --i-table q2/ZOTUtable-no-mitoch-no-chloro.qza --o-visualization q2/ZOTUtable-no-mitoch-no-chloro.qzv
qiime feature-table tabulate-seqs --i-data q2/rep-seqs-no-mitoch-no-chloro.qza --o-visualization q2/rep-seqs-no-mitoch-no-chloro.qzv
```

## 7) Align the representative sequences and construct a phylogenetic tree from the alignment using MAFFT (Multiple Alignment using Fast Fourier Transform)
```
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs-no-mitoch-no-chloro.qza \
--o-alignment aligned-rep-seqs-no-mitoch-no-chloro.qza \
--o-masked-alignment masked-aligned-rep-seqs-no-mitoch-no-chloro.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

### output ###
aligned-rep-seqs-no-mitoch-no-chloro.qza
masked-aligned-rep-seqs-no-mitoch-no-chloro.qza
unrooted-tree.qza
rooted-tree.qza
```

## 8) Export table to biom

```
mkdir export_file

cp ZOTU_table.qza export_file/
cp taxonomy.qza export_file/ 

qiime tools export \
    --input-path ZOTU_table.qza \
    --output-path Biom/

mv feature-table.biom ZOTU_table.biom

qiime tools export \
    --input-path ZOTUtable-no-mitoch-no-chloro.qza \
    --output-path Biom/

mv feature-table.biom ZOTUtable-no-mitoch-no-chloro.biom

biom convert -i Biom/ZOTU_table.biom -o Biom/ZOTU_table.tsv --to-tsv
biom convert -i Biom/ZOTUtable-no-mitoch-no-chloro.biom -o Biom/ZOTUtable-no-mitoch-no-chloro.tsv --to-tsv


qiime tools export --input-path taxonomy.qza --output-path exported_tax

cp exported_tax/taxonomy.tsv biom-taxonomy.tsv
```

## 9) Add taxonomy to the ZOTU table
```
biom add-metadata -i Biom/ZOTU_table.biom -o ZOTU_table_tax.biom --observation-metadata-fp biom-taxonomy.tsv --sc-separated taxonomy

biom add-metadata -i Biom/ZOTU_table.biom -o ZOTU_table_tax.biom --observation-metadata-fp=biom-taxonomy.tsv --sc-separated=taxonomy --observation-header=OTUID,taxonomy,confidence


biom add-metadata -i Biom/ZOTUtable-no-mitoch-no-chloro.biom -o ZOTUtable-no-mitoch-no-chloro_tax.biom --observation-metadata-fp biom-taxonomy.tsv --sc-separated taxonomy

biom convert -i ZOTU_table_tax.biom -o TSV/ZOTU_table_tax.tsv --header-key taxonomy --to-tsv
biom convert -i ZOTUtable-no-mitoch-no-chloro_tax.biom -o TSV/ZOTUtable-no-mitoch-no-chloro_tax.tsv --header-key taxonomy --to-tsv
```


















