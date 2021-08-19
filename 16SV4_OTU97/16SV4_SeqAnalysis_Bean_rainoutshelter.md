# Raw Sequences Analysis of V4 region of 16S rRNA  gene from common bean seed under rainoutshelter treatment (field experiment at MSU Agraonomy & UPRECH Farms)

## Analysis of 16S Miseq Data
Here we used the USEARCH pipeline (v10.0.240_i86linux64) for pre-processing raw sequences data and UPARSE method for OTU clustering (Edgar 2013). Additional analyses were conducted using QIIME1

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

"/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/16S_Bean_rainoutshelter/20210601_16SV4_PE250_renamed"

Sequencing Run 2:

"/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_rainoutshelter/16S_Bean_rainoutshelter/20210604_16SV4_PE250_renamed"

# Part I: Clustering

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

09:46 701Mb   100.0% 91.5% merged
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

03:11 627Mb   100.0% Filtering, 99.4% passed
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
Min size 1, median 1, max 3081804, avg 26.98

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

## 6) Precluster Sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_fast nosigs_uniques_filtered_merged.fastq -centroids_fastq denoised_nosigs_uniques_filtered_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size

### output ###
127936 seqs (tot.size 7983979), 127936 uniques, 0 singletons (0.0%)
Seqs  127936 (127.9k)
  Clusters  2281
  Max size  5731343 (5.7M)
  Avg size  3500.2
  Min size  2
Singletons  0, 0.0% of seqs, 0.0% of clusters
   Max mem  739Mb
      Time  4.00s
Throughput  32.0k seqs/sec.


#### for r2: ####


```

## 7) Closed Reference-based OTU Picking Using SILVA_132 Database
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global denoised_nosigs_uniques_filtered_merged.fastq -id 0.97 -db /mnt/home/bintarti/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -strand plus -uc ref_seqs.uc -dbmatched SILVA97_closed_reference.fasta -notmatchedfq failed_closed.fq

### output ###
100.0% Searching, 40.2% matched

Closed reference OTU picking:
Pick OTUs based on the Silva database in the home directory.
Produce some output files - ref_seqs.uc (pre-clustered), SILVA97_closed_reference.fasta will be the matched ones, and failed_closed.fq will be used in de novo OTU picking
```

## 8) De novo OTU picking
```
# sort by size
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize failed_closed.fq -fastaout sorted_failed_closed.fq

# cluster de novo
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus sorted_failed_closed.fq -minsize 2 -otus denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up

### output ###
100.0% 122 OTUs, 450 chimeras




```

## 9) Combine the Rep Sets Between De novo and SILVA Reference-based OTU Picking
```
cat SILVA97_closed_reference.fasta denovo_otus.fasta > FULL_REP_SET.fna
```

## 10) Map 'FULL_REP_SET.fna' Back to Pre-dereplicated Sequences and Make OTU Tables
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -usearch_global mergedfastq/merged.fq -db FULL_REP_SET.fna -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom

### output ###
8193798 / 8220965 mapped to OTUs (99.7%)
OTU_table.txt
OTU_jsn.biom
```

# Part II: Switch to QIIME 1.9.1

## 1) Convert OTU_table.txt to OTU_table.from_txt_json.biom
```
biom convert -i OTU_table.txt -o OTU_table.biom --table-type="OTU table" --to-json

### output ###
OTU_table.biom
```
## 2) Align sequences to the SILVA_132_QIIME_release database with PyNAST 
```
align_seqs.py -i FULL_REP_SET.fna -o alignment -t /mnt/home/bintarti/SILVA_132_QIIME_release/core_alignment/80_core_alignment.fna

### output ###
alignment/FULL_REP_SET_aligned.fasta
alignment/FULL_REP_SET_failures.fasta
alignment/FULL_REP_SET_log.txt
```
## 3) Filter failed alignment from OTU table
```
#Discard all OTUs listed in FULL_REP_SET_failures.fasta from OTU table

filter_otus_from_otu_table.py -i OTU_table.biom -o OTU_filteredfailedalignments.biom -e alignment/FULL_REP_SET_failures.fasta

#from FULL_REP_SET.fna file

filter_fasta.py -f FULL_REP_SET.fna -o FULL_REP_SET_filteredfailedalignments.fa -a alignment/FULL_REP_SET_aligned.fasta

### output ###
OTU_filteredfailedalignments.biom
FULL_REP_SET_filteredfailedalignments.fa
```
## 4) Assign taxonomy to the SILVA_132_QIIME_release database with UCLUST
```
assign_taxonomy.py -i FULL_REP_SET_filteredfailedalignments.fa -o taxonomy -r /mnt/home/bintarti/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -t /mnt/home/bintarti/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt
#-r reference -> path to silva db 
# -t taxonomy

### output ###
taxonomy/FULL_REP_SET_filteredfailedalignments_tax_assignments.txt
taxonomy/FULL_REP_SET_filteredfailedalignments_tax_assignments.log
```
## 5) Add taxonomy to the OTU table
```
echo "#OTUID"$'\t'"taxonomy"$'\t'"confidence" > templine.txt

cat templine.txt taxonomy/FULL_REP_SET_filteredfailedalignments_tax_assignments.txt >> taxonomy/FULL_REP_SET_filteredfailedalignments_tax_assignments_header.txt

biom add-metadata -i OTU_filteredfailedalignments.biom -o OTU_table_tax.biom --observation-metadata-fp=taxonomy/FULL_REP_SET_filteredfailedalignments_tax_assignments_header.txt  --sc-separated=taxonomy --observation-header=OTUID,taxonomy

### output ###
OTU_table_tax.biom
```
## 6) Filter non-bacteria/archaea
```
filter_taxa_from_otu_table.py -i OTU_table_tax.biom -o OTU_table_tax_filt.biom -n D_4__Mitochondria,D_3__Chloroplast,Chlorophyta,Unassigned

#remove same Mito and Chloro sequences from RepSeqs file
filter_fasta.py -f FULL_REP_SET_filteredfailedalignments.fa -o FULL_REP_SET_filteredfailedalignments_rmCM.fa -b OTU_table_tax_filt.biom 

#summarize OTU table

#1.original

biom summarize-table -i OTU_table_tax.biom -o OTU_table_tax_sum.txt

#2.filtered

biom summarize-table -i OTU_table_tax_filt.biom -o OTU_table_tax_filt_sum.txt

#optional sanity check:  count seqs in new fasta, and check that it has fewer than original

#orginal

count_seqs.py -i FULL_REP_SET_filteredfailedalignments.fa

#filtered

count_seqs.py -i FULL_REP_SET_filteredfailedalignments_rmCM.fa

### output ###
OTU_table_tax_filt.biom
FULL_REP_SET_filteredfailedalignments_rmCM.fa
OTU_table_tax_sum.txt
OTU_table_tax_filt_sum.txt
```
## 7) Rarefaction – will not be performed in this study. We will use reads normalization using CSS method from the MetagenomeSeq package on R.
```
single_rarefaction.py -d 11137 -o single_rare.biom -i OTU_table_tax_filt.biom

biom summarize-table -i single_rare.biom -o single_rare_sum.txt

### output ###
single_rare.biom
single_rare_sum.txt
```
## 8) Summarize global taxonomic data
```
summarize_taxa.py -i OTU_table_tax_filt.biom -o taxa_sum

### output ###
taxa_sum/
```
## 9) Make phylogeny with FastTree
```
#First, clean alignment by omitting highly variable regions before tree building - will make tree building more efficient

filter_alignment.py -i alignment/FULL_REP_SET_aligned.fasta -o alignment/filtered_alignment

#make phylogeny and root tree

make_phylogeny.py -i alignment/filtered_alignment/FULL_REP_SET_aligned_pfiltered.fasta -o rep_set.tre -r tree_method_default

### output ###
alignment/filtered_alignment/FULL_REP_SET_aligned_pfiltered.fasta
rep_set.tre
```
## 10) Convert and add taxonomy
```
# 1. filtered OTU table
biom convert -i OTU_table_tax_filt.biom -o OTU_table_tax_filt.txt --header-key taxonomy --to-tsv

### output ###
OTU_table_tax_filt.txt

# 2. original OTU table
biom convert -i OTU_table_tax.biom -o OTU_table_tax.txt --header-key taxonomy --to-tsv


```
# Part III: Switch to R

## 10) Removing DNA contaminant from the OTU_table_tax_filt.txt file was conducted on R using microDecon package (McKnight et al. 2019)
```
### output ###
removing 18 contaminant OTUs = 255 OTUs
```
## 11) Removing OTUs that do not present in sample
```
### output ###
removing 44 zero OTUs = 211 OTUs
```
## 12) Reads normalization using CSS method from the metagenomeSeq package (Paulson et al. 2013)
```
### output ###
file = otu_norm.txt 
# move the file from the local computer to the hpcc (/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_variability/16S_Bean_variability/test)
```
# Part IV: Switch to QIIME on hpcc

## 13) Convert "otu_norm.txt" to "otu_norm.biom"
```
biom convert -i otu_norm.txt -o otu_norm.biom --table-type="OTU table" --to-json

### output ###
otu_norm.biom
```
## 14) Prune the "rep_set.tre" based on the set of tip names from the "otu_norm.txt" (tips_to_keep.txt)
```
filter_tree.py -i rep_set.tre -t tips_to_keep.txt -o pruned.tre
```
## 15) Calculate PD_whole_tree (alpha diversity)
```
alpha_diversity.py -m PD_whole_tree -i otu_norm.biom -o PD -t pruned.tre
```








# Part II: Switch to QIIME2 (qiime2-2021.4)

## 1) Load Qiime2
```
conda activate qiime2-2021.4
```
## 2) Convert rep-seqs to q2 artifacts
```
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

## 5) Assign taxonomy using pre-trained (515F-806R) SILVA 138 referance database downloaded from qiime2 website: https://docs.qiime2.org/2021.4/data-resources/. Reference: Bokulich et al. 2018 (https://doi.org/10.1186/s40168-018-0470-z, http://doi.org/10.5281/zenodo.3891931). Classify-sklearn is a machine learning based  classification method with a Naive Bayes classifier
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























