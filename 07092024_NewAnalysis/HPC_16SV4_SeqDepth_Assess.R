##################################################################################################################################
# Bean seed microbiome drought 2024 (rain-shelter field experiment): Sequencing Depth Assessment & Decontamination  #########################
##################################################################################################################################

# Date: July 19th 2024
# By : A. Fina Bintarti
# CNRS, Lyon, FR

# INSTALL PACKAGES
library(multcomp)
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)


######################################################################################################################################
# SEQUENCING DEPTH ASSESSMENT: PLANT CONTAMINANT PROPORTION
######################################################################################################################################

# SET THE WORKING DIRECTORY
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07')
wd <- print(getwd())

# READ PROPORTION OF CHLOROPLAST AND MITOCHONDRIA

# read the unfiltered ASV table 
asv.tax.unfil <- read.table(file = 'ASVtable_unfiltered_tax_07.tsv', sep = '\t', header = TRUE, row.names = "ASV.ID", check.names = FALSE)
dim(asv.tax.unfil) #2051
# separate taxonomy from the ASV table
tax.unfil <- asv.tax.unfil[,'taxonomy',drop = FALSE]
head(tax.unfil)
#write.csv(tax.unfil, file = "taxonomy_unfil.csv") Edit the taxonomy format locally
#read taxonomy
tax.unfil.ed = read.csv("taxonomy_unfil_edited.csv", header=T)
head(tax.unfil.ed)
dim(tax.unfil.ed) # 2051
# ASV table
colnames(asv.tax.unfil)
asv.unfil <- asv.tax.unfil[,1:84] # unselect the taxonomy
colnames(asv.unfil)
head(asv.unfil) #2051 84
asv.unfil <- rownames_to_column(asv.unfil, var = "ASV.ID")
# merge taxonomy to the unfiltered otu table
asv.tax.unfil.ed <- merge(asv.unfil, tax.unfil.ed, by="ASV.ID")
# Rename the Domain from Bacteria to Plant when it contains chloroplast or mitochondria
asv.tax.unfil.ed2 <- asv.tax.unfil.ed %>%
  mutate(Domain = case_when(Family == 'Chloroplast' ~ 'Plant',
                            Family == 'Mitochondria' ~ 'Plant', TRUE ~ Domain)) 
as.factor(asv.tax.unfil.ed2$Domain)
str(asv.tax.unfil.ed2)
dim(asv.tax.unfil.ed2) # 2051

colnames(asv.tax.unfil.ed2)
asv.tax.unfil.ed2.bio <- asv.tax.unfil.ed2[,c(1:64,86:92)] # unselect Mock and negative controls from the asv table
colnames(asv.tax.unfil.ed2.bio)
asv.tax.unfil.ed2.bio <- column_to_rownames(asv.tax.unfil.ed2.bio, var = "ASV.ID")
dim(asv.tax.unfil.ed2.bio) #2051


asv.tax.unfil.bio.noassig <- asv.tax.unfil.ed2.bio %>% 
  filter(Domain != "Unassigned")
dim(asv.tax.unfil.bio.noassig) #2038 --> 13 unassigned ASVs in the experimental sample


# Remove the unassigned ASV
asv.tax.unfil.noassig <- asv.tax.unfil.ed2 %>% 
  filter(Domain != "Unassigned")
dim(asv.tax.unfil.noassig) #2038
as.factor(asv.tax.unfil.noassig$Domain)

# 1. Plant read proportion in the true samples

# select only biological sample from the asv table
colnames(asv.tax.unfil.noassig)
asv.unfil.bio <- asv.tax.unfil.noassig[,1:64] # unselect Mock and negative controls from the asv table
asv.unfil.bio <- column_to_rownames(asv.unfil.bio, var = "ASV.ID")
sort(rowSums(asv.unfil.bio, na.rm = FALSE, dims = 1), decreasing = F)
sum(asv.unfil.bio)
# remove ASVs that do not present in biological sample
asv.unfil.bio_1 <- asv.unfil.bio[which(rowSums(asv.unfil.bio) > 0),]
dim(asv.unfil.bio_1) # [1] 1154  63, ASV table before plant contaminant removal and normalization (metagenomeSeq), and decontamination
sort(rowSums(asv.unfil.bio_1, na.rm = FALSE, dims = 1), decreasing = F)
sum(asv.unfil.bio_1)
# select the taxonomy from the asv table
colnames(asv.tax.unfil.noassig)
tax.unfil.bio <- asv.tax.unfil.noassig[,c(85:91)] #
head(tax.unfil.bio)
# merge the taxonomy with biological asv table
tax.unfil.bio <- rownames_to_column(tax.unfil.bio, var = "ASV.ID")
asv.unfil.bio_1 <- rownames_to_column(asv.unfil.bio_1, var = "ASV.ID")
asv.unfil.bio.tax <- merge(asv.unfil.bio_1, tax.unfil.bio, by="ASV.ID")
dim(asv.unfil.bio_1)
# remove OTUs with unassigned taxonomy
#otu.unfil.bio.tax.ed <- otu.unfil.bio.tax[!(otu.unfil.bio.tax$Domain=="Unassigned"),] # there are no unassigned OTU in the biological samples
# gather into long data for plotting
long.dat <- gather(asv.unfil.bio.tax, Sample, Read, 2:64, factor_key = T)
# calculate the proportion 
df.unfil.bio <- long.dat %>%
  group_by(Sample, Domain) %>%
  summarise(read.number = sum(Read))
df.unfil.bio.prop <- df.unfil.bio %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)
df.unfil.bio.prop # plant reads ~99 %, bacteria reads ~ 0,25 %
view(df.unfil.bio.prop)

df.unfil.bio.prop.group <- df.unfil.bio.prop %>% group_by(Domain) %>%
  summarise(read.number = sum(read.number))
view(df.unfil.bio.prop.group)

# make plot of the bacteria vs plant plot in the experimental samples
install.packages("ggtext")
library(ggplot2)
library(ggtext)
plot.unfil.bio.prop <- ggplot(df.unfil.bio.prop, aes(x=Domain, y=percent, fill=Domain))+
  geom_violin(trim = F, scale="width") +
  scale_fill_manual(labels = c("Bacteria","Plant"),values=c("#CC79A7","#009E73"))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
  theme_bw()+
  expand_limits(y = 0)+
  labs(title = "A. Experimental Sample")+
  ylab("Read Proportion (%)")+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text= element_text(size = 14),
        strip.text = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = 14),
        axis.title.y = element_markdown(size=15,face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_summary(fun="median",geom="point", size=7, color="red", shape=95)
plot.unfil.bio.prop

# save the plot
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("230724_PlantPropBioSample.tiff",
       plot.unfil.bio.prop, device = "tiff",
       width = 5, height =5, 
       units= "in", dpi = 300, 
       compression="lzw", bg= "white")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.Plant read proportion in the negative controls

# select only negative controls from ASV table
colnames(asv.tax.unfil.noassig)
nc <- asv.tax.unfil.noassig[, c(1,78:85)] # select only the negative control of the extraction day
# remove ASVs that do not present in the negative controls
nc <- column_to_rownames(nc, var = "ASV.ID")
sort(rowSums(nc, na.rm = FALSE, dims = 1), decreasing = F) # check the presence of ASV in the negative controls
nc_1 <- nc[which(rowSums(nc) > 0),] # select only ASVs that present in at least one of the negative control
nc_1 <- rownames_to_column(nc_1, var = "ASV.ID")
# merge with the taxonomy
nc_1_tax <- merge(nc_1, tax.unfil.bio, by ="ASV.ID")
as.factor(nc_1_tax$Domain)
# remove ASVs with unassigned taxonomy
#nc_1_tax_ed <- nc_1_tax[!(nc_1_tax$Domain=="Unassigned"),] # removed 3 unassigned OTUs
#nc_1_ed <- nc_1_tax_ed[,c(1:9)]
#rownames(nc_1_ed) <- NULL
#nc_1_ed <- column_to_rownames(nc_1_ed, var = "OTUID")
# gather into long data for plotting
nc_1_tax_ed <- nc_1_tax[!(nc_1_tax$Domain=="Archaea"),]
colnames(nc_1_tax_ed)
nc.long.dat <- gather(nc_1_tax_ed, Sample, Read, 2:9, factor_key = T)
# calculate the proportion 
df.nc <- nc.long.dat %>%
  group_by(Sample, Domain) %>%
  summarise(read.number = sum(Read))
df.nc.prop <- df.nc %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)
df.nc.prop

# make plot of the bacteria vs plant plot in the negative controls
plot.nc.prop <- ggplot(df.nc.prop, aes(x=Domain, y=percent, fill=Domain))+
  geom_violin(trim = F, scale="width") +
  scale_fill_manual(labels = c("Bacteria", "Plant"),values=c("#CC79A7","#009E73"))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
  theme_bw()+
  expand_limits(y = 0)+
  labs(title = "B. Negative Control")+
  ylab("Read Proportion (%)")+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text= element_text(size = 14),
        strip.text = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = 14),
        axis.title.y = element_markdown(size=15,face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_summary(fun="median",geom="point", size=7, color="red", shape=95)
plot.nc.prop

# save the plot
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("240724_PlantPropNegCont.tiff",
       plot.nc.prop, device = "tiff",
       width = 5, height =5, 
       units= "in", dpi = 300, 
       compression="lzw", bg= "white")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

# 3.Plant read proportion in the mock

# select only mock from asv table
mock <- asv.tax.unfil.noassig[, c(1,65:72)]
mock <- column_to_rownames(mock, var = "ASV.ID")
# remove ASVs that do not present in the mocks
sort(rowSums(mock, na.rm = FALSE, dims = 1), decreasing = F) # check the presence of ASV in the mock samples
mock_1 <- mock[which(rowSums(mock) > 0),] # select only ASVs that present in at least one of the mock samples
mock_1 <- rownames_to_column(mock_1, var = "ASV.ID")
# merge with the taxonomy
mock_1_tax <- merge(mock_1, tax.unfil.bio, by ="ASV.ID")
# remove ASVs with unassigned taxonomy
#mock_1_tax_ed <- mock_1_tax[!(mock_1_tax$Domain=="Unassigned"),] # removed 1 unassigned OTUs
#mock_1_ed <- mock_1_tax_ed[,c(1:9)]
#rownames(mock_1_ed) <- NULL
colnames(mock_1_tax)
mock.long.dat <- gather(mock_1_tax, Sample, Read, 2:9, factor_key = T)
# calculate the proportion 
df.mock <- mock.long.dat %>%
  group_by(Sample, Domain) %>%
  summarise(read.number = sum(Read))
df.mock.prop <- df.mock %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)
df.mock.prop

# make plot of the bacteria vs plant plot in the mock samples
plot.mock.prop <- ggplot(df.mock.prop, aes(x=Domain, y=percent, fill=Domain))+
  geom_violin(trim = F, scale="width") +
  scale_fill_manual(labels = c("Bacteria", "Plant"),values=c("#CC79A7","#009E73"))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
  theme_bw()+
  expand_limits(y = 0)+
  labs(title = "B. Negative Control")+
  ylab("Read Proportion (%)")+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text= element_text(size = 14),
        strip.text = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = 14),
        axis.title.y = element_markdown(size=15,face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_summary(fun="median",geom="point", size=7, color="red", shape=95)
plot.mock.prop

# save the plot
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("240724_PlantPropMock.tiff",
       plot.mock.prop, device = "tiff",
       width = 5, height =5, 
       units= "in", dpi = 300, 
       compression="lzw", bg= "white")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

# PF Reads before and after plant contaminants removal

setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/')
#wd <- print(getwd())
seq_data <- read.table("220706_Fina_16sFWD_220705.txt", sep='\t', header=T)
seq_data <- rownames_to_column(seq_data,var = "ID")
seq_data
# calculate plant-unfiltered reads from each sample
colnames(asv.tax.unfil.noassig)
asv.unfil.noassig <- asv.tax.unfil.noassig[,1:84] # unselect the taxonomy
dim(asv.unfil.noassig) # 2038
head(asv.unfil.noassig)
#asv.tax.unfil.noassig <- column_to_rownames(asv.tax.unfil.noassig, var="ASV.ID")
unfil.read <- setNames(nm=c('SampleID','UnfilteredRead'),stack(colSums(asv.unfil.noassig))[2:1])
unfil.read
# merge plant-unfiltered reads to the seq_data table
reads_data <- merge(seq_data, unfil.read, by = "SampleID")
reads_data
# filter plant contaminant and unassigned ASVs from the ASV table
asv.tax.unfil.bio.noassig <- asv.tax.unfil.ed2.bio %>% 
  filter(Domain != "Unassigned")
dim(asv.tax.unfil.bio.noassig) #2038 --> 13 unassigned ASVs in the experimental sample
colnames(asv.tax.unfil.bio.noassig)
asv.unfil.bio.noassig <- asv.tax.unfil.bio.noassig[,1:63]
asv.unfil.bio.noassig1 <- asv.unfil.bio.noassig[which(rowSums(asv.unfil.bio.noassig) > 0),]
sort(rowSums(asv.unfil.bio.noassig1, na.rm = FALSE, dims = 1), decreasing = F)
dim(asv.unfil.bio.noassig1) #1154 ASVs including plant's ASVs
asv.unfil.bio.noassig1 <- rownames_to_column(asv.unfil.bio.noassig1, var = 'ASV.ID')
asv.tax.unfil.bio.noassig1 <- merge(tax.unfil.bio,asv.unfil.bio.noassig1, by="ASV.ID")
dim(asv.tax.unfil.bio.noassig1)

asv.tax.fil.bio.noassig <- asv.tax.unfil.bio.noassig1 %>% 
  filter(Domain != "Plant")
dim(asv.tax.fil.bio.noassig) #813 --> REMOVED= 1154-813= 341 plant ASVs -->29.54%

asv.fil.tax <- asv.tax.unfil.noassig[!(asv.tax.unfil.noassig$Domain=="Plant" | 
                                         asv.tax.unfil.noassig$Domain=="Unassigned"),]#removed 355 plant ASVs 
dim(asv.tax.unfil.noassig)
dim(asv.fil.tax)
as.factor(asv.fil.tax$Domain)
# separate the ASV table and tax table from the filtered ASV table
colnames(asv.fil.tax)
#asv.fil_1 <- asv.fil[which(rowSums(asv.fil) > 0),]
asv.fil <- asv.fil.tax[,1:84]
head(asv.fil)
sort(rowSums(asv.fil, na.rm = FALSE, dims = 1), decreasing = F) # no zero
#rownames(asv.fil) <- NULL
# calculate plant-filtered read in true samples
colnames(asv.fil)
asv.fil.bio <- asv.fil[,1:63]
dim(asv.fil.bio)
asv.fil.bio_1 <- asv.fil.bio[which(rowSums(asv.fil.bio) > 0),]
colnames(asv.fil.bio_1)
dim(asv.fil.bio_1) #813
sum(asv.fil.bio_1)

colnames(asv.unfil.bio_1)
dim(asv.unfil.bio_1) # 1154

prop.removed.ASV <- ((1154-813)/1154)*100
prop.removed.ASV # 29.5%

# calculate plant-filtered read in negative controls
asv.fil.nc <- asv.fil[,77:84]
asv.fil.nc_1 <- asv.fil.nc[which(rowSums(asv.fil.nc) > 0),]
# calculate plant-filtered read in mocks
asv.fil.mock <- asv.fil[,64:71]
asv.fil.mock_1 <- asv.fil.mock[which(rowSums(asv.fil.mock) > 0),]
sum(asv.fil.mock_1)
#calculate plant-filtered reads from each sample
#asv.fil <- column_to_rownames(asv.fil, var = "ASV.ID")
fil.read <- setNames(nm=c('SampleID','FilteredRead'),stack(colSums(asv.fil))[2:1])
fil.read
# merge plant-filtered reads to the reads_data table
reads_data2 <- merge(reads_data, fil.read, by = "SampleID")
# re-order the data based on the ID
reads_data2$ID <- as.numeric(reads_data2$ID)
reads_data.ord <- reads_data2[order(reads_data2$ID),] #order the data based on the ID
reads_data.ord
# plot PF_reads by well position (A1, B1, C1 - H1, A2, B2, etc.)
library(ggplot2)

# subset a new data frame of the negative/positive controls
# neg controls first
ID <- c("10","20","30","39","48","57","69", "79")
ID <- as.numeric(ID)
neg_controls <- data.frame(ID)
neg_controls<- inner_join(neg_controls, reads_data.ord)
neg_controls

# NC
ID <- c("80","81", "82", "83", "84")
ID <- as.numeric(ID)
NC <- data.frame(ID)
NC<- inner_join(NC, reads_data.ord)
NC

# resave sample_ID data frame as positive controls
ID <- c("9","19","29","38","47","56","68","78")
ID <- as.numeric(ID)
pos_controls <- data.frame(ID)
pos_controls <- inner_join(pos_controls, reads_data.ord)


# make unfiltered read plot
reads_data.ord.df <- reads_data.ord %>% mutate(new_column = case_when(
  SampleID == "Mock1"~ "Positive Control",
  SampleID == "Mock2"~ "Positive Control",
  SampleID == "Mock3"~ "Positive Control",
  SampleID == "Mock4"~ "Positive Control",
  SampleID == "Mock5"~ "Positive Control",
  SampleID == "Mock6"~ "Positive Control",
  SampleID == "Mock7"~ "Positive Control",
  SampleID == "Mock8"~ "Positive Control",
  SampleID == "NegCon1"~ "Negative Control",
  SampleID == "NegCon2"~ "Negative Control",
  SampleID == "NegCon3"~ "Negative Control",
  SampleID == "NegCon4"~ "Negative Control",
  SampleID == "NegCon5"~ "Negative Control",
  SampleID == "NegCon6"~ "Negative Control",
  SampleID == "NegCon7"~ "Negative Control",
  SampleID == "NegCon8"~ "Negative Control",
  TRUE ~ "True Sample"))
reads_data.ord.df
  
values = c("True Sample" = "black", "Positive Control" = "blue", "Negative Control" = "red")
#seq_data_1$sample_ID <- as.factor(seq_data_1$sample_ID)
#seq_data_1$PF_reads <- as.numeric(seq_data_1$PF_reads)
pf.reads.unfil <- ggplot(data=reads_data.ord.df, aes(x=ID, y=UnfilteredRead, colour = new_column)) + 
                     theme_bw()+
                     geom_point(size=2.5)+
                     scale_y_continuous(labels = scales::comma)+
                     scale_color_manual(values = values)+
  labs(title = "PF Reads: Plant-unfiltered")+
  ylab("Read number (plant-unfiltered)")+ xlab("ID (well-ordered)")+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size=13),
        #legend.key.size = unit(3,"line"),
        axis.title.x = element_text(size=15),
        axis.text= element_text(size = 14),
        plot.title = element_text(size = 14),
        axis.title.y = element_markdown(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pf.reads.unfil

# save the plot
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("250724_PFReads_Unfiltered.tiff",
       pf.reads.unfil, device = "tiff",
       width = 8, height =5, 
       units= "in", dpi = 300, 
       compression="lzw", bg= "white")

# make filtered read plot
pf.reads.fil <- ggplot(data=reads_data.ord.df, aes(x=ID, y=FilteredRead, colour = new_column)) + 
  theme_bw()+
  geom_point(size=2.5)+
  scale_y_continuous(labels = scales::comma)+
  scale_color_manual(values = values)+
  labs(title = "PF Reads: Plant-filtered")+
  ylab("Read number (plant-filtered)")+ xlab("ID (well-ordered)")+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size=13),
        #legend.key.size = unit(3,"line"),
        axis.title.x = element_text(size=15),
        axis.text= element_text(size = 14),
        plot.title = element_text(size = 14),
        axis.title.y = element_markdown(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pf.reads.fil

# save the plot
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("250724_PFReads_Filtered.tiff",
       pf.reads.fil, device = "tiff",
       width = 8, height =5, 
       units= "in", dpi = 300, 
       compression="lzw", bg= "white")

# select only true samples and negative controls
reads_data_samp_nc <- reads_data.ord %>%
  dplyr::filter(!ID %in% c("9","19","29","38","47","56","68","78"))

# make filtered read plot
pf.reads.fil2 <- (ggplot(data=reads_data_samp_nc, aes(x=ID, y=FilteredRead)) + 
                    geom_point()+
                    scale_y_continuous(labels = scales::comma)+
                    geom_point(data=neg_controls, 
                               aes(x=ID,y=FilteredRead), 
                               color='red',
                               size=2) +
                    geom_point(data=NC, 
                               aes(x=ID,y=FilteredRead), 
                               color='green',
                               size=2))
pf.reads.fil2

######################################################################################################################################
# DECONTAMINATION OF THE POTENTIAL MICROBIAL CONTAMINANT: DECONTAM PACKAGE (METHOD=PREVALENCE)
######################################################################################################################################

#BiocManager::install("phyloseq")
#BiocManager::install("decontam")
#install.packages("DT")
#BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
#install.packages("microViz",repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos")))
library(DT)
library(phyloseq)
library(decontam)
library(microbiome)
library(ComplexHeatmap)
library(microViz)

# Remove contaminant reads/ASVs from the experimental sample's ASV table using the decontam package (method=prevalence)
# Decontamination analysis of potential microbial contaminant was performed using the decontam package using the prevalence method.
# The method will compare the prevalence of ASVs found in the either negative or positive control with ASVs in the experimental samples.
# The decontamination analysis was performed using the negative and positive controls for each DNA extraction batch.

# load the plant and unassigned-filtered ASV table
asv.fil.no.NC.tax <- asv.fil.tax[,c(1:71,77:91)]
colnames(asv.fil.no.NC.tax)
# separate the ASV and tax table
asv.fil.no.NC <- asv.fil.no.NC.tax[,1:79]
colnames(asv.fil.no.NC)
sort(rowSums(asv.fil.no.NC, na.rm = FALSE, dims = 1), decreasing = F) # some zeros
asv.fil.no.NC_1 <- asv.fil.no.NC[which(rowSums(asv.fil.no.NC) > 0),]
sort(rowSums(asv.fil.no.NC_1, na.rm = FALSE, dims = 1), decreasing = F) # no ASV with zero read
dim(asv.fil.no.NC_1) # 1499 79
asv.fil.no.NC_1 <- rownames_to_column(asv.fil.no.NC_1, var="ASV.ID")
# check if there are samples that have zero read
sort(colSums(asv.fil.no.NC_1, na.rm = FALSE, dims = 1), decreasing = F) # no sample with 0 read
# load the plant and unassigned-filtered taxonomy table
tax.fil.no.NC <- asv.fil.no.NC.tax[,80:86]
head(tax.fil.no.NC)
dim(tax.fil.no.NC) # 1683 
tax.fil.no.NC <- rownames_to_column(tax.fil.no.NC, var="ASV.ID")
# merge back to the ASV table which some ASVs with 0 read have been removed previously
asv.fil.no.NC.tax2 <- merge(asv.fil.no.NC_1, tax.fil.no.NC, by="ASV.ID")
# separate the tax 
tax.fil.no.NC2 <- asv.fil.no.NC.tax2[,c(1,81:87)]
dim(tax.fil.no.NC2) # 1499 8
tax.fil.no.NC2 <- column_to_rownames(tax.fil.no.NC2, var="ASV.ID")
# load the metadata which contain the type of sample and the extraction date
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/')
meta.decon <- read.csv(file = "extraction_date.csv", check.names = F)
dim(meta.decon) # 79 8
head(meta.decon)
# make phyloseq object for the decontamination
asv.fil.no.NC_1 <- column_to_rownames(asv.fil.no.NC_1, var="ASV.ID")
asv.fil.no.NC.seq <- otu_table(asv.fil.no.NC_1, taxa_are_rows = TRUE)
tax.fil.no.NC.seq <- tax_table(as.matrix(tax.fil.no.NC2))
meta.decon$SampleID <- as.factor(meta.decon$SampleID)
rownames(meta.decon) <- meta.decon$SampleID
meta.decon.seq <- sample_data(meta.decon)
# make phyloseq object
physeq.decon <- merge_phyloseq(asv.fil.no.NC.seq,tax.fil.no.NC.seq,meta.decon.seq)
physeq.decon
head(otu_table(physeq.decon))
sample_data(physeq.decon) #the order of sample data in the previous metadata will follow the order of the otu TABLE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting the library size/reads for each sample after plant contaminant removal
df <- as.data.frame(sample_data(physeq.decon))
df$LibrarySize <- sample_sums(physeq.decon)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
df
# full library size after plant-filtering
libsize.fil.all <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_type)) + 
  geom_point()+
  theme_bw()
libsize.fil.all
# full library size after plant-filtering with removed Mock
libsize.fil.noMock <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_type)) + 
  geom_point()+
  theme_bw()+
ylim (0,70000)
libsize.fil.noMock
# full library size after plant-filtering with removed Mock and removed 1 sample with high number of read (208: 69663 reads)
libsize.fil.noMock2 <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_type)) + 
  geom_point()+
  theme_bw()+
  ylim (0,30000)
libsize.fil.noMock2
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("260724_LibSize_Filtered_All.tiff",
       libsize.fil.all, device = "tiff",
       width = 8, height =5, 
       units= "in", dpi = 300, 
       compression="lzw", bg= "white")
ggsave("260724_LibSize_Filtered_noMock.tiff",
       libsize.fil.noMock, device = "tiff",
       width = 8, height =5, 
       units= "in", dpi = 300, 
       compression="lzw", bg= "white")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

################ Detecting contaminant with decontam package per extraction batch ############################

# Decontam by extraction group :


# 1. Batch T1
library(microViz)

physeq_t1 <- physeq.decon %>% 
  ps_filter(batch == "1", .keep_all_taxa = TRUE)
sample_data(physeq_t1)
# default threshold (0,1)
sample_data(physeq_t1)$is.neg <- sample_data(physeq_t1)$sample_type == "negative"
contamdf.prev.01_t1 <- isContaminant(physeq_t1, method="prevalence", neg="is.neg")
table(contamdf.prev.01_t1$contaminant) # FALSE: 1497 TRUE:2 
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t1 <- transform_sample_counts(physeq_t1, function(abund) 1*(abund>0))
ps.pa.neg_t1 <- prune_samples(sample_data(ps.pa_t1)$is.neg == TRUE, ps.pa_t1)
ps.pa.pos_t1 <- prune_samples(sample_data(ps.pa_t1)$is.neg == FALSE, ps.pa_t1)
# Make data.frame of prevalence in positive and negative samples
df.pa_t1 <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t1), pa.neg=taxa_sums(ps.pa.neg_t1),
                    contaminant=contamdf.prev.01_t1$contaminant)
df.pa_t1
t1_decontam_plot <- ggplot(data=df.pa_t1, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t1 Decontam")
t1_decontam_plot
# grab the contaminant rows
row_indices_t1 <- which(contamdf.prev.01_t1$contaminant)
row_indices_t1
# extract the taxonomic information of the contaminants
webshot::install_phantomjs(force = TRUE)
taxonomy_table <- tibble()

for (i in row_indices_t1){
  loc <-  contamdf.prev.01_t1[i, 0]
  tax_key <- row.names(loc)
  tax_value <- tax.fil[tax_key, ]
  taxonomy_table <- rbind(taxonomy_table, tax_value)
}

cont.t1.neg <- datatable(taxonomy_table) 
as.data.frame(cont.t1.neg)
# save the contaminant table in the local computer
write.csv(df.pa_t1, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/contaminant-table-t1-neg.csv")
#removing contaminants from phyloseq object
colnames(df.pa_t1)
df.pa_t1 <- rownames_to_column(df.pa_t1, var = "ASV.ID")
subset.df_t1 <- subset(df.pa_t1, contaminant== "FALSE")
keep.taxa_t1 <- as.vector(subset.df_t1$ASV.ID)
length(keep.taxa_t1)

subset.df.remove_t1 <- subset(df.pa_t1, contaminant== "TRUE")
remove.taxa.list_t1 <- as.vector(subset.df.remove_t1)
remove.taxa.df_t1 <- data.frame(remove.taxa.list_t1)

decontam_phyloseq_t1 <- prune_taxa(keep.taxa_t1, physeq_t1)
decontam_phyloseq_t1
#filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t1<- subset(sample_data(decontam_phyloseq_t1), is.neg ==FALSE)
keep.samples_t1 <- as.vector(subset.metadata_t1$SampleID)
keep.samples_t1
decontam_phyloseq_t1 <- prune_samples(keep.samples_t1, decontam_phyloseq_t1)

### decontam with positive control

sample_data(decontam_phyloseq_t1)$is.neg <- sample_data(decontam_phyloseq_t1)$sample_type == "positive"
contamdf.prev0.1_t1pos <- isContaminant(decontam_phyloseq_t1, method="prevalence", neg="is.neg")
table(contamdf.prev0.1_t1pos$contaminant)# FALSE=1497  TRUE=0
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t1pos <- transform_sample_counts(decontam_phyloseq_t1, function(abund) 1*(abund>0))
ps.pa.neg_t1pos <- prune_samples(sample_data(ps.pa_t1pos)$is.neg == TRUE, ps.pa_t1pos)
ps.pa.pos_t1pos <- prune_samples(sample_data(ps.pa_t1pos)$is.neg == FALSE, ps.pa_t1pos)
# Make data.frame of prevalence in positive and negative samples
df.pa_t1pos <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t1pos), pa.neg=taxa_sums(ps.pa.neg_t1pos),
                    contaminant=contamdf.prev0.1_t1pos$contaminant)
t1pos_decontam_plot <- ggplot(data=df.pa_t1pos, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t1 Decontam")
t1pos_decontam_plot
# save the contaminant table in the local computer
write.csv(df.pa_t1pos, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/contaminant-table-t1-pos.csv")
##removing contaminants from phyloseq object
colnames(df.pa_t1pos)
df.pa_t1pos <- rownames_to_column(df.pa_t1pos, var = "ASV.ID")
subset.df_t1pos <- subset(df.pa_t1pos, contaminant== "FALSE")
keep.taxa_t1pos <- as.vector(subset.df_t1pos$ASV.ID)
length(keep.taxa_t1pos)

subset.df.remove_t1pos <- subset(df.pa_t1pos, contaminant== "TRUE")
remove.taxa.list_t1pos <- as.vector(subset.df.remove_t1pos)
remove.taxa.df_t1pos <- data.frame(remove.taxa.list_t1pos)

decontam_phyloseq_t1 <- prune_taxa(keep.taxa_t1pos, decontam_phyloseq_t1)
decontam_phyloseq_t1
##filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t1<- subset(sample_data(decontam_phyloseq_t1), is.neg ==FALSE)
keep.samples_t1 <- as.vector(subset.metadata_t1$SampleID)
keep.samples_t1
decontam_phyloseq_t1 <- prune_samples(keep.samples_t1, decontam_phyloseq_t1)
decontam_phyloseq_t1 
otu_table(decontam_phyloseq_t1)
#save
write.csv(otu_table(decontam_phyloseq_t1), "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam_ASV_t1pos.csv")

#______________________________________________________________________________________________________________________________________________________________________


# 2. Batch T2


physeq_t2 <- physeq.decon %>% 
  ps_filter(batch == "2", .keep_all_taxa = TRUE)
sample_data(physeq_t2)
physeq_t2
# default threshold (0,1)
sample_data(physeq_t2)$is.neg <- sample_data(physeq_t2)$sample_type == "negative"
contamdf.prev.01_t2 <- isContaminant(physeq_t2, method="prevalence", neg="is.neg")
table(contamdf.prev.01_t2$contaminant) # FALSE: 1496 TRUE: 3
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t2 <- transform_sample_counts(physeq_t2, function(abund) 1*(abund>0))
ps.pa.neg_t2 <- prune_samples(sample_data(ps.pa_t2)$is.neg == TRUE, ps.pa_t2)
ps.pa.pos_t2 <- prune_samples(sample_data(ps.pa_t2)$is.neg == FALSE, ps.pa_t2)
# Make data.frame of prevalence in positive and negative samples
df.pa_t2 <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t2), pa.neg=taxa_sums(ps.pa.neg_t2),
                       contaminant=contamdf.prev.01_t2$contaminant)
df.pa_t2
t2_decontam_plot <- ggplot(data=df.pa_t2, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t2 Decontam")
t2_decontam_plot
# grab the contaminant rows
row_indices_t2 <- which(contamdf.prev.01_t2$contaminant)
row_indices_t2
# extract the taxonomic information of the contaminants
#webshot::install_phantomjs()
taxonomy_table <- tibble()

for (i in row_indices_t2){
  loc <-  contamdf.prev.01_t2[i, 0]
  tax_key <- row.names(loc)
  tax_value <- tax.fil[tax_key, ]
  taxonomy_table <- rbind(taxonomy_table, tax_value)
}

datatable(taxonomy_table) 

# save the contaminant table in the local computer
write.csv(df.pa_t2, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t2-neg.csv")
#removing contaminants from phyloseq object
colnames(df.pa_t2)
df.pa_t2 <- rownames_to_column(df.pa_t2, var = "ASV.ID")
subset.df_t2 <- subset(df.pa_t2, contaminant== "FALSE")
keep.taxa_t2 <- as.vector(subset.df_t2$ASV.ID)
length(keep.taxa_t2)

subset.df.remove_t2 <- subset(df.pa_t2, contaminant== "TRUE")
remove.taxa.list_t2 <- as.vector(subset.df.remove_t2)
remove.taxa.df_t2 <- data.frame(remove.taxa.list_t2)

decontam_phyloseq_t2 <- prune_taxa(keep.taxa_t2, physeq_t2)
decontam_phyloseq_t2
#filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t2<- subset(sample_data(decontam_phyloseq_t2), is.neg ==FALSE)
keep.samples_t2 <- as.vector(subset.metadata_t2$SampleID)
keep.samples_t2
decontam_phyloseq_t2 <- prune_samples(keep.samples_t2, decontam_phyloseq_t2)

### decontam with positive control

sample_data(decontam_phyloseq_t2)$is.neg <- sample_data(decontam_phyloseq_t2)$sample_type == "positive"
contamdf.prev0.1_t2pos <- isContaminant(decontam_phyloseq_t2, method="prevalence", neg="is.neg")
table(contamdf.prev0.1_t2pos$contaminant)# FALSE=1496  TRUE=0
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t2pos <- transform_sample_counts(decontam_phyloseq_t2, function(abund) 1*(abund>0))
ps.pa.neg_t2pos <- prune_samples(sample_data(ps.pa_t2pos)$is.neg == TRUE, ps.pa_t2pos)
ps.pa.pos_t2pos <- prune_samples(sample_data(ps.pa_t2pos)$is.neg == FALSE, ps.pa_t2pos)
# Make data.frame of prevalence in positive and negative samples
df.pa_t2pos <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t2pos), pa.neg=taxa_sums(ps.pa.neg_t2pos),
                          contaminant=contamdf.prev0.1_t2pos$contaminant)
t2pos_decontam_plot <- ggplot(data=df.pa_t2pos, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t2 Decontam")
t2pos_decontam_plot
# save the contaminant table in the local computer
write.csv(df.pa_t2pos, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t2-pos.csv")
##removing contaminants from phyloseq object
colnames(df.pa_t2pos)
df.pa_t2pos <- rownames_to_column(df.pa_t2pos, var = "ASV.ID")
subset.df_t2pos <- subset(df.pa_t2pos, contaminant== "FALSE")
keep.taxa_t2pos <- as.vector(subset.df_t2pos$ASV.ID)
length(keep.taxa_t2pos)

subset.df.remove_t2pos <- subset(df.pa_t2pos, contaminant== "TRUE")
remove.taxa.list_t2pos <- as.vector(subset.df.remove_t2pos)
remove.taxa.df_t2pos <- data.frame(remove.taxa.list_t2pos)

decontam_phyloseq_t2 <- prune_taxa(keep.taxa_t2pos, decontam_phyloseq_t2)
decontam_phyloseq_t2
##filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t2<- subset(sample_data(decontam_phyloseq_t2), is.neg ==FALSE)
keep.samples_t2 <- as.vector(subset.metadata_t2$SampleID)
keep.samples_t2
decontam_phyloseq_t2 <- prune_samples(keep.samples_t2, decontam_phyloseq_t2)
otu_table(decontam_phyloseq_t2)
#save
write.csv(otu_table(decontam_phyloseq_t2), "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/decontam_ASV_t2pos.csv")

#______________________________________________________________________________________________________________________________________________________________________


# 3. Batch T3


physeq_t3 <- physeq.decon %>% 
  ps_filter(batch == "3", .keep_all_taxa = TRUE)
sample_data(physeq_t3)
physeq_t3
# default threshold (0,1)
sample_data(physeq_t3)$is.neg <- sample_data(physeq_t3)$sample_type == "negative"
contamdf.prev.01_t3 <- isContaminant(physeq_t3, method="prevalence", neg="is.neg")
table(contamdf.prev.01_t3$contaminant) # FALSE: 1495 TRUE: 4
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t3 <- transform_sample_counts(physeq_t3, function(abund) 1*(abund>0))
ps.pa.neg_t3 <- prune_samples(sample_data(ps.pa_t3)$is.neg == TRUE, ps.pa_t3)
ps.pa.pos_t3 <- prune_samples(sample_data(ps.pa_t3)$is.neg == FALSE, ps.pa_t3)
# Make data.frame of prevalence in positive and negative samples
df.pa_t3 <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t3), pa.neg=taxa_sums(ps.pa.neg_t3),
                       contaminant=contamdf.prev.01_t3$contaminant)
df.pa_t3
t3_decontam_plot <- ggplot(data=df.pa_t3, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t3 Decontam")
t3_decontam_plot
# grab the contaminant rows
row_indices_t3 <- which(contamdf.prev.01_t3$contaminant)
row_indices_t3
# extract the taxonomic information of the contaminants
#webshot::install_phantomjs()
taxonomy_table <- tibble()

for (i in row_indices_t3){
  loc <-  contamdf.prev.01_t3[i, 0]
  tax_key <- row.names(loc)
  tax_value <- tax.fil[tax_key, ]
  taxonomy_table <- rbind(taxonomy_table, tax_value)
}

datatable(taxonomy_table) 

# save the contaminant table in the local computer
write.csv(df.pa_t3, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t3-neg.csv")
#removing contaminants from phyloseq object
colnames(df.pa_t3)
df.pa_t3 <- rownames_to_column(df.pa_t3, var = "ASV.ID")
subset.df_t3 <- subset(df.pa_t3, contaminant== "FALSE")
keep.taxa_t3 <- as.vector(subset.df_t3$ASV.ID)
length(keep.taxa_t3)

subset.df.remove_t3 <- subset(df.pa_t3, contaminant== "TRUE")
remove.taxa.list_t3 <- as.vector(subset.df.remove_t3)
remove.taxa.df_t3 <- data.frame(remove.taxa.list_t3)

decontam_phyloseq_t3 <- prune_taxa(keep.taxa_t3, physeq_t3)
decontam_phyloseq_t3
#filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t3<- subset(sample_data(decontam_phyloseq_t3), is.neg ==FALSE)
keep.samples_t3 <- as.vector(subset.metadata_t3$SampleID)
keep.samples_t3
decontam_phyloseq_t3 <- prune_samples(keep.samples_t3, decontam_phyloseq_t3)

### decontam with positive control

sample_data(decontam_phyloseq_t3)$is.neg <- sample_data(decontam_phyloseq_t3)$sample_type == "positive"
contamdf.prev0.1_t3pos <- isContaminant(decontam_phyloseq_t3, method="prevalence", neg="is.neg")
table(contamdf.prev0.1_t3pos$contaminant)# FALSE=1495  TRUE=0
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t3pos <- transform_sample_counts(decontam_phyloseq_t3, function(abund) 1*(abund>0))
ps.pa.neg_t3pos <- prune_samples(sample_data(ps.pa_t3pos)$is.neg == TRUE, ps.pa_t3pos)
ps.pa.pos_t3pos <- prune_samples(sample_data(ps.pa_t3pos)$is.neg == FALSE, ps.pa_t3pos)
# Make data.frame of prevalence in positive and negative samples
df.pa_t3pos <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t3pos), pa.neg=taxa_sums(ps.pa.neg_t3pos),
                          contaminant=contamdf.prev0.1_t3pos$contaminant)
t3pos_decontam_plot <- ggplot(data=df.pa_t3pos, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t3 Decontam")
t3pos_decontam_plot
# save the contaminant table in the local computer
write.csv(df.pa_t3pos, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t3-pos.csv")
##removing contaminants from phyloseq object
colnames(df.pa_t3pos)
df.pa_t3pos <- rownames_to_column(df.pa_t3pos, var = "ASV.ID")
subset.df_t3pos <- subset(df.pa_t3pos, contaminant== "FALSE")
keep.taxa_t3pos <- as.vector(subset.df_t3pos$ASV.ID)
length(keep.taxa_t3pos)

subset.df.remove_t3pos <- subset(df.pa_t3pos, contaminant== "TRUE")
remove.taxa.list_t3pos <- as.vector(subset.df.remove_t3pos)
remove.taxa.df_t3pos <- data.frame(remove.taxa.list_t3pos)

decontam_phyloseq_t3 <- prune_taxa(keep.taxa_t3pos, decontam_phyloseq_t3)
decontam_phyloseq_t3
##filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t3 <- subset(sample_data(decontam_phyloseq_t3), is.neg ==FALSE)
keep.samples_t3 <- as.vector(subset.metadata_t3$SampleID)
keep.samples_t3
decontam_phyloseq_t3 <- prune_samples(keep.samples_t3, decontam_phyloseq_t3)
otu_table(decontam_phyloseq_t3)
#save
write.csv(otu_table(decontam_phyloseq_t3), "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/decontam_ASV_t3pos.csv")


#______________________________________________________________________________________________________________________________________________________________________


# 4. Batch T4


physeq_t4 <- physeq.decon %>% 
  ps_filter(batch == "4", .keep_all_taxa = TRUE)
sample_data(physeq_t4)
physeq_t4
# default threshold (0,1)
sample_data(physeq_t4)$is.neg <- sample_data(physeq_t4)$sample_type == "negative"
contamdf.prev.01_t4 <- isContaminant(physeq_t4, method="prevalence", neg="is.neg")
table(contamdf.prev.01_t4$contaminant) # FALSE: 1499
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t4 <- transform_sample_counts(physeq_t4, function(abund) 1*(abund>0))
ps.pa.neg_t4 <- prune_samples(sample_data(ps.pa_t4)$is.neg == TRUE, ps.pa_t4)
ps.pa.pos_t4 <- prune_samples(sample_data(ps.pa_t4)$is.neg == FALSE, ps.pa_t4)
# Make data.frame of prevalence in positive and negative samples
df.pa_t4 <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t4), pa.neg=taxa_sums(ps.pa.neg_t4),
                       contaminant=contamdf.prev.01_t4$contaminant)
df.pa_t3
t4_decontam_plot <- ggplot(data=df.pa_t4, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t4 Decontam")
t4_decontam_plot
# grab the contaminant rows
row_indices_t4 <- which(contamdf.prev.01_t4$contaminant)
row_indices_t4
# extract the taxonomic information of the contaminants
#webshot::install_phantomjs()
taxonomy_table <- tibble()

for (i in row_indices_t4){
  loc <-  contamdf.prev.01_t4[i, 0]
  tax_key <- row.names(loc)
  tax_value <- tax.fil[tax_key, ]
  taxonomy_table <- rbind(taxonomy_table, tax_value)
}

datatable(taxonomy_table) 

# save the contaminant table in the local computer
write.csv(df.pa_t4, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t4-neg.csv")
#removing contaminants from phyloseq object
colnames(df.pa_t4)
df.pa_t4 <- rownames_to_column(df.pa_t4, var = "ASV.ID")
subset.df_t4 <- subset(df.pa_t4, contaminant== "FALSE")
keep.taxa_t4 <- as.vector(subset.df_t4$ASV.ID)
length(keep.taxa_t4)

subset.df.remove_t4 <- subset(df.pa_t4, contaminant== "TRUE")
remove.taxa.list_t4 <- as.vector(subset.df.remove_t4)
remove.taxa.df_t4 <- data.frame(remove.taxa.list_t4)

decontam_phyloseq_t4 <- prune_taxa(keep.taxa_t4, physeq_t4)
decontam_phyloseq_t4
#filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t4 <- subset(sample_data(decontam_phyloseq_t4), is.neg ==FALSE)
keep.samples_t4 <- as.vector(subset.metadata_t4$SampleID)
keep.samples_t4
decontam_phyloseq_t4 <- prune_samples(keep.samples_t4, decontam_phyloseq_t4)

### decontam with positive control

sample_data(decontam_phyloseq_t4)$is.neg <- sample_data(decontam_phyloseq_t4)$sample_type == "positive"
contamdf.prev0.1_t4pos <- isContaminant(decontam_phyloseq_t4, method="prevalence", neg="is.neg")
table(contamdf.prev0.1_t4pos$contaminant)# FALSE=1499
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t4pos <- transform_sample_counts(decontam_phyloseq_t4, function(abund) 1*(abund>0))
ps.pa.neg_t4pos <- prune_samples(sample_data(ps.pa_t4pos)$is.neg == TRUE, ps.pa_t4pos)
ps.pa.pos_t4pos <- prune_samples(sample_data(ps.pa_t4pos)$is.neg == FALSE, ps.pa_t4pos)
# Make data.frame of prevalence in positive and negative samples
df.pa_t4pos <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t4pos), pa.neg=taxa_sums(ps.pa.neg_t4pos),
                          contaminant=contamdf.prev0.1_t4pos$contaminant)
t4pos_decontam_plot <- ggplot(data=df.pa_t4pos, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t4 Decontam")
t4pos_decontam_plot
# save the contaminant table in the local computer
write.csv(df.pa_t4pos, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t4-pos.csv")
##removing contaminants from phyloseq object
colnames(df.pa_t4pos)
df.pa_t4pos <- rownames_to_column(df.pa_t4pos, var = "ASV.ID")
subset.df_t4pos <- subset(df.pa_t4pos, contaminant== "FALSE")
keep.taxa_t4pos <- as.vector(subset.df_t4pos$ASV.ID)
length(keep.taxa_t4pos)

subset.df.remove_t4pos <- subset(df.pa_t4pos, contaminant== "TRUE")
remove.taxa.list_t4pos <- as.vector(subset.df.remove_t4pos)
remove.taxa.df_t4pos <- data.frame(remove.taxa.list_t4pos)

decontam_phyloseq_t4 <- prune_taxa(keep.taxa_t4pos, decontam_phyloseq_t4)
decontam_phyloseq_t4
##filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t4 <- subset(sample_data(decontam_phyloseq_t4), is.neg ==FALSE)
keep.samples_t4 <- as.vector(subset.metadata_t4$SampleID)
keep.samples_t4
decontam_phyloseq_t4 <- prune_samples(keep.samples_t4, decontam_phyloseq_t4)
otu_table(decontam_phyloseq_t4)
#save
write.csv(otu_table(decontam_phyloseq_t4), "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/decontam_ASV_t4pos.csv")


#______________________________________________________________________________________________________________________________________________________________________


# 5. Batch T5


physeq_t5 <- physeq.decon %>% 
  ps_filter(batch == "5", .keep_all_taxa = TRUE)
sample_data(physeq_t5)
physeq_t5
# default threshold (0,1)
sample_data(physeq_t5)$is.neg <- sample_data(physeq_t5)$sample_type == "negative"
contamdf.prev.01_t5 <- isContaminant(physeq_t5, method="prevalence", neg="is.neg")
table(contamdf.prev.01_t5$contaminant) # FALSE: 1499
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t5 <- transform_sample_counts(physeq_t5, function(abund) 1*(abund>0))
ps.pa.neg_t5 <- prune_samples(sample_data(ps.pa_t5)$is.neg == TRUE, ps.pa_t5)
ps.pa.pos_t5 <- prune_samples(sample_data(ps.pa_t5)$is.neg == FALSE, ps.pa_t5)
# Make data.frame of prevalence in positive and negative samples
df.pa_t5 <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t5), pa.neg=taxa_sums(ps.pa.neg_t5),
                       contaminant=contamdf.prev.01_t5$contaminant)
df.pa_t5
t5_decontam_plot <- ggplot(data=df.pa_t5, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t5 Decontam")
t5_decontam_plot
# grab the contaminant rows
row_indices_t5 <- which(contamdf.prev.01_t5$contaminant)
row_indices_t5
# extract the taxonomic information of the contaminants
#webshot::install_phantomjs()
taxonomy_table <- tibble()

for (i in row_indices_t5){
  loc <-  contamdf.prev.01_t5[i, 0]
  tax_key <- row.names(loc)
  tax_value <- tax.fil[tax_key, ]
  taxonomy_table <- rbind(taxonomy_table, tax_value)
}

datatable(taxonomy_table) 

# save the contaminant table in the local computer
write.csv(df.pa_t5, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t5-neg.csv")
#removing contaminants from phyloseq object
colnames(df.pa_t5)
df.pa_t5 <- rownames_to_column(df.pa_t5, var = "ASV.ID")
subset.df_t5 <- subset(df.pa_t5, contaminant== "FALSE")
keep.taxa_t5 <- as.vector(subset.df_t5$ASV.ID)
length(keep.taxa_t5)

subset.df.remove_t5 <- subset(df.pa_t5, contaminant== "TRUE")
remove.taxa.list_t5 <- as.vector(subset.df.remove_t5)
remove.taxa.df_t5 <- data.frame(remove.taxa.list_t5)

decontam_phyloseq_t5 <- prune_taxa(keep.taxa_t5, physeq_t5)
decontam_phyloseq_t5
#filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t5<- subset(sample_data(decontam_phyloseq_t5), is.neg ==FALSE)
keep.samples_t5 <- as.vector(subset.metadata_t5$SampleID)
keep.samples_t5
decontam_phyloseq_t5 <- prune_samples(keep.samples_t5, decontam_phyloseq_t5)

### decontam with positive control

sample_data(decontam_phyloseq_t5)$is.neg <- sample_data(decontam_phyloseq_t5)$sample_type == "positive"
contamdf.prev0.1_t5pos <- isContaminant(decontam_phyloseq_t5, method="prevalence", neg="is.neg")
table(contamdf.prev0.1_t5pos$contaminant)# FALSE=1499
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t5pos <- transform_sample_counts(decontam_phyloseq_t5, function(abund) 1*(abund>0))
ps.pa.neg_t5pos <- prune_samples(sample_data(ps.pa_t5pos)$is.neg == TRUE, ps.pa_t5pos)
ps.pa.pos_t5pos <- prune_samples(sample_data(ps.pa_t5pos)$is.neg == FALSE, ps.pa_t5pos)
# Make data.frame of prevalence in positive and negative samples
df.pa_t5pos <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t5pos), pa.neg=taxa_sums(ps.pa.neg_t5pos),
                          contaminant=contamdf.prev0.1_t5pos$contaminant)
t5pos_decontam_plot <- ggplot(data=df.pa_t5pos, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t5 Decontam")
t5pos_decontam_plot
# save the contaminant table in the local computer
write.csv(df.pa_t5pos, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t5-pos.csv")
##removing contaminants from phyloseq object
colnames(df.pa_t5pos)
df.pa_t5pos <- rownames_to_column(df.pa_t5pos, var = "ASV.ID")
subset.df_t5pos <- subset(df.pa_t5pos, contaminant== "FALSE")
keep.taxa_t5pos <- as.vector(subset.df_t5pos$ASV.ID)
length(keep.taxa_t5pos)

subset.df.remove_t5pos <- subset(df.pa_t5pos, contaminant== "TRUE")
remove.taxa.list_t5pos <- as.vector(subset.df.remove_t5pos)
remove.taxa.df_t5pos <- data.frame(remove.taxa.list_t5pos)

decontam_phyloseq_t5 <- prune_taxa(keep.taxa_t5pos, decontam_phyloseq_t5)
decontam_phyloseq_t5
##filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t5 <- subset(sample_data(decontam_phyloseq_t5), is.neg ==FALSE)
keep.samples_t5 <- as.vector(subset.metadata_t5$SampleID)
keep.samples_t5
decontam_phyloseq_t5 <- prune_samples(keep.samples_t5, decontam_phyloseq_t5)
otu_table(decontam_phyloseq_t5)
#save
write.csv(otu_table(decontam_phyloseq_t5), "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/decontam_ASV_t5pos.csv")

#______________________________________________________________________________________________________________________________________________________________________


# 6. Batch T6


physeq_t6 <- physeq.decon %>% 
  ps_filter(batch == "6", .keep_all_taxa = TRUE)
sample_data(physeq_t6)
physeq_t6
# default threshold (0,1)
sample_data(physeq_t6)$is.neg <- sample_data(physeq_t6)$sample_type == "negative"
contamdf.prev.01_t6 <- isContaminant(physeq_t6, method="prevalence", neg="is.neg")
table(contamdf.prev.01_t6$contaminant) # FALSE: 1499
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t6 <- transform_sample_counts(physeq_t6, function(abund) 1*(abund>0))
ps.pa.neg_t6 <- prune_samples(sample_data(ps.pa_t6)$is.neg == TRUE, ps.pa_t6)
ps.pa.pos_t6 <- prune_samples(sample_data(ps.pa_t6)$is.neg == FALSE, ps.pa_t6)
# Make data.frame of prevalence in positive and negative samples
df.pa_t6 <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t6), pa.neg=taxa_sums(ps.pa.neg_t6),
                       contaminant=contamdf.prev.01_t6$contaminant)
df.pa_t6
t6_decontam_plot <- ggplot(data=df.pa_t6, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t6 Decontam")
t6_decontam_plot
# grab the contaminant rows
row_indices_t6 <- which(contamdf.prev.01_t6$contaminant)
row_indices_t6
# extract the taxonomic information of the contaminants
#webshot::install_phantomjs()
taxonomy_table <- tibble()

for (i in row_indices_t6){
  loc <-  contamdf.prev.01_t6[i, 0]
  tax_key <- row.names(loc)
  tax_value <- tax.fil[tax_key, ]
  taxonomy_table <- rbind(taxonomy_table, tax_value)
}

datatable(taxonomy_table) 

# save the contaminant table in the local computer
write.csv(df.pa_t6, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t6-neg.csv")
#removing contaminants from phyloseq object
colnames(df.pa_t6)
df.pa_t6 <- rownames_to_column(df.pa_t6, var = "ASV.ID")
subset.df_t6 <- subset(df.pa_t6, contaminant== "FALSE")
keep.taxa_t6 <- as.vector(subset.df_t6$ASV.ID)
length(keep.taxa_t6)

subset.df.remove_t6 <- subset(df.pa_t6, contaminant== "TRUE")
remove.taxa.list_t6 <- as.vector(subset.df.remove_t6)
remove.taxa.df_t6 <- data.frame(remove.taxa.list_t6)

decontam_phyloseq_t6 <- prune_taxa(keep.taxa_t6, physeq_t6)
decontam_phyloseq_t6
#filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t6 <- subset(sample_data(decontam_phyloseq_t6), is.neg ==FALSE)
keep.samples_t6 <- as.vector(subset.metadata_t6$SampleID)
keep.samples_t6
decontam_phyloseq_t6 <- prune_samples(keep.samples_t6, decontam_phyloseq_t6)

### decontam with positive control

sample_data(decontam_phyloseq_t6)$is.neg <- sample_data(decontam_phyloseq_t6)$sample_type == "positive"
contamdf.prev0.1_t6pos <- isContaminant(decontam_phyloseq_t6, method="prevalence", neg="is.neg")
table(contamdf.prev0.1_t6pos$contaminant)# FALSE=1499  TRUE=0
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t6pos <- transform_sample_counts(decontam_phyloseq_t6, function(abund) 1*(abund>0))
ps.pa.neg_t6pos <- prune_samples(sample_data(ps.pa_t6pos)$is.neg == TRUE, ps.pa_t6pos)
ps.pa.pos_t6pos <- prune_samples(sample_data(ps.pa_t6pos)$is.neg == FALSE, ps.pa_t6pos)
# Make data.frame of prevalence in positive and negative samples
df.pa_t6pos <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t6pos), pa.neg=taxa_sums(ps.pa.neg_t6pos),
                          contaminant=contamdf.prev0.1_t6pos$contaminant)
t6pos_decontam_plot <- ggplot(data=df.pa_t6pos, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t6 Decontam")
t6pos_decontam_plot
# save the contaminant table in the local computer
write.csv(df.pa_t6pos, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t6-pos.csv")
##removing contaminants from phyloseq object
colnames(df.pa_t6pos)
df.pa_t6pos <- rownames_to_column(df.pa_t6pos, var = "ASV.ID")
subset.df_t6pos <- subset(df.pa_t6pos, contaminant== "FALSE")
keep.taxa_t6pos <- as.vector(subset.df_t6pos$ASV.ID)
length(keep.taxa_t6pos)

subset.df.remove_t6pos <- subset(df.pa_t6pos, contaminant== "TRUE")
remove.taxa.list_t6pos <- as.vector(subset.df.remove_t6pos)
remove.taxa.df_t6pos <- data.frame(remove.taxa.list_t6pos)

decontam_phyloseq_t6 <- prune_taxa(keep.taxa_t6pos, decontam_phyloseq_t6)
decontam_phyloseq_t6
##filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t6 <- subset(sample_data(decontam_phyloseq_t6), is.neg ==FALSE)
keep.samples_t6 <- as.vector(subset.metadata_t6$SampleID)
keep.samples_t6
decontam_phyloseq_t6 <- prune_samples(keep.samples_t6, decontam_phyloseq_t6)
otu_table(decontam_phyloseq_t6)
#save
write.csv(otu_table(decontam_phyloseq_t6), "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/decontam_ASV_t6pos.csv")

#______________________________________________________________________________________________________________________________________________________________________


# 7. Batch T7


physeq_t7 <- physeq.decon %>% 
  ps_filter(batch == "7", .keep_all_taxa = TRUE)
sample_data(physeq_t7)
physeq_t7
# default threshold (0,1)
sample_data(physeq_t7)$is.neg <- sample_data(physeq_t7)$sample_type == "negative"
contamdf.prev.01_t7 <- isContaminant(physeq_t7, method="prevalence", neg="is.neg")
table(contamdf.prev.01_t7$contaminant) # FALSE: 1488 TRUE: 11
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t7 <- transform_sample_counts(physeq_t7, function(abund) 1*(abund>0))
ps.pa.neg_t7 <- prune_samples(sample_data(ps.pa_t7)$is.neg == TRUE, ps.pa_t7)
ps.pa.pos_t7 <- prune_samples(sample_data(ps.pa_t7)$is.neg == FALSE, ps.pa_t7)
# Make data.frame of prevalence in positive and negative samples
df.pa_t7 <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t7), pa.neg=taxa_sums(ps.pa.neg_t7),
                       contaminant=contamdf.prev.01_t7$contaminant)
df.pa_t7
t7_decontam_plot <- ggplot(data=df.pa_t7, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t7 Decontam")
t7_decontam_plot
# grab the contaminant rows
row_indices_t7 <- which(contamdf.prev.01_t7$contaminant)
row_indices_t7
# extract the taxonomic information of the contaminants
#webshot::install_phantomjs()
taxonomy_table <- tibble()

for (i in row_indices_t7){
  loc <-  contamdf.prev.01_t7[i, 0]
  tax_key <- row.names(loc)
  tax_value <- tax.fil[tax_key, ]
  taxonomy_table <- rbind(taxonomy_table, tax_value)
}

datatable(taxonomy_table) 

# save the contaminant table in the local computer
write.csv(df.pa_t7, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t7-neg.csv")
#removing contaminants from phyloseq object
colnames(df.pa_t7)
df.pa_t7 <- rownames_to_column(df.pa_t7, var = "ASV.ID")
subset.df_t7 <- subset(df.pa_t7, contaminant== "FALSE")
keep.taxa_t7 <- as.vector(subset.df_t7$ASV.ID)
length(keep.taxa_t7)

subset.df.remove_t7 <- subset(df.pa_t7, contaminant== "TRUE")
remove.taxa.list_t7 <- as.vector(subset.df.remove_t7)
remove.taxa.df_t7 <- data.frame(remove.taxa.list_t7)

decontam_phyloseq_t7 <- prune_taxa(keep.taxa_t7, physeq_t7)
decontam_phyloseq_t7
#filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t7 <- subset(sample_data(decontam_phyloseq_t7), is.neg ==FALSE)
keep.samples_t7 <- as.vector(subset.metadata_t7$SampleID)
keep.samples_t7
decontam_phyloseq_t7 <- prune_samples(keep.samples_t7, decontam_phyloseq_t7)

### decontam with positive control

sample_data(decontam_phyloseq_t7)$is.neg <- sample_data(decontam_phyloseq_t7)$sample_type == "positive"
contamdf.prev0.1_t7pos <- isContaminant(decontam_phyloseq_t7, method="prevalence", neg="is.neg")
table(contamdf.prev0.1_t7pos$contaminant)# FALSE=1471  TRUE=17
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t7pos <- transform_sample_counts(decontam_phyloseq_t7, function(abund) 1*(abund>0))
ps.pa.neg_t7pos <- prune_samples(sample_data(ps.pa_t7pos)$is.neg == TRUE, ps.pa_t7pos)
ps.pa.pos_t7pos <- prune_samples(sample_data(ps.pa_t7pos)$is.neg == FALSE, ps.pa_t7pos)
# Make data.frame of prevalence in positive and negative samples
df.pa_t7pos <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t7pos), pa.neg=taxa_sums(ps.pa.neg_t7pos),
                          contaminant=contamdf.prev0.1_t7pos$contaminant)
t7pos_decontam_plot <- ggplot(data=df.pa_t7pos, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t7 Decontam")
t7pos_decontam_plot
# save the contaminant table in the local computer
write.csv(df.pa_t7pos, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t7-pos.csv")
##removing contaminants from phyloseq object
colnames(df.pa_t7pos)
df.pa_t7pos <- rownames_to_column(df.pa_t7pos, var = "ASV.ID")
subset.df_t7pos <- subset(df.pa_t7pos, contaminant== "FALSE")
keep.taxa_t7pos <- as.vector(subset.df_t7pos$ASV.ID)
length(keep.taxa_t7pos)

subset.df.remove_t7pos <- subset(df.pa_t7pos, contaminant== "TRUE")
remove.taxa.list_t7pos <- as.vector(subset.df.remove_t7pos)
remove.taxa.df_t7pos <- data.frame(remove.taxa.list_t7pos)

decontam_phyloseq_t7 <- prune_taxa(keep.taxa_t7pos, decontam_phyloseq_t7)
decontam_phyloseq_t7
##filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t7 <- subset(sample_data(decontam_phyloseq_t7), is.neg ==FALSE)
keep.samples_t7 <- as.vector(subset.metadata_t7$SampleID)
keep.samples_t7
decontam_phyloseq_t7 <- prune_samples(keep.samples_t7, decontam_phyloseq_t7)
otu_table(decontam_phyloseq_t7)
#save
write.csv(otu_table(decontam_phyloseq_t7), "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/decontam_ASV_t7pos.csv")

#______________________________________________________________________________________________________________________________________________________________________


# 8. Batch T8


physeq_t8 <- physeq.decon %>% 
  ps_filter(batch == "8", .keep_all_taxa = TRUE)
sample_data(physeq_t8)
physeq_t8
# default threshold (0,1)
sample_data(physeq_t8)$is.neg <- sample_data(physeq_t8)$sample_type == "negative"
contamdf.prev.01_t8 <- isContaminant(physeq_t8, method="prevalence", neg="is.neg")
table(contamdf.prev.01_t8$contaminant) # FALSE: 1490 TRUE: 9
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t8 <- transform_sample_counts(physeq_t8, function(abund) 1*(abund>0))
ps.pa.neg_t8 <- prune_samples(sample_data(ps.pa_t8)$is.neg == TRUE, ps.pa_t8)
ps.pa.pos_t8 <- prune_samples(sample_data(ps.pa_t8)$is.neg == FALSE, ps.pa_t8)
# Make data.frame of prevalence in positive and negative samples
df.pa_t8 <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t8), pa.neg=taxa_sums(ps.pa.neg_t8),
                       contaminant=contamdf.prev.01_t8$contaminant)
df.pa_t8
t8_decontam_plot <- ggplot(data=df.pa_t8, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t8 Decontam")
t8_decontam_plot
# grab the contaminant rows
row_indices_t8 <- which(contamdf.prev.01_t8$contaminant)
row_indices_t8
# extract the taxonomic information of the contaminants
#webshot::install_phantomjs()
taxonomy_table <- tibble()

for (i in row_indices_t8){
  loc <-  contamdf.prev.01_t8[i, 0]
  tax_key <- row.names(loc)
  tax_value <- tax.fil[tax_key, ]
  taxonomy_table <- rbind(taxonomy_table, tax_value)
}

datatable(taxonomy_table) 

# save the contaminant table in the local computer
write.csv(df.pa_t8, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t8-neg.csv")
#removing contaminants from phyloseq object
colnames(df.pa_t8)
df.pa_t8 <- rownames_to_column(df.pa_t8, var = "ASV.ID")
subset.df_t8 <- subset(df.pa_t8, contaminant== "FALSE")
keep.taxa_t8 <- as.vector(subset.df_t8$ASV.ID)
length(keep.taxa_t8)

subset.df.remove_t8 <- subset(df.pa_t8, contaminant== "TRUE")
remove.taxa.list_t8 <- as.vector(subset.df.remove_t8)
remove.taxa.df_t8 <- data.frame(remove.taxa.list_t8)

decontam_phyloseq_t8 <- prune_taxa(keep.taxa_t8, physeq_t8)
decontam_phyloseq_t8
#filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t8 <- subset(sample_data(decontam_phyloseq_t8), is.neg ==FALSE)
keep.samples_t8 <- as.vector(subset.metadata_t8$SampleID)
keep.samples_t8
decontam_phyloseq_t8 <- prune_samples(keep.samples_t8, decontam_phyloseq_t8)

### decontam with positive control

sample_data(decontam_phyloseq_t8)$is.neg <- sample_data(decontam_phyloseq_t8)$sample_type == "positive"
contamdf.prev0.1_t8pos <- isContaminant(decontam_phyloseq_t8, method="prevalence", neg="is.neg")
table(contamdf.prev0.1_t8pos$contaminant)# FALSE=1490  TRUE=0
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_t8pos <- transform_sample_counts(decontam_phyloseq_t8, function(abund) 1*(abund>0))
ps.pa.neg_t8pos <- prune_samples(sample_data(ps.pa_t8pos)$is.neg == TRUE, ps.pa_t8pos)
ps.pa.pos_t8pos <- prune_samples(sample_data(ps.pa_t8pos)$is.neg == FALSE, ps.pa_t8pos)
# Make data.frame of prevalence in positive and negative samples
df.pa_t8pos <- data.frame(pa.pos=taxa_sums(ps.pa.pos_t8pos), pa.neg=taxa_sums(ps.pa.neg_t8pos),
                          contaminant=contamdf.prev0.1_t8pos$contaminant)
t8pos_decontam_plot <- ggplot(data=df.pa_t8pos, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + labs(x="Prevalence (Negative Controls)", y="Prevalence (True Samples)", title="t8 Decontam")
t8pos_decontam_plot
# save the contaminant table in the local computer
write.csv(df.pa_t8pos, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/contaminant-table-t8-pos.csv")
##removing contaminants from phyloseq object
colnames(df.pa_t8pos)
df.pa_t8pos <- rownames_to_column(df.pa_t8pos, var = "ASV.ID")
subset.df_t8pos <- subset(df.pa_t8pos, contaminant== "FALSE")
keep.taxa_t8pos <- as.vector(subset.df_t8pos$ASV.ID)
length(keep.taxa_t8pos)

subset.df.remove_t8pos <- subset(df.pa_t8pos, contaminant== "TRUE")
remove.taxa.list_t8pos <- as.vector(subset.df.remove_t8pos)
remove.taxa.df_t8pos <- data.frame(remove.taxa.list_t8pos)

decontam_phyloseq_t8 <- prune_taxa(keep.taxa_t8pos, decontam_phyloseq_t8)
decontam_phyloseq_t8
##filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
subset.metadata_t8 <- subset(sample_data(decontam_phyloseq_t8), is.neg ==FALSE)
keep.samples_t8 <- as.vector(subset.metadata_t8$SampleID)
keep.samples_t8
decontam_phyloseq_t8 <- prune_samples(keep.samples_t8, decontam_phyloseq_t8)
otu_table(decontam_phyloseq_t8)
#save
write.csv(otu_table(decontam_phyloseq_t8), "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/decontam_ASV_t8pos.csv")

## Merge decontaminated phyloseq object together
phyloseq_decontam_all <- merge_phyloseq(decontam_phyloseq_t1, decontam_phyloseq_t2, decontam_phyloseq_t3, decontam_phyloseq_t4,
                                    decontam_phyloseq_t5, decontam_phyloseq_t6, decontam_phyloseq_t7, decontam_phyloseq_t8)
colnames(otu_table(phyloseq_decontam_all))
# remove ASVs/taxa with zero reads
phyloseq_decontam_all.1 <- prune_taxa(taxa_sums(phyloseq_decontam_all)>=1, phyloseq_decontam_all)
phyloseq_decontam_all.1 # 806 decontaminated taxa 63 true samples
sum((otu_table(phyloseq_decontam_all.1)))
# checking
sort(colSums(otu_table(phyloseq_decontam_all.1), na.rm = FALSE, dims = 1), decreasing = F)

saveRDS(phyloseq_decontam_all.1, file="/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/phyloseq_decontam_all1.rds", compress = TRUE)

# Remove potential contaminant ASV/taxa Phyromonas

# Read in decontaminated phyloseq object
decontaminated_physeq <- readRDS(file="/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/phyloseq_decontam_all1.rds")
decontaminated_physeq
sample_data(decontaminated_physeq)

decontam_asv_table <- as.data.frame(otu_table(decontaminated_physeq))
head(decontam_asv_table)
decontam_asv_table <- rownames_to_column(decontam_asv_table, var = "ASV.ID")
decontam_tax_table <- as.data.frame(tax_table(decontaminated_physeq))
decontam_tax_table <- rownames_to_column(decontam_tax_table, var = "ASV.ID")
head(decontam_tax_table)
decontam_asv_tax <- merge(decontam_asv_table, decontam_tax_table, by="ASV.ID", all = T)
#save to manually check
write.csv(decontam_asv_tax, "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/decontam_asv_tax.csv")

## Build a new phyloseq object with the phylogeny and new metadata

# Adding the phylogeny information (rooted tree)
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07')
rooted_tree <- ape::read.tree("rooted-tree.nwk")
rooted_tree
# Adding metadata
meta.rainout <- read.csv("metadata_rainout.csv")
view(meta.rainout)
meta.rainout # still contain sample number 402
# remove the sample 402 from the meta data
meta.rain.ed <- subset(meta.rainout, SampleID!= 402)
view(meta.rain.ed)
meta.rain.ed <- meta.rain.ed %>% 
  remove_rownames %>% 
  column_to_rownames(var="SampleID")
# match the order of the column names of the decontaminated ASV table with the row names of the metadata (sample identity)
head(decontam_asv_table)
decontam_asv_table <- column_to_rownames(decontam_asv_table, var = "ASV.ID")
dim(decontam_asv_table)
re_order <- match(rownames(meta.rain.ed), colnames(decontam_asv_table))
decontam_asv_table.ord  <- decontam_asv_table[ ,re_order]
# make phyloseq object
decontam.asv.seq <- otu_table(decontam_asv_table.ord, taxa_are_rows = TRUE) # asv table
head(decontam_tax_table)
decontam_tax_table <- column_to_rownames(decontam_tax_table, var = "ASV.ID")
decontam.tax.seq <- tax_table(as.matrix(decontam_tax_table))
meta.rain.ed <- rownames_to_column(meta.rain.ed, var = "SampleID")
meta.rain.ed$SampleID <- as.factor(meta.rain.ed$SampleID)
rownames(meta.rain.ed) <- meta.rain.ed$SampleID
meta.rain.seq <- sample_data(meta.rain.ed)
# merge phyloseq object
decontaminated_physeq.ed <- merge_phyloseq(decontam.asv.seq, decontam.tax.seq, meta.rain.seq, rooted_tree)
sample_data(decontaminated_physeq.ed)
# 1. unassigned, plant, and bacterial contaminants-filtered phyloseq object (not rarefied and not normalized)
saveRDS(decontaminated_physeq.ed, file="/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/phyloseqobj/decontaminated_unrare_notnormal_physeq.rds", compress = TRUE)
write.csv(otu_table(decontaminated_physeq.ed), "/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/decontam_ASV_with_phyro.csv")


## Removing Porphyromonas ASV

head(decontam_tax_table)
tax_df <- rownames_to_column(decontam_tax_table, var="ASV.ID")
tax_df$Keep <- tax_df$ASV.ID != "0232b56c81bb04fb94bdcfdf8b073ae1"
table(tax_df$Keep)
tax_df_keep <- filter(tax_df, tax_df$Keep == TRUE)
tax_list <- tax_df_keep$ASV.ID
decontaminated_physeq.ed2 <- prune_taxa(tax_list, decontaminated_physeq.ed)
decontaminated_physeq.ed2 # 805 taxa left and 63 sample
sort(colSums(otu_table(decontaminated_physeq.ed2), na.rm = FALSE, dims = 1), decreasing = F)
sum(otu_table(decontaminated_physeq.ed2))
sample_data(decontaminated_physeq.ed2)

## Save in the local computer:

# 1. unassigned, plant, and bacterial contaminants-filtered phyloseq object (not rarefied and not normalized)
saveRDS(decontaminated_physeq.ed2, file="/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/phyloseqobj/decontaminated_unrare_notnormal_physeq.rds", compress = TRUE)
# 2. unassigned, plant, and bacterial contaminants-filtered phyloseq object (rarefied to 1165 reads (remove sample with reads < 1000))
saveRDS(decon.rare.1165.seq, file="/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/phyloseqobj/decontaminated_rare_physeq.rds", compress = TRUE)
# 3. unassigned, plant, and bacterial contaminants-filtered phyloseq object (normalized using CSS method)
saveRDS(normalized.physeq, file="/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/phyloseqobj/decontaminated_normal_physeq.rds", compress = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# checking with the phyloseq object before decontamination
sample_data(physeq.decon)
Samples_toRemove <- c("Mock1", "Mock2", "Mock3", "Mock4", 
                      "Mock5", "Mock6", "Mock7", "Mock8", 
                      "NegCon1", "NegCon2", "NegCon3", "NegCon4",
                      "NegCon5", "NegCon6", "NegCon7", "NegCon8")
#To see what samples get removed, run the following; note, I have a column called "SampleID"
subset <- subset_samples(physeq.decon, SampleID %in% Samples_toRemove)
subset
#To remove those from your phyloseq object
physeq.decon.bio <- subset_samples(physeq.decon, !(SampleID %in% Samples_toRemove)) #This will return a ps object with the samples removed
physeq.decon.bio #only true samples
# checking
sort(rowSums(otu_table(physeq.decon.bio), na.rm = FALSE, dims = 1), decreasing = F) # contain zero reads ASVs
# remove ASVs/taxa with zero reads
physeq.decon.bio.1 <- prune_taxa(taxa_sums(physeq.decon.bio)>=1, physeq.decon.bio)
physeq.decon.bio.1 # 813 taxa, 63 true samples
# there are 7 ASVs (+1 Phyromonas ASV) that are removed after decontam
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


###################################################################################################
# rarefaction to 1165 reads (remove sample with reads < 1000)
###################################################################################################
# ASV Table
set.seed(333)
decon.rare.1165.seq <- rarefy_even_depth(decontaminated_physeq.ed2, sample.size = 1165,
                                       rngseed = 333, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
decon.rare.1165.seq # 1 samples removed (4001), 137 ASVs were removed, 668 ASVs and 62 samples remaining
sort(rowSums(otu_table(decon.rare.1165.seq), na.rm = FALSE, dims = 1), decreasing = F)
set.seed(13)
ggrare(decon.rare.1165.seq, step = 1, color = "treatment", label = "SampleID", se = FALSE)

# rarefaction plot

# run the ggrare function attached in the file "generating_rarecurve.r"
set.seed(13)
decontam.rare <- ggrare(phyloseq_decontam_all.1, step = 1,se = FALSE) #color = "treatment", label = "SampleID", se = FALSE)

#set up your own color palette
#install.packages("colorBlindness")
library(colorBlindness)
displayAvailablePalette(color="white")
#Palette <- c("#1F968BFF","#FDE725FF")
PairedColor12Steps
Brown2Blue10Steps
Palette <- c("#FF7F00", "#662F00")
names(Palette) <- levels(sample_data(aob.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.rare <- decontam.rare + 
  theme_bw()+
  scale_size_manual(values = 60)+
  xlim(0, 7500)+
  theme( strip.text.x = element_text(size=14, face='bold'),
         axis.text.x=element_text(size = 14),
         axis.text.y = element_text(size = 14),
         strip.text.y = element_text(size=18, face = 'bold'),
         plot.title = element_text(size =20 ,face='bold'),
         axis.title.y = element_text(size=15,face="bold"),
         axis.title.x = element_text(size=15,face="bold"),
         legend.position = "right",
         legend.title = element_text(size=15),
         legend.text = element_text(size = 13),
         plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
  ylab("Number of ASVs")+xlab("Reads")
plot.rare
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')

ggsave("rarecurve.tiff",
       plot.rare, device = "tiff",
       width = 10, height = 7, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")

############################################################################################################################################################################## 

# Reads normalization using Cumulative Sum Scaling (CSS) method in metagenomeSeq package

#install and load metagenomeSeq package
#BiocManager::install("metagenomeSeq")
library(metagenomeSeq)

#loading metadata
meta.df <- as.data.frame(sample_data(decontaminated_physeq.ed2))
view(meta.df)
# create MRexperiment object of the metadata
meta.df.ed <- meta.df %>% 
  remove_rownames %>% 
  column_to_rownames(var="SampleID")
phenotypeData <- AnnotatedDataFrame(meta.df.ed)
#loading taxonomy
tax.df <- as.data.frame(tax_table(decontaminated_physeq.ed2))
head(tax.df)
# create MRexperiment object of the tax table
TAXdata <- AnnotatedDataFrame(tax.df)
# loading ASV table
asv.df <- as.data.frame(otu_table(decontaminated_physeq.ed2))
# Creating MRexperiment object
model.obj <- newMRexperiment(asv.df, phenoData = phenotypeData, featureData = TAXdata)
model.obj
#sort(colSums(MRcounts(model.obj), na.rm = FALSE, dims = 1))

## Normalization

# 1. calculate the proper percentile by which to normalize counts
p <- cumNormStat(model.obj)
p
# 2. To calculate the scaling factors we simply run cumNorm
     #cumNorm=Calculates each column’s quantile and calculates the sum up to and including that quantile.
obj = cumNorm(model.obj, p = cumNormStat(model.obj))
normFactors(obj)
# 3. Exporting normalized count matrix
normalized.asv <- MRcounts(obj, norm = T, log = F)
normalized.asv <- as.data.frame(normalized.asv)
head(normalized.asv)
#sort(rowSums(normalized.asv, na.rm = FALSE, dims = 1), decreasing = F)
dim(normalized.asv) #there are 805 ASVs and 63 samples
# Make phyloseq object for the normalized ASV table
norm.seq <- otu_table(normalized.asv, taxa_are_rows = TRUE)
norm.seq
# merge phyloseq object
normalized.physeq <- merge_phyloseq(norm.seq, decontam.tax.seq, meta.rain.seq, rooted_tree)
normalized.physeq
otu_table(normalized.physeq)






# dissimilarity indices for community ecologist to make a distance structure (Jaccard distance between samples)
head(normalized.otu)
normalized.otu <- column_to_rownames(normalized.otu, var="OTU.ID")
normalized.otu_PA <- 1*(normalized.otu>0)
normalized.otu_PA
otu_dist <- vegdist(t(normalized.otu_PA), binary = TRUE, method = "bray") #Sorensen
otu_dist

set.seed(13)
permu_scheme <- how(within = Within(type = "free"),
                    plots = Plots(strata=beta.rainmap$var, type = "none"),
                    blocks = beta.rainmap$rep_block,
                    nperm = 999,
                    observed = TRUE)
perm.otu = adonis2(otu_dist ~ var+trt+loc,
                   permutations = permu_scheme, 
                   method="bray", 
                   by ="margin",
                   data = beta.rainmap)
perm.otu






