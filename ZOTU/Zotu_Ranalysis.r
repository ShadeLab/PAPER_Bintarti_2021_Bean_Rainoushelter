#######################################################################################################################
############### Bean seed microbiome analysis for the rain out shelter experiment: ZOTU #####################################
#######################################################################################################################
# Date: July 26th 2021
# By : Ari Fina Bintarti
# INSTALL PACKAGES
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
install.packages("ggpubr")
install.packages("car")
install.packages("agricolae")
install.packages("multcompView")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("sjmisc") 
install.packages("sjPlot")
install.packages("MASS")
install.packages("FSA")
install.packages('mvtnorm', dep = TRUE)
install.packages("rcompanion")
install.packages("onewaytests")
install.packages("PerformanceAnalytics")
install.packages("gvlma")
install.packages("userfriendlyscience")
install.packages("ggpmisc")
install.packages("fitdistrplus")
install.packages('BiocManager')
#install.packages("cowplot")
install.packages("dplyr")
install.packages("lme4")
install.packages("nlme")
install.packages("car")
install.packages("multcomp")
library(multcomp)
library(car)
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
#library(cowplot)
library(ggplot2)
library(reshape)
library(ggpubr)
library(car)
library(agricolae)
library(multcompView)
library(grid)
library(gridExtra)
library(sjmisc)
library(sjPlot)
library(MASS)
library(FSA)
library(rcompanion)
library(onewaytests)
library(ggsignif)
library(PerformanceAnalytics)
library(gvlma)
library(userfriendlyscience)
library(ggpmisc)
library(tibble)
library(fitdistrplus)
library(lme4)
library(nlme)

# SET THE WORKING DIRECTORY

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/ZOTU/20210601_16SV4')
wd <- print(getwd())

# READ PROPORTION OF CHLOROPLAST AND MITOCHONDRIA

#read the unfiltered zotu table 
zotu.unfil <- read.table(file = 'ZOTU_table_tax.tsv', sep = '\t', header = TRUE,check.names = FALSE)
zotu.unfil
tax.unfil <- zotu.unfil[,'taxonomy']
tax.unfil

#write.csv(tax.unfil, file = "tax.unfil.csv")
dim(zotu.unfil)
zotu.unfil <- zotu.unfil[,-81]
dim(zotu.unfil) # zotu= 122, zotu table still has Mock, NC, and PC in the sample
zotu.unfil <- column_to_rownames(zotu.unfil,var = "ZOTUID")
sort(rowSums(zotu.unfil, na.rm = FALSE, dims = 1), decreasing = F)

#read taxonomy
tax.unfil.ed = read.csv("tax.unfil.ed.csv", header=T)
rownames(tax.unfil.ed) <- rownames(zotu.unfil)
dim(tax.unfil.ed)

zotu.unfil <- rownames_to_column(zotu.unfil,var = "ZOTUID")
tax.unfil.ed <- rownames_to_column(tax.unfil.ed,var = "ZOTUID")
zotu.tax.unfiltered <- merge(zotu.unfil, tax.unfil.ed, by="ZOTUID")
write.csv(zotu.tax.unfiltered, file = "zotu.tax.unfiltered.csv")
#read the metadata


#select only biological sample from zotu table
zotu.bio.unfil <- zotu.unfil[,1:64] #unselect Mock, NC, and PC from the zotu table
dim(zotu.bio.unfil)
sort(rowSums(zotu.bio.unfil, na.rm = FALSE, dims = 1), decreasing = F)

#otu table of the negative control
NC.unfiltered <- zotu.unfil[,c(72:78)]#only negative control
NC.unfiltered
#NC.unfiltered <- column_to_rownames(NC.unfiltered,var="ZOTUID")
sort(rowSums(NC.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
NC1.unfiltered=NC.unfiltered[which(rowSums(NC.unfiltered) > 0),]
sort(rowSums(NC1.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
NC1.unfiltered <- rownames_to_column(NC1.unfiltered,var="ZOTUID")
NC1.tax.unfiltered <- merge(NC1.unfiltered, tax.unfil.ed, by="ZOTUID")
write.csv(NC1.tax.unfiltered, file = "NC1.tax.unfiltered.csv")


# remove zOTUs that do not present in sample
zotu.bio1.unfil <- zotu.bio.unfil[which(rowSums(zotu.bio.unfil) > 0),]
dim(zotu.bio1.unfil) # zotu= 99, zotu table before normalization using metagenomeSeq package and before decontamination
sort(rowSums(zotu.bio1.unfil, na.rm = FALSE, dims = 1), decreasing = F)

# load the zotu table
head(zotu.bio1.unfil)
zotu.bio1.unfil <- rownames_to_column(zotu.bio1.unfil, var = "ZOTUID")

# merge the taxonomy with otu table
head(tax.unfil.ed)
tax.unfil.ed <- rownames_to_column(tax.unfil.ed, var = "ZOTUID")
zotu.tax.unfil <- merge(zotu.bio1.unfil, tax.unfil.ed, by="ZOTUID")
dim(zotu.tax.unfil)

#select only the otu table and "Order"  & "Family"
#otu.tax.unfil.ed <- otu.tax.unfil[,c(1:48,52,53)]
#colnames(otu.tax.unfil.ed)

#edit the taxonomy
zotu.tax.unfil.ed <- zotu.tax.unfil %>%
    mutate(Taxonomy = case_when(Family == "Chloroplast" ~ 'Chloroplast',
                                  Family == "Mitochondria" ~ 'Mitochondria',
                                  Family == "Magnoliophyta" ~ 'Magnoliophyta',
                                  TRUE ~ 'Bacteria')) %>%
    mutate(Kingdom = case_when(Family == "Chloroplast" ~ 'Plant',
                                  Family == "Mitochondria" ~ 'Plant',
                                  Family == "Magnoliophyta" ~ 'Plant',
                                  TRUE ~ 'Bacteria'))

tail(zotu.tax.unfil.ed)
zotu.tax.unfil.ed1 <- zotu.tax.unfil.ed[,c(1:65,73:74)]
colnames(zotu.tax.unfil.ed1)
tail(zotu.tax.unfil.ed1)

long.dat <- gather(zotu.tax.unfil.ed1, Sample, Read, 2:65, factor_key = T)
long.dat


df.unfil <- long.dat %>%
  group_by(Sample, Kingdom) %>%
  summarize(read.number = sum(Read))
df.unfil1 <- df.unfil %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)

with(df.unfil1, sum(percent[Sample ==  "1001"]))

library(ggbeeswarm)
library(ggtext)
plot.unfil.king <- ggplot(df.unfil1, aes(x=Kingdom, y=percent, fill=Kingdom))+
                    geom_violin(trim = F, scale="width") +
                    #geom_beeswarm(dodge.width = 1, alpha = 0.3)+
                    scale_fill_manual(labels = c("Bacteria","Plant"),values=c("#88CCEE", "#117733"))+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    #geom_text(data=sum_rich_plant_new, aes(x=Plant,y=2+max.rich,label=Letter), vjust=0)+
                    labs(title = "A")+
                    ylab("Read Proportion (%)")+
                    theme(legend.position="none",
                          #axis.text.x=element_blank(),
                          #axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text= element_text(size = 14),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 14, face = 'bold'),
                          #axis.title.y=element_text(size=13,face="bold"),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          #plot.margin = unit(c(0, 0, 0, 0), "cm"))
                          stat_summary(fun="median",geom="point", size=7, color="red", shape=95)
                          #width=1, position=position_dodge(),show.legend = FALSE)

plot.unfil.king
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/ZOTU/20210601_16SV4/Figures')
ggsave("20210601_plant_proportion.eps",
       plot.unfil.king, device=cairo_ps,
       width = 5, height =5, 
       units= "in", dpi = 600)


df.unfil.tax <- long.dat %>%
  group_by(Sample, Taxonomy) %>%
  summarize(read.number = sum(Read))
df.unfil.tax1 <- df.unfil.tax %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)

plot.unfil.tax <- ggplot(df.unfil.tax1, aes(x=Taxonomy, y=percent, fill=Taxonomy))+
                    geom_violin(trim = F, scale="width") +
                    #geom_beeswarm(dodge.width = 1, alpha = 0.3)+
                    #scale_fill_manual(labels = c("A1","A2", "A3","B1","B2","B3","B4","B5","B6","C5","C6","C7"),values=c("#440154FF", "#482677FF","#3F4788FF","#238A8DFF","#1F968BFF","#20A386FF","#29AF7FF","#3CBC75F","#56C667FF","#B8DE29FF","#DCE318FF","#FDE725FF"))+
                    #scale_fill_viridis(discrete = T)+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    #geom_text(data=sum_rich_plant_new, aes(x=Plant,y=2+max.rich,label=Letter), vjust=0)+
                    labs(title = "A")+
                    ylab("Read Proportion (%)")+
                    theme(legend.position="none",
                          #axis.text.x=element_blank(),
                          #axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text= element_text(size = 14),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 14, face = 'bold'),
                          #axis.title.y=element_text(size=13,face="bold"),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          #plot.margin = unit(c(0, 0, 0, 0), "cm"))
                          stat_summary(fun="median",geom="point", size=7, color="red", shape=95)
                          #width=1, position=position_dodge(),show.legend = FALSE)

plot.unfil.tax
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/ZOTU/20210601_16SV4/Figures')
ggsave("20210601_chloromito_proportion.eps",
       plot.unfil.tax, device=cairo_ps,
       width = 7, height =5, 
       units= "in", dpi = 600)


# ANALYSIS OF READS AFTER CHLOROPLAST AND MITOCHONDRIA REMOVAL

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/ZOTU/20210601_16SV4')
zotu <- read.table('ZOTUtable-no-mitoch-no-chloro_tax.tsv', sep='\t', header=T, row.names = 1, check.names = FALSE)
zotu
tax <- zotu[,'taxonomy']
str(tax)
#write.csv(tax, file = "tax.fil.csv")
dim(zotu)
zotu <- zotu[,-80]
dim(zotu) # otu= 103, otu table still has Mock, NC, and PC in the sample
sort(rowSums(zotu, na.rm = FALSE, dims = 1), decreasing = F)

#read taxonomy
tax.ed = read.csv("tax.fil.ed.csv", header=T)
rownames(tax.ed) <- rownames(zotu)
dim(tax.ed)
#read the metadata

#select only biological sample from zotu table
zotu.bio <- zotu[,1:64] #unselect Mock, NC, and PC from the zotu table
dim(zotu.bio)
zotu.bio
sort(rowSums(zotu.bio, na.rm = FALSE, dims = 1), decreasing = F)
# remove zOTUs that do not present in sample
zotu.bio1=zotu.bio[which(rowSums(zotu.bio) > 0),]
dim(zotu.bio1) # otu= 80, zotu table before normalization using metagenomeSeq package and before decontamination
sort(rowSums(zotu.bio1, na.rm = FALSE, dims = 1), decreasing = F)
# merge zotu.bio1 with taxonomy to have match taxonomy table
head(zotu.bio1)
zotu.bio1 <- rownames_to_column(zotu.bio1,var = "ZOTUID")
#tax.ed <- rownames_to_column(tax.ed,var = "ZOTUID")
zotu.bio1.tax <- merge(zotu.bio1, tax.ed, by="ZOTUID")
dim(zotu.bio1.tax) 
#  separate  the sample 
# otu table
zotu.bac.fil <- zotu.bio1.tax[,c(1:65)]
head(zotu.bac.fil)
zotu.bac.fil <- column_to_rownames(zotu.bac.fil,var="ZOTUID")
sum(zotu.bac.fil)
dim(zotu.bac.fil)

#otu table of the negative control
NC <- zotu[,c(72:78)]#only negative control
NC
NC <- column_to_rownames(NC,var="ZOTUID")
sort(rowSums(NC, na.rm = FALSE, dims = 1), decreasing = F)
NC1=NC[which(rowSums(NC) > 0),]
sort(rowSums(NC1, na.rm = FALSE, dims = 1), decreasing = F)
NC1 <- rownames_to_column(NC1,var="ZOTUID")
NC1.tax <- merge(NC1, tax.ed, by="ZOTUID")
write.csv(NC1.tax, file = "NC1.tax.csv")

######################################################################################################################################
######################################################################################################################################

############## remove contaminant reads/otus from otu bac using microDecon package ##################

#install and load microDecon package
#install.packages("devtools")
#devtools::install_github("donaldtmcknight/microDecon")
library(microDecon)

#merge otu.NC to otu bac and put otu.NC as the first column
head(NC)
NC <- rownames_to_column(NC, var = "ZOTUID")
head(zotu.bio)
zotu.bio <- rownames_to_column(zotu.bio, var = "ZOTUID")
zotu.NC.bac <- merge(NC, zotu.bio, by="ZOTUID")
dim(zotu.NC.bac)

#decontamination
zotu.decon <- decon(data = zotu.NC.bac, numb.blanks = 7, numb.ind = 64, taxa = F)
zotu.decon$OTUs.removed # remove 49 OTUs
zotu.tab <- as.data.frame(zotu.decon$decon.table)

#remove NC from the otu.dec.table 
zotu.tab <- zotu.tab[,c(1,3:66)]
zotu.tab
dim(zotu.tab) #there are 54 otus and 64 samples

#merge otu.dec.table with taxonomy to have match taxonomy table
head(tax.ed)
head(zotu.tab)
#taxonomy <- rownames_to_column(taxonomy, var = "OTU.ID")
zotu.tab.tax <- merge(zotu.tab, tax.ed, by = "ZOTUID")
zotu.tab.tax

#separate the otu table and taxonomy table
#otu table
decon.zotu <- zotu.tab.tax[,c(1:65)]
dim(decon.zotu) #54 otus, 47samples
tax <- zotu.tab.tax[,c(1,66:72)]
dim(tax) #54 otus, 7columns/levels
head(decon.zotu)
decon.zotu <- column_to_rownames(decon.zotu, var = "ZOTUID")
#sort(colSums(decon.zotu, na.rm = FALSE, dims = 1), decreasing = F)
sort(rowSums(decon.zotu, na.rm = FALSE, dims = 1), decreasing = F)

#remove ZOTUs that do not present in sample: 
decon.zotu1=decon.zotu[which(rowSums(decon.zotu) > 0),]
sort(rowSums(decon.zotu1, na.rm = FALSE, dims = 1), decreasing = F)
#dim(otu1) #211 otus, 47samples, clean otu table after decontamination before reads normalization

######################################################################################################################################
######################################################################################################################################

############## Reads normalization using metagnomeSeq ####################

#install and load metagenomeSeq package
#BiocManager::install("metagenomeSeq", dependencies=TRUE)
library(metagenomeSeq)

#Creating a MRexperiment object

#loading count data (otu)
zotu.bac.fil
t(zotu.bac.fil)
#loading metadata


#create MRexperiment object 
phenotypeData <- AnnotatedDataFrame(map)
phenotypeData

OTUdata <- AnnotatedDataFrame(tax)
OTUdata

#create model
model <- newMRexperiment(decon.otu, phenoData = phenotypeData, featureData = OTUdata)
model

sort(rowSums(MRcounts(model), na.rm = FALSE, dims = 1), decreasing = T)

#normalising the data

#normalise the data to account for differences due to uneven sequencing depth
#metagenomeSeq uses Cumulative Sum Scaling (CSS) normalisation instead of rarefaction
#cumNormStatFast=Calculates the percentile for which to sum counts up to and scale by.

p <- cumNormStatFast(model, pFlag = T)

#cumNorm=Calculates each column’s quantile and calculates the sum up to and including that quantile.

bac.norm <- cumNorm(model, p = p)
bac.norm

#export count matrix
otu.norm <- MRcounts(bac.norm, norm = TRUE, log = F)

otu.norm <- as.data.frame(otu.norm)

head(sort(colSums(otu.norm, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums(otu.norm, na.rm = FALSE, dims = 1), decreasing = FALSE))

head(otu.norm)
otu.norm <- rownames_to_column(otu.norm, var="OTU.ID")
#write.table(otu.norm, "otu_norm.txt", sep = "\t", quote = F, row.names = F)
dim(otu.norm) #there are 211 otus and 47 samples




