#######################################################################################################################
############### Bean seed microbiome analysis for the rain out shelter experiment: OTU 97% ############################
#######################################################################################################################
# Date: August 18th 2021
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
library(dplyr)
library(plyr)
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

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
wd <- print(getwd())

# READ PROPORTION OF CHLOROPLAST AND MITOCHONDRIA

#read the unfiltered otu table 
otu.unfil <- read.table(file = 'OTU_table_tax.txt', sep = '\t', header = TRUE,check.names = FALSE)
otu.unfil
tax.unfil <- otu.unfil[,'taxonomy']
tax.unfil

#write.csv(tax.unfil, file = "tax.unfil.csv")
dim(otu.unfil) #[1] 325  81
colnames(otu.unfil)
otu.unfil <- otu.unfil[,-82]
dim(otu.unfil)# otu= 325, otu table still has Mock, NC, and PC in the sample
otu.unfil <- column_to_rownames(otu.unfil,var = "OTUID")
sort(rowSums(otu.unfil, na.rm = FALSE, dims = 1), decreasing = F)

#read taxonomy
tax.unfil.ed = read.csv("tax.unfil.ed.csv", header=T)
rownames(tax.unfil.ed) <- rownames(otu.unfil)
dim(tax.unfil.ed) #[1] 325   7

otu.unfil <- rownames_to_column(otu.unfil,var = "OTUID")
tax.unfil.ed <- rownames_to_column(tax.unfil.ed,var = "OTUID")
otu.tax.unfiltered <- merge(otu.unfil, tax.unfil.ed, by="OTUID")
View(otu.tax.unfiltered)
colnames(otu.tax.unfiltered)
#write.csv(otu.tax.unfiltered, file = "otu.tax.unfiltered.csv")
#read the metadata

#############################################################################################################################################################
#READ PROPORTION OF CHLOROPLAST AND MITOCHONDRIA OF EXPERIMENTAL SAMPLES

#select only biological sample from otu table
otu.bio.unfil <- otu.unfil[,1:65] #unselect Mock, NC, and PC from the otu table
dim(otu.bio.unfil)
colnames(otu.bio.unfil)
otu.bio.unfil <- column_to_rownames(otu.bio.unfil, var = "OTUID")
sort(rowSums(otu.bio.unfil, na.rm = FALSE, dims = 1), decreasing = F)
# remove OTUs that do not present in biological sample
otu.bio1.unfil <- otu.bio.unfil[which(rowSums(otu.bio.unfil) > 0),]
dim(otu.bio1.unfil) # [1] 244  64, otu table before plant contaminant removal  and normalization using metagenomeSeq package and before decontamination
sort(rowSums(otu.bio1.unfil, na.rm = FALSE, dims = 1), decreasing = F)
sum(otu.bio1.unfil)
# load the otu table
head(otu.bio1.unfil)
otu.bio1.unfil <- rownames_to_column(otu.bio1.unfil, var = "OTUID")
# merge the taxonomy with otu table
head(tax.unfil.ed)
#tax.unfil.ed <- rownames_to_column(tax.unfil.ed, var = "OTUID")
otu.tax.unfil <- merge(otu.bio1.unfil, tax.unfil.ed, by="OTUID")
dim(otu.tax.unfil)
colnames(otu.tax.unfil)

#select only the otu table and "Order"  & "Family"
#otu.tax.unfil.ed <- otu.tax.unfil[,c(1:48,52,53)]
#colnames(otu.tax.unfil.ed)

#edit the taxonomy
colnames(otu.tax.unfil)

otu.tax.unfil.ed <- otu.tax.unfil %>%
    mutate(Taxonomy = case_when(Order == "Chloroplast" ~ 'Chloroplast',
                                Phylum  == "Cyanobacteria"~ 'Chloroplast',
                                  Family == "Mitochondria" ~ 'Mitochondria',
                                  #Family == "Magnoliophyta" ~ 'Magnoliophyta',
                                  TRUE ~ 'Bacteria')) %>%
    mutate(Domain = case_when(Order == "Chloroplast" ~ 'Plant',
                              Phylum  == "Cyanobacteria"~ 'Plant',
                                  Family == "Mitochondria" ~ 'Plant',
                                  #Family == "Magnoliophyta" ~ 'Plant',
                                  TRUE ~ 'Bacteria'))

tail(otu.tax.unfil.ed)
otu.tax.unfil.ed
colnames(otu.tax.unfil.ed)
otu.tax.unfil.ed1 <- otu.tax.unfil.ed[,c(1:66,75)]
View(otu.tax.unfil.ed1)
colnames(otu.tax.unfil.ed1)
tail(otu.tax.unfil.ed1)

long.dat <- gather(otu.tax.unfil.ed1, Sample, Read, 2:65, factor_key = T)
long.dat

### 1. Plant contaminant proportion
detach(package:plyr)
df.unfil <- long.dat %>%
  group_by(Sample, Domain) %>%
  summarise(read.number = sum(Read))
df.unfil1 <- df.unfil %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)

#with(df.unfil1, sum(percent[Sample ==  "1001"]))

library(ggbeeswarm)
library(ggtext)
plot.unfil.dom <- ggplot(df.unfil1, aes(x=Domain, y=percent, fill=Domain))+
                    geom_violin(trim = F, scale="width") +
                    scale_fill_manual(labels = c("Bacteria","Plant"),values=c("#CC79A7", "#009E73"))+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(title = "A. Experimental Sample")+
                    ylab("Read Proportion (%)")+
                    theme(legend.position="none",
                          axis.title.x = element_blank(),
                          axis.text= element_text(size = 12),
                          strip.text = element_text(size=12),
                          plot.title = element_text(size = 14),
                          axis.title.y = element_markdown(size=13),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=7, color="red", shape=95)
                        

plot.unfil.dom
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_plant_proportion.eps",
       plot.unfil.dom, device=cairo_ps,
       width = 5, height =5, 
       units= "in", dpi = 600)


### 2. Chloroplast and Mitochondria contaminant proportion
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
                    labs(title = "B")+
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
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_chloromito_proportion.eps",
       plot.unfil.tax, device=cairo_ps,
       width = 7, height =5, 
       units= "in", dpi = 600)

#############################################################################################################################################################
#READ PROPORTION OF PLANT CONTAMINANTS OF NEGATIVE CONTROLS

# otu table of the negative control
colnames(otu.unfil)
NC.unfiltered <- otu.unfil[,c(1,73:79)]#only negative control
colnames(NC.unfiltered)
NC.unfiltered <- column_to_rownames(NC.unfiltered,var="OTUID")
sort(rowSums(NC.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
NC1.unfiltered=NC.unfiltered[which(rowSums(NC.unfiltered) > 0),]
sort(rowSums(NC1.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
NC1.unfiltered <- rownames_to_column(NC1.unfiltered,var="OTUID")
NC1.tax.unfiltered <- merge(NC1.unfiltered, tax.unfil.ed, by="OTUID")
NC1.unfiltered <- column_to_rownames(NC1.unfiltered,var="OTUID")
#write.csv(NC1.tax.unfiltered, file = "NC1.tax.unfiltered.csv")
head(NC1.unfiltered)
colnames(NC1.unfiltered)

#edit the taxonomy
colnames(NC1.tax.unfiltered)

NC1.tax.unfil.ed <- NC1.tax.unfiltered %>%
    mutate(Domain = case_when(Order == "Chloroplast" ~ 'Plant',
                                  Family == "Mitochondria" ~ 'Plant',
                                  TRUE ~ 'Bacteria'))


colnames(NC1.tax.unfil.ed)
NC1.tax.unfil.ed1 <- NC1.tax.unfil.ed[,c(1:9)]
colnames(NC1.tax.unfil.ed1)
#NC1.tax.unfil.ed1 <- column_to_rownames(NC1.tax.unfil.ed1, var = "OTUID")
#nc1.tax.unfil.ed1$Taxonomy <- as.factor(nc1.tax.unfil.ed)
tail(NC1.tax.unfil.ed1)
str(NC1.tax.unfil.ed1)
library(tidyr)
long.dat.nc.unfil <- gather(NC1.tax.unfil.ed1, Sample, Read, NC1r2:NC7r2, factor_key = T)
long.dat.nc.unfil

detach(package:plyr)
df.nc.unfil <- long.dat.nc.unfil %>%
  group_by(Sample, Domain) %>%
  summarise(read.number = sum(Read))
df.nc.unfil1 <- df.nc.unfil %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)

#with(df.nc.unfil1, sum(percent[Sample ==  "NC1r2"]))

library(ggbeeswarm)
library(ggtext)
plot.nc.unfil.dom <- ggplot(df.nc.unfil1, aes(x=Domain, y=percent, fill=Domain))+
                    geom_violin(trim = F, scale="width") +
                    scale_fill_manual(labels = c("Bacteria","Plant"),values=c("#CC79A7", "#009E73"))+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    ylab("Read Proportion (%)")+
                    theme(legend.position="none",
                          axis.title.x = element_blank(),
                          axis.text= element_text(size = 14),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 14, face = 'bold'),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=10, color="red", shape=95)
                          

plot.nc.unfil.dom
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_nc_plant_proportion.eps",
       plot.nc.unfil.dom, device=cairo_ps,
       width = 5, height =5, 
       units= "in", dpi = 600)

#############################################################################################################################################################
#READ PROPORTION OF PLANT CONTAMINANTS OF THE POSITIVE CONTROLS

# otu table of the positive control
colnames(otu.unfil)
PC.unfiltered <- otu.unfil[,c(1,66:72)]#only positive control
PC.unfiltered
PC.unfiltered <- column_to_rownames(PC.unfiltered,var="OTUID")
sort(rowSums(PC.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
PC1.unfiltered <- PC.unfiltered[which(rowSums(PC.unfiltered) > 0),]
sort(rowSums(PC1.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
PC1.unfiltered <- rownames_to_column(PC1.unfiltered,var="OTUID")
PC1.tax.unfiltered <- merge(PC1.unfiltered, tax.unfil.ed, by="OTUID")
PC1.unfiltered <- column_to_rownames(PC1.unfiltered,var="OTUID")
#write.csv(NC1.tax.unfiltered, file = "NC1.tax.unfiltered.csv")
sum(PC1.unfiltered)
dim(PC1.unfiltered)

#edit the taxonomy
colnames(PC1.tax.unfiltered)

PC1.tax.unfil.ed <- PC1.tax.unfiltered %>%
    mutate(Domain = case_when(Order == "Chloroplast" ~ 'Plant',
                                  Family == "Mitochondria" ~ 'Plant',
                                  TRUE ~ 'Bacteria'))
colnames(PC1.tax.unfil.ed)
PC1.tax.unfil.ed1 <- PC1.tax.unfil.ed[,c(1:9)]
colnames(PC1.tax.unfil.ed1)

library(tidyr)
long.dat.pc.unfil <- gather(PC1.tax.unfil.ed1, Sample, Read, Mock1r2:Mock7r2, factor_key = T)
long.dat.pc.unfil

detach(package:plyr)
df.pc.unfil <- long.dat.pc.unfil %>%
  group_by(Sample, Domain) %>%
  summarise(read.number = sum(Read))
df.pc.unfil1 <- df.pc.unfil %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)

library(ggbeeswarm)
library(ggtext)
plot.pc.unfil.dom <- ggplot(df.pc.unfil1, aes(x=Domain, y=percent, fill=Domain))+
                    geom_violin(trim = F, scale="width") +
                    scale_fill_manual(labels = c("Bacteria","Plant"),values=c("#CC79A7", "#009E73"))+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    ylab("Read Proportion (%)")+
                    theme(legend.position="none",
                          axis.title.x = element_blank(),
                          axis.text= element_text(size = 14),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 14, face = 'bold'),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=10, color="red", shape=95)
                          

plot.pc.unfil.dom
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_pc_plant_proportion.eps",
       plot.pc.unfil.dom, device=cairo_ps,
       width = 5, height =5, 
       units= "in", dpi = 600)

#############################################################################################################################################################
#READ PROPORTION OF PLANT CONTAMINANTS OF THE RTSF POSITIVE CONTROL

# otu table of the RTSF Zymo
colnames(otu.unfil)
otu.unfil <- column_to_rownames(otu.unfil, var = "OTUID")
zymo.unfiltered <- otu.unfil[,"ZymoMockDNAr2", drop=F]
zymo.unfiltered
#zymo.unfiltered <- column_to_rownames(zymo.unfiltered,var="OTUID")
sort(rowSums(zymo.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
zymo.unfiltered
zymo1.unfiltered <- subset(zymo.unfiltered,rowSums(zymo.unfiltered["ZymoMockDNAr2"]) > 0)
zymo1.unfiltered
sort(rowSums(zymo1.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
zymo1.unfiltered <- rownames_to_column(zymo1.unfiltered,var="OTUID")
zymo1.tax.unfiltered <- merge(zymo1.unfiltered, tax.unfil.ed, by="OTUID")
zymo1.unfiltered <- column_to_rownames(zymo1.unfiltered,var="OTUID")
#write.csv(zymo1.tax.unfiltered, file = "zymo1.tax.unfiltered.csv")
sum(zymo1.unfiltered)
dim(zymo1.unfiltered)

#edit the taxonomy
colnames(zymo1.tax.unfiltered)

zymo1.tax.unfil.ed <- zymo1.tax.unfiltered %>%
    mutate(Domain = case_when(Order == "Chloroplast" ~ 'Plant',
                                  Family == "Mitochondria" ~ 'Plant',
                                  TRUE ~ 'Bacteria'))
colnames(zymo1.tax.unfil.ed)
zymo1.tax.unfil.ed1 <- zymo1.tax.unfil.ed[,c(1:3)]
colnames(zymo1.tax.unfil.ed1)

library(tidyr)
long.dat.zymo.unfil <- zymo1.tax.unfil.ed1
long.dat.zymo.unfil$Read <- long.dat.zymo.unfil$ZymoMockDNAr2
long.dat.zymo.unfil

detach(package:plyr)
df.zymo.unfil <- long.dat.zymo.unfil %>%
  group_by(Domain) %>%
  summarise(read.number = sum(Read))
df.zymo.unfil1 <- df.zymo.unfil %>% 
  #group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)

library(ggbeeswarm)
library(ggtext)
plot.zymo.unfil.dom <- ggplot(df.zymo.unfil1, aes(x=Domain, y=percent, fill=Domain))+
                    geom_bar(stat='identity') +
                    #scale_fill_manual(labels = c("Bacteria","Plant"), values=as.vector(stepped(n=2)))+
                    scale_fill_manual(labels = c("Bacteria","Plant"),values=c("#CC79A7", "#009E73"))+
                    #geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    ylab("Read Proportion (%)")+
                    theme(legend.position="none",
                          axis.title.x = element_blank(),
                          axis.text= element_text(size = 14),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 14, face = 'bold'),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
                          #stat_summary(fun="median",geom="point", size=10, color="red", shape=95)
                          
plot.zymo.unfil.dom
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_zymo_plant_proportion.eps",
       plot.zymo.unfil.dom, device=cairo_ps,
       width = 5, height =5, 
       units= "in", dpi = 600)

#############################################################################################################################################################
#READ PROPORTION OF PLANT CONTAMINANTS OF THE RTSF NEGATIVE CONTROL

# otu table of the RTSF NC
colnames(otu.unfil)
#otu.unfil <- column_to_rownames(otu.unfil, var = "OTUID")
RTNC.unfiltered <- otu.unfil[,"RTSFNTCr2", drop=F]
RTNC.unfiltered
sort(rowSums(RTNC.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
RTNC1.unfiltered <- subset(RTNC.unfiltered,rowSums(RTNC.unfiltered["RTSFNTCr2"]) > 0)
RTNC1.unfiltered
sort(rowSums(RTNC1.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
RTNC1.unfiltered <- rownames_to_column(RTNC1.unfiltered,var="OTUID")
RTNC1.tax.unfiltered <- merge(RTNC1.unfiltered, tax.unfil.ed, by="OTUID")
RTNC1.unfiltered <- column_to_rownames(RTNC1.unfiltered,var="OTUID")
#write.csv(RTNC1.tax.unfiltered, file = "RTNC1.tax.unfiltered.csv")
sum(RTNC1.unfiltered)
dim(RTNC1.unfiltered)

#edit the taxonomy
colnames(RTNC1.tax.unfiltered)

RTNC1.tax.unfil.ed <- RTNC1.tax.unfiltered %>%
    mutate(Domain = case_when(Order == "Chloroplast" ~ 'Plant',
                                  Family == "Mitochondria" ~ 'Plant',
                                  TRUE ~ 'Bacteria'))
colnames(RTNC1.tax.unfil.ed)
RTNC1.tax.unfil.ed1 <- RTNC1.tax.unfil.ed[,c(1:3)]
colnames(RTNC1.tax.unfil.ed1)

library(tidyr)
long.dat.rtnc.unfil <- RTNC1.tax.unfil.ed1
long.dat.rtnc.unfil$Read <- long.dat.rtnc.unfil$RTSFNTCr2
long.dat.rtnc.unfil

detach(package:plyr)
df.rtnc.unfil <- long.dat.rtnc.unfil %>%
  group_by(Domain) %>%
  summarise(read.number = sum(Read))
df.rtnc.unfil1 <- df.rtnc.unfil %>% 
  #group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)

library(ggbeeswarm)
library(ggtext)
plot.rtnc.unfil.dom <- ggplot(df.rtnc.unfil1, aes(x=Domain, y=percent, fill=Domain))+
                    geom_bar(stat='identity') +
                    #scale_fill_manual(labels = c("Bacteria","Plant"), values=as.vector(stepped(n=2)))+
                    scale_fill_manual(labels = c("Bacteria","Plant"),values=c("#CC79A7", "#009E73"))+
                    #geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    ylab("Read Proportion (%)")+
                    theme(legend.position="none",
                          axis.title.x = element_blank(),
                          axis.text= element_text(size = 14),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 14, face = 'bold'),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
                          #stat_summary(fun="median",geom="point", size=10, color="red", shape=95)
                          

plot.rtnc.unfil.dom
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_rtnc_plant_proportion.eps",
       plot.rtnc.unfil.dom, device=cairo_ps,
       width = 5, height =5, 
       units= "in", dpi = 600)


#############################################################################################################################################################
# COMPILE ALL READ PROPORTION OF PLANT CONTAMINANTS FIGURES

library(patchwork)
plot.unfil.dom
plot.nc.unfil.dom
plot.pc.unfil.dom
plot.zymo.unfil.dom
plot.rtnc.unfil.dom
PlantContProp <- (plot.unfil.dom | plot.nc.unfil.dom | plot.pc.unfil.dom) / (plot.rtnc.unfil.dom | plot.zymo.unfil.dom)
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_rPlantContProp.eps",
       PlantContProp, device=cairo_ps,
       width = 8, height =5, 
       units= "in", dpi = 600)

# ANALYSIS OF READS AFTER CHLOROPLAST AND MITOCHONDRIA REMOVAL

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
wd <- print(getwd())
otu <- read.table('OTU_table_tax_filt.txt', sep='\t', header=T, row.names = 1, check.names = FALSE)
otu
head(otu)
colnames(otu)
tax <- otu[,'taxonomy']
str(tax)
#write.csv(tax, file = "tax.fil.csv")
dim(otu)
colnames(otu)
otu <- otu[,-81]
dim(otu) # [1] 298  79, otu table still has Mock, NC, and PC in the sample
sort(rowSums(otu, na.rm = FALSE, dims = 1), decreasing = F)
otu <- rownames_to_column(otu, var = "OTUID")

#read taxonomy
tax.ed = read.csv("tax.fil.ed.csv", header=T)
head(tax.ed)
colnames(otu)
otu <- column_to_rownames(otu, var = "OTUID")
rownames(tax.ed) <- rownames(otu)
dim(tax.ed)
#read the metadata

#select only biological sample from otu table
colnames(otu)
otu.bio <- otu[,1:64] #unselect Mock, NC, and PC from the otu table
colnames(otu.bio)
dim(otu.bio)
#otu.bio <- column_to_rownames(otu.bio,var = "OTUID")
sort(rowSums(otu.bio, na.rm = FALSE, dims = 1), decreasing = F)
# remove OTUs that do not present in sample
otu.bio1=otu.bio[which(rowSums(otu.bio) > 0),]
dim(otu.bio1) # otu= 218, otu table before normalization using metagenomeSeq package and before decontamination
sort(rowSums(otu.bio1, na.rm = FALSE, dims = 1), decreasing = F)
# merge otu.bio1 with taxonomy to have match taxonomy table
head(otu.bio1)
#otu.bio1 <- rownames_to_column(otu.bio1,var = "OTUID")
head(tax.ed)
tax.ed <- rownames_to_column(tax.ed,var = "OTUID")
otu.bio1 <- rownames_to_column(otu.bio1,var = "OTUID")
otu.bio1.tax <- merge(otu.bio1, tax.ed, by="OTUID")
dim(otu.bio1.tax) 
#  separate  the sample 
# otu table
otu.bac.fil <- otu.bio1.tax[,c(1:65)]
head(otu.bac.fil)
otu.bac.fil <- column_to_rownames(otu.bac.fil,var="OTUID")
sum(otu.bac.fil)
dim(otu.bac.fil)

#otu table of the negative control
NC <- otu[,c(72:78)]#only negative control
NC
#NC <- column_to_rownames(NC,var="OTUID")
sort(rowSums(NC, na.rm = FALSE, dims = 1), decreasing = F)
NC1=NC[which(rowSums(NC) > 0),]
sort(rowSums(NC1, na.rm = FALSE, dims = 1), decreasing = F)
NC1
NC1 <- rownames_to_column(NC1,var="OTUID")
tax.ed
NC1.tax <- merge(NC1, tax.ed, by="OTUID")
#write.csv(NC1.tax, file = "NC1.tax.csv")
dim(NC1)
NC1 <- column_to_rownames(NC1,var="OTUID")
sum(NC1)

#otu table of the positive control
colnames(otu)
PC <- otu[,c(65:71)]#only positive control
PC
#PC <- column_to_rownames(PC,var="OTUID")
sort(rowSums(PC, na.rm = FALSE, dims = 1), decreasing = F)
PC1=PC[which(rowSums(PC) > 0),]
sort(rowSums(PC1, na.rm = FALSE, dims = 1), decreasing = F)
PC1
PC1 <- rownames_to_column(PC1,var="OTUID")
tax.ed
PC1.tax <- merge(PC1, tax.ed, by="OTUID")
#write.csv(PC1.tax, file = "PC1.tax.csv")
dim(PC1)
PC1 <- column_to_rownames(PC1,var="OTUID")
sum(PC1)

# otu table of the RTSF Zymo
colnames(otu)
zymo.fil <- otu[,"ZymoMockDNAr2", drop=F]
zymo.fil
zymo.fil <- column_to_rownames(zymo.fil,var="OTUID")
sort(rowSums(zymo.fil, na.rm = FALSE, dims = 1), decreasing = F)
zymo.fil
zymo1.fil <- subset(zymo.fil,rowSums(zymo.fil["ZymoMockDNAr2"]) > 0)
zymo1.fil
sort(rowSums(zymo1.fil, na.rm = FALSE, dims = 1), decreasing = F)
zymo1.fil <- rownames_to_column(zymo1.fil,var="OTUID")
zymo1.tax.fil <- merge(zymo1.fil, tax.ed, by="OTUID")
zymo1.fil <- column_to_rownames(zymo1.fil,var="OTUID")
#write.csv(zymo1.tax.fil, file = "zymo1.tax.fil.csv")
sum(zymo1.fil)
dim(zymo1.fil)

# otu table of the RTSF NC
colnames(otu)
RTNC.fil <- otu[,"RTSFNTCr2", drop=F]
RTNC.fil
sort(rowSums(RTNC.fil, na.rm = FALSE, dims = 1), decreasing = F)
RTNC1.fil <- subset(RTNC.fil,rowSums(RTNC.fil["RTSFNTCr2"]) > 0)
RTNC1.fil
sort(rowSums(RTNC1.fil, na.rm = FALSE, dims = 1), decreasing = F)
RTNC1.fil <- rownames_to_column(RTNC1.fil,var="OTUID")
RTNC1.tax.fil <- merge(RTNC1.fil, tax.ed, by="OTUID")
RTNC1.fil <- column_to_rownames(RTNC1.fil,var="OTUID")
#write.csv(RTNC1.tax.fil, file = "RTNC1.tax.fil.csv")
sum(RTNC1.fil)
dim(RTNC1.fil)
#####################################################################################################################################
######################################################################################################################################

### Rarefaction curves ######
# using GlobalPatterns

library(phyloseq)

# 1. rarefaction curve for otu table after plant contaminant removal before microbial decontamination and normalization

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
wd <- print(getwd())
otu <- read.table('OTU_table_tax_filt.txt', sep='\t', header=T, row.names = 1, check.names = FALSE)
otu
otu #otu table after plant contaminant removal
colnames(otu)
head(otu)
otu <- otu[,-81]
dim(otu) # [1] 298  79, otu table still has Mock, NC, and PC in the sample
colnames(otu)
sort(rowSums(otu, na.rm = FALSE, dims = 1), decreasing = F)
# change name of ZymoMockDNAr2 to RTSF_ZymoMockDNAr2
library(dplyr)
is.data.frame(otu)
R.utils::detachPackage("plyr")

otu <- otu %>%
  dplyr::rename(RTSF_ZymoMockDNAr2=ZymoMockDNAr2)

colnames(otu)

# make phyloseq otu table and taxonomy
otu.phyl = otu_table(otu, taxa_are_rows = TRUE)
head(tax.ed)
tax.ed <- column_to_rownames(tax.ed, var = "OTUID")
tax.phyl = tax_table(as.matrix(tax.ed))

# make phyloseq map
map <- read.csv("metadata_part.csv")
head(map)
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id
map.phyl <- sample_data(map)

# make phyloseq object

phyl.obj <- merge_phyloseq(otu.phyl,tax.phyl,map.phyl)
phyl.obj
otu_table(phyl.obj)

#set seed
set.seed(42)

#rarefy the data
# make sure to run ggrare function in the "generating_rarecurfe.r" file
# data = phyloseq object of decontaminated non normalized otu table

# run the ggrare function attached in the file "generating_rarecurve.r"
p.rare <- ggrare(phyl.obj, step = 1, color = "sample_type", label = "sample_type", se = FALSE)

#set up your own color palette
#Palette <- c("#440154FF","#1F968BFF","#FDE725FF",)
#names(Palette) <- levels(sample_data(phyl.obj)$sample_type)
#Palette

#plot the rarecurve
#p <- ggrare(psdata, step = 1000, color = "SampleType", label = "Sample", se = FALSE)
library(ggtext)
rare <- p.rare + 
 #facet_wrap(~Plant, labeller = label_both)+
 theme_bw()+
 #scale_color_manual(values = Palette)+
 scale_size_manual(values = 60)+
 ylab("Number of OTUs")+
 xlab("Number of Reads")+
 #labs(title = "(a) Bacteria/archaea")+
 theme( strip.text.x = element_text(size=14, face='bold'),
        axis.text.x=element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size =20 ,face='bold'),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x = element_text(size=15,face="bold"),
        legend.position =  "right",
        legend.title = element_text(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
    
plot(rare)
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604rarefactioncurve.pdf",
       rare, device= "pdf",
       width = 7, height = 5, 
       units= "in", dpi = 600)

#####################################################################################################################################
######################################################################################################################################

### bacterial taxa composition of all samples (after plant contaminant removal and before microbial decontamination and normalization)

# make phyloseq object
otu #otu table after plant contaminant removal
colnames(otu)
sort(rowSums(otu, na.rm = FALSE, dims = 1), decreasing = F)

# make phyloseq otu table and taxonomy
head(otu)
otu.phyl = otu_table(otu, taxa_are_rows = TRUE)
head(tax.ed)
tax.ed <- column_to_rownames(tax.ed, var = "OTUID")
tax.phyl = tax_table(as.matrix(tax.ed))

# make phyloseq map
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
map <- read.csv("metadata_part.csv")
head(map)
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id
map.phyl <- sample_data(map)

# make phyloseq object

phyl.obj <- merge_phyloseq(otu.phyl,tax.phyl,map.phyl)
phyl.obj

# merge taxa by class

# 1. class - Bacteria
bac.cl <- tax_glom(phyl.obj, taxrank = "Class", NArm = F)
bac.cl.ra <- transform_sample_counts(bac.cl, function(x) x/sum(x))
bac.cl.ra

df.cl <- psmelt(bac.cl.ra) %>%
  group_by(sample_type, Class) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.cl$Class <- as.character(df.cl$Class)
#df.cl$Class[df.cl$Mean < 0.1] <- "Other"

# barplot of bacterial/archaeal composition across pods at Phylum level
#library(rcartocolor)
#display_carto_all(colorblind_friendly = TRUE)
#my_colors = carto_pal(12, "Safe")
#my_colors

# New facet label names for plant variable
#plant.labs <- c("Plant: A", "Plant: B", "Plant: C")
#names(plant.labs) <- c("A", "B", "C")

# Create the plot

#install.packages("pals")
library(pals)

cl <- ggplot(data=df.cl, aes(x=sample_type, y=Mean, fill=Class))
plot.cl <- cl + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(values=as.vector(stepped(n=24))) +
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     #scale_fill_manual(values=c("#44AA99", "#332288", "#117733","#CC6677","#DDCC77", "#88CCEE","#661100","#AA4499" ,"#888888"))+
                     theme(legend.position="right") + 
                     guides(fill=guide_legend(nrow=5))+
                     #labs(y= "Mean Relative Abundance", x="Plant")+
                     labs(y= "Mean Relative Abundance", x="Sample Type")+
                     theme(plot.title = element_text(size = 20, face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=14),
                           #axis.line.x = element_blank(),
                           #axis.text.x = element_blank(),
                           #axis.ticks.x = element_blank(),
                           #axis.title.x = element_blank(),
                           axis.title =  element_markdown(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           #strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           guides(fill=guide_legend(ncol=1,bycol=TRUE))
                           
plot.cl


setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_barplot_class.eps",
      plot.cl, device = "eps",
       width = 9.5, height =6.5, 
       units= "in", dpi = 600)


# merge taxa by genus

# 2. genus - Bacteria
bac.gen <- tax_glom(phyl.obj, taxrank = "Genus.ed", NArm = F)
bac.gen.ra <- transform_sample_counts(bac.gen, function(x) x/sum(x))
bac.gen.ra #197 taxa

df.gen <- psmelt(bac.gen.ra) %>%
  group_by(sample_type, Genus.ed) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.gen$Genus.ed <- as.character(df.gen$Genus.ed)
df.gen$Genus.ed[df.gen$Mean < 0.005] <- "Other (less than 0.5%)"

gen <- ggplot(data=df.gen, aes(x=sample_type, y=Mean, fill=Genus.ed))
plot.gen <- gen + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(name="Genus",values=as.vector(stepped(n=24))) +
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     #scale_fill_manual(values=c("#44AA99", "#332288", "#117733","#CC6677","#DDCC77", "#88CCEE","#661100","#AA4499" ,"#888888"))+
                     theme(legend.position="right") + 
                     guides(fill=guide_legend(nrow=5))+
                     #labs(y= "Mean Relative Abundance", x="Plant")+
                     labs(y= "Mean Relative Abundance", x="Sample Type")+
                     theme(plot.title = element_text(size = 20, face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=14),
                           #axis.line.x = element_blank(),
                           #axis.text.x = element_blank(),
                           #axis.ticks.x = element_blank(),
                           #axis.title.x = element_blank(),
                           axis.title =  element_markdown(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           #strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           guides(fill=guide_legend(ncol=1,bycol=TRUE))
                       
plot.gen
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_barplot_genus.eps",
      plot.gen, device = "eps",
       width = 13, height =7.5, 
       units= "in", dpi = 600)

#####################################################################################################################################
#####################################################################################################################################

##  2. bacterial taxa found in the negative control

# make phyloseq object

# otu table of negative control only
dim(NC)
NC <- rownames_to_column(NC, var = "OTUID")
head(NC)
dim(RTNC.fil)
head(RTNC.fil)
RTNC.fil <- rownames_to_column(RTNC.fil, var = "OTUID")
colnames(RTNC.fil)
#colnames(RTNC.fil)[2] <- "RTSF_NC"
ncrtnc <- merge(NC, RTNC.fil)
head(ncrtnc)
colnames(ncrtnc)

ncrtnc <- column_to_rownames(ncrtnc, var = "OTUID")
sort(rowSums(ncrtnc, na.rm = FALSE, dims = 1), decreasing = F)
ncrtnc1 <- ncrtnc[which(rowSums(ncrtnc) > 0),]
sort(rowSums(ncrtnc1, na.rm = FALSE, dims = 1), decreasing = F)
# taxonomy negative control
head(ncrtnc1)
ncrtnc1 <- rownames_to_column(ncrtnc1, var = "OTUID")
head(tax.ed)
tax.ed <- rownames_to_column(tax.ed, var = "OTUID")
ncrtnc1.tax <- merge(ncrtnc1, tax.ed, by="OTUID")
colnames(ncrtnc1.tax)
tax.ncrtnc <- ncrtnc1.tax[,c(1,10:18)]
head(tax.ncrtnc)
# make phyloseq otu table and taxonomy
ncrtnc1 <- column_to_rownames(ncrtnc1, var = "OTUID")
ncrtnc.phyl = otu_table(ncrtnc1, taxa_are_rows = TRUE)
tax.ncrtnc <- column_to_rownames(tax.ncrtnc, var = "OTUID")
tax.ncrtnc.phyl = tax_table(as.matrix(tax.ncrtnc))

# make phyloseq map
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
map <- read.csv("metadata_part.csv")
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id
map.phyl <- sample_data(map)

# make phyloseq object
ncrtnc.phyl.obj <- merge_phyloseq(ncrtnc.phyl,tax.ncrtnc.phyl,map.phyl)
ncrtnc.phyl.obj

# 1. genus - Bacteria
ncrtnc.gen <- tax_glom(ncrtnc.phyl.obj, taxrank = "Genus.ed", NArm = F)
ncrtnc.gen.ra <- transform_sample_counts(ncrtnc.gen, function(x) x/sum(x))
ncrtnc.gen.ra #61 taxa

df.ncrtnc.gen <- psmelt(ncrtnc.gen.ra) %>%
  group_by(Sample,Genus.ed) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.ncrtnc.gen$Genus.ed <- as.character(df.ncrtnc.gen$Genus.ed)
df.ncrtnc.gen$percent.mean <- df.ncrtnc.gen$Mean*100

ncrtnc.bubble.plot <- ggplot(data=df.ncrtnc.gen, aes(x=Sample, y=Genus.ed)) + 
                     geom_point(aes(size=percent.mean), alpha = 0.75, shape = 21) + 
                     scale_size_continuous(limits = c(0.0000000000000000000001, 100), range = c(1,10), breaks = c(0.1,1,10,50)) +
                     labs(size = "Mean Relative Abundance (%)", x ="Negative Controls", y="Taxa")+
                     theme(legend.key=element_blank(),
                     axis.title =  element_markdown(size=15,face="bold"),
                     axis.text.x = element_text(colour = "black", size = 12, face = "bold", vjust = 0.95, hjust = 1, angle=45), 
                     axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
                     legend.text = element_text(size = 10, face ="bold", colour ="black"), 
                     legend.title = element_text(size = 12, face = "bold"), 
                     panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
                     legend.position = "right") +  
                     scale_fill_manual(values = colours, guide = "none")
                           
ncrtnc.bubble.plot

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_NC_RTSFNC.bubble.plot.tiff",
      ncrtnc.bubble.plot, device = "tiff",
       width = 13.8, height =7.5, 
       units= "in", dpi = 600)

##  3. bacterial taxa found in the positive control

# make phyloseq object

# otu table of positive control and RTSF_Zymo mock
dim(PC)
colnames(PC)
PC <- rownames_to_column(PC, var = "OTUID")
dim(zymo.fil)
colnames(zymo.fil)
zymo.fil <- rownames_to_column(zymo.fil, var = "OTUID")
colnames(zymo.fil)[2] <- "RTSF_ZymoMockDNAr2"
colnames(zymo.fil)
#zymo.fil <- rownames_to_column(zymo.fil, var = "OTUID")
PC.zymo <- merge(PC, zymo.fil)
PC.zymo <- column_to_rownames(PC.zymo, var = "OTUID")
sort(rowSums(PC.zymo, na.rm = FALSE, dims = 1), decreasing = F)
PC.zymo1 <- PC.zymo[which(rowSums(PC.zymo) > 0),]
sort(rowSums(PC.zymo1, na.rm = FALSE, dims = 1), decreasing = F)
colnames(PC.zymo1)
# taxonomy positive control
head(PC.zymo1)
PC.zymo1 <- rownames_to_column(PC.zymo1, var = "OTUID")
head(tax.ed)
tax.ed <- rownames_to_column(tax.ed, var = "OTUID")
PC.zymo1.tax <- merge(PC.zymo1, tax.ed, by="OTUID")
colnames(PC.zymo1.tax)
tax.PC.zymo <- PC.zymo1.tax[,c(1,10:18)]
head(tax.PC.zymo)
# make phyloseq otu table and taxonomy
PC.zymo1 <- column_to_rownames(PC.zymo1, var = "OTUID")
PC.zymo.phyl = otu_table(PC.zymo1, taxa_are_rows = TRUE)
tax.PC.zymo <- column_to_rownames(tax.PC.zymo, var = "OTUID")
tax.PC.zymo.phyl = tax_table(as.matrix(tax.PC.zymo))
# make phyloseq map
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
map <- read.csv("metadata_part.csv")
colnames(map)
head(map)
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id
map.phyl <- sample_data(map)

# make phyloseq object
PC.zymo.phyl.obj <- merge_phyloseq(PC.zymo.phyl,tax.PC.zymo.phyl,map.phyl)
PC.zymo.phyl.obj  #121 taxa

# 1. genus - Bacteria
PC.zymo.gen <- tax_glom(PC.zymo.phyl.obj, taxrank = "Genus.ed", NArm = F)
PC.zymo.gen.ra <- transform_sample_counts(PC.zymo.gen, function(x) x/sum(x))
PC.zymo.gen.ra #61 taxa

df.PC.zymo.gen <- psmelt(PC.zymo.gen.ra) %>%
  group_by(Sample,Genus.ed) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.PC.zymo.gen$Genus.ed <- as.character(df.PC.zymo.gen$Genus.ed)
df.PC.zymo.gen$percent.mean <- df.PC.zymo.gen$Mean*100

PC.zymo.bubble.plot <- ggplot(data=df.PC.zymo.gen, aes(x=Sample, y=Genus.ed)) + 
                     geom_point(aes(size=percent.mean), alpha = 0.75, shape = 21) + 
                     scale_size_continuous(limits = c(0.0000000000000000000001, 100), range = c(1,10), breaks = c(0.1,1,10,50)) +
                     labs(size = "Mean Relative Abundance (%)", y="Taxa")+
                     theme(legend.key=element_blank(),
                     axis.title.y =  element_markdown(size=15,face="bold"),
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(colour = "black", size = 12, face = "bold", vjust = 0.95,  angle=45, hjust = 1), 
                     axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
                     legend.text = element_text(size = 10, face ="bold", colour ="black"), 
                     legend.title = element_text(size = 12, face = "bold"), 
                     panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
                     legend.position = "right") +  
                     scale_fill_manual(values = colours, guide = "none")
                           
PC.zymo.bubble.plot

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_PC.zymo.bubble.plot.tiff",
      PC.zymo.bubble.plot, device = "tiff",
       width = 12.5, height =7, 
       units= "in", dpi = 600)



#####################################################################################################################################
######################################################################################################################################

### bacterial taxa composition of all samples (before plant contaminant removal and before microbial decontamination and normalization)

# make phyloseq object
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
wd <- print(getwd())
# unfiltered otu table
otu.unfil
colnames(otu.unfil)
head(otu.unfil)
colnames(otu.unfil)[80] <- "RTSF_ZymoMockDNAr2"
otu.unfil <- column_to_rownames(otu.unfil, var = "OTUID")
sort(rowSums(otu.unfil, na.rm = FALSE, dims = 1), decreasing = F)

# make phyloseq otu table and taxonomy
otu.unfil.phyl = otu_table(otu.unfil, taxa_are_rows = TRUE)
head(tax.unfil.ed)
tax.unfil.ed <- column_to_rownames(tax.unfil.ed, var = "OTUID")
tax.unfil.phyl = tax_table(as.matrix(tax.unfil.ed))

# make phyloseq map
map <- read.csv("metadata_part.csv")
head(map)
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id
map.phyl <- sample_data(map)

# make phyloseq object

phyl.unfil.obj <- merge_phyloseq(otu.unfil.phyl,tax.unfil.phyl,map.phyl)
phyl.unfil.obj
otu_table(phyl.unfil.obj)


# merge taxa by class

# 1. class - Bacteria
bac.unfil.cl <- tax_glom(phyl.unfil.obj, taxrank = "Class", NArm = F)
bac.unfil.cl.ra <- transform_sample_counts(bac.unfil.cl, function(x) x/sum(x))
bac.unfil.cl.ra #23 taxa
otu_table(bac.unfil.cl.ra)

df.unfil.cl <- psmelt(bac.unfil.cl.ra) %>%
  group_by(sample_type, Class) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.unfil.cl$Class <- as.character(df.unfil.cl$Class)

#df.cl$Class[df.cl$Mean < 0.1] <- "Other"

# Create the plot

#install.packages("pals")
library(pals)

unfil.cl <- ggplot(data=df.unfil.cl, aes(x=sample_type, y=Mean, fill=Class))
plot.unfil.cl <- unfil.cl + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(values=as.vector(stepped(n=24))) +
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     #scale_fill_manual(values=c("#44AA99", "#332288", "#117733","#CC6677","#DDCC77", "#88CCEE","#661100","#AA4499" ,"#888888"))+
                     theme(legend.position="right") + 
                     guides(fill=guide_legend(nrow=5))+
                     #labs(y= "Mean Relative Abundance", x="Plant")+
                     labs(y= "Mean Relative Abundance", x="Sample Type")+
                     theme(plot.title = element_text(size = 20, face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=14),
                           #axis.line.x = element_blank(),
                           #axis.text.x = element_blank(),
                           #axis.ticks.x = element_blank(),
                           #axis.title.x = element_blank(),
                           axis.title =  element_markdown(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           #strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           guides(fill=guide_legend(ncol=1,bycol=TRUE))
                           
plot.unfil.cl


setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_barplot_class.unfiltered.eps",
      plot.unfil.cl, device = "eps",
       width = 9.5, height =6.5, 
       units= "in", dpi = 600)


# merge taxa by genus

# 2. genus - Bacteria
bac.unfil.gen <- tax_glom(phyl.unfil.obj, taxrank = "Genus.ed", NArm = F)
bac.unfil.gen.ra <- transform_sample_counts(bac.unfil.gen, function(x) x/sum(x))
bac.unfil.gen.ra #209 taxa

df.unfil.gen <- psmelt(bac.unfil.gen.ra) %>%
  group_by(sample_type, Genus.ed) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.unfil.gen$Genus.ed <- as.character(df.unfil.gen$Genus.ed)
df.unfil.gen$Genus.ed[df.unfil.gen$Mean < 0.005] <- "Other (less than 0.5%)"

unfil.gen <- ggplot(data=df.unfil.gen, aes(x=sample_type, y=Mean, fill=Genus.ed))
plot.unfil.gen <- unfil.gen + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(name="Genus", values=as.vector(stepped(n=24))) +
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     #scale_fill_manual(values=c("#44AA99", "#332288", "#117733","#CC6677","#DDCC77", "#88CCEE","#661100","#AA4499" ,"#888888"))+
                     theme(legend.position="right") + 
                     guides(fill=guide_legend(nrow=5))+
                     #labs(y= "Mean Relative Abundance", x="Plant")+
                     labs(y= "Mean Relative Abundance", x="Sample Type")+
                     theme(plot.title = element_text(size = 20, face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=14),
                           #axis.line.x = element_blank(),
                           #axis.text.x = element_blank(),
                           #axis.ticks.x = element_blank(),
                           #axis.title.x = element_blank(),
                           axis.title =  element_markdown(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           #strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           guides(fill=guide_legend(ncol=1,bycol=TRUE))
                           
plot.unfil.gen


setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_barplot_genus.unfiltered.eps",
      plot.unfil.gen, device = "eps",
       width = 13, height =7.5, 
       units= "in", dpi = 600)

#####################################################################################################################################
#####################################################################################################################################

##  2. bacterial taxa found in the negative control before plant contamination

# make phyloseq object

# otu table of negative control only
colnames(NC.unfiltered)
head(NC.unfiltered)
sort(rowSums(NC.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
NC.unfiltered1 <- NC.unfiltered[which(rowSums(NC.unfiltered) > 0),]

# taxonomy negative control
head(NC.unfiltered1)
NC.unfiltered1 <- rownames_to_column(NC.unfiltered1, var = "OTUID")
head(tax.unfil.ed)
tax.unfil.ed <- rownames_to_column(tax.ed, var = "OTUID")
colnames(tax.unfil.ed)
NC.unfiltered1.tax <- merge(NC.unfiltered1, tax.unfil.ed, by="OTUID")
colnames(NC.unfiltered1.tax)
tax.NC.unfiltered1 <- NC.unfiltered1.tax[,c(1,10:18)]
head(tax.NC.unfiltered1)
# make phyloseq otu table and taxonomy
NC.unfiltered1 <- column_to_rownames(NC.unfiltered1, var = "OTUID")
NC.unfiltered1.phyl = otu_table(NC.unfiltered1, taxa_are_rows = TRUE)
tax.NC.unfiltered1 <- column_to_rownames(tax.NC.unfiltered1, var = "OTUID")
tax.NC.unfiltered1.phyl = tax_table(as.matrix(tax.NC.unfiltered1))

# make phyloseq map
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
map <- read.csv("metadata_part.csv")
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id
map.phyl <- sample_data(map)

# make phyloseq object
NC.unfiltered1.phyl.obj <- merge_phyloseq(NC.unfiltered1.phyl,tax.NC.unfiltered1.phyl,map.phyl)
NC.unfiltered1.phyl.obj

# 1. genus - Bacteria
NC.unfiltered1.gen <- tax_glom(NC.unfiltered1.phyl.obj, taxrank = "Genus.ed", NArm = F)
NC.unfiltered1.gen.ra <- transform_sample_counts(NC.unfiltered1.gen, function(x) x/sum(x))
NC.unfiltered1.gen.ra #52 taxa

df.NC.unfiltered1.gen <- psmelt(NC.unfiltered1.gen.ra) %>%
  group_by(Sample,Genus.ed) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.NC.unfiltered1.gen$Genus.ed <- as.character(df.NC.unfiltered1.gen$Genus.ed)
df.NC.unfiltered1.gen$percent.mean <- df.NC.unfiltered1.gen$Mean*100

NC.unfiltered1.bubble.plot <- ggplot(data=df.NC.unfiltered1.gen, aes(x=Sample, y=Genus.ed)) + 
                     geom_point(aes(size=percent.mean), alpha = 0.75, shape = 21) + 
                     scale_size_continuous(limits = c(0.0000000000000000000001, 100), range = c(1,10), breaks = c(0.1,1,10,50)) +
                     labs(size = "Mean Relative Abundance (%)", x ="Negative Controls", y="Taxa")+
                     theme(legend.key=element_blank(),
                     axis.title =  element_markdown(size=15,face="bold"),
                     axis.text.x = element_text(colour = "black", size = 12, face = "bold", vjust = 0.95, hjust = 1, angle=45), 
                     axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
                     legend.text = element_text(size = 10, face ="bold", colour ="black"), 
                     legend.title = element_text(size = 12, face = "bold"), 
                     panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
                     legend.position = "right") +  
                     scale_fill_manual(values = colours, guide = "none")
                           
NC.unfiltered1.bubble.plot

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_NC.unfiltered1.bubble.plot.tiff",
      NC.unfiltered1.bubble.plot, device = "tiff",
       width = 13.8, height =7.5, 
       units= "in", dpi = 600)

##  3. bacterial taxa found in the positive control

# make phyloseq object

# otu table of positive control and RTSF_Zymo mock
dim(PC)
colnames(PC)
PC <- rownames_to_column(PC, var = "OTUID")
dim(zymo.fil)
colnames(zymo.fil)
zymo.fil <- rownames_to_column(zymo.fil, var = "OTUID")
colnames(zymo.fil)[2] <- "RTSF_ZymoMockDNAr2"
colnames(zymo.fil)
#zymo.fil <- rownames_to_column(zymo.fil, var = "OTUID")
PC.zymo <- merge(PC, zymo.fil)
PC.zymo <- column_to_rownames(PC.zymo, var = "OTUID")
sort(rowSums(PC.zymo, na.rm = FALSE, dims = 1), decreasing = F)
PC.zymo1 <- PC.zymo[which(rowSums(PC.zymo) > 0),]
sort(rowSums(PC.zymo1, na.rm = FALSE, dims = 1), decreasing = F)
colnames(PC.zymo1)
# taxonomy positive control
head(PC.zymo1)
PC.zymo1 <- rownames_to_column(PC.zymo1, var = "OTUID")
head(tax.ed)
tax.ed <- rownames_to_column(tax.ed, var = "OTUID")
PC.zymo1.tax <- merge(PC.zymo1, tax.ed, by="OTUID")
colnames(PC.zymo1.tax)
tax.PC.zymo <- PC.zymo1.tax[,c(1,10:18)]
head(tax.PC.zymo)
# make phyloseq otu table and taxonomy
PC.zymo1 <- column_to_rownames(PC.zymo1, var = "OTUID")
PC.zymo.phyl = otu_table(PC.zymo1, taxa_are_rows = TRUE)
tax.PC.zymo <- column_to_rownames(tax.PC.zymo, var = "OTUID")
tax.PC.zymo.phyl = tax_table(as.matrix(tax.PC.zymo))
# make phyloseq map
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
map <- read.csv("metadata_part.csv")
colnames(map)
head(map)
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id
map.phyl <- sample_data(map)

# make phyloseq object
PC.zymo.phyl.obj <- merge_phyloseq(PC.zymo.phyl,tax.PC.zymo.phyl,map.phyl)
PC.zymo.phyl.obj  #121 taxa

# 1. genus - Bacteria
PC.zymo.gen <- tax_glom(PC.zymo.phyl.obj, taxrank = "Genus.ed", NArm = F)
PC.zymo.gen.ra <- transform_sample_counts(PC.zymo.gen, function(x) x/sum(x))
PC.zymo.gen.ra #61 taxa

df.PC.zymo.gen <- psmelt(PC.zymo.gen.ra) %>%
  group_by(Sample,Genus.ed) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.PC.zymo.gen$Genus.ed <- as.character(df.PC.zymo.gen$Genus.ed)
df.PC.zymo.gen$percent.mean <- df.PC.zymo.gen$Mean*100

PC.zymo.bubble.plot <- ggplot(data=df.PC.zymo.gen, aes(x=Sample, y=Genus.ed)) + 
                     geom_point(aes(size=percent.mean), alpha = 0.75, shape = 21) + 
                     scale_size_continuous(limits = c(0.0000000000000000000001, 100), range = c(1,10), breaks = c(0.1,1,10,50)) +
                     labs(size = "Mean Relative Abundance (%)", y="Taxa")+
                     theme(legend.key=element_blank(),
                     axis.title.y =  element_markdown(size=15,face="bold"),
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(colour = "black", size = 12, face = "bold", vjust = 0.95,  angle=45, hjust = 1), 
                     axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
                     legend.text = element_text(size = 10, face ="bold", colour ="black"), 
                     legend.title = element_text(size = 12, face = "bold"), 
                     panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
                     legend.position = "right") +  
                     scale_fill_manual(values = colours, guide = "none")
                           
PC.zymo.bubble.plot

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_PC.zymo.bubble.plot.tiff",
      PC.zymo.bubble.plot, device = "tiff",
       width = 12.5, height =7, 
       units= "in", dpi = 600)



#####################################################################################################################################
######################################################################################################################################

### Shared taxa among all total samples (before plant contaminants removal)

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
## 1.calculate the occupancy of each OTUID across all samples

# unfiltered otu
# unfiltered otu table
otu.unfil
colnames(otu.unfil)
head(otu.unfil)
otu.unfil <- column_to_rownames(otu.unfil, var = "OTUID")
sort(rowSums(otu.unfil, na.rm = FALSE, dims = 1), decreasing = F)

# unfiltered taxonomy
head(tax.unfil.ed)
#tax.unfil.ed <- column_to_rownames(tax.unfil.ed, var = "OTUID")
tax.unfil.ed <- rownames_to_column(tax.unfil.ed, var = "OTUID")

# read map
map <- read.csv("metadata_part.csv")
head(map)
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id

##build a long data frame joining unfiltered otu table, map, and taxonomy
longdf.unfil <- data.frame(OTUID=as.factor(rownames(otu.unfil)), otu.unfil, check.names = F) %>%
  gather(sample_id, abun, -OTUID) %>%  #keep same column nameing as in mapping file, calling counts as "abun" (abundance)
  left_join(map) %>%  #will add the info form mapping file (grouped by the 'sample_id' column)
  left_join(tax.unfil.ed)  %>% #adding the taxonomy info (grouped by the 'OTUID' column)
  group_by(OTUID, sample_id) %>%
  summarise(n=sum(abun))

##build the new table: OTUID as rownames and sample_id as colnames
widedf.unfil <- as.data.frame(spread(longdf.unfil, OTUID, n, fill=0))
rownames(widedf.unfil) <-  widedf.unfil[,1]
widedf.unfil <- widedf.unfil[,-1]
widedf.unfil <- t(widedf.unfil)

## calculate the occupancy of each OTUID across all samples
widedf.unfil.PA <- 1*((widedf.unfil>0)==1)
Occ.unfil <- rowSums(widedf.unfil.PA)/ncol(widedf.unfil.PA)
df.Occ.unfil <- as.data.frame(Occ.unfil)
df.Occ.unfil <- rownames_to_column(df.Occ.unfil, var = "OTUID")
df.Occ.unfil.tax <- merge(df.Occ.unfil, tax.unfil.ed, by="OTUID")
sort.df.Occ.unfil.tax <- df.Occ.unfil.tax[order(df.Occ.unfil.tax$Occ.unfil, decreasing = TRUE),]

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
#write.csv(sort.df.Occ.unfil.tax, file = "sort.df.Occ.unfil.tax_all.csv")

##calculate the mean relative abundance of each OTUID across all samples 
widedf.unfil.RA  <- decostand(widedf.unfil, method="total", MARGIN=2)
widedf.unfil.RA
relabund.unfil <- rowSums(widedf.unfil.RA)
df.relabund.unfil <- as.data.frame(relabund.unfil)
df.relabund.unfil$meanRelAbund <- df.relabund.unfil$relabund.unfil/ncol(widedf.unfil.RA)
df.relabund.unfil = rownames_to_column(df.relabund.unfil, var = "OTUID")
sum(df.relabund.unfil$meanRelAbund)
sort.relabund.unfil <- df.relabund.unfil[order(df.relabund.unfil$meanRelAbund, decreasing = TRUE),]

##merge OCC table and mean relative abundance table
df.Occ.ra.unfil <- merge(df.Occ.unfil, df.relabund.unfil, by.x =c("OTUID"), by.y = c("OTUID"))
df.Occ.ra.unfil.tax <- merge(df.Occ.ra.unfil, tax.unfil.ed, by="OTUID")
sort.df.Occ.ra.unfil.tax <- df.Occ.ra.unfil.tax[order(df.Occ.ra.unfil.tax$Occ.unfil, decreasing = TRUE),]
#select OTUID with occ more than and equal to 50 %
Occ50.unfil <- subset(sort.df.Occ.ra.unfil.tax , sort.df.Occ.ra.unfil.tax$Occ.unfil>= 0.5)
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
#write.csv(Occ50.unfil, file = "Occ50.unfil.csv")
Occ50.unfil.ed <- read.csv("Occ50.unfil.ed.csv")

### Occupancy-mean relative abundance across all total samples before plant contaminants removal

Occ50.unfil.plot <- ggplot(Occ50.unfil.ed,aes(x=fct_reorder(OTUID.genus, Occ.unfil, .desc=T), y=Occ.unfil))+
 geom_bar(aes(), stat="identity")+
 #coord_flip()+
 #scale_fill_manual(values = palette)+
 labs(y= "Occupancy", x="OTU.ID")+
 theme_bw()+
 coord_flip()+
 theme(plot.title = element_text(size=16, face="bold"),
       axis.text=element_text(size=12, hjust = 0.5), 
       axis.title=element_text(size=14,face="bold"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       #legend.position = "right",
       legend.position = "none",
       panel.background = element_blank(),
       panel.grid = element_blank(),
       panel.border = element_blank(),
       axis.line.x = element_line(colour = "black"),
       axis.line.y = element_line(colour = "black"),
       plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))
Occ50.unfil.plot

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_Occ50.unfil.eps",
      Occ50.unfil.plot, device = "eps",
       width = 9, height =6.5, 
       units= "in", dpi = 600)
#####################################################################################################################################
######################################################################################################################################

### Shared taxa among samples (after plant contaminants removal)

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
## 1.calculate the occupancy of each OTUID across all samples

# plant filtered otu
otu
colnames(otu)
#otu <- column_to_rownames(otu, var = "OTUID")
sort(rowSums(otu, na.rm = FALSE, dims = 1), decreasing = F)

# filtered taxonomy
head(tax.ed)
#tax.ed <- column_to_rownames(tax.ed, var = "OTUID")
tax.ed <- rownames_to_column(tax.ed, var = "OTUID")

# read map
map <- read.csv("metadata_part.csv")
head(map)
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id

##build a long data frame joining unfiltered otu table, map, and taxonomy
longdf.fil <- data.frame(OTUID=as.factor(rownames(otu)), otu, check.names = F) %>%
  gather(sample_id, abun, -OTUID) %>%  #keep same column nameing as in mapping file, calling counts as "abun" (abundance)
  left_join(map) %>%  #will add the info form mapping file (grouped by the 'sample_id' column)
  left_join(tax.ed)  %>% #adding the taxonomy info (grouped by the 'OTUID' column)
  group_by(OTUID, sample_id) %>%
  summarise(n=sum(abun))

##build the new table: OTUID as rownames and sample_id as colnames
widedf.fil <- as.data.frame(spread(longdf.fil, OTUID, n, fill=0))
rownames(widedf.fil) <-  widedf.fil[,1]
widedf.fil <- widedf.fil[,-1]
widedf.fil <- t(widedf.fil)

## calculate the occupancy of each OTUID across all samples
widedf.fil.PA <- 1*((widedf.fil>0)==1)
Occ.fil <- rowSums(widedf.fil.PA)/ncol(widedf.fil.PA)
df.Occ.fil <- as.data.frame(Occ.fil)
df.Occ.fil <- rownames_to_column(df.Occ.fil, var = "OTUID")
df.Occ.fil.tax <- merge(df.Occ.fil, tax.ed, by="OTUID")
sort.df.Occ.fil.tax <- df.Occ.fil.tax[order(df.Occ.fil.tax$Occ.fil, decreasing = TRUE),]

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
write.csv(sort.df.Occ.fil.tax, file = "sort.df.Occ.fil.tax_all.csv")
#####################################################################################################################################
######################################################################################################################################

## 2.calculate the occupancy of each OTUID across all biological samples and all negative controls before plant contaminants removal

# subset otu only biological samples and negative controls
colnames(otu.unfil)
otu.bio.nc.unfil <- data.frame(otu.unfil[,c(1:64,72:78)], check.names = F)
colnames(otu.bio.nc.unfil)

##build a long data frame joining unfiltered otu table, map, and taxonomy
longdf.bio.nc.unfil <- data.frame(OTUID=as.factor(rownames(otu.bio.nc.unfil)), otu.bio.nc.unfil, check.names = F) %>%
  gather(sample_id, abun, -OTUID) %>%  #keep same column nameing as in mapping file, calling counts as "abun" (abundance)
  left_join(map) %>%  #will add the info form mapping file (grouped by the 'sample_id' column)
  left_join(tax.unfil.ed)  %>% #adding the taxonomy info (grouped by the 'OTUID' column)
  group_by(OTUID, sample_id) %>%
  summarise(n=sum(abun))

##build the new table: OTUID as rownames and sample_id as colnames
widedf.bio.nc.unfil <- as.data.frame(spread(longdf.bio.nc.unfil, OTUID, n, fill=0))
rownames(widedf.bio.nc.unfil) <-  widedf.bio.nc.unfil[,1]
widedf.bio.nc.unfil <- widedf.bio.nc.unfil[,-1]
widedf.bio.nc.unfil <- t(widedf.bio.nc.unfil)
colnames(widedf.bio.nc.unfil)

## calculate the occupancy of each OTUID across all biological samples and all negative controls
widedf.bio.nc.unfil.PA <- 1*((widedf.bio.nc.unfil>0)==1)
Occ.bio.nc.unfil <- rowSums(widedf.bio.nc.unfil.PA)/ncol(widedf.bio.nc.unfil.PA)
df.Occ.bio.nc.unfil <- as.data.frame(Occ.bio.nc.unfil)
df.Occ.bio.nc.unfil <- rownames_to_column(df.Occ.bio.nc.unfil, var = "OTUID")
df.Occ.bio.nc.unfil.tax <- merge(df.Occ.bio.nc.unfil, tax.unfil.ed, by="OTUID")
sort.df.Occ.bio.nc.unfil.tax <- df.Occ.bio.nc.unfil.tax[order(df.Occ.bio.nc.unfil.tax$Occ.bio.nc.unfil, decreasing = TRUE),]
View(sort.df.Occ.bio.nc.unfil.tax)

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
write.csv(sort.df.Occ.bio.nc.unfil.tax, file = "sort.df.Occ.unfil.tax_BioNc.csv")

#####################################################################################################################################
######################################################################################################################################

## calculate the occupancy of each OTUID across all biological samples and all negative controls after plant contaminants removal
## what taxa are shared among experimental samples and the negative controls

# subset otu only biological samples and negative controls
colnames(otu)
otu.bio.nc.fil <- data.frame(otu[,c(1:64,72:78)], check.names = F)
colnames(otu.bio.nc.fil)

##build a long data frame joining filtered otu table, map, and taxonomy
longdf.bio.nc.fil2 <- data.frame(OTUID=as.factor(rownames(otu.bio.nc.fil)), otu.bio.nc.fil, check.names = F) %>%
  gather(sample_id, abun, -OTUID) %>%  #keep same column nameing as in mapping file, calling counts as "abun" (abundance)
  left_join(map) %>%  #will add the info form mapping file (grouped by the 'sample_id' column)
  left_join(tax.ed)  %>% #adding the taxonomy info (grouped by the 'OTUID' column)
  group_by(Genus.ed,sample_id) %>%
  summarise(n=sum(abun))

##build the new table: Genus as rownames and sample_id as colnames
widedf.bio.nc.fil2 <- as.data.frame(spread(longdf.bio.nc.fil2, Genus.ed, n, fill=0))
rownames(widedf.bio.nc.fil2) <-  widedf.bio.nc.fil2[,1]
widedf.bio.nc.fil2 <- widedf.bio.nc.fil2[,-1]
widedf.bio.nc.fil2 <- t(widedf.bio.nc.fil2)
colnames(widedf.bio.nc.fil2)

## calculate the occupancy of each Genus across all biological samples and all negative controls
widedf.bio.nc.fil.PA2 <- 1*((widedf.bio.nc.fil2>0)==1)
Occ.bio.nc.fil2 <- rowSums(widedf.bio.nc.fil.PA2)/ncol(widedf.bio.nc.fil.PA2)
df.Occ.bio.nc.fil2 <- as.data.frame(Occ.bio.nc.fil2)
df.Occ.bio.nc.fil2 <- rownames_to_column(df.Occ.bio.nc.fil2, var = "Genus")
sort.df.Occ.bio.nc.fil2 <- df.Occ.bio.nc.fil2[order(df.Occ.bio.nc.fil2$Occ.bio.nc.fil2, decreasing = TRUE),]

##calculate the mean relative abundance of each Genus across experimental samples and the negative controls
widedf.bio.nc.fil2.RA  <- decostand(widedf.bio.nc.fil2, method="total", MARGIN=2)
widedf.bio.nc.fil2.RA
relabund <- rowSums(widedf.bio.nc.fil2.RA)
df.relabund <- as.data.frame(relabund)
df.relabund$meanRelAbund <- df.relabund$relabund/ncol(widedf.bio.nc.fil2.RA)
df.relabund = rownames_to_column(df.relabund, var = "Genus")
sum(df.relabund$meanRelAbund)
sort.relabund <- df.relabund[order(df.relabund$meanRelAbund, decreasing = TRUE),]

##merge OCC table and mean relative abundance table
df.Occ.ra <- merge(df.Occ.bio.nc.fil2, df.relabund, by.x =c("Genus"), by.y = c("Genus"))
sort.df.Occ.ra <- df.Occ.ra[order(df.Occ.ra$Occ.bio.nc.fil2, decreasing = TRUE),]
#select Genus with occ more than and equal to 2 %
Occ0.02 <- subset(sort.df.Occ.ra, sort.df.Occ.ra$Occ.bio.nc.fil2 >= 0.02)
#Occ1.pf
##sort the mean relative abundance
#sort_Occ1.pf <- Occ1.pf[order(Occ1.pf$meanRelAbund, decreasing = TRUE),]

### Occupancy-mean relative abundance across calculate the occupancy of each OTUID across all biological samples and all negative controls after plant contaminants removal

Occ.bio.nc.fil.plot <- ggplot(Occ0.02,aes(x=fct_reorder(Genus, Occ.bio.nc.fil2, .desc=T), y=Occ.bio.nc.fil2))+
 geom_bar(aes(), stat="identity")+
 #coord_flip()+
 #scale_fill_manual(values = palette)+
 labs(y= "Occupancy", x="Genus")+
 theme_bw()+
 coord_flip()+
 theme(plot.title = element_text(size=16, face="bold"),
       axis.text.x=element_text(size=10,vjust = 0.5, hjust = 1), 
       axis.title=element_text(size=12,face="bold"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       #legend.position = "right",
       legend.position = "none",
       panel.background = element_blank(),
       panel.grid = element_blank(),
       panel.border = element_blank(),
       axis.line.x = element_line(colour = "black"),
       axis.line.y = element_line(colour = "black"),
       plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))
Occ.bio.nc.fil.plot

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_expe_nc_0.02.eps",
      Occ.bio.nc.fil.plot, device = "eps",
       width = 5.5, height =6, 
       units= "in", dpi = 600)

##################################################################################################################

#  Subset  OTU that present only in the negative control(not present in the biological samples)

colnames(widedf.bio.nc.unfil.PA)
unique.nc.unfil <- as.data.frame(subset(widedf.bio.nc.unfil.PA, rowSums(widedf.bio.nc.unfil.PA[,1:64]) == 0))
colnames(unique.nc.unfil)
unique.nc.unfil2 <- as.data.frame(subset(unique.nc.unfil, rowSums(unique.nc.unfil[,65:71]) > 0))
unique.nc.unfil2 <- rownames_to_column(unique.nc.unfil2, var = "OTUID")
dim(unique.nc.unfil2) #  22 OTU present only in the negative control

unique.nc.unfil.tax <- merge(unique.nc.unfil2, tax.unfil.ed, by="OTUID")
dim(unique.nc.unfil.tax)

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
write.csv(unique.nc.unfil.tax, file = "unique.nc.unfil.tax.csv")

##################################################################################################################
## Making plot for the DNA cocentration

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
dna.con = read.csv("dnaconc.csv", header=T)
library(viridis)
library(grid)
dna.con$SampleID <- as.factor(dna.con$SampleID)
dna.con$batch <- as.factor(dna.con$batch)

#create list of dna. conc. plots
dna.conc.plot <- lapply(split(dna.con,dna.con$batch), function(x){
  #relevel factor partei by wert inside this subset
  x$SampleID <- factor(x$SampleID, levels=x$SampleID[order(x$DNA_conc_ng_per_ul,decreasing=F)])

  #make the plot
  p <- ggplot(x, aes(x = SampleID, y = DNA_conc_ng_per_ul, fill = batch, width=0.75)) +
    geom_bar(stat = "identity") +
    scale_fill_discrete(drop=F)+ #to force all levels to be considered, and thus different colors
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
    theme(legend.position="none")+
    labs(y="DNA concentration (ng/ul)", x="", title=unique(x$batch))+
    coord_flip()
})

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
setEPS()
postscript("dna_conct.eps", height = 7, width = 8)
do.call(grid.arrange,(c(dna.conc.plot, ncol=3)))
dev.off()
graphics.off()

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_barplot_genus.unfiltered.eps",
      plot.unfil.gen, device = "eps",
       width = 12, height =7.5, 
       units= "in", dpi = 600)
