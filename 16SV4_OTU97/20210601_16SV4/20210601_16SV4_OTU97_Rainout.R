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

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
wd <- print(getwd())

# READ PROPORTION OF CHLOROPLAST AND MITOCHONDRIA

#read the unfiltered zotu table 
otu.unfil <- read.table(file = 'OTU_table_tax.txt', sep = '\t', header = TRUE,check.names = FALSE)
otu.unfil
tax.unfil <- otu.unfil[,'taxonomy']
tax.unfil

#write.csv(tax.unfil, file = "tax.unfil.csv")
dim(otu.unfil) #[1] 356  81
colnames(otu.unfil)
otu.unfil <- otu.unfil[,-81]
dim(otu.unfil)# otu= 356, otu table still has Mock, NC, and PC in the sample
otu.unfil <- column_to_rownames(otu.unfil,var = "OTUID")
sort(rowSums(otu.unfil, na.rm = FALSE, dims = 1), decreasing = F)

#read taxonomy
tax.unfil.ed = read.csv("tax.unfil.ed copy.csv", header=T)
rownames(tax.unfil.ed) <- rownames(otu.unfil)
dim(tax.unfil.ed) #[1] 356   7

otu.unfil <- rownames_to_column(otu.unfil,var = "OTUID")
tax.unfil.ed <- rownames_to_column(tax.unfil.ed,var = "OTUID")
otu.tax.unfiltered <- merge(otu.unfil, tax.unfil.ed, by="OTUID")
View(otu.tax.unfiltered)
colnames(otu.tax.unfiltered)
#write.csv(otu.tax.unfiltered, file = "otu.tax.unfiltered.csv")
#read the metadata

#############################################################################################################################################################

#1. select only biological sample from otu table
otu.bio.unfil <- otu.unfil[,1:65] #unselect Mock, NC, and PC from the otu table
dim(otu.bio.unfil)
colnames(otu.bio.unfil)
otu.bio.unfil <- column_to_rownames(otu.bio.unfil, var = "OTUID")
sort(rowSums(otu.bio.unfil, na.rm = FALSE, dims = 1), decreasing = F)
# remove OTUs that do not present in biological sample
otu.bio1.unfil <- otu.bio.unfil[which(rowSums(otu.bio.unfil) > 0),]
dim(otu.bio1.unfil) # [1] 280  64, otu table before plant contaminant removal  and normalization using metagenomeSeq package and before decontamination
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

plot.unfil.dom
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_plant_proportion.eps",
       plot.unfil.dom, device=cairo_ps,
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
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_chloromito_proportion.eps",
       plot.unfil.tax, device=cairo_ps,
       width = 7, height =5, 
       units= "in", dpi = 600)

# otu table of the negative control
colnames(otu.unfil)
NC.unfiltered <- otu.unfil[,c(72:78)]#only negative control
NC.unfiltered
NC.unfiltered <- column_to_rownames(NC.unfiltered,var="OTUID")
sort(rowSums(NC.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
NC1.unfiltered=NC.unfiltered[which(rowSums(NC.unfiltered) > 0),]
sort(rowSums(NC1.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
NC1.unfiltered <- rownames_to_column(NC1.unfiltered,var="OTUID")
NC1.tax.unfiltered <- merge(NC1.unfiltered, tax.unfil.ed, by="OTUID")
NC1.unfiltered <- column_to_rownames(NC1.unfiltered,var="OTUID")
#write.csv(NC1.tax.unfiltered, file = "NC1.tax.unfiltered.csv")
sum(NC1.unfiltered)
dim(NC1.unfiltered)

# otu table of the positive control
colnames(otu.unfil)
PC.unfiltered <- otu.unfil[,c(65:71)]#only positive control
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

# otu table of the RTSF Zymo
colnames(otu.unfil)
otu.unfil
zymo.unfiltered <- otu.unfil[,"ZymoMockDNA", drop=F]
zymo.unfiltered
zymo.unfiltered <- column_to_rownames(zymo.unfiltered,var="OTUID")
sort(rowSums(zymo.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
zymo.unfiltered
zymo1.unfiltered <- subset(zymo.unfiltered,rowSums(zymo.unfiltered["ZymoMockDNA"]) > 0)
zymo1.unfiltered
sort(rowSums(zymo1.unfiltered, na.rm = FALSE, dims = 1), decreasing = F)
zymo1.unfiltered <- rownames_to_column(zymo1.unfiltered,var="OTUID")
zymo1.tax.unfiltered <- merge(zymo1.unfiltered, tax.unfil.ed, by="OTUID")
zymo1.unfiltered <- column_to_rownames(zymo1.unfiltered,var="OTUID")
#write.csv(zymo1.tax.unfiltered, file = "zymo1.tax.unfiltered.csv")
sum(zymo1.unfiltered)
dim(zymo1.unfiltered)


# ANALYSIS OF READS AFTER CHLOROPLAST AND MITOCHONDRIA REMOVAL

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
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
otu <- otu[,-80]
dim(otu) # [1] 320  79, otu table still has Mock, NC, and PC in the sample
sort(rowSums(otu, na.rm = FALSE, dims = 1), decreasing = F)
otu <- rownames_to_column(otu, var = "OTUID")

#read taxonomy
tax.ed = read.csv("tax.fil.ed.csv", header=T)
head(tax.ed)

rownames(tax.ed) <- rownames(otu)
dim(tax.ed)
#read the metadata

#select only biological sample from otu table
colnames(otu)
otu.bio <- otu[,1:64] #unselect Mock, NC, and PC from the otu table
dim(otu.bio)
otu.bio
sort(rowSums(otu.bio, na.rm = FALSE, dims = 1), decreasing = F)
# remove OTUs that do not present in sample
otu.bio1=otu.bio[which(rowSums(otu.bio) > 0),]
dim(otu.bio1) # otu= 245, otu table before normalization using metagenomeSeq package and before decontamination
sort(rowSums(otu.bio1, na.rm = FALSE, dims = 1), decreasing = F)
# merge otu.bio1 with taxonomy to have match taxonomy table
head(otu.bio1)
otu.bio1 <- rownames_to_column(otu.bio1,var = "OTUID")
tax.ed <- rownames_to_column(tax.ed,var = "OTUID")
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
zymo.fil <- otu[,"ZymoMockDNA", drop=F]
zymo.fil
zymo.fil <- column_to_rownames(zymo.fil,var="OTUID")
sort(rowSums(zymo.fil, na.rm = FALSE, dims = 1), decreasing = F)
zymo.fil
zymo1.fil <- subset(zymo.fil,rowSums(zymo.fil["ZymoMockDNA"]) > 0)
zymo1.fil
sort(rowSums(zymo1.fil, na.rm = FALSE, dims = 1), decreasing = F)
zymo1.fil <- rownames_to_column(zymo1.fil,var="OTUID")
zymo1.tax.fil <- merge(zymo1.fil, tax.ed, by="OTUID")
zymo1.fil <- column_to_rownames(zymo1.fil,var="OTUID")
#write.csv(zymo1.tax.fil, file = "zymo1.tax.fil.csv")
sum(zymo1.fil)
dim(zymo1.fil)

#####################################################################################################################################
######################################################################################################################################

### Rarefaction curves ######
# using GlobalPatterns

library(phyloseq)

# 1. rarefaction curve for otu table after plant contaminant removal before microbial decontamination and normalization

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
wd <- print(getwd())
otu <- read.table('OTU_table_tax_filt.txt', sep='\t', header=T, row.names = 1, check.names = FALSE)
otu
otu #otu table after plant contaminant removal
colnames(otu)
head(otu)
otu <- otu[,-80]
dim(otu) # [1] 324  79, otu table still has Mock, NC, and PC in the sample
sort(rowSums(otu, na.rm = FALSE, dims = 1), decreasing = F)

otu <- otu %>%
  dplyr::rename(RTSF_ZymoMockDNA=ZymoMockDNA)

# make phyloseq otu table and taxonomy
otu.phyl = otu_table(otu, taxa_are_rows = TRUE)
#read taxonomy
head(tax.ed)
rownames(tax.ed) <- rownames(otu)
#tax.ed <- column_to_rownames(tax.ed, var = "OTUID")
tax.phyl = tax_table(as.matrix(tax.ed))

# make phyloseq map
map <- read.csv("metadata_part.csv")
colnames(map)
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id
map.phyl <- sample_data(map)

# make phyloseq object

phyl.obj <- merge_phyloseq(otu.phyl,tax.phyl,map.phyl)
phyl.obj

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
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601rarefactioncurve.pdf",
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
tax.ed = read.csv("tax.fil.ed.csv", header=T)
head(tax.ed)
rownames(tax.ed) = rownames(otu)
tax.ed <- column_to_rownames(tax.ed, var = "OTUID")
tax.phyl = tax_table(as.matrix(tax.ed))

# make phyloseq map
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
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


setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_barplot_class.eps",
      plot.cl, device = "eps",
       width = 9.5, height =6.5, 
       units= "in", dpi = 600)


# merge taxa by genus

# 2. genus - Bacteria
bac.gen <- tax_glom(phyl.obj, taxrank = "Genus.ed", NArm = F)
bac.gen.ra <- transform_sample_counts(bac.gen, function(x) x/sum(x))
bac.gen.ra #202 taxa

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


setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_barplot_genus.eps",
      plot.gen, device = "eps",
       width = 12, height =7.5, 
       units= "in", dpi = 600)

#####################################################################################################################################
#####################################################################################################################################

##  2. bacterial taxa found in the negative control

# make phyloseq object

# otu table of negative control only
NC1
# taxonomy NC1
colnames(NC1.tax)
tax.negative <- NC1.tax[,c(1,9:17)]
dim(tax.negative)
# make phyloseq otu table and taxonomy
NC1.phyl = otu_table(NC1, taxa_are_rows = TRUE)
tax.negative <- column_to_rownames(tax.negative, var = "OTUID")
tax.negative.phyl = tax_table(as.matrix(tax.negative))
# make phyloseq map
map <- read.csv("metadata_part.csv")
head(map)
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id
map.phyl <- sample_data(map)

# make phyloseq object
NC1.phyl.obj <- merge_phyloseq(NC1.phyl,tax.negative.phyl,map.phyl)
NC1.phyl.obj

# 1. genus - Bacteria
NC1.gen <- tax_glom(NC1.phyl.obj, taxrank = "Genus.ed", NArm = F)
NC1.gen.ra <- transform_sample_counts(NC1.gen, function(x) x/sum(x))
NC1.gen.ra #67 taxa

df.NC1.gen <- psmelt(NC1.gen.ra) %>%
  group_by(Sample,Genus.ed) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.NC1.gen$Genus.ed <- as.character(df.NC1.gen$Genus.ed)
df.NC1.gen$percent.mean <- df.NC1.gen$Mean*100

NC1.bubble.plot <- ggplot(data=df.NC1.gen, aes(x=Sample, y=Genus.ed)) + 
                     geom_point(aes(size=percent.mean), alpha = 0.75, shape = 21) + 
                     scale_size_continuous(limits = c(0.0000000000000000000001, 100), range = c(1,10), breaks = c(0.1,1,10,50)) +
                     labs(size = "Mean Relative Abundance (%)", x ="Negative Controls", y="Taxa")+
                     theme(legend.key=element_blank(),
                     axis.title =  element_markdown(size=15,face="bold"),
                     axis.text.x = element_text(colour = "black", size = 12, face = "bold", vjust = 0.5, hjust = 0.5), 
                     axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
                     legend.text = element_text(size = 10, face ="bold", colour ="black"), 
                     legend.title = element_text(size = 12, face = "bold"), 
                     panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
                     legend.position = "right") +  
                     scale_fill_manual(values = colours, guide = "none")
                           
NC1.bubble.plot

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_NC1.bubble.plot.tiff",
      NC1.bubble.plot, device = "tiff",
       width = 12.5, height =7.5, 
       units= "in", dpi = 600)

##  3. bacterial taxa found in the positive control

# make phyloseq object

# otu table of positive control and RTSF_Zymo mock
dim(PC)
colnames(PC)
PC <- rownames_to_column(PC, var = "OTUID")
dim(zymo.fil)
colnames(zymo.fil)
colnames(zymo.fil)[2] <- "RTSF_ZymoMockDNA"
colnames(zymo.fil)
zymo.fil <- rownames_to_column(zymo.fil, var = "OTUID")
PC.zymo <- merge(PC, zymo.fil)
PC.zymo <- column_to_rownames(PC.zymo, var = "OTUID")
sort(rowSums(PC.zymo, na.rm = FALSE, dims = 1), decreasing = F)
PC.zymo1 <- PC.zymo[which(rowSums(PC.zymo) > 0),]
sort(rowSums(PC.zymo1, na.rm = FALSE, dims = 1), decreasing = F)
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
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
map <- read.csv("metadata_part.csv")

head(map)
map$sample_id <- as.factor(map$sample_id)
rownames(map) <- map$sample_id
map.phyl <- sample_data(map)

# make phyloseq object
PC.zymo.phyl.obj <- merge_phyloseq(PC.zymo.phyl,tax.PC.zymo.phyl,map.phyl)
PC.zymo.phyl.obj

# 1. genus - Bacteria
PC.zymo.gen <- tax_glom(PC.zymo.phyl.obj, taxrank = "Genus.ed", NArm = F)
PC.zymo.gen.ra <- transform_sample_counts(PC.zymo.gen, function(x) x/sum(x))
PC.zymo.gen.ra #60 taxa

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

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_PC.zymo.bubble.plot.tiff",
      PC.zymo.bubble.plot, device = "tiff",
       width = 13.5, height =7, 
       units= "in", dpi = 600)
#####################################################################################################################################
######################################################################################################################################

### bacterial taxa composition of all samples (before plant contaminant removal and before microbial decontamination and normalization)

# make phyloseq object
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
wd <- print(getwd())
# unfiltered otu table
otu.unfil
colnames(otu.unfil)
otu.unfil <- otu.unfil %>%
  dplyr::rename(RTSF_ZymoMockDNA=ZymoMockDNA)
head(otu.unfil)
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

# merge taxa by class

# 1. class - Bacteria
bac.unfil.cl <- tax_glom(phyl.unfil.obj, taxrank = "Class", NArm = F)
bac.unfil.cl.ra <- transform_sample_counts(bac.unfil.cl, function(x) x/sum(x))
bac.unfil.cl.ra #22 taxa

df.unfil.cl <- psmelt(bac.unfil.cl.ra) %>%
  group_by(sample_type, Class) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.unfil.cl$Class <- as.character(df.unfil.cl$Class)

#df.cl$Class[df.cl$Mean < 0.1] <- "Other"

# Create the plot

install.packages("pals")
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


setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_barplot_class.unfiltered.eps",
      plot.unfil.cl, device = "eps",
       width = 9.5, height =6.5, 
       units= "in", dpi = 600)


# merge taxa by genus

# 2. genus - Bacteria
bac.unfil.gen <- tax_glom(phyl.unfil.obj, taxrank = "Genus.ed", NArm = F)
bac.unfil.gen.ra <- transform_sample_counts(bac.unfil.gen, function(x) x/sum(x))
bac.unfil.gen.ra #220 taxa

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


setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_barplot_genus.unfiltered.eps",
      plot.unfil.gen, device = "eps",
       width = 12, height =7.5, 
       units= "in", dpi = 600)


#####################################################################################################################################
######################################################################################################################################
### Shared taxa among ALL samples (before plant contaminants removal)

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

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
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
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
#write.csv(Occ50.unfil, file = "Occ50.unfil.csv")
Occ50.unfil.ed <- read.csv("Occ50.unfil.ed.csv")

### Occupancy-mean relative abundance across all total samples before plant contaminants removal

Occ50.unfil.plot <- ggplot(Occ50.unfil.ed,aes(x=fct_reorder(OTUID.Genus, Occ.unfil, .desc=T), y=Occ.unfil))+
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

setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_Occ50.unfil.eps",
      Occ50.unfil.plot, device = "eps",
       width = 9, height =6.5, 
       units= "in", dpi = 600)


















## 2.calculate the occupancy of each OTUID across all biological samples and all negative controls

# subset otu only biological samples and negative controls
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

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
write.csv(sort.df.Occ.bio.nc.unfil.tax, file = "sort.df.Occ.unfil.tax_BioNc.csv")

#  Subset  OTU that present only in the negative control(not present in the biological samples)

colnames(widedf.bio.nc.unfil.PA)
unique.nc.unfil <- as.data.frame(subset(widedf.bio.nc.unfil.PA, rowSums(widedf.bio.nc.unfil.PA[,1:64]) == 0))
colnames(unique.nc.unfil)
unique.nc.unfil2 <- as.data.frame(subset(unique.nc.unfil, rowSums(unique.nc.unfil[,65:71]) > 0))
unique.nc.unfil2 <- rownames_to_column(unique.nc.unfil2, var = "OTUID")
dim(unique.nc.unfil2) #  23 OTU present only in the negative control

unique.nc.unfil.tax <- merge(unique.nc.unfil2, tax.unfil.ed, by="OTUID")
dim(unique.nc.unfil.tax)

setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
write.csv(unique.nc.unfil.tax, file = "unique.nc.unfil.tax.csv")









