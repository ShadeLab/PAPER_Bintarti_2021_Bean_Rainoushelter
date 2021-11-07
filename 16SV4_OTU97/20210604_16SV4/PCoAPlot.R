############## Reads normalization using metagnomeSeq ####################

# performed on the OTU table after plant contaminants removal.
# included all controls (positive, negative, RTSF's controls).
# skipped the decontamination step using microDecon or decontam packages.

#install and load metagenomeSeq package
#BiocManager::install("metagenomeSeq", dependencies=TRUE)
library(metagenomeSeq)

#Creating a MRexperiment object


#loading count data (otu)
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
wd <- print(getwd())
otu <- read.table('OTU_table_tax_filt.txt', sep='\t', header=T, row.names = 1, check.names = FALSE)
otu#otu table after plant contaminant removal
colnames(otu)
head(otu)
otu <- otu[,-81]
dim(otu) # [1] 298  79, otu table still has Mock, NC, and PC in the sample
colnames(otu)
sort(rowSums(otu, na.rm = FALSE, dims = 1), decreasing = F)
otu.ed <- otu %>%
  dplyr::rename(RTSF_ZymoMockDNAr2=ZymoMockDNAr2)
rownames(otu.ed)
colnames(otu.ed)
dim(otu.ed)

#subset only experimental samples from the otu table
colnames(otu)
otu_exp <- otu[,1:64]
colnames(otu_exp)
dim(otu_exp)
sort(rowSums(otu_exp, na.rm = FALSE, dims = 1), decreasing = F)
otu_exp1 <- otu_exp[which(rowSums(otu_exp) > 0),]
sort(rowSums(otu_exp1, na.rm = FALSE, dims = 1), decreasing = F)
dim(otu_exp1) #[1] 218  64
otu_exp1 <- column_to_rownames(otu_exp1, var = "OTUID")
head(otu_exp1)

#loading metadata
map.ed <- read.csv("metadata_part.csv")
head(map.ed)
map.ed$sample_id <- as.factor(map.ed$sample_id)
map.ed <- column_to_rownames(map.ed, var="sample_id")
View(map.ed)
rownames(map.ed)
dim(map.ed)

# loading metadata only experimental samples
map_exp <- read.csv("metadata_part _nocontrols.csv")
head(map_exp)
map_exp$sample_id <- as.factor(map_exp$sample_id)
map_exp <- column_to_rownames(map_exp, var="sample_id")
dim(map_exp)

#loading taxonomy
tax.ed = read.csv("tax.fil.ed.csv", header=T)
head(tax.ed)
dim(tax.ed)
rownames(tax.ed) <- rownames(otu)
colnames(tax.ed)
rownames(tax.ed)
dim(tax.ed)


# loading taxonomy only experimental samples
head(otu_exp1)
otu_exp1 <- rownames_to_column(otu_exp1, var = "OTUID")
head(tax.ed)
tax.ed <- rownames_to_column(tax.ed, var = "OTUID")
otu_exp_tax <- merge(otu_exp1, tax.ed, by="OTUID")
colnames(otu_exp_tax)
head(otu_exp_tax)
tax_exp <- otu_exp_tax[,c(1,66:74)]
head(tax_exp)
tax_exp <- column_to_rownames(tax_exp, var = "OTUID")
dim(tax_exp)
otu_exp1 <- column_to_rownames(otu_exp1,var = "OTUID")
tax_exp <- tax_exp[match(rownames(otu_exp1), rownames(tax_exp)), ]
rownames(tax_exp)

#create MRexperiment object 
library(metagenomeSeq)
phenotypeData.ed <- AnnotatedDataFrame(map.ed)
phenotypeData.ed

tax.ed <- column_to_rownames(tax.ed, var = "OTUID")
OTUdata.ed <- AnnotatedDataFrame(tax.ed)
OTUdata.ed

#create MRexperiment object only for experimental samples
phenotypeData.exp <- AnnotatedDataFrame(map_exp)
phenotypeData.exp

OTUdata.exp <- AnnotatedDataFrame(tax_exp)
OTUdata.exp

#create model
model.ed <- newMRexperiment(otu.ed, phenoData = phenotypeData.ed, featureData = OTUdata.ed)
model.ed

sort(rowSums(MRcounts(model.ed), na.rm = FALSE, dims = 1), decreasing = T)

#create model only for experimental samples
model.exp <- newMRexperiment(otu_exp1, phenoData = phenotypeData.exp, featureData = OTUdata.exp)
model.exp

sort(rowSums(MRcounts(model.exp), na.rm = FALSE, dims = 1), decreasing = T)

#normalising the data

#normalise the data to account for differences due to uneven sequencing depth
#metagenomeSeq uses Cumulative Sum Scaling (CSS) normalisation instead of rarefaction
#cumNormStatFast=Calculates the percentile for which to sum counts up to and scale by.

p <- cumNormStatFast(model.ed, pFlag = T)

p.exp <- cumNormStatFast(model.exp, pFlag = T)
#cumNorm=Calculates each column’s quantile and calculates the sum up to and including that quantile.

bac.norm <- cumNorm(model.ed, p = p)
bac.norm

bac.norm.exp <- cumNorm(model.exp, p = p.exp)

#export count matrix
otu.norm <- MRcounts(bac.norm, norm = TRUE, log = F)

otu.norm <- as.data.frame(otu.norm)

head(sort(colSums(otu.norm, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums(otu.norm, na.rm = FALSE, dims = 1), decreasing = FALSE))

head(otu.norm)

dim(otu.norm) #there are 298 otus and 80 samples

#export count matrix only experimental samples
otu.norm.exp <- MRcounts(bac.norm.exp, norm = TRUE, log = F)

otu.norm.exp <- as.data.frame(otu.norm.exp)

head(sort(colSums(otu.norm.exp, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums(otu.norm.exp, na.rm = FALSE, dims = 1), decreasing = FALSE))

head(otu.norm.exp)

dim(otu.norm.exp)#there are 218 otus and 64 samples


# 1. CALCULATE BETA DIVERSITY (PCoA PLOT) 

# dissimilarity indices for community ecologist to make a distance structure (Jaccard distance between samples)
bacnorm_PA <- 1*(otu.norm>0)
bacnorm_PA
otu_dist <- vegdist(t(bacnorm_PA), binary = TRUE, method = "jaccard") #Sorensen

bacnorm_PA.exp <- 1*(otu.norm.exp>0)
bacnorm_PA.exp
otu_dist.exp <- vegdist(t(bacnorm_PA.exp), binary = TRUE, method = "jaccard") #Sorensen
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa <- cmdscale(otu_dist, eig=T)
otu_pcoa

otu_pcoa.exp <- cmdscale(otu_dist.exp, eig=T)
otu_pcoa.exp
#env <- map.ed[,c(1:6)]

# scores of PC1 and PC2
ax1.scores=otu_pcoa$points[,1]
ax2.scores=otu_pcoa$points[,2] 

# scores of PC1 and PC2 for experimental samples only
ax1.scores.exp=otu_pcoa.exp$points[,1]
ax2.scores.exp=otu_pcoa.exp$points[,2] 
#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

# calculate percent variance explained, then add to plot
ax1 <- otu_pcoa$eig[1]/sum(otu_pcoa$eig)
ax2 <- otu_pcoa$eig[2]/sum(otu_pcoa$eig)
head(map.ed)

ax1.exp <- otu_pcoa.exp$eig[1]/sum(otu_pcoa.exp$eig)
ax2.exp <- otu_pcoa.exp$eig[2]/sum(otu_pcoa.exp$eig)
head(map_exp)

#map.ed <- column_to_rownames(map.ed, var = "sample_id")
map.pcoa <- cbind(map.ed,ax1.scores,ax2.scores)

map.pcoa.exp <- cbind(map_exp, ax1.scores.exp, ax2.scores.exp)

# simple plot
pcoa_plot <- plot(ax1.scores, ax2.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))
#plot(env_fit, p.max=0.05, col="red1")
# PCoA Plot 
require("ggrepel")
library(ggrepel)
library(viridis)
cbPalette <- c("#E69F00","#009E73", "#56B4E9", "#999999", "#0072B2", "#000000")


set.seed(13)
# 1. Grouping by treatment
#define as factor in the desired order
map.pcoa$treatment <-factor(map.pcoa$treatment, levels=c("shelter","open","positive_control","negative_control","RTSF_mock","RTSF_NC"))
set.seed(13)
tret.pcoa <- ggplot(data = map.pcoa)+
            theme_bw()+
            geom_point(aes(x = ax1.scores, y = ax2.scores,col=factor(treatment)),size=5, alpha =0.9)+
            #scale_colour_manual(values=as.vector(stepped(n=24)),labels=c("with rainout sheleter","without rainout shelter","positive control", "negative control", "RTSF positive control", "RTSF negative control"))+
            scale_colour_manual(values=cbPalette,labels=c("With rainout sheleter","Without rainout shelter","Positive control", "Negative control", "RTSF positive control", "RTSF negative control")) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2,3)*100,"% var. explained", sep=""))+
            coord_fixed() +
            labs(colour = "Treatment")+
            theme(legend.position="right",
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y  =element_text(size=14), 
            axis.title.y =element_text(size=15,face="bold"),
            legend.text=element_text(size=12),
            legend.title = element_text(size = 12),
            legend.spacing.x = unit(0.05, 'cm'))
tret.pcoa

# grouping by treatment for experimental samples only
map.pcoa.exp$treatment <-factor(map.pcoa.exp$treatment, levels=c("shelter","open"))
set.seed(13)
tret.pcoa.exp <- ggplot(data = map.pcoa.exp)+
            theme_bw()+
            geom_point(aes(x = ax1.scores.exp, y = ax2.scores.exp,col=factor(treatment)),size=5, alpha =0.9)+
            #scale_colour_manual(values=as.vector(stepped(n=24)),labels=c("with rainout sheleter","without rainout shelter","positive control", "negative control", "RTSF positive control", "RTSF negative control"))+
            scale_colour_manual(values=cbPalette,labels=c("With rainout sheleter","Without rainout shelter")) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1.exp,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2.exp,3)*100,"% var. explained", sep=""))+
            coord_fixed() +
            labs(colour = "Treatment")+
            theme(legend.justification = "left",
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y  =element_text(size=14), 
            axis.title.y =element_text(size=15,face="bold"),
            legend.text=element_text(size=12),
            legend.title = element_text(size = 12),
            legend.spacing.x = unit(0.05, 'cm'))
tret.pcoa.exp


# 2. Grouping by cultivar
map.pcoa$cultivar <-factor(map.pcoa$cultivar, levels=c("R99","B18504","Cayenne", "Rosetta","positive_control","negative_control","RTSF_mock","RTSF_NC"))
cbPalette2 <- c("#CC6666","#F0E442", "#66CC99", "#9999CC", "#56B4E9", "#999999", "#0072B2", "#000000")
library(pals)
set.seed(13)
cult.pcoa <- ggplot(data = map.pcoa)+
            theme_bw()+
            geom_point(aes(x = ax1.scores, y = ax2.scores,col=factor(cultivar)),size=5, alpha =0.9)+
            #scale_colour_manual(values=as.vector(stepped(n=24)),labels=c("with rainout sheleter","without rainout shelter","positive control", "negative control", "RTSF positive control", "RTSF negative control"))+
            scale_colour_manual(values=cbPalette2,labels=c("R99","B18504","Cayenne", "Rosetta","Positive control","Negative control", "RTSF positive control", "RTSF negative control")) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2,3)*100,"% var. explained", sep=""))+
            coord_fixed() +
            labs(colour = "Cultivar")+
            theme(legend.position="right",
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            axis.text=element_text(size=14), 
            axis.title=element_text(size=15,face="bold"),
            legend.text=element_text(size=12),
            legend.title = element_text(size = 12),
            legend.spacing.x = unit(0.05, 'cm'))
cult.pcoa

# grouping by cultivar for experimental samples only
map.pcoa.exp$cultivar <-factor(map.pcoa.exp$cultivar, levels=c("R99","B18504","Cayenne", "Rosetta"))
set.seed(13)
cult.pcoa.exp <- ggplot(data = map.pcoa.exp)+
            theme_bw()+
            geom_point(aes(x = ax1.scores.exp, y = ax2.scores.exp,col=factor(cultivar)),size=5, alpha =0.9)+
            #scale_colour_manual(values=as.vector(stepped(n=24)),labels=c("with rainout sheleter","without rainout shelter","positive control", "negative control", "RTSF positive control", "RTSF negative control"))+
            scale_colour_manual(values=cbPalette2) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1.exp,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2.exp,3)*100,"% var. explained", sep=""))+
            coord_fixed() +
            labs(colour = "Cultivar")+
            theme(legend.justification = "left",
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            axis.text=element_text(size=14), 
            axis.title=element_text(size=15,face="bold"),
            legend.text=element_text(size=12),
            legend.title = element_text(size = 12),
            legend.spacing.x = unit(0.05, 'cm'))
cult.pcoa.exp

# 3. Grouping by location
cbPalette3 <- c("#CC6699", "#009999","#56B4E9", "#999999", "#0072B2", "#000000")
#define as factor in the desired order
map.pcoa$location <-factor(map.pcoa$location, levels=c("MSU","UPREC","positive_control","negative_control","RTSF_mock","RTSF_NC"))
set.seed(13)
loc.pcoa <- ggplot(data = map.pcoa)+
            theme_bw()+
            geom_point(aes(x = ax1.scores, y = ax2.scores,col=factor(location)),size=5, alpha =0.9)+
            scale_colour_manual(values=cbPalette3,labels=c("MSU","UPREC","Positive control", "Negative control", "RTSF positive control", "RTSF negative control")) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2,3)*100,"% var. explained", sep=""))+
            coord_fixed() +
            labs(colour = "Location")+
            theme(legend.position="right",
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y  =element_text(size=14), 
            axis.title.y =element_text(size=15,face="bold"), 
            legend.text=element_text(size=12),
            legend.title = element_text(size = 12),
            legend.spacing.x = unit(0.05, 'cm'))
loc.pcoa

# grouping by location for experimental samples only
map.pcoa.exp$location <-factor(map.pcoa.exp$location, levels=c("MSU","UPREC"))
set.seed(13)
loc.pcoa.exp <- ggplot(data = map.pcoa.exp)+
            theme_bw()+
            geom_point(aes(x = ax1.scores.exp, y = ax2.scores.exp, col=factor(location)),size=5, alpha =0.9)+
            #scale_colour_manual(values=as.vector(stepped(n=24)),labels=c("with rainout sheleter","without rainout shelter","positive control", "negative control", "RTSF positive control", "RTSF negative control"))+
            scale_colour_manual(values=cbPalette3)+
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1.exp,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2.exp,3)*100,"% var. explained", sep=""))+
            coord_fixed() +
            labs(colour = "Location")+
            theme(legend.justification = "left",
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y  =element_text(size=14), 
            axis.title.y =element_text(size=15,face="bold"),
            legend.text=element_text(size=12),
            legend.title = element_text(size = 12),
            legend.spacing.x = unit(0.05, 'cm'))
loc.pcoa.exp
######## Calculated the statistical analysis of beta diversity using nested permanova #########

#1. experimental samples only
set.seed(13)
otu_dist.exp #ignore location and cultivar
adonis.treat <- adonis(otu_dist.exp ~ map_exp$treatment, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
adonis.treat # p-val= 0.471

adonis.treat.loc <- adonis(otu_dist.exp ~ map_exp$treatment*map_exp$location, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
adonis.treat.loc

adonis(otu_dist.exp ~ map_exp$treatment+map_exp$location, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)


adonis.loc <- adonis(otu_dist.exp ~ map_exp$location, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
adonis.loc # p-val= 0.169

adonis.cul <- adonis(otu_dist.exp ~ map_exp$cultivar, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
adonis.cul # p-val= 0.028*

adonis.rep <- adonis(otu_dist.exp ~ map_exp$rep, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
adonis.rep # p-val= 0.004**

adonis.batch <- adonis(otu_dist.exp.nc ~ map_exp_nc$batch, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
adonis.batch # p-val= 0.001 ***


map_exp$rep <- as.factor(map_exp$rep)
map_exp$batch <- as.factor(map_exp$batch)
map_exp$rep_loc <- as.factor(map_exp$rep_loc)
map_exp$location <- as.factor(map_exp$location)
map_exp$treatment <- as.factor(map_exp$treatment)
map_exp$subplot <- as.factor(map_exp$subplot)
map_exp$plot <- as.factor(map_exp$plot)
map_exp$plot_loc <- as.factor(map_exp$plot_loc)
map_exp$cultivar <- as.factor(map_exp$cultivar)

str(map_exp)

#within
set.seed(13)
adonis(otu_dist.exp ~ rep_loc+location*cultivar*treatment,
                      strata=map_exp$rep_loc, data = map_exp,
                      permutations = 999,
                      method="jaccard")
#main 
set.seed(13)
nested.npmanova(otu_dist.exp ~ location + rep_loc,
                data = map_exp,
                method = "jac", 
                permutations = 999)



####### Compile all PCoA Plots #########
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
library(patchwork)


loc.pcoa
tret.pcoa
cult.pcoa
loc.pcoa.exp
tret.pcoa.exp
cult.pcoa.exp

PCoA.Plot <- (loc.pcoa|loc.pcoa.exp)/(tret.pcoa|tret.pcoa.exp)/(cult.pcoa|cult.pcoa.exp)
PCoA.Plot

ggsave("PCoA.Plot.png",
       PCoA.Plot, device="png",
       width = 15, height = 12, 
       units= "in", dpi = 600)

################################################################################################################################
############################### PCoA plot for the experimental and negative controls only grouped by batch ########################

#subset only experimental samples and the negative controls from the otu table
colnames(otu)
otu_exp_nc <- otu[,c(1:64,72:78)]
colnames(otu_exp_nc)
dim(otu_exp_nc)
sort(rowSums(otu_exp_nc, na.rm = FALSE, dims = 1), decreasing = F)
otu_exp_nc1 <- otu_exp_nc[which(rowSums(otu_exp_nc) > 0),]
sort(rowSums(otu_exp_nc1, na.rm = FALSE, dims = 1), decreasing = F)
dim(otu_exp_nc1) #[1] 239  71
head(otu_exp_nc1)
#otu_exp_nc1 <- column_to_rownames(otu_exp_nc1, var = "OTUID")

# loading taxonomy only experimental samples
otu_exp_nc1 <- rownames_to_column(otu_exp_nc1, var = "OTUID")
head(tax.ed)
tax.ed <- rownames_to_column(tax.ed, var = "OTUID")
otu_exp_nc_tax <- merge(otu_exp_nc1, tax.ed, by="OTUID")
colnames(otu_exp_nc_tax)
head(otu_exp_nc_tax)
tax_exp_nc <- otu_exp_nc_tax[,c(1,73:81)]
head(tax_exp_nc)
tax_exp_nc <- column_to_rownames(tax_exp_nc, var = "OTUID")
dim(tax_exp_nc)
otu_exp_nc1 <- column_to_rownames(otu_exp_nc1,var = "OTUID")
tax_exp_nc <- tax_exp_nc[match(rownames(otu_exp_nc1), rownames(tax_exp_nc)), ]
rownames(tax_exp_nc)

# metadata
colnames(map.ed)
map.ed <- rownames_to_column(map.ed, var = "sample_id")
map.ed.tib <- as_tibble(map.ed)
map.ed.tib
map_exp_nc <- map.ed.tib %>%
 slice(1:64,72:78)
map_exp_nc$sample_id <- as.factor(map_exp_nc$sample_id)
map_exp_nc <- column_to_rownames(map_exp_nc, var="sample_id")

#create MRexperiment object 
library(metagenomeSeq)
phenotypeData.exp.nc <- AnnotatedDataFrame(map_exp_nc)
phenotypeData.exp.nc

OTUdata.exp.nc <- AnnotatedDataFrame(tax_exp_nc)
OTUdata.exp.nc

#create model
model.exp.nc <- newMRexperiment(otu_exp_nc1, phenoData = phenotypeData.exp.nc, featureData = OTUdata.exp.nc)
model.exp.nc
sort(rowSums(MRcounts(model.exp.nc), na.rm = FALSE, dims = 1), decreasing = T)

#normalising the data
#normalise the data to account for differences due to uneven sequencing depth
#metagenomeSeq uses Cumulative Sum Scaling (CSS) normalisation instead of rarefaction
#cumNormStatFast=Calculates the percentile for which to sum counts up to and scale by.
p.exp.nc <- cumNormStatFast(model.exp.nc, pFlag = T)

#cumNorm=Calculates each column’s quantile and calculates the sum up to and including that quantile.
bac.norm.exp.nc <- cumNorm(model.exp.nc, p = p.exp.nc)
bac.norm.exp.nc

#export count matrix
otu.norm.exp.nc <- MRcounts(bac.norm.exp.nc, norm = TRUE, log = F)
otu.norm.exp.nc <- as.data.frame(otu.norm.exp.nc)
head(sort(colSums(otu.norm.exp.nc, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums(otu.norm.exp.nc, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(otu.norm.exp.nc)
dim(otu.norm.exp.nc)#there are 239 otus and 71 samples


# 1.CALCULATE BETA DIVERSITY (PCoA PLOT) 

# dissimilarity indices for community ecologist to make a distance structure (Jaccard distance between samples)
bacnorm_PA.exp.nc <- 1*(otu.norm.exp.nc>0)
bacnorm_PA.exp.nc
otu_dist.exp.nc <- vegdist(t(bacnorm_PA.exp.nc), binary = TRUE, method = "jaccard") #Sorensen
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa.exp.nc <- cmdscale(otu_dist.exp.nc, eig=T)
otu_pcoa.exp.nc
# scores of PC1 and PC2 
ax1.scores.exp.nc=otu_pcoa.exp.nc$points[,1]
ax2.scores.exp.nc=otu_pcoa.exp.nc$points[,2] 
# calculate percent variance explained, then add to plot
ax1.exp.nc <- otu_pcoa.exp.nc$eig[1]/sum(otu_pcoa.exp.nc$eig)
ax2.exp.nc <- otu_pcoa.exp.nc$eig[2]/sum(otu_pcoa.exp.nc$eig)
map.pcoa.exp.nc <- cbind(map_exp_nc, ax1.scores.exp.nc, ax2.scores.exp.nc)
# simple plot
pcoa_plot_exp_nc <- plot(ax1.scores.exp.nc, ax2.scores.exp.nc, xlab=paste("PCoA1: ",round(ax1.exp.nc,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2.exp.nc,3)*100,"% var. explained", sep=""))
# PCoA Plot 
require("ggrepel")
library(ggrepel)
library(viridis)
cbPalette.batch <- palette.colors(palette = "Set 1")
set.seed(13)
palette.colors(palette = "Pastel 1")
nccol <- c("1" = "#FBB4AE", "2" = "#B3CDE3", "3" = "#CCEBC5", "4" = "#DECBE4", "5"="#FED9A6", "6" ="#FFFFCC", "7"="#E5D8BD")
# Grouping by batch
#define as factor in the desired order
map.pcoa.exp.nc$batch <-factor(map.pcoa.exp.nc$batch, levels=c("1","2","3","4","5","6","7"))
str(map.pcoa.exp.nc)
set.seed(13)
# neg controls 
sample_id <- c("NC1r2","NC2r2","NC3r2","NC4r2","NC5r2","NC6r2","NC7r2")
neg_controls <- data.frame(sample_id)
map.pcoa.exp.nc <- rownames_to_column(map.pcoa.exp.nc, var = "sample_id")
neg_controls <- inner_join(neg_controls, map.pcoa.exp.nc)
neg_controls$batch <-factor(neg_controls$batch, levels=c("1","2","3","4","5","6","7"))

map.pcoa.exp.exp <- map.pcoa.exp.nc %>%
  slice(1:64)


batch.pcoa <- ggplot(data = map.pcoa.exp.exp,aes(x = ax1.scores.exp.nc, y = ax2.scores.exp.nc))+
            theme_bw()+
            geom_point(aes(colour=factor(batch)), size=5)+
            #scale_colour_manual(values=as.vector(stepped(n=24)),labels=c("with rainout sheleter","without rainout shelter","positive control", "negative control", "RTSF positive control", "RTSF negative control"))+
            scale_colour_manual(values=cbPalette.batch)+
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1.exp.nc,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2.exp.nc,3)*100,"% var. explained", sep=""))+
            #coord_fixed() +
            labs(colour = "Batch")+
            theme(legend.position="right",
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            axis.text  =element_text(size=14), 
            axis.title =element_text(size=15,face="bold"),
            legend.text=element_text(size=12),
            legend.title = element_text(size = 12))+
 geom_point(data=neg_controls, 
            aes(x = ax1.scores.exp.nc, y = ax2.scores.exp.nc, shape=factor(batch)),size=5,alpha =2)+
            
            scale_shape_manual(values=seq(0,15))+
            labs(shape = "Negative control")
batch.pcoa
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("PCoA.Plot.Batch.png",
       batch.pcoa, device="png",
       width = 12, height = 8, 
       units= "in", dpi = 600)      






