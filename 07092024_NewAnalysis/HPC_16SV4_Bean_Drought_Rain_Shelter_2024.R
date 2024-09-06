
##################################################################################################################################
########### Bean seed microbiome drought 2024 (rain-shelter field experiment): Community (Statistical) Analyses  ###############################
##################################################################################################################################

# Date: July 19th 2024
# By : A. Fina Bintarti
# CNRS, Lyon, FR

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
install.packages("ggpmisc")
install.packages("fitdistrplus")
install.packages('BiocManager')
install.packages("dplyr")
install.packages("lme4")
install.packages("nlme")
install.packages("car")
install.packages("multcomp")
install.packages("DHARMa")
library(multcomp)
library(car)
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
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
library(ggpmisc)
library(tibble)
library(fitdistrplus)
library(lme4)
library(nlme)
library(DHARMa)

# Read back the phyloseq object:
# 1. unassigned, plant, and bacterial contaminants-filtered phyloseq object (not rarefied and not normalized)
decontaminated_unrare_notnormal_physeq <- readRDS(file="/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/phyloseqobj/decontaminated_unrare_notnormal_physeq.rds")
# 2. unassigned, plant, and bacterial contaminants-filtered phyloseq object (rarefied to 1165 reads (remove sample with reads < 1000))
decontaminated_rare_physeq <- readRDS(file="/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/phyloseqobj/decontaminated_rare_physeq.rds")
# 3. unassigned, plant, and bacterial contaminants-filtered phyloseq object (normalized using CSS method)
decontaminated_normal_physeq <- readRDS(file="/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/decontam/phyloseqobj/decontaminated_normal_physeq.rds")

###############################################################################################################################################
# ALPHA DIVERSITY 
###############################################################################################################################################
#BiocManager::install("DECIPHER")
library(DECIPHER)
#devtools::install_github("RIVM-IIV-Microbiome/biomeUtils", dependencies = T)
library(biomeUtils)


# Calculate alpha diversity indices:
# 1. not rarefied and not normalized
decon.asv.df <-  as.data.frame(otu_table(decontaminated_unrare_notnormal_physeq))
decon.asv.df_pa <- 1*(decon.asv.df>0)
sum(decon.asv.df)
# Richness
decon.s <- specnumber(decon.asv.df, MARGIN = 2) # richness
decon.richness <- as.data.frame(decon.s) 
# Shannon
decon.h <- vegan::diversity(t(decon.asv.df), index = 'shannon') # Shannon index
decon.shannon <- as.data.frame(decon.h)
# Pielou's Evenness
decon.pielou <- decon.h/log(decon.s) # Pielou's evenness
decon.evenness <- as.data.frame(decon.pielou)
# Simpson
decon.d <- vegan::diversity(t(decon.asv.df), index = 'simpson') # Simpson index
decon.simpson <- as.data.frame(decon.d)
# Inverse Simpson
decon.inv.d <- vegan::diversity(t(decon.asv.df), index = 'invsimpson')
# Faith's PD
decon.pd <- calculatePD(decontaminated_unrare_notnormal_physeq, justDF = TRUE, include_root = TRUE)
decon.pd.df <- decon.pd
str(decon.pd.df)
# merge
decon.pd.df$SHA <- decon.h
decon.pd.df$EVE <- decon.pielou
decon.pd.df$SIM <- decon.d
decon.pd.df$INSIM <- decon.inv.d
view(decon.pd.df)

# 2. rarefied
rare.asv.df <-  as.data.frame(otu_table(decontaminated_rare_physeq))
rare.asv.df_pa <- 1*(rare.asv.df>0)
# Richness
rare.s <- specnumber(rare.asv.df, MARGIN = 2) # richness
rare.richness <- as.data.frame(rare.s) 
# Shannon
rare.h <- vegan::diversity(t(rare.asv.df), index = 'shannon') # Shannon index
rare.shannon <- as.data.frame(rare.h)
# Pielou's Evenness
rare.pielou <- rare.h/log(rare.s) # Pielou's evenness
rare.evenness <- as.data.frame(rare.pielou)
# Simpson
rare.d <- vegan::diversity(t(rare.asv.df), index = 'simpson') # Simpson index
rare.simpson <- as.data.frame(rare.d)
# Inverse Simpson
rare.inv.d <- vegan::diversity(t(rare.asv.df), index = 'invsimpson')
# Faith's PD
rare.pd <- calculatePD(decontaminated_rare_physeq, justDF = TRUE, include_root = TRUE)
rare.pd.df <- rare.pd %>%
  dplyr::rename(PD_rare = PD,
                SR_rare = SR)
view(rare.pd.df)
# merge
rare.pd.df$SHA_rare <- rare.h
rare.pd.df$EVE_rare <- rare.pielou
rare.pd.df$SIM_rare <- rare.d
rare.pd.df$INSIM_rare <- rare.inv.d
view(rare.pd.df)

# 2. CSS normalized
norm.asv.df <-  as.data.frame(otu_table(decontaminated_normal_physeq))
norm.asv.df_pa <- 1*(norm.asv.df>0)

# Richness
norm.s <- specnumber(norm.asv.df, MARGIN = 2) # richness
norm.richness <- as.data.frame(norm.s) 
# Shannon
norm.h <- vegan::diversity(t(norm.asv.df), index = 'shannon') # Shannon index
norm.shannon <- as.data.frame(norm.h)
# Pielou's Evenness
norm.pielou <- norm.h/log(norm.s) # Pielou's evenness
norm.evenness <- as.data.frame(norm.pielou)
# Simpson
norm.d <- vegan::diversity(t(norm.asv.df), index = 'simpson') # Simpson index
norm.simpson <- as.data.frame(norm.d)
# Inverse Simpson
norm.inv.d <- vegan::diversity(t(norm.asv.df), index = 'invsimpson')
# Faith's PD
norm.pd <- calculatePD(decontaminated_normal_physeq, justDF = TRUE, include_root = TRUE)
norm.pd.df <- norm.pd %>%
  dplyr::rename(PD_norm = PD,
                SR_norm = SR)
view(norm.pd.df)
# merge
norm.pd.df$SHA_norm <- norm.h
norm.pd.df$EVE_norm <- norm.pielou
norm.pd.df$SIM_norm <- norm.d
norm.pd.df$INSIM_norm <- norm.inv.d
view(norm.pd.df)

## Merge together:
decon.alpha <- decon.pd.df %>%
  select(c(1, 15:20))
rare.alpha <- rare.pd.df %>%
  select(c(1, 15:20))
join.alpha <- left_join(decon.alpha, rare.alpha, by = "SampleID")
join.meta.data <- left_join(norm.pd.df, join.alpha, by = "SampleID")
view(join.meta.data)
str(join.meta.data)
join.meta.data$Number_of_plants_harvested <- as.numeric(join.meta.data$Number_of_plants_harvested)
int_cols <- sapply(join.meta.data, is.integer)
join.meta.data[int_cols] <- lapply(join.meta.data[int_cols], as.factor)
chr_cols <- sapply(join.meta.data, is.character)
join.meta.data[chr_cols] <- lapply(join.meta.data[chr_cols], as.factor)
# save
setwd("/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07")
#write.csv(join.meta.data, file = "join.meta.data.csv")
# Read back
join.meta.data.ed <- read.csv("join.meta.data.ed.csv")
str(join.meta.data.ed)
join.meta.data.ed$Number_of_plants_harvested <- as.numeric(join.meta.data.ed$Number_of_plants_harvested)
join.meta.data.ed$SR_norm <- as.numeric(join.meta.data.ed$SR_norm)
join.meta.data.ed$SR_rare <- as.numeric(join.meta.data.ed$SR_rare)
join.meta.data.ed$SR <- as.numeric(join.meta.data.ed$SR)
int_cols <- sapply(join.meta.data.ed, is.integer)
join.meta.data.ed[int_cols] <- lapply(join.meta.data.ed[int_cols], as.factor)
chr_cols <- sapply(join.meta.data.ed, is.character)
join.meta.data.ed[chr_cols] <- lapply(join.meta.data.ed[chr_cols], as.factor)
#_______________________________________________________________________________________________________________________
# STATS ALPHA
#_______________________________________________________________________________________________________________________
#install.packages("datarium")
#install.packages("rstatix")
#install.packages("lmerTest")
#install.packages("emmeans")
library(datarium)
library(rstatix)
library(lmerTest)
library(emmeans)
library(bestNormalize)

# the importance of the ‘drought’ effect, the ‘variety’ effect and their interaction. 
# it is also interesting to investigate the interaction of those treatment effect with location.



#(response variable ~ location ´ cultivar ´ treatment + (1|block), data = data).
# Yield ~ loc * var * treatment + loc:block + loc:block:var


# Richness:

#box plot
rich.box <- ggboxplot(
  join.meta.data.ed, x = "var", y = "SR", 
  color = "treatment", palette = "jco", facet.by = "loc")
rich.box
#linear mixed model
rich.lmer <- lme4::lmer(SR ~ treatment*var*loc + (1|block:loc),
             data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(rich.lmer, test="F", type="III") 
#Assumption tests:
rich.resid <- resid(rich.lmer)
shapiro.test(rich.resid) # not normal
hist(rich.resid) # right skewed
plot(simulateResiduals(rich.lmer)) # not okay
leveneTest(residuals(rich.lmer) ~ join.meta.data.ed$treatment*join.meta.data.ed$var*join.meta.data.ed$loc) # not homogen
rich.Lave <- join.meta.data.ed %>% levene_test(SR ~ treatment*var*loc) # not homogen
rich.Lave
# Transform the data
BNobject.rich <- bestNormalize(join.meta.data.ed$SR)
BNobject.rich #Standardized sqrt(x + a), a=0
SR.sqrt <- sqrt_x(join.meta.data.ed$SR, a = 0)
join.meta.data.ed$SR.sqrt <- SR.sqrt$x.t
# Data transformation
join.meta.data.ed$log10.SR <- log10(join.meta.data.ed$SR)
# try LMM
rich.lmer2 <- lmerTest::lmer(SR.sqrt ~ treatment*var*loc + (1|loc:block),
                        data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(rich.lmer2, test="F", type="III") 
#nested ANOVA
rich.lm <- lm(SR.sqrt ~ treatment*var*loc + loc/var/treatment,
              data = join.meta.data.ed)
anova(rich.lm)
#Assumption tests:
rich.resid2 <- resid(rich.lmer2)
shapiro.test(rich.resid2) # normal
hist(rich.resid2) # fine
plot(simulateResiduals(rich.lmer2)) # okay
leveneTest(residuals(rich.lmer2) ~ join.meta.data.ed$treatment*join.meta.data.ed$var*join.meta.data.ed$loc) # good
rich.Lave2 <- join.meta.data.ed %>% levene_test(SR.sqrt ~ treatment*var*loc) # good
rich.Lave2
#pairwise comparison between drought treatment
pwc.rich.trt <- join.meta.data.ed %>%
  group_by(var,loc) %>%
  emmeans_test(SR.sqrt ~ treatment, p.adjust.method = "BH", model = rich.lmer2)
pwc.rich.trt
pwc.rich.var <- join.meta.data.ed %>%
  group_by(loc,treatment) %>%
  emmeans_test(SR ~ var, p.adjust.method = "BH")
view(pwc.rich.var)
# pairwise t-test
#newdata <- join.meta.data.ed[-57,]
#newdata.df <- droplevels(newdata)
#str(newdata.df)
pwc_t.rich <- join.meta.data.ed %>%
  dplyr::group_by(var,loc) %>%
  pairwise_t_test(SR.sqrt ~ treatment, p.adjust.method = "BH")
#select(-df, -statistic)# Remove details
pwc_t.rich
# 1. between cultivars:
pwc.rich.cult.em <- join.meta.data.ed %>%
  group_by(loc,treatment) %>%
  emmeans_test(SR.sqrt ~ var, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model=rich.lm)
View(pwc.rich.cult.em) # nothing significant
pwc.rich.cult.t <- join.meta.data.ed %>%
  group_by(loc,treatment) %>%
  pairwise_t_test(SR.sqrt ~ var, p.adjust.method = "BH") 
View(pwc.rich.cult.t)


# Richness Plot

#install.packages("rcartocolor")
library(rcartocolor)
library(RColorBrewer)
library(ggtext)
carto_pal(n = NULL, 'Safe')
palette.colors(n = NULL, "Polychrome 36")

colfunc1 <- colorRampPalette(c("#8DA0CB","white"))
colfunc1(10)
plot(rep(1,10),col=colfunc1(10),pch=19,cex=3)

colfunc2 <- colorRampPalette(c("#FC8D62","white"))
colfunc2(10)
plot(rep(1,10),col=colfunc2(10),pch=19,cex=3)

colfunc3 <- colorRampPalette(c("#66C2A5","white"))
colfunc3(10)
plot(rep(1,10),col=colfunc3(10),pch=19,cex=3)

colfunc4 <- colorRampPalette(c("#E78AC3","white"))
colfunc4(10)
plot(rep(1,10),col=colfunc4(10),pch=19,cex=3)

join.meta.data.ed$name <- factor(join.meta.data.ed$name, 
                            levels = c('B18504', 'B18504_drought', 'Cayenne', 
                                       'Cayenne_drought', 'R99', 'R99_drought',
                                       'Rosetta', 'Rosetta_drought'),
                            labels = c('B18504_control', 'B18504_drought', 'Cayenne_control', 
                                       'Cayenne_drought', 'R99_control', 'R99_drought',
                                       'Rosetta_control', 'Rosetta_drought'))
join.meta.data.ed$loc <- factor(join.meta.data.ed$loc, levels = c("LP", "UP"),
                  labels = c("Lower Peninsula", "Upper Peninsula"))


RS.plot <- ggplot(join.meta.data.ed , aes(x=var, y=SR)) +
  #geom_violin(position = position_dodge(width = .75),size=0.5, trim = F) +
  #geom_point(shape = 21,size=2, position = position_jitterdodge(),alpha=1)+
  geom_boxplot(aes(fill=name))+
  theme_bw() +
  scale_fill_manual(values = c("#8DA0CB","#E5E9F3","#FC8D62","#FEE5DC","#66C2A5","#DCF1EB","#E78AC3","#F9E4F1"),
                    labels=c('B18504_control', 'B18504_drought', 'Cayenne_control', 
                             'Cayenne_drought', 'R99_control', 'R99_drought', 'Rosetta_control', 'Rosetta_drought'))+
  labs(y="Richness", subtitle = "C x L *")+
  facet_wrap(~ loc)+
  scale_y_continuous(limits = c(0, 150))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.text = element_text(size=25),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                                               size = 23,
                                               lineheight = 1,
                                               padding = margin(5.5, 5.5, 5.5, 5.5),
                                               margin = margin(0, 0, 5.5, 0),
                                               linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))
  #geom_label(data = aob.stat_text.BS.shan,label=aob.stat_text.BS.shan$label, hjust=0,colour="black", size=6, fontface="bold")
RS.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
# pairwise t-test
pwc_t.rich <- join.meta.data.ed %>%
  dplyr::group_by(var,loc) %>%
  pairwise_t_test(log10.SR ~ treatment, p.adjust.method = "BH")
# add x y position
rich.xy <- pwc_t.rich  %>% 
  add_xy_position(x = "var", dodge = 0.8) # 
# plotting the pairwise comparisons among treatments 
RS.plot2 <- RS.plot + 
  stat_pvalue_manual(rich.xy, x = "var", y.position = 143,
                     label = "{p.adj}{p.adj.signif}",size=6, hide.ns = T)
                     #tip.length = 0.01, hide.ns = F)
RS.plot2
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("Richness.tiff",
       RS.plot2, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")







################################################################################################

# 1. not rarefied and not normalized

#--------Faith's PD---------

#box plot
pd.decon.box <- ggboxplot(
  join.meta.data.ed, x = "var", y = "PD", 
  color = "treatment", palette = "jco", facet.by = "loc")
pd.decon.box
#linear mixed model
pd.decon.lmer <- lmerTest::lmer(PD ~ treatment*var*loc +(1|plot), 
                            data = join.meta.data.ed) #contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
anova(pd.decon.lmer, type = 3)
car::Anova(pd.decon.lmer, test="F", type="III") 
#anova
pd.decon.aov <- join.meta.data.ed %>% anova_test(PD ~ treatment*var*loc, type = 2)
pd.decon.aov
#pairwise comparisons
pwc.pd.decon <- join.meta.data.ed %>%
  group_by(loc,var) %>%
  emmeans_test(PD ~ treatment, p.adjust.method = "BH", model=pd.decon.lmer)
#select(-df, -statistic, -p)
view(pwc.pd.decon)
# pairwise t-test
newdata.df
pwc_t.pd.decon <- newdata.df %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(PD ~ treatment, paired = T, p.adjust.method = "BH")
#select(-df, -statistic)# Remove details
pwc_t.pd.decon #R99     LP    PD    open   shelter     0.00868 **


#--------SHANNON---------

#box plot
sha.decon.box <- ggboxplot(
  join.meta.data, x = "var", y = "SHA", 
  color = "treatment", palette = "jco", facet.by = "loc")
sha.decon.box
#linear mixed model
sha.decon.lmer <- lmerTest::lmer(SHA ~ treatment*var*loc + (1|block), 
                                data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
anova(sha.decon.lmer, type = 3)
car::Anova(sha.decon.lmer, test="F", type="III") 
#anova
sha.decon.aov <- join.meta.data %>% anova_test(SHA ~ treatment*var*loc, type = 3)
sha.decon.aov # var   3  47 3.563 0.021     * 1.85e-01
              # loc   1  47 7.188 0.010     * 1.33e-01
    # treatment:var   3  47 3.033 0.038     * 1.62e-01
#pairwise comparisons
pwc.sha.decon <- join.meta.data %>%
  group_by(var, loc) %>%
  emmeans_test(SHA ~ treatment, p.adjust.method = "BH", model=sha.decon.lmer)
#select(-df, -statistic, -p)
pwc.sha.decon # Cayenne UP    treatment SHA   open   shelter 0.00759 **
# pairwise t-test
pwc_t.sha.decon <- join.meta.data %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(SHA ~ treatment, paired = F, p.adjust.method = "BH")
#select(-df, -statistic)# Remove details
pwc_t.sha.decon


test.sha.decon <- lmerTest::lmer(SHA ~ treatment*var*loc +
                                   (1|plot.ed:loc),
                                data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(test.sha.decon, test="F", type="III") 


#--------EVENNESS---------

#box plot
even.decon.box <- ggboxplot(
  join.meta.data, x = "var", y = "EVE", 
  color = "treatment", palette = "jco", facet.by = "loc")
even.decon.box
#linear mixed model
even.decon.lmer <- lmerTest::lmer(EVE ~ treatment*var*loc + (1|block), 
                                 data = join.meta.data, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
anova(even.decon.lmer, type = 3)
car::Anova(even.decon.lmer, test="F", type="III")
#anova
even.decon.aov <- join.meta.data %>% anova_test(EVE ~ treatment*var*loc, type = 2)
even.decon.aov
#pairwise comparisons
pwc.even.decon <- join.meta.data %>%
  group_by(var, loc) %>%
  emmeans_test(EVE ~ treatment, p.adjust.method = "BH") #model=even.decon.lmer)
#select(-df, -statistic, -p)
pwc.even.decon #Cayenne UP    treatment EVE   open   shelter  0.00302 **
# pairwise t-test
pwc_t.even.decon <- join.meta.data %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(EVE ~ treatment, paired = F, p.adjust.method = "BH")
#select(-df, -statistic)# Remove details
pwc_t.even.decon




#--------SIMPSON---------

#box plot
simp.decon.box <- ggboxplot(
  join.meta.data, x = "var", y = "SIM", 
  color = "treatment", palette = "jco", facet.by = "loc")
simp.decon.box
#linear mixed model
simp.decon.lmer <- lmerTest::lmer(SIM ~ treatment*var*loc + (1|block), 
                                 data = join.meta.data, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
anova(simp.decon.lmer, type = 3)
car::Anova(simp.decon.lmer, test="F", type="III")
#anova
simp.decon.aov <- join.meta.data %>% anova_test(SIM ~ treatment*var*loc, type = 3)
simp.decon.aov
#pairwise comparisons
pwc.simp.decon <- join.meta.data %>%
  group_by(var, loc) %>%
  emmeans_test(SIM ~ treatment, p.adjust.method = "BH", model=simp.decon.lmer)
#select(-df, -statistic, -p)
pwc.simp.decon #Cayenne UP    treatment SIM   open   shelter  0.0000311 ****


#--------INVERSE SIMPSON---------

#box plot
invsimp.decon.box <- ggboxplot(
  join.meta.data, x = "var", y = "INSIM", 
  color = "treatment", palette = "jco", facet.by = "loc")
invsimp.decon.box
#linear mixed model
invsimp.decon.lmer <- lmerTest::lmer(INSIM ~ treatment*var*loc + (1|block), 
                                 data = join.meta.data, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
anova(invsimp.decon.lmer, type = 3)
car::Anova(invsimp.decon.lmer, test="F", type="III")
#anova
invsimp.decon.aov <- join.meta.data %>% anova_test(INSIM ~ treatment*var*loc, type = 3)
invsimp.decon.aov
#pairwise comparisons
pwc.invsimp.decon <- join.meta.data %>%
  group_by(var, loc) %>%
  emmeans_test(INSIM ~ treatment, p.adjust.method = "BH",  model=invsimp.decon.lmer)
#select(-df, -statistic, -p)
pwc.invsimp.decon #Cayenne UP    treatment INSIM open   shelter  0.0246 *


#############################################################################################
# 2. Rarefied


#--------Faith's PD---------

#box plot
pd.rare.box <- ggboxplot(
  join.meta.data.ed, x = "var", y = "PD_rare", 
  color = "treatment", palette = "jco", facet.by = "loc")
pd.rare.box
#linear mixed model
pd.rare.lmer <- lmerTest::lmer(PD_rare ~ treatment*var*loc + (1|loc:block), 
                                data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(pd.rare.lmer, test="F", type="III") 
# Anova
pd.rare.lm <- lm(PD_rare ~ treatment*var*loc + loc/var/treatment,
              data = join.meta.data.ed)
anova(pd.rare.lm)
#Assumption tests:
pd.rare.resid <- resid(pd.rare.lmer)
shapiro.test(pd.rare.resid) # normal, p-value = 0.11
hist(pd.rare.resid) # normal distributed
plot(simulateResiduals(pd.rare.lmer)) # looks okay
leveneTest(residuals(pd.rare.lmer) ~ newdata.df$treatment*newdata.df$var*newdata.df$loc) #okay 0.56
#pairwise comparison
pwc.pd.rare <- join.meta.data.ed %>%
  group_by(var, loc) %>%
  emmeans_test(PD_rare ~ treatment, p.adjust.method = "BH", model=pd.rare.lmer)
pwc.pd.rare #not signif
# pairwise t-test
pwc_t.pd.rare <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(PD_rare ~ treatment, p.adjust.method = "BH")
pwc_t.pd.rare #significant
# 1. between cultivars:
pwc.pd.rare.cult.em <- join.meta.data.ed %>%
  group_by(loc,treatment) %>%
  emmeans_test(PD_rare ~ var, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model=pd.rare.lmer)
View(pwc.pd.rare.cult.em) # nothing significant
pwc.pd.rare.cult.t <- join.meta.data.ed %>%
  group_by(loc,treatment) %>%
  pairwise_t_test(PD_rare ~ var, p.adjust.method = "BH") 
View(pwc.pd.rare.cult.t) # nothing significant


# PD RARE Plot

PD.rare.plot <- ggplot(join.meta.data.ed , aes(x=var, y=PD_rare)) +
  geom_boxplot(aes(fill=name))+
  theme_bw() +
  scale_fill_manual(values = c("#8DA0CB","#E5E9F3","#FC8D62","#FEE5DC","#66C2A5","#DCF1EB","#E78AC3","#F9E4F1"),
                    labels=c('B18504_control', 'B18504_drought', 'Cayenne_control', 
                             'Cayenne_drought', 'R99_control', 'R99_drought', 'Rosetta_control', 'Rosetta_drought'))+
  labs(y="Faith's PD", subtitle = "W x C *")+
  facet_wrap(~ loc)+
  scale_y_continuous(limits = c(0, 13))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.text = element_text(size=25),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                                               size = 23,
                                               lineheight = 1,
                                               padding = margin(5.5, 5.5, 5.5, 5.5),
                                               margin = margin(0, 0, 5.5, 0),
                                               linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))
#geom_label(data = aob.stat_text.BS.shan,label=aob.stat_text.BS.shan$label, hjust=0,colour="black", size=6, fontface="bold")
PD.rare.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
# pairwise t-test
pwc_t.pd.rare <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(PD_rare ~ treatment, p.adjust.method = "BH")
# add x y position
pd.rare.xy <- pwc_t.pd.rare  %>% 
  add_xy_position(x = "var", dodge = 0.8) # 
# plotting the pairwise comparisons among treatments 
PD.rare.plot2 <- PD.rare.plot + 
  stat_pvalue_manual(pd.rare.xy, x = "var", y.position = 12.3,
                     label = "{p.adj}{p.adj.signif}",size=6, hide.ns = T)
#tip.length = 0.01, hide.ns = F)
PD.rare.plot2
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("PD_RARE.tiff",
       PD.rare.plot2, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")


#--------SHANNON---------

#box plot
sha.rare.box <- ggboxplot(
  join.meta.data.ed, x = "var", y = "SHA_rare", 
  color = "treatment", palette = "jco", facet.by = "loc")
sha.rare.box
#linear mixed model
sha.rare.lmer <- lmerTest::lmer(SHA_rare ~ treatment*var*loc + (1|block:loc), 
                                 data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(sha.rare.lmer, test="F", type="III") 
#Assumption tests:
sha.rare.resid <- resid(sha.rare.lmer)
shapiro.test(sha.rare.resid) # not normal, p-value = 0.001
hist(sha.rare.resid) # slightly normal distributed
plot(simulateResiduals(sha.rare.lmer)) # looks okay
leveneTest(residuals(sha.rare.lmer) ~ newdata.df$treatment*newdata.df$var*newdata.df$loc) #okay 0.152
join.meta.data.ed <- join.meta.data.ed %>% select(-SHA.rare.orderNorm)
str(join.meta.data.ed)
# transform the data
#install.packages('bestNormalize')
library(bestNormalize)
BNobject <- bestNormalize(join.meta.data.ed$SHA_rare)
BNobject
SHA.rare.orderNorm <- orderNorm(join.meta.data.ed$SHA_rare)
join.meta.data.ed$SHA.rare.orderNorm <- SHA.rare.orderNorm$x.t

# try again linear mixed model
sha.rare.lmer2 <- lmerTest::lmer(SHA.rare.orderNorm ~ treatment*var*loc + (1|block:loc), 
                                data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(sha.rare.lmer2, test="F", type="III") 
#Assumption tests:
sha.rare.resid2 <- resid(sha.rare.lmer2)
shapiro.test(sha.rare.resid2) # not normal, p-value = 0.4
hist(sha.rare.resid2) # slightly normal distributed
plot(simulateResiduals(sha.rare.lmer2)) # looks okay
leveneTest(residuals(sha.rare.lmer2) ~ newdata.df$treatment*newdata.df$var*newdata.df$loc) #okay 0.47
#pairwise comparison
pwc.sha.rare <- join.meta.data.ed %>%
  group_by(var, loc) %>%
  emmeans_test(SHA.rare.orderNorm ~ treatment, p.adjust.method = "BH", model=sha.rare.lmer2)
pwc.sha.rare #signif
# pairwise t-test
pwc_t.sha.rare <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(SHA.rare.orderNorm ~ treatment, p.adjust.method = "BH")
pwc_t.sha.rare #significant
# 1. between cultivars:
pwc.sha.rare.cult.em <- join.meta.data.ed %>%
  group_by(loc,treatment) %>%
  emmeans_test(SHA.rare.orderNorm ~ var, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model=sha.rare.lmer2)
View(pwc.sha.rare.cult.em) # there are significance
pwc.sha.rare.cult.t <- join.meta.data.ed %>%
  group_by(loc,treatment) %>%
  pairwise_t_test(SHA.rare.orderNorm ~ var, p.adjust.method = "BH") 
View(pwc.sha.rare.cult.t)

# SHANNON RARE Plot

SHA.rare.plot <- ggplot(join.meta.data.ed , aes(x=var, y=SHA_rare)) +
  geom_boxplot(aes(fill=name))+
  theme_bw() +
  scale_fill_manual(values = c("#8DA0CB","#E5E9F3","#FC8D62","#FEE5DC","#66C2A5","#DCF1EB","#E78AC3","#F9E4F1"),
                    labels=c('B18504_control', 'B18504_drought', 'Cayenne_control', 
                             'Cayenne_drought', 'R99_control', 'R99_drought', 'Rosetta_control', 'Rosetta_drought'))+
  labs(y="Shannon Index", subtitle = "L *, W x C **")+
  facet_wrap(~ loc)+
  scale_y_continuous(limits = c(0, 4))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.text = element_text(size=25),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                                               size = 23,
                                               lineheight = 1,
                                               padding = margin(5.5, 5.5, 5.5, 5.5),
                                               margin = margin(0, 0, 5.5, 0),
                                               linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))
#geom_label(data = aob.stat_text.BS.shan,label=aob.stat_text.BS.shan$label, hjust=0,colour="black", size=6, fontface="bold")
SHA.rare.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
# pairwise t-test
pwc_t.sha.rare <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(SHA.rare.orderNorm ~ treatment, p.adjust.method = "BH")
# add x y position
sha.rare.xy <- pwc_t.sha.rare  %>% 
  add_xy_position(x = "var", dodge = 0.8) # 
# plotting the pairwise comparisons among treatments 
SHA.rare.plot2 <- SHA.rare.plot + 
  stat_pvalue_manual(sha.rare.xy, x = "var", y.position = 3.8,
                     label = "{p.adj}{p.adj.signif}",size=6, hide.ns = T)
#tip.length = 0.01, hide.ns = F)
SHA.rare.plot2
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("SHA_RARE.tiff",
       SHA.rare.plot2, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")


#--------EVENNESS---------

#box plot
eve.rare.box <- ggboxplot(
  join.meta.data.ed, x = "var", y = "EVE_rare", 
  color = "treatment", palette = "jco", facet.by = "loc")
eve.rare.box
#linear mixed model
even.rare.lmer <- lmerTest::lmer(EVE_rare ~ treatment*var*loc + (1|block:loc), 
                                  data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(even.rare.lmer, test="F", type="III")
#Assumption tests:
eve.rare.resid <- resid(even.rare.lmer)
shapiro.test(eve.rare.resid) # not normal
hist(eve.rare.resid) # not normal distributed
plot(simulateResiduals(even.rare.lmer)) # not good
leveneTest(residuals(even.rare.lmer) ~ newdata.df$treatment*newdata.df$var*newdata.df$loc) #not good
# transform data
BNobject.even.rare <- bestNormalize(join.meta.data.ed$EVE_rare)
BNobject.even.rare
EVE.rare.orderNorm <- orderNorm(join.meta.data.ed$EVE_rare)
join.meta.data.ed$EVE.rare.orderNorm <- EVE.rare.orderNorm$x.t
# try again linear mixed model
eve.rare.lmer2 <- lmerTest::lmer(EVE.rare.orderNorm ~ treatment*var*loc + (1|block:loc), 
                                 data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(eve.rare.lmer2, test="F", type="III")
#Assumption tests:
eve.rare.resid2 <- resid(eve.rare.lmer2)
shapiro.test(eve.rare.resid2) # normal, p-value = 0.4
hist(eve.rare.resid2) # normal distributed
plot(simulateResiduals(eve.rare.lmer2)) # looks okay
leveneTest(residuals(eve.rare.lmer2) ~ newdata.df$treatment*newdata.df$var*newdata.df$loc) #okay 0.9
#pairwise comparison
pwc.eve.rare <- join.meta.data.ed %>%
  group_by(var, loc) %>%
  emmeans_test(EVE.rare.orderNorm ~ treatment, p.adjust.method = "BH", model=eve.rare.lmer2)
pwc.eve.rare #signif
# pairwise t-test
pwc_t.eve.rare <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(EVE.rare.orderNorm ~ treatment, p.adjust.method = "BH")
pwc_t.eve.rare #significant

# EVENNESS RARE Plot

EVE.rare.plot <- ggplot(join.meta.data.ed , aes(x=var, y=EVE_rare)) +
  geom_boxplot(aes(fill=name))+
  theme_bw() +
  scale_fill_manual(values = c("#8DA0CB","#E5E9F3","#FC8D62","#FEE5DC","#66C2A5","#DCF1EB","#E78AC3","#F9E4F1"),
                    labels=c('B18504_control', 'B18504_drought', 'Cayenne_control', 
                             'Cayenne_drought', 'R99_control', 'R99_drought', 'Rosetta_control', 'Rosetta_drought'))+
  labs(y="Pielou's Evenness")+
  facet_wrap(~ loc)+
  scale_y_continuous(limits = c(0, 0.85))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.text = element_text(size=25),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                                               size = 23,
                                               lineheight = 1,
                                               padding = margin(5.5, 5.5, 5.5, 5.5),
                                               margin = margin(0, 0, 5.5, 0),
                                               linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))
EVE.rare.plot
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("EVE_RARE.tiff",
       EVE.rare.plot, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")






simp.rare.lmer <- lmerTest::lmer(SIM_rare ~ treatment*var*loc + (1|block), 
                                  data = join.meta.data, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
anova(simp.rare.lmer, type = 3)
car::Anova(simp.rare.lmer, test="F", type="III")
simp.rare.aov <- join.meta.data %>% anova_test(SIM_rare ~ treatment*var*loc, type = 3)
simp.rare.aov


invsimp.rare.lmer <- lmerTest::lmer(INSIM_rare ~ treatment*var*loc + (1|block), 
                                     data = join.meta.data, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
anova(invsimp.rare.lmer, type = 3)
car::Anova(invsimp.rare.lmer, test="F", type="III")
invsimp.rare.aov <- join.meta.data %>% anova_test(INSIM_rare ~ treatment*var*loc, type = 3)
invsimp.rare.aov


#############################################################################################
# 3. Normalized


#--------Faith's PD---------

# 3. 
#box plot
pd.norm.box <- ggboxplot(
  join.meta.data.ed, x = "var", y = "PD_norm", 
  color = "treatment", palette = "jco", facet.by = "loc")
pd.norm.box
#linear mixed model
pd.norm.lmer <- lmerTest::lmer(PD_norm ~ treatment*var*loc + (1|block:loc),
                               data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(pd.norm.lmer, test="F", type="III") 
#Assumption tests:
pd.norm.resid <- resid(pd.norm.lmer)
shapiro.test(pd.norm.resid) # slightly normal, p-value = 0.05165
hist(pd.norm.resid) # normal distributed
plot(simulateResiduals(pd.norm.lmer)) # looks okay
leveneTest(residuals(pd.norm.lmer) ~ join.meta.data.ed$treatment*join.meta.data.ed$var*join.meta.data.ed$loc) #okay 0.0677
pd.norm.Lave <- join.meta.data.ed %>% levene_test(PD_norm ~ treatment*var*loc) # fine
pd.norm.Lave
#transform data
BNobject.pd.norm <- bestNormalize(join.meta.data.ed$PD_norm)
BNobject.pd.norm #Standardized Log_b(x + a) Transformation, a=0, b=10
PD.norm.logx <- log_x(join.meta.data.ed$PD_norm, a = 0, b = 10)
join.meta.data.ed$PD.norm.logx <- PD.norm.logx$x.t
# log transform
join.meta.data.ed$log10.PD_norm <- log10(join.meta.data.ed$PD_norm)
str(join.meta.data.ed)
#try again linear mixed model
pd.norm.lmer2 <- lmerTest::lmer(log10.PD_norm ~ treatment*var*loc + (1|block:loc) ,
                               data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(pd.norm.lmer2, test="F", type="III") 
pd.norm.lmer3 <- lmerTest::lmer(PD.norm.logx ~ treatment*var*loc + (1|block:loc) ,
                                data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(pd.norm.lmer3, test="F", type="III") 
# Anova
pd.norm.lm <- lm(PD.norm.logx ~ treatment*var*loc + loc/var/treatment,
                 data = join.meta.data.ed)
anova(pd.norm.lm)
#Assumption tests:
pd.norm.resid2 <- resid(pd.norm.lmer2)
shapiro.test(pd.norm.resid2) # normal, p-value = 0.32
hist(pd.norm.resid2) # normal distributed
plot(simulateResiduals(pd.norm.lmer2)) # looks okay
leveneTest(residuals(pd.norm.lmer2) ~ join.meta.data.ed$treatment*join.meta.data.ed$var*join.meta.data.ed$loc) #okay 0.0677
pd.norm.Lave2 <- join.meta.data.ed %>%levene_test(log10.PD_norm ~ treatment*var*loc) # fine
pd.norm.Lave2
#pairwise comparison
pwc.pd.norm <- join.meta.data.ed %>%
  group_by(var, loc) %>%
  emmeans_test(PD.norm.logx ~ treatment, p.adjust.method = "BH", model=pd.norm.lmer3)
pwc.pd.norm #not signif
# pairwise t-test
pwc_t.pd.norm <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(log10.PD_norm ~ treatment, p.adjust.method = "BH")
pwc_t.pd.norm #significant R99  LP  PD_norm open  shelter 0.00743 **

# PD Plot

PD.norm.plot <- ggplot(join.meta.data.ed , aes(x=var, y=PD_norm)) +
  geom_boxplot(aes(fill=name))+
  theme_bw() +
  scale_fill_manual(values = c("#8DA0CB","#E5E9F3","#FC8D62","#FEE5DC","#66C2A5","#DCF1EB","#E78AC3","#F9E4F1"),
                    labels=c('B18504_control', 'B18504_drought', 'Cayenne_control', 
                             'Cayenne_drought', 'R99_control', 'R99_drought', 'Rosetta_control', 'Rosetta_drought'))+
  labs(y="Faith's PD")+
  facet_wrap(~ loc)+
  scale_y_continuous(limits = c(0, 13))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.text = element_text(size=25),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                                               size = 23,
                                               lineheight = 1,
                                               padding = margin(5.5, 5.5, 5.5, 5.5),
                                               margin = margin(0, 0, 5.5, 0),
                                               linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))
#geom_label(data = aob.stat_text.BS.shan,label=aob.stat_text.BS.shan$label, hjust=0,colour="black", size=6, fontface="bold")
PD.norm.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
# pairwise t-test
pwc_t.pd.norm <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(log10.PD_norm ~ treatment, p.adjust.method = "BH")
# add x y position
pd.xy <- pwc_t.pd.norm  %>% 
  add_xy_position(x = "var", dodge = 0.8) # 
# plotting the pairwise comparisons among treatments 
PD.norm.plot2 <- PD.norm.plot + 
  stat_pvalue_manual(pd.xy, x = "var", y.position = 12.3,
                     label = "{p.adj}{p.adj.signif}",size=6, hide.ns = T)
#tip.length = 0.01, hide.ns = F)
PD.norm.plot2
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("PD.tiff",
       PD.norm.plot2, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")



#--------SHANNON---------

#box plot
sha.norm.box <- ggboxplot(
  join.meta.data.ed, x = "var", y = "SHA_norm", 
  color = "treatment", palette = "jco", facet.by = "loc")
sha.norm.box
#linear mixed model
sha.norm.lmer <- lmerTest::lmer(SHA_norm ~ treatment*var*loc + (1|block:loc) , 
                                data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(sha.norm.lmer, test="F", type="III") #
#Assumption tests:
sha.norm.resid <- resid(sha.norm.lmer)
shapiro.test(sha.norm.resid) # not normal, p-value = 0.02
hist(sha.norm.resid) # slightly not normal distributed - negative skew (left)
plot(simulateResiduals(sha.norm.lmer)) # violated
leveneTest(residuals(sha.norm.lmer) ~ join.meta.data.ed$treatment*join.meta.data.ed$var*join.meta.data.ed$loc) #okay 0.0677
sha.norm.Lave <- join.meta.data.ed %>%
  #group_by(loc) %>%
  levene_test(SHA_norm ~ treatment*var*loc) # fine
sha.norm.Lave
#transform data
BNobject.sha.norm <- bestNormalize(join.meta.data.ed$SHA_norm)
BNobject.sha.norm #orderNorm Transformation
SHA.norm.orderNorm <- orderNorm(join.meta.data.ed$SHA_norm)
join.meta.data.ed$SHA.norm.orderNorm <- SHA.norm.orderNorm$x.t
# try lmm again
sha.norm.lmer2 <- lmerTest::lmer(SHA.norm.orderNorm ~ treatment*var*loc + (1|block:loc) , 
                                data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(sha.norm.lmer2, test="F", type="III") #
# try Assumption tests again
sha.norm.resid2 <- resid(sha.norm.lmer2)
shapiro.test(sha.norm.resid2) # normal, p-value = 0.19
hist(sha.norm.resid2) # normal
plot(simulateResiduals(sha.norm.lmer2)) # 
leveneTest(residuals(sha.norm.lmer2) ~ join.meta.data.ed$treatment*join.meta.data.ed$var*join.meta.data.ed$loc) #okay 0.7
sha.norm.Lave <- join.meta.data.ed %>%
  levene_test(SHA.norm.orderNorm ~ treatment*var*loc) # fine
sha.norm.Lave
#pairwise comparison
pwc.sha.norm <- join.meta.data.ed %>%
  group_by(var, loc) %>%
  emmeans_test(SHA.norm.orderNorm ~ treatment, p.adjust.method = "BH", model=sha.norm.lmer2)
pwc.sha.norm #Cayenne UP treatment box.SHA_norm open  shelter 0.0165 * 
# pairwise t-test
#newdata <- join.meta.data.ed[-57,]
#newdata.df <- droplevels(newdata)
#view(newdata.df)
pwc_t.sha.norm <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(SHA.norm.orderNorm ~ treatment, p.adjust.method = "BH")
#select(-df, -statistic)# Remove details
pwc_t.sha.norm #

# SHA Plot

SHA.norm.plot <- ggplot(join.meta.data.ed , aes(x=var, y=SHA_norm)) +
  geom_boxplot(aes(fill=name))+
  theme_bw() +
  scale_fill_manual(values = c("#8DA0CB","#E5E9F3","#FC8D62","#FEE5DC","#66C2A5","#DCF1EB","#E78AC3","#F9E4F1"),
                    labels=c('B18504_control', 'B18504_drought', 'Cayenne_control', 
                             'Cayenne_drought', 'R99_control', 'R99_drought', 'Rosetta_control', 'Rosetta_drought'))+
  labs(y="Shannon Index", subtitle = "C *, L *, W x C *")+
  facet_wrap(~ loc)+
  scale_y_continuous(limits = c(0, 4))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.text = element_text(size=25),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                                               size = 23,
                                               lineheight = 1,
                                               padding = margin(5.5, 5.5, 5.5, 5.5),
                                               margin = margin(0, 0, 5.5, 0),
                                               linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))
SHA.norm.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
# pairwise t-test
pwc_t.sha.norm <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(SHA.norm.orderNorm ~ treatment, p.adjust.method = "BH")
# add x y position
sha.norm.xy <- pwc_t.sha.norm  %>% 
  add_xy_position(x = "var", dodge = 0.8) # 
# plotting the pairwise comparisons among treatments 
SHA.norm.plot2 <- SHA.norm.plot + 
  stat_pvalue_manual(sha.norm.xy, x = "var", y.position = 3.8,
                     label = "{p.adj}{p.adj.signif}",size=6, hide.ns = T)
#tip.length = 0.01, hide.ns = F)
SHA.norm.plot2
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("SHA.tiff",
       SHA.norm.plot2, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")


#--------EVENNESS--------

#box plot
even.norm.box <- ggboxplot(
  join.meta.data.ed, x = "var", y = "EVE_norm", 
  color = "treatment", palette = "jco", facet.by = "loc")
even.norm.box
#linear mixed model
even.norm.lmer <- lmerTest::lmer(EVE_norm ~ treatment*var*loc + (1|block:loc) , 
                                 data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(even.norm.lmer, test="F", type="III") # 
# Assumption tests 
eve.norm.resid <- resid(even.norm.lmer)
shapiro.test(eve.norm.resid) # not normal
hist(eve.norm.resid) # not normal
plot(simulateResiduals(even.norm.lmer)) # 
leveneTest(residuals(even.norm.lmer) ~ join.meta.data.ed$treatment*join.meta.data.ed$var*join.meta.data.ed$loc) #okay 0.0677
eve.norm.Lave <- join.meta.data.ed %>% levene_test(EVE_norm ~ treatment*var*loc) 
eve.norm.Lave # not homogen
# Data transformation
#transform data
BNobject.eve.norm <- bestNormalize(join.meta.data.ed$EVE_norm)
BNobject.eve.norm #orderNorm Transformation
EVE.norm.orderNorm <- orderNorm(join.meta.data.ed$EVE_norm)
join.meta.data.ed$EVE.norm.orderNorm <- EVE.norm.orderNorm$x.t
# use box cox transform
#join.meta.data.ed <- join.meta.data.ed %>% select(-log10.SHA_norm)
tran.eve <- car::powerTransform(join.meta.data.ed$EVE_norm)
box.EVE_norm <- bcPower(join.meta.data.ed$EVE_norm, coef(tran.eve))
join.meta.data.ed$box.EVE_norm <- box.EVE_norm
str(join.meta.data.ed)
#linear mixed model
even.norm.lmer2 <- lmerTest::lmer(EVE.norm.orderNorm ~ treatment*var*loc + (1|block:loc) , 
                                 data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(even.norm.lmer2, test="F", type="III") # 
# Assumption tests 
eve.norm.resid2 <- resid(even.norm.lmer2)
shapiro.test(eve.norm.resid2) # normal 0.38
hist(eve.norm.resid2) # not normal
plot(simulateResiduals(even.norm.lmer2)) # 
leveneTest(residuals(even.norm.lmer2) ~ join.meta.data.ed$treatment*join.meta.data.ed$var*join.meta.data.ed$loc) #okay 0.257
eve.norm.Lave <- join.meta.data.ed %>% levene_test(EVE.norm.orderNorm ~ treatment*var*loc) 
eve.norm.Lave # good 0.257
# pairwise t-test
pwc_t.eve.norm <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(EVE.norm.orderNorm ~ treatment, p.adjust.method = "BH")
#select(-df, -statistic)# Remove details
pwc_t.eve.norm #

# EVENNESS Plot

EVE.norm.plot <- ggplot(join.meta.data.ed , aes(x=var, y=EVE_norm)) +
  geom_boxplot(aes(fill=name))+
  theme_bw() +
  scale_fill_manual(values = c("#8DA0CB","#E5E9F3","#FC8D62","#FEE5DC","#66C2A5","#DCF1EB","#E78AC3","#F9E4F1"),
                    labels=c('B18504_control', 'B18504_drought', 'Cayenne_control', 
                             'Cayenne_drought', 'R99_control', 'R99_drought', 'Rosetta_control', 'Rosetta_drought'))+
  labs(y="Pielou's Evenness")+
  facet_wrap(~ loc)+
  scale_y_continuous(limits = c(0, 0.85))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.text = element_text(size=25),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                                               size = 23,
                                               lineheight = 1,
                                               padding = margin(5.5, 5.5, 5.5, 5.5),
                                               margin = margin(0, 0, 5.5, 0),
                                               linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))
EVE.norm.plot
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("EVE.tiff",
       EVE.norm.plot, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")

#--------SIMPSON--------


#box plot
simp.norm.box <- ggboxplot(
  join.meta.data.ed, x = "var", y = "SIM_norm", 
  color = "treatment", palette = "jco", facet.by = "loc")
simp.norm.box
#linear mixed model
simp.norm.lmer <- lmerTest::lmer(SIM_norm ~ treatment*var*loc + (1|block:loc), 
                                 data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(simp.norm.lmer, test="F", type="III")
#Assumption tests:
simp.norm.resid <- resid(simp.norm.lmer)
shapiro.test(simp.norm.resid) # not normal
hist(simp.norm.resid) #not normal distributed
plot(simulateResiduals(simp.norm.lmer)) # violated
leveneTest(residuals(simp.norm.lmer) ~ join.meta.data.ed$treatment*join.meta.data.ed$var*join.meta.data.ed$loc) #okay 0.0688
simp.norm.Lave <- join.meta.data.ed %>% levene_test(SIM_norm ~ treatment*var*loc) # fine
simp.norm.Lave


#--------INVERSE SIMPSON--------


#box plot
inv.norm.box <- ggboxplot(
  join.meta.data.ed, x = "var", y = "INSIM_norm", 
  color = "treatment", palette = "jco", facet.by = "loc")
inv.norm.box
# linear mixed model
invsimp.norm.lmer <- lmerTest::lmer(INSIM_norm ~ treatment*var*loc + (1|block:loc), 
                                    data = join.meta.data.ed, contrasts = list(treatment="contr.sum",var="contr.sum",loc="contr.sum"))
car::Anova(invsimp.norm.lmer, test="F", type="III")
#Assumption tests:
inv.norm.resid <- resid(invsimp.norm.lmer)
shapiro.test(inv.norm.resid) # not normal
hist(inv.norm.resid) #not normal distributed
plot(simulateResiduals(invsimp.norm.lmer)) # violated
leveneTest(residuals(invsimp.norm.lmer) ~ join.meta.data.ed$treatment*join.meta.data.ed$var*join.meta.data.ed$loc) #okay 0.0688
inv.norm.Lave <- join.meta.data.ed %>% levene_test(INSIM_norm ~ treatment*var*loc) # fine
inv.norm.Lave
#transform data
#join.meta.data.ed <- join.meta.data.ed %>% select(-INSIM_box)
# NO TRANSFORMATIONS ARE SUCCESFULLY IMPLEMENTED
# USE NON PARAMETRIC TEST
#install.packages("ARTool")
library("ARTool")
insim.norm.lmer.m <- art(INSIM_norm ~ treatment*var*loc + (1|block:loc), data = join.meta.data.ed)
anova(insim.norm.lmer.m)
# pairwise t-test
pwc_t.insim.norm <- join.meta.data.ed %>%
  dplyr::group_by(var, loc) %>%
  pairwise_t_test(sqrt.INSIM ~ treatment, p.adjust.method = "BH")
#select(-df, -statistic)# Remove details
pwc_t.insim.norm #


#######################################################################################################################
## BETA DIVERSITY WHOLE DATASET
#######################################################################################################################

# Calculating dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)

# 1. Rarefied

rare.asv.df <-  as.data.frame(otu_table(decontaminated_rare_physeq))
sort(rowSums(rare.asv.df, na.rm = FALSE, dims = 1), decreasing = FALSE)
dim(rare.asv.df)
colnames(rare.asv.df)
# a. Bray-Curtis 
rare.asv_dist_bc <- vegdist(t(rare.asv.df), method = "bray")
# b. Jaccard 
rare.asv_dist_jc <- vegdist(t(rare.asv.df), binary = TRUE, method = "jaccard")
# c. Weighted Unifrac
rare.asv_dist_wUF <- UniFrac(decontaminated_rare_physeq, weighted=TRUE, normalized = TRUE)
# d. Unweighted UniFrac 
rare.asv_dist_uwUF <- UniFrac(decontaminated_rare_physeq, weighted=FALSE, normalized = TRUE)

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis 
rare.asv_pcoa_bc <- cmdscale(rare.asv_dist_bc, eig=T)
# Jaccard 
rare.asv_pcoa_jc <- cmdscale(rare.asv_dist_jc, eig=T)
# Weighted UniFrac 
rare.asv_pcoa_wUF <- cmdscale(rare.asv_dist_wUF, eig=T)
# Unweighted UniFrac 
rare.asv_pcoa.uwUF <- cmdscale(rare.asv_dist_uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis 
ax1.scores.rare <- rare.asv_pcoa_bc$points[,1]
ax2.scores.rare <- rare.asv_pcoa_bc$points[,2] 
# Jaccard 
ax1.scores.j.rare <- rare.asv_pcoa_jc$points[,1]
ax2.scores.j.rare <- rare.asv_pcoa_jc$points[,2]
# Weighted UniFrac 
ax1.scores.wUF.rare <- rare.asv_pcoa_wUF$points[,1]
ax2.scores.wUF.rare <- rare.asv_pcoa_wUF$points[,2]
# Unweighted UniFrac 
ax1.scores.uwUF.rare <- rare.asv_pcoa.uwUF$points[,1]
ax2.scores.uwUF.rare <- rare.asv_pcoa.uwUF$points[,2]

# 4. calculate percent variance explained, then add to plot
map.pcoa.rare <- join.meta.data.ed[-25,]
map.pcoa.rare <- droplevels(map.pcoa.rare)
str(map.pcoa.rare)
# Bray-curtis 
ax1.rare.BC <- rare.asv_pcoa_bc$eig[1]/sum(rare.asv_pcoa_bc$eig)
ax2.rare.BC <- rare.asv_pcoa_bc$eig[2]/sum(rare.asv_pcoa_bc$eig)
rare.BC.PCOA.map <- cbind(map.pcoa.rare,ax1.rare.BC,ax2.rare.BC)
# Jaccard 
ax1.rare.JC <- rare.asv_pcoa_jc$eig[1]/sum(rare.asv_pcoa_jc$eig)
ax2.rare.JC <- rare.asv_pcoa_jc$eig[2]/sum(rare.asv_pcoa_jc$eig)
rare.JC.PCOA.map <- cbind(map.pcoa.rare,ax1.rare.JC,ax2.rare.JC)
# Weighted UniFrac 
ax1.rare.wUF <- rare.asv_pcoa_wUF$eig[1]/sum(rare.asv_pcoa_wUF$eig)
ax2.rare.wUF <- rare.asv_pcoa_wUF$eig[2]/sum(rare.asv_pcoa_wUF$eig)
rare.wUF.PCOA.map <- cbind(map.pcoa.rare,ax1.rare.wUF,ax2.rare.wUF)
# Unweighted UniFrac 
ax1.rare.uwUF <- rare.asv_pcoa.uwUF$eig[1]/sum(rare.asv_pcoa.uwUF$eig)
ax2.rare.uwUF <- rare.asv_pcoa.uwUF$eig[2]/sum(rare.asv_pcoa.uwUF$eig)
rare.uwUF.PCOA.map <- cbind(map.pcoa.rare,ax1.rare.uwUF,ax2.rare.uwUF)

###############################################################################
# 6. PCoA Plot Rarefied  
require("ggrepel")
library(ggrepel)

# a. Bray-Curtis:
set.seed(33)

rare.PCOA.BC.plot <- ggplot(data = rare.BC.PCOA.map, aes(x=ax1.scores.rare, y=ax2.scores.rare, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.rare.BC,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.rare.BC,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "A. Bray Curtis (rarefied)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
  #stat_ellipse(aes(colour = var))
rare.PCOA.BC.plot

# b. Jaccard:
set.seed(33)

rare.PCOA.JC.plot <- ggplot(data = rare.JC.PCOA.map, aes(x=ax1.scores.j.rare, y=ax2.scores.j.rare, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.rare.JC,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.rare.JC,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "B. Jaccard (rarefied)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
rare.PCOA.JC.plot

# C. Weighted UniFrac:
set.seed(33)

rare.PCOA.WuF.plot <- ggplot(data = rare.wUF.PCOA.map, aes(x=ax1.scores.wUF.rare, y=ax2.scores.wUF.rare, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.rare.wUF,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.rare.wUF,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "C. Weighted UniFrac (rarefied)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
rare.PCOA.WuF.plot

# D. UnWeighted UniFrac:
set.seed(33)

rare.PCOA.UWuF.plot <- ggplot(data = rare.uwUF.PCOA.map, aes(x=ax1.scores.uwUF.rare, y=ax2.scores.uwUF.rare, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.rare.uwUF,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.rare.uwUF,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "C. Weighted UniFrac (rarefied)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
rare.PCOA.UWuF.plot

# Using Adonis2 package 

set.seed(13)
rare.adonis.BC <- adonis2(rare.asv_dist_bc ~ treatment*var*loc,
                              data=map.pcoa.rare, 
                              permutations = 9999) # significant
rare.adonis.BC 
# Df SumOfSqs      R2      F Pr(>F)   
# treatment          1  0.10579 0.03438 2.4300 0.0047 **
# var                3  0.18697 0.06075 1.4316 0.1429   
# loc                1  0.08373 0.02721 1.9233 0.0588 . 
# treatment:var      3  0.18655 0.06062 1.4284 0.1065   
# treatment:loc      1  0.06547 0.02127 1.5039 0.2126   
# var:loc            3  0.22998 0.07473 1.7609 0.0246 * 
# treatment:var:loc  3  0.21637 0.07031 1.6567 0.0403 * 

set.seed(13)
rare.adonis.JC <- adonis2(rare.asv_dist_jc ~ treatment*var*loc,
                          data=map.pcoa.rare, 
                          permutations = 9999) # just location is significant
rare.adonis.JC 
# Df SumOfSqs      R2      F Pr(>F)  
# treatment          1   0.2991 0.01966 1.2070 0.0808 .
# var                3   0.7770 0.05108 1.0452 0.2736  
# loc                1   0.3102 0.02039 1.2517 0.0493 *
# treatment:var      3   0.7528 0.04949 1.0126 0.4206  
# treatment:loc      1   0.2631 0.01730 1.0619 0.2970  
# var:loc            3   0.7395 0.04862 0.9948 0.4966  
# treatment:var:loc  3   0.6716 0.04415 0.9034 0.8833 

set.seed(13)
rare.adonis.WuF <- adonis2(rare.asv_dist_wUF ~ treatment*var*loc,
                          data=map.pcoa.rare, 
                          permutations = 9999) # 
rare.adonis.WuF 
# Df SumOfSqs      R2      F Pr(>F)   
# treatment          1  0.00632 0.01580 1.1431 0.3610   
# var                3  0.02860 0.07148 1.7236 0.0487 * 
# loc                1  0.01620 0.04049 2.9289 0.0040 **
# treatment:var      3  0.02934 0.07333 1.7680 0.0317 * 
# treatment:loc      1  0.01013 0.02532 1.8312 0.1037   
# var:loc            3  0.02687 0.06717 1.6196 0.0542 . 
# treatment:var:loc  3  0.02819 0.07046 1.6988 0.0419 * 

set.seed(13)
rare.adonis.UWuF <- adonis2(rare.asv_dist_uwUF ~ treatment*var*loc,
                           data=map.pcoa.rare, 
                           permutations = 9999) #
rare.adonis.UWuF
# Df SumOfSqs      R2      F Pr(>F)  
# treatment          1   0.1500 0.01838 1.1402 0.2588  
# var                3   0.3936 0.04822 0.9972 0.4776  
# loc                1   0.1805 0.02211 1.3717 0.0825 .
# treatment:var      3   0.4951 0.06064 1.2540 0.0522 .
# treatment:loc      1   0.1714 0.02099 1.3022 0.1177  
# var:loc            3   0.4143 0.05075 1.0496 0.3428  
# treatment:var:loc  3   0.3061 0.03749 0.7753 0.9499 

############################################################################################################################################

# 2. Normalized
library(phyloseq)
norm.asv.df <-  as.data.frame(otu_table(decontaminated_normal_physeq))
sort(rowSums(norm.asv.df, na.rm = FALSE, dims = 1), decreasing = FALSE)
dim(norm.asv.df) # 805, 63
colnames(norm.asv.df)
# a. Bray-Curtis 
norm.asv_dist_bc <- vegdist(t(norm.asv.df), method = "bray")
# b. Jaccard 
norm.asv_dist_jc <- vegdist(t(norm.asv.df), binary = TRUE, method = "jaccard")
# c. Weighted Unifrac
norm.asv_dist_wUF <- UniFrac(decontaminated_normal_physeq, weighted=TRUE, normalized = TRUE)
# d. Unweighted UniFrac 
norm.asv_dist_uwUF <- UniFrac(decontaminated_normal_physeq, weighted=FALSE, normalized = TRUE)

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis 
norm.asv_pcoa_bc <- cmdscale(norm.asv_dist_bc, eig=T)
# Jaccard 
norm.asv_pcoa_jc <- cmdscale(norm.asv_dist_jc, eig=T)
# Weighted UniFrac 
norm.asv_pcoa_wUF <- cmdscale(norm.asv_dist_wUF, eig=T)
# Unweighted UniFrac 
norm.asv_pcoa.uwUF <- cmdscale(norm.asv_dist_uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis 
ax1.scores.norm <- norm.asv_pcoa_bc$points[,1]
ax2.scores.norm <- norm.asv_pcoa_bc$points[,2] 
# Jaccard 
ax1.scores.j.norm <- norm.asv_pcoa_jc$points[,1]
ax2.scores.j.norm <- norm.asv_pcoa_jc$points[,2]
# Weighted UniFrac 
ax1.scores.wUF.norm <- norm.asv_pcoa_wUF$points[,1]
ax2.scores.wUF.norm <- norm.asv_pcoa_wUF$points[,2]
# Unweighted UniFrac 
ax1.scores.uwUF.norm <- norm.asv_pcoa.uwUF$points[,1]
ax2.scores.uwUF.norm <- norm.asv_pcoa.uwUF$points[,2]

# 4. calculate percent variance explained, then add to plot
map.pcoa.norm <- join.meta.data.ed

# Bray-curtis 
ax1.norm.BC <- norm.asv_pcoa_bc$eig[1]/sum(norm.asv_pcoa_bc$eig)
ax2.norm.BC <- norm.asv_pcoa_bc$eig[2]/sum(norm.asv_pcoa_bc$eig)
norm.BC.PCOA.map <- cbind(map.pcoa.norm,ax1.norm.BC,ax2.norm.BC)
# Jaccard 
ax1.norm.JC <- norm.asv_pcoa_jc$eig[1]/sum(norm.asv_pcoa_jc$eig)
ax2.norm.JC <- norm.asv_pcoa_jc$eig[2]/sum(norm.asv_pcoa_jc$eig)
norm.JC.PCOA.map <- cbind(map.pcoa.norm,ax1.norm.JC,ax2.norm.JC)
# Weighted UniFrac 
ax1.norm.wUF <- norm.asv_pcoa_wUF$eig[1]/sum(norm.asv_pcoa_wUF$eig)
ax2.norm.wUF <- norm.asv_pcoa_wUF$eig[2]/sum(norm.asv_pcoa_wUF$eig)
norm.wUF.PCOA.map <- cbind(map.pcoa.norm,ax1.norm.wUF,ax2.norm.wUF)
# Unweighted UniFrac 
ax1.norm.uwUF <- norm.asv_pcoa.uwUF$eig[1]/sum(norm.asv_pcoa.uwUF$eig)
ax2.norm.uwUF <- norm.asv_pcoa.uwUF$eig[2]/sum(norm.asv_pcoa.uwUF$eig)
norm.uwUF.PCOA.map <- cbind(map.pcoa.norm,ax1.norm.uwUF,ax2.norm.uwUF)


# Using Adonis2 package 

set.seed(13)
norm.adonis.BC <- adonis2(norm.asv_dist_bc ~ treatment*var*loc,
                          data=map.pcoa.norm, 
                          permutations = 9999) # significant
norm.adonis.BC 
# Df SumOfSqs      R2      F Pr(>F)
# treatment          1   0.0948 0.02373 1.5539 0.1262
# var                3   0.2482 0.06211 1.3557 0.1277
# loc                1   0.0652 0.01631 1.0679 0.3573
# treatment:var      3   0.1854 0.04640 1.0129 0.4275
# treatment:loc      1   0.0828 0.02072 1.3571 0.1972
# var:loc            3   0.2163 0.05411 1.1812 0.2364
# treatment:var:loc  3   0.2353 0.05889 1.2854 0.1662 

set.seed(13)
norm.adonis.JC <- adonis2(norm.asv_dist_jc ~ treatment*var*loc,
                          data=map.pcoa.norm, permutations = 9999) # just location is significant
norm.adonis.JC 
# Df SumOfSqs      R2      F Pr(>F)  
# treatment          1   0.3125 0.01929 1.2085 0.0822 .
# var                3   0.8178 0.05048 1.0540 0.2349  
# loc                1   0.2945 0.01818 1.1388 0.1465  
# treatment:var      3   0.8157 0.05035 1.0513 0.2570  
# treatment:loc      1   0.2636 0.01627 1.0192 0.4053  
# var:loc            3   0.8122 0.05013 1.0468 0.2730  
# treatment:var:loc  3   0.7292 0.04501 0.9398 0.7516 

set.seed(13)
norm.adonis.WuF <- adonis2(norm.asv_dist_wUF ~ treatment*var*loc,
                           data=map.pcoa.norm, 
                           permutations = 9999) # 
summary(norm.adonis.WuF)
# Df SumOfSqs      R2      F Pr(>F)   
# treatment          1  0.00681 0.01700 1.2716 0.2944   
# var                3  0.02995 0.07474 1.8638 0.0297 * 
# loc                1  0.01578 0.03936 2.9448 0.0031 **
# treatment:var      3  0.03024 0.07545 1.8815 0.0165 * 
# treatment:loc      1  0.01097 0.02738 2.0481 0.0580 . 
# var:loc            3  0.02650 0.06614 1.6492 0.0520 . 
# treatment:var:loc  3  0.02872 0.07168 1.7874 0.0268 * 
#______________________

# subset



set.seed(13)
norm.adonis.WuF2 <- adonis2(norm.asv_dist_wUF ~ name,
                           data=map.pcoa.norm, 
                           permutations = 9999) # 
norm.adonis.WuF2



set.seed(13)
norm.adonis.UWuF <- adonis2(norm.asv_dist_uwUF ~ treatment*var*loc,
                            data=map.pcoa.norm, 
                            permutations = 9999) #
norm.adonis.UWuF
# Df SumOfSqs      R2      F Pr(>F)  
# treatment          1   0.1394 0.01563 0.9992 0.4560  
# var                3   0.4430 0.04966 1.0584 0.3202  
# loc                1   0.2000 0.02242 1.4334 0.0639 .
# treatment:var      3   0.5492 0.06157 1.3121 0.0359 *
# treatment:loc      1   0.1669 0.01871 1.1961 0.2003  
# var:loc            3   0.5457 0.06117 1.3036 0.0384 *
# treatment:var:loc  3   0.3182 0.03567 0.7601 0.9576   

###############################################################################
#  PCoA Plot Normalized 
require("ggrepel")
library(ggrepel)

# a. Bray-Curtis:
set.seed(33)

norm.PCOA.BC.plot <- ggplot(data = norm.BC.PCOA.map, aes(x=ax1.scores.norm, y=ax2.scores.norm, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.norm.BC,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.norm.BC,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "A. Bray Curtis (CSS-normalized)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = loc))
norm.PCOA.BC.plot

# b. Jaccard:
set.seed(33)

norm.PCOA.JC.plot <- ggplot(data = norm.JC.PCOA.map, aes(x=ax1.scores.j.norm, y=ax2.scores.j.norm, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                    labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.norm.JC,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.norm.JC,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "B. Jaccard (CSS-normalized)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
norm.PCOA.JC.plot

# C. Weighted UniFrac:
set.seed(33)

norm.PCOA.WuF.plot <- ggplot(data = norm.wUF.PCOA.map, aes(x=ax1.scores.wUF.norm, y=ax2.scores.wUF.norm, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 3,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar:",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime:",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#8DA0CB","#8DA0CB","#FC8D62", "#FC8D62", "#66C2A5", "#66C2A5","#E78AC3","#E78AC3")) +
  geom_mark_ellipse(aes(fill = name), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.norm.wUF,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.norm.wUF,3)*100,"% var. explained", sep=""))+
 # xlim(c(-0.3,0.3)) +
 # ylim(c(-0.15,0.15)) +
  # scale_x_continuous(limits = symmetric_limits) +
 # scale_y_continuous(limits = symmetric_limits) +
 #labs(subtitle = "Weighted UniFrac (CSS-normalized)")+

  theme(legend.position="bottom",
        legend.text = element_text(size=20),
        legend.title = element_text(size=22),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
norm.PCOA.WuF.plot

setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("PCOA.WuF.normal.tiff",
       norm.PCOA.WuF.plot, device = "tiff",
       width = 13, height = 10, 
       units= "in", dpi =300,compression = "lzw",  bg= "white")




# D. UnWeighted UniFrac:
set.seed(33)

norm.PCOA.UWuF.plot <- ggplot(data = norm.uwUF.PCOA.map, aes(x=ax1.scores.uwUF.norm, y=ax2.scores.uwUF.norm, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.norm.uwUF,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.norm.uwUF,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "C. Unweighted UniFrac (CSS-normalized)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
norm.PCOA.UWuF.plot


# code for background theme for ggplot so I don't have to have it in every plot
my_theme <- theme(panel.background = element_rect(fill = "white", colour = "white"), 
                  panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid', colour = "light gray"),
                  panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', colour = "light gray"))
seed_WU_ord <- ordinate(physeq = decontaminated_normal_physeq, method="PCoA", distance=norm.asv_dist_wUF)
seed_WU_ord 
# plot ordination
seed_WU_PCoA <- plot_ordination(decontaminated_normal_physeq, seed_WU_ord, color="loc", shape="treatment", ) + theme(aspect.ratio=1) + geom_point(size=3) + my_theme 
plot(seed_WU_PCoA)


# 3D PCoA Plot
#install.packages("dartR")
BiocManager::install("SNPRelate")
library(SNPRelate)
library(dartR)
gl.pcoa.plot.3d(
  norm.asv_pcoa_wUF,
  norm.asv_dist_wUF,
  xaxis = 1,
  yaxis = 2,
  zaxis = 3,
  radius = 8,
  verbose = NULL
)











#######################################################################################################################
## BETA DIVERSITY SUBSET BY LOCATION
#######################################################################################################################

# Calculating dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)

# 1. Rarefied

rare.asv.df <-  as.data.frame(otu_table(decontaminated_rare_physeq))
colnames(rare.asv.df)
dim(rare.asv.df)
rare.asv.LP <- rare.asv.df[,1:31]
sort(rowSums(rare.asv.LP, na.rm = FALSE, dims = 1), decreasing = FALSE)
rare.asv.LP1 <- rare.asv.LP[rowSums(rare.asv.LP)>0,]
sort(rowSums(rare.asv.LP1, na.rm = FALSE, dims = 1), decreasing = FALSE)

decon_rare_seq.LP <- subset_samples(decontaminated_rare_physeq, loc=="LP")
sort(rowSums(otu_table(decon_rare_seq.LP), na.rm = FALSE, dims = 1), decreasing = FALSE)
decon_rare_seq.LP1 <- prune_taxa(taxa_sums(decon_rare_seq.LP)>0, decon_rare_seq.LP)
sort(rowSums(otu_table(decon_rare_seq.LP1), na.rm = FALSE, dims = 1), decreasing = FALSE)

# a. Bray-Curtis 
rare.dist.LP.bc <- vegdist(t(rare.asv.LP1), method = "bray")
# b. Jaccard 
rare.dist.LP.jc <- vegdist(t(rare.asv.LP1), binary = TRUE, method = "jaccard")
# c. Weighted Unifrac
rare.dist.LP.WuF <- UniFrac(decon_rare_seq.LP1, weighted=TRUE, normalized = TRUE)
# d. Unweighted UniFrac 
rare.dist.LP.uwUF <- UniFrac(decon_rare_seq.LP1, weighted=FALSE, normalized = TRUE)

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis 
rare.pcoa.bc.LP <- cmdscale(rare.dist.LP.bc, eig=T)
# Jaccard 
rare.pcoa.jc.LP <- cmdscale(rare.dist.LP.jc, eig=T)
# Weighted UniFrac 
rare.pcoa.wUF.LP <- cmdscale(rare.dist.LP.WuF, eig=T)
# Unweighted UniFrac 
rare.pcoa.uwUF.LP <- cmdscale(rare.dist.LP.uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis 
ax1.scor.bc.rare.LP <- rare.pcoa.bc.LP$points[,1]
ax2.scor.bc.rare.LP <- rare.pcoa.bc.LP$points[,2] 
# Jaccard 
ax1.scor.j.rare.LP <- rare.pcoa.jc.LP$points[,1]
ax2.scor.j.rare.LP <- rare.pcoa.jc.LP$points[,2]
# Weighted UniFrac 
ax1.scor.wUF.rare.LP <- rare.pcoa.wUF.LP$points[,1]
ax2.scor.wUF.rare.LP <- rare.pcoa.wUF.LP$points[,2]
# Unweighted UniFrac 
ax1.scor.uwUF.rare.LP <- rare.pcoa.uwUF.LP$points[,1]
ax2.scor.uwUF.rare.LP <- rare.pcoa.uwUF.LP$points[,2]

# 4. calculate percent variance explained, then add to plot
map.pcoa.rare <- join.meta.data.ed[-25,]
map.pcoa.rare <- droplevels(map.pcoa.rare)
map.pcoa.rare.LP <- map.pcoa.rare[1:31,]
map.pcoa.rare.LP <- droplevels(map.pcoa.rare.LP)
view(map.pcoa.rare.LP)
# Bray-curtis 
ax1.rare.BC.LP <- rare.pcoa.bc.LP$eig[1]/sum(rare.pcoa.bc.LP$eig)
ax2.rare.BC.LP <- rare.pcoa.bc.LP$eig[2]/sum(rare.pcoa.bc.LP$eig)
rare.BC.PCOA.map.LP <- cbind(map.pcoa.rare.LP,ax1.rare.BC.LP,ax2.rare.BC.LP)
# Jaccard 
ax1.rare.JC.LP <- rare.pcoa.jc.LP$eig[1]/sum(rare.pcoa.jc.LP$eig)
ax2.rare.JC.LP <- rare.pcoa.jc.LP$eig[2]/sum(rare.pcoa.jc.LP$eig)
rare.JC.PCOA.map.LP <- cbind(map.pcoa.rare.LP,ax1.rare.JC.LP,ax2.rare.JC.LP)
# Weighted UniFrac 
ax1.rare.wUF.LP <- rare.pcoa.wUF.LP$eig[1]/sum(rare.pcoa.wUF.LP$eig)
ax2.rare.wUF.LP <- rare.pcoa.wUF.LP$eig[2]/sum(rare.pcoa.wUF.LP$eig)
rare.wUF.PCOA.map.LP <- cbind(map.pcoa.rare.LP,ax1.rare.wUF.LP,ax2.rare.wUF.LP)
# Unweighted UniFrac 
ax1.rare.uwUF.LP <- rare.pcoa.uwUF.LP$eig[1]/sum(rare.pcoa.uwUF.LP$eig)
ax2.rare.uwUF.LP <- rare.pcoa.uwUF.LP$eig[2]/sum(rare.pcoa.uwUF.LP$eig)
rare.uwUF.PCOA.map.LP <- cbind(map.pcoa.rare.LP,ax1.rare.uwUF.LP,ax2.rare.uwUF.LP)


# PCOA PLOT RAREFIED - LP

# a. Bray-Curtis:
set.seed(33)

rare.PCOA.BC.plot.LP <- ggplot(data = rare.BC.PCOA.map.LP, aes(x=ax1.scor.bc.rare.LP, y=ax2.scor.bc.rare.LP, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#8DA0CB","#8DA0CB","#FC8D62", "#FC8D62", "#66C2A5", "#66C2A5","#E78AC3","#E78AC3")) +
  geom_mark_ellipse(aes(fill = name), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.rare.BC.LP,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.rare.BC.LP,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "A. Bray Curtis - Lower Peninsula (rarefied)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
rare.PCOA.BC.plot.LP

# b. Jaccard:
set.seed(33)

rare.PCOA.JC.plot.LP <- ggplot(data = rare.JC.PCOA.map.LP, aes(x=ax1.scor.j.rare.LP, y=ax2.scor.j.rare.LP, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#8DA0CB","#8DA0CB","#FC8D62", "#FC8D62", "#66C2A5", "#66C2A5","#E78AC3","#E78AC3")) +
  geom_mark_ellipse(aes(fill = name), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.rare.JC.LP,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.rare.JC.LP,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "B. Jaccard - Lower Peninsula (rarefied)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
rare.PCOA.JC.plot.LP

# C. Weighted UniFrac:
set.seed(33)

rare.PCOA.WuF.plot.LP <- ggplot(data = rare.wUF.PCOA.map.LP, aes(x=ax1.scores.wUF.rare, y=ax2.scores.wUF.rare, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.rare.wUF,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.rare.wUF,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "C. Weighted UniFrac (rarefied)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
rare.PCOA.WuF.plot.LP

# D. UnWeighted UniFrac:
set.seed(33)

rare.PCOA.UWuF.plot <- ggplot(data = rare.uwUF.PCOA.map, aes(x=ax1.scores.uwUF.rare, y=ax2.scores.uwUF.rare, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.rare.uwUF,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.rare.uwUF,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "C. Weighted UniFrac (rarefied)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
rare.PCOA.UWuF.plot

# 2. Normalized

# LOWER PENINSULA
norm.asv.df <-  as.data.frame(otu_table(decontaminated_normal_physeq))
norm.asv.LP <- norm.asv.df[,1:32]
sort(rowSums(norm.asv.LP, na.rm = FALSE, dims = 1), decreasing = FALSE)
norm.asv.LP1 <- norm.asv.LP[rowSums(norm.asv.LP)>0,]
sort(rowSums(norm.asv.LP1, na.rm = FALSE, dims = 1), decreasing = FALSE)

decon_norm_seq.LP <- subset_samples(decontaminated_normal_physeq, loc=="LP")
sort(rowSums(otu_table(decon_norm_seq.LP), na.rm = FALSE, dims = 1), decreasing = FALSE)
decon_norm_seq.LP1 <- prune_taxa(taxa_sums(decon_norm_seq.LP)>0, decon_norm_seq.LP)
sort(rowSums(otu_table(decon_norm_seq.LP1), na.rm = FALSE, dims = 1), decreasing = FALSE)

# a. Bray-Curtis 
norm.dist.LP.bc <- vegan::vegdist(t(norm.asv.LP1), method = "bray")
norm.dist.LP.bc
# b. Jaccard 
norm.dist.LP.jc <- vegan::vegdist(t(norm.asv.LP1), binary = TRUE, method = "jaccard")
norm.dist.LP.jc
# c. Weighted Unifrac
norm.dist.LP.WuF <- UniFrac(decon_norm_seq.LP1, weighted=TRUE, normalized = TRUE)
# d. Unweighted UniFrac 
norm.dist.LP.uwUF <- UniFrac(decon_norm_seq.LP1, weighted=FALSE, normalized = TRUE)

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis 
norm.pcoa.bc.LP <- cmdscale(norm.dist.LP.bc, eig=T)
# Jaccard 
norm.pcoa.jc.LP <- cmdscale(norm.dist.LP.jc, eig=T)
# Weighted UniFrac 
norm.pcoa.wUF.LP <- cmdscale(norm.dist.LP.WuF, eig=T)
# Unweighted UniFrac 
norm.pcoa.uwUF.LP <- cmdscale(norm.dist.LP.uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis 
ax1.scor.bc.norm.LP <- norm.pcoa.bc.LP$points[,1]
ax2.scor.bc.norm.LP <- norm.pcoa.bc.LP$points[,2] 
# Jaccard 
ax1.scor.j.norm.LP <- norm.pcoa.jc.LP$points[,1]
ax2.scor.j.norm.LP <- norm.pcoa.jc.LP$points[,2]
# Weighted UniFrac 
ax1.scor.wUF.norm.LP <- norm.pcoa.wUF.LP$points[,1]
ax2.scor.wUF.norm.LP <- norm.pcoa.wUF.LP$points[,2]
# Unweighted UniFrac 
ax1.scor.uwUF.norm.LP <- norm.pcoa.uwUF.LP$points[,1]
ax2.scor.uwUF.norm.LP <- norm.pcoa.uwUF.LP$points[,2]

# 4. calculate percent variance explained, then add to plot
map.pcoa.norm <- join.meta.data.ed
map.pcoa.norm.LP <- map.pcoa.norm[1:32,]
view(map.pcoa.norm.LP)
map.pcoa.norm.LP <- droplevels(map.pcoa.norm.LP)
# Bray-curtis 
ax1.norm.BC.LP <- norm.pcoa.bc.LP$eig[1]/sum(norm.pcoa.bc.LP$eig)
ax2.norm.BC.LP <- norm.pcoa.bc.LP$eig[2]/sum(norm.pcoa.bc.LP$eig)
norm.BC.PCOA.map.LP <- cbind(map.pcoa.norm.LP,ax1.norm.BC.LP,ax2.norm.BC.LP)
# Jaccard 
ax1.norm.JC.LP <- norm.pcoa.jc.LP$eig[1]/sum(norm.pcoa.jc.LP$eig)
ax2.norm.JC.LP <- norm.pcoa.jc.LP$eig[2]/sum(norm.pcoa.jc.LP$eig)
norm.JC.PCOA.map.LP <- cbind(map.pcoa.norm.LP,ax1.norm.JC.LP,ax2.norm.JC.LP)
# Weighted UniFrac 
ax1.norm.wUF.LP <- norm.pcoa.wUF.LP$eig[1]/sum(norm.pcoa.wUF.LP$eig)
ax2.norm.wUF.LP <- norm.pcoa.wUF.LP$eig[2]/sum(norm.pcoa.wUF.LP$eig)
norm.wUF.PCOA.map.LP <- cbind(map.pcoa.norm.LP,ax1.norm.wUF.LP,ax2.norm.wUF.LP)
# Unweighted UniFrac 
ax1.norm.uwUF.LP <- norm.pcoa.uwUF.LP$eig[1]/sum(norm.pcoa.uwUF.LP$eig)
ax2.norm.uwUF.LP <- norm.pcoa.uwUF.LP$eig[2]/sum(norm.pcoa.uwUF.LP$eig)
norm.uwUF.PCOA.map.LP <- cbind(map.pcoa.norm.LP,ax1.norm.uwUF.LP,ax2.norm.uwUF.LP)


set.seed(33)

norm.PCOA.WuF.plot.LP <- ggplot(data = norm.wUF.PCOA.map.LP, aes(x=ax1.scor.wUF.norm.LP, y=ax2.scor.wUF.norm.LP, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.norm.wUF.LP,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.norm.wUF.LP,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "Weighted UniFrac: Lower Peninsula (CSS-normalized)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
norm.PCOA.WuF.plot.LP

# UPPER PENINSULA

norm.asv.df <-  as.data.frame(otu_table(decontaminated_normal_physeq))
norm.asv.UP <- norm.asv.df[,33:63]
sort(rowSums(norm.asv.UP, na.rm = FALSE, dims = 1), decreasing = FALSE)
norm.asv.UP1 <- norm.asv.UP[rowSums(norm.asv.UP)>0,]
sort(rowSums(norm.asv.UP1, na.rm = FALSE, dims = 1), decreasing = FALSE)

decon_norm_seq.UP <- subset_samples(decontaminated_normal_physeq, loc=="UP")
sort(rowSums(otu_table(decon_norm_seq.UP), na.rm = FALSE, dims = 1), decreasing = FALSE)
decon_norm_seq.UP1 <- prune_taxa(taxa_sums(decon_norm_seq.UP)>0, decon_norm_seq.UP)
sort(rowSums(otu_table(decon_norm_seq.UP1), na.rm = FALSE, dims = 1), decreasing = FALSE)

# a. Bray-Curtis 
norm.dist.UP.bc <- vegan::vegdist(t(norm.asv.UP1), method = "bray")
norm.dist.UP.bc
# b. Jaccard 
norm.dist.UP.jc <- vegan::vegdist(t(norm.asv.UP1), binary = TRUE, method = "jaccard")
norm.dist.UP.jc
# c. Weighted Unifrac
norm.dist.UP.WuF <- UniFrac(decon_norm_seq.UP1, weighted=TRUE, normalized = TRUE)
# d. Unweighted UniFrac 
norm.dist.UP.uwUF <- UniFrac(decon_norm_seq.UP1, weighted=FALSE, normalized = TRUE)

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis 
norm.pcoa.bc.UP <- cmdscale(norm.dist.UP.bc, eig=T)
# Jaccard 
norm.pcoa.jc.UP <- cmdscale(norm.dist.UP.jc, eig=T)
# Weighted UniFrac 
norm.pcoa.wUF.UP <- cmdscale(norm.dist.UP.WuF, eig=T)
# Unweighted UniFrac 
norm.pcoa.uwUF.UP <- cmdscale(norm.dist.UP.uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis 
ax1.scor.bc.norm.UP <- norm.pcoa.bc.UP$points[,1]
ax2.scor.bc.norm.UP <- norm.pcoa.bc.UP$points[,2] 
# Jaccard 
ax1.scor.j.norm.UP <- norm.pcoa.jc.UP$points[,1]
ax2.scor.j.norm.UP <- norm.pcoa.jc.UP$points[,2]
# Weighted UniFrac 
ax1.scor.wUF.norm.UP <- norm.pcoa.wUF.UP$points[,1]
ax2.scor.wUF.norm.UP <- norm.pcoa.wUF.UP$points[,2]
# Unweighted UniFrac 
ax1.scor.uwUF.norm.UP <- norm.pcoa.uwUF.UP$points[,1]
ax2.scor.uwUF.norm.UP <- norm.pcoa.uwUF.UP$points[,2]

# 4. calculate percent variance explained, then add to plot
map.pcoa.norm <- join.meta.data.ed
map.pcoa.norm.UP <- map.pcoa.norm[33:63,]
map.pcoa.norm.UP <- droplevels(map.pcoa.norm.UP)
str(map.pcoa.norm.UP)
view(map.pcoa.norm.UP)
# Bray-curtis 
ax1.norm.BC.UP <- norm.pcoa.bc.UP$eig[1]/sum(norm.pcoa.bc.UP$eig)
ax2.norm.BC.UP <- norm.pcoa.bc.UP$eig[2]/sum(norm.pcoa.bc.UP$eig)
norm.BC.PCOA.map.UP <- cbind(map.pcoa.norm.UP,ax1.norm.BC.UP,ax2.norm.BC.UP)
# Jaccard 
ax1.norm.JC.UP <- norm.pcoa.jc.UP$eig[1]/sum(norm.pcoa.jc.UP$eig)
ax2.norm.JC.UP <- norm.pcoa.jc.UP$eig[2]/sum(norm.pcoa.jc.UP$eig)
norm.JC.PCOA.map.UP <- cbind(map.pcoa.norm.UP,ax1.norm.JC.UP,ax2.norm.JC.UP)
# Weighted UniFrac 
ax1.norm.wUF.UP <- norm.pcoa.wUF.UP$eig[1]/sum(norm.pcoa.wUF.UP$eig)
ax2.norm.wUF.UP <- norm.pcoa.wUF.UP$eig[2]/sum(norm.pcoa.wUF.UP$eig)
norm.wUF.PCOA.map.UP <- cbind(map.pcoa.norm.UP,ax1.norm.wUF.UP,ax2.norm.wUF.UP)
# Unweighted UniFrac 
ax1.norm.uwUF.UP <- norm.pcoa.uwUF.UP$eig[1]/sum(norm.pcoa.uwUF.UP$eig)
ax2.norm.uwUF.UP <- norm.pcoa.uwUF.UP$eig[2]/sum(norm.pcoa.uwUF.UP$eig)
norm.uwUF.PCOA.map.UP <- cbind(map.pcoa.norm.UP,ax1.norm.uwUF.UP,ax2.norm.uwUF.UP)


set.seed(33)

norm.PCOA.WuF.plot.UP <- ggplot(data = norm.wUF.PCOA.map.UP, aes(x=ax1.scor.wUF.norm.UP, y=ax2.scor.wUF.norm.UP, colour=var))+
  theme_bw()+
  geom_point(aes(color = var, shape = treatment), size = 1.5,stroke=1) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.norm.wUF.UP,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.norm.wUF.UP,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "Weighted UniFrac: Upper Peninsula (CSS-normalized)")+
  theme(legend.position="right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))
#stat_ellipse(aes(colour = var))
norm.PCOA.WuF.plot.UP

library(BiodiversityR)
set.seed(13)
norm.WuF.npmanova <- nested.npmanova(norm.dist.UP.WuF~var+treatment, data=map.pcoa.norm.UP, 
                permutations=999)
norm.WuF.npmanova #Nested anova for treatment nested within var 


perm <- how(within = Within(type = "free"),
                   plots = Plots(strata = map.pcoa.norm.UP$block, type = "none"),
                   nperm = 9999,
                   observed = TRUE)
#they specify that plots are to be freely permuted within blocks but that blocks are not allowed to permute
set.seed(13)
norm.WuF.adon.strata.UP <- adonis2(norm.dist.UP.WuF ~ var*treatment, data=map.pcoa.norm.UP, 
                                      permutations = perm)
norm.WuF.adon.strata.UP

norm.WuF.adon.UP <- adonis2(norm.dist.UP.WuF ~ var*treatment, data=map.pcoa.norm.UP,permutations=999)
norm.WuF.adon.UP




set.seed(13)
norm.WuF.npmanova <- nested.npmanova(norm.dist.UP.WuF~var+treatment, data=map.pcoa.norm.UP, 
                                     permutations=9999)
norm.WuF.npmanova













##########################################################################################################################################
# CAP ANALYSIS

#install.packages("parallel")
#install.packages("ggforce")
#install.packages("BiodiversityR")
library(parallel)
library(BiodiversityR) # ALWAYS LOAD FROM THE R CONSOLE!!!!
library(ggforce)

### 1. Rarefied

### 1 A. Bray-Curtis
rare.asv_dist_bc 
# metadata
str(map.pcoa.rare)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(rare.asv_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(rare.asv_dist_bc ~ treatment, data = map.pcoa.rare, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}
par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)
success$class.success
#here for example you would choose 9 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
cap.rare.BC <- CAPdiscrim(rare.asv_dist_bc ~ treatment, data = map.pcoa.rare, m = 9, permutations = 9999, add = TRUE) # 70.97% ,Significance of this percentage was 0.0007

success <- cbind(data.frame(cap.rare.BC$group), data.frame(cap.rare.BC$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(cap.rare.BC$PCoA)
success <- success[order(success$source), ]
success

rare.cap1.bc <- paste("CAP1 (", round((100/sum(cap.rare.BC$lda.other$svd^2) * cap.rare.BC$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
rare.cap2.bc <- paste("CAP2 (", round((100/sum(cap.rare.BC$lda.other$svd^2) * cap.rare.BC$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with ggplot2

rare.cap.bc.plot <- 
  ggplot(as.data.frame(cap.rare.BC$x), aes(x = cap.rare.BC$x[,1], y = cap.rare.BC$x[,2])) +
  geom_point(aes(color = map.pcoa.rare$var, shape = map.pcoa.rare$treatment), size = 3, stroke=1) +
  xlab(rare.cap1.bc) + ylab(rare.cap2.bc) +
  scale_color_manual(values = c("#8DA0CB","#FC8D62","#66C2A5","#E78AC3"),
                     name = "Cultivar",
                     labels = c("B18504", "Cayenne", "R99","Rosetta")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Water regime",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#8DA0CB","#8DA0CB","#FC8D62", "#FC8D62", "#66C2A5", "#66C2A5","#E78AC3","#E78AC3")) +
  geom_mark_ellipse(aes(fill = map.pcoa.rare$name), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  #labs(subtitle = "C. AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.subtitle = element_text(size=20, face="bold")) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=20)) +
  #annotate("text",x=-9.5,y=-1,label= "90%", hjust = 0, size = 4,color="#009E73") +
  #annotate("text",x=-7.5,y=3,label= "95%", hjust = 0, size = 4,color="#009E73") +
  #annotate("text",x=-1.7,y=-7,label= "100%", hjust = 0, size = 4, color="#FF618C")+
  #annotate("text",x=-1.7,y=-1.5,label= "100%", hjust = 0, size = 4,color="#FF618C")+
  #annotate("text",x=2,y=6,label= "95%", hjust = 0, size = 4, color="#E69F00")+
  #annotate("text",x=6.5,y=0.4,label= "85%", hjust = 0, size = 4,color="#E69F00")+
  #annotate("text",x=-14,y=-10,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 6) +
  #annotate("text", x=-14, y=-11.5, label= "Pillai's test=4.1***", hjust = 0, size = 6)+
  guides(colour=guide_legend(override.aes = list(size=7)),shape=guide_legend(override.aes = list(size=7)))
rare.cap.bc.plot
setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOA_CAP_bulk_bray.tiff",
       aoa.cap.plot, device = "tiff",
       width = 4, height =3, 
       units= "in", dpi = 600)

####### Relative Abundance

# Normalized phyloseq object
tax_table(decontaminated_normal_physeq)
# Order
norm.order <- tax_glom(decontaminated_normal_physeq, taxrank = "Order", NArm = F)
norm.order.ra <- transform_sample_counts(norm.order, function(x) x/sum(x))
norm.order.df <- as.data.frame(psmelt(norm.order.ra) %>%
  group_by(loc, var, treatment, Order) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean))
#norm.order.df$Order <- as.character(norm.order.df$Order)
norm.order.df$Order[norm.order.df$Mean < 0.005] <- "Other (<0.5%)"
norm.order.df$Order[is.na(norm.order.df$Order)] <- "Other (<0.5%)"
str(norm.order.df)
view(norm.order.df)
norm.order.df$var <- factor(norm.order.df$var)
norm.order.df$Order <- factor(norm.order.df$Order)
norm.order.df$treatment <- factor(norm.order.df$treatment, levels = c("open", "shelter"),
                                  labels = c("Control", "Drought"))
norm.order.df$loc <- factor(norm.order.df$loc, levels = c("LP", "UP"),
                                  labels = c("Lower Peninsula", "Upper Peninsula"))
# Genus
norm.gen <- tax_glom(decontaminated_normal_physeq, taxrank = "Genus", NArm = F)
norm.gen.ra <- transform_sample_counts(norm.gen, function(x) x/sum(x))
norm.gen.df <- as.data.frame(psmelt(norm.gen.ra) %>%
                               group_by(loc, var, treatment, Genus) %>%
                               summarize(Mean = mean(Abundance)) %>%
                               arrange(-Mean))
view(norm.gen.df)
norm.gen.df$Genus <- as.character(norm.gen.df$Genus)
norm.gen.df$Genus[norm.gen.df$Mean < 0.005] <- "Other<0.5%"
norm.gen.df$Genus[is.na(norm.gen.df$Genus)] <- "Other<0.5%"
norm.gen.df$var <- factor(norm.gen.df$var)
norm.gen.df$Genus <- factor(norm.gen.df$Genus)
norm.gen.df$treatment <- factor(norm.gen.df$treatment, levels = c("open", "shelter"),
                                  labels = c("Control", "Drought"))
norm.gen.df$loc <- factor(norm.gen.df$loc, levels = c("LP", "UP"),
                            labels = c("Lower Peninsula", "Upper Peninsula"))




library(RColorBrewer)
n <- 19
colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
col <- sample(col_vec, n)
area <- rep(1,n)
pie(area, col = col)
#brewer.pal(12, "Paired")
#display.brewer.pal(n)
#display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, colorblindFriendly=T)
#install.packages('pals')
library(pals)
pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome, 
          stepped, tol, watlington,
          show.names=F)
col.order <- kelly(n=22)


#install.packages("ggh4x")
library(ggh4x)
set.seed(13)
norm.order.plot <- ggplot(norm.order.df, aes(x=treatment, y=Mean, fill=Order)) + 
  geom_bar(aes(), stat="identity", position="fill") + 
  theme_bw()+
  scale_fill_manual(values=col.order)+
  #facet_wrap(~ loc, strip.position="top", nrow = 1)+
  facet_nested(~loc+var,
               nest_line = element_line(linetype = 1, linewidth = 0.5), 
               scales="free")+
  guides(fill=guide_legend(ncol=1))+
  labs(y= "Mean Relative Abundance")+
  theme(
    legend.direction = "vertical",
    legend.position="right",
    plot.title = element_text(size = 25, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=20),
    #axis.line.x = element_blank(),
    #axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=18),
    axis.title.x = element_blank(),
    #axis.title.x = element_blank(),
    axis.title.y =element_text(size=22),
    legend.text=element_text(size = 20),
    legend.title = element_text(size=22, face="bold"),
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    #strip.background = element_blank(),
    strip.background=element_rect(color="grey30", fill="white"),
    strip.text.x = element_text(size = 20),
    #panel.border = element_rect(colour = "black", fill = NA,linewidth= 0.2),
    panel.border=element_rect(color="grey30",fill = NA))+
  scale_y_continuous(expand = c(0,0))+
  guides(x="axis_nested")
#geom_vline(xintercept = 5, linetype="dotted", 
#color = "blue", linewidth=1.5)
norm.order.plot 

setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("RA_norm_order.tiff",
       norm.order.plot, device = "tiff",
       width = 14, height =8, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")



col.gen <- c("#F2F3F4","#222222","#F3C300","#875692","#F38400","#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856",
             "#E68FAC","#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600",
             "#654522", "#E25822", "#2B3D26","#332288")

set.seed(13)
norm.gen.plot <- ggplot(norm.gen.df, aes(x=treatment, y=Mean, fill=Genus)) + 
  geom_bar(aes(), stat="identity", position="fill") + 
  theme_bw()+
  scale_fill_manual(values=col.gen)+
  facet_nested(~loc+var,
               nest_line = element_line(linetype = 1, linewidth = 0.5), 
               scales="free")+
  guides(fill=guide_legend(ncol=1))+
  labs(y= "Mean Relative Abundance")+
  theme(
    legend.direction = "vertical",
    legend.position="right",
    plot.title = element_text(size = 25, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=20),
    #axis.line.x = element_blank(),
    #axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=18),
    axis.title.x = element_blank(),
    #axis.title.x = element_blank(),
    axis.title.y =element_text(size=22),
    legend.text=element_text(size = 20),
    legend.title = element_text(size=22, face="bold"),
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    #strip.background = element_blank(),
    strip.background=element_rect(color="grey30", fill="white"),
    strip.text.x = element_text(size = 20),
    #panel.border = element_rect(colour = "black", fill = NA,linewidth= 0.2),
    panel.border=element_rect(color="grey30",fill = NA))+
  scale_y_continuous(expand = c(0,0))+
  guides(x="axis_nested")
#geom_vline(xintercept = 5, linetype="dotted", 
#color = "blue", linewidth=1.5)
norm.gen.plot 
setwd('/Users/emiliedehon/Documents/Documents_EmilieMacAir/Project_Bean_Drought_01072024/HPCC/Biom_07/Figure')
ggsave("RA_norm_genus.tiff",
       norm.gen.plot , device = "tiff",
       width = 14, height =8, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")


#################################################################################################
#DIFFERENTIAL ABUNDANCE ANALYSIS TEST
#################################################################################################

# ANCOM - BC
BiocManager::install("ANCOMBC")
library(ANCOMBC)
library(mia)
remove.packages("ANCOMBC")
# make tse data from the not rarefied and not normalized phyloseq object

sample_data(decontaminated_unrare_notnormal_physeq)$id_com <- join.meta.data.ed$id_comb
sample_data(decontaminated_unrare_notnormal_physeq)
tse.asv = mia::makeTreeSummarizedExperimentFromPhyloseq(decontaminated_unrare_notnormal_physeq)
tse.asv 

raw.gen <- agglomerateByRank(tse.asv, "Genus")
raw.gen
raw.gen.pseq = makePhyloseqFromTreeSummarizedExperiment(raw.gen)
tax_table(raw.gen.pseq)

out.gen = ancombc2(data = tse.asv, assay_name = "counts",
                             tax_level = "Genus",
                             fix_formula = "id_com",
                             rand_formula =NULL,
                             p_adj_method = "BH", prv_cut = 0.10, lib_cut = 0,
                             group = "id_com", pairwise = T,
                             struc_zero = TRUE, neg_lb = F, 
                             iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
                             alpha = 0.05, global = TRUE,
                             n_cl = 1, verbose = TRUE)
View(out.gen$res_pair)
out.gen$res_pair$p_id_comB18504UP_id_comB18504_droughtUP
trace(ancombc2, edit = TRUE)


# DAA using glmmTMB

# 1. Subset data by location and cultivar

# 1. Lower Peninsula

# 1. B18504 - LP
B18504.LP.seq<- subset_samples(decontaminated_unrare_notnormal_physeq, loc=="LP" & var=="B18504")
B18504.LP.seq1 <- prune_taxa(taxa_sums(B18504.LP.seq)>0, B18504.LP.seq) # 299 ASVs
B18504.LP.seq1
# 2. Cayenne - LP
Cayen.LP.seq<- subset_samples(decontaminated_unrare_notnormal_physeq, loc=="LP" & var=="Cayenne")
Cayen.LP.seq1 <- prune_taxa(taxa_sums(Cayen.LP.seq)>0, Cayen.LP.seq) # 144 ASVs
Cayen.LP.seq1
# 3. R99 - LP
R99.LP.seq<- subset_samples(decontaminated_unrare_notnormal_physeq, loc=="LP" & var=="R99")
R99.LP.seq1 <- prune_taxa(taxa_sums(R99.LP.seq)>0, R99.LP.seq) # 223 ASVs
R99.LP.seq1
# 4. Rosetta - LP
Roset.LP.seq<- subset_samples(decontaminated_unrare_notnormal_physeq, loc=="LP" & var=="Rosetta")
Roset.LP.seq1 <- prune_taxa(taxa_sums(Roset.LP.seq)>0, Roset.LP.seq) # 214 ASVs
Roset.LP.seq1

###############################################################################
# Filter low-abundant taxa
# keeping ASVs with at least 0.01 % relative abundance across all samples
physeq.subset <-B18504.LP.seq1
physeq.subset 
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0001)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 266

################################################################################
#Lets generate a prevalence table (number of samples each taxa occurs in) for each taxa.
prevalencedf = apply(X = otu_table(physeq.subset),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevalencedf = data.frame(Prevalence = prevalencedf,
                          TotalAbundance = taxa_sums(physeq.subset))
prevalencedf[1:10,]
dim(prevalencedf)
# calculate prevalence
ps = physeq.subset
df_tmp <- psmelt(ps)
df_tmp$sample <- 0
df_tmp$sample[df_tmp$Abundance > 0] <- 1 #E: DON'T UNDERSTAND WHY THIS IS DONE
df_otu_prev_ttt <- data.frame(matrix(ncol=nlevels(as.factor(df_tmp$id_com)),
                                     nrow=nlevels(as.factor(df_tmp$OTU)), 
                                     dimnames=list(levels(as.factor(df_tmp$OTU)),
                                                   levels(as.factor(df_tmp$id_com)))))
view(df_otu_prev_ttt)
#attention il ya Sample et sample
for (i in unique(df_tmp$OTU)) {
  for (j in unique(df_tmp$id_com)) {
    df_otu_prev_ttt[i,j] <- sum(df_tmp$sample[df_tmp$OTU == i & df_tmp$id_com == j],na.rm = T) / nrow(df_tmp[df_tmp$OTU == i & df_tmp$id_com == j,]) *100
    print(paste(i,j,df_otu_prev_ttt[i,j]),sep="\t")
    #print(df_otu_prev_ttt[i,j])
  }
  
}

df_otu_prev_ttt$max_prev <- apply(df_otu_prev_ttt,MARGIN=1, FUN=max)
view(df_otu_prev_ttt)

# filter otu par prevalence
otu_table(physeq.subset)
ps =  physeq.subset 
df_prev = df_otu_prev_ttt
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 50,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75  # B18504.LP.seq1=54 taxa
otu_table(physeq.subset.75)

# Run the model
#install.packages("glmmTMB", type="source")
library(glmmTMB)
library(lme4)
#.rs.restartR()

otu_table(physeq.subset.75)
sample_data(physeq.subset.75)
tmp_T3s <- physeq.subset.75
str(tmp_T3s)
#  treatment
a = tibble("sample"= tmp_T3s@sam_data$SampleID,
           "treatment"= as.character(tmp_T3s@sam_data$treatment))
# force control as intercept
a[a == "open"] <- "1a"
a = as.factor(a$treatment)
# offset
o = log(sample_sums(physeq.subset.75)) # using filtered data
#Since library sizes are different between groups, accounting for library size results in the gene no longer being DE at the 5% significance level. Not correcting for sequencing depth would thus result in spurious results.
# random effect
z <- as.factor(tmp_T3s@sam_data$SampleID)

# model with pairwise comparison ---------------------------------------------------------------------------------
glmT3s.sum.global = data.frame()
glmT3s.pairwise.global = data.frame()
#fam=genpois(link = "log")

for (i in 1:length(taxa_names(tmp_T3s))) {
  
  OTU = taxa_names(tmp_T3s)[i] 
  
  # response variable
  y = as.vector(tmp_T3s@otu_table[OTU,]@.Data)
  
  tryCatch({
    ### model
    glmT3s <- glmmTMB(y ~ a + (1|z) , family="poisson", offset = o)
    #glmT3s <- glm(y ~ a, family='poisson')
    glmT3s.sum = summary(glmT3s)$coefficients$cond
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum))
    glmT3s.sum
    glmT3s.sum.global = rbind(glmT3s.sum.global,glmT3s.sum)
    ### multiple comparison
    glmT3s.pairwise = emmeans(glmT3s,pairwise~a)
    # select p value
    glmT3s.pairwise.sum = summary(glmT3s.pairwise)
    glmT3s.pairwise.sum = glmT3s.pairwise.sum[["contrasts"]]
    # extract summary
    tmp_df = glmT3s.pairwise.sum
    # keep only comparisons of interest
    tmp = unlist(strsplit(as.character(tmp_df$contrast)," - "))
    tmp_df[,"a"] <- tmp[seq(1,length(tmp),by=2)]
    tmp_df[,"b"] <- tmp[seq(2,length(tmp),by=2)]
    #tmp_df = tmp_df[grep("Ni",tmp_df$b), ]
    tmp_df = cbind("OTU"=OTU,tmp_df)
    # extract results in data frame
    glmT3s.pairwise.global = rbind(glmT3s.pairwise.global,tmp_df)
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  rm(OTU,y,glmT3s,glmT3s.sum)
}

glmT3s.model.global = glmT3s.sum.global
glmT3s.pairwise.global = glmT3s.pairwise.global
#When there is only one thing to test, there is no multiplicity issue, and hence no multiplicity adjustment to the P values.
glmT3s.pairwise.global$p.adjust <- p.adjust(glmT3s.pairwise.global$p.value, method = "fdr")
view(glmT3s.pairwise.global)






