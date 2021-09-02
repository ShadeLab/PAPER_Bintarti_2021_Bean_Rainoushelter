# PF read data by well position script
# Abby Sulesky-Grieb


### make csv file of sample name and PF reads column with well position/sample_ID,
### can't have commas in numbers - use data type in excel to switch to "number"
### read in csv file
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210601_16SV4')
wd <- print(getwd())
seq_data_1 <- read.csv("20210601_PF_read_template.csv")


### plot PF_reads by well position (A1, B1, C1 - H1, A2, B2, etc.)

library(ggplot2)
ggplot(data=seq_data_1, aes(x=sample_ID, y=PF_reads)) + geom_point()

### subset a new dataframe of the negative/positive controls

library(dplyr)

# neg controls first
sample_ID <- c("72","73","74","75","76","77","78", "80")
sample_ID <- as.numeric(sample_ID)
neg_controls <- data.frame(sample_ID)
neg_controls<- inner_join(neg_controls, seq_data_1)

# resave sample_ID data frame as positive controls
sample_ID <- c("65","66","67","68","69","70","71","79")
sample_ID <- as.numeric(sample_ID)
pos_controls <- data.frame(sample_ID)
pos_controls <- inner_join(pos_controls, seq_data_1)

#seq_data_1$sample_ID <- as.factor(seq_data_1$sample_ID)
str(seq_data_1)
#seq_data_1$PF_reads <- as.numeric(seq_data_1$PF_reads)
pf.reads.unfil <- (ggplot(data=seq_data_1, aes(x=sample_ID, y=PF_reads)) + 
        geom_point()+
geom_point(data=neg_controls, 
           aes(x=sample_ID,y=PF_reads), 
           color='red',
           size=2) +
geom_point(data=pos_controls, 
             aes(x=sample_ID,y=PF_reads), 
             color='blue',
             size=2))
 
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_PFReads_unfiltered.eps",
      pf.reads.unfil, device = "eps",
       width = 8, height =5, 
       units= "in", dpi = 600)



## arrange data frame by PF_reads

seq_data_ordered <- seq_data_1 %>% dplyr::arrange(PF_reads)
order_ID <- as.numeric(c(1:80))
seq_data_ordered <- cbind(seq_data_ordered, order_ID)

### subset new dataframes of the negative/positive controls
# neg controls first
order_ID <- c("1","2","3","4","5","6","7","8")
order_ID <- as.numeric(order_ID)
neg_controls_ordered <- data.frame(order_ID)
neg_controls_ordered <- inner_join(neg_controls_ordered, seq_data_ordered)

# pos controls
order_ID <- c("20","21","29","36","46","64","69","73")
order_ID <- as.numeric(order_ID)
pos_controls_ordered <- data.frame(order_ID)
pos_controls_ordered <- inner_join(pos_controls_ordered, seq_data_ordered)


pf.reads.unfil.ordered <- (ggplot(data=seq_data_ordered, aes(x=order_ID, y=PF_reads)) + geom_point() + 
  geom_point(data=neg_controls_ordered, 
             aes(x=order_ID,y=PF_reads), 
             color='red',
             size=2) +
  geom_point(data=pos_controls_ordered, 
             aes(x=order_ID,y=PF_reads), 
             color='blue',
             size=2))
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210601_16SV4')
ggsave("20210601_PFReads_ordered_unfiltered.eps",
      pf.reads.unfil.ordered, device = "eps",
       width = 8, height =5, 
       units= "in", dpi = 600)

##################################################################################################################################
##################################################################################################################################

### read in csv file
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
wd <- print(getwd())
seq_data_2 <- read.csv("20210604_PF_read_template.csv")


### plot PF_reads by well position (A1, B1, C1 - H1, A2, B2, etc.)

library(ggplot2)
ggplot(data=seq_data_2, aes(x=sample_ID, y=PF_reads)) + geom_point()

### subset a new dataframe of the negative/positive controls

library(dplyr)

# neg controls first
sample_ID <- c("72","73","74","75","76","77","78", "80")
sample_ID <- as.numeric(sample_ID)
neg_controls <- data.frame(sample_ID)
neg_controls<- inner_join(neg_controls, seq_data_2)

# resave sample_ID data frame as positive controls
sample_ID <- c("65","66","67","68","69","70","71","79")
sample_ID <- as.numeric(sample_ID)
pos_controls <- data.frame(sample_ID)
pos_controls <- inner_join(pos_controls, seq_data_2)

#seq_data_1$sample_ID <- as.factor(seq_data_1$sample_ID)
str(seq_data_2)
#seq_data_1$PF_reads <- as.numeric(seq_data_1$PF_reads)
pf.reads.unfil.2 <- (ggplot(data=seq_data_2, aes(x=sample_ID, y=PF_reads)) + 
        geom_point()+
geom_point(data=neg_controls, 
           aes(x=sample_ID,y=PF_reads), 
           color='red',
           size=2) +
geom_point(data=pos_controls, 
             aes(x=sample_ID,y=PF_reads), 
             color='blue',
             size=2))
pf.reads.unfil.2
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_PFReads_unfiltered.eps",
      pf.reads.unfil.2, device = "eps",
       width = 8, height =5, 
       units= "in", dpi = 600)

### arrange data frame by PF_reads

seq_data_ordered2 <- seq_data_2 %>% dplyr::arrange(PF_reads)
order_ID <- as.numeric(c(1:80))
seq_data_ordered2 <- cbind(seq_data_ordered2, order_ID)

### subset new dataframes of the negative/positive controls
# neg controls first
order_ID <- c("1","2","3","4","5","6","7","8")
order_ID <- as.numeric(order_ID)
neg_controls_ordered2 <- data.frame(order_ID)
neg_controls_ordered2 <- inner_join(neg_controls_ordered2, seq_data_ordered2)

# pos controls
order_ID <- c("20","22","30","34","46","65","68","72")
order_ID <- as.numeric(order_ID)
pos_controls_ordered2 <- data.frame(order_ID)
pos_controls_ordered2 <- inner_join(pos_controls_ordered2, seq_data_ordered2)


pf.reads.unfil.ordered2 <- (ggplot(data=seq_data_ordered2, aes(x=order_ID, y=PF_reads)) + geom_point() + 
  geom_point(data=neg_controls_ordered2, 
             aes(x=order_ID,y=PF_reads), 
             color='red',
             size=2) +
  geom_point(data=pos_controls_ordered2, 
             aes(x=order_ID,y=PF_reads), 
             color='blue',
             size=2))
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_PFReads_ordered_unfiltered.eps",
      pf.reads.unfil.ordered2, device = "eps",
       width = 8, height =5, 
       units= "in", dpi = 600)

## 1. filtered reads

### read in csv file
setwd('/Users/arifinabintarti/Documents/PAPER/PAPER_Bintarti_2021_Bean_Rainoutshelter/16SV4_OTU97/20210604_16SV4')
wd <- print(getwd())
seq_data_fil2 <- read.csv("20210604_PF_read_template_filtered.csv")


### plot PF_reads by well position (A1, B1, C1 - H1, A2, B2, etc.)

library(ggplot2)
ggplot(data=seq_data_fil2, aes(x=sample_ID, y=PF_reads)) + geom_point()

### subset a new dataframe of the negative/positive controls

library(dplyr)

# neg controls first
sample_ID <- c("72","73","74","75","76","77","78", "80")
sample_ID <- as.numeric(sample_ID)
neg_controls <- data.frame(sample_ID)
neg_controls<- inner_join(neg_controls, seq_data_fil2)

# resave sample_ID data frame as positive controls
sample_ID <- c("65","66","67","68","69","70","71","79")
sample_ID <- as.numeric(sample_ID)
pos_controls <- data.frame(sample_ID)
pos_controls <- inner_join(pos_controls, seq_data_fil2)

#seq_data_1$sample_ID <- as.factor(seq_data_1$sample_ID)
str(seq_data_fil2)
#seq_data_1$PF_reads <- as.numeric(seq_data_1$PF_reads)
pf.reads.fil.2 <- (ggplot(data=seq_data_fil2, aes(x=sample_ID, y=PF_reads)) + 
        geom_point()+
geom_point(data=neg_controls, 
           aes(x=sample_ID,y=PF_reads), 
           color='red',
           size=2) +
geom_point(data=pos_controls, 
             aes(x=sample_ID,y=PF_reads), 
             color='blue',
             size=2))
pf.reads.fil.2
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_PFReads_fil.eps",
      pf.reads.fil.2, device = "eps",
       width = 8, height =5, 
       units= "in", dpi = 600)

### arrange data frame by PF_reads

seq_data_ordered_fil2 <- seq_data_fil2 %>% dplyr::arrange(PF_reads)
order_ID <- as.numeric(c(1:80))
seq_data_ordered_fil2 <- cbind(seq_data_ordered_fil2, order_ID)

### subset new dataframes of the negative/positive controls
# neg controls first
order_ID <- c("1","59","60","67","68","70","71","72")
order_ID <- as.numeric(order_ID)
neg_controls_ordered_fil2 <- data.frame(order_ID)
neg_controls_ordered_fil2 <- inner_join(neg_controls_ordered_fil2, seq_data_ordered_fil2)

# pos controls
order_ID <- c("73","74","75","76","77","78","79","80")
order_ID <- as.numeric(order_ID)
pos_controls_ordered_fil2 <- data.frame(order_ID)
pos_controls_ordered_fil2 <- inner_join(pos_controls_ordered_fil2, seq_data_ordered_fil2)


pf.reads.fil.ordered2 <- (ggplot(data=seq_data_ordered_fil2, aes(x=order_ID, y=PF_reads)) + geom_point() + 
  geom_point(data=neg_controls_ordered_fil2, 
             aes(x=order_ID,y=PF_reads), 
             color='red',
             size=2) +
  geom_point(data=pos_controls_ordered_fil2, 
             aes(x=order_ID,y=PF_reads), 
             color='blue',
             size=2))
pf.reads.fil.ordered2
setwd('/Users/arifinabintarti/Documents/Research/Seeds_microbiome/Rainoutshelter/16SV4_OTU97/20210604_16SV4')
ggsave("20210604_PFReads_ordered_filtered.eps",
      pf.reads.fil.ordered2, device = "eps",
       width = 8, height =5, 
       units= "in", dpi = 600)

