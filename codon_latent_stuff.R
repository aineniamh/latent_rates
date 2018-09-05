library(ggplot2)
library(grid)
library(gridExtra)


setwd("~/Google Drive/phd/latent_rates/")
cbbPalette <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#CC79A7")
d<-  read.csv('spliced_cds_work/identify_and_classify_snps_results.csv', header = T, sep=',')

data<-  read.csv('clustering_info_EBOV_2017.csv', header = T, sep=',')
data$'Y'<- 'A'
f <- ggplot(data, aes(y=Cluster,x=Position))+ 
  geom_vline(aes(xintercept=Position), colour="grey")+
  geom_point(aes(color=Cluster), size =4)+
  scale_colour_manual(values=cbbPalette)+
  theme_classic() +
  theme(text = element_text(size=15)) 
f

f <- ggplot(data, aes(x=Cluster))+ 
  geom_bar()+
  theme_classic() +
  theme(text = element_text(size=15)) 
f

d$'Y'<- 1





f <- ggplot(d, aes(x=Type))+ 
  geom_bar(aes(fill=Type))+
  labs(x = "SNP Position",y="") +
  theme_classic() +
  theme(text = element_text(size=15), legend.position="none") 
f

f <- ggplot(d, aes(y=Type,x=Position, group = Type, color=as.factor(Position_in_Codon)))+ 
  geom_jitter( width = 0, height=0.1)+
  theme_classic() +
  theme(text = element_text(size=15), legend.position="none") 
f

h<- ggplot(d, aes(x=SNP)) + 
  geom_bar()+
  theme_classic()
h

m<-d[which(d$Type=='not_monophyly'), ]
m$'Y2'<- 2
m$'Y3'<-1.5
f <- ggplot(m, aes(y=Y,x=Position, color=as.factor(Position_in_Codon)))+ 
  geom_point()+
  theme_classic() +
  geom_point(aes(x=Position, y=Y2,color=Gene))+
  geom_point(aes(x=Position, y=Y3,color=SNP))+
  theme(text = element_text(size=15), legend.position="none") 
f


dl <- d[which(d$N_diffs<200), ]
dt <- dl[which(dl$Type_diff=='transition'|dl$Type_diff=='transversion'), ]
f <- ggplot(dt, aes(x=log10(dt$N_diffs)))+ 
  geom_histogram(aes(fill= dt$Type_diff))+
  labs(x = "",y="Number of differences") +
  theme_classic() +
  theme(text = element_text(size=15), legend.position="none") +
  ggtitle('Types of differences')
f
f <- ggplot(dt, aes(x=reorder(dt$Difference,
                              dt$Difference,
                              function(x)-length(x))))+ 
  geom_bar(aes(fill= dt$Type_diff))+
  labs(x = "",y="Number of differences") +
  theme_classic() +
  coord_flip() +
  theme(text = element_text(size=15), legend.position="top", legend.title = element_blank()) +
  ggtitle('Types of differences')
f
x=reorder(Position,Position,
          function(x)-length(x))


