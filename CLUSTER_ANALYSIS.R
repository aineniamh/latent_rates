install.packages("ggplot2")
install.packages("ggExtra")
install.packages("grid")
install.packages("gridExtra")
install.packages("plotly")
install.packages("ggrepel")
install.packages("scales")
install.packages("viridis")
install.packages("treeio")
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
library(ggtree)
library(ggplot2)
library(ggExtra)
library(grid)
library(treeio)
library(viridis)
library(gridExtra)
library("ggrepel") # for spreading text labels on the plot, you can replace with `geom_text` if you want
library("scales")
install.packages("gtable")
library(gtable)
library(reshape2)
setwd("~/Google Drive/phd/latent_rates/")


d<-  read.csv('spliced_cds_work/identify_and_classify_snps_results.csv', header = T, sep=',')
d<-  read.csv('wout17/identify_and_classify_snps_results.csv', header = T, sep=',')
d<-  read.csv('wout17/cluster_vs_cluster_mutations.csv', header = T, sep=',')


g<- read.csv('gene_start_stop.csv',header = T)
g$y1<-0
g$y2<-1

d$Position_in_Codon<- as.factor(d$Position_in_Codon)
d_nona<- subset(d, Cluster!='NA')
d_nona<- subset(d_nona, Cluster!='M')
d_nona<- subset(d_nona, Cluster!='A:B:C:D:E:M')
d_nona<- subset(d_nona, Cluster!='A:B:C:D:M')
d_nona<- subset(d_nona, Cluster!='A:B:C:E:M')
d_nona<- subset(d_nona, Cluster!='A:B:D:E:M')
d_nona<- subset(d_nona, Cluster!='A:B:D:M')
d_nona<- subset(d_nona, Cluster!='A:B:E:M')
d_nona<- subset(d_nona, Cluster!='A:D:E:M')
d_nona<- subset(d_nona, Cluster!='A:C:D:E:M')
d_nona<- subset(d_nona, Cluster!='A:C:D:M')
d_nona<- subset(d_nona, Cluster!='A:C:E:M')
d_nona<- subset(d_nona, Cluster!='B:C:D:E:M')
d_nona<- subset(d_nona, Cluster!='B:C:E:M')
d_nona<- subset(d_nona, Cluster!='B:D:M')

not_with_tree <- subset(d, Type=='not_monophyly')

p1<- ggplot(data=g)+
  geom_rect(data=g, aes(xmin=Start,xmax=Stop,ymin=y1,ymax=y2),fill="white", alpha=0.5)+
  theme_void()+
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0),"cm"))+
  annotate("text", x=(g$Start+g$Stop)/2, y=0.2, label=g$Gene, size = 4)
p1
p2 <- ggplot(not_with_tree, aes(y=Cluster,x=Position))+ 
  #geom_vline(data=g,xintercept=g$Stop, linetype="dashed", colour="grey")+
  geom_rect(data=not_with_tree, aes(xmin=2221, xmax=3244), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_rect(data=not_with_tree, aes(xmin=4225, xmax=6256), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_rect(data=not_with_tree, aes(xmin=7123, xmax=7879), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_hline(color="grey", size =0.5,yintercept = c(5,10,15,20,25))+
  geom_text(data=g, inherit.aes = FALSE, aes(x=(g$Start+g$Stop)/2,y=Inf,label=g$Gene),
            vjust=1.2, size=4)+
  geom_point(size =2)+
  theme_classic()+
  ggtitle("EBOV SNPs by Mutation")+
  scale_y_discrete(expand = c(0.05,0))+
  theme(legend.position = "bottom", legend.title = element_blank(),plot.title = element_text(hjust=0.5))+
  theme(text = element_text(size=15))
  #scale_x_continuous(limits = c(0, 13584))+
  #scale_x_continuous(breaks = g$Stop)+
  #scale_color_viridis(discrete=TRUE)
  
            #vjust=1.2, size=4)+
  #scale_colour_manual(values=c("#3F226B","#ADEA75"))
p2

p3 <- ggplot(not_with_tree, aes(x=SNP))+
  theme_classic()+
  xlab("Cases of Homoplasy")+
  ylab("Count")+
  theme(text = element_text(size=35))+
  geom_histogram(stat = "count", fill="#D55E00")
p3


g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

gA=ggplot_gtable(ggplot_build(p1))
gB=ggplot_gtable(ggplot_build(p2))
maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
grid.newpage()

gl<-list(gA, gB, tree)
grid.arrange(
  grobs = gl,
  widths = c(1, 1),
  heights = c(.1,.9),
  layout_matrix = rbind(c(1, NA),
                        c(2, 3))
)

nwk <- read.newick(file="wout17/RAxML_bipartitions.raxml_rapid_spliced_cds_tree")
nwk_nodes <- read.newick(file="wout17/RAxML_bipartitions.raxml_rapid_spliced_cds_tree.node_labelled.nw")
nwk$node.label<-nwk_nodes$node.label
annotations<- read.csv("wout17/clustering_annotations_just_the_tips.txt", header = TRUE, sep="\t")
annotations$Year<-as.factor(annotations$Year)
cbPalette <- c("#D55E00","#56B4E9","#56B4E9", "#E69F00", "#009E73","#56B4E9", "#CC79A7")

tree<- ggtree(tr = nwk) +
  geom_tree(color="#999999")
treea<- tree %<+% annotations + 
  geom_tippoint(aes(color=Outbreak),size=4)+
  scale_color_manual(values=cbPalette)+
  theme(legend.position = "right",
        legend.title = element_blank(),
        text = element_text(size=20),
        legend.text = element_text(size=18, colour = "#999999"))
treea

tree<- ggtree(tr = nwk) +
  #geom_text2(aes(subset=!isTip, label=label), hjust=-.3)+
  #geom_text2(aes(subset=!isTip, label=node), vjust=-.3)+
  geom_hilight(node=102, fill="#CC6666",alpha=0.7,extend = 0.001)  +
  geom_hilight(node=56, fill="#83CC66",alpha=0.7,extend = 0.001) + 
  geom_hilight(node=60, fill="#66CCA0",alpha=0.7,extend = 0.001) + 
  geom_hilight(node=61, fill="#CCBD66",alpha=0.7,extend = 0.001)+
  geom_hilight(node=77, fill="#66A0CC",alpha=0.7,extend = 0.001)  + 
  geom_hilight(node=69, fill="#8366CC",alpha=0.7,extend = 0.001)  + 
  geom_hilight(node=67, fill="#CC79A7",alpha=0.7,extend = 0.001)  +
  geom_cladelabel(node=102, label="A", align = TRUE, offset=0.002) +
  geom_cladelabel(node=56, label="B2", align = FALSE) +
  geom_cladelabel(node=55, label="B", align = TRUE, offset=0.002) +
  geom_cladelabel(node=60, label="B3", align = FALSE) +
  geom_cladelabel(node=61, label="B1", align = FALSE) +
  geom_cladelabel(node=77, label="C", align = TRUE, offset=0.002) +
  geom_cladelabel(node=69, label="D", align = TRUE, offset=0.002)+
  geom_cladelabel(node=67, label="E", align = TRUE, offset=0.002) +
  geom_tree()
  #geom_text(aes(label=label))
  #geom_text2(aes(label="M", subset=isTip & f), nudge_x = 0.0005, size=7)#+
  #geom_text2(aes(label="2014", subset=isTip & f), nudge_x = 0.0005, size=7)
tree<-flip(tree, 68, 67)
tree<-flip(tree, 77, 69)
tree

metadata <- read.csv("wout17/pos_uracil.txt",header = TRUE,sep='\t')
metadata<-subset(metadata, Country!='?')
mcc <- read.beast(file="wout17/1610_makona.mcc.tree")
tree<- ggtree(tr = mcc)+
  geom_tree(colour="grey")

tree1<- tree %<+% metadata + 
  geom_tippoint(aes(color=One),size=1.2)+
  scale_color_continuous(low="darkblue",high="orange", breaks= c(29,29.5,30))+
  #scale_x_continuous(breaks = c(29,29.5,30))+
  ggtitle("Codon Position One")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        #plot.title = element_text(size=20),
        text = element_text(size=20),
        legend.text = element_text(size=12))
tree1

tree2<- tree %<+% metadata + 
  geom_tippoint(aes(color=Two),size=1.2)+
  scale_color_continuous(low="darkblue",high="orange")+
  ggtitle("Codon Position Two")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text(size=20),
        legend.text = element_text(size=12))
tree2

tree3<- tree %<+% metadata + 
  geom_tippoint(aes(color=Three),size=1.2)+
  ggtitle("Codon Position Three")+
  scale_color_continuous(low="darkblue",high="orange", 
                         breaks=c(31.5,32,32.5,33), 
                         labels=c("31.5","32","32.5","33"))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text(size=20),
        legend.text = element_text(size=12))
tree3

tree4<- tree %<+% metadata + 
  geom_tippoint(aes(color=Intergenic),size=1.2)+
  scale_color_continuous(low="darkblue",high="orange")+
  ggtitle("Intergenic")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text(size=20),
        legend.text = element_text(size=12))
tree4

grid.arrange(tree1, tree2, tree3, tree4, ncol=2)

full<-read.csv('full_snp_record.csv', header = T, sep=',')
full
melted_full <- melt(full, id.vars=c("Position"))
melted_full_yes<- subset(melted_full, value=="Yes")
melted_full_no<- subset(melted_full, value=="No")
colnames(flipped_full) <- full$SNP
flipped_full<-flipped_full[2:78,]

p2 <- ggplot(melted_full_no, aes(y=variable,x=Position, color=variable))+ 
  geom_vline(data=g,xintercept=g$Stop, linetype="dashed", colour="grey")+
  #geom_line(color="black", size =0.5)+
  geom_point(size =0.1)+
  theme_classic()+
  theme(legend.position = "none", legend.title = element_blank())+
  theme(text = element_text(size=15),
        axis.text.y=element_blank())
  
  #scale_x_continuous(breaks = g$Stop)+
  #scale_color_viridis(discrete=TRUE)
p2
############


#grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
#ggExtra::ggMarginal(f, type = "violin", margins = "x", size = 4, marginCol = "red")
adar<-  read.csv('wout17/adar_signals_tips.csv', header = T, sep=',')
cB2 <-c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7",
        "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
        "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
        "#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7",
        "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
tc<-  read.csv('wout17/t_to_c_mutations_tips.csv', header = T, sep=',')
tc$Date<- as.Date(tc$Date, "%Y-%m-%d")
adar$Date<- as.Date(adar$Date, "%Y-%m-%d")
#tc<-tc[tc$Sequence %in% adar$Sequence, ]
#g$position<-(g$Start+g$Stop)/2
#tc$Sequence2 <- reorder(tc$Sequence, tc$Date)
#adar$Sequence2 <- reorder(adar$Sequence, adar$Date)
adar$Signal<-as.factor(adar$Signal)
tc[order(tc$Date, decreasing = TRUE),]
clus<- unique(tc$Cluster)
adar1 <- ggplot(tc, aes(y=reorder(Sequence, tc$Date),x=Position))+ 
  geom_rect(data=tc, aes(xmin=2221, xmax=3244), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_rect(data=tc, aes(xmin=4225, xmax=6256), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_rect(data=tc, aes(xmin=7123, xmax=7879), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_point(size =1, color="#e3d8f3")+
  geom_point(data=adar, aes(y=reorder(Sequence, adar$Date), x=Position),
              size =1, shape=3,color="darkblue")+
  theme_classic()+
  xlab("T -> C Mutations")+
  ylab("")+
  #scale_y_discrete(limits=reorder(tc$Sequence, tc$Date))+
  #ggtitle("ADAR Signals")+
  scale_y_discrete(expand = c(0.05,0), labels=element_blank())+
  geom_text(data=g,aes(x=(g$Start+g$Stop)/2,y=Inf,label=Gene),
            vjust=1.2, size=5.5)+
  theme(text = element_text(size=20),
        plot.title = element_text(size = 20, hjust=0.5),
        axis.ticks.y = element_blank(),
        legend.position = "none")
adar1

tcsub<-subset(tc, grepl("2014",tc$Date))
asub<-subset(adar, grepl("2014",adar$Date))
tcsub<-subset(tcsub, grepl("Lomela",tcsub$Sequence)==FALSE)
asub<-subset(asub, grepl("Lomela",asub$Sequence)==FALSE)
adar2 <- ggplot(tcsub, aes(y=reorder(Sequence, tcsub$Date),x=Position))+ 
  geom_rect(data=tcsub, aes(xmin=2221, xmax=3244), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_rect(data=tcsub, aes(xmin=4225, xmax=6256), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_rect(data=tcsub, aes(xmin=7123, xmax=7879), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_point(size =1, color="#e3d8f3")+
  geom_point(data=asub, aes(y=reorder(Sequence, asub$Date), x=Position),
             size =1, shape=3,color="darkblue")+
  theme_classic()+
  xlab("T -> C Mutations")+
  ylab("")+
  #scale_y_discrete(limits=reorder(tc$Sequence, tc$Date))+
  #ggtitle("ADAR Signals")+
  scale_y_discrete(expand = c(0.05,0))+
  geom_text(data=g,aes(x=(g$Start+g$Stop)/2,y=Inf,label=Gene),
            vjust=1.2, size=4)+
  theme(text = element_text(size=12),
        plot.title = element_text(size = 20, hjust=0.5),
        legend.position = "none")
adar2


pa <- ggplot_gtable(ggplot_build(p1))
p2 <- ggplot_gtable(ggplot_build(adar1))
maxWidth = unit.pmax(pa$widths[2:3], p2$widths[2:3])
pa$widths[2:3] <- maxWidth
p2$widths[2:3] <- maxWidth
grid.arrange(pa,p2,heights = c(.05,.9))


tree2<- ggtree(tr = nwk, layout="unrooted") +
  geom_hilight(node=43, fill="#A7DBC0",alpha=0.7,extend = 0.001) +
  geom_cladelabel(node=53, label="2014\nW.Af", align = TRUE, offset=0.002) +
  geom_cladelabel(node=43, label="2003\nDRC", align = TRUE, offset=0.002)+
  geom_hilight(node=53, fill="#A7DBC0",alpha=0.7,extend = 0.001) 
tree2<-flip(tree2, 43, 44)
tree2<-flip(tree2, 53, 45)
tree2

pa <- ggplot_gtable(ggplot_build(p1))
p2 <- ggplot_gtable(ggplot_build(adar1))
maxWidth = unit.pmax(pa$widths[2:3], p2$widths[2:3])
pa$widths[2:3] <- maxWidth
p2$widths[2:3] <- maxWidth
grid.arrange(pa,p2,heights = c(.1,.9))

#########
ordering <- c('A','B1','B2','B3','C', 'D','E',
              'B1:B3','B1:B2:B3', 'B1|B1:B2:B3',  
              'C:D','C|D','C:D|D', 'C|C:D',
              'C:D:E', 'C|C:D:E','C:D:E|D', 'C|E',
              'A|C:D', 'A|D', 'A|E', 'B1:B2:B3|C', 'B1:B2:B3|D', 'B1|C', 
              'B1|D', 'B2|C', 'B2|D', 'B2|E', 'B3|C', 'B3|D', 'B3|E')

ordering<- c('A', 'B1','B2','B3','C', 'D', 'E',
             'B1:B3','B1:B2:B3','C:D','C:D:E',
            #'C|C:D','C:D|D',
             'C:D:E|D','C|C:D:E',
             'A|C:D', 'A|D', 'A|E',  
             'B1|B1:B2:B3','B1|C', 'B1|D', 'B1:B2:B3|C', 'B1:B2:B3|D',    
             'B2|C', 'B2|D', 'B2|E',  'B3|C', 'B3|D', 'B3|E', 'C|E')
branch<-read.csv('wout17/cluster_vs_cluster_mutations.csv', header = T, sep=',')
p2 <- ggplot(d, aes(x=Position))+
  theme_classic()+
  geom_histogram(stat = "count",fill="#D55E00")+
  xlab("Cluster")+
  ylab("Count")
p2

p2 <- ggplot(branch, aes(y=Clusters,x=Position))+ 
  #geom_vline(data=g,xintercept=g$Stop, linetype="dashed", colour="grey")+
  geom_rect(data=branch, aes(xmin=2221, xmax=3244), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_rect(data=branch, aes(xmin=4225, xmax=6256), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_rect(data=branch, aes(xmin=7123, xmax=7879), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_hline(color="grey", size =0.5,yintercept = c(17.5), linetype="dashed")+
  geom_text(data=g, inherit.aes = FALSE, aes(x=(g$Start+g$Stop)/2,y=Inf,label=g$Gene),
            vjust=1.2, size=4)+
  geom_point(size =2,color="darkblue")+
  theme_classic()+
  ggtitle("EBOV Substitutions by Clusters")+
  scale_y_discrete(limits = rev(ordering),expand = c(0.05,0))+
  theme(legend.position = "bottom", legend.title = element_blank(),plot.title = element_text(hjust=0.5))+
  theme(text = element_text(size=15))

p2


write.table(wide_u, file = "wout17/wide_uracil.txt", sep='\t',row.names=FALSE)

b<-ggplot(branch, aes(x=Position, y=Clusters))+
  geom_rect(data=branch, aes(xmin=2221, xmax=3244), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_rect(data=branch, aes(xmin=4225, xmax=6256), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_rect(data=branch, aes(xmin=7123, xmax=7879), ymin=-Inf, ymax=Inf,fill="#e9e9e9",color="#e9e9e9", alpha=0.5) +
  geom_point(color="D55E00")+
  ggtitle("Shared substitutions")+
  theme_classic()+
  geom_text(data=g,aes(x=(g$Start+g$Stop)/2,y=Inf,label=Gene),
            vjust=1.2, size=4)+
  theme(text = element_text(size=13),
        plot.title = element_text(size=20,hjust=0.5))+
  scale_y_discrete(limits = rev(ordering),expand = c(0.05,0))
b

cluster<-  read.csv('wout17/branch_mutations_counts.csv', header = T, sep=',')
wide_cluster<- dcast(cluster, Cluster1~Cluster2)
row.names(wide_cluster)<-wide_cluster$Cluster1
wide_cluster<-wide_cluster[,-1]
get_upper_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
upper_tri <- get_upper_tri(wide_cluster)
upper_tri$Cluster1<- row.names(upper_tri)

melted_tri<-melt(upper_tri, na.rm=TRUE)
melted_tri$logged<-log10(melted_tri$value+1)
sub<- subset(melted_tri, Cluster1==variable)

p <- ggplot(melted_tri, aes(Cluster1,variable)) + 
  geom_tile(aes(fill=melted_tri$logged), colour="white")+
  theme_classic()+
  scale_fill_gradient(low = 'white',high = 'darkblue') +
  theme(axis.ticks=element_blank(),
        axis.line=element_blank(),
        plot.title = element_text(size=20,hjust=0.5),
        text = element_text(size=15),
        legend.position = "none")+
  labs(x='',y='',title='Shared Substitutions')+
  geom_text(aes(label=value), size= 5)+
  geom_text(data=sub, aes(label=value), color="white", size=5)+
  scale_y_discrete(limits = rev(levels(cluster$Cluster2)))+
  scale_x_discrete(position = "top")
p
cB <- c("#D55E00","#E69F00","#CC79A7","#999999")
###########

base_content<-  read.csv('wout17/1610_base_content.csv', header = T, sep=',')
base_content$Date<-as.Date(base_content$Date, '%Y-%m-%d')
base_content<-subset(base_content, base_content$Base %in% c('A','T','C','G'))
base_content$ordered = factor(base_content$Position, levels=c('One','Two','Three','Intergenic'))
baseinfo<-ggplot(base_content, aes(x=Date, y=Percentage,  color=Base))+
  geom_point(size=1)+
  geom_smooth(color="black", size=0.5, method = "loess")+
  #ggtitle("Position 1")+
  theme_classic()+
  theme(text =element_text(size=10), 
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 20),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "3 month")
  #scale_colour_gradient(low="darkblue",high="orange")
baseinfo2<-baseinfo + 
  facet_wrap(~interaction(ordered,Base),scales = "free_y") +
  theme(strip.background = element_blank())+
        #strip.text.x = c("One","Two","Three","Intersection"))+
  scale_color_viridis(discrete = TRUE)
baseinfo2

base_content<-  read.csv('spliced_cds_work/39_ref_base_content.csv', header = T, sep=',')
base_content$Date<-as.Date(base_content$Date, '%Y-%m-%d')
base_content<-subset(base_content, base_content$Base %in% c('A','T','C','G'))
base_content$ordered = factor(base_content$Position, levels=c('One','Two','Three'))
baseinfo<-ggplot(base_content, aes(x=Date, y=Percentage,  color=Base))+
  geom_point(size=1)+
  geom_smooth(color="black", size=0.5, method = "loess")+
  #ggtitle("Position 1")+
  theme_classic()+
  theme(text =element_text(size=10), 
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 20),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "3 years")
baseinfo
#scale_colour_gradient(low="darkblue",high="orange")
baseinfo2<-baseinfo + 
  facet_grid(~interaction(ordered,Base),scales = "free_y") +
  theme(strip.background = element_blank())+
  #strip.text.x = c("One","Two","Three","Intersection"))+
  scale_color_viridis(discrete = TRUE)
baseinfo2

#########
u<-  read.csv('wout17/1610_uracil_content.csv', header = T, sep=',')
u$Date<-as.Date(u$Date, '%Y-%m-%d')
u$log_pcent<-u$Percentage
u$Position <- factor(u$Position, levels = c("One","Two",'Three',"Intergenic"))
wide_u<- dcast(u, Sequence~Position)

u1<-subset(u, Position=="One")
u2<-subset(u, Position=="Two")
u3<-subset(u, Position=="Three")
u4<-subset(u, Position=="Intergenic")

u_w_info <- read.csv('wout17/pos_uracil.txt',sep = '\t',header = T)
u_w_info$Date<-as.Date(u_w_info$Date, '%Y-%m-%d')
uinfo1<-ggplot(u_w_info, aes(x=Date, y=Country,  color=u_w_info$One))+
  #geom_smooth()+
  geom_jitter(size=1)+
  ggtitle("Position 1")+
  theme_classic()+
  theme(text =element_text(size=10), 
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "3 month")+
  scale_colour_gradient(low="darkblue",high="orange")
uinfo1
uinfo2<-ggplot(u_w_info, aes(x=Date, y=Country,  color=u_w_info$Two))+
  #geom_smooth()+
  geom_jitter(size=1)+
  ggtitle("Position 2")+
  theme_classic()+
  theme(text =element_text(size=10), 
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "3 month")+
  scale_colour_gradient(low="darkblue",high="orange")
uinfo2
uinfo3<-ggplot(u_w_info, aes(x=Date, y=Country,  color=u_w_info$Three))+
  #geom_smooth()+
  geom_jitter(size=1)+
  ggtitle("Position 3")+
  theme_classic()+
  theme(text =element_text(size=10), 
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "3 month")+
  scale_colour_gradient(low="darkblue",high="orange")
uinfo3
uinfo4<-ggplot(u_w_info, aes(x=Date, y=Country,  color=u_w_info$Intergenic))+
  #geom_smooth()+
  geom_jitter(size=1)+
  ggtitle("Intergenic")+
  theme_classic()+
  theme(text =element_text(size=10), 
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "3 month")+
  scale_colour_gradient(low="darkblue",high="orange")
uinfo4

grid.arrange(uinfo1,uinfo2,uinfo3,uinfo4,ncol=2)

uplot1<-ggplot(u1, aes(x=Date, y=Percentage,  color=Percentage))+
  #geom_smooth()+
  geom_point()+
  scale_x_date(date_breaks = "3 month")+
  theme_classic()+
  ggtitle("Position 1")+
  #scale_y_continuous()
  theme(text =element_text(size=10), 
        legend.position = "bottom",
        plot.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  #scale_colour_gradient()
  scale_colour_gradient(low="darkblue",high="orange")
uplot1

uplot2<-ggplot(u2, aes(x=Date, y=Percentage, color=Percentage))+
  #geom_smooth()+
  geom_point()+
  scale_x_date(date_breaks = "3 month")+
  theme_classic()+
  ggtitle("Position 2")+
  #scale_y_continuous()
  theme(text =element_text(size=10), 
        legend.position = "none",
        plot.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_colour_gradient(low="darkblue",high="orange")
  
uplot2
uplot3<-ggplot(u3, aes(x=Date, y=Percentage, color=Percentage))+
  #geom_smooth()+
  geom_point()+
  scale_x_date(date_breaks = "3 month")+
  theme_classic()+
  ggtitle("Position 3")+
  #scale_y_continuous()
  theme(text =element_text(size=10), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(size = 20) )+
  scale_colour_gradient(low="darkblue",high="orange")

uplot3
u4<- subset(u4, u4$Percentage<36.5)
uplot4<-ggplot(u4, aes(x=Date, y=Percentage, color=Percentage))+
  #geom_smooth()+
  geom_point()+
  scale_x_date(date_breaks = "3 month")+
  theme_classic()+
  ggtitle("Intergenic")+
  #scale_y_continuous(limits = c(34.3, 36.5))+
  theme(text =element_text(size=10), 
        legend.position = "none",
        plot.title = element_text(size = 20),axis.text.x = element_text(angle = 90, hjust = 1) )+
  scale_colour_gradient(low="darkblue",high="orange")

uplot4

grid.arrange(uplot1, uplot2, uplot3, uplot4, nrow=2)


c<-  read.csv('wout17/1610_codon_usage_summary_with_dates.csv', header = T, sep=',')
c$Date<-as.Date(c$Date, '%Y-%m-%d')

c$Hum_dist <- log2(c$Frequency/c$Human)

cearly$Frequency<- round(cearly$Frequency, 2)

wide_c<- dcast(c, Sequence + Date ~ Codon, value.var="Frequency")
cearly<-subset(c, Date=='2014-03-17')
cearly$mylabel<- paste(cearly$Codon, cearly$Frequency, sep = ' : ')
ggplot(cearly, aes(x = 1, y = Three, fill = Frequency, label =  mylabel)) + 
  geom_tile() + facet_grid(One ~ Two, labeller = label_both) + 
  scale_fill_gradient(low='grey',high='#D55E00')+
  geom_text(colour = "black") +
  theme_classic()+
  
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("")+
  labs(title = "Codon Frequency Earliest Sample")


codonplot<-ggplot(c, aes(x=Date, y=Hum_dist))+
  #geom_smooth()+
  #geom_smooth()+
  geom_point(size = 0.5)+
  scale_x_date(date_breaks = "3 month")+
  theme_classic()+
  ggtitle("Codon Usage")+
  theme(text =element_text(size=15), plot.title = element_text(size = 20, hjust = 0.5) )
  #scale_color_manual(values = cB)
#scale_color_viridis(discrete = TRUE)
codonplot

disth<-  read.csv('wout17/euclidean_distance_wrt_human.csv', header = T, sep=',')
diste<-  read.csv('wout17/euclidean_distance_wrt_earliest.csv', header = T, sep=',')

disth$Date<-as.Date(disth$Date, '%Y-%m-%d')
diste$Date<-as.Date(diste$Date, '%Y-%m-%d')
distploth<-ggplot(disth, aes(x=Date, y=Distance))+
  #geom_smooth()+
  #geom_smooth()+
  
  geom_point(size = 1, colour="red")+
  geom_smooth(method = "loess",color="black")+
  scale_x_date(date_breaks = "3 month")+
  theme_classic()+
  ylab("Distance from Human Codon Usage")+
  #ggtitle("EBOV Codon Usage:\nEuclidean Distance Relative to Human Codon Usage")+
  theme(text =element_text(size=15), 
        axis.text.x = element_text(size =10,angle = 90, hjust = 1),
        plot.title = element_text(size = 20, hjust = 0.5) )
#scale_color_manual(values = cB)
distploth
distplot<-ggplot(diste, aes(x=Date, y=Distance))+
  #geom_smooth()+
  #geom_smooth()+
  
  geom_point(size = 1, colour="darkblue")+
  geom_smooth(method = "loess",color="black")+
  scale_x_date(date_breaks = "3 month")+
  theme_classic()+
  ylab("Distance from Earliest Virus Sequence")+
  #ggtitle("EBOV Codon Usage:\nEuclidean Distance Relative to Earliest Virus Sequence")+
  theme(text =element_text(size=15), 
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        plot.title = element_text(size = 20, hjust = 0.5) )
#scale_color_manual(values = cB)
distplot
grid.arrange(distplot, distploth, ncol=2,
             top=("EBOV Codon Usage:\nEuclidean Distance"))
