#### script for Fig4 and Fig S12
#### NM nov 2023


library(ggpubr)
library(cowplot)
library(ggplot2)

#### loqd ICE profile
load('~/breadth.all.RData')

### species composition
load('~/df.species.otu_vc_icp1.RData')
df.species.otu_vc_icp1$sample=rownames(df.species.otu_vc_icp1)

df.all=merge(breadth.all,df.species.otu_vc_icp1[,c('sample','ICP1','ICP2','ICP3','Vc')],by='sample')
df.all$ICP1_VC_ratio=df.all$ICP1/df.all$Vc

df.all$ICE=as.factor(df.all$ICE)
levels (df.all$ICE)[levels (df.all$ICE) == 'No_ICE'] <- 'ICE-'

icp1=ggplot(df.all , aes(x=ICE,y=log10(ICP1_VC_ratio+1))) +
  geom_point(fill='#858279ff', size=0.9, shape=21, colour="black",
             position=position_jitter(width=0.21, height=0)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
  xlab('Dehydration')+
  ylab('Log10(ICP1:Vc ratio)')+
  theme_bw()+
  theme(axis.title=element_text(size=8),
        legend.position = 'right',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
  stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
icp1

df.all$ICP2_VC_ratio=df.all$ICP2/df.all$Vc

icp2=ggplot(df.all , aes(x=ICE,y=log10(ICP2_VC_ratio+1))) +
  geom_point(fill='#858279ff', size=0.9, shape=21, colour="black",
             position=position_jitter(width=0.21, height=0)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
  xlab('Dehydration')+
  ylab('Log10(ICP2:Vc ratio)')+
  theme_bw()+
  theme(axis.title=element_text(size=8),
        legend.position = 'right',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
  stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
icp2

df.all$ICP3_VC_ratio=df.all$ICP3/df.all$Vc
icp3=ggplot(df.all , aes(x=ICE,y=log10(ICP3_VC_ratio+1))) +
  geom_point(fill='#858279ff', size=0.9, shape=21, colour="black",
             position=position_jitter(width=0.21, height=0)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
  xlab('Dehydration')+
  ylab('Log10(ICP3:Vc ratio)')+
  theme_bw()+
  theme(axis.title=element_text(size=8),
        legend.position = 'right',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
  stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
icp3

df.all$all_phages=df.all$ICP1+df.all$ICP2+df.all$ICP3
df.all$phage_VC_ratio=df.all$all_phages/df.all$Vc

phages=ggplot(df.all , aes(x=ICE,y=log10(phage_VC_ratio+1))) +
  geom_point(fill='#858279ff', size=0.9, shape=21, colour="black",
             position=position_jitter(width=0.21, height=0)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
  xlab('Dehydration')+
  ylab('Log10(phages:Vc ratio)')+
  theme_classic()+
  theme_bw()+
  theme(axis.title=element_text(size=8),
        legend.position = 'right',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
  stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
phages

### Statistic tests

## ICP1
kruskal.test(ICP1_VC_ratio~ICE ,
             data=df.all)

PT = dunnTest(ICP1_VC_ratio~ICE ,
              data=df.all,
              method='bh')

##icp2
kruskal.test(ICP2_VC_ratio~ICE ,
             data=df.all)

PT = dunnTest(ICP2_VC_ratio~ICE ,
              data=df.all,
              method='bh')
##icp3
kruskal.test(ICP3_VC_ratio~ICE ,
             data=df.all)

PT = dunnTest(ICP3_VC_ratio~ICE ,
              data=df.all,
              method='bh')

##all
kruskal.test(phage_VC_ratio~ICE ,
             data=df.all)

PT = dunnTest(phage_VC_ratio~ICE ,
              data=df.all,
              method='bh')

####
########
## add dehydration
load('~/my_metadata1.antibio.RData')
colnames(df.all)[1]='Sample'

df.all.dehyd=merge(df.all,my_metadata1.antibio,by='Sample')

levels (df.all.dehyd$ICE)[levels (df.all.dehyd$ICE) == 'No_ICE'] <- 'ICE-'
levels (df.all.dehyd$Dehydration_Status)[levels (df.all.dehyd$Dehydration_Status) == '1'] <- 'Mild'
levels (df.all.dehyd$Dehydration_Status)[levels (df.all.dehyd$Dehydration_Status) == '2'] <- 'Moderate'
levels (df.all.dehyd$Dehydration_Status)[levels (df.all.dehyd$Dehydration_Status) == '3'] <- 'Severe'

### change levels order
df.all.dehyd$Dehydration_Status <- factor(df.all.dehyd$Dehydration_Status, levels=c('Severe','Moderate', 'Mild'))


icp1_dehydr1=ggplot(df.all.dehyd , aes(x=ICE,y=log10(ICP1_VC_ratio+1))) +
   geom_boxplot(outlier.colour='black',outlier.fill = '#858279ff',outlier.shape = 21, 
               aes(fill=Dehydration_Status), colour="grey20",outlier.size = 0.7) +
  xlab('Dehydration')+
  ylab('Log10(ICP1:Vc ratio)')+
  theme_bw()+
  scale_fill_manual(values = c( "#FF0000","#F2AD00","#5BBCD6","#00A08A" ))+
  theme(axis.title=element_text(size=7),
        legend.position = c(0.35,0.9),
        legend.text =element_text(size=5),
        legend.title =element_blank(),
        legend.box.spacing = unit(0.1, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill='transparent'),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
icp1_dehydr1

icp2_dehydr1=ggplot(df.all.dehyd , aes(x=ICE,y=log10(ICP2_VC_ratio+1))) +
  geom_boxplot(outlier.colour='black',outlier.fill = '#858279ff',outlier.shape = 21, 
               aes(fill=Dehydration_Status), colour="grey20",outlier.size = 0.7) +
  xlab('Dehydration')+
  ylab('Log10(ICP2:Vc ratio)')+
  theme_bw()+
  scale_fill_manual(values = c( "#FF0000","#F2AD00","#5BBCD6","#00A08A" ))+
  theme(axis.title=element_text(size=7),
        legend.position = 'none',
        legend.text =element_text(size=5),
        legend.title =element_text(size=6),
        legend.box.spacing = unit(0.1, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill='transparent'),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
icp2_dehydr1


icp3_dehydr1=ggplot(df.all.dehyd , aes(x=ICE,y=log10(ICP3_VC_ratio+1))) +
  geom_boxplot(outlier.colour='black',outlier.fill = '#858279ff',outlier.shape = 21, 
               aes(fill=Dehydration_Status), colour="grey20",outlier.size = 0.7) +
  xlab('Dehydration')+
  ylab('Log10(ICP3:Vc ratio)')+
  theme_bw()+
  scale_fill_manual(values = c( "#FF0000","#F2AD00","#5BBCD6","#00A08A" ))+
  theme(axis.title=element_text(size=7),
        legend.position = c(0.35,0.9),
        legend.text =element_text(size=5),
        legend.title =element_blank(),
        legend.box.spacing = unit(0.1, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill='transparent'),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
icp3_dehydr1

phage_dehydr1=ggplot(df.all.dehyd , aes(x=ICE,y=log10(phage_VC_ratio+1))) +
  geom_boxplot(outlier.colour='black',outlier.fill = '#858279ff',outlier.shape = 21, 
               aes(fill=Dehydration_Status), colour="grey20",outlier.size = 0.7) +
  xlab('Dehydration')+
  ylab('Log10(phages:Vc ratio)')+
  theme_bw()+
  scale_fill_manual(values = c( "#FF0000","#F2AD00","#5BBCD6","#00A08A" ))+
  theme(axis.title=element_text(size=7),
        legend.position = 'none',
        legend.text =element_text(size=5),
        legend.title =element_text(size=6),
        legend.box.spacing = unit(0.1, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill='transparent'),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
phage_dehydr1

pdf('Fig4.pdf',width = 3,height = 4.3,pointsize = 0.5)
plot_grid(icp1,phages,icp1_dehydr1,phage_dehydr1, ncol=2,nrow=2,labels = c('A','B','C','D'),label_size = 7)######, rel_heights=c(5,2),rel_widths =c(2,0.5))
dev.off()

pdf('FigS12.pdf',width = 3,height = 4.3,pointsize = 0.5)
plot_grid(icp2,icp3,icp2_dehydr1,icp3_dehydr1, ncol=2,nrow=2,labels = c('A','B','C','D'),label_size = 7)######, rel_heights=c(5,2),rel_widths =c(2,0.5))
dev.off()
