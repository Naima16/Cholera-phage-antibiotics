
###  Script for fig S13
###  NM nov 2023


library(readxl)
library(ggpubr)
library(FSA)
library(cowplot)
library(ggplot2)


breadth.all=read.csv('~/breadth.all.csv')

## qPCR data
ICP_VC_qpcr=read_excel('~/mHDM_Aggregate_9_12_23_Draft.xlsx',sheet='mHDM_VC_qPCR')

df_ratio=ICP_VC_qpcr[,c('Sample_ID','qPCR_ICP1_PFUml_CT28','qPCR_tcpA_CFUml_CT28','qPCR_Ratio_ICP1_VC_CT28','qPCR_ICP2_PFUml_CT28','qPCR_ICP3_PFUmlCT28')]

df_ratio[df_ratio=='NEG'] = '0'

df_ratio$qPCR_ICP1_PFUml_CT28=as.numeric(df_ratio$qPCR_ICP1_PFUml_CT28)
df_ratio$qPCR_tcpA_CFUml_CT28=as.numeric(df_ratio$qPCR_tcpA_CFUml_CT28)
df_ratio$qPCR_ICP2_PFUml_CT28=as.numeric(df_ratio$qPCR_ICP2_PFUml_CT28)
df_ratio$qPCR_ICP3_PFUmlCT28=as.numeric(df_ratio$qPCR_ICP3_PFUmlCT28)
df_ratio$qPCR_Ratio_ICP1_VC_CT28=as.numeric(df_ratio$qPCR_Ratio_ICP1_VC_CT28)


mydf_ratio_vcplus= df_ratio[df_ratio $qPCR_tcpA_CFUml_CT28>0 | df_ratio $qPCR_ICP1_PFUml_CT28>0 | df_ratio $qPCR_ICP2_PFUml_CT28>0 | df_ratio $qPCR_ICP3_PFUmlCT28>0,]
colnames(mydf_ratio_vcplus)[1]='Sample'

df.all=merge(mydf_ratio_vcplus,breadth.all,by='Sample')

df.all$ICE=as.factor(df.all$ICE)
levels (df.all$ICE)[levels (df.all$ICE) == 'No_ICE'] <- 'ICE-'


icp1=ggplot(df.all , aes(x=ICE,y=log10(qPCR_Ratio_ICP1_VC_CT28+1))) +
  geom_point(fill='#858279ff', size=0.7, shape=21, colour="black",
             position=position_jitter(width=0.21, height=0)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
  xlab('Dehydration')+
  ylab('Log10(ICP1:Vc ratio from qPCR)')+
  #theme_classic()+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'right',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
icp1


## icp2
df.all.vcplus=df.all[df.all$qPCR_tcpA_CFUml_CT28>0,]

df.all.vcplus$ICP2_VC_ratio_qPCR=df.all.vcplus$qPCR_ICP2_PFUml_CT28 /df.all.vcplus$qPCR_tcpA_CFUml_CT28

icp2=ggplot(df.all.vcplus , aes(x=ICE,y=log10(ICP2_VC_ratio_qPCR+1))) +
  geom_point(fill='#858279ff', size=0.7, shape=21, colour="black",
             position=position_jitter(width=0.21, height=0)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
  xlab('Dehydration')+
  ylab('Log10(ICP2:Vc ratio from qPCR)')+
  #theme_classic()+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'right',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
icp2

### icp3
df.all.vcplus$ICP3_VC_ratio_qPCR=df.all.vcplus$qPCR_ICP3_PFUmlCT28 /df.all.vcplus$qPCR_tcpA_CFUml_CT28

icp3=ggplot(df.all.vcplus , aes(x=ICE,y=log10(ICP3_VC_ratio_qPCR+1))) +
  geom_point(fill='#858279ff', size=0.7, shape=21, colour="black",
             position=position_jitter(width=0.21, height=0)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
  xlab('Dehydration')+
  ylab('Log10(ICP3:Vc ratio from qPCR)')+
  #theme_classic()+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'right',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
icp3

## all phages

df.all.vcplus$all_phages=df.all.vcplus$qPCR_ICP1_PFUml_CT28+df.all.vcplus$qPCR_ICP2_PFUml_CT28+df.all.vcplus$qPCR_ICP3_PFUmlCT28
df.all.vcplus$phage_VC_ratio_qPCR=df.all.vcplus$all_phages/df.all.vcplus$qPCR_tcpA_CFUml_CT28

phages=ggplot(df.all.vcplus , aes(x=ICE,y=log10(phage_VC_ratio_qPCR+1))) +
  geom_point(fill='#858279ff', size=0.7, shape=21, colour="black",
             position=position_jitter(width=0.21, height=0)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
  xlab('Dehydration')+
  ylab('Log10(phages:Vc ratio from qPCR)')+
  #theme_classic()+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'right',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
phages


##########.  plot by dehydration

my_metadata=read.csv('~/metadata.csv')
df.all.met=merge(df.all, my_metadata,by='Sample')

levels (df.all.met$Dehydration_Status)[levels (df.all.met$Dehydration_Status) == '1'] <- 'Mild'
levels (df.all.met$Dehydration_Status)[levels (df.all.met$Dehydration_Status) == '2'] <- 'Moderate'
levels (df.all.met$Dehydration_Status)[levels (df.all.met$Dehydration_Status) == '3'] <- 'Severe'


## change the ref level in prder to plot the dehydration=1 (end of infection) in the rda, the reference is not plotted
df.all.met$Dehydration_Status <- factor(df.all.met$Dehydration_Status, levels=c('Severe','Moderate', 'Mild'))


icp1_dehydr1=ggplot(df.all.met , aes(x=ICE,y=log10(qPCR_Ratio_ICP1_VC_CT28+1))) +
  geom_boxplot(outlier.colour='black',outlier.fill = '#858279ff',outlier.shape = 21, 
               aes(fill=Dehydration_Status), colour="grey20",outlier.size = 0.7) +
  xlab('Dehydration')+
  ylab('Log10(ICP1:Vc ratio from qPCR)')+
  #theme_classic()+
  theme_bw()+
  scale_fill_manual(values = c( "#FF0000","#F2AD00","#5BBCD6","#00A08A" ))+
  theme(axis.title=element_text(size=7),
        legend.position = c(0.45,0.8),
        legend.text =element_text(size=5),
        legend.title =element_text(size=6),
        legend.box.spacing = unit(0.1, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill='transparent'),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1))
     
        #legend.title = element_text(angle=0,size=7),
        #legend.text  = element_text(angle=0,size=7),
        
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
icp1_dehydr1


### remove vc=0 because ratio = inf
df.all.met2=merge(df.all.vcplus,my_metadata.antibio,by='sample')

levels (df.all.met2$Dehydration_Status)[levels (df.all.met2$Dehydration_Status) == '1'] <- 'Mild'
levels (df.all.met2$Dehydration_Status)[levels (df.all.met2$Dehydration_Status) == '2'] <- 'Moderate'
levels (df.all.met2$Dehydration_Status)[levels (df.all.met2$Dehydration_Status) == '3'] <- 'Severe'


##change the ref level in prder to plot the dehydration=1 (end of infection) in the rda, the reference is not plotted
df.all.met2$Dehydration_Status <- factor(df.all.met2$Dehydration_Status, levels=c('Severe','Moderate', 'Mild'))

icp2_dehydr2=ggplot(df.all.met2 , aes(x=ICE,y=log10(ICP2_VC_ratio_qPCR+1))) +
  geom_boxplot( outlier.colour='black',outlier.fill = '#858279ff',outlier.shape = 21,
                aes(fill=Dehydration_Status), colour="grey20",outlier.size = 0.7)  +
  scale_fill_manual(values = c( "#FF0000","#F2AD00","#5BBCD6","#00A08A" ))+
  xlab('Dehydration')+
  ylab('Log10(ICP2:Vc ratio from qPCR)')+
  theme_classic()+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        #legend.position = 'right',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1),
        legend.position = 'none') 
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
icp2_dehydr2


icp3_dehydr2=ggplot(df.all.met2 , aes(x=ICE,y=log10(ICP3_VC_ratio_qPCR+1))) +
  geom_boxplot(outlier.colour='black',outlier.fill = '#858279ff',outlier.size = 0.7,outlier.shape = 21,aes(fill=Dehydration_Status), colour="grey20") +
  scale_fill_manual(values = c( "#FF0000","#F2AD00","#5BBCD6","#00A08A" ))+
  xlab('Dehydration')+
  ylab('Log10(ICP3:Vc ratio from qPCR)')+
  theme_classic()+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1),
        legend.position = 'none') 
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
icp3_dehydr2

all_dehydr2=ggplot(df.all.met2 , aes(x=ICE,y=log10(phage_VC_ratio_qPCR+1))) +
  geom_boxplot(outlier.colour='black',outlier.fill = '#858279ff',outlier.shape = 21,  
               aes(fill=Dehydration_Status), colour="grey20",outlier.size = 0.7)  +
  scale_fill_manual(values = c( "#FF0000","#F2AD00","#5BBCD6","#00A08A" ))+
  xlab('Dehydration')+
  ylab('Log10(phages:Vc ratio from qPCR)')+
  theme_classic()+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        #legend.position = 'right',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=6,hjust = 1),
        legend.position = 'none') 
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
all_dehydr2

pdf('FigS14.pdf',width = 5.5,height = 5,pointsize = 0.5)
plot_grid(icp1,icp2,icp3,phages,
          icp1_dehydr1,icp2_dehydr2,icp3_dehydr2,all_dehydr2,
          labels = c('A','B','C','D','E','F','G','H'),label_size = 8,
          ncol=4,nrow=2)
dev.off()
