## script for S7 and S8
## NM nov 2023



library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)


#### ARGs : antibiotics resistance genes table as predicted by deepARG
all.arg=read.csv('~/merge_quant_subtype_deepARG.csv',check.names = F)

names(all.arg)[names(all.arg)=="16s-NormalizedReadCount"] <- "NormalizedReadCount_16s"
names(all.arg)[names(all.arg)=="ARG-group"] <- "ARG_group"

all.arg=all.arg[,colnames(all.arg) != 'ReadCount']

ARG_df <- spread(all.arg, ARG_group, NormalizedReadCount_16s)

ARG_df[is.na(ARG_df)] <- 0
rownames(ARG_df)=ARG_df$Sample

ARG_df=ARG_df[,colnames(ARG_df) != 'Sample']

#### antibio
antibio=read.csv('~/github_data/metadata_all.csv')

df.arg.antib=merge(all.arg,antibio,by='Sample')
df.arg.antib$CIP=as.numeric(df.arg.antib$CIP)
range(df.arg.antib$CIP)

### - outliers
df.arg.antib=df.arg.antib[!df.arg.antib$Sample %in% c('1394550','D17181815','1319650','D02182446','D17181817','D05181982'), ]


########## normalized args
df.species.otu_vc_icp1=read.table('~/df.species.otu_vc_icp1.csv',sep=',',header=T)

ratio_vc_icp1=df.species.otu_vc_icp1[,c('Sample','Vc','ICP1')]
colnames(ratio_vc_icp1)[2]='Vc_reads'
colnames(ratio_vc_icp1)[3]='ICP1_reads'

df.all1=merge(df.arg.antib,ratio_vc_icp1[c('Sample','Vc_reads','ICP1_reads')],by='Sample')
df.all1$CIP=as.numeric(df.all1$CIP)
range(df.all1$CIP)

df.arg.antib.clust=df.all1[df.all1$ARG_group %in% c('MEXI','TET34','TET35','DFRA1',
                                                    'DFRA16','APH_3____I','APH_6__I',
                                                    'CHLORAMPHENICOL_EXPORTER','VIBRIO_CHOLERAE_OMPT',
                                                    'VIBRIO_CHOLERAE_OMPU',
                                                    'VIBRIO_CHOLERAE_VARG','FLOR'),]

df.arg.antib.clust=df.arg.antib.clust %>%
  mutate(CIP_fact = case_when(
    CIP< 0.063 ~ "ND",
    CIP >= 0.063 ~ "D",
  ))

df.arg.antib.clust.vc=df.arg.antib.clust[df.arg.antib.clust$Vc_reads>0,]
df.arg.antib.clust.vc$Vc_Normalized_read=df.arg.antib.clust.vc$NormalizedReadCount_16s/df.arg.antib.clust.vc$Vc_reads

cip_anaerobic=ggplot(df.arg.antib.clust.vc , aes(x = CIP_fact,y=log10(Vc_Normalized_read))) +
  geom_boxplot(aes(fill=CIP_fact),outlier.size = 0.8,outlier.alpha = 0.3,width = 0.6)+ 
  facet_wrap(~ARG_group,scales = "free",ncol=4)+
  xlab('CIP anaerobic')+
  scale_fill_manual(values = c( "#00A08A",'#F2AD00'))+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'none',
        legend.text = element_text(size=6),
        legend.title  = element_text(size=7),
        axis.text=element_text(size=5),
        axis.text.x=element_text(angle=0,size=5,hjust = 1),
        strip.text = element_text(size = 5)) #+
cip_anaerobic

## cip aerobic
df.arg.antib.clust.vc=df.arg.antib.clust.vc %>%
  mutate(CIP_fact_aerobic = case_when(
    CIP< 0.016 ~ "ND",
    CIP >= 0.016 ~ "D",
  ))

cip_erobic=ggplot(df.arg.antib.clust.vc , aes(x = CIP_fact_aerobic,y=log10(Vc_Normalized_read))) +
  geom_boxplot(aes(fill=CIP_fact_aerobic),outlier.size = 0.8,outlier.alpha = 0.3,width=0.6)+#,col=mut_type))+aes(fill=mutation_type)
  facet_wrap(~ARG_group,scales = "free",ncol=4)+
  xlab('CIP aerobic')+
  scale_fill_manual(values = c( "#00A08A",'#F2AD00',"#F2AD00","#FF0000","#00A08A" ))+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'none',
        legend.text = element_text(size=6),
        legend.title  = element_text(size=7),
        axis.text=element_text(size=5),
        axis.text.x=element_text(angle=0,size=5,hjust = 1),
        strip.text = element_text(size = 5))#+
cip_erobic

pdf('S8.pdf',width = 5.5,height = 6,pointsize = 0.5)
plot_grid(cip_anaerobic,cip_erobic,ncol=2,nrow=1,labels = c('A','B'),label_size = 8)######, rel_heights=c(5,2),rel_widths =c(2,0.5))
dev.off()

##

##   aerobic
df.arg.antib.clust.vc=df.arg.antib.clust.vc %>%
  mutate(AZI_fact = case_when(
    AZI < 1 ~ "ND",
    AZI >= 1 ~ "D",
  ))


azi_erobic=ggplot(df.arg.antib.clust.vc , aes(x = AZI_fact,y=log10(Vc_Normalized_read))) +
  geom_boxplot(aes(fill=AZI_fact),outlier.size = 0.8,outlier.alpha = 0.3,width = 0.6)+#,col=mut_type))+aes(fill=mutation_type)
  facet_wrap(~ARG_group,scales = "free",ncol=4)+
  xlab('AZI (aerobic)')+
  scale_fill_manual(values = c( "#00A08A",'#F2AD00',"#F2AD00","#FF0000","#00A08A" ))+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'none',
        legend.text = element_text(size=6),
        legend.title  = element_text(size=7),
        axis.text=element_text(size=5),
        axis.text.x=element_text(angle=0,size=5,hjust = 1),
        strip.text = element_text(size = 5))#+
azi_erobic

## anaerobic
df.arg.antib.clust.vc=df.arg.antib.clust.vc %>%
  mutate(AZI_fact_anaerobic = case_when(
    AZI < 8 ~ "ND",
    AZI >= 8 ~ "D",
  ))

azi_anaerobic=ggplot(df.arg.antib.clust.vc , aes(x = AZI_fact_anaerobic,y=log10(Vc_Normalized_read))) +
  geom_boxplot(aes(fill=AZI_fact_anaerobic),outlier.size = 0.8,outlier.alpha = 0.3,width = 0.6)+#,col=mut_type))+aes(fill=mutation_type)
  facet_wrap(~ARG_group,scales = "free",ncol=4)+
  xlab('AZI (anaerobic)')+
  scale_fill_manual(values = c( "#00A08A",'#F2AD00',"#F2AD00","#FF0000","#00A08A" ))+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'none',
        legend.text = element_text(size=6),
        legend.title  = element_text(size=7),
        axis.text=element_text(size=5),
        axis.text.x=element_text(angle=0,size=5,hjust = 1),
        strip.text = element_text(size = 5))#+
azi_anaerobic


pdf('S7.pdf',width = 5.5,height = 6,pointsize = 0.5)
plot_grid(azi_anaerobic,azi_erobic,ncol=2,nrow=1,labels = c('A','B'),label_size = 8)######, rel_heights=c(5,2),rel_widths =c(2,0.5))
dev.off()

##  HB correction
p=0.0047 #0.056
p=0.0023 ##0.03
p=0.00011 #0.0013
p=0.0063 # 0.075
p=0.0039  #0.05
p.adjust(p, method = 'BH', n = 12)


