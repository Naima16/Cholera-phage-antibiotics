### script for figure S17
### NM nov 2023


library(plyr)
library(readxl)
library(dplyr)
library(ggplot2)

#### SNVs from inStrain
snv.filter3.biallelic=read.csv('~/snv.filter3.biallelic.csv')

### mean mutation freq/sample
freq.mut <- ddply(snv.filter3.biallelic, .(snv.filter3.biallelic$Sample, snv.filter3.biallelic$mutation_type), summarize,mean_freq=mean(var_freq))
colnames(freq.mut)=c('Sample','mutation_type','mean_freq')


##### per SXT+ or SXT-
breadth.all=read.csv('~/breadth.all.csv')

df.1=merge(freq.mut,breadth.all[,c('Sample','ICE')],by='Sample')
df.1$ICE=as.factor(df.1$ICE)

########
df.2=df.1[df.1$mutation_type %in% c('N','S','I'),]
df.2$mutation_type <- factor(df.2$mutation_type, levels=c("N", "S","I"))

####qPCR
ICP_VC_qpcr=read_excel('~/mHDM_Aggregate_9_12_23_Draft.xlsx',sheet='mHDM_VC_qPCR')
df_ratio=ICP_VC_qpcr[,c('Sample_ID','qPCR_ICP1_PFUml_CT28','qPCR_tcpA_CFUml_CT28')]

df_ratio[df_ratio=='NEG'] = '0'

df_ratio$qPCR_ICP1_PFUml_CT28=as.numeric(df_ratio$qPCR_ICP1_PFUml_CT28)
df_ratio$qPCR_tcpA_CFUml_CT28=as.numeric(df_ratio$qPCR_tcpA_CFUml_CT28)
colnames(df_ratio)[1]='Sample'

df.2.qpcr=merge(df.2,df_ratio,by='Sample')

df.2.qpcr=df.2.qpcr %>%
  mutate(icp1_abund = case_when(
    qPCR_ICP1_PFUml_CT28 > 0 ~ 'plus',
    qPCR_ICP1_PFUml_CT28 == 0  ~ 'min',
  ))

df.2.qpcr$ICE <- factor(df.2.qpcr$ICE, levels=c('No_ICE', 'ind6', 'ind5'))

qPCR_plot=ggplot(df.2.qpcr , aes(x = icp1_abund,y=log10(mean_freq)) )+
  geom_boxplot(aes(fill=ICE))+#
  ylab('Log10(mean mutation freq)')+
  facet_wrap(mutation_type~ICE,scales = "free",ncol=3,labeller = label_both)+
  #scale_fill_manual(values = c( "#5BBCD6","#F2AD00","#FF0000","#00A08A",'pink' ))+
  scale_fill_manual(values = c(  "#a739caff","#858279ff" , "#afbeafff")) +
  theme_bw()+
  theme(axis.title.y=element_text(size=8),
        strip.text = element_text(size=6),
        legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=5,hjust=0.5),
        axis.text.y=element_text(size=5),
        strip.text.x = element_blank()
        ) +
  scale_x_discrete(labels=c('ICP1-', 'ICP1+')) #+
qPCR_plot

ggsave('FigS17.pdf',width = 2.3,height = 3.8)

