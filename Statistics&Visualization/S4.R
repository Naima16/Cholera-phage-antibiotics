library(readxl)
library(ggplot2)
library(FSA)


# ##qPCR
ICP_VC_qpcr=read_excel('~/mHDM_Aggregate_9_12_23_Draft.xlsx',sheet='mHDM_VC_qPCR')
df_ratio=ICP_VC_qpcr[,c('Sample_ID','qPCR_Ratio_ICP1_VC_CT28','qPCR_tcpA_CFUml_CT28')]
colnames(df_ratio_POS)[1]='Sample'
df_ratio[df_ratio=='NEG']='0'
df_ratio$qPCR_Ratio_ICP1_VC_CT28_num=as.numeric(df_ratio$qPCR_Ratio_ICP1_VC_CT28)
df_ratio$qPCR_tcpA_CFUml_CT28_num=as.numeric(df_ratio$qPCR_tcpA_CFUml_CT28)

#### only vc>0 as in MG manuscript
df_ratio.vc=df_ratio[df_ratio$qPCR_tcpA_CFUml_CT28_num>0,]
colnames(df_ratio.vc)[1]='Sample'

### patient metadata
metadata=read.csv('~/metadata.csv')

df_dehyd=merge(df_ratio.vc,metadata[,c('Sample','Dehydration_Status')])

df_dehyd$Dehydration_Status=as.factor(df_dehyd$Dehydration_Status)
levels (df_dehyd$Dehydration_Status)[levels (df_dehyd$Dehydration_Status) == '1'] <- 'Mild'
levels (df_dehyd$Dehydration_Status)[levels (df_dehyd$Dehydration_Status) == '2'] <- 'Moderate'
levels (df_dehyd$Dehydration_Status)[levels (df_dehyd$Dehydration_Status) == '3'] <- 'Severe'

df_dehyd$Dehydration_Status <- factor(df_dehyd$Dehydration_Status, levels=c('Severe', 'Moderate', 'Mild'))

ggplot(df_dehyd, aes(x=Dehydration_Status, y=log10(qPCR_Ratio_ICP1_VC_CT28_num+1))) +
  geom_point(aes(fill=Dehydration_Status),color='black', size=1, shape=21,
             position=position_jitter(width=0.21, height=0)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20",width = 0.8 ) +
  xlab('Dehydration')+
  ylab('Log10 (ICP1:VC ratio from qPCR)')+
  theme_bw()+
  theme(axis.title=element_text(size=3),
        legend.position = 'none',
        axis.text.y=element_text(size=5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6),
        axis.text.x=element_text(angle=0,size=6,hjust = 0.5))+
  scale_fill_manual(values = c( "#FF0000" ,"#F2AD00","#34b1f0ff"))
ggsave('FigS4.pdf',width = 1.8,height = 2.5)


kruskal.test(qPCR_Ratio_ICP1_VC_CT28_num~Dehydration_Status ,
             data=df_dehyd)

PT = dunnTest(qPCR_Ratio_ICP1_VC_CT28_num~Dehydration_Status ,
              data=df_dehyd,
              method='bh')
PT
