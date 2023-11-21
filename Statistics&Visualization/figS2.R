## script for fig S2

library(readxl)
library(ggplot2)
library(cowplot)


###qPCR
ICP_VC_qpcr=read_excel('~/mHDM_Aggregate_9_12_23_Draft.xlsx',sheet='mHDM_VC_qPCR')

df_ratio=ICP_VC_qpcr[,c('Sample_ID','qPCR_tcpa_AvgCt','qPCR_tcpA_CFUml_CT28','qPCR_ICP1_AvgCt','qPCR_ICP1_PFUml_CT28')]

### species abundances
df.species.otu.comp.otu_vc_icp1=read.csv('~/df.species.otu.comp.otu_vc_icp1.csv')

colnames(df_ratio)[1]='Sample'
df.species.otu.comp.otu_vc_icp1$Sample=rownames(df.species.otu.comp.otu_vc_icp1)

### plot
comon_data=merge(df_ratio,df.species.otu.comp.otu_vc_icp1[,c('Sample','Vc','ICP1','ICP2','ICP3')],by='Sample')

vc_1=ggplot(comon_data , aes(x = Vc,y=qPCR_tcpa_AvgCt)) +
  geom_point(fill='#00A08A',size=1,pch=21)+#,col=mut_type))+aes(fill=mutation_type)
  xlab('Vc (from MG)')+
  ylab('Vc (Ct from qPCR)')+
  scale_fill_manual(values = c( "#5BBCD6","#F2AD00","#FF0000","#00A08A",'pink' ,'black','yellow','purple','blue','cyan','orange'))+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'bottom',
        legend.text =element_text(size=3),
        axis.text.x=element_text(angle=0,size=5,hjust = 1),
        axis.text.y=element_text(size=5,hjust = 1),
        strip.text.x = element_text(size = 4.5))

############
############ replace NEG with 0
comon_data[comon_data=='NEG']='0'
comon_data$qPCR_tcpA_CFUml_CT28=as.numeric(comon_data$qPCR_tcpA_CFUml_CT28)
comon_data$qPCR_ICP1_PFUml_CT28=as.numeric(comon_data$qPCR_ICP1_PFUml_CT28)

vc_ct28=ggplot(comon_data , aes(x = Vc,y=qPCR_tcpA_CFUml_CT28)) +
  geom_point(fill='#00A08A',size=1,pch=21)+#,col=mut_type))+aes(fill=mutation_type)
  xlab('Vc (from MG)')+
  ylab('Vc (CFU from qPCR)')+
  scale_fill_manual(values = c( "#5BBCD6","#F2AD00","#FF0000","#00A08A",'pink' ,'black','yellow','purple','blue','cyan','orange'))+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        axis.text.x=element_text(angle=0,size=5,hjust = 1),
        axis.text.y=element_text(size=5,hjust = 1),
  geom_smooth(method='glm', se=F,color="#FF0000",size=0.6,alpha=0.6) + ##avant ct glm
  stat_cor(size=1.8,label.x.npc = "left",
           label.y.npc = "top",method='pearson')

#####icp1

icp1_1=ggplot(comon_data , aes(x = ICP1,y=qPCR_ICP1_AvgCt)) +
  geom_point(fill='#00A08A',size=1,pch=21)+#,col=mut_type))+aes(fill=mutation_type)
  xlab('ICP1 (from MG)')+
  ylab('ICP1 (Ct from qPCR)')+
  scale_fill_manual(values = c( "#5BBCD6","#F2AD00","#FF0000","#00A08A",'pink' ,'black','yellow','purple','blue','cyan','orange'))+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'bottom',
        legend.text =element_text(size=3),
        axis.text.x=element_text(angle=0,size=5,hjust = 1),
        axis.text.y=element_text(size=5,hjust = 1),
        strip.text.x = element_text(size = 4.5))
 
 icp1_ct28=ggplot(comon_data , aes(x = ICP1,y=qPCR_ICP1_PFUml_CT28)) +
  geom_point(fill='#00A08A',size=1,pch=21)+#,col=mut_type))+aes(fill=mutation_type)
  xlab('ICP1 (from MG)')+
  ylab('ICP1 (PFU from qPCR)')+
  scale_fill_manual(values = c( "#5BBCD6","#F2AD00","#FF0000","#00A08A",'pink' ,'black','yellow','purple','blue','cyan','orange'))+
  theme_bw()+
   theme(axis.title=element_text(size=7),
        legend.position = 'bottom',
        legend.text =element_text(size=3),
        axis.text.x=element_text(angle=0,size=5,hjust = 1),
        axis.text.y=element_text(size=5,hjust = 1),
        strip.text.x = element_text(size = 4.5))+
geom_smooth(method='glm', se=F,color="#FF0000",size=0.6,alpha=0.6)+  ##avant ct glm
stat_cor(size=1.8,label.x.npc = "left",
  label.y.npc = "top",method='pearson')


pdf('FigS2.pdf',width = 4,height = 4,pointsize = 0.5)
plot_grid(vc_1,icp1_1,vc_ct28,icp1_ct28,ncol=2,labels = c('A','B','C','D'), nrow=2,label_size = 8)
dev.off()

