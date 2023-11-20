#################################### NM nov 2023
#################################### fig 1B from 'Phage predation is a biomarker for disease severity and shapes pathogen genetic diversity in cholera patients', madi et al. 2023

library(ggpubr)
library(ggplot2)
library(FSA)

##reads number to calculate icp1/vc ratio
df.species.otu_vc_icp1=read.csv('~/df.species.otu_vc_icp1.csv')
## because of the ratio, discard null denumeratr
df.species.otu_vc_icp1.Vcpos=df.species.otu_vc_icp1[df.species.otu_vc_icp1$Vc>0,]

df.species.otu_vc_icp1.Vcpos$ICP1_VC_ratio=df.species.otu_vc_icp1.Vcpos$ICP1/df.species.otu_vc_icp1.Vcpos$Vc
df.species.otu_vc_icp1.Vcpos$Sample=rownames(df.species.otu_vc_icp1.Vcpos)

##rel abundance
df.species.otu.comp.otu_vc_icp1=read.csv('~/df.species.otu.comp.otu_vc_icp1.csv')
df.species.otu.comp.otu_vc_icp1$Sample=rownames(df.species.otu.comp.otu_vc_icp1)

df.all.vcplus.ratio=merge(df.species.otu.comp.otu_vc_icp1[,c('Sample','Vc','ICP1')],df.species.otu_vc_icp1.Vcpos[,c('Sample','ICP1_VC_ratio')],by='Sample')
dim(df.all.vcplus.ratio)  ## 323  4

## patient metadata
patient_metadata=read.csv('~/metadata.csv')
df.all.vcplus.ratio.met=merge(df.all.vcplus.ratio,patient_metadata,by='Sample')


levels (df.all.vcplus.ratio.met$Dehydration_Status)[levels (df.all.vcplus.ratio.met$Dehydration_Status) == '1'] <- 'Mild'
levels (df.all.vcplus.ratio.met$Dehydration_Status)[levels (df.all.vcplus.ratio.met$Dehydration_Status) == '2'] <- 'Moderate'
levels (df.all.vcplus.ratio.met$Dehydration_Status)[levels (df.all.vcplus.ratio.met$Dehydration_Status) == '3'] <- 'Severe'

#####
df.all.vcplus.ratio.met$Dehydration_Status <- factor(df.all.vcplus.ratio.met$Dehydration_Status, levels=c('Severe', 'Moderate', 'Mild'))



ggplot(df.all.vcplus.ratio.met, aes(x=Dehydration_Status, y=log10(ICP1_VC_ratio+1))) +
  geom_point(aes(fill=Dehydration_Status),color='black', size=1.8, shape=21,
             position=position_jitter(width=0.21, height=0)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
  xlab('Dehydration')+
  ylab('Log10 (ICP1/cholerae+1)')+
  theme_bw()+
  theme(axis.title=element_text(size=3),
        legend.position = 'none',
        axis.text.y=element_text(size=7),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x=element_text(angle=0,size=7,hjust = 1))+
  scale_fill_manual(values = c( "#FF0000" ,"#F2AD00","#34b1f0ff"))
  #stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
ggsave('Fig1B.pdf',width = 2.7,height = 3)

#### statistical test
##   Kruskal-wallis test + Dunn's post-hoc test
kruskal.test(ICP1_VC_ratio~Dehydration_Status ,
             data=df.all.vcplus.ratio.met)

##   post-hoc test
PT = dunnTest(ICP1_VC_ratio~Dehydration_Status ,
              data=df.all.vcplus.ratio.met,
              method='bh')
PT


