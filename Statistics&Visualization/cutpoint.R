## cutpoint analysis
## NM janv 2024

rm(list=ls())
library(dplyr)
library(cutpointr)

df.species.otu_vc_icp1=read.csv('~/df.species.otu_vc_icp1.csv')

df.species.otu_vc_icp1.Vcpos=df.species.otu_vc_icp1[df.species.otu_vc_icp1$Vc>0,]

df.species.otu_vc_icp1.Vcpos$ICP1_VC_ratio=df.species.otu_vc_icp1.Vcpos$ICP1/df.species.otu_vc_icp1.Vcpos$Vc

df.species.otu_vc_icp1.Vcpos$Sample=rownames(df.species.otu_vc_icp1.Vcpos)

##rel abundance
df.species.otu.comp.otu_vc_icp1=read.csv('~/df.species.otu.comp.otu_vc_icp1.csv')

df.species.otu.comp.otu_vc_icp1$Sample=rownames(df.species.otu.comp.otu_vc_icp1)

df.all.vcplus.ratio1=merge(df.species.otu.comp.otu_vc_icp1[,c('Sample','Vc','ICP1')],df.species.otu_vc_icp1.Vcpos[,c('Sample','ICP1_VC_ratio')],by='Sample')

my_metadata=read.csv('~/metadata.csv')

df1=merge(df.all.vcplus.ratio1[,c('Sample','Vc','ICP1','ICP1_VC_ratio')], my_metadata[,c('Sample','Dehydration_Status')],by='Sample')

levels (df1$Dehydration_Status)[levels (df1$Dehydration_Status) == '1'] <- 'Mild'
levels (df1$Dehydration_Status)[levels (df1$Dehydration_Status) == '2'] <- 'Moderate'
levels (df1$Dehydration_Status)[levels (df1$Dehydration_Status) == '3'] <- 'Severe'

#####

df1$Dehydration_Status <- factor(df1$Dehydration_Status, levels=c('Severe', 'Moderate', 'Mild'))

df1=df1 %>%
  mutate(dehydration = case_when(
    Dehydration_Status=='Mild' ~ 'mild',
    Dehydration_Status=='Moderate' | Dehydration_Status == 'Severe' ~ "not_mild"
  ))

cp <- cutpointr(df1, ICP1_VC_ratio, dehydration, 
                method = maximize_metric, pos_class = 'not_mild', neg_class = "mild", 
                metric = sum_sens_spec,boot_runs = 1000) 

summary(cp)