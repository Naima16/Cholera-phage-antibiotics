
## Script for S3
## NM nov 2023

library(cowplot)
library(dplyr)
library(FSA)
library(ggplot2)


my_path='~/naima/cholera_project'
df.species.otu_vc_icp1=('my_path/df.species.otu_vc_icp1.csv',sep=',',header=T)

## because of the ratio, discard null denumeratr
df.species.otu_vc_icp1.Vcpos=df.species.otu_vc_icp1[df.species.otu_vc_icp1$Vc>0,]

df.species.otu_vc_icp1.Vcpos$ICP1_VC_ratio=df.species.otu_vc_icp1.Vcpos$ICP1/df.species.otu_vc_icp1.Vcpos$Vc

## rel abundance
df.species.otu.comp.otu_vc_icp1=read.table('my_path/df.species.otu.comp.otu_vc_icp1.csv',sep=',',header=T)

df.all.vcplus.ratio=merge(df.species.otu.comp.otu_vc_icp1[,c('Sample','Vc','ICP1')],df.species.otu_vc_icp1.Vcpos[,c('Sample','ICP1_VC_ratio')],by='Sample')

## patient metadata
patient_metadata=read.csv('my_path/metadata.csv')

df.all.vcplus.ratio.met=merge(df.all.vcplus.ratio,patient_metadata,by='Sample')

df.all.vcplus.ratio.met$Dehydration_Status=as.factor(df.all.vcplus.ratio.met$Dehydration_Status)

levels(df.all.vcplus.ratio.met$Dehydration_Status)
levels (df.all.vcplus.ratio.met$Dehydration_Status)[levels (df.all.vcplus.ratio.met$Dehydration_Status) == '1'] <- 'Mild'
levels (df.all.vcplus.ratio.met$Dehydration_Status)[levels (df.all.vcplus.ratio.met$Dehydration_Status) == '2'] <- 'Moderate'
levels (df.all.vcplus.ratio.met$Dehydration_Status)[levels (df.all.vcplus.ratio.met$Dehydration_Status) == '3'] <- 'Severe'

#####
df.all.vcplus.ratio.met$Dehydration_Status <- factor(df.all.vcplus.ratio.met$Dehydration_Status, levels=c('Severe', 'Moderate', 'Mild'))


#### 24h  dirrhoea treshold

df.all.vcplus.ratio.met=df.all.vcplus.ratio.met %>%
  mutate(Diarrhea = case_when(
    Duration_of_Dirrhoea_in_Hrs>=24  ~ ">=24h",
    #Duration_of_Dirrhoea_in_Hrs>48 &  Duration_of_Dirrhoea_in_Hrs<72 ~ "mid",
    Duration_of_Dirrhoea_in_Hrs<24~ "<24h"
  ))
df.all.vcplus.ratio.met$Diarrhea=as.factor(df.all.vcplus.ratio.met$Diarrhea)


vc_24=ggplot(df.all.vcplus.ratio.met , aes(x=Dehydration_Status,y=Vc)) +  #Dehydration_Status Duration_of_Dirrhoea_in_Hrs
  geom_boxplot(aes(fill=Dehydration_Status),width=0.5)+#,col=mut_type))+
  #scale_fill_manual(values = c( "#FF0000","#F2AD00","#00A08A" ))+
  theme_bw()+
  facet_wrap(~Diarrhea,scales = "free",ncol=3,labeller = label_both)+  #labeller = label_both,
  stat_cor(size=3,label.x.npc = 'left',
           label.y.npc = 'top')+
  scale_fill_manual(values = c( "#FF0000" ,"#F2AD00","#34b1f0ff"))+
  
  theme(axis.title=element_text(size=9),
        strip.text=element_text( size = 6),
        legend.position = 'none',#####'bottom',
        axis.text.x=element_text(angle=0,size=6,hjust = 0.4),
        axis.title.x=element_blank(),
        axis.title.y= element_text(color = "black", size = 7.5),
        axis.text.y=element_text(size=5,hjust = 1))
vc_24

###. stat tests
## > 24h
df1=df.all.vcplus.ratio.met[df.all.vcplus.ratio.met$Diarrhea==">=24h",]

kruskal.test(Vc~Dehydration_Status ,
             data=df1)

PT = dunnTest(Vc~Dehydration_Status ,
              data=df1,
              method="bh")
PT

##<24h
df2=df.all.vcplus.ratio.met[df.all.vcplus.ratio.met$Diarrhea=="<24h",]


kruskal.test(Vc~Dehydration_Status ,
             data=df2)

#########
######### ICP1
icp1_24=ggplot(df.all.vcplus.ratio.met , aes(x=Dehydration_Status,y=ICP1)) +  #Dehydration_Status Duration_of_Dirrhoea_in_Hrs
  geom_boxplot(aes(fill=Dehydration_Status),width=0.5)+#,col=mut_type))+
  #scale_fill_manual(values = c( "#FF0000","#F2AD00","#00A08A" ))+
  theme_bw()+
  facet_wrap(~Diarrhea,scales = "free",ncol=3,labeller = label_both)+  #labeller = label_both,
  stat_cor(size=3,label.x.npc = 'left',
           label.y.npc = 'top')+
  scale_fill_manual(values = c( "#FF0000" ,"#F2AD00","#34b1f0ff"))+
  
  theme(axis.title=element_text(size=9),
        strip.text=element_text( size = 6),
        legend.position = 'none',#####'bottom',
        axis.text.x=element_text(angle=0,size=6,hjust = 0.4),
        axis.title.x=element_blank(),
        axis.title.y= element_text(color = "black", size = 7.5),
        axis.text.y=element_text(size=5,hjust = 1))
icp1_24

##### 72h  dirrhoea treshold

df.all.vcplus.ratio.met=df.all.vcplus.ratio.met %>%
  mutate(Diarrhea = case_when(
    Duration_of_Dirrhoea_in_Hrs>=72  ~ ">=72h",
    #Duration_of_Dirrhoea_in_Hrs>48 &  Duration_of_Dirrhoea_in_Hrs<72 ~ "mid",
    Duration_of_Dirrhoea_in_Hrs<72~ "<72h"
  ))
df.all.vcplus.ratio.met$Diarrhea=as.factor(df.all.vcplus.ratio.met$Diarrhea)


vc_72=ggplot(df.all.vcplus.ratio.met , aes(x=Dehydration_Status,y=Vc)) +  
  geom_boxplot(aes(fill=Dehydration_Status),width=0.5)+#,col=mut_type))+
  theme_bw()+
  facet_wrap(~Diarrhea,scales = "free",ncol=3,labeller = label_both)+  
  stat_cor(size=3,label.x.npc = 'left',
           label.y.npc = 'top')+
  scale_fill_manual(values = c( "#FF0000" ,"#F2AD00","#34b1f0ff"))+
  theme(axis.title=element_text(size=9),
        strip.text=element_text( size = 6),
        legend.position = 'none',#####'bottom',
        axis.text.x=element_text(angle=0,size=6,hjust = 0.4),
        axis.title.x=element_blank(),
        axis.title.y= element_text(color = "black", size = 7.5),
        axis.text.y=element_text(size=5,hjust = 1))
vc_72


icp1_72=ggplot(df.all.vcplus.ratio.met , aes(x=Dehydration_Status,y=ICP1)) +  
  geom_boxplot(aes(fill=Dehydration_Status),width=0.5)+#,col=mut_type))+
  theme_bw()+
  facet_wrap(~Diarrhea,scales = "free",ncol=3,labeller = label_both)+  
  stat_cor(size=3,label.x.npc = 'left',
           label.y.npc = 'top')+
  scale_fill_manual(values = c( "#FF0000" ,"#F2AD00","#34b1f0ff"))+
  theme(axis.title=element_text(size=9),
        strip.text=element_text( size = 6),
        legend.position = 'none',#####'bottom',
        axis.text.x=element_text(angle=0,size=6,hjust = 0.4),
        axis.title.x=element_blank(),
        axis.title.y= element_text(color = "black", size = 7.5),
        axis.text.y=element_text(size=5,hjust = 1))
icp1_72

##post-hoc test
## > 72
df3=df.all.vcplus.ratio.met[df.all.vcplus.ratio.met$Diarrhea==">=72h",]

kruskal.test(Vc~Dehydration_Status ,
             data=df3)
dunnTest(Vc~Dehydration_Status ,
         data=df3,
         method="bh")


##<72h
df4=df.all.vcplus.ratio.met[df.all.vcplus.ratio.met$Diarrhea=="<72h",]

kruskal.test(Vc~Dehydration_Status ,
             data=df4)
PT = dunnTest(Vc~Dehydration_Status ,
              data=df4,
              method="bh")
PT


kruskal.test(ICP1~Dehydration_Status ,
             data=df4)
PT = dunnTest(ICP1~Dehydration_Status ,
              data=df4,
              method="bh")
PT

###
pdf('S3.pdf',width = 4.5,height = 5,pointsize = 0.5)
plot_grid(vc_24,icp1_24,vc_72,icp1_72,ncol=2,nrow=2,labels = c('A','B','C','D'),label_size = 7)#####, rel_heights=c(5,2),rel_widths =c(2,0.5))
dev.off()


