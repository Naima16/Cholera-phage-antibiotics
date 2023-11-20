### GLMM with SNV count as a function of phage relative abundance and antibiotics and their interaction
### This script to produce table S8 and fig 4A
### NM nov 2023



library(plyr)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(jtools)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(interactions)

### SNVs file from inStrain
snv.filter3.biallelic=read.csv('~/snv.filter3.biallelic.csv')

##load antibio
patient_antibio=read.csv('~/metadata_all.csv')

# merge frames together
df.all=df_list %>% reduce(inner_join, by='Sample') 
df.all=merge(snv.filter3.biallelic,patient_antibio,by='Sample')

## remove antibio outliers
df.all.out=df.all[!df.all$Sample %in% c('1394550','D17181815','1319650','D02182446','D17181817','D05181982'), ]

N_mut=df.all.out[df.all.out$mutation_type=='N' & df.all.out$var_freq>0.1,]

counts.mut <- ddply(N_mut, .(N_mut$Sample), nrow)
colnames(counts.mut)=c('Sample','SNV_nb')

## add sxt profile
breadth.all=read.csv('~/breadth.all.csv')
names(breadth.all)[names(breadth.all)=="sample"] <- "Sample"

df.patient=unique(df.all.out[,c('Sample',"CIP" , "AZI" , "DOX" )])
                  
df.2=merge(counts.mut,df.patient,by='Sample')
df.3=merge(df.2,breadth.all[,c('Sample','ICE')],by='Sample')
df.3$ICE=as.factor(df.3$ICE)

### rel abundances of the dominant species 
df.species.otu.comp.otu_vc_icp1=read.csv('df.species.otu.comp.otu_vc_icp1.csv')
##########

df.species.otu.comp.otu_vc_icp1$Sample=rownames(df.species.otu.comp.otu_vc_icp1)

df.all=merge(df.species.otu.comp.otu_vc_icp1,df.3,by='Sample')

df.all$CIP=as.numeric(df.all$CIP)

### GLMMs

datsc=df.all
pvar1='ICP1'
pvar2='Vc'
pvar3='CIP'
pvar4='AZI'
pvar5='DOX'

datsc[pvar1] <- lapply(df.all[pvar1],scale)
datsc[pvar2] <- lapply(df.all[pvar2],scale)
datsc[pvar3] <- lapply(df.all[pvar3],scale)
datsc[pvar4] <- lapply(df.all[pvar4],scale)
datsc[pvar5] <- lapply(df.all[pvar5],scale)

## select which distribution between poisson, nbinom1 and nbinom2 is the best
mod_all=glmmTMB(SNV_nb ~ ICP1  + Vc + CIP + AZI + DOX+
                        Vc:ICP1:CIP + Vc:ICP1:AZI+ Vc:ICP1:DOX,
                data=datsc,
                disp=~1,
                list(family="poisson",link="log"))

summary(mod_all)


mod_all.nb1=glmmTMB(SNV_nb ~ ICP1  + Vc + CIP + AZI+DOX+
                        Vc:ICP1:CIP + Vc:ICP1:AZI+ Vc:ICP1:DOX,
                data=datsc,
                disp=~1,
                family=nbinom1)

mod_all.nb2=glmmTMB(SNV_nb ~ ICP1  + Vc + CIP + AZI+DOX+
                            Vc:ICP1:CIP + Vc:ICP1:AZI+ Vc:ICP1:DOX,
                    data=datsc,
                    disp=~1,
                    family=nbinom2)

AICtab(mod_all,mod_all.nb1,mod_all.nb2)

#### 
mod.1 <- glmmTMB(SNV_nb ~ ICP1  + Vc + CIP + DOX+AZI+
                          ICP1:CIP+
                          ICP1:AZI+
                          Vc:ICP1,
                  data=datsc,
                  disp=~1,
                 family=nbinom2)
summary(mod.1)


mod.2 <- glmmTMB(SNV_nb ~ ICP1  + Vc + CIP + DOX+AZI+ ICP1:AZI+ Vc:ICP1,
                 data=datsc,
                 disp=~1,
                 family=nbinom2)


mod.3 <- glmmTMB(SNV_nb ~ ICP1  + Vc + CIP + DOX+AZI+ Vc:ICP1,
                 data=datsc,
                 disp=~1,
                 family=nbinom2)

mod.4 <- glmmTMB(SNV_nb ~ ICP1  + Vc + CIP + AZI+ Vc:ICP1,
                 data=datsc,
                 disp=~1,
                 family=nbinom2)

mod.5 <- glmmTMB(SNV_nb ~ ICP1  + Vc +  AZI + Vc:ICP1,
                 data=datsc,
                 disp=~1,
                 family=nbinom2)

mod.6<- glmmTMB(SNV_nb ~ ICP1  + Vc + Vc:ICP1,
                 data=datsc,
                 disp=~1,
                family=nbinom2)

mod.7<- glmmTMB(SNV_nb ~  Vc + Vc:ICP1,
                data=datsc,
                disp=~1,
                family=nbinom2)
summary(mod.7)
mod.8<- glmmTMB(SNV_nb ~  Vc , data=datsc,
                disp=~1,
                family=nbinom2)

mod.9<- glmmTMB(SNV_nb ~  ICP1 ,
                data=datsc,
                disp=~1,
                family=nbinom2)


AICtab(mod_all.nb2,mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9)

summary(mod.7)

save(mod.7,file='mod.7.RData')


r2(mod.7)




sims_0 <- simulateResiduals(mod.7)

png('mod.7_dharma.png')
plot(sims_0,quantreg = T)
dev.off()

mod.7.null<- glmmTMB(SNV_nb ~  1,
                data=datsc,
                disp=~1,
                family=nbinom2)

anova(mod.7,mod.7.null)


## plot the glmm

p1=as_grob(effect_plot(mod.7,'Vc',main.title= NULL,
                       rug=FALSE,
                       line.thickness=0.8,
                       y.label='Number of SNVs',
                       x.label='Standardized Vc',
                       interval = F,colors ='gray35' ) +
                   theme_classic2()+
                   theme(axis.text.x = element_text(angle=0,size=7),
                         axis.text.y = element_text(size=7),#
                         legend.title = element_text(color = "black", size = 5),
                         legend.text = element_text(color = "black", size = 5),
                         legend.position = c(0.2,0.8),
                         legend.background =  element_blank(),
                         axis.title=element_text(angle=0,size=8)))



p4=interact_plot(
        mod.7,vary.lty = F,
        colors=c('lightskyblue2','deepskyblue2','blue4'),
        ICP1,
        Vc,
        x.label='Standardized ICP1',
        y.label='Number of SNVs',
        line.thickness = 0.8)+
        ylim(1, 93)+
        xlim(-0.3,4)+
        theme_classic2()+
        theme(axis.text.x = element_text(angle=0,size=7),
              axis.text.y = element_text(size=7),
              legend.position = c(0.2,0.8),
              legend.text =element_text(size=6),
              legend.title =element_text(size=7),
              legend.box.spacing = unit(0.07, "cm"),
              legend.key.size = unit(0.4, "cm"),
              legend.background = element_rect(fill='transparent'),
              axis.title=element_text(angle=0,size=8))

p4


pdf('mod.7.pdf',width = 1.4,height = 3.5,pointsize = 8)
plot_grid(p1,p4,ncol=1,nrow=2, rel_heights=c(1,1),rel_widths =c(2.2,3))
dev.off()
