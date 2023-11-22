
## RDA plotted in Figure 1C and S6 as well as indicator species analysis
###


library(wesanderson)
library(purrr)
library(forcats)
library(tidyverse)
library(vegan)
library(dplyr)
library(ggrepel)
library(indicspecies)

pa = wes_palettes %>% 
  names()
pal = wes_palette(name = pa[7], n = 5, type = "discrete")
pal2 = wes_palette(name = pa[8], n = 5, type = 'discrete')
pal3 = wes_palette(name = pa[1], n = 5, type = 'discrete')

cols = c(rev(pal2), rev(pal),rev(pal3))
cols = c(pal3,pal2)

mypalette <- colorRampPalette(pal)(5)

#### 
df.all.subset=read.csv('~/df.all.subset.csv')

patient_antibio=read.csv('~/metadata_all.csv')

df_list <- list(df.all.subset, patient_antibio)      

df.rda1=df_list %>% reduce(inner_join, by='Sample')  #full_join

## remove antbiotics outliers
df.rda=df.rda1[!df.rda1$Sample %in% c('1331650','1394550','D17181815','1319650','D02182446','D17181817','D05181982'),]


df.rda=df.rda[,colnames(df.rda)%in% c('Vibrio.cholerae','Escherichia.coli','Bacteroides.vulgatus','Bifidobacterium.longum','Shigella.flexneri',
                                      'Bifidobacterium.breve','Streptococcus.mitis','ICP1','ICP2','ICP3',
                                      'Dehydration_Status','AZI','DOX','CIP','Collection_Month','Area_Code','Age_in_Years',
                                      'Duration_of_Dirrhoea_in_Hrs',"Sex",'Nature_of_Stool','Vomiting')]

df.rda$CIP=as.numeric(df.rda$CIP)
rownames(df.rda)=df.rda$Sample

metadata.df$Dehydration_Status=as.factor(metadata.df$Dehydration_Status)

levels (df.rda$Dehydration_Status)[levels (df.rda$Dehydration_Status) == '1'] <- 'Mild'
levels (df.rda$Dehydration_Status)[levels (df.rda$Dehydration_Status) == '2'] <- 'Moderate'
levels (df.rda$Dehydration_Status)[levels (df.rda$Dehydration_Status) == '3'] <- 'Severe'


sp.df=df.rda[,1:7] 
rownames(sp.df)=df.rda$Sample


metadata.df=df.rda[,8:21]
rownames(metadata.df)=df.rda$Sample

sp.df[1:7] <- sapply(sp.df[1:7],as.numeric)
sapply(sp.df[1:7],class)


##log chord
log.df=decostand(sp.df,method='log')
physeq.hel=decostand(log.df,method='norm')

rda_1=vegan::rda(physeq.hel~.,metadata.df)#RDA Analysis

RsquareAdj(rda_1)

anova.cca(rda_1, step = 10000, by = "term")

## forward selection 
fwd.sel <- ordiR2step(rda(physeq.hel ~ CIP+DOX+AZI+ICP1+ICP2+ICP3,
                          data = metadata.df), # lower model limit (simple!)
                      scope = formula(rda_1), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = F, # can't surpass the "full" model's R2
                      pstep = 1000000000000000000000,
                      trace = FALSE) # change to TRUE to see the selection process!


##
fwd.sel$call

rda_1.rda.signif$call

RsquareAdj(rda_1.rda.signif) 
anova.cca(rda_1.rda.signif, step = 10000)
anova.cca(rda_1.rda.signif, step = 10000, by = "term")
anova.cca(rda_1.rda.signif, step = 1000, by = "axis")

## plot
ii=summary(rda_1.rda.signif)  
sp=as.data.frame(ii$species[,1:2])*5
st=as.data.frame(ii$sites[,1:2])*3
yz=as.data.frame(ii$biplot[,1:2])*4
yz1=yz
yz1=yz[rownames(yz) %in% c('CIP','DOX','AZI','ICP1','ICP2','ICP3','Dehydration_StatusMild','Dehydration_StatusSevere','Vomiting2','Age_in_Years'),]

options(ggrepel.max.overlaps = Inf)

p=ggplot() +
  geom_point(data = st,aes(RDA1,RDA2,color=metadata.df$Dehydration_Status),size=1.2,alpha=0.8) +#,,color=mymetadata1$VC,shape=metad$Antibotic
  scale_color_manual(values = c( "#F2AD00", '#34b1f0ff',"#FF0000"))
  p=p+geom_segment(data = sp,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                 arrow = arrow(angle=10,length = unit(0.2,"cm"),
                 type = "closed"),linetype=1, size=0.25,colour = "#000078ff",alpha=0.8)  
p=p+ 
  geom_text_repel(data = sp,aes(RDA1,RDA2,label=row.names(sp)),size=2,
                  color="#000078ff",segment.color = 'grey',fontface ='italic',alpha=1) 

p=p+geom_segment(data = yz1,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                 arrow = arrow(angle=10,length = unit(0.2,"cm"),
                               type = "closed"),linetype=1, size=0.25,colour = "#006000ff",alpha=1)  ##scale=0.4 ds num1 version

p=p+geom_text_repel(data = yz1,aes(RDA1,RDA2,label=row.names(yz1)),size=2,color='#006000ff',segment.color = 'grey',fontface ='plain')+ #bold
  labs(x=paste("RDA 1 (", format(100 *ii$cont[[1]][2,1], digits=4), "%)", sep=""),
       y=paste("RDA 2 (", format(100 *ii$cont[[1]][2,2], digits=4), "%)", sep=""))+
  geom_hline(yintercept=0,linetype=3,size=0.5,color='gray') + 
  geom_vline(xintercept=0,linetype=3,size=0.5,color='gray')+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill = F)+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'), 
        legend.key.height = unit(0.5, 'cm'), 
        legend.key.width = unit(0.5, 'cm'),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6),
        legend.position = 'top',
        axis.title = element_text(size=7),
        axis.text = element_text(size=5))+
  labs(colour="Dehydration status")
p
ggsave('rda_whole.pdf',height = 5,width = 5)

###plot only species, phages, dehydartion, vomit, age

ii=summary(rda_1.rda.signif)  #View analysis results
sp=as.data.frame(ii$species[,1:2])*5
st=as.data.frame(ii$sites[,1:2])*3
yz=as.data.frame(ii$biplot[,1:2])*4


options(ggrepel.max.overlaps = Inf)

p=ggplot() + 
  geom_point(data = st,aes(RDA1,RDA2,color=metadata.df$Dehydration_Status),size=1.8,alpha=0.6) 
  scale_color_manual(values = c( "#F2AD00", "#FF0000", "#00A08A",'darkorchid1')) 
p=p+geom_segment(data = sp,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                 arrow = arrow(angle=10,length = unit(0.2,"cm"),
                 type = "closed"),linetype=1, size=0.25,colour = "#000078ff",alpha=0.8) 

p=p+ 
  geom_text_repel(data = sp,aes(RDA1,RDA2,label=row.names(sp)),size=3,
                  color="#000078ff",segment.color = 'grey',fontface ='italic',alpha=0.9) ##size=2.5 ds version1

yz1=yz[rownames(yz) %in% c('CIP','DOX','AZI','ICP1','ICP2','ICP3','Dehydration_Status1','Dehydration_Status3','Vomiting2','Age_in_Years'),]
p=p+geom_segment(data = yz1,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                 arrow = arrow(angle=10,length = unit(0.2,"cm"),
                 type = "closed"),linetype=1, size=0.25,colour = "#006000ff",alpha=0.8)  ##scale=0.4 ds num1 version

p=p+geom_text_repel(data = yz1,aes(RDA1,RDA2,label=row.names(yz1)),size=3,color='#006000ff',segment.color = 'grey',fontface ='plain')+ #bold
  labs(x=paste("RDA 1 (", format(100 *ii$cont[[1]][2,1], digits=4), "%)", sep=""),
       y=paste("RDA 2 (", format(100 *ii$cont[[1]][2,2], digits=4), "%)", sep=""))+
  
  geom_hline(yintercept=0,linetype=3,size=0.5,color='gray') + 
  geom_vline(xintercept=0,linetype=3,size=0.5,color='gray')+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill = F)+
  theme_bw()+
  theme(legend.key.size = unit(0.1, 'cm'), #change legend key size
        legend.key.height = unit(0.1, 'cm'), #change legend key height
        legend.key.width = unit(0.1, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6),
        legend.position = 'top')+
  labs(colour="Dehydration")
p

ggsave('rda.pdf',height = 7,width = 7)


##### indicator species

df.rda=df.rda  %>%
  relocate(Dehydration_Status, .after = ICP2)

df.indicspecies=df.rda[,1:38]
groups=df.indicspecies$Dehydration_Status

## no combination of groups
indval.sing = multipatt(df.indicspecies[,-38], groups,  func='r.g',duleg = TRUE, control = how(nperm=9999)) 
summary(indval.sing)
summary(indval.sing,indvalcomp=TRUE)
indval.sing$sign
