library(tidyverse)
library(tidyr)
library(ggplot2)
library(purrr)
library(wesanderson)
library(forcats)


pa = wes_palettes %>% 
  names()
pal = wes_palette(name = pa[7], n = 5, type = "discrete")
pal2 = wes_palette(name = pa[8], n = 5, type = 'discrete')
pal3 = wes_palette(name = pa[1], n = 5, type = 'discrete')

cols = c(rev(pal2), rev(pal),rev(pal3))
cols = c(pal3,pal2)

mypalette <- colorRampPalette(pal)(5)

### species abundances
df.all.subset=read.csv('~/df.all.subset.csv')
df.all.subset$Sample=rownames(df.all.subset)

### patient metadata and antibiotics concentrations
patient_antibio=read.csv('~/metadata_all.csv')

df_list <- list(df.all.subset, patient_antibio)      

df1=df_list %>% reduce(inner_join, by='Sample')  #full_join

## remove antbx outliers
df2=df1[!df.MFA1$Sample %in% c('1331650','1394550','D17181815','1319650','D02182446','D17181817','D05181982'),]


df.rda=df2[,colnames(df.MFA)%in% c('Vibrio.cholerae','Escherichia.coli','Bacteroides.vulgatus','Bifidobacterium.longum','Shigella.flexneri','Bifidobacterium.breve','Streptococcus.mitis',
                                      'ICP1','ICP2','ICP3','Dehydration_Status','AZI','DOX','CIP','Collection_Month','Area_Code','Age_in_Years',
                                      'Duration_of_Dirrhoea_in_Hrs',"Sex",'Nature_of_Stool')]

rownames(df.rda)=df.MFA$Sample

## change the ref level in prder to plot the dehydration=1 (end of infection) in the rda, the reference is not plotted
levels (df.rda$Dehydration_Status)[levels (df.rda$Dehydration_Status) == '1'] <- 'Mild'
levels (df.rda$Dehydration_Status)[levels (df.rda$Dehydration_Status) == '2'] <- 'Moderate'
levels (df.rda$Dehydration_Status)[levels (df.rda$Dehydration_Status) == '3'] <- 'Severe'

# change the ref level in prder to plot the dehydration=1 (end of infection) in the rda, the reference is not plotted
df.rda$Dehydration_Status <- factor(df.rda$Dehydration_Status, levels=c('Severe','Moderate', 'Mild'))


my.df=df.rda[,c('Dehydration_Status','Vibrio.cholerae','Escherichia.coli',"Bacteroides.vulgatus" ,"Bifidobacterium.longum",'ICP1','ICP2','ICP3',"Bifidobacterium.breve")]
abund1.long <- gather(my.df, species, abund, 'Vibrio.cholerae':'Bifidobacterium.breve', factor_key=TRUE)
abund1.long$species=as.factor(abund1.long$species)


abund1.long$Dehydration_Status=as.factor(abund1.long$Dehydration_Status)


ggplot(abund1.long , aes(x = species,y=abund)) +
  geom_boxplot(aes(fill=Dehydration_Status),color='black')+#,col=mut_type))+,col=BREX_fac#FF0000
  ylab('Relative abundance')+
  facet_wrap(~Dehydration_Status ,ncol=3)+  #scales = "free",,labeller = label_both
   theme_bw()+
  theme(axis.title=element_text(size=3),
        legend.position = 'none',
        legend.text = element_text(size=6),
        legend.title  = element_text(size=7),
        axis.text=element_text(size=7),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x=element_text(angle=60,size=7,hjust = 1,face ='italic'),
        strip.text.x = element_text(size = 7))+
  scale_fill_manual(values = c( "#FF0000" ,"gray","#34b1f0ff"))
#stat_compare_means()
ggsave('Fig1_A.pdf',width = 4,height = 3.5)

## add indicator species (script in RDA.r line 187) with inkscape
