### script for figure S18 and S19
### NM nov 20023


library(plyr)
library(purrr)  #reduce
library(dplyr)  #%>%
library(ggplot2)

ICP1_snv.filter3.biallelic=read.csv('~/ICP1_snv.filter3.biallelic.csv')
head(ICP1_snv.filter3.biallelic)

## add ICE
breadth.all=read.table('~/breadth.all.csv',sep=',',header=T)

## high freq mutations
ICP1_snv.filter3.biallelic_high=ICP1_snv.filter3.biallelic[ICP1_snv.filter3.biallelic$var_freq>0.1,]

### count snvs per sample per mut type
counts.mut <- ddply(ICP1_snv.filter3.biallelic, .(ICP1_snv.filter3.biallelic$Sample, ICP1_snv.filter3.biallelic$mutation_type), nrow)
colnames(counts.mut)=c('Sample','mutation_type','nb_snv')

df.species.otu.comp.otu_vc_icp1=read.table('~/df.species.otu.comp.otu_vc_icp1.csv',sep=',',header=T)

df_list <- list(counts.mut,df.species.otu_vc_icp1[,c('Sample','Vc','ICP1')], breadth.all[,c('Sample','ICE')])      
df.2=df_list %>% reduce(inner_join, by='Sample')  ## purr library

####
df.2$ICP1_Vc_ratio=df.2$ICP1/df.2$Vc
df.2$mutation_type=factor(df.2$mutation_type,levels=c('N','S','I'))

##ICP1
df.2$ICE <- factor(df.2$ICE, levels=c('No_ICE', 'ind6', 'ind5'))

ggplot(df.2[df.2$mutation_type %in% c('N','S','I'),] , aes(x=log10(ICP1_Vc_ratio+1), y=nb_snv,fill=ICE)) + 
  geom_point(pch=21,size=1.3) +
  ylab('SNV number')+
  xlab('log10(ICP1:Vc ratio)')+
  scale_fill_manual(values = c("#a739caff","#858279ff", "#afbeafff","#F2AD00" ))+
  theme_bw()+
  geom_smooth(method = 'lm',se=F,col='gray',alpha=0.5)+ #method = 'gam'
  facet_wrap(mutation_type~ICE,nrow=3)+ #,labeller = label_both
  theme(axis.title=element_text(size=7),
        legend.position = 'none',
        axis.text.x=element_text(size=4.5,hjust=0.5),
        axis.text.y=element_text(size=4.5),
        strip.text.x = element_blank())
# ggpubr::stat_cor(size=1.7,label.x.npc = "left",
#                  label.y.npc = "top",method='spearman',
#                  color = 'black',r.accuracy=0.01,p.accuracy=0.001,
#                  show.legend=T)
ggsave('S18.pdf',width = 2.9,height = 3.3)


### S19
freq.mut <- ddply(ICP1_snv.filter3.biallelic_high, .(ICP1_snv.filter3.biallelic_high$Sample, ICP1_snv.filter3.biallelic_high$mutation_type), summarize,mean_freq=mean(var_freq))
colnames(freq.mut)=c('Sample','mutation_type','mean_freq')

df_list <- list(freq.mut,df.species.otu_vc_icp1[,c('Sample','Vc','ICP1')], breadth.all[,c('Sample','ICE')])      
df.1=df_list %>% reduce(inner_join, by='Sample')  ## purr library

####
df.1$ICP1_Vc_ratio=df.1$ICP1/df.1$Vc
df.1$mutation_type=factor(df.1$mutation_type,levels=c('N','S','I'))

##ICP1
df.1$ICE <- factor(df.1$ICE, levels=c('No_ICE', 'ind6', 'ind5'))
ggplot(df.1[df.1$mutation_type %in% c('N','S','I'),] , aes(x=log10(ICP1_Vc_ratio+1), y=log10(mean_freq),fill=ICE)) + #
  geom_point(pch=21,size=1.3) +
  ylab('Log10 mean freq (high freq)')+
  xlab('log10(ICP1:Vc ratio)')+
  scale_fill_manual(values = c("#a739caff","#858279ff", "#afbeafff","#F2AD00" ))+
  theme_bw()+
  geom_smooth(method = 'lm',se=F,col='gray',alpha=0.5)+ #method = 'gam'
  facet_wrap(mutation_type~ICE,nrow=3)+ #,labeller = label_both
  theme(axis.title=element_text(size=7),
        legend.position = 'none',
        axis.text.x=element_text(size=4.5,hjust=0.5),
        axis.text.y=element_text(size=4.5),
        strip.text.x = element_blank())
 # ggpubr::stat_cor(size=1.7,label.x.npc = "left",
                  # label.y.npc = "top",method='spearman',
                  # color = 'black',r.accuracy=0.01,p.accuracy=0.001,
                  # show.legend=T)
ggsave('S19.pdf',width = 2.9,height = 3.3)

