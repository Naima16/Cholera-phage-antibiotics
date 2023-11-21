
### fig S5C
### NM nov 2023

library(dplyr)
library(plyr)
library(purrr)
library(ggplot2)

## SNVs from inStrain
snv.filter3.biallelic=read.csv('~/snv.filter3.biallelic.csv')

### species composition
df.all=read.csv('~/df.species.otu.comp.otu_vc_icp1.csv')

## SXT profile
breadth.all=read.csv('~/breadth.all.csv')

### mean mutation freq/sample
freq.mut <- ddply(snv.filter3.biallelic, .(snv.filter3.biallelic$Sample, snv.filter3.biallelic$mutation_type), summarize,mean_freq=mean(var_freq))
colnames(freq.mut)=c('Sample','mutation_type','mean_freq')

df_list <- list(freq.mut,df.all[,c('Sample','ICP1')], breadth.all[,c('Sample','ICE')])      
df.1=df_list %>% reduce(inner_join, by='Sample')  ## purr library

df.1$ICE=as.factor(df.1$ICE)
df.1$ICE <- factor(df.1$ICE, levels=c("No_ICE", "ind6","ind5"))

df.2=df.1[df.1$mutation_type %in% c('N','S','I'),]
df.2$mutation_type <- factor(df.2$mutation_type, levels=c("N", "S","I"))

df.2=df.2 %>%
  mutate(icp1_abund = case_when(
    ICP1>=0.001 ~ 'plus',
    ICP1<0.001 ~ 'min',
  ))

ggplot(df.2 , aes(x = icp1_abund,y=log10(mean_freq)) )+
  geom_boxplot(aes(fill=ICE))+#
  ylab('Log10(mean mutation freq)')+
  facet_wrap(mutation_type~ICE,scales = "free",ncol=3,labeller = label_both)+
  scale_fill_manual(values = c(  "#a739caff", "#858279ff", "#afbeafff")) +
  theme_bw()+
  theme(axis.title.y=element_text(size=8),
        legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=6,hjust=0.25),
        axis.text.y=element_text(size=6),
        strip.text.x = element_blank()) +
  scale_x_discrete(labels=c('ICP1-', 'ICP1+')) +
  stat_compare_means(size=1.7,aes(label = paste0("p = ", after_stat(p.format))),vjust=2,hjust=0.3)
ggsave('Fig5C.pdf',width = 2.5,height = 5)




