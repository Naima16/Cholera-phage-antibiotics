### scripts for S15 and S16
### NM nov 2023

library(dplyr)
library(tidyr)
library(cowplot)
library(FSA)

snv.filter3.biallelic=read.csv('~/snv.filter3.biallelic.csv')

### mutator-non mutator table
mutator_status=read.table('~/mutatorstatus.csv',header=T,sep=',')

transversion_transition_df <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(transversion_transition_df)=c('Sample','A->G','A->C','A->T','G->T','G->C','G->A','C->T','C->G','C->A','T->C','T->G','T->A')

sample_set=unique(snv.filter3.biallelic$Sample)

for (i in sample_set)
{ 
  df_sample=snv.filter3.biallelic[snv.filter3.biallelic$Sample==i,]  
  
  nb_AG=nrow(df_sample[df_sample$ref_base =='A' & df_sample$var_base=='G',])
  nb_AC=nrow(df_sample[df_sample$ref_base =='A' & df_sample$var_base=='C',])
  nb_AT=nrow(df_sample[df_sample$ref_base =='A' & df_sample$var_base=='T',])
  
  nb_GT=nrow(df_sample[df_sample$ref_base =='G' & df_sample$var_base=='T',])
  nb_GC=nrow(df_sample[df_sample$ref_base =='G' & df_sample$var_base=='C',])
  nb_GA=nrow(df_sample[df_sample$ref_base =='G' & df_sample$var_base=='A',])
  
  nb_CT=nrow(df_sample[df_sample$ref_base =='C' & df_sample$var_base=='T',])
  nb_CG=nrow(df_sample[df_sample$ref_base =='C' & df_sample$var_base=='G',])
  nb_CA=nrow(df_sample[df_sample$ref_base =='C' & df_sample$var_base=='A',])
  
  nb_TC=nrow(df_sample[df_sample$ref_base =='T' & df_sample$var_base=='C',])
  nb_TG=nrow(df_sample[df_sample$ref_base =='T' & df_sample$var_base=='G',])
  nb_TA=nrow(df_sample[df_sample$ref_base =='T' & df_sample$var_base=='A',])
  
  row_to_add=c(i,nb_AG,nb_AC,nb_AT,nb_GT,nb_GC,nb_GA,nb_CT,nb_CG,nb_CA,nb_TC,nb_TG,nb_TA)
  transversion_transition_df[nrow(transversion_transition_df) + 1, ] <- row_to_add 

}

mut=c('A->G','A->C','A->T','G->T','G->C','G->A','C->T','C->G','C->A','T->C','T->G','T->A')
mut=as.data.frame(mut)


mut=mut %>% 
  mutate (mut_type = case_when (
    mut %in% c('A->G','G->A','T->C','C->T')  ~ 'Transition',
    !mut %in% c('A->G','G->A','T->C','C->T') ~ 'Transversion'
  ))


transversion_transition_df.long <- gather(transversion_transition_df, SNV_sort, SNV_count, 'A->G':'T->A', factor_key=TRUE)
transversion_transition_df.long$SNV_count=as.numeric(transversion_transition_df.long$SNV_count)
transversion_transition_df.long.1=merge(transversion_transition_df.long,mutator_status,by='Sample')

colnames(mut)[1]='SNV_sort'
transversion_transition_df.long2=merge(transversion_transition_df.long.1,mut,by='SNV_sort')
transversion_transition_df.long2$SNV_sort<- factor(transversion_transition_df.long2$SNV_sort,levels = c('A->G','G->A','C->T','T->C', 'G->T','C->A', 'A->C','A->T','G->C','C->G','T->G','T->A'))

transversion_transition_df.long2$hyper_Mutator <- factor(transversion_transition_df.long2$hyper_Mutator, levels=c("yes","no","high_snv_noDnaRep","low_snv_DnaRep" ))

data.1=unique(transversion_transition_df.long2[,c('Sample','hyper_Mutator','SNV_sort','SNV_sort')])

C=ggplot(transversion_transition_df.long2[transversion_transition_df.long2$hyper_Mutator %in% c("no","yes" ),] , aes(x = SNV_sort,y=SNV_count)) +
  geom_boxplot(aes(fill=hyper_Mutator))+
  xlab('')+
  ylab('SNV count')+
  facet_wrap(~hyper_Mutator,ncol=2,scales = "free")+ 
  scale_fill_manual(values = c( "#FF0000","#F2AD00","#00A08A" ,"#5BBCD6"))+
  theme_bw()+
  theme(axis.title=element_text(size=8),
        legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=60,size=6,vjust = 0.5,colour = 'black'),
        axis.text.y=element_text(size=5,hjust = 1),
        strip.text.x =element_blank())

C


#######
df.species.otu.comp.otu_vc_icp1=read.table('~/df.species.otu.comp.otu_vc_icp1.csv',sep=',',header=T)

data.2=unique(transversion_transition_df.long2[,c('Sample','hyper_Mutator')])

df.1=merge(df.species.otu.comp.otu_vc_icp1,data.2,by='Sample')      

df.1$hyper_Mutator <- factor(df.1$hyper_Mutator, levels=c("yes","no","high_snv_noDnaRep","low_snv_DnaRep" ))

A=ggplot(df.1 , aes(x=hyper_Mutator,y=Vc)) +
  geom_boxplot(outlier.colour=NA, aes(fill=hyper_Mutator), colour="black") +
  scale_fill_manual(values = c( '#FF0000','#F2AD00','#858279ff','#afbeafff'))+
  ylab('Vc')+
  theme_bw()+
  theme(axis.title=element_text(size=9),
        legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=5,hjust = 1),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
A

B=ggplot(df.1 , aes(x=hyper_Mutator,y=ICP1)) +
  geom_boxplot(outlier.colour=NA, aes(fill=hyper_Mutator), colour="black") +
  scale_fill_manual(values = c( '#FF0000','#F2AD00','#858279ff','#afbeafff'))+
  ylab('ICP1')+
  theme_bw()+
  theme(axis.title=element_text(size=9),
        legend.position = 'none',
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=7,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=5,hjust = 1))
  stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
B


A_B=cowplot::plot_grid(A,B,ncol=1,nrow=2,labels = c('A','B'),label_size = 8)

pdf('S15.pdf',width = 5,height = 3.2,pointsize = 0.5)
cowplot::plot_grid(A_B,C,ncol=2,nrow=1,rel_widths =c(1,2),labels=c('','C'),label_size = 8)
dev.off()

## stat test
kruskal.test(Vc~hyper_Mutator ,
             data=df.1)

# Kruskal-Wallis rank sum test
# 
# data:  Vc by hyper_Mutator
# Kruskal-Wallis chi-squared = 31.001, df = 3, p-value = 8.496e-07

PT = dunnTest(Vc~hyper_Mutator ,
              data=df.1,
              method='bh')
PT

# Comparison          Z      P.unadj        P.adj
# 1 high_snv_noDnaRep - low_snv_DnaRep  2.4249327 1.531123e-02 3.062245e-02
# 2             high_snv_noDnaRep - no  0.5051969 6.134206e-01 6.134206e-01
# 3                low_snv_DnaRep - no -2.3498577 1.878059e-02 2.817088e-02
# 4            high_snv_noDnaRep - yes  3.5163199 4.375734e-04 1.312720e-03
# 5               low_snv_DnaRep - yes -1.0560248 2.909569e-01 3.491483e-01
# 6                           no - yes  4.8930279 9.929634e-07 5.957781e-06

### icp1
kruskal.test(ICP1~hyper_Mutator ,
             data=df.1)

# Kruskal-Wallis rank sum test
# 
# data:  ICP1 by hyper_Mutator
# Kruskal-Wallis chi-squared = 6.9979, df = 3, p-value = 0.07197


#####
### sxt breadth
breadth.all=read.csv('~/breadth.all.csv')
breadth.all$ICE =as.factor(breadth.all$ICE )
levels(breadth.all$ICE)[levels(breadth.all$ICE)=='No_ICE']='ICE-'
breadth.all$ICE <- factor(breadth.all$ICE, levels=c("ICE-","ind6",'ind5' ))

##  juste les samples qui ont passé inStrain filter
snv.filter3.biallelic=read.csv('/Users/naimamadi/cholera_nov2023/new_Figs_16nov2023/github_data/snv.filter3.biallelic.csv')

load('~/snv.filter3.biallelic.RData')
length(unique(snv.filter3.biallelic$Sample)) 

breadth.all.instrain=breadth.all[breadth.all$Sample %in% snv.filter3.biallelic$Sample,]
dim(breadth.all.instrain) #131

ind5=ggplot(breadth.all.instrain , aes(x = ICE,y=ind5_breadth) )+
  geom_boxplot(outlier.colour='gray41', aes(fill=ICE), colour="black") +
  scale_fill_manual(values = c('#a739caff','#858279ff','#afbeafff'))+
   ylab('ind5 breadth of coverage')+
  theme_bw()+
  theme(axis.title=element_text(size=9),
        legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=8,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=5,hjust = 1))
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
ind5

ind6=ggplot(breadth.all.instrain , aes(x = ICE,y=breadth.ind6) )+
  geom_boxplot(outlier.colour='gray41', aes(fill=ICE), colour="black") +
  scale_fill_manual(values = c('#a739caff','#858279ff','#afbeafff'))+
   ylab('ind6 breadth of coverage')+
  theme_bw()+
  theme(axis.title=element_text(size=9),
        legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,size=8,hjust = 0.8,colour = 'black'),
        axis.text.y=element_text(size=5,hjust = 1))
stat_compare_means(size=2,aes(label = paste0("p = ", after_stat(p.format))))
ind6

pdf('S16.pdf',width =4,height = 2,pointsize = 3)
cowplot::plot_grid(ind5,ind6,nrow = 1,labels = c('A','B'),label_size = 8)
dev.off()


####   statistical test ind6
kruskal.test(breadth.ind6~ICE ,
             data=breadth.all.instrain)

# Kruskal-Wallis rank sum test
# 
# data:  ICP1 by hyper_Mutator
# Kruskal-Wallis chi-squared = 6.9979, df = 3, p-value = 0.07197

PT = dunnTest(breadth.ind6~ICE ,
              data=breadth.all.instrain,
              method='bh')
PT
# Comparison         Z      P.unadj        P.adj
# 1 ICE- - ind5 -6.474614 9.505465e-11 1.425820e-10
# 2 ICE- - ind6 -9.298522 1.424109e-20 4.272328e-20
# 3 ind5 - ind6 -5.456287 4.861948e-08 4.861948e-08

### ind5
kruskal.test(ind5_breadth~ICE ,
             data=breadth.all.instrain)


PT = dunnTest(ind5_breadth~ICE ,
              data=breadth.all.instrain,
              method='bh')

PT
# Comparison         Z      P.unadj        P.adj
# 1 ICE- - ind5 -8.565140 1.079445e-17 3.238335e-17
# 2 ICE- - ind6 -1.780957 7.491955e-02 7.491955e-02
# 3 ind5 - ind6  5.474065 4.398277e-08 6.597416e-08
