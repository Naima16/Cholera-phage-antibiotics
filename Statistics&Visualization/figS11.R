### script for S11
### NM nov 2023

library(ggpubr)
library(ggplot2)
library(sjPlot)
library(cowplot)

breadth.all=read.table('~/breadth.all.csv',sep=',',header=T)
breadth.all$ICE=as.factor(breadth.all$ICE)
levels(breadth.all$ICE)[levels(breadth.all$ICE)=='No_ICE']='ICE-'

plot_A=ggplot(breadth.all, aes(x =ind5_breadth ,y=breadth.ind6,fill=ICE)) + #
  geom_point(size=1.4,pch=21)+
  geom_abline(slope=1, intercept=0,color='#FF0000',lwd=0.8,alpha=0.6)+
  geom_hline(yintercept=90,color='#FF0000',linetype="dashed",alpha=0.5)+
  geom_vline(xintercept=90,color='#FF0000',linetype="dashed",alpha=0.5)+
  ylab('Coverage breadth ind6')+
  xlab('Coverage breadth ind5')+
  scale_fill_manual(values = c( '#afbeafff','#858279ff','#a739caff' ))+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = c(0.26,0.95),
        legend.direction = 'horizontal',
        legend.text = element_text(size=5),
        legend.title  = element_blank(),
        axis.text=element_text(size=4),
        axis.text.x=element_text(size=6,hjust = 0.5),
        strip.text = element_text(size =5))+
  scale_x_continuous(breaks=c(0, 10,20,30,40,50,60,70,80,90,100)) +
  scale_y_continuous(breaks=c(0, 10,20,30,40,50,60,70,80,90,100))
plot_A


### 
df.species.otu.comp.otu_vc_icp1=read.table('~/df.species.otu.comp.otu_vc_icp1.csv',sep=',',header=T)
df.all1=merge(breadth.all,df.species.otu.comp.otu_vc_icp1[,c('Sample','Vc')],by='Sample')

## instrain output
ICP1_snv.filter3.biallelic=read.csv('~/snv.filter3.biallelic.csv')
df.all2=df.all1[df.all1$Sample %in% snv.filter3.biallelic$Sample,]

inbox_icp1=ggplot(df.all2 , aes(x=ICE,y=Vc)) +
  geom_boxplot(aes(fill=ICE),col='black',outlier.alpha = 0.7,width=0.6)+
  xlab('')+
  ylab('Vc')+
  scale_fill_manual(values = c( '#afbeafff','#858279ff','#a739caff' ))+
  theme_bw()+
  theme(axis.title=element_text(size=7),
        legend.position = 'none',
        axis.text.x=element_text(size=6,hjust = 1),
        axis.text.y=element_text(size=4,hjust = 1),
        plot.margin=margin(t = 30, r = 10, b = 30, l = 0, unit = 'pt')
        )+
  
  stat_compare_means(size=1.6,aes(label = paste0("P = ", after_stat(p.format))),label.x=3,label.y=1)
inbox_icp1

kruskal.test( Vc~ICE, data=df.all2)  ## p=0.36


pdf('S11.pdf',width = 4,height = 2.8,pointsize = 0.5)
cowplot::plot_grid(plot_A,inbox_icp1,ncol=2,nrow=1,rel_widths =c(2,1.3))
dev.off()



