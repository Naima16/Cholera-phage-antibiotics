### GAM with the average frequency mutation in Vc as a function of ICP1, antibitics and the ICE profile.
### these results are showed in Figure 4B and Table S9.
### NM nov 2023


library(plyr)
library(mgcv)
library(itsadug)



snv.filter3.biallelic=read.csv('~/snv.filter3.biallelic.csv',sep=',',header=T)
freq.mut <- ddply(snv.filter3.biallelic, .(snv.filter3.biallelic$Sample, snv.filter3.biallelic$mutation_type), summarize,mean_freq=mean(var_freq))

colnames(freq.mut)=c('Sample','mutation_type','mean_freq')

patient_antibio=read.csv('~/metadata_all.csv')

df.all=merge(freq.mut,patient_antibio,by='Sample')

##remove antibio outliers
df.all.out=df.all[!df.all$Sample %in% c('1331650','1394550','D17181815','1319650','D02182446','D17181817','D05181982'), ]

### dominant species with rel abundances
df.species.otu.comp.otu_vc_icp1=read.table('df.species.otu.comp.otu_vc_icp1.csv',sep=',',header=T)
##########

df.all=merge(df.species.otu.comp.otu_vc_icp1[,c('Sample','Vc','ICP1','ICP2','ICP3')],df.all.out,by='Sample')

### add SXT status
breadth.all=read.csv('~/breadth.all.csv')

df.all.sxt=merge(df.all,breadth.all[,c('Sample','ICE')],by='Sample')
df.all.sxt$ICE=as.factor(df.all.sxt$ICE)
df.all.sxt$mutation_type=as.factor(df.all.sxt$mutation_type)

df.all.sxt1=df.all.sxt[df.all.sxt$mutation_type %in% c('I','N','S'),]

df.all.sxt1$ICE.by.mut <- with(df.all.sxt1, interaction(mutation_type,ICE,
                                                         drop=TRUE)) 

### GAMs
mod1 <- gam(mean_freq ~ s(ICP1,by=ICE.by.mut) ,
              data=df.all.sxt1,
              family=betar(link='logit'))
summary(mod1)

save(mod1,file='mod1.RData')

mod4=gam(mean_freq ~ s(ICP1,by=ICE.by.mut)+
                   s(AZI,by=mutation_type)+ 
                   s(CIP,by=mutation_type),
                    data=df.all.sxt1,
                    family=betar(link='logit'))
summary(mod4)

mod2=gam(mean_freq ~ s(ICP1,by=ICE.by.mut)+
         s(AZI,by=mutation_type)+ 
         s(CIP,by=mutation_type)+
         DOX*mutation_type,
       data=df.all.sxt1,
       family=betar(link='logit'))
summary(mod2)

mod3=gam(mean_freq ~ s(ICP1,by=ICE.by.mut)+
         s(AZI,by=mutation_type)+ 
         s(CIP,by=mutation_type)+
         DOX*mutation_type+
         te(ICP1,AZI,by=mutation_type)+
         te(ICP1,CIP,by=mutation_type)+
         DOX:ICP1:mutation_type,
       data=df.all.sxt1,
       family=betar(link='logit'))
summary(mod3)



AICtab(mod4,mod2,mod3,mod1)

gam.check(mod1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)

### plot
pdf('mod1_plot.pdf',width = 4.7,height = 3,pointsize = 8)
par(mfrow=c(1,3),mar = c(4,4,1,0), cex = 0.7)
itsadug::plot_smooth(mod1, view='ICP1', plot_all='ICE.by.mut',cond=list(ICE.by.mut=c('N.No_ICE','N.ind6','N.ind5')), 
                     cex.main=0.8,rug=F,
                     col=c( "#a739caff", "#858279ff", "#afbeafff"),
                     transform=plogis,ylim=c(0,0.5),
                     se=0,lwd=3)
itsadug::plot_smooth(mod1, view='ICP1', plot_all='ICE.by.mut',cond=list(ICE.by.mut=c('S.No_ICE','S.ind6','S.ind5')), 
                     cex.main=0.8,rug=F,
                     col=c( "#a739caff", "#858279ff", "#afbeafff"),
                     transform=plogis,ylim=c(0,0.5),
                     se=0,lwd=3)
itsadug::plot_smooth(mod1, view='ICP1', plot_all='ICE.by.mut',cond=list(ICE.by.mut=c('I.No_ICE','I.ind6','I.ind5')), 
                     cex.main=0.8,rug=F,
                     col=c( "#a739caff", "#858279ff", "#afbeafff"),
                     transform=plogis,ylim=c(0,0.5),
                     se=0,lwd=3)
dev.off()

