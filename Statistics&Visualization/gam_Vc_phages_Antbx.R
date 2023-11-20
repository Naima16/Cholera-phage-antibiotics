
### script for fig2 and table S6
####
library(purrr)
library(mgcv)
library(performance)
library(bbmle)
library(itsadug)

####. rel abundances
load('~/df.species.otu.comp.otu_vc_icp1.RData')

### metadata+antibio+AMR
patient_antibio=read.csv('~/metadata_all.csv')

##########

df.species.otu.comp.otu_vc_icp1$Sample=rownames(df.species.otu.comp.otu_vc_icp1)

df.all=merge(df.species.otu.comp.otu_vc_icp1,patient_antibio,by='Sample')
colnames(df.all)
dim(df.all)

df.all$CIP=as.numeric(df.all$CIP)

## remove antbx outliers 
df.all.out=df.all[!df.all$Sample %in% c('1331650','1394550','D17181815','1319650','D02182446','D17181817','D05181982'), ]


df.gam=df.all.out [,c('Sample','Vc','ICP1','ICP3','ICP2','CIP','AZI','DOX','Dehydration_Status')]

## GAMs with Vc as a function of ICP1+antibiotics (abundances from MG)
library(mgcv)
##df.gam[df.gam$Sample=='1331650',], il n'est pas ds la data !!!!!

gam1 <- gam(Vc ~ s(ICP1)+s(AZI)+ te(ICP1,AZI)+
            s(Dehydration_Status,bs="re"),
            data=df.gam,
            family=betar(link='logit'))

summary(gam1)     
r2(gam1)
gam.check(gam1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)


gam2 <- gam(Vc ~ s(ICP1)+s(AZI)+
                    s(Dehydration_Status,bs="re"),
                    data=df.gam,
                    family=betar(link='logit'))

gam3 <- gam(Vc ~ s(ICP1)+s(AZI)+ s(CIP)+ s(DOX)+
                    te(ICP1,AZI)+
                    te(ICP1,CIP)+
                    te(ICP1,DOX)+
                    s(Dehydration_Status,bs="re"),
                    data=df.gam,
                    family=betar(link='logit'))

gam4 <- gam(Vc ~ s(ICP1)+s(AZI)+s(CIP)+s(DOX)+
                    s(Dehydration_Status,bs="re"),
            data=df.gam,
            family=betar(link='logit'))

gam5 <- gam(Vc ~ s(ICP1)+s(DOX)+te(ICP1,DOX)+
                    s(Dehydration_Status,bs="re"),
            data=df.gam,
            family=betar(link='logit'))

gam6 <- gam(Vc ~ s(ICP1)+s(CIP)+te(ICP1,CIP)+
                    s(Dehydration_Status,bs="re"),
            data=df.gam,
            family=betar(link='logit'))
gam7 <- gam(Vc ~ s(ICP1)+s(AZI)+te(ICP1,AZI)+
                    s(Dehydration_Status,bs="re"),
            data=df.gam,
            family=betar(link='logit'))

gam8 <- gam(Vc ~ s(ICP1)+s(AZI)+ te(ICP1,AZI),
            data=df.gam,
            family=betar(link='logit'))


AICtab(gam1,gam2,gam3,gam4,gam5,gam6,gam7,gam8)

save(gam1,file='gam1.RData')
r2(gam1) 


###plot gam1 

library(itsadug)

pdf('fig2.pdf',width = 3.5,height = 2.7,pointsize = 8)
par(mfrow=c(1,2),mar = c(4,4,1,0), cex = 0.7)

plot_smooth(gam1, view="AZI", #10,20,30,40,50,60,
            col='black',lwd=2.3,sim.ci=T,shade=T,se=F,ylim=c(0,0.4),
            legend_plot_all = T,ann=F,bty='l',print.summary = F,rug = F,hide.label = T,
            transform=plogis) #,transform=exp)#,transform=plogis)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='Vc', cex.lab=1.3,line=2.3)
title(xlab='AZI', cex.lab=1.3,line=2.3)


plot_smooth(gam1, view="ICP1", plot_all='AZI',cond=list(AZI=c(0 ,10 ,30, 50, 70)), #10,20,30,40,50,60,
            col=c('#4b0000ff','#970400ff',"#f20400ff",'#f28d00ff','#f2c500ff'),
           ### col=c('#FF0000','darkorchid4','#00A08A','#F2AD00','#5BBCD6'),
            lwd=2.3,sim.ci=T,shade=T,se=F,
            legend_plot_all = T,ann=F,bty='l',print.summary = F,rug = F,
            hide.label = T,transform=plogis,ylim=c(0,0.4)) #,transform=exp)#,transform=plogis)
title(main='AZI', cex.main = 1,line=0)

axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='Vc', cex.lab=1.3,line=2.3)
title(xlab='ICP1', cex.lab=1.3,line=2.3)


dev.off()


