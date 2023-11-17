## GAM with qPCR Vc as a function of qPCR ICP1 and antibiotices
###  this script produces table S7 and figS10

library(readxl)
library(mgcv)
library(performance)


###qPCR data
ICP_VC_qpcr=read_excel('~/mHDM_Aggregate_9_12_23_Draft.xlsx',sheet='mHDM_VC_qPCR') 
df_ratio=ICP_VC_qpcr[,c('Sample_ID','qPCR_ICP1_PFUml_CT28','qPCR_tcpA_CFUml_CT28','qPCR_Ratio_ICP1_VC_CT28')]

colnames(df_ratio)[1]='Sample'

df.all=merge(df_ratio,metadata_all,by='Sample')
df.all$CIP=as.numeric(df.all$CIP)

##remove antibiotics outliers
df.all.out=df.all[!df.all$Sample %in% c('1331650','1394550','D17181815','1319650','D02182446','D17181817','D05181982'), ]

### replace NEG values by 0
df.all.out[df.all.out$qPCR_tcpA_CFUml_CT28=='NEG','qPCR_tcpA_CFUml_CT28'] = '0'
df.all.out$qPCR_tcpA_CFUml_CT28=as.numeric(df.all.out$qPCR_tcpA_CFUml_CT28)


df.all.out[df.all.out$qPCR_ICP1_PFUml_CT28=='NEG','qPCR_ICP1_PFUml_CT28'] = '0'
df.all.out$qPCR_ICP1_PFUml_CT28=as.numeric(df.all.out$qPCR_ICP1_PFUml_CT28)

 
df.all.out$log_qpcr_vc=log10(1+(df.all.out$qPCR_tcpA_CFUml_CT28))
df.all.out$Dehydration_Status=as.factor(df.all.out$Dehydration_Status)


gam1.log=gam(log_qpcr_vc ~ s(qPCR_ICP1_PFUml_CT28) + s(AZI) + te(qPCR_ICP1_PFUml_CT28, AZI) + 
           s(Dehydration_Status, bs = "re"),
         data=df.all.out)

summary(gam1.log)


gam2=gam(log_qpcr_vc ~ s(qPCR_ICP1_PFUml_CT28) + s(AZI) + s(Dehydration_Status, bs = "re"),
         data=df.all.out)

gam3=gam(log_qpcr_vc ~ s(qPCR_ICP1_PFUml_CT28) + s(AZI) + s(CIP) + s(DOX) + 
           te(qPCR_ICP1_PFUml_CT28, AZI) + 
           te(qPCR_ICP1_PFUml_CT28, CIP) + 
           te(qPCR_ICP1_PFUml_CT28, DOX) +
           s(Dehydration_Status, bs = "re"),
         data=df.all.out)

gam4 = gam(log_qpcr_vc ~ s(qPCR_ICP1_PFUml_CT28) + s(AZI) + s(CIP) + 
           s(DOX) + s(Dehydration_Status, bs = "re"),
         data=df.all.out)

gam5 = gam(log_qpcr_vc ~ s(qPCR_ICP1_PFUml_CT28)+ s(DOX) + te(qPCR_ICP1_PFUml_CT28, DOX) + 
             s(Dehydration_Status, bs = "re"),
           data=df.all.out)


gam6= gam(log_qpcr_vc ~ s(qPCR_ICP1_PFUml_CT28)+ s(CIP)+ te(qPCR_ICP1_PFUml_CT28, CIP) +
            s(Dehydration_Status, bs = "re"),
          
          data=df.all.out)

gam7= gam(log_qpcr_vc ~ s(qPCR_ICP1_PFUml_CT28)+ s(AZI) + 
            te(qPCR_ICP1_PFUml_CT28, AZI, by = Dehydration_Status),
          data=df.all.out)


gam8= gam(log_qpcr_vc ~ s(qPCR_ICP1_PFUml_CT28) + s(AZI) + te(qPCR_ICP1_PFUml_CT28, AZI),
          data=df.all.out)


AICtab(gam1.log,gam2,gam3,gam4,gam5,gam6,gam7,gam8)
# dAIC df  
# gam3      0.0 10.1
# gam4      0.4 10.4
# gam1.log  0.9 7.2 
# gam2      1.3 7.6 
# gam7      5.3 14.2
# gam5     11.8 7.8 
# gam6     13.0 7.4 
# gam8     30.6 5.5 

save(gam3,file='gam3.RData')


##Plot
gam.check(gam3,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)


pdf('FigS10.pdf',width = 3.5,height = 3,pointsize = 8)
par(mfrow=c(1,2),mar = c(4,4,1,0), cex = 0.7)

plot_smooth(gam3, view="qPCR_ICP1_PFUml_CT28",  
            plot_all='AZI',
            cond=list(AZI=c(0 ,10 ,30,60,70)),#10,20,30,40,50,60,
            lwd=2.5,sim.ci=T,shade=T,se=F,
            col=c('black','#5BBCD6','#00A08A','#FF0000','darkorchid4'),
            legend_plot_all = T,ann=F,bty='l',
            print.summary = F,rug = F,hide.label = T) #,transfo
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='qPCR Vc (log10)', cex.lab=1,line=2.5)
title(xlab='qPCR ICP1', cex.lab=1,line=2.5)

plot_smooth(gam3, view="AZI",  
            #plot_all='AZI',
            #cond=list(AZI=c(0 ,10 ,30,60,70)),#10,20,30,40,50,60,
            lwd=2.5,sim.ci=T,shade=T,se=F,
            col=c('black','#5BBCD6','#00A08A','#FF0000','darkorchid4'),
            legend_plot_all = T,ann=F,bty='l',
            print.summary = F,rug = F,hide.label = T) #,transfo

axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='qPCR Vc (log10)', cex.lab=1,line=2.5)
title(xlab='AZI', cex.lab=1,line=2.5)
dev.off()
