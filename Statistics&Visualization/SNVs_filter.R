
## SNVs filtering
## create snv.filter3.biallelic.csv based on inStrain output
## NM nov 2023

library(plyr)

### inStrain output
snv.all=read.csv('~/snv.all.csv',header = T)

for(i in 1:dim(snv.all)[1])
{
  nn=unlist(strsplit(as.character(snv.all[i,'scaffold']),"_"))[1]
  if (nn != 'phage') {
    snv.all[i,'genome']='1458300'
  } else {
    snv.all[i,'genome']=unlist(strsplit(as.character(snv.all[i,'scaffold']),"_"))[2]
  }
}

VC.snv.NoFilt=snv.all[snv.all$genome=='1458300',]
ICP1.snv.Nofiltr=snv.all[snv.all$genome=='ICP1',]
colnames(VC.snv.NoFilt)[20]='Sample' 

#### filter out covergae 
## coverage profile from inStrain
cov.all=read.csv('~/cov.all.csv',header = T)

sample.median <- ddply( VC.snv.NoFilt, .(Sample), function(x) median(x$position_coverage) )
colnames(sample.median)=c('Sample','median_cov')

##scaffold length filter
for (i in 1:dim(VC.snv.NoFilt)[1])
{
  scaf_name=VC.snv.NoFilt[i,'scaffold']
  scaf_length=unique(cov.all[cov.all$scaffold==scaf_name,'length'])
  if(VC.snv.NoFilt[i,'position'] >100 & VC.snv.NoFilt[i,'position']< scaf_length-100) {
    VC.snv.NoFilt[i,'scaff_filter']='yes'
  } else {
    VC.snv.NoFilt[i,'scaff_filter']='no'
  }  
}

## nb of sites with cov filter and scaffold filter (avec python:58563)
for (i in 1:dim(VC.snv.NoFilt)[1])
{
  samp=VC.snv.NoFilt[i,'Sample']
  med.sample=sample.median[sample.median$Sample==samp,'median_cov']
  inf_limit=0.3*med.sample
  sup_limit=3*med.sample
  if (VC.snv.NoFilt[i,'position_coverage'] > inf_limit & VC.snv.NoFilt[i,'position_coverage'] < sup_limit){
    VC.snv.NoFilt[i,'cov_filter']='yes'
  }  else   {VC.snv.NoFilt[i,'cov_filter']='no'}
}

## filter out position with cov < 20
VC.snv.NoFilt.20=VC.snv.NoFilt[VC.snv.NoFilt$position_coverage>=20,]

## nb of sites with cov filter and scaffold filter 
snv.vc.filt=VC.snv.NoFilt.20[VC.snv.NoFilt.20$scaff_filter=='yes' & VC.snv.NoFilt.20$cov_filter=='yes',]
dim(snv.vc.filt)  #107313     23
length(unique(snv.vc.filt$Sample)) ##198

### filter out sites that did not pass coverage filter in more than 2 samples
sites.cov.not.passed=ddply(VC.snv.NoFilt.20, .(scaffold,position), summarize, samples_not_ok = sum(cov_filter=='no'))
ss=sites.cov.not.passed[sites.cov.not.passed$samples_not_ok>2,]
snv.filter1=merge(snv.vc.filt,sites.cov.not.passed,by=c('scaffold','position'))
snv.filter2=snv.filter1[snv.filter1$samples_not_ok<2,]

#### filter out positions (scaff+position) from the positive control

t1=snv.filter2[snv.filter2$Sample=='Positivecontrol3',] 
posControl_filtered=VC.snv.NoFilt[VC.snv.NoFilt$Sample=='Positivecontrol3',]  
dim(posControl_filtered)

positiv_scaffold=unique(posControl_filtered$scaffold)
length(positiv_scaffold)  

snv.filter2$scaff_pos <- paste(snv.filter2$scaffold, snv.filter2$position)

for (i in 1:dim(snv.filter2)[1])
{
  if(snv.filter2[i,'scaff_pos'] %in% posControl_filtered$scaff_pos)  
    snv.filter2[i,'control3Filter']='no' 
  
  else 
    snv.filter2[i,'control3Filter']='yes'
}

snv.filter3.biallelic=snv.filter3[snv.filter2$control3Filter=='yes' & snv.filter2$allele_count==2,]
save(snv.filter3.biallelic,file='snv.filter3.biallelic.RData')
