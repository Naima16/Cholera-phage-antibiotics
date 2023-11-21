### Figure S5 in supplementary figures 
### 

library(phyloseq)
library(vegan)
library (dplyr)
library(ggrepel)


biomfilename1 ='~/bracken_v1.biom'
data1 <- import_biom(biomfilename1, parseFunction=parse_taxonomy_default)

##change rank names from rank1,.. 
colnames(tax_table(data1)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                                o = "Order", f = "Family", g = "Genus", s = "Species")

OTU1=otu_table(data1)
for(i in 1:dim(OTU1)[2])
  colnames(OTU1)[i]=unlist(strsplit(as.character(colnames(OTU1)[i]),"[_]"))[1]

TAX1 = tax_table(data1)
physeq1 = phyloseq(OTU1, TAX1)

physeq1.comp <- microbiome::transform(physeq1, 'compositional')
d1=otu_table(physeq1.comp)

df_myseq1=as.data.frame(t(d1))


####data 2
biomfilename2 ='~/bracken_kraken_v2.biom'
data2 <- import_biom(biomfilename2, parseFunction=parse_taxonomy_default)

##change rank names from rank1,.. 
colnames(tax_table(data2)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                                o = "Order", f = "Family", g = "Genus", s = "Species")

for(i in 1:dim(OTU2)[2])
  colnames(OTU2)[i]=unlist(strsplit(as.character(colnames(OTU2)[i]),"[_]"))[1]

TAX2 = tax_table(data2)
physeq2 = phyloseq(OTU2, TAX2)

physeq2.comp <- microbiome::transform(physeq2, 'compositional')
d2=otu_table(physeq2.comp)

df_myseq2=as.data.frame(t(d2))

df_myseq2$Sample=rownames(df_myseq2)
df_myseq1$Sample=rownames(df_myseq1)

################ 

merged_physeq=merge_phyloseq(physeq1, physeq2)
ps0.comp <- microbiome::transform(merged_physeq, 'compositional')

d3=otu_table(ps0.comp)
df.3=as.data.frame(t(d3))

#remove positive control and human reads
df.4=df.3[rownames(df.3)!="Positivecontrol3",colnames(df.3) != '9606']

## hellinger transform species composition

physeq.hel=decostand(df.4,method='hellinger')
pca_res2=rda(physeq.hel)
sc=scores(pca_res2, display = "sp", scaling = 0)

sc_cp=sc
for (i in 1:nrow(sc_cp))
{
  sc_cp[i,'PC1']=abs( sc_cp[i,'PC1'])
  sc_cp[i,'PC2']=abs( sc_cp[i,'PC2'])
}

dim_1=sc_cp[order(-sc_cp[,1]),1,drop=F]
dim_2=sc_cp[order(-sc_cp[,2]),2,drop=F]

colnames(dim_1)[1]='score'
colnames(dim_2)[1]='score'

high_scor_sp=unique(c(rownames(dim_1)[1:20],rownames(dim_2)[1:20]))

################# PCA to determine dominant species
#############. 
high_scor_sp2=c(high_scor_sp,'979525','979535','979533'). ## add phages
high_scor_sp_tax=merge(as.data.frame(high_scor_sp2),TAX1[rownames(TAX1)%in% high_scor_sp2],by.x='high_scor_sp2',by.y='row.names')

## high_scor_sp_tax
# "Vibrio.cholerae"                           "Escherichia.coli"                         
# [3] "Escherichia.albertii"                      "Klebsiella.pneumoniae"                    
# [5] "Shigella.dysenteriae"                      "Shigella.flexneri"                        
# [7] "Shigella.boydii"                           "Salmonella.enterica"                      
# [9] "Brachyspira.pilosicoli"                    "Faecalibacterium.prausnitzii"             
# [11] "Roseburia.Lachnospiraceae bacterium GAM79" "Roseburia.[Eubacterium] rectale"          
# [13] "Roseburia.intestinalis"                    "Streptococcus.pasteurianus"               
# [15] "Streptococcus.gallolyticus"                "Streptococcus.equinus"                    
# [17] "Streptococcus.mitis"                       "Streptococcus.infantarius"                
# [19] "Streptococcus.lutetiensis"                 "Lactobacillus.ruminis"                    
# [21] "Enterococcus.faecium"                      "Bifidobacterium.longum"                   
# [23] "Bifidobacterium.breve"                     "Bifidobacterium.kashiwanohense"           
# [25] "Collinsella.aerofaciens"                   "Fusobacterium.mortiferum"                 
# [27] "Bacteroides.thetaiotaomicron"              "Bacteroides.vulgatus"                     
# [29] "Bacteroides.fragilis"                      "Prevotella.dentalis"                      
# [31] "Prevotella.intermedia"                     "Prevotella.ruminicola"                    
# [33] "ICP3"                                      "Bifidobacterium.pseudocatenulatum"        
# [35] "Bifidobacterium.adolescentis"              "ICP1"                                     
# [37] "ICP2" 

## only dominant species 34 +  3 phages ICP1, ICP2. and ICP3
# 979525 = icp1
# 979535 = icp3
# 979533 =icp2
# 666 = Vc

df.all.subset=df.4[,colnames(df.4)%in% c(high_scor_sp,'979525','979535','979533')]  
### save(df.all.subset,file='df.all.subset.RData')
dim(df.all.subset)  ##344,37
write.csv(df.all.subset,file='df.all.subset.csv',quote=F,row.names = F)

##change species name
for (i in 1:dim(df.all.subset)[2])
{
  nom_col=colnames(df.all.subset)[i]
  sp=high_scor_sp_tax[high_scor_sp_tax$high_scor_sp2==nom_col,'Species']
  gen=high_scor_sp_tax[high_scor_sp_tax$high_scor_sp2==nom_col,'Genus']
  new_col=paste(unlist(strsplit(as.character(gen),"[__]"))[3],unlist(strsplit(as.character(sp),"[__]"))[3],sep='.')
  new_col
  colnames(df.all.subset)[i]=new_col
}
df.all.subset=as.data.frame(df.all.subset)

names(df.all.subset)[names(df.all.subset)=='Cepunavirus.Vibrio phage ICP2'] <- 'ICP2'
names(df.all.subset)[names(df.all.subset)=='NA.Vibrio phage ICP1'] <- 'ICP1'
names(df.all.subset)[names(df.all.subset)=='Teseptimavirus.Vibrio phage ICP3'] <- 'ICP3'

### define VC_phage_abund factor to color samples based on their Vc-phage abundance
names(df.4)[names(df.4)=='Cepunavirus.Vibrio phage ICP2'] <- 'ICP2'
names(df.4)[names(df.4)=='NA.Vibrio phage ICP1'] <- 'ICP1'
names(df.4)[names(df.4)=='Teseptimavirus.Vibrio phage ICP3'] <- 'ICP3'

names(df.4)[names(df.4)=='666'] <- 'Vc'
names(df.4)[names(df.4)=='979525'] <- 'ICP1'

df.4=df.4%>%
  mutate(VC_phage_abund = case_when(
    ICP1>0.001 & Vibrio.cholerae>0.005 ~ 'VCplus_ICP1plus',
    ICP1>0.001 & Vibrio.cholerae<=0.005 ~ 'VCmin_ICP1plus',
    ICP1<=0.001 & Vibrio.cholerae>0.005 ~ 'VCplus_ICP1min',
    ICP1<=0.001 & Vibrio.cholerae<=0.005 ~ 'VCmin_ICP1min'
  ))
##

##plot pca
ii=summary(pca_res2)  
sp=as.data.frame(ii$species[,1:2])*2
st=as.data.frame(ii$sites[,1:2])  
sp_dominant_pc=sp[rownames(sp)%in%high_scor_sp2,]

for (i in 1:dim(sp_dominant_pc)[1])
{
  nom=rownames(sp_dominant_pc)[i]
  if(any(row.names(TAX1) == nom))
  {
    sp=TAX1[rownames(TAX1)==nom,'Species']
    gen=TAX1[rownames(TAX1)==nom,'Genus']
    new_col=paste(unlist(strsplit(as.character(gen),"[__]"))[3],unlist(strsplit(as.character(sp),"[__]"))[3],sep='.')
    new_col
    rownames(sp_dominant_pc)[i]=new_col
  } else
  {
    sp=TAX1[rownames(TAX2)==nom,'Species']
    gen=TAX1[rownames(TAX2)==nom,'Genus']
    new_col=paste(unlist(strsplit(as.character(gen),"[__]"))[3],unlist(strsplit(as.character(sp),"[__]"))[3],sep='.')
    new_col
    rownames(sp_dominant_pc)[i]=new_col
  }
}

rownames(sp_dominant_pc)[rownames(sp_dominant_pc)=='Cepunavirus.Vibrio phage ICP2'] <- 'ICP2'
rownames(sp_dominant_pc)[rownames(sp_dominant_pc)=='NA.Vibrio phage ICP1'] <- 'ICP1'
rownames(sp_dominant_pc)[rownames(sp_dominant_pc)=='Teseptimavirus.Vibrio phage ICP3'] <- 'ICP3'

options(ggrepel.max.overlaps = Inf)

p=ggplot() +
  geom_point(data = st,aes(PC1,PC2,color=df.4$VC_phage_abund),size=1.5,alpha=0.7) +
  scale_color_manual(values = c( "#F2AD00", "#FF0000", "#00A08A","darkorchid1"))+
  scale_shape_manual(values = c(3:17))
p=p+geom_segment(data = sp_dominant_pc,aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(angle=10,length = unit(0.2,"cm"),
                               type = "closed"),linetype=1, size=0.25,colour = "#00A08A",alpha=0.3)  
p=p+geom_text_repel(data = sp_dominant_pc,aes(PC1,PC2,label=row.names(sp_dominant_pc)),
                  size=2.5,color='#000078ff',segment.color = 'grey',fontface ='italic') 
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
        legend.position = 'bottom')+
  labs(colour="VC+ (metagenomic)")
p
ggsave('figS5.pdf',height = 7,width = 7)
