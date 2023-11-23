
## script for S9
## NM nov 2023

library(FactoMineR)
library(factoextra)
library(tidyr)
library(dplyr)
library(cowplot)

#### ARGs from deepARG
all.arg=read.csv('/Users/naimamadi/cholera_sept_2023/data_tables/deeparg/merge_quant_subtype_deepARG.csv',check.names = F)

names(all.arg)[names(all.arg)=="16s-NormalizedReadCount"] <- "NormalizedReadCount_16s"
names(all.arg)[names(all.arg)=="ARG-group"] <- "ARG_group"

all.arg=all.arg[,colnames(all.arg) != 'ReadCount']

ARG_df <- spread(all.arg, ARG_group, NormalizedReadCount_16s)

ARG_df[is.na(ARG_df)] <- 0
rownames(ARG_df)=ARG_df$Sample

ARG_df=ARG_df[,colnames(ARG_df) != 'Sample']

### species composition for the most dominant species
df.all.subset=read.csv('~/df.all.subset.csv')
rownames(df.all.subset)=df.all.subset$Sample
df.all.subset=df.all.subset[,colnames(df.all.subset) != 'Sample']

mfa.df=merge(df.all.subset,ARG_df,by='row.names')
rownames(mfa.df)=mfa.df$Row.names
mfa.df=mfa.df[,colnames(mfa.df) !='Row.names']

mfa.df=mfa.df  %>%
  relocate(ICP3, .after = ICP1)

##MFA with bacteri and ARGs
mfa.df_withoutphages=mfa.df[,!colnames(mfa.df) %in% c('ICP1','ICP2','ICP3')]

mfa.out=FactoMineR::MFA(mfa.df_withoutphages,group=c(34,634),type=c('s','s'),name.group=c('Bacteria',"ARG"))
#####
summary(mfa.out)

tab_contrib=get_mfa_quanti_var(mfa.out)$contrib[,c(1,2)]
tab.contrib=tab_contrib[1:20,]

options(ggrepel.max.overlaps = 60)

pdf('FigS9.pdf',width = 7,height = 7)
fviz_mfa_var(mfa.out, "quanti.var", 
             palette = c( "darkorchid4", "turquoise4"), 
             col.var.sup = "violet", repel = TRUE,
             labelsize = 2,alpha.var='contrib',select.var = list(contrib=200))
dev.off()


## contribution to each axis
ax1=as_grob(fviz_contrib(mfa.out, choice = "quanti.var", axes = 1, top = 20,
                         palette = "jco")  +
              theme_minimal() +
              theme(axis.text = element_text(size=5,angle=45,hjust=1),
                    legend.position='none',
                    axis.title=element_text(size=7),
                    axis.title.x=element_blank(),
                    text = element_text(size = 7)))
ax2=as_grob(fviz_contrib(mfa.out, choice = "quanti.var", axes = 2, top = 20,
                         palette = "jco") +
              theme_minimal() +
              theme(axis.text = element_text(size=5,angle=45,hjust=1),
                    legend.position='right',
                    axis.title=element_text(size=7),
                    axis.title.x=element_blank(),
                    text = element_text(size = 7)))
 
pdf('Contrib_20Top.pdf',width = 6,height = 3,pointsize = 0.5)
plot_grid(ax1,ax2,ncol=2,nrow=1,rel_widths =c(1,1.3))
dev.off()




