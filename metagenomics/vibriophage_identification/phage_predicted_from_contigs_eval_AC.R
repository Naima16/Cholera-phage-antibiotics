# November 2023
# Aline Cuénod

# load packages
library('tidyr')
library('dplyr')
library('ggplot2')
library('data.table')
library('patchwork')
library('scales')
library('readr')

setwd('./01_data/data')

# import virus output 
genomad_virus <- read.table('genomad/genomad_virus_all.csv', sep = '\t', header = T)
# remove positive control
genomad_virus <- genomad_virus[!grepl('control', genomad_virus$seq_name),]

# number of viral contigs
length(unique(genomad_virus$seq_name))

# add sample and virus family
genomad_virus['sample'] <- gsub('\\_.*', '',genomad_virus$file)
genomad_virus['virus_family'] <- gsub('.*\\;', '',genomad_virus$taxonomy)
# add whether the 'virus' contig was assigned a prophage or not
genomad_virus['provirus_or_not'] <- ifelse(genomad_virus$topology == 'Provirus', 'Predicted Pro-phage', 'Predicted Phage')

#import iphop output
iphop <- read.table('iphop/iphop_all.csv', sep = '\t', header = T)
# remove control sample
iphop <- iphop[!grepl('control',iphop$Virus),] 
# the iphop output includes all prediction with a score above 90%, this can be more than one per sample. Therefore, supset to only one per sample (the one with the highest score) 
iphop <- iphop %>% group_by(Virus) %>% slice_max(Confidence.score)
# add host genus assignment
iphop['host_genus_short'] <- gsub('(.*\\;g\\_\\_)', '', iphop$Host.genus)
iphop['vibrio_host'] <- ifelse(grepl('Vibrio',iphop$host_genus_short), 'Vibrio', 'Other')
# some have ties as score, check these
iphop_check <- iphop[duplicated(iphop$Virus),]
table(iphop_check$vibrio_host) # all have host 'Other' (no Vibrio). These do not matter for out analysis (as we focus on Vibrio phages), so just choose one
# select one sample with highest score per sample
iphop <- iphop %>% group_by(Virus) %>% slice_max(Confidence.score, with_ties = FALSE)

# merge to genomad output (some do not have an iphop output, iphop only outputs predictions with a score > 90%)
genomad_virus <- merge(genomad_virus, iphop, by.x = c('seq_name', 'sample'), by.y = c('Virus','sample'), all = T)
# color by if vibrio host assigned
basic_virus_occurence_3 <- ggplot(genomad_virus, aes(x=sample, fill = vibrio_host)) + 
  geom_bar(stat="count") + facet_grid(rows = vars(provirus_or_not), scales = "free_y")  + ylab('Nr. of unique viral contigs') + scale_fill_manual(values = c('skyblue', 'red'), na.value = 'lightgrey') + theme_light() + theme(axis.text.x=element_blank())

# For further analysis export all contigs which iphop assigned as 'Vibrio phages'
genomad_virus_vibrio <- genomad_virus[genomad_virus$host_genus_short == 'Vibrio',]
genomad_virus_vibrio <- genomad_virus_vibrio[!is.na(genomad_virus_vibrio$seq_name),]
# count the number of 'Vibrio' phage contigs
length(unique(genomad_virus_vibrio$seq_name))
# count the number of samples
length(unique(genomad_virus_vibrio$sample))

# I blasted (blastn) all the contigs which where identified as 'Vibrio' phages against all Vibrio phages included in the INPHARED database and the PLE and TLR satellitephages
# I realised that the phages in the inphared database are not named consistently (not all ICP1 are named ICP1). Therefore, I included all inphared genomes (and PLE and TLC) to run in vContact2, which includes reference phages (RefSeq)
blastout_phages <- read.table('blastn/blastout_all.csv', sep = '\t', header = T)
blastout_phages_all <- blastout_phages
# check why PLE does not come up
blastout_phages_all['PLE'] <- ifelse(grepl('KC152960', 'KC152961', 'MF176135',  blastout_phages_all$sseqid), 'PLE', 'no')
table(blastout_phages_all['PLE']) # I think these might not be in because they where maybe missed by geNomad (?)

# shorten query name
blastout_phages['sample'] <- gsub('\\_.*', '', blastout_phages$query)

# filter for 80% alignment length and 90% coverage
blastout_phages <- blastout_phages[blastout_phages$pident >= 80, ]
#blastout_phages <- blastout_phages[blastout_phages$length >= (blastout_phages$qlen*0.6),]
blastout_phages <- blastout_phages[blastout_phages$qcovs >= 90,]
# consider only best blast hit per query (there can be duplicates in the db, eg multiple CTXphi)
# sort by bitscore and only consider highest
blastout_phages <- blastout_phages %>% group_by(qseqid) %>% arrange(desc(bitscore)) %>% slice_head(n=1)
blastout_phages['inphared_accession'] <- gsub('\\..*', '', blastout_phages$sseqid)
# add meta-data (phage names) for db entries
inphred_meta <- read.csv2('blastn/1Oct2023_data_excluding_refseq.tsv', sep = '\t')
setdiff(blastout_phages$inphared_accession, inphred_meta$Accession) # this is the TLC reference its not included in the inphared database
# subset to data of interest (only inphared entries which are matching one of our viral contigs)
inphred_meta <- inphred_meta[inphred_meta$Accession %in% blastout_phages$inphared_accession, c('Accession','Description','Classification','Genome.Length..bp.')]
inphred_meta['Description'] <- gsub('Vibrio phage ICP1_2012_A', 'Vibrio phage ICP1', inphred_meta$Description)
# add tag, so that its clear that these information relate to the inphared db entries
colnames(inphred_meta) <- paste0('blasthit_', colnames(inphred_meta)) # add this to ake clear that the information refers to the blasthit and not to the query
# merge
blastout_phages <- merge(blastout_phages, inphred_meta, by.x = 'inphared_accession', by.y = "blasthit_Accession", all.x = T)
colnames(blastout_phages) <- paste0('bl', colnames(blastout_phages))

# plot the contigs where iphop identified a vibrio host
setdiff(blastout_phages$blqseqid, genomad_virus_vibrio$seq_name)
genomad_virus_vibrio_b <- merge(genomad_virus_vibrio, blastout_phages, by.x = 'seq_name', by.y = 'blqseqid', all.x = T)

overview <- ggplot(genomad_virus_vibrio_b) +
  geom_bar(aes(x=provirus_or_not, fill = blblasthit_Description), stat="count") 
# clarify the TLC description. If no blasthit description available, no blast hit was found
genomad_virus_vibrio_b$blblasthit_Description <- ifelse(grepl('contig', genomad_virus_vibrio_b$blsseqid), 'Vibrio satellite phage TLC', 
                                                              ifelse(!is.na(genomad_virus_vibrio_b$blblasthit_Description), genomad_virus_vibrio_b$blblasthit_Description, 'No blast hit'))
# sort by frequency
genomad_virus_vibrio_b$blblasthit_Description <- factor(genomad_virus_vibrio_b$blblasthit_Description, levels = names(sort(table(genomad_virus_vibrio_b$blblasthit_Description), decreasing = TRUE)))

blasthits <- ggplot(genomad_virus_vibrio_b) +
  geom_bar(aes(blblasthit_Description, fill = 'red'), stat="count") +
  facet_grid(rows = vars(provirus_or_not)) + 
  ylab('Number of viral contigs') + 
  xlab('descriprion blasthit') +
  theme_light() +   theme(axis.text.x=element_text(angle=60,hjust=1), legend.position = 'none') 

# remove duplicates, so that I can plot number of samples in which at least one contig was assigned to
genomad_virus_vibrio_b_per_sample <- genomad_virus_vibrio_b %>% group_by(sample, blblasthit_Description) %>% slice_max(Confidence.score, with_ties = FALSE)

blasthits_per_sample <- ggplot(genomad_virus_vibrio_b_per_sample) +
  geom_bar(aes(blblasthit_Description, fill = 'red'), stat="count") +
  facet_grid(rows = vars(provirus_or_not)) + 
  ylab('Number of samples with at least one contig matching') + 
  xlab('descriprion blasthit') +
  theme_light() +   theme(axis.text.x=element_text(angle=60,hjust=1), legend.position = 'none') 
blasthits_per_sample

# I realised that the phages in the inphared database are not named consistently (not all ICP1 are named ICP1). Therefore, I included all inphared and all our identified Vc phage genomes to run in vContact2
# import inphared metadata
inphred_meta_all <- read.csv2('blastn/1Oct2023_data_excluding_refseq.tsv', sep = '\t')
inphred_meta_Vibrio <- inphred_meta_all[inphred_meta_all$Host == 'Vibrio',]

# import vContact2 output  of the inphared data (I also run our metagenome contigs together with the inphared references in vContact2 which gives a somehow unexpercted results, I think mainly because the viral contigs I put in are not necessarily complete viruses, so it eg. clusters all metagenome sequences which have 99% pident to ICP1 in another VC than OCP itself, I think because of incompleteness)
vContact2_inphared <- read.csv('vContact2/genome_by_genome_overview.csv')

# add column stating VC of well known Vibrio phages
vContact2_inphared['VC_no_subcluster'] <- gsub('(VC\\_\\d+)(\\_\\d+$)', '\\1',vContact2_inphared$VC) #here: https://bitbucket.org/MAVERICLab/vcontact2/src/master/ it states that everything after the '_' refers to the subcluster and the VC is defined by the number before this (eg VC_30)
# add a 'desctiption' section to all VC which include a well known reference
vContact2_inphared['VC_description'] <- ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('Vibrio.*ICP1', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_ICP1', 
                                               ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('Vibrio.*ICP2', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_ICP2', 
                                                      ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('Vibrio.*ICP3', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_ICP3', 
                                                             ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('Vibrio.*CTXphi', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_CTXphi', 
                                                                    ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('Vibrio.*K139', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_K139',
                                                                           ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('Vibrio.*12A4', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_12A4',
                                                                                  ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('Vibrio.*fs1', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_fs1', 
                                                                                         ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('KC152960', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_PLE1', 
                                                                                                ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('KC152961', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_PLE2', 
                                                                                                       ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('MF176135', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_PLE3',  
                                                                                                              ifelse(vContact2_inphared$VC_no_subcluster == unique(vContact2_inphared[grepl('_contig', vContact2_inphared$Genome),]$VC_no_subcluster), 'Same_VC_as_TLC', vContact2_inphared$VC_no_subcluster)))))))))))
# add information whether the seqs come from our sequences ('Metagenomes') from the INPHARED database (incl. satellite phages) or from the refseq database included in vContact2
vContact2_inphared['dataset'] <- ifelse(vContact2_inphared$Genome %in% genomad_virus$seq_name, 'Metagenomes', 
                                      ifelse(vContact2_inphared$Genome %in% inphred_meta_Vibrio$Accession, 'inphared', 
                                             ifelse(vContact2_inphared$Genome %in% c('136Vc08_contig_2', 'KC152960', 'KC152961', 'MF176135'), 'Satellites', 'RefSeq')))

# check which VC are assigned to samples which have ctxphi blasthit
vContact2_inphared[vContact2_inphared$Genome %in% blastout_phages[blastout_phages$blblasthit_Description == 'Vibrio phage CTXphi',]$blqseqid,]$VC_no_subcluster
# add this information to the blastn results
genomad_virus_vibrio_b_per_sample_VC <- merge(genomad_virus_vibrio_b_per_sample, vContact2_inphared[,c('Genome', 'VC_description')], by.x='blinphared_accession', by.y='Genome', all.x = TRUE, all.y=FALSE)

# add column which includes info when there was no blasthit 
genomad_virus_vibrio_b_per_sample_VC['VC_description_blast'] <- ifelse(genomad_virus_vibrio_b_per_sample_VC$blblasthit_Description == 'No blast hit', 'No blast hit', genomad_virus_vibrio_b_per_sample_VC$VC_description)
genomad_virus_vibrio_b_per_sample_VC$VC_description_blast <- gsub('Same_VC_as_', '',genomad_virus_vibrio_b_per_sample_VC$VC_description_blast)
genomad_virus_vibrio_b_per_sample_VC$VC_description_blast <- factor(genomad_virus_vibrio_b_per_sample_VC$VC_description_blast, levels= c( "CTXphi","ICP1","ICP3","ICP2","TLC","K139","No blast hit"))

blast_vC_per_sample <- ggplot(genomad_virus_vibrio_b_per_sample_VC) +
  geom_bar(aes(VC_description_blast), fill = 'darkblue', color = 'darkblue', stat="count") +
  facet_grid(rows = vars(provirus_or_not)) + 
  ylab('Number of samples with at least one contig matching') + 
  xlab('VC of blastn match') +
  theme_light() +   theme(axis.text.x=element_text(angle=60,hjust=1), legend.position = 'none') 
blast_vC_per_sample

# count samples per viruses
length(unique(genomad_virus_vibrio_b_per_sample_VC[genomad_virus_vibrio_b_per_sample_VC$VC_description_blast == 'No blast hit', ]$sample))
length(unique(genomad_virus_vibrio_b_per_sample_VC[genomad_virus_vibrio_b_per_sample_VC$VC_description_blast == 'CTXphi', ]$sample))
length(unique(genomad_virus_vibrio_b_per_sample_VC[genomad_virus_vibrio_b_per_sample_VC$VC_description_blast == 'ICP1', ]$sample))
length(unique(genomad_virus_vibrio_b_per_sample_VC[genomad_virus_vibrio_b_per_sample_VC$VC_description_blast == 'ICP2', ]$sample))
length(unique(genomad_virus_vibrio_b_per_sample_VC[genomad_virus_vibrio_b_per_sample_VC$VC_description_blast == 'ICP3', ]$sample))
length(unique(genomad_virus_vibrio_b_per_sample_VC[genomad_virus_vibrio_b_per_sample_VC$VC_description_blast == 'K139', ]$sample))
length(unique(genomad_virus_vibrio_b_per_sample_VC[genomad_virus_vibrio_b_per_sample_VC$VC_description_blast == 'TLC', ]$sample))


#pdf('./blast_vC_per_sample.pdf', height = 4.5, width = 3)
#blast_vC_per_sample
#dev.off()

# look at the ones with no blast hit more closely
# import the network
c1 <- read_table("/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/03_AMR_meta/01_data/09_vContact2/output_inphared_and_own_and_sat/c1.ntw", 
                 col_names = c("source","target","weight"))

c1_subset <- c1[c1$source %in% genomad_virus$seq_name | c1$target %in% genomad_virus$seq_name,]

# check how many of the ones with no blast hit are singletos
check_no_blast_network <- vContact2_inphared[vContact2_inphared$Genome %in% genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit',]$seq_name,]
length(unique(gsub('\\_.*','',check_no_blast_network$Genome)))
table(check_no_blast_network$VC.Status) # 245 Singletons and 22 outlier, check the remaining ones more closely
check_no_blast_network<-check_no_blast_network[!check_no_blast_network$VC.Status %in% c('Outlier','Singleton'),]
length(unique(gsub('\\_.*', '',check_no_blast_network$Genome))) # these are 57 phages of 38 samples which cluster in the network, but where no blasthit was found. check what references are in the cluster
table(check_no_blast_network$VC_no_subcluster)
check_no_blast_clusters <- vContact2_inphared[vContact2_inphared$VC_no_subcluster %in% check_no_blast_network$VC_no_subcluster,]
check_no_blast_clusters <- check_no_blast_clusters[!is.na(check_no_blast_clusters$VC_no_subcluster),]
check_no_blast_clusters <- check_no_blast_clusters[!check_no_blast_clusters$VC_no_subcluster == '',]
table(check_no_blast_clusters$VC_no_subcluster, check_no_blast_clusters$dataset) # all but three of the clusters are metagenome only clusters. Check the three which also include refseq / inphared db seq
table(check_no_blast_network$VC_no_subcluster) # these are the number of no-blast-hit-contigs per VC
check_no_blast_clusters_ref <- check_no_blast_clusters[check_no_blast_clusters$VC_no_subcluster %in% check_no_blast_clusters[check_no_blast_clusters$dataset %in% c('inphared', 'RefSeq', 'Satellites'),]$VC_no_subcluster,]
# VC_no_subcluster_285 are mainly Erwinia phages (ok that there is no blast hit to a vibrio phage)
# there is only one genome from VC3 which is a VC including CTXphi. The blast hit of this particular genome is low qcov (<40%), so its ok its not in.
# eighteen cluster within VC_6 which is the VC which includes the satellite TLC

# check size of metagenome only cluster
check_no_blast_clusters_meta_only_clusters <- check_no_blast_clusters[!(check_no_blast_clusters$VC_no_subcluster %in% check_no_blast_clusters[check_no_blast_clusters$dataset %in% c('inphared','RefSeq', 'Satellites'),]$VC_no_subcluster),]
check_no_blast_clusters_meta_only_clusters_unique <- check_no_blast_clusters_meta_only_clusters[!duplicated(check_no_blast_clusters_meta_only_clusters$VC_no_subcluster),]

mean(check_no_blast_clusters_meta_only_clusters_unique$VC.Size)
range(check_no_blast_clusters_meta_only_clusters_unique$VC.Size)

# add VC of the blasthit first
vContact2_info_blastref <- vContact2_inphared[,c('Genome', 'VC_description')]
colnames(vContact2_info_blastref) <- c('blastref_accession', 'VC_description_blastref')

genomad_virus_vibrio_b_cVo <- merge(genomad_virus_vibrio_b, vContact2_info_blastref, by.x='blinphared_accession', by.y='blastref_accession', all.x = TRUE, all.y=FALSE)

sum_vContact_blast <- merge(genomad_virus_vibrio_b_cVo[,c("seq_name","sample","length",'VC_description_blastref',"blblasthit_Classification","blblasthit_Genome.Length..bp.")], vContact2_inphared, by.x = 'seq_name', by.y = 'Genome')
sum_vContact_blast$VC_description_blastref <- ifelse(is.na(sum_vContact_blast$VC_description_blastref), 'No blast hit', sum_vContact_blast$VC_description_blastref)

# to reduce confusion, remove 'same_VC_as_' from the blasthit desctiption
sum_vContact_blast$VC_description_blastref <- gsub('Same_VC_as_', '', sum_vContact_blast$VC_description_blastref)
sum_vContact_blast$VC_description <- ifelse(grepl('Singleton|Outlier|Overlap', sum_vContact_blast$VC.Status), sum_vContact_blast$VC.Status, sum_vContact_blast$VC_description)
# plot the contigs with and without blasthits
sum_vContact_blast_no_blasthit <- sum_vContact_blast[sum_vContact_blast$VC_description_blastref == 'No blast hit',]
sum_vContact_blast$VC_description <- ifelse(sum_vContact_blast$VC_description %in% check_no_blast_clusters_meta_only_clusters_unique$VC_no_subcluster, 'VC with no References\n(metagenomes from this study only)', sum_vContact_blast$VC_description)

sum(sum_vContact_blast$VC_description_blastref == 'No blast hit')
table(sum_vContact_blast[sum_vContact_blast$VC_description_blastref == 'No blast hit',]$VC_description)
table(sum_vContact_blast[sum_vContact_blast$VC_description_blastref == 'ICP1',]$VC_description)

# Count contigs with which blasthit clustered in which VC
sum_vContact_blast_sum  <- sum_vContact_blast %>%  
  group_by(VC_description_blastref, VC_description) %>%  
  summarize(Count = n()) 

sum_vContact_blast_sum$VC_description_blastref <- ifelse(is.na(sum_vContact_blast_sum$VC_description_blastref), 'No blast hit', sum_vContact_blast_sum$VC_description_blastref)
sum_vContact_blast_sum$VC_description_blastref <- factor(sum_vContact_blast_sum$VC_description_blastref, levels= c( "CTXphi","ICP1","ICP3","ICP2", "TLC","K139","No blast hit"))
sum_vContact_blast_sum$VC_description <- factor(sum_vContact_blast_sum$VC_description, levels = c("Same_VC_as_CTXphi", "VC_4", "Same_VC_as_ICP2", "Same_VC_as_ICP3", "Same_VC_as_TLC","Same_VC_as_K139", "VC with no References\n(metagenomes from this study only)", "VC_285", "Overlap (VC_28/VC_29)", "Singleton", "Outlier"))

blast_vCo_comp <- ggplot(sum_vContact_blast_sum, aes(x=VC_description_blastref, y=Count, fill=VC_description)) + 
  geom_bar(stat='identity', position = position_dodge2(width = 0.9, preserve = "single"), width = 1) + 
  scale_fill_brewer(palette = 'Set3') + 
  theme_light() + theme(legend.position = 'bottom', axis.text.x=element_text(angle=60,hjust=1)) + ylab('Number of contigs') + xlab('BLASTn hit') 

sum_vContact_blast_sum_onlyhits <- sum_vContact_blast_sum[sum_vContact_blast_sum$VC_description_blastref != 'No blast hit',]
sum_vContact_blast_sum_onlyhits$VC_description <- factor(sum_vContact_blast_sum_onlyhits$VC_description, levels = c("Same_VC_as_CTXphi", "VC_4", "Same_VC_as_ICP2", "Same_VC_as_ICP3", "Same_VC_as_TLC","Same_VC_as_K139", "Singleton", "Outlier"))

blast_vCo_comp_only_blasthits <- ggplot(sum_vContact_blast_sum_onlyhits, aes(x=VC_description_blastref, y=Count, fill=VC_description)) + 
  geom_bar(stat='identity', position = position_dodge2(preserve = "single"), width = 1) + 
  scale_fill_manual(values = c(brewer_pal(palette = "Set3")(length(levels(sum_vContact_blast_sum_onlyhits$VC_description))-2), 'lightgrey', 'darkgrey')) + 
  theme_light() + theme(legend.position = 'bottom', axis.text.x=element_text(angle=60,hjust=1)) + ylab('Number of contigs') + xlab('BLASTn hit') 
blast_vCo_comp_only_blasthits

#pdf('./blast_vCo_comp_only_blasthits.pdf', height = 3.8, width = 3.6)
#blast_vCo_comp_only_blasthits + theme(legend.position = 'bottom', axis.text.x=element_text(angle=60,hjust=1))
#dev.off()


# focus on the 'No blast hit'
sum_vContact_blast_no_blasthit$VC_description <- gsub('VC_285', 'VC including Erwinia spp. phages', sum_vContact_blast_no_blasthit$VC_description)

sum_vContact_blast_no_blasthit_sum  <- sum_vContact_blast_no_blasthit %>%  
  group_by(VC_description_blastref, VC_description) %>%  
  summarize(Count = n()) 

sum_vContact_blast_no_blasthit_sum$VC_description_blastref <- ifelse(is.na(sum_vContact_blast_no_blasthit_sum$VC_description_blastref), 'No blast hit', sum_vContact_blast_no_blasthit_sum$VC_description_blastref)
sum_vContact_blast_no_blasthit_sum$VC_description <- factor(sum_vContact_blast_no_blasthit_sum$VC_description, levels = names(sort(table(sum_vContact_blast_no_blasthit$VC_description), decreasing = T)))

blast_vCo_comp_no_blasthit <- ggplot(sum_vContact_blast_no_blasthit_sum, aes(x=VC_description, y=Count, fill=VC_description)) + 
  geom_bar(stat='identity') + 
  scale_fill_manual(values = c('lightgrey', 'darkgrey', brewer_pal(palette = "Set3")(5)[5], rep('black', 14), brewer_pal(palette = "Set3")(1), 'darkmagenta')) + 
  theme_light() + theme(legend.position = 'bottom', axis.text.x=element_text(angle=60,hjust=1)) + ylab('Number of contigs') + xlab('Viral Cluster Assignment') 
blast_vCo_comp_no_blasthit

# export with no legend and finx in illustrator
#pdf('./blast_vCo_comp_no_blasthits_no_legend_1.pdf', height = 3.8, width = 3.6)
#blast_vCo_comp_no_blasthit + theme(legend.position = 'none', axis.text.x=element_text(angle=60,hjust=1))
#dev.off()


# check if the ones with no blast hits have a lower virus score, higher fdr
wilcox.test(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit',]$virus_score, 
            genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description != 'No blast hit',]$virus_score, 
            paired = F)
quantile(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit',]$virus_score)
mean(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit',]$virus_score)
sd(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit',]$virus_score)
quantile(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description != 'No blast hit',]$virus_score)
mean(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description != 'No blast hit',]$virus_score)
sd(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description != 'No blast hit',]$virus_score)

wilcox.test(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit',]$fdr, 
            genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description != 'No blast hit',]$fdr, 
            paired = F)
quantile(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit',]$fdr)
mean(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit',]$fdr)
sd(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit',]$fdr)

quantile(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description != 'No blast hit',]$fdr)
mean(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description != 'No blast hit',]$fdr)
sd(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description != 'No blast hit',]$fdr)

quantile(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit',]$n_hallmarks)
quantile(genomad_virus_vibrio_b[genomad_virus_vibrio_b$blblasthit_Description != 'No blast hit',]$n_hallmarks)

# plot 
genomad_virus_vibrio_b['binary_if_blasthit'] <- ifelse(genomad_virus_vibrio_b$blblasthit_Description == 'No blast hit', 'No blast hit', 'Blast hit')
genomad_virus_vibrio_b_long <- genomad_virus_vibrio_b %>% pivot_longer(cols = c('virus_score', 'fdr'), names_to = 'geNomad_metric', values_to = 'geNomad_value' )

ggplot(genomad_virus_vibrio_b_long) + geom_boxplot(aes(x=binary_if_blasthit, y=geNomad_value)) +facet_grid(geNomad_metric~., scales = 'free_y')



