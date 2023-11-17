
#####
###### Bash scripts for metagenomics reads analyses

### Reads taxonomic annotation
### kraken on raw reads
kraken2-build --db=~/kraken2_db --threads=10
cd  ~/fastq_reads

for i in *_R1.fastq.gz
do
   prefix=$(basename $i _R1.fastq.gz)    #D05181978
    kraken2   --use-names --threads 48  --db ~/kraken2_db  --report  kraken_out/${prefix}.kraken   --gzip-compressed  --paired  ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz > kraken_out/${prefix}_tax
done

###
######
###########
## bracken on kraken's output
## build kraken database
bracken-build -d  ~/kraken2_db -t 10

cd  kraken_out
for k in *.kraken
do
prefix=$(basename $k .kraken)
bracken -d ~/kraken2_db -i ${prefix}.kraken -o Braken_out/${prefix}.bracken -r 150 -l S -t 32
done

## merge bracken output
kraken-biom  kraken_out/*_bracken_species.kraken  -o  bracken_kraken_out.biom --fmt json

## import bracken merged output in R
biomfilename ='bracken_kraken_out.biom'
data <- import_biom(biomfilename, parseFunction=parse_taxonomy_default)

###
####
#######
###### metagenomics reads cleaning

## remove PhiX genome and the human genome, the decontamined fastq are in 'decontamined_cholera_fastq' folder
bowtie2  -x  human_phix_ref/ref/GRCh38_Phix -1 ${prefix}_1.fastq.gz  -2  ${prefix}_2.fastq.gz  --un-conc  decontamined_cholera_fastq/${prefix}_%.fastq.gz  -S  bowtie_output/${prefix}.sam 


##
###
iu-gen-configs samples_v1.txt  -o config_out_v1 ##avant c t config_out_v1  
for ini in  config_out_v1/*.ini; do iu-filter-quality-minoche --ignore-deflines  $ini; done


####### Read assembling
#### MEGAHIT
megahit -1 ${prefix}-QUALITY_PASSED_R1.fastq  -2 ${prefix}-QUALITY_PASSED_R2.fastq     --min-contig-len  1000  -m 0.9  -o  assemblies/${prefix}/  -t 40

## mapping reads against all the assemblies
1. build index for the assemblies

bowtie2-build $i  mapping/Bowtie_build_output/$prefix/contigs  --threads 32

2. mapping reads from each sample $i against the assemblies

cd /home/naima17/projects/def-shapiro/naima17/cholera_v2/Instrain_pipe/mapping/Bowtie_build_output


for i in *

do
   prefix=$(basename ${i}) ##1274300
   mkdir /home/naima17/scratch/cholera_nov2020/mapping/Bam_files/${prefix}
   shuf -n 20 /home/naima17/projects/def-shapiro/naima17/cholera_v2/Instrain_pipe/mapping/sample_list_filt.csv > /home/naima17/projects/def-shapiro/naima17/cholera_v2/Instrain_pipe/mapping/liste.csv
   ##add the name of the sample being mapped to the list
   ### si le sample lui meme n’est pas ds la liste, ajouter le sample a la liste
   if ! (grep -Fxq $prefix mapping/liste.csv) ; then echo $prefix >> mapping/liste.csv; fi
   for j in `cat mapping/liste.csv`
   do

        prefix2=$j
        bowtie2 --threads 32 -x mapping/Bowtie_build_output/$i/contigs -1 filtred_to_analyse/${prefix2}-QUALITY_PASSED_R1.fastq -2 filtred_to_analyse/${prefix2}-QUALITY_PASSED_R2.fastq -S mapping/Bam_files/${prefix}/${prefix2}.sam
        samtools view -S -b mapping/Bam_files/${prefix}/${prefix2}.sam > mapping/Bam_files/${prefix}/${prefix2}-RAW.bam
        samtools sort -@ 32 mapping/Bam_files/${prefix}/${prefix2}-RAW.bam -o mapping/Bam_files/${prefix}/${prefix2}.bam
        rm mapping/Bam_files/${prefix}/${prefix2}.sam mapping/Bam_files/${prefix}/${prefix2}-RAW.bam
        samtools index mapping/Bam_files/${prefix}/${prefix2}.bam

   done
done


### metabat binning

## create depth file (metabat and Concoct pipeline)
for i in mapping/Bam_files/*
do
prefix=$(basename ${i})   ##D17181835
jgi_summarize_bam_contig_depths --outputDepth    mapping/depth_files/${prefix}_metabat_depth.txt  $i/*.bam
cut_up_fasta.py renamed_assemblies/${prefix}_final_contigs.fa -c 10000 -o 0 --merge_last -b  ${prefix}_contigs_10K.bed  > ${prefix}_contigs_10K.fa
concoct_coverage_table.py  ${prefix}_contigs_10K.bed  $i/*.bam >  mapping/depth_files/${prefix}_concoct_coverage_table.tsv
done

###
#######
### run metabat
cd  mapping/depth_files_metabat
for i in   *.txt
do
prefix=$(basename $i _metabat_depth.txt)
mkdir  Metabat/${prefix}
metabat -i  renamed_assemblies/${prefix}_final_contigs.fa   -a ${prefix}_metabat_depth.txt  -o   Metabat/${prefix}/bin  -t 32
done

######
#########
#### run concoct
cd   mapping/depth_files_concoct
for i in   *.tsv
do
prefix=$(basename $i _concoct_coverage_table.tsv)
mkdirConcoct/${prefix}
concoct --composition_file mapping/bed_files/${prefix}_contigs_10K.fa  --coverage_file mapping/depth_files_concoct/${prefix}_concoct_coverage_table.tsv -b   Concoct/${prefix}  --threads 32
merge_cutup_clustering.py Concoct/${prefix}/clustering_gt1000.csv > Concoct/${prefix}/clustering_merged.csv
mkdir Concoct/${prefix}/fasta_bins
extract_fasta_bins.py renamed_assemblies/${prefix}_final_contigs.fa   Concoct/${prefix}/clustering_merged.csv --output_path  Concoct/${prefix}/fasta_bins
done

## Dastool
###### DAS Tool is an automated method that integrates the results of a flexible number of binning algorithms 
###### to calculate an optimized, non-redundant set of bins from a single assembly

for  d  in /Metabat/*
do
prefix=$(basename $d )
Fasta_to_Scaffolds2Bin.sh  -e fa -i  Metabat/${prefix}  > dastool/sample_data/${prefix}_scaffold2bin_metabat.tsv
perl -pe "s/,/\t/g;" Concoct/${prefix}/clustering_merged.csv > dastool/sample_data/${prefix}_scaffold2bin_concoct.tsv
mkdir dastool_output/${prefix}
DAS_Tool  -i sample_data/${prefix}_scaffold2bin_metabat.tsv,sample_data/${prefix}_scaffold2bin_concoct.tsv -l metabat,concoct -c renamed_assemblies/${prefix}_final_contigs.fa  -o  dastool_output/${prefix}/DASToolRun1  --write_bins 1 -t 32 
done

#### checkM 
for  d  in dastool_output/*
do

prefix=$(basename $d )
cd $d/DASToolRun1_DASTool_bins
mkdir checkm/${prefix}
checkm lineage_wf -x fa .  checkm/${prefix}
done

## rename all bins from dastools
for  d  in dastool_output/*
do
prefix=$(basename $d )
 for FILENAME in $d/DASToolRun1_DASTool_bins/* ; do
            prefix2=$(basename $FILENAME )
            NEWNAME=$prefix'_'$prefix2
            cp $FILENAME  dastool_allBins/$NEWNAME
 done
done

## drep
mkdir complete_only_MAGs
dRep dereplicate complete_only_MAGs  -g  dastool_allBins/*.fa    -p 32

## ##GTDBK, taxonomic annotation of the MAGs from drep
gtdbtk identify --genome_dir  complete_only_MAGs/dereplicated_genomes  --out_dir GTDBK/gtdbtk_id  -x fa --cpus 32
gtdbtk align --identify_dir GTDBK/gtdbtk_id --out_dir GTDBK/gtdbtk_align --cpus 32
mkdir gtbk_mmap
gtdbtk classify --genome_dir  complete_only_MAGs/dereplicated_genomes  --align_dir GTDBK/gtdbtk_align --out_dir GTDBK/gtdbtk_classify -x fa --cpus 32 --scratch_dir  gtbk_mmap/mmap-file

## run inStrain 
## concatenate all MAGs and ICP1
## all MAGs + NCBI ICP1(NC_015157.1 Vibrio phage ICP1, complete genome)
## all_MAGs_phage.fa

cd dereplicated_genomes

cat *.fa > allGenomes_v1/all_MAGs.fa
cp allGenomes_v1/all_MAGs.fa   ~/allGenomes_v1/all_MAGs_phage.fa
cat ICP1_ncbi.fasta >>  ~/allGenomes_phages_v1/all_MAGs_phage.fa

## run inStrain, mapp reads against the concatenated MAGs+ICP1 file
mkdir  ~/allGenomes_v1/all_genomes_index
bowtie2-build    ~/allGenomes_v1/all_MAGs_phage.fa   ~/allGenomes_v1/all_genomes_index/contigs  --threads 32

cd  filtred_to_analyse

for i in *_R1.fastq
do
   prefix=$(basename $i _R1.fastq)
   bowtie2 --threads 32 -x ~/allGenomes_v1/all_genomes_index/contigs  -1  ${prefix}_R1.fastq  -2  ${prefix}_R2.fastq -S /sam_Instrain_v1/${prefix}.sam
done


## create stb file
parse_stb.py --reverse -f complete_only_MAGs/dereplicated_genomes/* -o representitve_genomes.stb

## Prodigal
prodigal -i allGenomes_v1/all_MAGs.fa   -d all_MAGs_genes_v1.fna -o prodigal_bact.out
prodigal -i allGenomes_v1/phages.fa -d  phages_v1.fna -o prodigal_phages.out -p meta
prodigal -i allGenomes_v1/all_MAGs_phage.fa   -d all_MAGs_phage_genes_v1.fna -o prodigal_bact_phage.out

## run inStrain
cd sam_Instrain_v1
for i in *sam
prefix=$(basename $i sam)
do
inStrain profile  $i  allGenomes_v1/all_MAGs_phage.fa  -o  Instrain_profiles_v1/IS_${prefix} -p 6  -g prodigal/all_MAGs_phage_genes_v1.fna  -s  representitve_genomes_phage.stb  --database_mode
done