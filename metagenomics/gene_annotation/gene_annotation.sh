## Gene prediction with prodigal and gene annotation with eggNOG
## NM nov 2023

## Prodigal
prodigal -i ~/allGenomes_v1/all_MAGs_phage.fa  -d  all_MAGs_phage_genes.fna -o prodigal_bact_phage.out

## Eggnog on the result of prodigal avec cette cde : 
emapper.py -i all_MAGs_phage_genes.fna -o ~/output_folder --itype CDS —translate

