#!/bin/bash
#SBATCH -t 4-00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

###########################################
## GATK 4.1.4.1 Haplotypecaller launcher ##
###########################################

# This program runs GATK 4.1.4.1 Haplotypecaller with standard settings to generate
# a per-chromosome VCF on the finis terrae II server of CESGA.

# The chromosome must be defined while launching the script as such:

# ./cesga_CombineGVCFs.sh <chromosome>

module load gatk/4.1.4.1

###################################
## VARIABLE and PATHS definition ##
###################################

# Reference Genome:
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/Ref_Genome_LyCa/lc4.fa

# Input bams - c_ll_ca_0249 and c_ll_ca_0253 didn't pass QC and are excluded:
INbamARRAY=($(ls /mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_bams/*_indelrealigner.bam | grep -vE "c_ll_ca_0249|c_ll_ca_0253"))

# Output Files:
OUTvcf=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_vcfs/${1}_LyCa_ref.vcf

# chromosome BED file:
BED=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/Ref_Genome_LyCa/CHR_BEDs/${1}_CHR_coordinates.bed


###############################
## GATK 4.1.1.0 CombineGVCFs ##
###############################

gatk HaplotypeCaller \
   -R $REF \
   $(for bam in ${INbamARRAY[@]}; do echo "-I ${bam}";done) \
   -L $BED \
   -O $OUTvcf
