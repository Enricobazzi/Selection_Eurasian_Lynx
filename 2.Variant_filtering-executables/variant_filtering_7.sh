#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

# This script will use the following softwares:

# BEDtools 2.28.0
# BCFtools 1.9
# GATK 4.1.4.1

module load bedtools
module load bcftools
module load gatk/4.1.4.1

# # The VCF file name (without the extensions) must be defined while launching
# the script as such:

# ./Depth_Filter_7.sh <VCFfilename>

###################################
## VARIABLE and PATHS definition ##
###################################

# Reference Genome:
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/Ref_Genome_LyCa/lc4.fa

# AllIndividuals input VCF
INVCF=$LUSTRE/LL_selection/LyCaRef_vcfs/${1}.filter6.vcf

# VCF Directory
OUTdir=$LUSTRE/LL_selection/LyCaRef_vcfs

# Depth per sample Table generated in R
DPStable=$LUSTRE/LL_selection/SamTools_Depth/depth_per_sample.csv

# Create a copy of the VCF file with a new name that will be used to filter out
# variants with excessively low/high depth:
cp $INVCF $OUTdir/${1}.filter7.vcf
OUTVCF=$OUTdir/${1}.filter7.vcf

#############################
## Applying filters - LOOP ##
#############################

# List of all datasets in an array:
popARRAY=($(ls $LUSTRE/LL_selection/LyCaRef_bams/*.depth.bamlist | rev | cut -d'/' -f1 | rev | cut -d'.' -f1))

# For each dataset:
for pop in ${popARRAY[@]}
  do

    # (1) Extract VCF of individuals of dataset
    samplesARRAY=($(cat $LUSTRE/LL_selection/LyCaRef_bams/${pop}.depth.bamlist | rev | cut -d'/' -f1 | rev | cut -d'_' -f1-4))

    gatk SelectVariants \
    -R $REF \
    -V $INVCF \
    $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
    -O $OUTdir/${pop}.depth.vcf

    # Minimum and Maximum depth values
    max=$(grep -w ${pop} ${DPStable} | cut -d',' -f5 | cut -d'.' -f1)
    min=$(grep -w ${pop} ${DPStable} | cut -d',' -f6)
    echo "Maximum depth of ${pop} is ${max}, Minimum depth is ${min}"

    # (2) extract the excessive missingness variant with BCFtools filter
    echo "extracting excessively low/high depth variants from $pop VCF"
    bcftools filter -i "INFO/DP < ${min} || INFO/DP > ${max}" -Ov $OUTdir/${pop}.depth.vcf \
    > $OUTdir/${pop}.applydepthfilter.vcf

    # (3) filter the excessively missing variants from the new VCF file of all samples
    echo "subtracting excessively low/high depth variants of $pop from output VCF"
    bedtools subtract -a $OUTVCF \
    -b $OUTdir/${pop}.applydepthfilter.vcf -header \
    > tmp && mv tmp $OUTVCF

done
