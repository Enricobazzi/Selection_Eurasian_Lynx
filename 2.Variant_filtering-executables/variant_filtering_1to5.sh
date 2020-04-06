#!/bin/bash
#SBATCH -t 1-00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

#######################
## Variant Filtering ##
#######################

# With this script I want to apply all of the General Hard filters to my VCF dataset.
# These include:

# (1) Repetitive/Low mappability regions
# (2) Indels + Non-biallelic sites
# (3) Lynx lynx exclusive substitutions (frequency = 1 in all populations)
# (4,5) Hard quality filters, as of GATK standard practices

# It will generate an output hard-filtered.vcf file and a log file reporting the number
# of variants removed and the number of variants left at each step.

# As these filters are independent of species and sequencing technology they can be
# applied to the complete VCF file directly. A detailed explenation for each
# filtering step is available in the 2.Variant_Filtering.md MarkDown file.

# This script will run on the finis terrae II server of CESGA using the following softwares:

# GATK 4.1.1.4
# BEDtools 2.28.0
# BCFtools 1.9

# The VCF file name (without the .vcf extension) must be defined while launching
# the script as such:

# ./variant_filtering_1to5.sh <VCFfilename>

module load gatk/4.1.4.1
module load bedtools
module load bcftools

###################################
## VARIABLE and PATHS definition ##
###################################

# VCFfilename
filename=($(echo ${1}))

# Reference Genome:
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/Ref_Genome_LyCa/lc4.fa

# Input VCF File:
INvcf=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_vcfs/${filename}.vcf

# BED File of Masked regions:
MASKbed=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/Ref_Genome_LyCa/repetitive_regions/lc_rep_ALL_scaffold_coord.bed

# Step 1 output VCF
ST1out=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_vcfs/${filename}.filter1.vcf

# Step 2 output VCF
ST2out=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_vcfs/${filename}.filter2.vcf

# Step 3 output VCF
ST3out=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_vcfs/${filename}.filter3.vcf

# Step 4 output VCF
ST4out=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_vcfs/${filename}.filter4.vcf

# Step 5 output VCF
ST5out=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_vcfs/${filename}.filter5.vcf


############################################
## (1) Repetitive/Low mappability regions ##
############################################
echo "step 1"
# Apply the filter with BedTools subtract
bedtools subtract -a ${INvcf} -b ${MASKbed} -header | uniq > ${ST1out}

######################################
## (2) Indels + Non-biallelic sites ##
######################################
echo "step 2"
# Apply the filter with GATK SelectVariants
gatk SelectVariants \
  -select-type SNP \
  --restrict-alleles-to BIALLELIC \
  -R ${REF} \
  -V ${ST1out} \
  -O ${ST2out}

##################################################################
## (3) Lynx genus wide exclusive substitutions from Felis catus ##
##################################################################
echo "step 3"
# Apply the filter with BCFtools view
bcftools view -e 'INFO/AF=1.00' ${ST2out} > ${ST3out}

######################################
## (4) and (5) Hard quality filters ##
######################################
echo "step 4"
# Filter all except for the RanksSums:
gatk SelectVariants \
  --selectExpressions "QUAL >= 30 && QD >= 2.0 && FS <= 60.0 && MQ >= 40.0" \
  -R ${REF} \
  -V ${ST3out} \
  -O ${ST4out}

echo "step 5"
# Filter RankSums with bcftools view:
bcftools view -e 'INFO/MQRankSum<-12.5 | INFO/ReadPosRankSum<-8.0' ${ST4out} > ${ST5out}
