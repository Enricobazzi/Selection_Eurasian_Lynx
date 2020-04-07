#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

#############################################################
START=$(date)
echo "Missingness_Filter_6 SCRIPT for $1 starting : $START"
#############################################################

####################################
## Missingness Filtering launcher ##
####################################

# With this script I want to apply a filter based on data missingness to my
# dataset of Lynx lynx individuals, composed of 12 populations.
# I will keep only variants which are present in all popolations and in each population
# at least in four individuals (N_SAMPLES-N_MISSING is the number of non-missing samples).

# To do so I will use the following softwares:

# BEDtools 2.28.0
# BCFtools 1.9

module load bedtools
module load bcftools

# # The VCF file name (without the .filter5.vcf extension) must be defined while launching
# the script as such:

# ./variant_filtering_6.sh <VCFfilename>

#####################################
## Applying filters - Preparations ##
#####################################

# List populations from namelist files in an array (for loop):
popARRAY=($(ls $LUSTRE/LL_selection/LyCaRef_bams/*.namelist | rev | cut -d'/' -f1 | rev | grep -v "all-samples" | cut -d '.' -f1 | sort -u))

# Create a copy of the VCF file with a new name that will be used to filter out
# excessively missing variants:
cp $LUSTRE/LL_selection/LyCaRef_vcfs/${1}.filter5.vcf $LUSTRE/LL_selection/LyCaRef_vcfs/${1}.filter6.vcf

#######################################
## Applying filters - The great LOOP ##
#######################################

for pop in ${popARRAY[@]}
  do

# (2) use the namelist file to divide the VCF by pop
  echo "filtering ${pop} individuals from original VCF"
  bcftools view -S $LUSTRE/LL_selection/LyCaRef_bams/${pop}.namelist -Ov \
  $LUSTRE/LL_selection/LyCaRef_vcfs/${1}.filter5.vcf \
  > $LUSTRE/LL_selection/LyCaRef_vcfs/${pop}_LyCa_ref.filter5.subset.vcf

# (3) extract the excessive missingness variant with BCFtools filter and
# (4) filter the excessively missing variants from the new VCF file of all samples
  echo "extracting missing variants from ${pop} VCF and filtering them out"

  bcftools filter -i "N_SAMPLES-N_MISSING < 4" -Ov $LUSTRE/LL_selection/LyCaRef_vcfs/${pop}_LyCa_ref.filter5.subset.vcf \
  > $LUSTRE/LL_selection/LyCaRef_vcfs/${pop}_LyCa_ref.filter5.subset.missing.vcf

  bedtools subtract -a $LUSTRE/LL_selection/LyCaRef_vcfs/${1}.filter6.vcf \
  -b $LUSTRE/LL_selection/LyCaRef_vcfs/${pop}_LyCa_ref.filter5.subset.missing.vcf -header \
  > tmp && mv tmp $LUSTRE/LL_selection/LyCaRef_vcfs/${1}.filter6.vcf

done

###########################################################
END=$(date)
echo "Missingness_Filter_6 SCRIPT for $1 ended : $END"
###########################################################
