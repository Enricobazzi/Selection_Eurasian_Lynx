#!/bin/bash

##############################
## Generating Allele Counts ##
##############################

# This script will be used to generate Allele Count data from the starting "all-individuals" VCF.
# The VCF file name and absolute path, as well as the output directory, must be
# defined while launching the script as such:

# ./Allele_Count_generation.sh <VCFfile> <OutPutDirectory>

# First first have to divide my VCF with all the individuals, into different VCF for each population.
# GATK SelectVariants version 4.1.0.0 will be used for this step

# Then the number of reference and alternative alleles can be computed and joined into
# a final file to be used as input for BayPass. This file is called all.allelecounts
# and will be located in the defined output directory

###################################
## VARIABLE and PATHS definition ##
###################################

# Input VCF
INVCF="$1"

# Array of Populations of input VCF - remove unwanted populations (ba, og, no, po and cr)
# also change "ka" and "to" populations to Mongolia "mo"
popARRAY=($(grep -m1 "#CHROM" ${INVCF} | tr '\t' '\n' | grep "c_ll" | cut -d"_" -f3 | sort -u | grep -vE "ba|og|no|po|cr" | sed 's/ka/mo/' | sed 's/to/mo/' | sort -u))

# Output Directiory
OUTdir="$2"

# Reference Genome:
REF=/GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.fa

###################################################
## Divide all individuals VCF in Population VCFs ##
###################################################

# Allele count file - make sure to start with an empty one
touch $OUTdir/all.allelecounts > $OUTdir/all.allelecounts

echo " - dividing ${INVCF} by populations and writing allele counts -"

for pop in ${popARRAY[@]}
  do

    if [[ ${pop} == mo ]]
    then
    samplesARRAY=($(grep -m1 "#CHROM" ${INVCF} | tr '\t' '\n' | grep -E "c_ll_ka|c_ll_to" | grep -vE "c_ll_vl_0137|c_ll_tu_0154"))

    else
    samplesARRAY=($(grep -m1 "#CHROM" ${INVCF} | tr '\t' '\n' | grep "c_ll_${pop}" | grep -vE "c_ll_vl_0137|c_ll_tu_0154"))

    fi

    echo "extracting ${pop} population..."

    # select samples from pop i
    /opt/gatk-4.1.0.0/gatk SelectVariants \
    -R $REF \
    -V $INVCF \
    $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
    -O $OUTdir/${pop}.vcf

    echo "counting ${pop} alleles..."

  # Extract allele counts and write file (IMPORTANT - check VCF format for column selection)
    grep -v "#" $OUTdir/${pop}.vcf | \
    grep -o -E 'AC=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}.*AN=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' | \
    cut -d';' -f1,3 | sed 's/AC=//g' | sed 's/;AN=/ /g' | \
    awk -F' ' 'BEGIN { OFS = " " } {print $2-$1,$1}' \
    > $OUTdir/${pop}.allelecounts

  # this if step is to avoid that the first column of the pasted file is empty
    ncols=($(wc -l $OUTdir/all.allelecounts))
    if [[ $ncols -eq 0 ]]
      then
        # First iteration goes as it is
        cat $OUTdir/${pop}.allelecounts > $OUTdir/all.allelecounts
      else
        # else Paste to final Allele Count file
        paste -d' ' $OUTdir/all.allelecounts $OUTdir/${pop}.allelecounts > tmp \
        && mv tmp $OUTdir/all.allelecounts
    fi
done
