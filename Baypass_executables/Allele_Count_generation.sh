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

# Array of Populations of input VCF
popARRAY=($(grep "#CHROM" "$INVCF" | tr '\t' '\n' | grep "c_ll" | cut -d"_" -f3 | sort -u))

# Output Directiory
OUTdir="$2"

# Reference Genome:
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/test/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa

###################################################
## Divide all individuals VCF in Population VCFs ##
###################################################

# Allele count file
touch $OUTdir/all.allelecounts

echo " - dividing $INVCF by populations and writing allele counts -"

for i in ${popARRAY[@]}
  do
    samplesARRAY=($(grep "#CHROM" "$INVCF" | tr '\t' '\n' | grep "c_ll_$i"))

    echo "extracting $i population..."

    # select samples from pop i
    /opt/gatk-4.1.0.0/gatk SelectVariants \
    -R $REF \
    -V $INVCF \
    $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
    -O $OUTdir/${i}.vcf

    echo "counting $i alleles..."

  # Extract allele counts and write file
    grep -v "#" $OUTdir/${i}.vcf | cut -d';' -f1,3 | \
    grep -o -E 'AC=[[:digit:]]{1,3};AN=[[:digit:]]{1,3}' | \
    sed 's/AC=//g' | sed 's/;AN=/ /g' | awk -F' ' 'BEGIN { OFS = " " } {print $2-$1,$1}' \
    > $OUTdir/${i}.allelecounts

  # this if step is to avoid that the first column of the pasted file is empty
    ncols=($(wc -l $OUTdir/all.allelecounts))
    if [[ $ncols -eq 0 ]]
      then
        # First iteration goes as it is
        cat $OUTdir/${i}.allelecounts > $OUTdir/all.allelecounts
      else
        # else Paste to final Allele Count file
        paste -d' ' $OUTdir/all.allelecounts $OUTdir/${i}.allelecounts > tmp \
        && mv tmp $OUTdir/all.allelecounts
    fi
done
