---
title: "BayPass"
author: "Enrico"
date: "15 October 2019"
output: html_document
---

In this Markdown I will describe and perform all the steps to perform a genome-wide scan for signatures of selection in Eurasian Lynx populations, using the software BayPass 2.1

The analyses will be performed on the genomics-b server of EBD, unless specified otherwise.

This software uses Allele Frequency data and Bayesian Hierarchical Models to generate a distribution of differentiation coefficients of SNPs across the genome. First a scaled population covariance matrix based on population allele frequencies is generated, to evaluate the relationships between populations. Based on this covariance matrix, and the supposed ancestral allele frequency (inferred from weighted mean of the reference allele frequency), a CORE model is generated.

SNPs that present values of differentiation (XtX, a SNP-specific Fst explicitly corrected for the scaled covariance of population allele frequencies) that exceed the amount expected under the core model can be identified as candidate selection loci.

Furthermore, this software allows the evaluation of the association of particular SNPs to some environmental Covariate. With one measurement for population for each covariate, the software will evaluate the data under to additional models:

(1) The Standard Covariate model (STD): which adds an association covariable (given by the correlation coefficient between the covariate measurement and XtX) to the CORE model.
(2) The Auxiliary Variable Covariate model (AUX): which further builds on the STD model by attaching a binary variable (0 or 1 -> association or no association) to each locus' regression coefficient. The posterior mean of this variable will indicate the posterior probability of the association of that variable with a particular SNP.

In order to run the software we will need:

– Allele Count data for all the considered populations: in the form of a space delimited file, with one row for each SNP and two columns for each population, with the allele counts for the reference and alternative alleles.

– Environmental variable files: in the form of different files with one column per population and their respective measurement for that particular covariate.

The Allele Count file will have to be generated from the VCF file containing information on all the populations of interest.

The Environmental variable files have already been generated and can be found at /home/ebazzicalupo/BayPass/Covariate_Data in the genomics-b server of EBD.

## Generating Allele Count data file

In order to generate the Allele Count file I will first have to divide my VCF with all the individuals, into different VCF for each population. Then the number of reference and alternative alleles can be computed and joined into a final file to be used as input for BayPass.

The script to perform this step can be found at Baypass_executables/Allele_Count_generation.sh and was launched as such

```
screen -S allelecounts
script allelecounts_baypass.log

./Allele_Count_generation.sh \
/home/ebazzicalupo/BayPass/VCF/ll_perspecies.trimmed_filtered1.ann_wout_no_po_ba_no_fixed_max_2alleles_min_0.02.vcf \
/home/ebazzicalupo/BayPass/AlleleCounts
```

## CORE Model

Now I will run BayPass with the generated allele count data using the CORE model with no Co-Variate data.

```
screen -S baypass_CORE
script baypass_CORE.log

# Define working directory:
WD=/home/ebazzicalupo/BayPass

# Run BayPass
g_baypass -npop 10 -gfile $WD/AlleleCounts/all.allelecounts -outprefix $WD/OutPut/CORE

```
The results under the CORE model can be analyzed in R on my laptop. The results directory (OutPut) was copied in the project directory

```{R}

# Load Libraries and Functions #
require(corrplot) ; require(ape)
library(geigen)
source("/Users/enricobazzicalupo/Desktop/Pruebas_Baypass/utils/baypass_utils.R")

# Upload estimate of omega #
omega=as.matrix(read.table("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/OutPut/CORE_mat_omega.out"))
pop.names=c("CR","KA","KI","LA","OG","TO","TU","UR","VL","YA")
dimnames(omega)=list(pop.names,pop.names)

# Compute and visualize the correlation matrix #
cor.mat=cov2cor(omega)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))
# Visualize the correlation matrix as hierarchical clustering tree #
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
plot(bta14.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

# Estimates of the XtX differentiation measures #
anacore.snp.res=read.table("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/OutPut/CORE_summary_pi_xtx.out",h=T)
plot(anacore.snp.res$M_XtX)


```
