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
./Allele_Count_generation.sh \
/home/ebazzicalupo/BayPass/VCF/ll_perspecies.trimmed_filtered1.ann_wout_no_po_ba_no_fixed_max_2alleles_min_0.02.vcf \
/home/ebazzicalupo/BayPass/AlleleCounts

```
