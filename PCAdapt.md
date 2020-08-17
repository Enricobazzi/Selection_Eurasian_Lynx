---
title: "PCAdapt"
author: "Enrico"
date: "25 October 2019"
output: html_document
editor_options:
  chunk_output_type: console
---

In this markdown I will describe all the steps I took while analyzing my database with the software PCAdapt (v. XXX for R).

This software works like this XXXXXXX

The first thing I will need to do is divide the VCF into 50 subsets (like I did for the BayPass analysis), and convert them to PLINK format, which is the input format for PCAdapt.

```
cd ~/PCAdapt

VCF=ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc

# Divide in 50 datasets (1 SNP every 50)
for n in {0..49}
 do
  echo "extracting dataset ${n} vcf"
  grep "#" VCF/${VCF}.vcf > VCF/${VCF}.${n}.vcf
  grep -v "#" VCF/${VCF}.vcf | awk -v number="$n" 'NR % 50 == 0+number' >> VCF/${VCF}.${n}.vcf
done

# Change 0 to 50
mv VCF/${VCF}.0.vcf VCF/${VCF}.50.vcf
```
Now I can convert each one to PLINK format
```
for n in {1..50}
 do
  echo "converting dataset ${n} to plink"
  plink_1.9 --vcf VCF/${VCF}.${n}.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --make-bed --out PLINK/${VCF}.${n}
done
```
To analyze the data with PCAdapt I have copied the PLINK files on my laptop in the folder:
/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/PCAdapt_plink
And the I could use R to run the analyses.

First I loaded all the necessary packages and defined the path to my input files:
```{R}
library(pcadapt)
library(tidyverse)
path_to_files <- "/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/PCAdapt_plink/"
poplist.names <- c(rep("Caucasus", 12), rep("Mongolia", 4), rep("Kirov", 13),
                   rep("Latvia", 7), rep("Mongolia", 2), rep("Tuva", 6),
                   rep("Urals", 7), rep("Vladivostok", 9), rep("Yakutia", 9))
```
Firs thing I need to do is to find the right number of principal components (K) to use. To do this I will see, comparing scree plots of XXX, with a high K, if there is consistency in the K plateau across the 50 datasets.
```{R}
for (n in 1:50){
  file <- paste0(path_to_files,
               "ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.",n,".bed")

  filename <- read.pcadapt(file, type = "bed")

  x <- pcadapt(input = filename, K = 10)

  plot(x, option = "screeplot")
}
```
I can see that K=3 seems to be always the start of the plateau, hence the optimal value.

Now I can re-analyze every dataset with K=3 and build a final result table with the information of each SNP together with the calculated p-value.
```{R}
PCAdapt_results <- data.frame()

for (n in 1:50){
  file <- paste0(path_to_files,
               "ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.",n,".bed")
  filename <- read.pcadapt(file, type = "bed")

  x <- pcadapt(input = filename, K = 3)
  
  pvalues <- data.frame(pvalue = x$pvalues, SNPnum = seq(n, 2100553, by = 50))
  
  PCAdapt_results <- data.frame(rbind(PCAdapt_results, pvalues))

}

# order the table based on the SNP number
PCAdapt_results <- PCAdapt_results %>% arrange(SNPnum)

# Add SNP ID information
SNPIDs <- read_tsv("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/finalset.maf5pc.SNPIDs", col_names = F)[,2-3] %>%
    rename("scaffold" =  X2, "position" = X3)

# Add SNP IDs to the total snps table
PCAdapt_results <- data.frame(PCAdapt_results,SNPIDs)
 
# Add Odd/Even for colors in Manhplot
CHR <- data.frame()
for (n in 1:length(unique(PCAdapt_results$scaffold))){
   scaffold_lines <- PCAdapt_results %>% filter(scaffold==unique(PCAdapt_results$scaffold)[n])
   if((n %% 2) == 0) {
    num <- "Even"
   } else {
    num <- "Odd"
   } 
   number <- data.frame(colora = rep(num, nrow(scaffold_lines)))
   CHR <- data.frame(rbind(CHR, number))
}
PCAdapt_results <- data.frame(cbind(PCAdapt_results, CHR))

# Save the table
write.table(x = PCAdapt_results,
            file = paste0("PCAdapt_results/PCAdapt_results.tsv"),
            quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

# Outlier threshold (>99.9%)
pval.thresh=as.numeric(quantile(-log10(PCAdapt_results$pvalue),probs=0.999,na.rm=T))

# Manhattan plot of PCAdapt results for whole genome:
manhplot <- ggplot(PCAdapt_results, aes(x = SNPnum, y = -log10(pvalue), color=colora)) +
    geom_point(alpha = 0.75, stat = "identity", size=1) +
    scale_color_manual(values= c("Black","darkgray")) +
    # TITLE
    # LEGEND
    # X AXIS SCAFFOLD NAMES BELOW
    geom_hline(yintercept=pval.thresh,linetype="dashed", size=0.5, color="red")
manhplot
```
Now I can use the results table as an input file for candidate window analysis with GenWin.