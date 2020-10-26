---
title: "PCA_outliers"
author: "Enrico"
date: "13 April 2020"
output: html_document
editor_options:
  chunk_output_type: console
---

In this R markdown I will run PCA focusing only on candidate regions from the selection scan analyses.

First I will upload my tables of outlier SNPs to the EBD genomics-b server to convert them into BED files for filtering out non-candidate regions
```
# copy CORE+PCAdapt combined results:
scp /Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/GenWin_results/combined_differentiation_GenWin_windows_outliers.tsv \
ebazzicalupo@genomics-b.ebd.csic.es:~/GenWin_results/

# copy AUX results:
scp /Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/GenWin_results/*_WholeGenomeOutliers.tsv \
ebazzicalupo@genomics-b.ebd.csic.es:~/GenWin_results/
```
In order to convert the tables to BED format I ran the following:
```
cd ~/GenWin_results/

for table in $(ls *.tsv)
 do
  name=($(echo ${table} | cut -d'_' -f1))
  if [[ ${name} == combined ]]
   then
   name="differentiation"
  fi
  echo "${name}"
  tail -n +2 ${table} | awk '{print $6, $1, $2}' > ${name}_W_outliers.bed
done
```
Then I can run PCA with plink on the outlier windows using the BED files
```
VCF=/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.vcf

for bed in $(ls ~/GenWin_results/*_W_outliers.bed)
 do
  name=($(echo "${bed}" | rev | cut -d'/' -f1 | rev | cut -d'_' -f1))
  echo "${name}"

  # create a range file from BED for PLINK
  nlines=($(wc -l < ${bed}))
  paste ${bed} <(yes ${name} | head -${nlines}) > ~/GenWin_results/${name}_range

  # run plink
  plink_1.9 --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --extract range ~/GenWin_results/${name}_range \
  --make-bed --pca --out ~/PCA_outliers/${name}
done
```
I will then copy the eigenval and eigenvec files on my laptop to process them with R:
```
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/PCA_outliers/*eigen* \
Documents/Selection_Eurasian_Lynx/PCA_outliers/
```
In R the first thing to do will be to load the necessary libraries
```{R}
# load tidyverse package
library(tidyverse)
library(plotly)
```
Then I can proceed with analyzing the results in a loop
```{R}
variables <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "differentiation")

for (i in 1:length(variables)){

 var <- variables[i]
 
 # import eigen vec and val
 pca <- read_table2(paste0("PCA_outliers/",var,".eigenvec"),
                    col_names = FALSE)
 eigenval <- scan(paste0("PCA_outliers/",var,".eigenval"))
 
 # remove nuisance column
 pca <- pca[,-1]
 # set names
 names(pca)[1] <- "ind"
 names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
 
 # add population
loc <- rep(NA, length(pca$ind))
loc[grep("ca", pca$ind)] <- "caucasus"
loc[grep("ka", pca$ind)] <- "mongolia"
loc[grep("ki", pca$ind)] <- "kirov"
loc[grep("la", pca$ind)] <- "latvia"
loc[grep("to", pca$ind)] <- "mongolia"
loc[grep("tu", pca$ind)] <- "tuva"
loc[grep("ur", pca$ind)] <- "urals"
loc[grep("vl", pca$ind)] <- "vladivostok"
loc[grep("ya", pca$ind)] <- "yakutia"

# remake data.frame
pca <- as.tibble(data.frame(pca, loc))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
# then make a plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# plot PCA without labels
ggplot(pca, aes(PC1, PC2, col = loc, label=ind)) + geom_point(size = 3) +
  coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave(paste0(var,"_outliers.PCA.pdf"), path = "PCA_outliers/",
       width=25,height=25,units="cm")
}
```

```{R}
```