---
title: "Variant_filtering"
author: "Enrico"
date: "13 April 2020"
output: html_document
---

Once I have completed all of the filtering steps of my VCF, I can proceed to run a Principal Component Analysis, in order to explore the overall population structure and visualize my data. For this I will be using plink 1.9 which is not installed on CESGA's FT2 server. For this I will copy my VCF file to the EBD genomics-a server and run plink there (1.9 is installed):
```
cd ~/LL_selection/LyCaRef_vcfs/
scp csebdjg2@ft2.cesga.es:/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf .
```
To run this I have followed a tutorial found at: https://speciationgenomics.github.io/pca/

I will be working on the VCF from the final filtering step:
```
VCF=~/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf
```
First thing to do is to prune my dataset to eliminate possibly linked SNPs, as PCA assumes the data to be independent. With this command check for correlation between SNPs in 50kbp windows, with a 10bp step size, and the r2 threshold will be 0.1

```
cd ~/LL_selection/pca_plink

# perform linkage pruning - i.e. identify prune sites
plink_1.9 --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out ll_LyCa_ref
```
A threshold of 0.1 leaves us with ~900K SNPs, pruning out ~4M SNPs. I will later value if this is too restrictive or not.

Now I can run PCA on this subset
```
plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ll_LyCa_ref.prune.in \
--make-bed --pca --out ll_LyCa_ref
```
The output includes: eigenvalues (eigenval) and eigenvectors (eigenvec) from PCA, BED/BIM/FAM files (plink format, might be useful for other analyses).

I will copy the eigenval and eigenvec files on my laptop to process them with R
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
In R the first thing to do will be to load the necessary libraries
```{R}
# load tidyverse package
library(tidyverse)
```
Then we can import the data and set it up to be analyzed:
```{R}
# read in data
pca <- read_table2("plink_pca/ll_LyCa_ref.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# add population
loc <- rep(NA, length(pca$ind))
loc[grep("ba", pca$ind)] <- "balkans"
loc[grep("ca", pca$ind)] <- "caucasus"
loc[grep("cr", pca$ind)] <- "carpathians"
loc[grep("ka", pca$ind)] <- "mongolia"
loc[grep("ki", pca$ind)] <- "kirov"
loc[grep("la", pca$ind)] <- "latvia"
loc[grep("no", pca$ind)] <- "norway"
loc[grep("og", pca$ind)] <- "mongolia"
loc[grep("po", pca$ind)] <- "poland"
loc[grep("to", pca$ind)] <- "mongolia"
loc[grep("tu", pca$ind)] <- "tuva"
loc[grep("ur", pca$ind)] <- "urals"
loc[grep("vl", pca$ind)] <- "vladivostok"
loc[grep("ya", pca$ind)] <- "yakutia"

# remake data.frame
pca <- as.tibble(data.frame(pca, loc))
```
To plot the data we first calculate the percentage of variance explained by each PC and plot them
```{R}
# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
# then make a plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
```
Finally, to plot the PCA:
```{R}
b <- ggplot(pca, aes(PC1, PC2, col = loc)) + geom_point(size = 3)
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
```
It appears that the 3 Balkan individuals from LYNX_16 (h_ll_ba_0214,h_ll_ba_0215 and c_ll_ba_0216) are very big outliers in the analysis and are affecting the ability to detect other patterns of variance. I will try eliminating them from the analysis.
```
cd ~/LL_selection/pca_plink

grep -E "h_ll_ba_0214|h_ll_ba_0215|c_ll_ba_0216" ll_LyCa_ref.fam > balkans_samplestoremove.txt

VCF=~/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ll_LyCa_ref.prune.in --remove balkans_samplestoremove.txt \
--make-bed --pca --out ll_LyCa_ref.noba
```
Copy new eigenval and eigenvec files on my laptop
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.noba.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
Analyze with R:
```{R}
# load tidyverse package
library(tidyverse)

# read in data
pca <- read_table2("plink_pca/ll_LyCa_ref.noba.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.noba.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# add population
loc <- rep(NA, length(pca$ind))
loc[grep("ba", pca$ind)] <- "balkans"
loc[grep("ca", pca$ind)] <- "caucasus"
loc[grep("cr", pca$ind)] <- "carpathians"
loc[grep("ka", pca$ind)] <- "mongolia"
loc[grep("ki", pca$ind)] <- "kirov"
loc[grep("la", pca$ind)] <- "latvia"
loc[grep("no", pca$ind)] <- "norway"
loc[grep("og", pca$ind)] <- "mongolia"
loc[grep("po", pca$ind)] <- "poland"
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
ggplotly()

# plot PCA with labels
ggplot(pca, aes(PC1, PC2, col = loc, label=ind)) + geom_point(size = 3) + 
  coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + geom_text(hjust = 1, vjust = 1)
```

