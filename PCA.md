---
title: "PCA"
author: "Enrico"
date: "13 April 2020"
output: html_document
editor_options:
  chunk_output_type: console
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

It appears that the 3 Balkan individuals from LYNX_16 (h_ll_ba_0214,h_ll_ba_0215 and c_ll_ba_0216) are very big outliers in the analysis and are affecting the ability to detect other patterns of variance. I will try eliminating them from the analysis.
```
cd ~/LL_selection/pca_plink

grep -E "h_ll_ba_0214|h_ll_ba_0215|c_ll_ba_0216" ll_LyCa_ref.fam > balkans_samplestoremove.txt

VCF=~/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ll_LyCa_ref.prune.in --remove balkans_samplestoremove.txt \
--make-bed --pca --out ll_LyCa_ref.noba
```
I also want to see if patterns of missingness are affecting my results. For this I will use a different approach with plink: I will not eliminate odd samples (e.g. the 3 balkan individuals), but instead filter out all of the SNPs were information is missing for any individual (--geno 0.001):
```
cd ~/LL_selection/pca_plink

VCF=~/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ll_LyCa_ref.prune.in --geno 0.001 \
--make-bed --pca --out ll_LyCa_ref.nomiss
```
As it appears that missingness is not affecting, I will try to also remove the 3 outlier caucasus samples:
```
cd ~/LL_selection/pca_plink

grep -E "h_ll_ba_0214|h_ll_ba_0215|c_ll_ba_0216|c_ll_ca_0245|c_ll_ca_0248|c_ll_ca_0254" ll_LyCa_ref.fam > balkans-caucasus_samplestoremove.txt

VCF=~/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ll_LyCa_ref.prune.in --geno 0.001 --remove balkans-caucasus_samplestoremove.txt \
--threads 15 --make-bed --pca --out ll_LyCa_ref.nomiss.noba.noca
```
The co-authors also want a version of the PCA with only samples from the Caucasus population.
```
cd ~/LL_selection/pca_plink

# Get list of Caucasus samples:
grep "ll_ca" ll_LyCa_ref.fam > caucasus_samplelist.txt

VCF=~/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ll_LyCa_ref.prune.in --keep caucasus_samplelist.txt --geno 0.001 \
--make-bed --pca --out ll_LyCa_ref.caucasus
```
Also only Caucasus population but without "bad" samples
```
cd ~/LL_selection/pca_plink

# Get list of Caucasus samples:
grep "ll_ca" ll_LyCa_ref.fam | grep -vE "0245|0248|0254" > caucasus9_samplelist.txt

VCF=~/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ll_LyCa_ref.prune.in --keep caucasus9_samplelist.txt --geno 0.001 \
--make-bed --pca --out ll_LyCa_ref.caucasus9
```
For a version of the PCA with only the samples from the Western part of the distribution:
```
cd ~/LL_selection/pca_plink

# Get list of Western samples:
grep -E "ll_ba|ll_cr|ll_ki|ll_la|ll_no|ll_po|ll_ur" ll_LyCa_ref.fam | grep -vE "h_ll_ba_0214|h_ll_ba_0215|c_ll_ba_0216" > western_samplelist.txt

VCF=~/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ll_LyCa_ref.prune.in --keep western_samplelist.txt \
--make-bed --pca --out ll_LyCa_ref.western
```
For a version of the PCA with only the samples from the Eastern part of the distribution:
```
cd ~/LL_selection/pca_plink

# Get list of Eastern samples:
grep -vE "ll_ba|ll_cr|ll_ki|ll_la|ll_no|ll_po|ll_ur|ll_ca" ll_LyCa_ref.fam > eastern_samplelist.txt

VCF=~/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ll_LyCa_ref.prune.in --keep eastern_samplelist.txt \
--make-bed --pca --out ll_LyCa_ref.eastern
```
For a version of the PCA with only the samples from the Eastern part of the distribution, but not the 2 og samples (mongolia):
```
cd ~/LL_selection/pca_plink

# Get list of Eastern samples:
grep -vE "ll_ba|ll_cr|ll_ki|ll_la|ll_no|ll_po|ll_ur|ll_ca|ll_og" ll_LyCa_ref.fam > eastern_noog_samplelist.txt

VCF=~/LL_selection/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ll_LyCa_ref.prune.in --keep eastern_noog_samplelist.txt \
--make-bed --pca --out ll_LyCa_ref.eastern.noog
```
I will copy the eigenval and eigenvec files on my laptop to process them with R:
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
Copy eigenval and eigenvec files of the version without 3 Balkan individuals on my laptop:
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.noba.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
Copy eigenval and eigenvec files of the version without missing individuals on my laptop:
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.nomiss.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
Copy eigenval and eigenvec files of the version without missing and without 3 Balkan and 3 Caucasus individuals on my laptop:
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.nomiss.noba.noca.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
Copy eigenval and eigenvec files of the version with only Caucasus samples on my laptop:
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.caucasus.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
Copy eigenval and eigenvec files of the version with only Caucasus 9 samples on my laptop:
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.caucasus9.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
Copy eigenval and eigenvec files of the version with only Western samples on my laptop:
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.western.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
Copy eigenval and eigenvec files of the version with only Eastern samples on my laptop:
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.eastern.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
Copy eigenval and eigenvec files of the version with only Eastern samples without og on my laptop:
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LL_selection/pca_plink/ll_LyCa_ref.eastern.noog.eigen* Documents/Selection_Eurasian_Lynx/plink_pca/
```
In R the first thing to do will be to load the necessary libraries
```{R}
# load tidyverse package
library(tidyverse)
library(plotly)
```
Then we can import the data and set it up to be analyzed:
```{R}
# read in data - choose the dataset -> ONE ONLY:
# all - no filter
pca <- read_table2("plink_pca/ll_LyCa_ref.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.eigenval")
# no balkans
pca <- read_table2("plink_pca/ll_LyCa_ref.noba.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.noba.eigenval")
# no missingness
pca <- read_table2("plink_pca/ll_LyCa_ref.nomiss.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.nomiss.eigenval")
# no missingness - no balkan no caucasus
pca <- read_table2("plink_pca/ll_LyCa_ref.nomiss.noba.noca.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.nomiss.noba.noca.eigenval")
# Caucasus only:
pca <- read_table2("plink_pca/ll_LyCa_ref.caucasus.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.caucasus.eigenval")
# Caucasus 9 only:
pca <- read_table2("plink_pca/ll_LyCa_ref.caucasus9.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.caucasus9.eigenval")
# Intergenic:
pca <- read_table2("plink_pca/ll_intergenic_LyCa_ref.noba.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_intergenic_LyCa_ref.noba.eigenval")
# Western:
pca <- read_table2("plink_pca/ll_LyCa_ref.western.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.western.eigenval")
# Eastern:
pca <- read_table2("plink_pca/ll_LyCa_ref.eastern.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.eastern.eigenval")
# Eastern - no og:
pca <- read_table2("plink_pca/ll_LyCa_ref.eastern.noog.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca/ll_LyCa_ref.eastern.noog.eigenval")

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

# for caucasus only:
loc <- rep(NA, length(pca$ind))
loc[grep("0240", pca$ind)] <- "Lesser Caucasus"
loc[grep("0241", pca$ind)] <- "Dagestan"
loc[grep("0242", pca$ind)] <- "Lesser Caucasus"
loc[grep("0243", pca$ind)] <- "Dagestan"
loc[grep("0244", pca$ind)] <- "Dagestan"
loc[grep("0245", pca$ind)] <- "Dagestan"
loc[grep("0247", pca$ind)] <- "Lesser Caucasus"
loc[grep("0248", pca$ind)] <- "Western Greater Caucasus"
loc[grep("0252", pca$ind)] <- "Western Greater Caucasus"
loc[grep("0254", pca$ind)] <- "Lesser Caucasus"
loc[grep("0259", pca$ind)] <- "Dagestan"
loc[grep("0260", pca$ind)] <- "Dagestan"


# remake data.frame
pca <- as.tibble(data.frame(pca, loc))
```
To plot the data we first calculate the percentage of variance explained by each PC and plot them
```{R}
# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
pve <- data.frame(PC = 1:9, pve = eigenval/sum(eigenval)*100)
# then make a plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
```
Finally, to plot the PCA:
```{R}
# plot PCA without labels
ggplot(pca, aes(PC1, PC2, col = loc, label=ind)) + geom_point(size = 3) +
  coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# Save the plot - Choose appropriate name for version of data loaded
#ggsave("PCA.all.pdf", path = "plink_pca/", width=25,height=25,units="cm")
#ggsave("PCA.noba.pdf", path = "plink_pca/", width=25,height=25,units="cm")
#ggsave("PCA.nomiss.pdf", path = "plink_pca/", width=25,height=25,units="cm")
#ggsave("PCA.nomiss.noba.noca.pdf", path = "plink_pca/", width=25,height=25,units="cm")
#ggsave("PCA.western.pdf", path = "plink_pca/", width=25,height=25,units="cm")
#ggsave("PCA.eastern.pdf", path = "plink_pca/", width=25,height=25,units="cm")
#ggsave("PCA.eastern.noog.pdf", path = "plink_pca/", width=25,height=25,units="cm")
#ggsave("PCA.caucasus9.pdf", path = "plink_pca/", width=25,height=25,units="cm")

ggplotly()

# plot PCA with labels
library(ggrepel)

ggplot(pca, aes(PC1, PC2, col = loc, label=ind)) + geom_point(size = 3) +
  coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + geom_label_repel()
```
3D plot - (colors missing):
```{R}
library(threejs)
scatterplot3js(pca$PC1, pca$PC2, pca$PC3,
              labels=pca$ind,
               size=0.7, grid=FALSE)
```
