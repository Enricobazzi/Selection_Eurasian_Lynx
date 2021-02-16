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
# copy results:
scp /Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/GenWin_results/*_GenWin_windows_outliers.tsv \
ebazzicalupo@genomics-b.ebd.csic.es:~/GenWin_results/
```
In order to convert the tables to BED format I ran the following:
```
cd ~/GenWin_results/

# Make bed of windows for each variable
for table in $(ls *_outliers.tsv)
 do
  name=($(echo ${table} | sed 's/_GenWin_windows_outliers.tsv//g'))
  echo "${name}"
  tail -n +2 ${table} | awk '{print $6, $1, $2}' | tr ' ' '\t' > ${name}_GenWin_windows_outliers.bed
done
```
I want to try to remove the superoutlier region. I believe it lays between position 17M and 19M of scaffold_17_arrow_ctg1. I will try removing it with bedtools subtract
```
# Create bed of region
paste <(echo "scaffold_17_arrow_ctg1") <(echo "17000000") <(echo "19000000") > superoutlier_region.bed

# Bedtools subtract from all beds
for bed in $(ls ~/GenWin_results/*_GenWin_windows_outliers.bed)
 do
  name=($(echo "${bed}" | rev | cut -d'/' -f1 | rev | sed 's/_GenWin_windows_outliers.bed//g'))
  echo "${name}"
  bedtools subtract -a ${bed} -b superoutlier_region.bed > ${name}_GenWin_windows_outliers.nosuper.bed
done
```
Then I can run PCA with plink on the outlier windows using the BED files
```
VCF=/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.vcf

for bed in $(ls ~/GenWin_results/*.nosuper.bed)
 do
  name=($(echo "${bed}" | rev | cut -d'/' -f1 | rev | sed 's/_GenWin_windows_outliers.nosuper.bed//g'))
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
In R the first thing to do will be to load the necessary libraries and colors
```{R}
# load tidyverse package
library(tidyverse)
library(plotly)
library(RColorBrewer)
library (viridis)

cols <- c("Caucasus"="#B8860b",
          "NE-Poland"=viridis_pal()(5)[3],
          "Urals"="#0F4909",
          "Balkans"="#A035AF",
          "Carpathians"=brewer.pal(12,"Paired")[9],
          "Kirov"=viridis_pal()(5)[1],
          "Latvia"=brewer.pal(12,"Paired")[3],
          "Norway"=viridis_pal()(5)[2],
          "Tuva"=brewer.pal(12,"Paired")[8],
          "Mongolia"=brewer.pal(12,"Paired")[7],
          "Vladivostok"=brewer.pal(12,"Paired")[5],
          "Yakutia"=brewer.pal(12,"Paired")[6],
          "Sierra Morena"=brewer.pal(8, "Greys") [5],
          "Doñana"=brewer.pal(8, "Greys") [8],
          "Caucasus-Pseudodiploid"="#B8860b",
          "NE-Poland-Pseudodiploid"=viridis_pal()(5)[3],
          "Balkans-Pseudodiploid"="#A035AF",
          "Carpathians-Pseudodiploid"=brewer.pal(12,"Paired")[9],
          "Yakutia-Pseudodiploid"=brewer.pal(12,"Paired")[6])

```
Then I can proceed with analyzing the results in a loop
```{R}
variables <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "combined_differentiation", "snow_days", "jan_depth")

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
 loc[grep("ca", pca$ind)] <- "Caucasus"
 loc[grep("ka", pca$ind)] <- "Mongolia"
 loc[grep("ki", pca$ind)] <- "Kirov"
 loc[grep("la", pca$ind)] <- "Latvia"
 loc[grep("to", pca$ind)] <- "Mongolia"
 loc[grep("tu", pca$ind)] <- "Tuva"
 loc[grep("ur", pca$ind)] <- "Urals"
 loc[grep("vl", pca$ind)] <- "Vladivostok"
 loc[grep("ya", pca$ind)] <- "Yakutia"

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
   scale_color_manual(values=cols) +
   ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

 ggsave(paste0(var,"_outliers.PCA.pdf"), path = "PCA_outliers/",
        width=25,height=25,units="cm")
}
```
To do a PCA of the super outlier window
```
# create a range file from BED for PLINK
nlines=($(wc -l < superoutlier_region.bed))
paste superoutlier_region.bed <(yes superoutlier | head -${nlines}) > ~/GenWin_results/superoutlier_range
# run plink
plink_1.9 --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract range ~/GenWin_results/superoutlier_range \
--make-bed --pca --out ~/PCA_outliers/superoutlier
```
Copy the eigen files
```
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/PCA_outliers/superoutlier.eigen* \
Documents/Selection_Eurasian_Lynx/PCA_outliers/
```
Plot PCA in R
```{R}
var <- "superoutlier"
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
loc[grep("ca", pca$ind)] <- "Caucasus"
loc[grep("ka", pca$ind)] <- "Mongolia"
loc[grep("ki", pca$ind)] <- "Kirov"
loc[grep("la", pca$ind)] <- "Latvia"
loc[grep("to", pca$ind)] <- "Mongolia"
loc[grep("tu", pca$ind)] <- "Tuva"
loc[grep("ur", pca$ind)] <- "Urals"
loc[grep("vl", pca$ind)] <- "Vladivostok"
loc[grep("ya", pca$ind)] <- "Yakutia"

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
  scale_color_manual(values=cols) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave(paste0(var,"_outliers.PCA.pdf"), path = "PCA_outliers/",
       width=25,height=25,units="cm")
```
To run PCA on groups of outlier windows (see Window_analysis.md) on genomics-b
I copied the file balkans_samplestoremove.txt from my PCA analysis which was on genomics-a server
```
cd ~/GenWin_results/Window_analysis

VCF=/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.vcf
VCF=/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

for bed in $(ls ~/GenWin_results/Window_analysis/*_windows.bed)
 do
  name=($(echo "${bed}" | rev | cut -d'/' -f1 | rev | sed 's/_windows.bed//g'))
  echo "${name}"

  # create a range file from BED for PLINK
  nlines=($(wc -l < ${bed}))
  paste ${bed} <(yes ${name} | head -${nlines}) > ${name}_range

  # run plink
  plink_1.9 --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --extract range ${name}_range \
  --remove balkans_samplestoremove.txt \
  --make-bed --pca --out ${name}
done
```
And on the small windowset
```
VCF=/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.vcf

for bed in $(ls ~/GenWin_results/small_windowset_nocomb.bed)
 do
  name=($(echo "${bed}" | rev | cut -d'/' -f1 | rev | sed 's/_nocomb.bed//g'))
  echo "${name}"

  # create a range file from BED for PLINK
  nlines=($(wc -l < ${bed}))
  paste ${bed} <(yes ${name} | head -${nlines}) > ${name}_range

  # run plink
  plink_1.9 --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --extract range ${name}_range \
  --remove balkans_samplestoremove.txt \
  --make-bed --pca --out ${name}
done
```
I also run PCA on neutral set of windows (see RDA for range and samplestoremove files), the combined differentiation and the haploblock (superoutlier) using ALL samples:
```
# NEUTRAL
cd ~/GenWin_results/Window_analysis
name=12k10k_neutral
VCF=/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract range /home/ebazzicalupo/RDA/${name}.range \
--remove /home/ebazzicalupo/RDA/samplestoremove.txt --geno 0.001 \
--make-bed --pca --out ${name}

# COMBINED DIFFERENTIATION
cd ~/GenWin_results/Window_analysis
name=combined_differentiation
VCF=/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf
bed=/home/ebazzicalupo/GenWin_results/combined_differentiation_GenWin_windows_outliers.nosuper.bed
name=combined_differentiation
nlines=($(wc -l < ${bed}))
paste ${bed} <(yes ${name} | head -${nlines}) > ${name}_range

plink_1.9 --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract range ${name}_range \
--remove /home/ebazzicalupo/RDA/samplestoremove.txt --geno 0.001 \
--make-bed --pca --out ${name}_persample

# HAPLOBLOCK
cd ~/GenWin_results/Window_analysis
name=superoutlier
VCF=/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

plink_1.9 --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract range /home/ebazzicalupo/GenWin_results/${name}_range \
--remove /home/ebazzicalupo/RDA/samplestoremove.txt --geno 0.001 \
--make-bed --pca --out ${name}_persample

```
I will then copy the eigenval and eigenvec files on my laptop to process them with R:
```
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/GenWin_results/Window_analysis/*eigen* \
Documents/Selection_Eurasian_Lynx/Window_analysis/

# for the neutral PCA
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/GenWin_results/Window_analysis/12k10k_neutral.eigen* \
Documents/Selection_Eurasian_Lynx/Window_analysis/

# for combined diff
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/GenWin_results/Window_analysis/combined_differentiation_persample.eigen* \
Documents/Selection_Eurasian_Lynx/Window_analysis/

# for superoutlier
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/GenWin_results/Window_analysis/superoutlier_persample.eigen* \
Documents/Selection_Eurasian_Lynx/Window_analysis/
```
In R the first thing to do will be to load the necessary libraries and colors
```{R}
# load tidyverse package
library(tidyverse)
library(plotly)
library(RColorBrewer)
library (viridis)

cols <- c("Caucasus"="#B8860b",
          "NE-Poland"=viridis_pal()(5)[3],
          "Urals"="#0F4909",
          "Balkans"="#A035AF",
          "Carpathians"=brewer.pal(12,"Paired")[9],
          "Kirov"=viridis_pal()(5)[1],
          "Latvia"=brewer.pal(12,"Paired")[3],
          "Norway"=viridis_pal()(5)[2],
          "Tuva"=brewer.pal(12,"Paired")[8],
          "Mongolia"=brewer.pal(12,"Paired")[7],
          "Vladivostok"=brewer.pal(12,"Paired")[5],
          "Yakutia"=brewer.pal(12,"Paired")[6])

```
Then I can proceed with analyzing the results in a loop
```{R}
variables <- c("g1", "g2", "g3", "g4", "g5", "g6", "g7","small_windowset")
# for the other give var value and run rest of loop without loop
var <- "12k10k_neutral"
var <- "combined_differentiation_persample"
var <- "superoutlier_persample"
for (i in 1:length(variables)){

 var <- variables[i]

 # import eigen vec and val
 pca <- read_table2(paste0("Window_analysis/",var,".eigenvec"),
                    col_names = FALSE)
 eigenval <- scan(paste0("Window_analysis/",var,".eigenval"))

 # remove nuisance column
 pca <- pca[,-1]
 # set names
 names(pca)[1] <- "ind"
 names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

 # add population
 loc <- rep(NA, length(pca$ind))
 loc[grep("ca", pca$ind)] <- "Caucasus"
 loc[grep("po", pca$ind)] <- "NE-Poland"
 loc[grep("ba", pca$ind)] <- "Balkans"
 loc[grep("cr", pca$ind)] <- "Carpathians"
 loc[grep("no", pca$ind)] <- "Norway"
 loc[grep("ka", pca$ind)] <- "Mongolia"
 loc[grep("ki", pca$ind)] <- "Kirov"
 loc[grep("la", pca$ind)] <- "Latvia"
 loc[grep("to", pca$ind)] <- "Mongolia"
 loc[grep("og", pca$ind)] <- "Mongolia"
 loc[grep("tu", pca$ind)] <- "Tuva"
 loc[grep("ur", pca$ind)] <- "Urals"
 loc[grep("vl", pca$ind)] <- "Vladivostok"
 loc[grep("ya", pca$ind)] <- "Yakutia"

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
   scale_color_manual(values=cols) +
   ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

 ggsave(paste0(var,"_outliers.PCA.pdf"), path = "PCA_outliers/",
        width=25,height=25,units="cm")
}
ggplotly()
```
Investigate PC1 of superoutlier to find Haplotypes of samples in haploblock
```{R}
pctable <- data.frame(PC1=pca$PC1, sample=pca$ind, pop=loc)
#pctable %>% filter(pctable$pop == "Yakutia")

# Draw a pie chart for AF of each population
populations <- as.vector(unique(pctable$pop))
aftable <- data.frame()
for (i in 1:length(populations)){
  popu <- populations[i]
  poppctable <- pctable %>% filter(pctable$pop == popu)
  AF0=((2*length(which(poppctable$PC1<(-0.05))))+length(which(poppctable$PC1>-0.05&poppctable$PC1<0.05)))/(2*length(poppctable$PC1))
  AF1=((2*length(which(poppctable$PC1>0.05)))+length(which(poppctable$PC1>-0.05&poppctable$PC1<0.05)))/(2*length(poppctable$PC1))
  popafs <- data.frame(population=popu, af0=AF0, af1=AF1)
  aftable <- rbind(aftable, popafs)

  pietable <- data.frame(af=c(AF0, AF1), allele=c("0","1"))
  
  AFpie <- ggplot(pietable, aes(x="", y=af, fill=allele)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("lightgrey", "black")) +
  theme_void()+ theme(legend.position = "none")
  
  ggsave(paste0("~/Documents/Selection_Eurasian_Lynx/PCA_outliers/",
                       popu, "_superoutlier_AFpie.png"), AFpie)
  
}
write.table(aftable, "PCA_outliers/superoutlier_popafs.tsv", sep = "\t", row.names = F)
```
