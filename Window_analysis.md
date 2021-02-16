---
title: "Window_analysis"
author: "Enrico"
date: "10 November 2020"
output: html_document
editor_options:
  chunk_output_type: console
---
In this MD I will run exploratory analysis on my GenWin results. The objective is to find groups of variables (worldclim or snow) that can be grouped together for the fact that they show similar genomic association results.

The first idea is to run a Multiple Correspondence Analysis (MCA). Using a list of all possible outlier windows in the different variables, the fact of being an outlier (1) or not (0) will be used as a (binary) categorical variable for the MCA.

For this I will generate two lists of all possible outlier windows. One where I merge overlapping windows coming from different variables (small windowset), and one where I remove overlaps by dividing windows into non-overlapping fragments (big windowset).

Working on genomics-b server in the GenWin_results folder.

Get dataset of all windows and unite haploblock into one window.
NOTE: I also have added a version of ALL windows WITHOUT the combined_differentiation results
```
cd ~/GenWin_results

# Make BED files from outlier TSV files
for table in $(ls *_outliers.tsv)
 do
  name=($(echo ${table} | sed 's/_GenWin_windows_outliers.tsv//g'))
  echo "${name}"
  tail -n +2 ${table} | awk '{print $6, $1, $2}' | tr ' ' '\t' > ${name}_GenWin_windows_outliers.bed
done

# get all windows
cat *_GenWin_windows_outliers.bed | sort -k1,1 -k2,2n | uniq > all_windows.bed

# get windows in superoutlier region
bedtools intersect -a all_windows.bed -b superoutlier_region.bed -c | awk -F "\t" '{if($4 == 1){print}}' \
> superoutlier_region_windows.bed

# remove superoutlier region windows and add a single window for the whole region
bedtools subtract -a all_windows.bed -b superoutlier_region_windows.bed > all_windows_onesuper.bed
paste <(head -1 superoutlier_region_windows.bed | cut -f 1,2) \
<(tail -1 superoutlier_region_windows.bed | cut -f 3) >> all_windows_onesuper.bed
sort -k1,1 -k2,2n all_windows_onesuper.bed | uniq > tmp && mv tmp all_windows_onesuper.bed

## WITHOUT the combined_differentiation results

# get all windows WITHOUT combined differentiation
cat [^c]*_GenWin_windows_outliers.bed | sort -k1,1 -k2,2n | uniq > all_windows_nocomb.bed

# get windows in superoutlier region
bedtools intersect -a all_windows_nocomb.bed -b superoutlier_region.bed -c | awk -F "\t" '{if($4 == 1){print}}' \
> superoutlier_region_windows_nocomb.bed

# remove superoutlier region windows and add a single window for the whole region
bedtools subtract -a all_windows_nocomb.bed -b superoutlier_region_windows_nocomb.bed \
> all_windows_nocomb_onesuper.bed
paste <(head -1 superoutlier_region_windows_nocomb.bed | cut -f 1,2) \
<(tail -1 superoutlier_region_windows_nocomb.bed | cut -f 3) >> all_windows_nocomb_onesuper.bed
sort -k1,1 -k2,2n all_windows_nocomb_onesuper.bed | uniq > tmp && mv tmp all_windows_nocomb_onesuper.bed

```
Get small windowset by merging overlapping windows:
```
# Use bedtools merge (mergeBed) to merge overlapping features into single feature
mergeBed -i all_windows_onesuper.bed > small_windowset.bed

## WITHOUT the combined_differentiation results
mergeBed -i all_windows_nocomb_onesuper.bed > small_windowset_nocomb.bed

```
This has generated a set of 888 windows.

To get the big windowset using bedops --partition (I need to install it first):
```
conda create -n bedops
conda activate bedops
conda install -c bioconda bedops
bedops --partition all_windows_onesuper.bed > big_windowset.bed

## WITHOUT the combined_differentiation results
bedops --partition all_windows_nocomb_onesuper.bed > big_windowset_nocomb.bed
```
This has generated a set of 1458 windows.

Check which windows from the two windowsets are outliers for each variable using bedtools intersect
```
varlist=($(ls *_GenWin_windows_outliers.bed | sed 's/_GenWin_windows_outliers.bed//'))
for var in ${varlist[@]}
 do

  echo ${var}
  bedtools intersect -a big_windowset.bed -b ${var}_GenWin_windows_outliers.bed -c |
  awk -F "\t" -OFS="\t" '{if($4 >= 1){print $1, $2, $3, 1} else {print $0}}' > ${var}_big_intersect.bed

  bedtools intersect -a small_windowset.bed -b ${var}_GenWin_windows_outliers.bed -c |
  awk -F "\t" -OFS="\t" '{if($4 >= 1){print $1, $2, $3, 1} else {print $0}}' > ${var}_small_intersect.bed

done
```
Download these intersects on laptop to work with R
```
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/GenWin_results/*_intersect.bed ~/Documents/Selection_Eurasian_Lynx/Window_analysis
```
Prepare R:
```{R}
library(tidyverse)
variables <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "combined_differentiation", "snow_days", "jan_depth")
```
Upload the results into R creating a matrix
```{R}
# version with 0 and 1
for (i in 1:length(variables)){

  var=variables[i]

  intersect <- read_tsv(paste0("Window_analysis/", var, "_big_intersect.bed"),
           col_names = F) %>%
    rename("scaffold" =  X1, "start" = X2, "end" = X3)

  if (i==1){
    names(intersect)[names(intersect) == "X4"] <- var
    matrix <- unite(intersect, "window", scaffold, start, end, sep="-")
  } else {
    colu <- intersect$X4
    matrix <- cbind(matrix, colu)
    names(matrix)[names(matrix) == "colu"] <- var
  }
}

# version with var_yes/var_no
for (i in 1:length(variables)){

  var=variables[i]

  intersect <- read_tsv(paste0("Window_analysis/", var, "_big_intersect.bed"),
           col_names = F) %>%
    rename("scaffold" =  X1, "start" = X2, "end" = X3)

  if (i==1){
    intersect$X5 <- str_replace(intersect$X4, "1", paste0(var, "_yes")) %>%
      str_replace(., "0", paste0(var, "_no"))
    intersect <- intersect[,-4]
    names(intersect)[names(intersect) == "X5"] <- var
    matrix <- unite(intersect, "window", scaffold, start, end, sep="-")
  } else {
    colu <- intersect$X4 %>% str_replace(., "1", paste0(var, "_yes")) %>%
      str_replace(., "0", paste0(var, "_no")) %>% str_replace(., "bio1bio10_no_yes", "bio10_no")
    matrix <- cbind(matrix, colu)
    names(matrix)[names(matrix) == "colu"] <- var
  }
}

# create final matrix (ROWS <- if I want row names)
ROWS <- matrix$window
matrix1 <- as.matrix(data.frame(matrix[-1], row.names = ROWS))
matrix2 <- t(matrix1)
```
Correspondece analysis using the package CA
```{R}
library(ca)
library(ggrepel)

ca1 <- ca(matrix1)
ca2 <- ca(matrix2)

plot(ca1)
plot(ca2)

ca2_data_frame <- data.frame(var=row.names(ca2$rowcoord[,1:2]), ca2$rowcoord[,1:2])
ca1_data_frame <- data.frame(var=1:1530, ca1$rowcoord[,1:2])

# plot PCA with labels
winplot <- ggplot(data = ca1_data_frame, aes(x = Dim1, y = Dim2)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_vline(xintercept = 0, colour = "black") +
  geom_point(alpha = 0.8, size = 5) +
  labs(x = "Dim 1", y = "Dim 2") +
  ggtitle("CA plot of variables") +
  theme_bw()

varplot <- ggplot(data = ca2_data_frame, aes(x = Dim1, y = Dim2, label = var)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_vline(xintercept = 0, colour = "black") +
  geom_point(alpha = 0.8, size = 5) +
  labs(x = "Dim 1", y = "Dim 2") +
  ggtitle("CA plot of variables") +
  geom_label_repel() +
  theme_bw()

```
Modularity and Nestedness
```{R}
library(bipartite)
wnodf=networklevel(matrix1, index="weighted NODF")
visweb(matrix1, type="nested")

modreal <- computeModules(web = matrix2, method = "Beckett")
#modreal_max <- metaComputeModules(moduleObject = matrix1, method = "Beckett", N=10)

modreal@likelihood
plotModuleWeb(modreal)
printoutModuleInformation(modreal)
modinfo <- listModuleInformation(modreal)
```
Plot Unipartite matrix (number windows shared between variables)
```{R}
library(igraph)
unimat <- as.one.mode(matrix2, project = "higher")

igraph.mat <- graph_from_adjacency_matrix(unimat, diag = F, mode = "undirected", weighted = T)
E(igraph.mat)
plot(igraph.mat, edge.width = (E(igraph.mat)$weight)/30, layout = layout_nicely, vertex.label.color="black", vertex.color="lightblue", vertex.label.cex=0.5)

visweb(unimat)
```
MCA trial with FactoMineR
```{R}
#install.packages("FactoMineR")
library(FactoMineR)
library(plotly)

mca <- MCA(matrix1, graph=FALSE)
eig.val <- mca$eig
sum(eig.val[, 2])
barplot(eig.val[, 2],
        names.arg = 1:nrow(eig.val),
        main = "Variances Explained by Dimensions (%)",
        xlab = "Principal Dimensions",
        ylab = "Percentage of variances",
        col ="steelblue")
varplot <- plot(mca, invisible = c("ind", "quali.sup", "quanti.sup"),
     cex = 0.8,
     autoLab = "no")

winplot <- plot(mca,
     invisible = c("var", "quali.sup", "quanti.sup"),
     cex = 0.8,                                    
     autoLab = "yes")
ggplotly(winplot, tooltip = "all")
```

Having run different analyses with binary results for each variable, I want to try to get continuous values instead. I will do this by calculating the weighted mean of the Wstat for each window in the small windowset, which contains larger windows created by collapsing all windows which have been calculated as outlier for any variable (see above).

To calculate the weighted mean for each window, I also need to have a BED file of all windows for each variable, with its Wstat value as 4th column (after normal chr, start, end of a BED file). This is generated from the TSV of GenWin results for all windows:
```
cd ~/GenWin_results

for table in $(ls *_windows.tsv)
 do
  name=($(echo ${table} | sed 's/_GenWin_windows.tsv//g'))
  echo "${name}"
  tail -n +2 ${table} | awk '{print $6, $1, $2, $5}' | tr ' ' '\t' |
  sort -k1,1 -k2,2n | grep -v "NA" > ${name}_GenWin_windows.bed
done
```
Intersect the windowset with bed file calculating the wx and w values for the weighted mean formula for which the weighted mean is equal to:
∑ weight (w) * value (x) ÷ ∑ weight (w)
or the sum of values multiplied by their weight divided by the sum of weights:
```
varlist=($(ls *_GenWin_windows_outliers.bed | sed 's/_GenWin_windows_outliers.bed//'))

for var in ${varlist[@]}
 do
  echo ${var}

  bedtools intersect -a small_windowset.bed -b ${var}_GenWin_windows.bed -wo |
  awk '{print $1, $2, $3, $7*($8/10000), $8/10000}' | sort -k1,1 -k2,2n | tr ' ' '\t' \
  > ${var}_small_windowset_weights.bed
done

## WITHOUT the combined_differentiation results
varlist=($(ls *_GenWin_windows_outliers.bed | sed 's/_GenWin_windows_outliers.bed//' | grep -v "combined_differentiation"))

for var in ${varlist[@]}
 do
  echo ${var}

  bedtools intersect -a small_windowset_nocomb.bed -b ${var}_GenWin_windows.bed -wo |
  awk '{print $1, $2, $3, $7*($8/10000), $8/10000}' | sort -k1,1 -k2,2n | tr ' ' '\t' \
  > ${var}_small_windowset_weights_nocomb.bed
done
```
Now to calculate the weighted mean of Wstat values for each window in the windowset:
```
for var in ${varlist[@]}
 do
  echo ${var}
  windows=($(cut -f1-3 ${var}_small_windowset_weights.bed | tr '\t' '-' | uniq))
  touch ${var}_small_windowset_wmeans.bed
  for window in ${windows[@]}
    do
     win=($(echo ${window} | tr '-' '\t'))
     grep -f <(echo ${win[@]} | tr ' ' '\t') ${var}_small_windowset_weights.bed |
     awk '{ wx += $4; w += $5 } END { print $1, $2, $3, wx/w; }' | uniq | tr ' ' '\t' \
     >> ${var}_small_windowset_wmeans.bed
  done
done

## WITHOUT the combined_differentiation results
for var in ${varlist[@]}
 do
  echo ${var}
  windows=($(cut -f1-3 ${var}_small_windowset_weights_nocomb.bed | tr '\t' '-' | uniq))
  touch ${var}_small_windowset_wmeans_nocomb.bed
  for window in ${windows[@]}
    do
     win=($(echo ${window} | tr '-' '\t'))
     grep -f <(echo ${win[@]} | tr ' ' '\t') ${var}_small_windowset_weights_nocomb.bed |
     awk '{ wx += $4; w += $5 } END { print $1, $2, $3, wx/w; }' | uniq | tr ' ' '\t' \
     >> ${var}_small_windowset_wmeans_nocomb.bed
  done
done
```
Download data to laptop for analysis with R:
```
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/GenWin_results/*_small_windowset_wmeans.bed ~/Documents/Selection_Eurasian_Lynx/Window_analysis

## WITHOUT the combined_differentiation results
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/GenWin_results/*_small_windowset_wmeans_nocomb.bed ~/Documents/Selection_Eurasian_Lynx/Window_analysis

```
Prepare R:
```{R}
library(tidyverse)
library(stats)
library(pheatmap)
library(viridis)

variables <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "combined_differentiation", "snow_days", "jan_depth")

## WITHOUT the combined_differentiation results
variables <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "snow_days", "jan_depth")

```
Upload the results into R creating a matrix
```{R}
for (i in 1:length(variables)){

  var=variables[i]
  intersect <- read_tsv(paste0("Window_analysis/", var, "_small_windowset_wmeans.bed"),
           col_names = F) %>%
    rename("scaffold" =  X1, "start" = X2, "end" = X3)

  if (i==1){
    names(intersect)[names(intersect) == "X4"] <- var
    matrix <- unite(intersect, "window", scaffold, start, end, sep="-")
  } else {
    colu <- intersect$X4
    matrix <- cbind(matrix, colu)
    names(matrix)[names(matrix) == "colu"] <- var
  }
}

## WITHOUT the combined_differentiation results
for (i in 1:length(variables)){

  var=variables[i]
  intersect <- read_tsv(paste0("Window_analysis/", var, "_small_windowset_wmeans_nocomb.bed"),
           col_names = F) %>%
    rename("scaffold" =  X1, "start" = X2, "end" = X3)

  if (i==1){
    names(intersect)[names(intersect) == "X4"] <- var
    matrix <- unite(intersect, "window", scaffold, start, end, sep="-")
  } else {
    colu <- intersect$X4
    matrix <- cbind(matrix, colu)
    names(matrix)[names(matrix) == "colu"] <- var
  }
}

# create final matrix (ROWS <- if I want row names)
ROWS <- matrix$window
matrix1 <- as.matrix(data.frame(matrix[-1])) # columns = variables
matrix2 <- t(matrix1) # columns = windows
```
Calculate euclidian distance between windows from matrix1 which has them as rows and apply Ward's Hierarchical Clustering algorithm to find groups. Create a table for each group with its windows
```{R}
# Groups of windows based on association with variables
eu_dist1 <- dist(matrix1, method = "euclidean")
clust1 <- hclust(eu_dist1, method = "ward.D")
membs1 <- cutree(clust1, k=7)
plot(clust1, hang = -1, cex = 0.1, labels = FALSE, ylab = NULL)
# order of groups in dendrogram is:
# 1, 7, 2, 5, 3, 4, 6 with comb
# 4, 6, 3, 5, 2, 1, 7 without comb

# Also find groups of variables based on the GEA results
eu_dist2 <- dist(matrix2, method = "euclidean")
clust2 <- hclust(eu_dist2, method = "ward.D")
membs2 <- cutree(clust2, k=6)
plot(clust2, hang = -1, cex = 0.6, ylab = NULL)

# write a table for each window group:
w_members_all <- data.frame(cbind(matrix, membs1))
for (i in 1:7){
 w_members <- w_members_all %>% subset(., membs1 == i)

 write.table(x = w_members,
            file = paste0("Window_analysis/g", i ,"_windows.tsv"),
            quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")
}
```
Use pheatmap to draw a heatmap of the clusters of variables and windows
```{R}
pheatmap(matrix2, scale = "column", cluster_rows = clust2, cluster_cols = clust1,
         color = viridis(17, option =  "B"), cutree_cols = 8)
```
Copy the window groups tables to genomics-b
```
scp Documents/Selection_Eurasian_Lynx/Window_analysis/g*_windows.tsv ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/GenWin_results/Window_analysis
```
Transform the window group tables into BED files:
```
cd ~/GenWin_results/Window_analysis

for table in $(ls g*_windows.tsv)
 do
  name=($(echo ${table} | sed 's/_windows.tsv//g'))
  echo "${name}"
  tail -n +2 ${table} | cut -f1 | tr '-' '\t' |
  sort -k1,1 -k2,2n > ${name}_windows.bed
done
```
On genomics-b run a PCA on each group of windows (see PCA_outliers.md)

Intersect with the reference genome annotation to see list of genes in each group of windows:
```
cd ~/GenWin_results/Window_analysis

# Intersect beds with gff3 file:
for bed in $(ls *_windows.bed)
 do
  name=($(echo ${bed} | sed 's/_windows.bed//g'))
  echo ${name}
  bedtools intersect -wa -wb -a ${bed} \
  -b /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.gff3 \
  > ${name}.intersect
done

# get gene list for each variable:
for intersect in $(ls *.intersect)
 do
  name=($(echo ${intersect} | sed 's/.intersect//g'))
  echo ${name}
  grep -Eo "gene=.*;" ${intersect} | cut -d';' -f1 | sort -u > ${name}.genes
done
```
And for the whole small windowset without the combined differentiation
```
bedtools intersect -wa -wb -a ../small_windowset_nocomb.bed \
-b /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.gff3 \
> small_windowset_nocomb.intersect

grep -Eo "gene=.*;" small_windowset_nocomb.intersect | cut -d';' -f1 | sort -u \
> small_windowset_nocomb.genes
```
