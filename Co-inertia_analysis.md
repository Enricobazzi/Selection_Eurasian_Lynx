---
title: "Co-inertia Analysis"
author: "Enrico"
date: "1 December 2020"
output: html_document
editor_options:
  chunk_output_type: console
---
In this markdown I will describe the steps I took in order to run a co-inertia analysis on my window analysis results.

What I will need is a neutral set of windows to use as the baseline, so I can see how my window analysis results differ from this neutral pattern.

To get the set of neutral windows I will subtract the outlier windows from the intergenic bed (ask DANI how he generated it) and generate a random subset of 10 thousand windows of 12kbp length.
```
cd /home/ebazzicalupo/co-inertia

bedtools random -l 12000 -n 10000 -g <(bedtools subtract -a /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.intergenic.buffer5000.bed  -b /home/ebazzicalupo/GenWin_results/small_windowset.bed) | sort -k1,1 -k2,2n > 12k10k_neutral.bed
```
To run co-inertia I need a matrix of allele frequencies of each population for each SNP inside my windows.
```
# Input VCF
INVCF=/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf
# Reference
REF=/GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.fa
# Array of Populations of input VCF 
popARRAY=($(grep -m1 "#CHROM" ${INVCF} | tr '\t' '\n' | grep "c_ll" | cut -d"_" -f3 | sort -u | sed 's/ka/mo/' | sed 's/to/mo/' | sed 's/og/mo/' | sort -u))
# Array of Window Group beds
wgbedARRAY=($(ls /home/ebazzicalupo/GenWin_results/Window_analysis/g*_windows.bed))


# create VCF for each population
for pop in ${popARRAY[@]}
  do
    
    if [[ ${pop} == mo ]]
     then
      samplesARRAY=($(grep -m1 "#CHROM" ${INVCF} | tr '\t' '\n' | grep -E "c_ll_ka|c_ll_to"))

     else
      samplesARRAY=($(grep -m1 "#CHROM" ${INVCF} | tr '\t' '\n' | grep "c_ll_${pop}" | grep -vE "c_ll_vl_0137|c_ll_tu_0154|h_ll_ba_0214|h_ll_ba_0215|c_ll_ba_0216"))

    fi

    echo "extracting ${pop} population..."

    # select samples from pop i
    /opt/gatk-4.1.0.0/gatk SelectVariants \
    -R $REF \
    -V $INVCF \
    $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
    -O ${pop}.vcf
done

# Get VCF of windows only (neutral + 7 wgs)  
for pop in ${popARRAY[@]}
  do
  echo ${pop}
  # population VCF of neutral windows
  bedtools intersect -a ${pop}.vcf -b 12k10k_neutral.bed > ${pop}.neutral.vcf
  
  # population VCF of each windowgroup windows
  for bed in ${wgbedARRAY[@]}
   do
    wg=($(echo ${bed} | rev | cut -d'/' -f1 | rev | cut -d'_' -f1))
    echo ${wg}
    bedtools intersect -a ${pop}.vcf -b ${bed} > ${pop}.${wg}.vcf
  done
done

# Generate a matrix of allele frequency for each window group (neutral + 7 wgs)
for i in $(ls *.*.vcf | tr ' ' '\n' | cut -d'.' -f2 | sort -u)
 do
 echo ${i}
  echo ${popARRAY[@]} | tr ' ' '\t' > ${i}.afs
  paste <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' ba.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' ca.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' cr.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' ki.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' la.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' mo.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' no.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' po.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' tu.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' ur.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' vl.${i}.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' ya.${i}.vcf | sed 's/AF=//') \
  >> ${i}.afs
done
```
Download the matrixes of afs to laptop to analyze with R:
```
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/co-inertia/*.afs ~/Documents/Selection_Eurasian_Lynx/co-inertia/
```
Prepare R:
```{R}
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ade4)
cols <- c("ca"="#B8860b",
          "po"=viridis_pal()(5)[3],
          "ur"="#0F4909",
          "ba"="#A035AF",
          "cr"=brewer.pal(12,"Paired")[9],
          "ki"=viridis_pal()(5)[1],
          "la"=brewer.pal(12,"Paired")[3],
          "no"=viridis_pal()(5)[2],
          "tu"=brewer.pal(12,"Paired")[8],
          "mo"=brewer.pal(12,"Paired")[7],
          "vl"=brewer.pal(12,"Paired")[5],
          "ya"=brewer.pal(12,"Paired")[6])
```
Explore consistency of Neutral matrix results:
```{R}
# Load the neutral AF matrix:
neutral_matrix <- read_tsv("co-inertia/neutral.afs", 
           col_names = T)
neutral_matrix2 <- t(neutral_matrix)

ran_neutral <- t(sample_n(neutral_matrix, 10000))

pca <- prcomp(ran_neutral, scale. = F)
summary(pca)
plot(pca$x, col = cols)
```
Run co-inertia:
```{R}
wgs <- c("g1", "g2", "g3", "g4", "g5", "g6", "g7")
neutral_matrix <- read_tsv("co-inertia/neutral.afs", 
           col_names = T)

for(i in 1:length(wgs)){
  wg <- wgs[i]
  
  wg_matrix <- t(read_tsv(paste0("co-inertia/", wg,".afs"), 
           col_names = T)[1:1000,])

  n_rand <- ncol(wg_matrix)
  ran_neutral <- t(sample_n(neutral_matrix, n_rand))

  wg_pca <- dudi.pca(wg_matrix, scannf = FALSE, nf = 2)
  neutral_pca <- dudi.pca(ran_neutral, scannf = FALSE, nf = 2)
  
  coin <- coinertia(dudiX = neutral_pca, dudiY = wg_pca, scannf = FALSE, nf = 2)
  
  coin$lX
  plot(coin$lX$AxcX1, coin$lX$AxcX2, main="neutral", col = c(1:12))
  coin$lY
  plot(coin$lY$AxcY1, coin$lY$AxcY2, main="selected", col = c(1:12))
  plot (coin, xax = 1, yax = 2)
}
```


## RDA

I'll try a different approach using RDA.

I chose 1 representant from each variable group:
bio5, bio6, bio13, bio14, jan_depth and snow_days

To get a BED of outlier windows from any of these:
```{bash}
cd GenWin_results
# get all windows WITHOUT combined differentiation
for var in bio5 bio6 bio13 bio14 jan_depth snow_days
 do
  echo ${var}
  cat ${var}_GenWin_windows_outliers.bed >> tmp
done

sort -k1,1 -k2,2n tmp | uniq > sixvars.bed && rm tmp

# get windows in superoutlier region
bedtools intersect -a sixvars.bed -b superoutlier_region.bed -c | awk -F "\t" '{if($4 == 1){print}}' \
> superoutlier_sixvars.bed

# remove superoutlier region windows and add a single window for the whole region
bedtools subtract -a sixvars.bed -b superoutlier_sixvars.bed \
> sixvars_onesuper.bed
paste <(head -1 superoutlier_sixvars.bed | cut -f 1,2) \
<(tail -1 superoutlier_sixvars.bed | cut -f 3) >> sixvars_onesuper.bed
sort -k1,1 -k2,2n sixvars_onesuper.bed | uniq > tmp && mv tmp sixvars_onesuper.bed

# Get small windowset by merging
mergeBed -i sixvars_onesuper.bed > small_windowset_sixvars.bed

# go to co-inertia for af files
cd ~/co-inertia

# filter each vcf for the BED with the outlier windows for the six variables
# Input VCF
INVCF=/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf
# Array of Populations of input VCF 
popARRAY=($(grep -m1 "#CHROM" ${INVCF} | tr '\t' '\n' | grep "c_ll" | cut -d"_" -f3 | sort -u | sed 's/ka/mo/' | sed 's/to/mo/' | sed 's/og/mo/' | sort -u))
for pop in ${popARRAY[@]}
  do
  echo ${pop}
  # population VCF of sixvars windows
  bedtools intersect -a ${pop}.vcf -b ~/GenWin_results/small_windowset_sixvars.bed \
  > ${pop}.sixvars.vcf
done
echo ${popARRAY[@]} | tr ' ' '\t' > sixvars.afs
paste <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' ba.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' ca.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' cr.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' ki.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' la.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' mo.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' no.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' po.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' tu.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' ur.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' vl.sixvars.vcf | sed 's/AF=//') \
  <(grep -o -E 'AF=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' ya.sixvars.vcf | sed 's/AF=//') \
  >> sixvars.afs
```
Download AF to laptop for analysis with R
```
scp ebazzicalupo@genomics-b.ebd.csic.es:~/co-inertia/sixvars.afs Documents/Selection_Eurasian_Lynx/co-inertia/
```
Prepare R:
```{R}
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
```
In R
```{R}
# Load matrix
sixvars_matrix <- t(read_tsv("co-inertia/sixvars.afs", 
           col_names = T)[1:10000,])

# Load environmental data bio5, bio6, bio13, bio14
sixvar_env <- read_tsv("WorldClim_table.tsv", col_names = T)[ c(5,6,13,14) ,] %>% 
  column_to_rownames(., var="pop") %>% t(.)

sixvar_env <- cbind(sixvar_env, t(read_tsv("Snow_table.tsv") %>% column_to_rownames(., var="variable"))) %>% data.frame(.)

# Remove the missing populations [[SHOULD CALCULATE ENV DATA FOR THEM AS WELL]]
sixvars_matrix <- sixvars_matrix[ c(2,4,5,6,9,10,11,12) , ]

rda <- rda(sixvars_matrix ~ ., data=sixvar_env, scale=T)
RsquareAdj(rda)
summary(eigenvals(rda, model = "constrained"))
signif.full <- anova.cca(rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

plot(rda, scaling=3)          # default is axes 1 and 2
plot(rda, choices = c(1, 3), scaling=3)  # axes 1 and 3


ola <- data.frame(pop = c("Caucaus","Kirov","Latvia","Mongolia","Tuva","Urals", "Vladivostok", "Yakutia"), n = c(1:8))
levels(ola$n) <- c("Caucaus","Kirov","Latvia","Mongolia","Tuva","Urals", "Vladivostok", "Yakutia")
eco <- ola$n
bg <- c("#B8860b",viridis_pal()(5)[1],brewer.pal(12,"Paired")[3],brewer.pal(12,"Paired")[7],
        brewer.pal(12,"Paired")[8],"#0F4909",
        brewer.pal(12,"Paired")[5],brewer.pal(12,"Paired")[6]) 

plot(rda, type="n", scaling=3)
points(rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the wolves
text(rda, scaling=3, display="bp", col="#0868ac", cex=1)     # the predictors
legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
```
See RDA.md for details on the per-sample version of RDA