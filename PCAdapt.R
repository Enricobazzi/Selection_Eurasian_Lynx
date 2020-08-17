#############
## PCAdapt ##
#############

# In this script I will run the R package PCAdapt on my local computer, using the same starting VCF data as BayPass,
# which includes 10 populations of eurasian lynx.

# First I need to convert the VCF file to PLINK format with vcftools (on genomics EBD server):

# vcftools --vcf /home/ebazzicalupo/PCAdapt/VCF/ll_perspecies.trimmed_filtered1.ann_wout_no_po_ba_no_fixed_max_2alleles_min_0.02.vcf \
# --plink --out /home/ebazzicalupo/PCAdapt/PLINK/ll_perspecies.PLINK

# Then I convert it to plinkBED with PLINK

# plink_1.9 --file /home/ebazzicalupo/PCAdapt/PLINK/ll_perspecies.PLINK --make-bed --out /home/ebazzicalupo/PCAdapt/PLINK/ll_perspecies.PLINK

# I will follow the guide found at :

# https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html

# Load Package:
# install.packages("pcadapt")
library(pcadapt)

# Upload VCF data file
path_to_file <- "/Users/enricobazzicalupo/Documents/PCAdapt/PLINK/ll_perspecies.PLINK.bed"
filename <- read.pcadapt(path_to_file, type = "bed")

##############################################
## Principal Component Analysis of the Data ##
##############################################

# Run the pcadapt command
# computing test statistics and p-values based on the correlations between SNPs and the first K principal components (PCs)
x <- pcadapt(input = filename, K = 10)

# Scree plot -> Visualize in descending order the amount of variation explained by each PC
plot(x, option = "screeplot")

# Define populations
poplist.int <- c(rep(1, 6), rep(2, 4), rep(3, 13), rep(4, 6), rep(5, 2), rep(6, 2), rep(7, 6), rep(8, 6), rep(9, 8), rep(10, 8))
poplist.names <- c(rep("cr", 6), rep("ka", 4), rep("ki", 13), rep("la", 6), rep("og", 2), rep("to", 2), rep("tu", 6), rep("ur", 6), rep("vl", 8), rep("ya", 8))

# Choosing the number K of Principal Components - See when it looses structure?

# Plot PC1 vs PC2
plot(x, option = "scores", pop = poplist.names)
# Plot PC3 vs PC4
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
# Plot PC4 vs PC5
plot(x, option = "scores", i = 4, j = 5, pop = poplist.names)

###############################################
## Computing the test statistic based on PCA ##
###############################################

# The test statistic for detecting outlier SNPs is the Mahalanobis distance:
# a multi-dimensional approach that measures how distant is a point from the mean
# Optimized by selecting optimum number of PCs
lala <- data.frame(pvalue = x$pvalues)
summary(x)
# Manhattan plot of -log10 of p-values
plot(x , option = "manhattan")
# check uniform distribution of p-values with q-q plot
plot(x, option = "qqplot")
# Histogram of p-values
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x, option = "stat.distribution")
