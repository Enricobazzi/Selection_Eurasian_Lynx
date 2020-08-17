---
title: "GenWin"
author: "Enrico"
date: "6 August 2020"
output: html_document
editor_options:
  chunk_output_type: console
---

In this markdown I will use the R package GenWin to analyze the results from PCAdapt and BayPass. 

GenWin is a program that trys to detect inflection points in summary statistics calculated in a locus per locus manner (e.g. Fst). Using these inflection points it will divide the genome in windows which are not of arbitrary size, but based on properties of the summary statistic values. The program will output the number of SNPs inside each window and give a summary value based on the spline (Wstat). Wstat is a value that depends on the mean of the summary statistic, weighted against the mean and standard deviation of the dataset and the number of SNPs inside each window.

Wstat can be then used to calculate outlier windows with a quantile criteria (e.g. above 99.9% as done in the GenWin paper).

I will use it to divide the Lynx lynx genome in windows based on different summary statistics.

## Outliers based on levels of Differentiation

The combined p-value obtained from the combination of BayPass' CORE model and PCAdapt results, will be used to find outlier windows based on levels of differentiation between populations.
The p-value for the XtX values of CORE model will be calculated based on the distribution of the PODs simulated XtX values, which represents the null hypothesis of no selection.

Load the necessary packages into R
```{R}
library(tidyverse)
#library(metaRNASeq)
library(GenWin)
```
Upload the results from the CORE model, PODs and PCAdapt into R
```{R}
PCAdapt_results=read.table(paste0("PCAdapt_results/PCAdapt_results.tsv"),h=T)
XtX_results=read.table(paste0("BayPass_results/CORE_XtX_results.tsv"),h=T)
pod.xtx=read.table("BayPass_OutPut/PODS_summary_pi_xtx.out",h=T)$M_XtX
```
Create a combined dataframe of the 2 tables, with a column for scaffold, position, XtX value, p-value calculated from XtX distribution, p-value calculated by PCAdapt.
```{R}
# Function for z.test p-value
z.test = function(x,mu,popvar){
  one.tail.p <- NULL
  z.score <- (mean(x)-mu)/(popvar/sqrt(length(x)))
  one.tail.p <- pnorm(abs(z.score),lower.tail = FALSE)
}
# Mean and SD of POD XtX values for the z-test
mu <- mean(pod.xtx)
popvar <- sd(pod.xtx)

# Calculate p-values for all XtX values of CORE model
p.xtx <- c()
for (n in 1:nrow(XtX_results)){
  p=z.test(XtX_results$M_XtX[n],mu,popvar)
  p.xtx[n] <- p
}

# Combine into a single table all of the results
combined.table <- data.frame(scaffold = XtX_results$scaffold, position = XtX_results$position,
                             SNPnum = XtX_results$SNPnum, colora = XtX_results$colora, 
                             M_XtX = XtX_results$M_XtX,
                             p_XtX = p.xtx,
                             p_PCA = PCAdapt_results$pvalue)

# Remove NA values (PCAdapt has some)
combined.table <- na.omit(combined.table)
```
Then I can calculate the combined p-value from the 2 methods using Fisher's method
```{R}
# get a list of p-values
pvals <- list(p_XtX = combined.table$p_XtX, p_PCA = combined.table$p_PCA)

# create an empty vector to add combined p-values
p_combined <- c()

# for all p-values pairs calculate combined p-value with fisher's method
# NOTE: combined p-value of 2 points is too low and is approximated by R to 0. This causes
# problems with the log and is manually changed to be highest possible p-value (among the data)
for (n in 1:nrow(combined.table)){
  p2 <- c(pvals[[1]][n], pvals[[2]][n])
  p_combined2=pchisq((sum(log(p2))*-2), df=length(p2)*2, lower.tail=F)
  if (p_combined2 == 0.000000e+00) {
    p_combined2=0.1e-300
    p_combined[n] <- p_combined2
  } else {
    p_combined[n] <- p_combined2
  }
}

# add the raw p-value calculated with metaRNASeq package to the table
combined.table <- data.frame(combined.table, p_combined = p_combined,
                             log_p_comb = -log10(p_combined))


# Write the table with the results
write.table(x = combined.table,
            file = paste0("GenWin_results/combined_differentiation_results_nona.tsv"),
            quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")
```
Draw the manhattan plot of the combined p-values
```{R}
# Outlier threshold (>99.9%)
pval.thresh=as.numeric(quantile(combined.table$log_p_comb,probs=0.999,na.rm=T))

# Manhattan plot of combined p-values for whole genome:
manhplot <- ggplot(combined.table, aes(x = SNPnum, y = log_p_comb, color=colora)) +
    geom_point(alpha = 0.75, stat = "identity", size=1) +
    scale_color_manual(values= c("Black","darkgray")) +
    # TITLE
    # LEGEND
    # X AXIS SCAFFOLD NAMES BELOW
    geom_hline(yintercept=pval.thresh,linetype="dashed", size=0.5, color="red")
manhplot
```
Analyze the combined p-value dataset with GenWin
```{R}
# Load table if not already loaded
combined.table=data.frame(read.table("GenWin_results/combined_differentiation_results_nona.tsv",
                          h=T))

# Empty Spline output table to fill with GenWin results
all_spline <- data.frame()

# Loop through all scaffolds as GenWin works on one scaffold at the time
for (n in 1:length(unique(combined.table$scaffold))) {

 # Scaffold:
 chr=as.character(unique(combined.table$scaffold)[n])
 # Get SNPs table for the scaffold only
 data_chr <- combined.table %>% filter(scaffold==chr)

 # Analyze data with GenWin - Produces a table and a plot
 spline <- splineAnalyze(data_chr$log_p_comb, data_chr$position,
                         smoothness=100, plotRaw=F, plotWindows=F,  method=3)
 # Get the table with the results
 spline.data <- as.data.frame(spline$windowData)
 # Add the scaffold name to the results table
 spline.data$scaffold <- rep(chr,nrow(spline.data))

 # Alternate Odd and Even - for later plot of all scaffolds
 if((n %% 2) == 0) {
  spline.data$color <- rep("Even",nrow(spline.data))
 } else {
  spline.data$color <- rep("Odd",nrow(spline.data))
 } 
 # Fill the GenWin results table
 all_spline <- data.frame(rbind(all_spline, data.frame(spline.data)))
}

# Add the Window length
all_spline$WindowLength <- all_spline$WindowStop-all_spline$WindowStart
# Add the Window Number
all_spline$WindowNumber <- 1:nrow(all_spline)

# Write final table
write.table(x = all_spline, 
            file = paste0("GenWin_results/combined_differentiation_GenWin_windows.tsv"),
            quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

# Calculate Outlier value
all_spline.thresh=as.numeric(quantile(all_spline$Wstat,probs=0.999,na.rm=T))

# Write table of outlier windows
all_spline_total_outliers <- subset(all_spline, Wstat > all_spline.thresh)

write.table(x = all_spline_total_outliers, 
          file = paste0("GenWin_results/combined_differentiation_GenWin_windows_outliers.tsv"),
          quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

# Manhattan plot of window results for whole genome:
manhplot <- ggplot(all_spline, aes(x = WindowNumber, y = Wstat, color=color)) +
    geom_point(alpha = 0.75, stat = "identity", size=1) +
    scale_color_manual(values= c("Black","darkgray")) +
    # TITLE
    # LEGEND
    # X AXIS SCAFFOLD NAMES BELOW
    geom_hline(yintercept=all_spline.thresh,linetype="dashed", size=0.5, color="red")
manhplot
```
Create a manhattan plot for each scaffold of the combined p-values and of the resulting GenWin calculated windows
```{R}
# Loop through all scaffolds
for (n in 1:length(unique(combined.table$scaffold))) {
 # Scaffold:
 chr=unique(combined.table$scaffold)[n]
 # Get SNPs table for the scaffold only
 data_chr <- combined.table %>% filter(scaffold==chr)
 # Get Window table for the scaffold only 
 chr_spline <- all_spline %>% filter(scaffold==chr)
 
 # Manhattan plot of combined p-values for scaffold:
 manhplot <- ggplot(data_chr, aes(x = SNPnum, y = log_p_comb, color=colora)) +
    geom_point(alpha = 0.75, stat = "identity", size=1) +
    scale_color_manual(values= c("Black","darkgray")) +
    # TITLE
    # LEGEND
    geom_hline(yintercept=pval.thresh,linetype="dashed", size=0.5, color="red")
 manhplot
 # Save the scaffold plot
  filename <- paste("GenWin_results/", chr, "_combined_manhplot_pvalues.pdf", sep="")
  ggsave(filename)  


 # Manhattan plot of window results for scaffold bubbly:
 manhplot <- ggplot(chr_spline, aes(x = WindowNumber, y = Wstat)) +
     geom_point(alpha = 0.75, stat = "identity", size=chr_spline$SNPcount*0.1) +
     # TITLE
     # LEGEND
     geom_hline(yintercept=all_spline.thresh,linetype="dashed", size=0.5, color="red")
  manhplot
  # Save the scaffold plot
  filename <- paste("GenWin_results/", chr,"_bubble_WStat_manhplot.pdf", sep="")
  ggsave(filename)  
}
```

## Outliers based on levels of association to external variables

To detect windows of high association with the different bioclimatic variables, I will use the Bayes Factor, calculated by BayPass, as a summary statistic.

```{R}
# Define the variables to work in a loop
variables <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

for (k in 1:length(variables)){

 # Variable:
 var=variables[k]
 # Empty SNPs dataframe to fill
 snps.table <- data.frame()
 
 # Fill SNPs dataframe with all BayPass results for the variable
 for (n in 1:50){
  # upload the dataset table:
  var.snp=read.table(paste0("BayPass_OutPut/AUX_",  var,"_",n,"_summary_betai.out"),h=T)
  
  # calculate the SNP database number sequence (1 every 50) 
  # NOTE: 2100553 is the total number of SNPs in all datasets
  var.snp <- data.frame(var.snp, SNPnum = seq(n, 2100553, by = 50))
  
  # add the dataset rows to the total snps table
  snps.table <- rbind(snps.table, var.snp)
 }
 # order the table based on the SNP number
 snps.table <- snps.table %>% arrange(SNPnum)
 
 # Add SNP ID information from SNPIDs table
 SNPIDs <- read_tsv("BayPass_OutPut/finalset.maf5pc.SNPIDs", col_names = F)[,2-3] %>%
    rename("scaffold" =  X2, "position" = X3)
  
 # Add SNP IDs to the total snps table
 snps.table <- data.frame(snps.table,SNPIDs)
 
 # Empty Spline output table to fill with GenWin results
 all_spline <- data.frame()
 # Empty Spline output table to fill with the OUTLIERS of GenWin results
 all_spline_outliers <- data.frame()
 
 # Loop through all scaffolds as GenWin works on one scaffold at the time
 for (n in 1:length(unique(snps.table$scaffold))) {
 
  # Scaffold:
  chr=unique(snps.table$scaffold)[n]
  # Get SNPs table for the scaffold only
  data_chr <- snps.table %>% filter(scaffold==chr)
 
  # Analyze data with GenWin - Produces a table and a plot
  spline <- splineAnalyze(data_chr$BF.dB., data_chr$position,
                          smoothness=100, plotRaw=T, plotWindows=T,  method=3)

  # Get the table with the results
  spline.data <- as.data.frame(spline$windowData)
  # Add the scaffold name to the results table
  spline.data$scaffold <- rep(chr,nrow(spline.data))
  # Add the variable
  spline.data$var <- rep(var,nrow(all_spline))
  # Add the Window length
  spline.data$WindowLength <- all_spline$WindowStop-all_spline$WindowStart

  # Alternate Odd and Even - for later plot of all scaffolds
  if((n %% 2) == 0) {
   spline.data$color <- rep("Even",nrow(spline.data))
  } else {
   spline.data$color <- rep("Odd",nrow(spline.data))
  } 

  # Fill the GenWin results table
  all_spline <- data.frame(rbind(all_spline, data.frame(spline.data)))
 }
 
 # Add the Window Number to spline dataset
 all_spline$WindowNumber <- 1:nrow(all_spline)
 
 # Calculate threshold for all genome
 all_spline.thresh=as.numeric(quantile(all_spline$Wstat,probs=0.999,na.rm=T))
 # Get whole-genome outliers table
 all_spline_total_outliers <- subset(all_spline, Wstat > all_spline.thresh)
 
 # Save the tables generated
 # All windows:
 write.table(x = all_spline,file = paste0("GenWin_results/",var,"_GenWin_windows.tsv"),
             quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")
 # Whole-Genome threshold outlier windows:
 write.table(x = all_spline_total_outliers,
             file = paste0("GenWin_results/",var,"_GenWin_windows_WholeGenomeOutliers.tsv"),
             quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")
}
```
Plot the results of GenWin for each variable (whole genome and single scaffolds)
```{R}
variables <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

for (n in 1:length(variables)){

 var=variables[n]
 
 # Upload GenWin results table
 all_spline=read.table(paste0("GenWin_results/",var,"_GenWin_windows.tsv"),h=T)
 # Add the variable
 all_spline$var <- rep(var,nrow(all_spline))
 # Add the Window length
 all_spline$WindowLength <- all_spline$WindowStop-all_spline$WindowStart
 # Add the Window Number
 all_spline$WindowNumber <- 1:nrow(all_spline)
 
 # Calculate Outlier value
 all_spline.thresh=as.numeric(quantile(all_spline$Wstat,probs=0.999,na.rm=T))
  
 # Manhattan plot of window results for whole genome:
 manhplot <- ggplot(all_spline, aes(x = WindowNumber, y = Wstat, color=color)) +
     geom_point(alpha = 0.75, stat = "identity", size=1) +
     scale_color_manual(values= c("Black","darkgray")) +
     # TITLE
     # LEGEND
     # X AXIS SCAFFOLD NAMES BELOW
     geom_hline(yintercept=all_spline.thresh,linetype="dashed", size=0.5, color="red")
 manhplot
 # Save the plot
 filename <- paste("GenWin_results/allspline_", var, "_WStat_manhattan.pdf", sep="")
 ggsave(filename)
 
 # Extract table of outlier windows
 all_spline_total_outliers <- subset(all_spline, Wstat > all_spline.thresh)

 # List consecutive outlier windows
 consec.windows <- grep("-",unlist(strsplit((seqToHumanReadable(all_spline_total_outliers$WindowNumber)),", ")),value=TRUE)
 # Create empty table to add consecutive outlier windows
 consec.w.table <- data.frame()
 # Fill consecutive outlier windows table with entries from general table
 if (length(consec.windows) > 0) { 
   for (w in 1:length(consec.windows)) {
    win=consec.windows[w]
    consec.w.length <- length(unlist(strsplit(win,"-"))[1]:unlist(strsplit(win,"-"))[2])
    for (k in unlist(strsplit(win,"-"))[1]:unlist(strsplit(win,"-"))[2]) {
     c.w.t.entry <- subset(all_spline_total_outliers, WindowNumber == k)
     consec.w.table <- data.frame(rbind(consec.w.table, data.frame(c.w.t.entry)))
    }
   }
 }
 
 # Write consecutive W table
 write.table(x = consec.w.table,
             file = paste0("GenWin_results/",var,"_consec_windows.tsv"),
             quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

 # Scaffolds:
 scaffolds <- levels(unique(all_spline$scaffold))
 
 # Draw bubbly plot for each scaffold
 for (s in 1:length(scaffolds)){
 
  chr=scaffolds[s]
  chr_spline <- subset(all_spline, scaffold == chr)
 
  # Manhattan plot of window results for scaffold bubbly:
  manhplot <- ggplot(chr_spline, aes(x = WindowNumber, y = Wstat)) +
     geom_point(alpha = 0.75, stat = "identity", size=chr_spline$SNPcount*0.1) +
     # TITLE
     # LEGEND
     # X AXIS SCAFFOLD NAMES BELOW
     geom_hline(yintercept=all_spline.thresh,linetype="dashed", size=0.5, color="red")
  manhplot
  # Save the scaffold plot
  filename <- paste("GenWin_results/", var, "_", chr ,"_bubble_WStat_manhattan.pdf", sep="")
  ggsave(filename)
  
  step_plot <- ggplot(chr_spline, aes(x = WindowLength, y = Wstat, order = WindowNumber)) +
      geom_segment() + 
      geom_hline(yintercept=all_spline.thresh,linetype="dashed", size=0.5, color="red")

 }
}
```
See how many outlier SNPs (BF>20) are present in each outlier window & other exploratory stuff
```{R}
for (n in 1:length(variables)){

 var=variables[n]
 
 # Upload GenWin results outliers table for variable
 spline_outliers=read.table(paste0("BayPass_plots&tables/",var,
                              "_GenWin_windows_WholeGenomeOutliers.tsv"),h=T)
 # Upload SNPs outliers table for variable
 snps_outliers=read.table(paste0("BayPass_plots&tables/",var,
                              "_outliers_SNPs.tsv"),h=T)
 
 # Create empty dataframe for new outlier window table with n of outliers snps for each window
 spline_outliers_nsnps <- data.frame()
 # Loop all windows one at the time
 for (w in 1:nrow(spline_outliers)){
   # Get one window
   window <- data.frame(spline_outliers[w,])
   # Get outlier SNPs inside that window
   snps.in.window <- subset(snps_outliers,
          position >= window$WindowStart & position <= window$WindowStop & scaffold == as.character(window$scaffold))
   # Add the number of outlier SNPs to the window table
   window$NoutlierSNPs=nrow(snps.in.window)
   # Add window entry to dataframe
   spline_outliers_nsnps <- data.frame(rbind(spline_outliers_nsnps, window))

 }
 assign(paste0(var,"_winoutliers_nsnps"), data.frame(spline_outliers_nsnps))
}





table.WG=read.table(paste0("BayPass_plots&tables/",var,"_GenWin_windows_WholeGenomeOutliers.tsv"),h=T)
table.WG$var <- rep(var,nrow(table.WG))
assign(paste0(var,"_WGoutliers_data"), data.frame(table.WG))

table.CHR=read.table(paste0("BayPass_plots&tables/",var,"_GenWin_windows_PerChrOutliers.tsv"),h=T)
table.CHR$var <- rep(var,nrow(table.CHR))
assign(paste0(var,"_CHRoutliers_data"), data.frame(table.CHR))

# rows in common wg vs chr
nrow(inner_join(table.WG, table.CHR))

# rows exclusive to wg
anti_join(table.WG, table.CHR)
# rows exclusive to chr
anti_join(table.CHR, table.WG)



for (n in 1:length(levels(unique(table.CHR$scaffold)))){
  chr=levels(unique(table.CHR$scaffold))[n]
  data_chr <- table.CHR %>% filter(scaffold==chr)
  assign(paste0(chr,"_perChr_data"), data.frame(data_chr))
  
}

```
