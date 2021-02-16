---
title: "Climatic_Variables_biomes"
author: "Enrico"
date: "20 January 2021"
output: html_document
editor_options:
  chunk_output_type: console
---
In this markdown I will extract the Biome of each population from a shapefile downloaded from:
(https://ecoregions2017.appspot.com/)
(https://academic.oup.com/bioscience/article/67/6/534/3102935)

Prepare R
```{R}
library(tidyverse)
library(raster)
library(rgdal)
library(geobuffer)
library(grDevices)
library(rgeos)
library(sf)
```
Input and prepare the biomes layer and population coordinates
```{R}
biomes_shape <- readOGR(paste0("~/Downloads/Ecoregions2017/Ecoregions2017.shp"))
lala <- biomes_shape@data
biome_numbers_table <- data.frame(name=unique(lala$BIOME_NAME)[-15],num=unique(lala$BIOME_NUM))
biomes <- biomes_shape[3]
head(biomes, rows = 100)
r <- raster()
extent(r) <- extent(biomes)
biome_raster <- rasterize(biomes, r, 'BIOME_NUM')
plot(biome_raster)

coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/csv_LL_selection_coords_wholeset.csv", col_names = T, delim = ',') 
# define populations
samples <- coord_table$id %>% unique
elements_2_remove <- c("c_ll_cr_0211", "c_ll_ba_0233")
samples = samples[!(samples %in% elements_2_remove)]

table <- data.frame()
for (i in 1:length(samples)){
  # define sample
  sample <- samples[i]
  # get population coordinates
  coord_sample <- coord_table %>% filter(coord_table$id == sample)
  # in a dataframe with columns x and y
  coords <- data.frame(x=as.numeric(coord_sample$longitude),
                       y=as.numeric(coord_sample$latitude))
  # get buffer around coordinates
  buff <- geobuffer_pts(xy = coords, dist_m = 100000)
  buffer <- buff[1] # first layer only (though I guess with 1 sample it's useless)

  row <- data.frame(sample=coord_sample$id)
  # Get mean of biovariable value of all points within the polygon
  values <- extract(
      biome_raster,                  # raster layer
      buffer,               # SPDF with centroids for buffer
      na.rm = TRUE,         # remove NAs
      df=TRUE)              # return a dataframe
  
  b <- data.frame(biome_n=unique(values$layer))
  # add the mean to the population table
  row <- cbind(row, b)
  table <- rbind(table, row)
}
table_nona <- na.omit(table)
write.table(table_nona, "biomes_persample.tsv", sep = "\t", row.names = F)
write.table(biome_numbers_table, "biome_numbers_table.tsv", sep = "\t", row.names = F)

# Version of this table where only most present biome is kept for each sample
table2 <- data.frame()

for (i in 1:length(samples)){
  # define sample
  sample <- samples[i]
  # get population coordinates
  coord_sample <- coord_table %>% filter(coord_table$id == sample)
  # in a dataframe with columns x and y
  coords <- data.frame(x=as.numeric(coord_sample$longitude),
                       y=as.numeric(coord_sample$latitude))
  # get buffer around coordinates
  buff <- geobuffer_pts(xy = coords, dist_m = 100000)
  buffer <- buff[1] # first layer only (though I guess with 1 sample it's useless)

  row <- data.frame(sample=coord_sample$id)
  # Get mean of biovariable value of all points within the polygon
  values <- extract(
      biome_raster,                  # raster layer
      buffer,               # SPDF with centroids for buffer
      na.rm = TRUE,         # remove NAs
      df=TRUE)              # return a dataframe

  b <- data.frame(biome_n=tail(names(sort(table(values$layer))), 1))
  # add the mean to the population table
  row <- cbind(row, b)
  table2 <- rbind(table2, row)
}
write.table(table2, "biomes_persample_majority.tsv", sep = "\t", row.names = F)



limits <- raster(buff)@extent %>% as.list # get xmin, xmax, ymin , ymax
cropbox <- c(limits)
DEMcrop <- crop(biome_raster, cropbox)
plot(DEMcrop, legend=T, axes=FALSE, main=paste0(sample[i],"_",names(biome_raster)))
points(coords, pch = 19, add=T)
plot(buffer, add=T)
```
