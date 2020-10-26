---
title: "Climatic_Variables_snow"
author: "Enrico"
date: "28 September 2020"
output: html_document
editor_options:
  chunk_output_type: console
---

In this markdown I will download, process and analyze Snow data in order to obtain a value of snow depth and snow cover days for each Lynx population. These values will be then analyzed using BayPass to find association between environment and genotype.

Load the libraries to be used
```{R}
library(tidyverse)
library(raster)
library(rgdal)
library(geobuffer)
library(grDevices)
library(rgeos)
```
As the coordinate system of the snow data is different from the one of the samples, I convert the population polygons as shape files, import them into qGis and convert the coordinate system there (coulnd't find an easy way to do it in R although there probably is).
Create a shape file for each population's polygon
```{R}
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/samples_selection.csv", col_names = T, delim = ',')
# define populations
populations <- coord_table$pop %>% unique
# remove populations (tuva as it only has 1 coordinate for all samples)
elements_2_remove <- "Tuva"
populations = populations[!(populations %in% elements_2_remove)]

for (i in 1:length(populations)){
  # get population coordinates
  coord_pop <- coord_table %>% filter(coord_table$pop == populations[i])
  
  coords <- data.frame(x=coord_pop$longitude,y=coord_pop$latitude)
  
  # find convex hull (which points form a polygon enclosing the rest of the coordinates)
  con.hull.pos <- chull(coords) # find positions of convex hull
  
  # put it in a data frame to use
  con.hull <- rbind(coords[con.hull.pos,],coords[con.hull.pos[1],])
  
  # Create a Polygon object in Raster (spatial polygon dataframe) for convex hull
  P1 <- Polygon(con.hull)
  P1s <- Polygons(list(P1), 1)
  SPlgs <- SpatialPolygons(list(P1s))
  df <- data.frame(ID = 1)
  SPDF <- SpatialPolygonsDataFrame(SPlgs, df)
  shapefile(SPDF, filename=paste0("Snow_data_daily/",
                                  populations[i], "_SPDF_shape.shp"))
}

# ADD Tuva
coord_tuva <- coord_table %>% filter(coord_table$pop == "Tuva")
coords <- data.frame(x=coord_tuva$longitude,y=coord_tuva$latitude)
points <- SpatialPoints(coords, proj4string = worldclim@crs)
buff <- geobuffer_pts(xy = coords, dist_m = 100*10^3)
buffer <- buff[1] # this is because buff has 1 layer for each sample but they
gArea(buffer) # ~4
SPDF <- SpatialPolygonsDataFrame(buffer, df)
shapefile(SPDF, filename=paste0("Snow_data_daily/Tuva_SPDF_shape.shp"))
```
After converting the polygon coordinates I can import the new shapefiles into R and extract the avarage snow depth of each day of each year for each population
```{R}
world <- readOGR("/Users/enricobazzicalupo/Dropbox/LL_LC_LR_Databases/Snow_data_analysis/world_polar_stereographic/world_polar_stereographic.shp")

coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/samples_selection.csv", col_names = T, delim = ',')

# define populations
populations <- coord_table$pop %>% unique

# plot polygons to see if all is well
plot(world)
for (i in 1:length(populations)){
pop_poly <- readOGR(paste0("Snow_data_daily/",populations[i],"_nc.shp"))
plot(pop_poly, col = "red", add=T)
}
# all is well

# Loop through all years
for (i in 1:length(populations)){
 
 # get population polygon from shapefile
 pop_poly <- readOGR(paste0("Snow_data_daily/",populations[i],"_nc.shp"))
 
 TABLE <- data.frame()

 for (y in 1999:2019){
  # Get data for year
  data <- stack(paste0("Snow_data_daily/cmc_sdepth_dly_", y, "_v01.2.tif"))
 
  for (k in 1:nlayers(data)){
   # subset for just one month
   DEM <- data[[k]]
   # Get mean of value of all points within the polygon
   means <- raster::extract(
                   DEM,                  # raster layer
                   pop_poly,             # polygon of population
                   fun=mean,             # get the mean
                   na.rm = TRUE,         # remove NAs
                   df=TRUE)              # return a dataframe
 
   # add the mean to the population table
   snow_level <- means[[2]][1]
   day_data <- data.frame(pop=populations[i], year=y, day=k, mean_snow=snow_level)
   TABLE <- data.frame(rbind(TABLE,day_data))
     }
   }
   write.table(x = TABLE, file = paste0("Snow_data_daily/",populations[i],"_snow_daily_means.tsv"),
               quote=FALSE, col.names = T, row.names = FALSE, sep= "\t")
}
```
Looking at literature, Ratkiewicz et al. (2014) used avarage snow depth in january and days with snow cover as variables to correlate with genetic structure.

Extract a table with each population as a column, and two rows: (1) avarage snow depth in january; and (2) avarage yearly number of days with snow (depth>0)
```{R}
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/samples_selection.csv", col_names = T, delim = ',')
# define populations
populations <- coord_table$pop %>% unique

snow_data <- data.frame(variable=c("jan_depth","snow_days"))
for (i in 1:length(populations)){
  # define population
  population=populations[i]
  # load population daily snow data
  pop_snow_daily <- read_delim(paste0("Snow_data_daily/",population,"_snow_daily_means.tsv"), 
                               col_names = T, delim = '\t')
  pop_snow_yearly <- data.frame()
  for (y in 1999:2019){
    # january mean depth for year
    pop_snow_january <- subset(pop_snow_daily, day<32 & year==y)
    jan_mean_row <- data.frame(year = y, c1 = mean(pop_snow_january$mean_snow))
    # number of snow days for year
    pop_snow_days <- subset(pop_snow_daily, mean_snow>0 & year==y)
    snow_days_row <- data.frame(c2 = length(pop_snow_days$mean_snow))
    pop_snow_yearly <- data.frame(rbind(pop_snow_yearly, cbind(jan_mean_row, snow_days_row)))
  }
  snow_data <- data.frame(cbind(snow_data, 
                                data.frame(v1=(c(mean(pop_snow_yearly$c1), mean(pop_snow_yearly$c2))))))
  names(snow_data)[names(snow_data) == "v1"] <- population
}
write.table(x = snow_data,
            file = "Snow_table.tsv", quote=FALSE, col.names = T, row.names = FALSE, sep= "\t")
```
I exported it as a TSV for easier manipulation in excel. From there I will:

- change the decimals to only 2;
- change the order of the columns so that they are in alphabetical order (as in allele counts file: ca - ki - la - mo - tu - ur - vl - ya).

Finally I can copy the final table to the BayPass folder in the EBD genomics server for GEA analysis:
```
scp Snow_table.tsv ebazzicalupo@genomics-b.ebd.csic.es:~/BayPass/Covariate_Data/
```
