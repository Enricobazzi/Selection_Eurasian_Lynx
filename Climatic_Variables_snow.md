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
library(sf)
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
To chek how populations are similar to each other for the snow variables:
```{R}
library(RColorBrewer)
library (viridis)

cols <- c("ca"="#B8860b",
          "ur"="#0F4909", 
          "ki"=viridis_pal()(5)[1], 
          "la"=brewer.pal(12,"Paired")[3], 
          "tu"=brewer.pal(12,"Paired")[8], 
          "mo"=brewer.pal(12,"Paired")[7], 
          "vl"=brewer.pal(12,"Paired")[5], 
          "ya"=brewer.pal(12,"Paired")[6])


snow_data <- read_delim("Snow_table.tsv", col_names = T, delim = '\t')
snow_data_long <- snow_data%>%pivot_longer(cols=ca:ya, names_to="pop")

jan_depth <- snow_data_long %>% filter(variable=="jan_depth")%>%
  mutate(pop=factor(pop, levels=c("ca","ki","la","ur","mo","tu","vl","ya")))

ggplot(jan_depth, aes(x=pop, y=value, color=pop))+
  geom_point(size=10, stat='identity')+
  scale_color_manual(values=cols)+
  theme_bw()

ggsave("jan_depth_values_dots.pdf", path = "PCA_outliers/",
        width=25,height=25,units="cm")

snow_days <- snow_data_long %>% filter(variable=="snow_days")%>%
  mutate(pop=factor(pop, levels=c("ca","ki","la","ur","mo","tu","vl","ya")))
ggplot(snow_days, aes(x=pop, y=value, color=pop))+ 
  geom_point(size=10, stat='identity')+
  scale_color_manual(values=cols)+
  theme_bw()
ggsave("snow_days_values_dots.pdf", path = "PCA_outliers/",
        width=25,height=25,units="cm")
```

In order to run RDA using individual samples genotypes and I will generate a table with avarage values of each snow variable for a buffer around the coordinates of each sample.

To prepare R for this:
```{R}
library(raster)
library(tidyverse)
library(grDevices)
library(rgeos)
library(geobuffer)
library(rgdal)
library(sf)
```
I also prepared a table with the coordinates of each sample. Samples without coordinates (cr_0211 and ba_0233) will not be included in the RDA. I will also remove the "bad" Balkans samples which will also not be included in the RDA (c_ll_ba_0216, h_ll_ba_0214 and h_ll_ba_0215).

Using st_transform from the package sf I can convert my buffer's (spatial polygon) coordinate system to the snow data's coordinate system directly in R (no need to pass through qGis like above). I will generate a table with avarage snow depth of each sample for each day from 1999 to 2019 (my snow data).
```{R}
# Load coordinates
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/csv_LL_selection_coords_wholeset.csv", col_names = T, delim = ',') 
# define samples
samples <- coord_table$id %>% unique
# remove bad samples and ones without coordinates
elements_2_remove <- c("c_ll_cr_0211", "c_ll_ba_0233","c_ll_ba_0216", 
                       "h_ll_ba_0214", "h_ll_ba_0215")
samples = samples[!(samples %in% elements_2_remove)]

# Loop through all samples:
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
  # transform polygon to snow data coordinates system
  buffer2 <- st_transform(as(buffer, "sf"), crs = "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=10 +k=1 +x_0=0 +y_0=0 +a=6371200 +b=6371200 +units=m +no_defs")
 
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
                   buffer2,             # polygon of population
                   fun=mean,             # get the mean
                   na.rm = TRUE,         # remove NAs
                   df=TRUE)              # return a dataframe
 
   # add the mean to the population table
   snow_level <- means[[2]][1]
   day_data <- data.frame(sample=samples[i], year=y, day=k, mean_snow=snow_level)
   TABLE <- data.frame(rbind(TABLE,day_data))
     }
   }
   write.table(x = TABLE, 
               file = paste0("Snow_data_daily/",samples[i],"_snow_daily_means.tsv"),
               quote=FALSE, col.names = T, row.names = FALSE, sep= "\t")
}
```
Now to get per sample snow_days and jan_depth:
```{R}
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/csv_LL_selection_coords_wholeset.csv", col_names = T, delim = ',')
# define populations
samples <- coord_table$id %>% unique
# remove bad samples and ones without coordinates
elements_2_remove <- c("c_ll_cr_0211", "c_ll_ba_0233","c_ll_ba_0216", 
                       "h_ll_ba_0214", "h_ll_ba_0215")
samples = samples[!(samples %in% elements_2_remove)]

snow_data <- data.frame()
for (i in 1:length(samples)){
  # define population
  sample=samples[i]
  # load population daily snow data
  sam_snow_daily <- read_delim(paste0("Snow_data_daily/",sample,"_snow_daily_means.tsv"), 
                               col_names = T, delim = '\t')
  sam_snow_yearly <- data.frame()
  for (y in 1999:2019){
    # january mean depth for year
    sam_snow_january <- subset(sam_snow_daily, day<32 & year==y)
    jan_mean_row <- data.frame(year = y, c1 = mean(sam_snow_january$mean_snow))
    # number of snow days for year
    sam_snow_days <- subset(sam_snow_daily, mean_snow>0 & year==y)
    snow_days_row <- data.frame(c2 = length(sam_snow_days$mean_snow))
    sam_snow_yearly <- data.frame(rbind(sam_snow_yearly, cbind(jan_mean_row, snow_days_row)))
  }
  snow_data <- data.frame(rbind(snow_data, 
                                data.frame(sample=sample,jan_depth=mean(sam_snow_yearly$c1),
                                                 snow_days=mean(sam_snow_yearly$c2))))
}
write.table(x = snow_data,
            file = "Snow_table_persample.tsv", quote=FALSE,
            col.names = T, row.names = FALSE, sep= "\t")
```
