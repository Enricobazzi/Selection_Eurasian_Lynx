---
title: "Climatic_Variables"
author: "Enrico"
date: "31 March 2020"
output: html_document
editor_options:
  chunk_output_type: console
---

This is the steps I took in order to extract values of the bioclimatic variables of WorldClim for each Lynx lynx population.

The objective is to have a single value for each variable for each population in a table that can be used as input while running the software BayPass.

The strategy we adopted was to draw a polygon enclosing all of the coordinates of the samples in the population, and extract the mean value for all of the points in that area.

These are the libraries I used

```{R}
library(raster)
library(readr)
library(dplyr)
library(grDevices)
library(ggplot2)
library(rgeos)
#devtools::install_github("valentinitnelav/geobuffer") # installed with devtools for buffer creation
library(geobuffer)
```

I imported the data for all biovariables from the WorldClim database (resolution 10m -> ~340km2 cells), and my table with coordinates and population data for each individual.

```{R}
# first load all environmental data:
worldclim <- getData("worldclim", var = "bio", res = 10)
# and all coordinates data:
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/csv_LL_selection_coords_finalset.csv", col_names = T, delim = ',')
# define populations
populations <- coord_table$pop %>% unique
```

I have one population for which all individuals are registered to the same coordinates. I will eliminate it from the first extraction of values and deal with it later:

```{R}
# remove populations (tuva as it only has 1 coordinate for all samples)
elements_2_remove <- "tu"
populations = populations[!(populations %in% elements_2_remove)]
```

I will loop through the rest of the populations, creating the polygon enclosing its coordinates, extracting mean values for each bio-variable for cells within the polygon, and drawing a plot of the area with cell values and polygon borders.

The mean values will be stored in a data frame which is created empty before the loop. The graphs of each population are saved in the folder WorldClim_graphs

```{R}
# create a table with 1 column to bind to the rest of the calculated data
table <- data.frame(rbind(data.frame(C1="pop"),data.frame(C1=names(worldclim))))
# now lets loop though populations:
for (i in 1:length(populations)){
  # get population coordinates
  coord_pop <- coord_table %>% filter(coord_table$pop == populations[i])
  # in a dataframe with columns x and y
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

  # Now we loop through all biovariables to extract their values within our polygon:
  pop_table <- data.frame(pop=populations[i], stringsAsFactors = F)

  # set up divided plot with no margins:
  jpeg((paste0("WorldClim_graphs/graph_",populations[i],".pdf")))
  par(mfrow=c(5,4), mar = c(0, 2, 2, 0))

  for (k in 1:nlayers(worldclim)){

      # subset for just one biovariable
      DEM <- worldclim[[k]]
      names(DEM) # print the name of the biovariable

      # Get mean of biovariable value of all points within the polygon
      means <- raster::extract(
                  DEM,                  # raster layer
                  SPDF,                 # SPDF with centroids for buffer
                  fun=mean,             # get the mean
                  na.rm = TRUE,         # remove NAs
                  df=TRUE)              # return a dataframe

      # add the mean to the population table
      row <- data.frame(means[2])
      colnames(row) <- "pop"
      pop_table <- rbind(pop_table, row)

      # plot all the biovariables and the polygon for the population:
      limits <- raster(SPlgs)@extent %>% as.list # get xmin, xmax, ymin , ymax
      cropbox <- c(limits)
      DEMcrop <- crop(DEM, cropbox)
      plot(DEMcrop, legend=F, axes=FALSE, main=paste0(populations[i],"_",names(DEM)))
      points(coords, pch = 19, add=T)
      plot(SPDF, add=T)
  }
  # add the population table to the general table
  table <- cbind(table,pop_table)
  # dev off to reset plot
  dev.off()
}
```

For Tuva I need to do a buffer around the single location I have coordinates for.

I wanted to know the size of the polygons I have for the other populations to decide the final size of the buffer around Tuva.

```{R}
for (i in 1:length(populations)){
  # get population coordinates
  coord_pop <- coord_table %>% filter(coord_table$pop == populations[i])
  # in a dataframe with columns x and y
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
  print(paste0(populations[i]," ",gArea(SPDF)))

}
```

There is high variability of polygon area:
Ki has 0.19 with a very narrow triangle shape,while the biggest area is Yakutia (which has a far away point), with 18.39.
Others values are: 2-5 (mo,ur,la); 8-11 (vl,ca).

I import the data for Tuva only first:

```{R}
# I will work separetly on Tuva population only:
coord_tuva <- coord_table %>% filter(coord_table$pop == "tu")
coords <- data.frame(x=coord_tuva$longitude,y=coord_tuva$latitude)
points <- SpatialPoints(coords, proj4string = worldclim@crs)
```

I will try getting an area of around ~4 for Tuva -> 100km around the point

```{R}
# create a circular buffer of 100km around - coords is now Tuva (see above for import of Tuva only)
buff <- geobuffer_pts(xy = coords, dist_m = 100*10^3)
buffer <- buff[1] # this is because buff has 1 layer for each sample but they are all the same
gArea(buffer) # ~4
```

Now to extract the mean values for the different variables for that area I loop through the variables, extract mean in a data frame (again created empty before the loop) and the graph is saved in the folder WorldClim_graphs

```{R}
# Now we loop through all biovariables to extract their values within our buffer area:
tuva_table <- data.frame(pop="tu", stringsAsFactors = F)

# set up divided plot with no margins:
jpeg((paste0("WorldClim_graphs/graph_tu.pdf")))
par(mfrow=c(5,4), mar = c(0, 2, 2, 0))

for (k in 1:nlayers(worldclim)){

  # subset for just one biovariable
  DEM <- worldclim[[k]]
  names(DEM) # print the name of the biovariable

  # Get mean of biovariable value of all points within the polygon
  means <- raster::extract(
    DEM,                  # raster layer
    buffer,                 # SPDF with centroids for buffer
    fun=mean,             # get the mean
    na.rm = TRUE,         # remove NAs
    df=TRUE)              # return a dataframe

  # add the mean to the population table
  row <- data.frame(means[2])
  colnames(row) <- "pop"
  tuva_table <- rbind(tuva_table, row)

  # plot all the biovariables and the polygon for the population:
  limits <- raster(buffer)@extent %>% as.list # get xmin, xmax, ymin , ymax
  cropbox <- c(limits)
  DEMcrop <- crop(DEM, cropbox)
  plot(DEMcrop, legend=F, axes=FALSE, main=paste0("tu_",names(DEM)))
  points(coords, pch = 19, add=T)
  plot(buffer, add=T)

}
```

Now I to add Tuva to the table of the rest of populations and export the final table:

```{R}
# add Tuva column:
table <- cbind(table, tuva_table)

# export final table
write.table(x = table,file = "WorldClim_table.tsv",quote=FALSE, col.names = F, row.names = FALSE, sep= "\t")
```

I exported it as a TSV for easier manipulation in excel. From there I will:

- change the decimals to only 2;
- change the order of the columns so that they are in alphabetical order (as in allele counts file: ca - ki - la - mo - tu - ur - vl - ya).

Then I will add the additional environmental data sent by Krystoff. This will be done manually as he sent the data in another table. This data is regarding prey type ...

Finally I can copy the final table to the BayPass folder in the EBD genomics server:
```
scp WorldClim_table.tsv ebazzicalupo@genomics-b.ebd.csic.es:~/BayPass/Covariate_Data/
```
