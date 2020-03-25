# I will use this script to get environmental data for my different populations
# from the WorldClim database.

library(raster)
library(readr)
library(dplyr)

# Load data from WorldClim
worldclim <- getData("worldclim", var = "bio", res = 10)

#################################
## Example with bio1 and bio12 ##
#################################
bios <- worldclim[[c("bio1", "bio12")]]
names(bios) <- c("Temp","Prec")
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/CSV_LL_selection_coordinates_noNA.csv", col_names = T, delim = ',')
coords <- data.frame(x=coord_table$longitude,y=coord_table$latitude)
points <- SpatialPoints(coords, proj4string = bios@crs)
values <- extract(bios,points)
df <- cbind.data.frame(coordinates(points),values)
# Plot although it is weird (points get displaced if plot is changed of proportions)
plot(bios[[5]])
plot(points,add=T)

####################################
## Extracting population avarages ##
####################################
# with all Worldclim variables for YAKUTIA population ("ya") it would be
bios <- worldclim[[c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")]]
names(bios) <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/CSV_LL_selection_coordinates_noNA.csv", col_names = T, delim = ',')
coord_pop <- coord_table %>% filter(coord_table$pop == "ya")
coords <- data.frame(x=coord_pop$longitude,y=coord_pop$latitude)
points <- SpatialPoints(coords, proj4string = bios@crs)
values <- extract(bios,points)
df <- cbind.data.frame(coordinates(points),values)

# in a loop for all populations
bios <- worldclim[[c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")]]
names(bios) <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/CSV_LL_selection_coordinates_noNA.csv", col_names = T, delim = ',')
populations <- coord_table$pop %>% unique
table <- data.frame(WCvar=c("x","y","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"))
for (i in 1:length(populations)){
  coord_pop <- coord_table %>% filter(coord_table$pop == populations[i])
  coords <- data.frame(x=coord_pop$longitude,y=coord_pop$latitude)
  points <- SpatialPoints(coords, proj4string = bios@crs)
  values <- extract(bios,points)
  df <- cbind.data.frame(coordinates(points),values)
  means <- data.frame(mean=colMeans(df, na.rm = T), row.names = NULL)
  colnames(means) <- populations[i]
  table <- cbind(table,means)
}

write.table(x = table,file = "WorldClim_biovariables.txt",quote=FALSE, col.names = T, row.names = FALSE, sep= "\t")


## POLYGONS

# with all Worldclim variables for YAKUTIA population ("ya") it would be
bios <- worldclim[[c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")]]
names(bios) <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/CSV_LL_selection_coordinates_noNA.csv", col_names = T, delim = ',')
coord_pop <- coord_table %>% filter(coord_table$pop == "ya")
coords <- data.frame(x=coord_pop$longitude,y=coord_pop$latitude)
points <- SpatialPoints(coords, proj4string = bios@crs)
values <- extract(bios,points)
df <- cbind.data.frame(coordinates(points),values)
