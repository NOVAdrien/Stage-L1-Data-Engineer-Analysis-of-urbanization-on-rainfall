# RAIN GAUGES Selection

# CHARGE  LIBRARIES
library(magrittr)
library(rgeos)
library(rgdal)
library(mapview)
library(tidyr)
library(rnoaa)
library(raster)
devtools::install_github("ropensci/FedData")
library(FedData)
library(ggmap)
library(ncdf4)
library(progress)

# SET WORKING DIRECTORY
DIR = dirname(rstudioapi::getSourceEditorContext()$path) # Folder of the R scrit
setwd(DIR) # set working directory

# OPTIONS TO RUN
start_year = 1980
end_year = 2010
land = 'US'
coord_zone = c(32.3,-95,40.6,-83.4)
Pas = 10 # Tested before not to exceed ncdc_stations function limit 1000
urb_min <- 5
urb_max <- 95
Rmax <- 40
Rmin <- 0
THETA <- 45

# Map system
crs_WRF   <- CRS("+proj=lcc +lat_1=28 +lat_2=50 +lat_0=39.7000122070312 +lon_0=-98 +x_0=0 +y_0=0 +a=6370000 +b=6370000 +units=m +no_defs")
crs_WGS <- crs("+proj=longlat +datum=WGS84 +no_defs ")

# CHARGE THE raingauge DATA
# Put your data in dataframe format (...)

start_date = as.POSIXct(paste0(start_year,"-01-01 00:00:00"))
end_date = as.POSIXct(paste0(end_year,"-12-31 23:00:00"))

# Creating matrix and sequence for uploading
my_coord <- seq(from = coord_zone[1], to = coord_zone[3], length.out = Pas)
my_pluv <- matrix(ncol = 9)
colnames(my_pluv) <- c("elevation", "mindate", "maxdate", "latitude", "name", 
                       "datacoverage", "id", "elevationUnit", "longitude" )

# Uploading stations data following the Pas
for (x in 1 : length(my_coord)) {
  
  my_lat_min <- my_coord[x]
  my_lat_max <- my_coord[x + 1]
  stations <- ncdc_stations(datasetid = 'PRECIP_HLY', 
                            extent = c(my_lat_min, coord_zone[2], my_lat_max, coord_zone[4]), 
                            token = "ynaiHGszPBanriYXUYYVzhEdkuOVyJqe",
                            limit = 1000) # 1000 is upper limit
  
  # Adding data above in the matrix
  my_pluv <- rbind(my_pluv, stations$data)
  
}

# Removing first empty line
my_pluv <- my_pluv[-1,]

# Extracting year from start & end date
my_pluv$mindate <- as.POSIXct(my_pluv$mindate)
my_pluv$maxdate <- as.POSIXct(my_pluv$maxdate)

# Stations selection with data periode
my_pluv_2 <- my_pluv[my_pluv$mindate<=start_date,]
my_pluv_2 <- my_pluv_2[my_pluv_2$maxdate>=end_date,]

# Test length pluv number
length(my_pluv_2[,1])

# Adding col "urbain"
my_pluv_3 <- my_pluv_2
my_pluv_3$urbain = 0

# Set timer
pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = length(my_pluv_3$id),
  clear = FALSE
)

# Creating map for each station and completing the data my_pluv_3
for (i_station in 56:length(my_pluv_3$id)) {
  
  # Location of the raingauge
  point_WGS <- SpatialPointsDataFrame(coords = my_pluv_3[i_station,c("longitude", "latitude")],
                                      data = my_pluv_3[i_station,c("longitude", "latitude")],
                                      proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # project un planar coordinates
  point_planar  <- spTransform(point_WGS,crs_WRF)
  
  # add a radius of Rmax
  Circle_LCC <- gBuffer(point_planar,width = Rmax*1000)
  
  # project to WGS coordinates
  Circle_WGS <- spTransform(Circle_LCC,crs_WGS)
  
  # view your raingauge with its radius of Rmax
  mapview(Circle_WGS)
  
  # Get data for stations
  ras_nlcd <- get_nlcd(template = Circle_WGS,
                       label = substr(my_pluv_3$id[i_station], 6, 1000000L),
                       dataset = 'landcover',
                       force.redo = FALSE,
                       year = 2019)
  
  # Set 0 for rural and 1 for urban
  ras_nlcd[values(ras_nlcd) == 41| values(ras_nlcd) == 71| 
             values(ras_nlcd) == 43| values(ras_nlcd) == 52| 
             values(ras_nlcd) == 81| values(ras_nlcd) == 21| 
             values(ras_nlcd) == 11| values(ras_nlcd) == 82| 
             values(ras_nlcd) == 42| values(ras_nlcd) == 90|
             values(ras_nlcd) == 95| values(ras_nlcd) == 31] = 0
  ras_nlcd[values(ras_nlcd) == 22| values(ras_nlcd) == 23| 
             values(ras_nlcd) == 24] = 1
  
  # Create vector with data
  ras_nlcd <- raster(ras_nlcd)
  
  ras_nlcd_lcc  <- projectRaster(from = ras_nlcd,
                                 crs = crs(point_planar))
  
  ras_nlcd_lcc <- aggregate(ras_nlcd_lcc, fact=100)
  my_pluv_3$urbain[i_station] = cellStats(ras_nlcd_lcc, stat = 'mean')
  
  # Update the progress bar
  pb$tick() 
  gc()
  
}

# Saving table with data
write.table(my_pluv_3, file = paste0(land, '_', start_year, '_', end_year, '_', urb_min, '%_', urb_max, '%_', Rmin, '-', Rmax, 'km'), sep = ";")
