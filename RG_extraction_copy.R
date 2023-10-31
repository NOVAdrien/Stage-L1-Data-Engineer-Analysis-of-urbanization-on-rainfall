# RAIN GAUGES ANALYSIS #

#----------------------------------------------------------------------------------------------------------------------------------
# CHARGE  LIBRARIES #
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
#----------------------------------------------------------------------------------------------------------------------------------
# SET WORKING DIRECTORY #
DIR = dirname(rstudioapi::getSourceEditorContext()$path) # Folder of the R scrit
setwd(DIR) # set working directory

#----------------------------------------------------------------------------------------------------------------------------------
# OPTIONS TO RUN #
THETA <- 45 # angle defining the zone
pressure <- 850 # Air pressure
Rmin <- 0  # distance between the raingauge and the zone (recommended 0-30km)
Rmax <- 20  # distance between the raingauge and end of the zone (recommended 20-60km)
start_year = 2002
end_year = 2005
land = 'US'
#city = 'Atlanta'

#----------------------------------------------------------------------------------------------------------------------------------
# CHARGE THE raingauge DATA #
# Put your data in dataframe format (...)

# Here stations selection
stations <- ncdc_stations(datasetid = 'PRECIP_HLY', 
                          startdate = paste0(start_year, '0101'), 
                          enddate = paste0(end_year, '1231'), 
                          extent = c(32.3,-95,40.6,-83.4), 
                          token = "ynaiHGszPBanriYXUYYVzhEdkuOVyJqe",
                          limit = 1000) # 1000 is upper limit

#----------------------------------------------------------------------------------------------------------------------------------
# PLOT YOUR raingauge IN A MAP #
RG_set <- matrix(ncol = 3, nrow = length(stations$data$id))
RG_set[,1] <- stations$data$id
RG_set[,2] <- stations$data$longitude
RG_set[,3] <- stations$data$latitude
colnames(RG_set) <- c('id', 'longitude', 'latitude')
RG_set <- data.frame(RG_set)
RG_set$longitude <- as.numeric(RG_set$longitude)
RG_set$latitude <- as.numeric(RG_set$latitude)
RG_spatial <- SpatialPointsDataFrame(coords = RG_set[,c("longitude", "latitude")],
                                  data = RG_set[,c("longitude", "latitude")],
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))
mapview(RG_spatial)

#Import wnd data
u_all <- brick(paste0("D:/Documents/0Polytech Sorbonne/S3/Stage/Analyses/Données/", land, " Uall ", start_year, " ", end_year, ".nc"))
v_all <- brick(paste0("D:/Documents/0Polytech Sorbonne/S3/Stage/Analyses/Données/", land, " Vall ", start_year, " ", end_year, ".nc"))

mapview(u_all[[1]])
mapview(v_all[[1]])
#----------------------------------------------------------------------------------------------------------------------------------
# SET CRS (coordinate reference systems) FOR THE WIND BASED ANALYSIS
crs_WRF   <- CRS("+proj=lcc +lat_1=28 +lat_2=50 +lat_0=39.7000122070312 +lon_0=-98 +x_0=0 +y_0=0 +a=6370000 +b=6370000 +units=m +no_defs")
crs_WGS <- crs("+proj=longlat +datum=WGS84 +no_defs ")

#----------------------------------------------------------------------------------------------------------------------------------
# ALL WIND DIRECTIONS FOR THE ANALYSIS (360 is better precision but take more time to compute)
windir <- seq(from = 1, to = 360, by = 1)
#----------------------------------------------------------------------------------------------------------------------------------
# THE LOOP  
system.time(
for (i_station in 1:length(stations$data$id)) {
  
  # CREATE THE DATAFRAME TO BE FILLED WITH THE RESULTS
  raingauge <- matrix(nrow = 1, ncol = 8)
  colnames(raingauge) <- c('ID', 'Date', 'precip', 'wind', 'upwind_min', 'upwind_max', 'urban', 'rural')
  raingauge <- data.frame(raingauge)
  
  # Location of the raingauge
  point_WGS <- SpatialPointsDataFrame(coords = RG_set[i_station,c("longitude", "latitude")],
                                      data = RG_set[i_station,c("longitude", "latitude")],
                                      proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # project un planar coordinates
  point_planar  <- spTransform(point_WGS,crs_WRF)
  
  # add a radius of Rmax
  Circle_LCC <- gBuffer(point_planar,width = Rmax*1000)
  # project to WGS coordinates
  Circle_WGS <- spTransform(Circle_LCC,crs_WGS)
  
  # view your raingauge with its radius of Rmax
  mapview(Circle_WGS)
  
  # IMPORT WIND DATA 
  #(download here, by year, by component : https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=form)
  vec_year <- seq(as.POSIXct(paste0(start_year,"-01-01 00:00:00")),
                  as.POSIXct(paste0(end_year,"-12-31 23:00:00")),
                  by ="1 hour")
  e <- extent(Circle_WGS)
  u = crop(u_all,e)
  v = crop(v_all,e)
  vecU = cellStats(u, stat='mean', na.rm=TRUE)
  vecV = cellStats(v, stat='mean', na.rm=TRUE)
  source("plot.windrose.R")
  vecW = sqrt(vecU^2 + vecV^2) # vitesse du vent
  
  df_wind = data.frame(Date = as.POSIXct(vec_year, format='%Y-%m-%d %H',tz='UTC'),  
                       u = vecU,
                       v = vecV,
                       ws = vecW,
                       wdir_rose = mod(atan2(vecU,vecV)*180/pi + 180,360),
                       wdir_to = mod(atan2(vecV,vecU)*180/pi + 180,360),
                       wdir_from = mod(mod(atan2(vecV,vecU)*180/pi + 180,360) + 180,360))
  
  df_wind$upwind_min = mod(df_wind$wdir_from-THETA/2,360)
  df_wind$upwind_max = mod(df_wind$wdir_from+THETA/2,360)
  
  png(file=paste0('Images/',substring(text = stations$data$id[i_station], first = 6),'_windrose.png'),
      width = 20,
      height = 20,
      units = "cm",
      res = 600)
  
  windrose <- plot.windrose(data = df_wind,
                            spd = "ws",
                            dir = "wdir_rose")
  print(windrose)
  
  dev.off()
  
  # GET NLCD (or charge yours)  also test 'impervious'
  ras_nlcd <- get_nlcd(template = Circle_WGS,
                       label = substr(stations$data$id[i_station], 6, 1000000L),
                       dataset = 'landcover',
                       force.redo = TRUE,
                       year = 2019)
  
  # Change values, 1 for urban, 0 for others (if you want to study other land use than urban put it at 1 and urban at 0)
  ras_nlcd[values(ras_nlcd) == 41| values(ras_nlcd) == 71| 
             values(ras_nlcd) == 43| values(ras_nlcd) == 52| 
             values(ras_nlcd) == 81| values(ras_nlcd) == 21| 
             values(ras_nlcd) == 11| values(ras_nlcd) == 82| 
             values(ras_nlcd) == 42| values(ras_nlcd) == 90|
             values(ras_nlcd) == 95| values(ras_nlcd) == 31] = 0
  ras_nlcd[values(ras_nlcd) == 22| values(ras_nlcd) == 23| 
             values(ras_nlcd) == 24] = 1
  
  # For Radius and angle raster, work on planar coordinates then reproject on wgs
 
  
   ras_nlcd_lcc  <- projectRaster(from = ras_nlcd,
                                 crs = point_planar@proj4string@projargs)
  
  # aggregate for faster run
  ras_nlcd_lcc <- aggregate(ras_nlcd_lcc, fact=100)
  
  png(file=paste0('Images/',substring(text = stations$data$id[i_station], first = 6),'_landuse.png'), width = 20, height = 20, units = "cm", res = 600)
  my_plot <- plot(ras_nlcd)
  print(my_plot)
  dev.off()
  
  png(file=paste0('Images/',substring(text = stations$data$id[i_station], first = 6),'_landuse_agg.png'), width = 20, height = 20, units = "cm", res = 600)
  my_plot <- plot(ras_nlcd_lcc)
  print(my_plot)
  dev.off()
  
  # pixel distance from center
  ras_Temp <- ras_nlcd_lcc
  ras_R <- distanceFromPoints(ras_Temp, coordinates(point_planar))/1000 #in km
  x <- init(ras_Temp, 'x') - coordinates(point_planar)[1]
  y <- init(ras_Temp, 'y') - coordinates(point_planar)[2]
  t <- atan2(y,x)/pi*180 +180# Angles range from 0 to 360 (same as wdir_to of df_wind table)
  ras_NA = ras_Temp
  values(ras_NA) = NA

  # create matrix for the 360 wind directions
  wind_360 <- matrix(nrow = 360, ncol = 2)
  colnames(wind_360) <- c('upwind_min', 'upwind_max')
  wind_360[,1] <- seq(from = 1, to = 360, by = 1)
  wind_360[,2] <- wind_360[,1] + THETA
  wind_360[(360-THETA):360,2] <- seq(from= 0, to = THETA, by= 1)
  wind_360 <- data.frame(wind_360)
  
  # we create the 360 mask, one for each wind direction
  b_upwind = stack(replicate(ras_NA,n=360)) 
  tmin = b_upwind
  values(tmin) = rep(wind_360$upwind_min,each = dim(b_upwind)[1]*dim(b_upwind)[2])
  tmax = b_upwind
  values(tmax) = rep(wind_360$upwind_max,each = dim(b_upwind)[1]*dim(b_upwind)[2])
  b_upwind[ ((tmin < tmax) & (t > tmin) & (t < tmax)) | ((tmin > tmax) & (((t > tmin) & (t < 360)) | ((t > 0) & (t < tmax))))] = 1
  
  b_upwindR = b_upwind # I need to add option to reduce radius
  
  # landuse with wind mask
  Res = b_upwindR*ras_nlcd_lcc
  
  # we calculate the urban and rural % for each wind direction
  for (k in 1:360) {
    wind_360$urban[k] <-  mean(values(Res[[k]]), na.rm =T)
    wind_360$rural[k] <- 1 - mean(values(Res[[k]]), na.rm =T)
  }
  
# must run by year because it is not possible to download more that a year of precipitation at the same time  
for (i_year in start_year:end_year) {
 
  # download precipitation data
  precip <- ncdc(datasetid='PRECIP_HLY', 
                 stationid = stations$data$id[i_station], 
                 datatypeid='HPCP',
                 startdate = paste0(i_year, '-01-01'), 
                 enddate = paste0(i_year, '-12-31'), 
                 token = "ynaiHGszPBanriYXUYYVzhEdkuOVyJqe", limit = 1000)
  
    if(length(precip$data$value) != 0 ) {
    
    precip <- precip$data
    
    if(length(precip$value) != 0 ) {
    
   # remove missing and 0 values
    precip <- precip[precip$value != 99999,]
    precip <- precip[precip$value != 0,]
    
    if(length(precip$value) != 0) {
  
    rg_data <- matrix(nrow = length(precip$date), ncol = 8)
    colnames(rg_data) <- c('ID', 'Date', 'precip', 'wind', 'upwind_min', 'upwind_max', 'urban', 'rural')
    rg_data <- data.frame(rg_data)
    rg_data$ID <- precip$station
    rg_data$Date <- as.POSIXct(precip$date, format = '%Y-%m-%dT%H', tz='UTC') 
    rg_data$precip <- precip$value 

    rain_date <- match(x = rg_data$Date, table = df_wind$Date)
    df_wind <- df_wind[rain_date,]
    
    rg_data$wind <- round(df_wind$wdir_to, digits = 0)
    rg_data$upwind_min <- round(df_wind$upwind_min, digits = 0)
    rg_data$upwind_max <- round(df_wind$upwind_max, digits = 0)
    rg_data$upwind_min[rg_data$upwind_min == 0] <- 360
    
    for (i_wind in 1:length(rg_data$ID)) {
      
      y_wind <- rg_data$upwind_min[i_wind]
      rg_data$urban[i_wind] <- wind_360$urban[y_wind]
      rg_data$rural[i_wind] <- wind_360$rural[y_wind]
      
      }

    raingauge <- rbind(raingauge, rg_data)
    rm('precip')
    
    } else { print("no precipitation during this period")}
    } else { print("no precipitation during this period")}
    } else { print("no precipitation during this period")}
 }
  write.table(raingauge, file = paste0(land,'/',substring(text = rg_data$ID[1], first = 6),".txt"), sep = ";")
}
  
)

