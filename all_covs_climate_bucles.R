library(sp)
library(raster)
library(gstat)
library(automap)
library(dismo)

# Load data and eliminate erroneous values and correct the variable class Year and Month
dat <- read.delim2("Estaciones Meteorologicas.csv",header = TRUE)
dat$A?o <- factor(dat$A?o)
dat$Mes <- factor(dat$Mes)
dat$T.Maxima[order(dat$T.Maxima, decreasing = T)][1:25]
dat$T.Maxima[dat$T.Maxima > 45] <- NA

# Estimate the average monthly value of each meteorological station, for each climatic variable
dat_agg_tmean <- aggregate(x= dat$T.Media, by=list(dat$Estacion, dat$Mes),FUN=mean,na.rm=TRUE)
dat_agg_tmin <- aggregate(x= dat$T.Minima, by=list(dat$Estacion, dat$Mes),FUN=mean,na.rm=TRUE)
dat_agg_tmax <- aggregate(x= dat$T.Maxima, by=list(dat$Estacion, dat$Mes),FUN=mean,na.rm=TRUE)
dat_agg_rain <- aggregate(x= dat$Precipitacion, by=list(dat$Estacion, dat$Mes),FUN=mean,na.rm=TRUE)

# Create SpatialPointsDataFrame of the Weather Stations
dat_sp <- unique(dat[, c("Estacion", "Latitud", "Longitud")])
coordinates(dat_sp) <- ~ Longitud + Latitud
dat_sp@proj4string <- CRS("+init=epsg:4326")

# Load the Digital Elevation Model, coordinate transformation of the SpatialPointsDataFrame 
# and generate a grid for use in interpolation by kriging method
DEM_32719 <- raster("covs/DEM/DEM_32719.tif")
dat_sp<- spTransform(dat_sp, DEM_32719@crs)
dat_sp$DEM_32719<- extract(DEM_32719,dat_sp)
DEM.sp <- as(DEM_32719, "SpatialGridDataFrame")

# Each loop performs spatial interpolation using the kriging model for each climate variable: 
# Tmean, Tmax, Tmin and Rain, using the correlation of data with respect to height (DEM).
# The generated data are saved in lists
# Use the 3 column to extract the data of the month "i",
# adding it to the column "i + 1" with the name of "x",
# Done in kirging using the column "x" and then modify the name of that column by the name Mes_ "i"

###################              Loop   Tmean        ##################################

kriging_tmean <- list()
raster_tmean <- list()
start_time_tmean<- Sys.time()

for(i in 1:12){
  dat_sp@data[,i+2] <- dat_agg_tmean[dat_agg_tmean$Group.2 == i,][,3]
  names(dat_sp@data)[i+2] <- "x"
  print(paste0("iteration #", i))
  kriging_tmean[[i]] <- autoKrige(formula = x ~ DEM_32719, input_data = dat_sp, new_data = DEM.sp, fix.values = c(0,20191,0.59))
  names(dat_sp@data)[i+2] <- paste0("mes_", i)
  raster_tmean[[i]] <- raster(kriging_tmean[[i]]$krige_output["var1.pred"])
}


save(kriging_tmean, file="covs/weather_results/modelos_kriging_tmean.RData")
# do.call generates an stack of all raster layers in one raster
tmean <- do.call(stack, raster_tmean)
writeRaster(tmean,"covs/weather_results/tmean.tif",format="GTiff")
names(tmean) <- month.name
plot(tmean)

end_time_tmean <- Sys.time()

end_time_tmean - start_time_tmean

###################              Loop  tmax           ##################################

kriging_tmax <- list()
raster_tmax<- list()
start_time_tmax <- Sys.time()

for(i in 1:12){
  dat_sp@data[,i+2] <- dat_agg_tmax[dat_agg_tmax$Group.2 == i,][,3]
  names(dat_sp@data)[i+2] <- "x"
  print(paste0("iteration #", i))
  kriging_tmax[[i]] <- autoKrige(formula = x ~ DEM_32719, input_data = dat_sp, new_data = DEM.sp, fix.values = c(0,20191,0.59))
  names(dat_sp@data)[i+2] <- paste0("mes_", i)
  raster_tmax[[i]] <- raster(kriging_tmax[[i]]$krige_output["var1.pred"])
}


save(kriging_tmax, file="covs/weather_results/modelos_kriging_tmax.RData")
tmax <- do.call(stack, raster_tmax)
writeRaster(tmax,"covs/weather_results/tmax.tif",format="GTiff")
names(tmax) <- month.name
plot(tmax)

end_time_tmax<- Sys.time()

end_time_tmax - start_time_tmax

###################              Loop  tmin           ##################################

kriging_tmin <- list()
raster_tmin <- list()
start_time_tmin <- Sys.time()

for(i in 1:12){
  dat_sp@data[,i+2] <- dat_agg_tmin[dat_agg_tmin$Group.2 == i,][,3]
  names(dat_sp@data)[i+2] <- "x"
  print(paste0("iteration #", i))
  kriging_tmin[[i]] <- autoKrige(formula = x ~ DEM_32719, input_data = dat_sp, new_data = DEM.sp, fix.values = c(0,20191,0.59))
  names(dat_sp@data)[i+2] <- paste0("mes_", i)
  raster_tmin[[i]] <- raster(kriging_tmin[[i]]$krige_output["var1.pred"])
}


save(kriging_tmin, file="covs/weather_results/modelos_kriging_tmin.RData")
tmin <- do.call(stack, raster_tmin)
writeRaster(tmin,"covs/weather_results/tmin.tif",format="GTiff")
names(tmin) <- month.name
plot(tmin)


end_time_tmin <- Sys.time()

end_time_tmin - start_time_tmin

###################              Loop  Rain           ##################################

kriging_rain <- list()
raster_rain<- list()
start_time_rain<- Sys.time()

for(i in 1:12){
  dat_sp@data[,i+2] <- dat_agg_rain[dat_agg_rain$Group.2 == i,][,3]
  names(dat_sp@data)[i+2] <- "x"
  print(paste0("iteration #", i))
  kriging_rain[[i]] <- autoKrige(formula = x ~ DEM_32719, input_data = dat_sp, new_data = DEM.sp, fix.values = c(0,20191,0.59))
  names(dat_sp@data)[i+2] <- paste0("mes_", i)
  raster_rain[[i]] <- raster(kriging_rain[[i]]$krige_output["var1.pred"])
}

save(kriging_rain, file="covs/weather_results/modelos_kriging_rain.RData")
rain <- do.call(stack, raster_rain)
writeRaster(tmax,"covs/weather_results/rain.tif",format="GTiff")
names(rain) <- month.name
plot(rain)

end_time_rain<- Sys.time()

end_time_rain - start_time_rain

###################              END LOOPS            ##################################
end_time_rain - start_time_tmean 

# Load of stacks
rain_raster<-stack("covs/weather_results/rain.tif")
tmean_raster<-stack("covs/weather_results/tmean.tif")
tmax_raster<-stack("covs/weather_results/tmax.tif")
tmin_raster<-stack("covs/weather_results/tmin.tif")

# Calculation of biovars
biovars<-biovars(rain_raster,tmin_raster,tmax_raster)
biovars
writeRaster(biovars,"results/biovars.tif",format="GTiff")
names(biovars)<-c("biovars.1","biovars.2","biovars.3","biovars.4","biovars.5","biovars.6","biovars.7","biovars.8","biovars.9","biovars.10","biovars.11","biovars.12","biovars.13","biovars.14","biovars.15","biovars.16","biovars.17","biovars.18","biovars.19")
save(biovars,file="covs/biovars_results/biovars.RData")
