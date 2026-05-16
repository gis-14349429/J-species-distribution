#Set working directory
setwd("C:/Users/User/OneDrive - The University of Manchester/GEOG71922/Assessment-2/14349429_A2")

#Install packages
install.packages(c("terra","gdistance","raster","igraph"))
install.packages(c("sf"))
library(terra)
library(raster)
library(sf)
#Read in the species data
beetle<- read.csv("scot_beetle_community.csv")

#view the first six rows of the data
head(beetle) 

#Load land cover for the analysis
UK=rast("LCMUK_2000.tif")
#Plot the landcover data
plot(UK)
env<-read.csv("scot_beetle_env.csv")
head(env)
#make spatial points layer

#create crs object
env.latlong<-data.frame(X=env$X, Y=env$Y)

class(beetle)

#Use coordinates object to create our spatial points object
beetle.sp<-vect(env.latlong, geom = c("X","Y"))

#check that the points now have our desired crs. 
crs(beetle.sp)<-"epsg:27700"
#check all plot together - but first project the env data to epsg:4326
beetle.sp<-project(beetle.sp,crs(UK))
plot(UK)
plot(beetle.sp, col="red", add=TRUE)
names(env)
plot(env)
