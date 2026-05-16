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
env<-read.csv("scot_beetle_env.csv")
head(env)
class(beetle)
names(env)
#create crs object
env.latlong<-data.frame(X=env$X, Y=env$Y)
# Install necessary packages if they are not already installed
install.packages(c("Hmsc", "coda", "corrplot"))
library(Hmsc)
library(coda)
library(corrplot)
#Use coordinates object to create our spatial points object
beetle.sp<-vect(env.latlong, geom = c("X","Y"))
#check that the points now have our desired crs. 
crs(beetle.sp)<-"epsg:27700"
#check all plot together - but first project the env data to epsg:4326
beetle.sp<-project(beetle.sp,crs(UK))
plot(UK)
plot(beetle.sp, col="red", add=TRUE)

#Transforming beetles data to the same CRS as the environment data
beetleSF=st_as_sf(beetle.sp,coords=c("x","y"),crs="EPSG:27700")
#Extracting the environment data to the points data

