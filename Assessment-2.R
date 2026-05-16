#Set working directory
setwd("C:/Users/User/OneDrive - The University of Manchester/GEOG71922/Assessment-2/14349429_A2")

#Install packages
install.packages(c("terra","gdistance","raster","igraph"))
install.packages(c("sf", "Hmsc", "coda", "corrplot"))
library(terra)
library(raster)
library(sf)
library(Hmsc)
library(coda)
library(corrplot)
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
#Use coordinates object to create our spatial points object
beetle.sp<-vect(env.latlong, geom = c("X","Y"))
#check that the points now have our desired crs. 
crs(beetle.sp)<-"epsg:27700"
#check all plot together - but first project the env data to epsg:4326
beetle.sp<-project(beetle.sp,crs(UK))
plot(UK)
plot(beetle.sp, col="red", add=TRUE)

#Transforming beetles data to the same CRS as the environment data
beetleSF=st_as_sf(beetle.sp)
#Extracting the environment data to the points data
envData=terra::extract(UK,beetleSF)
crdsbeetleSF=st_coordinates(beetleSF)
#Get cell number at different locations
cn=cellFromXY(UK, crdsbeetleSF)
library(dplyr)
library(tidyr)
#merge data by the "cn" column.
allData<-cbind(beetle, env)
head(allData)
names(allData)
#keep only unqiue coordinates
allData=distinct(bradypus,geometry, .keep_all = TRUE)
