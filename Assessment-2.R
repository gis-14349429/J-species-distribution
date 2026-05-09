#Set working directory
setwd("C:/Users/User/OneDrive - The University of Manchester/GEOG71922/Assessment-2/14349429_A2")

#Install packages
install.packages(c("terra","gdistance","raster","igraph"))
library(terra)
library(raster)
#Read in the species data
beetle<- read.csv("scot_beetle_community.csv")

#view the first six rows of the data
head(beetle) 

#Load landcover for the analysis
UK=rast("LCMUK_2000.tif")
#Plot the landcover data
plot(UK)
env<-read.csv("scot_beetle_env.csv")
head(env)

