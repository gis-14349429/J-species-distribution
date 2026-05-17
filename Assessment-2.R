#Set working directory
setwd("C:/Users/s16188ss/OneDrive - The University of Manchester/GEOG71922/Assessment-2/14349429_A2")

#Install packages
install.packages(c("terra","gdistance","raster","igraph"))
install.packages(c("sf", "Hmsc", "coda", "corrplot"))
library(terra)
library(raster)
library(sf)
library(Hmsc)
library(coda)
library(corrplot)
library(dplyr)
library(tidyr)

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
#merge species and environmental data.
allData<-cbind(beetle, env)
head(allData)
names(allData)

#Remove duplicate columns
allData = allData[, !duplicated(names(allData))]
ncol(allData)
names(allData)

# Prepare species data matrix - only need columns 3 to 24
Y <- data.frame(allData[, c("sp1", "sp2", "sp3", "sp4", "sp5", "sp6",
                                      "sp7", "sp8", "sp9", "sp10", "sp11", "sp12", 
                                      "sp13", "sp14", "sp15", "sp16", "sp17", "sp18", 
                                      "sp19", "sp20", "sp21", "sp22")])
X <- data.frame(allData[, c("Org", "pH", "AvailP", "AvailK", "Moist", "Litter", 
                                      "Bryophyte", "Plants_m2", "CanopyHeight", "Stem.density", 
                                      "Biom0_5", "Repro.biom", "Elevation", "Management")])

#convert to data frames and ensure all columns are numeric
Y=data.frame(apply(Y,MARGIN=2,FUN=as.numeric))

#convert env data
X=data.frame(apply(X,MARGIN=2,FUN=as.numeric))
nrow(Y)==nrow(X)
nrow(Y)
nrow(X)
ncol(X)
names(X)
#Assigning site number to address sampling bias
Sites=as.factor(allData$Sites)
#this needs to be set up as a data frame - as the Hmsc( ) function only accepts data frames (not matrices)
studyDesign <- data.frame(sample = Sites)
#use the HmscRandomLevel() function to store the random level
rL <- HmscRandomLevel(units = studyDesign$sample)
#Write model formula
XFormula = ~ pH + Moist + CanopyHeight + Bryophyte + AvailK + AvailP + Stem.density
#Construct the model
m <- Hmsc(Y = Y, XData = X, XFormula = XFormula, distr = "probit",
          studyDesign = studyDesign,ranLevels = list(sample = rL))
# Set MCMC parameters
nChains <- 2
nIter <- 1000
nBurn <- 1000

set.seed(11)
# Fit the model
fit <- sampleMcmc(
  m,
  samples = nIter,
  transient = nBurn,
  nChains = nChains
)
# Check MCMC convergence
mcmc_chains <- convertToCodaObject(fit)

par(mfrow = c(4, 4)) #adjust plotting scheme
par(mar = c(2, 2, 2, 2))  # Adjust margin sizes (bottom, left, top, right)
plot(mcmc_chains$Beta)#plot

#Extract the posterior estimates of the beta parameters (model coefficients)
postBeta <- getPostEstimate(fit, "Beta")

# Adjust margin sizes
par(mar = c(3, 3, 0.5, 0.5))  

#Plot the response of species to env covariates
plotBeta(fit, postBeta, param = "Mean")
#Extract the potential scale reduction factor as a model diagnostic (should ideally be close to one for most runs)
mpost = convertToCodaObject(fit)

psrf.omega = gelman.diag(mpost$Omega[[1]], multivariate =
                           FALSE)$psrf

#plot values
hist(psrf.omega)
