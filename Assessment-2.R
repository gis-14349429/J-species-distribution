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
class(env)
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
XFormula = ~ pH + Moist + Plants_m2 + Biom0_5 + Org + Elevation + Management
#Construct the model
m <- Hmsc(Y = Y, XData = X, XFormula = XFormula, distr = "probit",
          studyDesign = studyDesign,ranLevels = list(sample = rL))
# Set MCMC parameters
nChains <- 3
nIter <- 2000
nBurn <- 2000

set.seed(11)
# Fit the model
fit <- sampleMcmc(
  m,
  samples = nIter,
  transient = nBurn,
  nChains = nChains
)
dev.off()
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
#Compute associations using the fitted model
OmegaCor <- computeAssociations(fit)

# Extract the first set of associations from OmegaCor - mean between-species associations
toPlot <- OmegaCor[[1]]$mean
# Open a new plot window
dev.off()
# Plot the correlation matrix (interspecific interactions) using corrplot 
corrplot(toPlot,
         method = "color",
         col = colorRampPalette(c("blue", "white", "red"))(20),
         title = NULL,
         type = "lower",
         tl.col = "black",
         tl.cex = 0.5,  # Adjust text label size
         mar = c(0, 0, 6, 0))  # Adjust margins
#compute gradient
Gradient = constructGradient(fit,focalVariable = "pH")

#make prediction
predY = predict(fit, Gradient=Gradient, expected = TRUE)

#plot
plotGradient(fit, Gradient, pred=predY, measure="S", index = 1, showData = FALSE)
# environmental predictors from the env layer 

#Set NAs to zero 
env[is.na(env)]<-0

#get coordinates of all raster cells
env_grid = read.csv("scot_beetle_env.csv")

#transform all raster covariates to a matrix of values
xy.grid = data.frame(x=env_grid$X, y=env_grid$Y)

#Extratcing variables only in the formula
XData.grid<- data.frame(
  Org = as.numeric(env_grid$Org),
  pH = as.numeric(env_grid$pH),
  Moist = as.numeric(env_grid$Moist),
  Plants_m2 = as.numeric(env_grid$Plants_m2),
  Biom0_5 = as.numeric(env_grid$Biom0_5),
  Elevation = as.numeric(env_grid$Elevation),
  Management = as.numeric(env_grid$Management)
)
nrow(XData.grid)==nrow(xy.grid)

# We next use the prepareGradient function to convert the environmental and spatial
# predictors into a format that can be used as input for the predict function

Gradient = prepareGradient(fit, XDataNew = XData.grid,sDataNew = xy.grid)



# We are now ready to compute the predicted distribution (takes a minute so please be patient)

# Note that we used expected = TRUE to predict occurrence probabilities (e.g. 0.2) instead of occurrences (0 or 1) 
predY = predict(fit, Gradient=Gradient,expected=TRUE)


EpredY=Reduce("+",predY)/length(predY)


# EpredY is a matrix of posterior mean occurrence probabilities
# The next step is to post-process the predictions to those community features
# that we wish to illustrate over the prediction space. With the script below,
# we derive from the predictions for species richness


SpeciesRichness=rowSums(EpredY)

#coordinates
xy =xy.grid[,1:2]

#prediction for Beetle
H = EpredY[,1]

mapData=data.frame(xy,SpeciesRichness,H)


# We will use the ggplot function from the ggplot2 package, 
library(ggplot2)

# We can plot variation in species richness.

ggplot(data = mapData, aes(x=x, y=y, color=SpeciesRichness))+geom_point(size=3)  + scale_color_gradient(low="blue",high="green") + coord_equal()
colnames(XData.grid)
