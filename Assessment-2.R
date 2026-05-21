#Set working directory
setwd("../../Assessment-2/14349429_A2")

#Install packages
install.packages(c("terra","gdistance","raster","igraph","dplyr","tidyr"))
install.packages(c("sf", "Hmsc", "coda", "corrplot", "ggplot2"))
library(terra)
library(raster)
library(sf)
library(Hmsc)
library(coda)
library(corrplot)
library(dplyr)
library(tidyr)
library(ggplot2)

#Read in the species data
beetle<- read.csv("scot_beetle_community.csv")

#Load land cover for the analysis
UK=rast("LCMUK_2000.tif")

#Load environment data 
env<-read.csv("scot_beetle_env.csv")

#create CRS object
env.latlong<-data.frame(X=env$X, Y=env$Y)

#Use coordinates object to create our spatial points object
beetle.sp<-vect(env.latlong, geom = c("X","Y"))

#check that the points now have our desired CRS. 
crs(beetle.sp)<-"epsg:27700"

#check all plot together - but first project the species points on the UK landcover map 
beetle.sp<-project(beetle.sp,crs(UK))
plot(UK)
plot(beetle.sp, col="red", add=TRUE)

#Transforming beetles data to the same CRS as the environment data
beetleSF=st_as_sf(beetle.sp)

#Merge species and environmental data.
allData<-cbind(beetle, env)

#Remove duplicate columns
allData = allData[, !duplicated(names(allData))]

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

#Assigning site number to address sampling bias
Sites=as.factor(allData$Sites)
#Setting up as a dataframe as HMSC function only accepts dataframe
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
nIter <- 1000
nBurn <- 1000

set.seed(20)
#Fit the model
fit <- sampleMcmc(
  m,
  samples = nIter,
  transient = nBurn,
  nChains = nChains
)
dev.off()
#Check MCMC convergence
mcmc_chains <- convertToCodaObject(fit)

#adjust plotting scheme
par(mfrow = c(4, 4)) 
#Adjust margin sizes (bottom, left, top, right)
par(mar = c(2, 2, 2, 2))  
#Plot
plot(mcmc_chains$Beta)

#Extract the posterior estimates of the beta parameters (model coefficients)
postBeta <- getPostEstimate(fit, "Beta")

#Adjust margin sizes
par(mar = c(3, 3, 0.5, 0.5))  

#Plot the response of species to env covariates
plotBeta(fit, postBeta, param = "Mean")
#Extract the potential scale reduction factor as a model diagnostic
mpost = convertToCodaObject(fit)

psrf.omega = gelman.diag(mpost$Omega[[1]], multivariate =
                           FALSE)$psrf

#plot values
hist(psrf.omega)
#Compute associations using the fitted model
OmegaCor <- computeAssociations(fit)

#Extract the first set of associations from OmegaCor - mean between-species associations
toPlot <- OmegaCor[[1]]$mean
#Open a new plot window
dev.off()
#Plot the correlation matrix (inter specific interactions) using corrplot 
corrplot(toPlot,
         method = "color",
         col = colorRampPalette(c("blue", "white", "red"))(20),
         title = NULL,
         type = "lower",
         tl.col = "black",
         tl.cex = 0.6,  # Adjust text label size
         mar = c(0, 0, 0, 0))  # Adjust margins
#compute gradient
Gradient = constructGradient(fit,focalVariable = "pH")

#make prediction
predY = predict(fit, Gradient=Gradient, expected = TRUE)

#plot
plotGradient(fit, Gradient, pred=predY, measure="S", index = 1, showData = FALSE)

#Create gradient for spatial prediction
xy.grid = data.frame(x=env$X, y=env$Y)
#Envrironmental data grid
XData.grid=data.frame(apply(X, MARGIN=2, FUN=as.numeric))
nrow(X)==nrow(xy.grid)

#Compute AUC for model validation
preds=computePredictedValues(fit)
modelEval=evaluateModelFit(hM=fit, predY=preds)
cat("Mean AUC:", round(mean(modelEval$AUC, na.rm=TRUE), 3), "\n")

Gradient = prepareGradient(fit, XDataNew = XData.grid,sDataNew = xy.grid)

predY = predict(fit, Gradient=Gradient,expected=TRUE)

EpredY=Reduce("+",predY)/length(predY)

SpeciesRichness=rowSums(EpredY)

#coordinates
xy =xy.grid[,1:2]

#prediction for Beetle
H = EpredY[,1]

mapData=data.frame(xy,SpeciesRichness,H)

# We can plot variation in species richness.

ggplot(data = mapData, aes(x=x, y=y, color=SpeciesRichness))+geom_point(size=3)  + scale_color_gradient(low="blue",high="green") + coord_equal()

