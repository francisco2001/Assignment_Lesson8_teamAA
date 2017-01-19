#Assignment Lesson 8
#Team AA
#Authors : Arnan Araza ; Francisco Arias 

#Cleaning Workspace
rm(list=ls())

#Importing Libraries
library(raster)

#Calling Functions
source ("R/re_extreme.R")
source("R/Subtract.R")

## Load data
load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")
load("data/trainingPoly.rda")


#Assigning a value of NA to these pixels with values above 100 representing others than forest

vcfGewata[vcfGewata > 100] <- NA

# Build a brick containing all data

gewata <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7,vcfGewata)
names(gewata) <- c("band1", "band2", "band3", "band4", "band5", "band7","VCF")

hist(gewata)

# Extract all data to a data.frame
df <- as.data.frame(getValues(gewata))

### Relationship between the Landsat bands and the VCF. 

pairs (gewata)

### Applying a Linear Model

model <- lm(VCF ~ band1 + band2 + band3 + band4 + band5 + band7, data = df)
summary(model)
squareroot <- summary(model)$r.square

### Prediction Model

predict_VCF <- predict(gewata, model=model, na.rm=TRUE)

### Plotting the predicted tree cover and comparing it with the original forest cover

opar <- par(mfrow=c(1,2))
plot(predict_VCF, main="Predicted forest cover", zlim=c(0,100))
plot(vcfGewata, main="Original forest cover")
par(opar)

# Subtraction between predicted tree cover and original VCF  raster values.
### Comparisson 

SubtVFC_predict <- vcfGewata - predict_VCF
par(mfrow = c(1,2))
plot(SubtVFC_predict, main = "Diference between predicted and original VCF")
hist(SubtVFC_predict)


### Compute the RMSE between predicted and the actual tree cover values

#First we convert rasterLayers to dataframe

predicted_df <- as.data.frame(predict_VCF)
vcfGewata_df <- as.data.frame(vcfGewata)

# Second Calculating the RMSE

RMSE <- sqrt( mean((predicted_df-vcfGewata_df)^2 , na.rm = TRUE))
print ( paste("The RMSE value between the predicted and the actual tree cover is", RMSE))

#### Calculating RMSE for trainingpolygons

#Converting using the as.numeric() function, 
# which takes the factor levels

trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)

#Rasterization of trainingPolygons

classes <- rasterize(trainingPoly, gewata, field='Code')

# Applying zonal statistics using mean value 

VCF_orig_zonal <- zonal(gewata$VCF, classes, fun=mean, na.rm=TRUE)
predict_zonal <- zonal(predict_VCF, classes, fun=mean, na.rm=TRUE)

#Tranforming to data.frame. to be use in RMSE calculations 

VCF_orig_zonal_df <- as.data.frame(VCF_orig_zonal)
predict_zonal_df <- as.data.frame(predict_zonal)

#RMSE using training polygons

RMSE_class1 <- sqrt(mean(VCF_orig_zonal_df[1,2]-predict_zonal_df[1,2])^2)
RMSE_class2 <- sqrt(mean(VCF_orig_zonal_df[2,2]-predict_zonal_df[2,2])^2)
RMSE_class3 <- sqrt(mean(VCF_orig_zonal_df[3,2]-predict_zonal_df[3,2])^2)

print ( paste("The RMSE value between the predicted and the actual tree cover is", RMSE_class1))
print ( paste("The RMSE value between the predicted and the actual tree cover is", RMSE_class2))
print ( paste("The RMSE value between the predicted and the actual tree cover is", RMSE_class3))

# Creating a data.frame of RMSE per zone.  

Class_RMSE_matrix = matrix(c(1, 2, 3, round(RMSE_class1, digits=2), round(RMSE_class2, digits=2), round(RMSE_class3, digits=2),"cropland","forest", "wetland"), nrow=3, ncol=3)
colnames(Class_RMSE_matrix) <- c("Class", "RMSE", "Description")
print (Class_RMSE_matrix)











