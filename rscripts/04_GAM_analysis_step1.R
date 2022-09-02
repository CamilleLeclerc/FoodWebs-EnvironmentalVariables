rm(list = ls()) #Removes all objects from the current workspace (R memory)
mypath <- rprojroot::find_package_root_file


##------------------------------------
##LOADING PACKAGES, FUNCTIONS AND DATA
##------------------------------------
##PACKAGES##
library(corrplot)
library(dplyr)
library(ggpubr)
library(mgcv)
library(PerformanceAnalytics)
library(randomForest)
library(spdep)
library(tibble)


##FUNCTIONS##
source("rfunctions/misc.R")
source("rfunctions/vif_func.R")


##DATA##
myload(envdata, list_lake, dir = mypath("data"))
myload(lake_pca_coord, dir = mypath("outputs"))




##----------------
# DATA PREPARATION
##----------------
geopoints <- as.matrix(cbind(list_lake$lat_pla, list_lake$long_pla))
knn <- spdep::knn2nb(spdep::knearneigh(geopoints, k = round(67/3), longlat = TRUE))
bw <- max(unlist(spdep::nbdists(knn, geopoints, longlat = TRUE)))


class(lake_pca_coord)
lake_pca_coord <- as.data.frame(lake_pca_coord)
class(lake_pca_coord)
ggdensity(lake_pca_coord$Dim.1, xlab = "Dim 1 coordinates")
ggdensity(lake_pca_coord$Dim.2, xlab = "Dim 2 coordinates")
ggqqplot(lake_pca_coord$Dim.1)
ggqqplot(lake_pca_coord$Dim.2)



rownames(envdata) <- envdata$cd.lac
envdata <- envdata%>% dplyr::select(-cd.lac)
chart.Correlation(envdata)
keep.dat <- vif_func(in_frame = envdata, thresh = 3, trace = T)




##----------------------------------
# GENERALIZED ADDITIVE MODELS (GAMs)
##----------------------------------
#Script based on https://github.com/JfvBraga/FoodwebSpace and https://doi.org/10.5281/zenodo.596810
colnames(envdata)
envdata$sqr_area <- (envdata$area)^2
envdata$sqr_max_depth <- (envdata$max_depth)^2
envdata$sqr_shorel_dvlpt_index <- (envdata$shorel_dvlpt_index)^2
envdata$sqr_Tmean <- (envdata$Tmean)^2
envdata$sqr_isothermality <- (envdata$isothermality)^2
envdata$sqr_Tseasonality <- (envdata$Tseasonality)^2
envdata$sqr_PT <- (envdata$PT)^2
envdata$sqr_NO3 <- (envdata$NO3)^2
envdata$sqr_COD <- (envdata$COD)^2


# 1.RESPONSE VARIABLE = TROPHIC DIVERSITY DESCRIPTOR (PCA AXIS 1)
#----------------------------------------------------------------
auto.dim1 <- autocov_dist(z = lake_pca_coord$Dim.1, xy = geopoints, longlat = TRUE, nbs = bw, style = "B")  #Create the autocavariate variable
colnames(envdata)
RF.dim1 <- randomForest(auto.dim1 ~ ., data = envdata, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.dim1.resid <- auto.dim1 -  RF.dim1$predicted     #Spatial residuals variable (obs minus fitted)
envdata$dim1 <- lake_pca_coord$Dim.1
envdata$sp_resid_dim1 <- RF.dim1.resid

gam.dim1 <- mgcv::gam(formula = dim1 ~  s(sp_resid_dim1, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = envdata,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  
summary(gam.dim1)
saveRDS(gam.dim1, "outputs/gam_dim1.rds")


gam.dim1.best <- mgcv::gam(formula = dim1 ~ 
                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                           data = envdata,
                           select = TRUE,
                           family = "gaussian",
                           method = "REML")    
summary(gam.dim1.best)
k.check(gam.dim1.best)
gam.check(gam.dim1.best) #https://noamross.github.io/gams-in-r-course/chapter2

type <- "deviance"  
resid <- residuals(gam.dim1.best, type = type)
linpred <- napredict(gam.dim1.best$na.action, gam.dim1.best$linear.predictors)
observed.y <- napredict(gam.dim1.best$na.action, gam.dim1.best$y)
par(mfrow = c(2,2))
qq.gam(gam.dim1.best, rep = 0, level = 0.9, type = type, rl.col = 2, rep.col = "gray80")
hist(resid, xlab = "Residuals", main = "Histogram of residuals")
plot(linpred, resid, main = "Resids vs. linear pred.", xlab = "linear predictor", ylab = "residuals")
plot(fitted(gam.dim1), observed.y, xlab = "Fitted Values", ylab = "Response", main = "Response vs. Fitted Values")

concurvity(gam.dim1.best, full = TRUE)

saveRDS(gam.dim1.best, "outputs/gam_dim1_reduce_model.rds")



    
# 2.RESPONSE VARIABLE = VERTICAL STRUCTURE DESCRIPTOR (PCA AXIS 2)
#-----------------------------------------------------------------
auto.dim2 <- autocov_dist(z = lake_pca_coord$Dim.2, xy = geopoints, longlat = TRUE, nbs = bw, style = "B")  #Create the autocavariate variable
colnames(envdata)
envdata <- select(envdata, -c(dim1, sp_resid_dim1))
RF.dim2 <- randomForest(auto.dim2 ~ ., data = envdata, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.dim2.resid <- auto.dim2 -  RF.dim2$predicted     #Spatial residuals variable (obs minus fitted)
envdata$dim2 <- lake_pca_coord$Dim.2
envdata$sp_resid_dim2 <- RF.dim2.resid


gam.dim2 <- mgcv::gam(formula = dim2 ~  s(sp_resid_dim2, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = envdata,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  
summary(gam.dim2)
saveRDS(gam.dim2, "outputs/gam_dim2.rds")


gam.dim2.best <- mgcv::gam(formula = dim2 ~ s(sp_resid_dim2, bs = "cr", k = 3) +
                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(Tmean, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3),
                           data = envdata,
                           select = TRUE,
                           family = "gaussian",
                           method = "REML")  
summary(gam.dim2.best)
k.check(gam.dim2.best)
gam.check(gam.dim2.best) #https://noamross.github.io/gams-in-r-course/chapter2

type <- "deviance"  
resid <- residuals(gam.dim2.best, type = type)
linpred <- napredict(gam.dim2.best$na.action, gam.dim2.best$linear.predictors)
observed.y <- napredict(gam.dim2.best$na.action, gam.dim2.best$y)
par(mfrow = c(2,2))
qq.gam(gam.dim2.best, rep = 0, level = 0.9, type = type, rl.col = 2, rep.col = "gray80")
hist(resid, xlab = "Residuals", main = "Histogram of residuals")
plot(linpred, resid, main = "Resids vs. linear pred.", xlab = "linear predictor", ylab = "residuals")
plot(fitted(gam.dim2), observed.y, xlab = "Fitted Values", ylab = "Response", main = "Response vs. Fitted Values")

concurvity(gam.dim2.best, full = TRUE)

saveRDS(gam.dim2.best, "outputs/gam_dim2_reduce_model.rds")
