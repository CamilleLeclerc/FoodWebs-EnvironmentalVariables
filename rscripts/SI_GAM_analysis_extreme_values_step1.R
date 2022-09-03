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
rownames(geopoints) <- list_lake$cd.lac


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


#1.WITHOUT EXTREM VALUES OF SHORELINE DEVELOPMENT INDEX
#------------------------------------------------------
#RESPONSE VARIABLE = TROPHIC DIVERSITY DESCRIPTOR (PCA AXIS 1)
SDI <- envdata %>% dplyr::filter(shorel_dvlpt_index < max(shorel_dvlpt_index))
SDIgeopoints <- geopoints %>% as.data.frame(.) %>% dplyr::filter(rownames(.) %in% rownames(SDI))
SDIknn <- spdep::knn2nb(spdep::knearneigh(as.matrix(SDIgeopoints), k = round(66/3), longlat = TRUE))
SDIbw <- max(unlist(spdep::nbdists(SDIknn, as.matrix(SDIgeopoints), longlat = TRUE)))

auto.dim1 <- autocov_dist(z = lake_pca_coord %>% dplyr::filter(rownames(.) %in% rownames(SDI)) %>% dplyr::select(Dim.1) %>% as.list(.) %>% unlist(.),
                          xy = as.matrix(SDIgeopoints),
                          longlat = TRUE,
                          nbs = SDIbw, 
                          style = "B")  #Create the autocavariate variable

colnames(SDI)
RF.dim1 <- randomForest(auto.dim1 ~ ., data = SDI, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.dim1.resid <- auto.dim1 -  RF.dim1$predicted     #Spatial residuals variable (obs minus fitted)
lakeSDI <- lake_pca_coord %>% dplyr::select(Dim.1) %>% dplyr::filter(rownames(.) %in% rownames(SDI))
SDI$dim1 <- lakeSDI$Dim.1
SDI$sp_resid_dim1 <- RF.dim1.resid

gam.dim1 <- mgcv::gam(formula = dim1 ~  s(sp_resid_dim1, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = SDI,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  
summary(gam.dim1)
saveRDS(gam.dim1, "outputs/gam_dim1_extremevalueSDI.rds")


gam.dim1.best <- mgcv::gam(formula = dim1 ~ 
                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                           data = SDI,
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

saveRDS(gam.dim1.best, "outputs/gam_dim1_extremevalueSDI_reduce_model.rds")
rm(auto.dim1, RF.dim1, RF.dim1.resid, lakeSDI, gam.dim1, gam.dim1.best, type, resid, linpred, observed.y)


#RESPONSE VARIABLE = VERTICAL STRUCTURE DESCRIPTOR (PCA AXIS 2)
auto.dim2 <- autocov_dist(z = lake_pca_coord %>% dplyr::filter(rownames(.) %in% rownames(SDI)) %>% dplyr::select(Dim.2) %>% as.list(.) %>% unlist(.),
                          xy = as.matrix(SDIgeopoints),
                          longlat = TRUE,
                          nbs = SDIbw, 
                          style = "B")  #Create the autocavariate variable

colnames(SDI)
RF.dim2 <- randomForest(auto.dim2 ~ ., data = SDI, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.dim2.resid <- auto.dim2 -  RF.dim2$predicted     #Spatial residuals variable (obs minus fitted)
lakeSDI <- lake_pca_coord %>% dplyr::select(Dim.2) %>% dplyr::filter(rownames(.) %in% rownames(SDI))
SDI$dim2 <- lakeSDI$Dim.2
SDI$sp_resid_dim2 <- RF.dim2.resid

gam.dim2 <- mgcv::gam(formula = dim2 ~  s(sp_resid_dim2, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = SDI,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  
summary(gam.dim2)
saveRDS(gam.dim2, "outputs/gam_dim2_extremevalueSDI.rds")


gam.dim2.best <- mgcv::gam(formula = dim2 ~ s(sp_resid_dim2, bs = "cr", k = 3) +
                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(Tmean, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3),
                           data = SDI,
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

saveRDS(gam.dim2.best, "outputs/gam_dim2_extremevalueSDI_reduce_model.rds")
rm(auto.dim2, RF.dim2, RF.dim2.resid, lakeSDI, gam.dim2, gam.dim2.best, type, resid, linpred, observed.y)
rm(SDI, SDIgeopoints, SDIknn, SDIbw)




#2.WITHOUT EXTREM VALUES OF DISSOLVED ORGANIC CARBON
#---------------------------------------------------
#RESPONSE VARIABLE = TROPHIC DIVERSITY DESCRIPTOR (PCA AXIS 1)
DOC <- envdata %>% dplyr::filter(COD < max(COD))
DOCgeopoints <- geopoints %>% as.data.frame(.) %>% dplyr::filter(rownames(.) %in% rownames(DOC))
DOCknn <- spdep::knn2nb(spdep::knearneigh(as.matrix(DOCgeopoints), k = round(66/3), longlat = TRUE))
DOCbw <- max(unlist(spdep::nbdists(DOCknn, as.matrix(DOCgeopoints), longlat = TRUE)))

auto.dim1 <- autocov_dist(z = lake_pca_coord %>% dplyr::filter(rownames(.) %in% rownames(DOC)) %>% dplyr::select(Dim.1) %>% as.list(.) %>% unlist(.),
                          xy = as.matrix(DOCgeopoints),
                          longlat = TRUE,
                          nbs = DOCbw, 
                          style = "B")  #Create the autocavariate variable

colnames(DOC)
RF.dim1 <- randomForest(auto.dim1 ~ ., data = DOC, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.dim1.resid <- auto.dim1 -  RF.dim1$predicted     #Spatial residuals variable (obs minus fitted)
lakeDOC <- lake_pca_coord %>% dplyr::select(Dim.1) %>% dplyr::filter(rownames(.) %in% rownames(DOC))
DOC$dim1 <- lakeDOC$Dim.1
DOC$sp_resid_dim1 <- RF.dim1.resid

gam.dim1 <- mgcv::gam(formula = dim1 ~  s(sp_resid_dim1, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = DOC,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  
summary(gam.dim1)
saveRDS(gam.dim1, "outputs/gam_dim1_extremevalueDOC.rds")


gam.dim1.best <- mgcv::gam(formula = dim1 ~ 
                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                           data = DOC,
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

saveRDS(gam.dim1.best, "outputs/gam_dim1_extremevalueDOC_reduce_model.rds")
rm(auto.dim1, RF.dim1, RF.dim1.resid, lakeDOC, gam.dim1, gam.dim1.best, type, resid, linpred, observed.y)


#RESPONSE VARIABLE = VERTICAL STRUCTURE DESCRIPTOR (PCA AXIS 2)
auto.dim2 <- autocov_dist(z = lake_pca_coord %>% dplyr::filter(rownames(.) %in% rownames(DOC)) %>% dplyr::select(Dim.2) %>% as.list(.) %>% unlist(.),
                          xy = as.matrix(DOCgeopoints),
                          longlat = TRUE,
                          nbs = DOCbw, 
                          style = "B")  #Create the autocavariate variable

colnames(DOC)
RF.dim2 <- randomForest(auto.dim2 ~ ., data = DOC, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.dim2.resid <- auto.dim2 -  RF.dim2$predicted     #Spatial residuals variable (obs minus fitted)
lakeDOC <- lake_pca_coord %>% dplyr::select(Dim.2) %>% dplyr::filter(rownames(.) %in% rownames(DOC))
DOC$dim2 <- lakeDOC$Dim.2
DOC$sp_resid_dim2 <- RF.dim2.resid

gam.dim2 <- mgcv::gam(formula = dim2 ~  s(sp_resid_dim2, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = DOC,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  
summary(gam.dim2)
saveRDS(gam.dim2, "outputs/gam_dim2_extremevalueDOC.rds")


gam.dim2.best <- mgcv::gam(formula = dim2 ~ s(sp_resid_dim2, bs = "cr", k = 3) +
                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(Tmean, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3),
                           data = DOC,
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

saveRDS(gam.dim2.best, "outputs/gam_dim2_extremevalueDOC_reduce_model.rds")
rm(auto.dim2, RF.dim2, RF.dim2.resid, lakeDOC, gam.dim2, gam.dim2.best, type, resid, linpred, observed.y)
rm(DOC, DOCgeopoints, DOCknn, DOCbw)




#3.WITHOUT EXTREM VALUES OF MEAN ANNUAL TEMPERATURE
#--------------------------------------------------
#RESPONSE VARIABLE = TROPHIC DIVERSITY DESCRIPTOR (PCA AXIS 1)
TempMean <- envdata %>% dplyr::filter(Tmean < max(Tmean))
TempMeangeopoints <- geopoints %>% as.data.frame(.) %>% dplyr::filter(rownames(.) %in% rownames(TempMean))
TempMeanknn <- spdep::knn2nb(spdep::knearneigh(as.matrix(TempMeangeopoints), k = round(66/3), longlat = TRUE))
TempMeanbw <- max(unlist(spdep::nbdists(TempMeanknn, as.matrix(TempMeangeopoints), longlat = TRUE)))

auto.dim1 <- autocov_dist(z = lake_pca_coord %>% dplyr::filter(rownames(.) %in% rownames(TempMean)) %>% dplyr::select(Dim.1) %>% as.list(.) %>% unlist(.),
                          xy = as.matrix(TempMeangeopoints),
                          longlat = TRUE,
                          nbs = TempMeanbw, 
                          style = "B")  #Create the autocavariate variable

colnames(TempMean)
RF.dim1 <- randomForest(auto.dim1 ~ ., data = TempMean, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.dim1.resid <- auto.dim1 -  RF.dim1$predicted     #Spatial residuals variable (obs minus fitted)
lakeTempMean <- lake_pca_coord %>% dplyr::select(Dim.1) %>% dplyr::filter(rownames(.) %in% rownames(TempMean))
TempMean$dim1 <- lakeTempMean$Dim.1
TempMean$sp_resid_dim1 <- RF.dim1.resid

gam.dim1 <- mgcv::gam(formula = dim1 ~  s(sp_resid_dim1, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = TempMean,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  
summary(gam.dim1)
saveRDS(gam.dim1, "outputs/gam_dim1_extremevalueTempMean.rds")


gam.dim1.best <- mgcv::gam(formula = dim1 ~ 
                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                           data = TempMean,
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

saveRDS(gam.dim1.best, "outputs/gam_dim1_extremevalueTempMean_reduce_model.rds")
rm(auto.dim1, RF.dim1, RF.dim1.resid, lakeTempMean, gam.dim1, gam.dim1.best, type, resid, linpred, observed.y)


#RESPONSE VARIABLE = VERTICAL STRUCTURE DESCRIPTOR (PCA AXIS 2)
auto.dim2 <- autocov_dist(z = lake_pca_coord %>% dplyr::filter(rownames(.) %in% rownames(TempMean)) %>% dplyr::select(Dim.2) %>% as.list(.) %>% unlist(.),
                          xy = as.matrix(TempMeangeopoints),
                          longlat = TRUE,
                          nbs = TempMeanbw, 
                          style = "B")  #Create the autocavariate variable

colnames(TempMean)
RF.dim2 <- randomForest(auto.dim2 ~ ., data = TempMean, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.dim2.resid <- auto.dim2 -  RF.dim2$predicted     #Spatial residuals variable (obs minus fitted)
lakeTempMean <- lake_pca_coord %>% dplyr::select(Dim.2) %>% dplyr::filter(rownames(.) %in% rownames(TempMean))
TempMean$dim2 <- lakeTempMean$Dim.2
TempMean$sp_resid_dim2 <- RF.dim2.resid

gam.dim2 <- mgcv::gam(formula = dim2 ~  s(sp_resid_dim2, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = TempMean,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  
summary(gam.dim2)
saveRDS(gam.dim2, "outputs/gam_dim2_extremevalueTempMean.rds")


gam.dim2.best <- mgcv::gam(formula = dim2 ~ s(sp_resid_dim2, bs = "cr", k = 3) +
                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(Tmean, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3),
                           data = TempMean,
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

saveRDS(gam.dim2.best, "outputs/gam_dim2_extremevalueTempMean_reduce_model.rds")
rm(auto.dim2, RF.dim2, RF.dim2.resid, lakeTempMean, gam.dim2, gam.dim2.best, type, resid, linpred, observed.y)
rm(TempMean, TempMeangeopoints, TempMeanknn, TempMeanbw)




#4.WITHOUT EXTREM VALUES OF TEMPERATURE SEASONALITY
#--------------------------------------------------
#RESPONSE VARIABLE = TROPHIC DIVERSITY DESCRIPTOR (PCA AXIS 1)
TempSeason <- envdata %>% dplyr::filter(Tseasonality < max(Tseasonality))
TempSeasongeopoints <- geopoints %>% as.data.frame(.) %>% dplyr::filter(rownames(.) %in% rownames(TempSeason))
TempSeasonknn <- spdep::knn2nb(spdep::knearneigh(as.matrix(TempSeasongeopoints), k = round(66/3), longlat = TRUE))
TempSeasonbw <- max(unlist(spdep::nbdists(TempSeasonknn, as.matrix(TempSeasongeopoints), longlat = TRUE)))

auto.dim1 <- autocov_dist(z = lake_pca_coord %>% dplyr::filter(rownames(.) %in% rownames(TempSeason)) %>% dplyr::select(Dim.1) %>% as.list(.) %>% unlist(.),
                          xy = as.matrix(TempSeasongeopoints),
                          longlat = TRUE,
                          nbs = TempSeasonbw, 
                          style = "B")  #Create the autocavariate variable

colnames(TempSeason)
RF.dim1 <- randomForest(auto.dim1 ~ ., data = TempSeason, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.dim1.resid <- auto.dim1 -  RF.dim1$predicted     #Spatial residuals variable (obs minus fitted)
lakeTempSeason <- lake_pca_coord %>% dplyr::select(Dim.1) %>% dplyr::filter(rownames(.) %in% rownames(TempSeason))
TempSeason$dim1 <- lakeTempSeason$Dim.1
TempSeason$sp_resid_dim1 <- RF.dim1.resid

gam.dim1 <- mgcv::gam(formula = dim1 ~  s(sp_resid_dim1, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = TempSeason,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  
summary(gam.dim1)
saveRDS(gam.dim1, "outputs/gam_dim1_extremevalueTempSeason.rds")


gam.dim1.best <- mgcv::gam(formula = dim1 ~ 
                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(isothermality, bs = "cr", k = 3) + 
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                           data = TempSeason,
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

saveRDS(gam.dim1.best, "outputs/gam_dim1_extremevalueTempSeason_reduce_model.rds")
rm(auto.dim1, RF.dim1, RF.dim1.resid, lakeTempSeason, gam.dim1, gam.dim1.best, type, resid, linpred, observed.y)


#RESPONSE VARIABLE = VERTICAL STRUCTURE DESCRIPTOR (PCA AXIS 2)
auto.dim2 <- autocov_dist(z = lake_pca_coord %>% dplyr::filter(rownames(.) %in% rownames(TempSeason)) %>% dplyr::select(Dim.2) %>% as.list(.) %>% unlist(.),
                          xy = as.matrix(TempSeasongeopoints),
                          longlat = TRUE,
                          nbs = TempSeasonbw, 
                          style = "B")  #Create the autocavariate variable

colnames(TempSeason)
RF.dim2 <- randomForest(auto.dim2 ~ ., data = TempSeason, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.dim2.resid <- auto.dim2 -  RF.dim2$predicted     #Spatial residuals variable (obs minus fitted)
lakeTempSeason <- lake_pca_coord %>% dplyr::select(Dim.2) %>% dplyr::filter(rownames(.) %in% rownames(TempSeason))
TempSeason$dim2 <- lakeTempSeason$Dim.2
TempSeason$sp_resid_dim2 <- RF.dim2.resid

gam.dim2 <- mgcv::gam(formula = dim2 ~  s(sp_resid_dim2, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = TempSeason,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  
summary(gam.dim2)
saveRDS(gam.dim2, "outputs/gam_dim2_extremevalueTempSeason.rds")


gam.dim2.best <- mgcv::gam(formula = dim2 ~ s(sp_resid_dim2, bs = "cr", k = 3) +
                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(Tmean, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3),
                           data = TempSeason,
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

saveRDS(gam.dim2.best, "outputs/gam_dim2_extremevalueTempSeason_reduce_model.rds")
rm(auto.dim2, RF.dim2, RF.dim2.resid, lakeTempSeason, gam.dim2, gam.dim2.best, type, resid, linpred, observed.y)
rm(TempSeason, TempSeasongeopoints, TempSeasonknn, TempSeasonbw)


