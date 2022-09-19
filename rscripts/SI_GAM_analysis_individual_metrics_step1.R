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
myload(lake_pca_coord, topological_metrics, dir = mypath("outputs"))




##----------------
# DATA PREPARATION
##----------------
geopoints <- as.matrix(cbind(list_lake$lat_pla, list_lake$long_pla))
knn <- spdep::knn2nb(spdep::knearneigh(geopoints, k = round(67/3), longlat = TRUE))
bw <- max(unlist(spdep::nbdists(knn, geopoints, longlat = TRUE)))


ggdensity(topological_metrics$SR.total, xlab = "SR.total")
ggdensity(topological_metrics$SR.invertebrate, xlab = "SR.invertebrate")
ggdensity(topological_metrics$Linkage.density, xlab = "Linkage.density")
ggdensity(topological_metrics$SR.vertebrate, xlab = "SR.vertebrate")
ggdensity(topological_metrics$Mean.food.chain.length, xlab = "Mean.food.chain.length")
ggdensity(topological_metrics$Mean.trophic.level, xlab = "Mean.trophic.level")
ggqqplot(topological_metrics$SR.total)
ggqqplot(topological_metrics$SR.invertebrate)
ggqqplot(topological_metrics$Linkage.density)
ggqqplot(topological_metrics$SR.vertebrate)
ggqqplot(topological_metrics$Mean.food.chain.length)
ggqqplot(topological_metrics$Mean.trophic.level)


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


# 1.RESPONSE VARIABLE = SR.TOTAL
#-------------------------------
auto.srtotal <-  autocov_dist(z = topological_metrics$SR.total, xy = geopoints, longlat = TRUE, nbs = bw, style = "B")     # create the autocavariate variable
colnames(envdata)
RF.srtotal <- randomForest(auto.srtotal ~ ., data = envdata, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.srtotal.resid <- auto.srtotal -  RF.srtotal$predicted     #Spatial residuals variable (obs minus fitted)
envdata$srtotal <- topological_metrics$SR.total
envdata$sp_resid_srtotal <- RF.srtotal.resid

gam.srtotal <- mgcv::gam(formula = srtotal ~  s(sp_resid_srtotal, bs = "cr", k = 3) +
                                              s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                              s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                              s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = envdata,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  

summary(gam.srtotal)
saveRDS(gam.srtotal, "outputs/SI_IndividualMetricsAnalysis/gam_srtotal.rds")


gam.srtotal.best <- mgcv::gam(formula = srtotal ~ 
                                                  s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                                  s(isothermality, bs = "cr", k = 3) + 
                                                  s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3),
                              data = envdata,
                              select = TRUE,
                              family = "gaussian",
                              method = "REML")  
summary(gam.srtotal.best)

saveRDS(gam.srtotal.best, "outputs/SI_IndividualMetricsAnalysis/gam_srtotal_reduce_model.rds")


# 2.RESPONSE VARIABLE = SR.INVERTEBRATE
#--------------------------------------
auto.srinvertebrate <-  autocov_dist(z = topological_metrics$SR.invertebrate, xy = geopoints, longlat = TRUE, nbs = bw, style = "B")     # create the autocavariate variable
colnames(envdata)
envdata <- select(envdata, -c(srtotal, sp_resid_srtotal))
RF.srinvertebrate <- randomForest(auto.srinvertebrate ~ ., data = envdata, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.srinvertebrate.resid <- auto.srinvertebrate -  RF.srinvertebrate$predicted     #Spatial residuals variable (obs minus fitted)
envdata$srinvertebrate <- topological_metrics$SR.invertebrate
envdata$sp_resid_srinvertebrate <- RF.srinvertebrate.resid

gam.srinvertebrate <- mgcv::gam(formula = srinvertebrate ~  s(sp_resid_srinvertebrate, bs = "cr", k = 3) +
                                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                                            s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                         data = envdata,
                         select = TRUE,
                         family = "gaussian",
                         method = "REML")  

summary(gam.srinvertebrate)
saveRDS(gam.srinvertebrate, "outputs/SI_IndividualMetricsAnalysis/gam_srinvertebrate.rds")


gam.srinvertebrate.best <- mgcv::gam(formula = srinvertebrate ~ 
                                                                s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                                                s(isothermality, bs = "cr", k = 3) + 
                                                                s(PT, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                              data = envdata,
                              select = TRUE,
                              family = "gaussian",
                              method = "REML")  
summary(gam.srinvertebrate.best)

saveRDS(gam.srinvertebrate.best, "outputs/SI_IndividualMetricsAnalysis/gam_srinvertebrate_reduce_model.rds")


# 3.RESPONSE VARIABLE = LINKAGE.DENSITY
#--------------------------------------
auto.LD <-  autocov_dist(z = topological_metrics$Linkage.density, xy = geopoints, longlat = TRUE, nbs = bw, style = "B")     # create the autocavariate variable
colnames(envdata)
envdata <- select(envdata, -c(srinvertebrate, sp_resid_srinvertebrate))
RF.LD <- randomForest(auto.LD ~ ., data = envdata, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.LD.resid <- auto.LD -  RF.LD$predicted     #Spatial residuals variable (obs minus fitted)
envdata$LD <- topological_metrics$Linkage.density
envdata$sp_resid_LD <- RF.LD.resid

gam.LD <- mgcv::gam(formula = LD ~  s(sp_resid_LD, bs = "cr", k = 3) +
                                    s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                    s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                    s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                                data = envdata,
                                select = TRUE,
                                family = "gaussian",
                                method = "REML")  

summary(gam.LD)
saveRDS(gam.LD, "outputs/SI_IndividualMetricsAnalysis/gam_LD.rds")


gam.LD.best <- mgcv::gam(formula = LD ~ s(sp_resid_LD, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                                     data = envdata,
                                     select = TRUE,
                                     family = "gaussian",
                                     method = "REML")  
summary(gam.LD.best)

saveRDS(gam.LD.best, "outputs/SI_IndividualMetricsAnalysis/gam_LD_reduce_model.rds")


# 4.RESPONSE VARIABLE = SR.VERTEBRATE
#------------------------------------
auto.srvertebrate <-  autocov_dist(z = topological_metrics$SR.vertebrate, xy = geopoints, longlat = TRUE, nbs = bw, style = "B")     # create the autocavariate variable
colnames(envdata)
envdata <- select(envdata, -c(LD, sp_resid_LD))
RF.srvertebrate <- randomForest(auto.srvertebrate ~ ., data = envdata, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.srvertebrate.resid <- auto.srvertebrate -  RF.srvertebrate$predicted     #Spatial residuals variable (obs minus fitted)
envdata$srvertebrate <- topological_metrics$SR.vertebrate
envdata$sp_resid_srvertebrate <- RF.srvertebrate.resid

gam.srvertebrate <- mgcv::gam(formula = srvertebrate ~  s(sp_resid_srvertebrate, bs = "cr", k = 3) +
                                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                                data = envdata,
                                select = TRUE,
                                family = "gaussian",
                                method = "REML")  

summary(gam.srvertebrate)
saveRDS(gam.srvertebrate, "outputs/SI_IndividualMetricsAnalysis/gam_srvertebrate.rds")


gam.srvertebrate.best <- mgcv::gam(formula = srvertebrate ~ s(sp_resid_srvertebrate, bs = "cr", k = 3) +
                                                            s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                                            s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3),
                                     data = envdata,
                                     select = TRUE,
                                     family = "gaussian",
                                     method = "REML")  
summary(gam.srvertebrate.best)

saveRDS(gam.srvertebrate.best, "outputs/SI_IndividualMetricsAnalysis/gam_srvertebrate_reduce_model.rds")


# 5.RESPONSE VARIABLE = MEAN.FOOD.CHAIN.LENGTH
#---------------------------------------------
auto.MFCL <-  autocov_dist(z = topological_metrics$Mean.food.chain.length, xy = geopoints, longlat = TRUE, nbs = bw, style = "B")     # create the autocavariate variable
colnames(envdata)
envdata <- select(envdata, -c(srvertebrate, sp_resid_srvertebrate))
RF.MFCL <- randomForest(auto.MFCL ~ ., data = envdata, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.MFCL.resid <- auto.MFCL -  RF.MFCL$predicted     #Spatial residuals variable (obs minus fitted)
envdata$MFCL <- topological_metrics$Mean.food.chain.length
envdata$sp_resid_MFCL <- RF.MFCL.resid

gam.MFCL <- mgcv::gam(formula = MFCL ~  s(sp_resid_MFCL, bs = "cr", k = 3) +
                                        s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                        s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                        s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                              data = envdata,
                              select = TRUE,
                              family = "gaussian",
                              method = "REML")  

summary(gam.MFCL)
saveRDS(gam.MFCL, "outputs/SI_IndividualMetricsAnalysis/gam_MFCL.rds")


gam.MFCL.best <- mgcv::gam(formula = MFCL ~ 
                                            s(area, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                            s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                            s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3),
                                   data = envdata,
                                   select = TRUE,
                                   family = "gaussian",
                                   method = "REML")  
summary(gam.MFCL.best)

saveRDS(gam.MFCL.best, "outputs/SI_IndividualMetricsAnalysis/gam_MFCL_reduce_model.rds")


# 6.RESPONSE VARIABLE = MEAN.TOPHIC.LEVEL
#-----------------------------------------
auto.MTL <-  autocov_dist(z = topological_metrics$Mean.trophic.level, xy = geopoints, longlat = TRUE, nbs = bw, style = "B")     # create the autocavariate variable
colnames(envdata)
envdata <- select(envdata, -c(MFCL, sp_resid_MFCL))
RF.MTL <- randomForest(auto.MTL ~ ., data = envdata, na.action = na.roughfix)    #Random forest model (PCA AXIS ~ predictors)
RF.MTL.resid <- auto.MTL -  RF.MTL$predicted     #Spatial residuals variable (obs minus fitted)
envdata$MTL <- topological_metrics$Mean.trophic.level
envdata$sp_resid_MTL <- RF.MTL.resid

gam.MTL <- mgcv::gam(formula = MTL ~  s(sp_resid_MTL, bs = "cr", k = 3) +
                                      s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                      s(Tmean, bs = "cr", k = 3) + s(isothermality, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                      s(PT, bs = "cr", k = 3) + s(NO3, bs = "cr", k = 3) + s(COD, bs = "cr", k = 3),
                      data = envdata,
                      select = TRUE,
                      family = "gaussian",
                      method = "REML")  

summary(gam.MTL)
saveRDS(gam.MTL, "outputs/SI_IndividualMetricsAnalysis/gam_MTL.rds")


gam.MTL.best <- mgcv::gam(formula = MTL ~ s(sp_resid_MTL, bs = "cr", k = 3) +
                                          s(area, bs = "cr", k = 3) + s(max_depth, bs = "cr", k = 3) + s(shorel_dvlpt_index, bs = "cr", k = 3) +
                                          s(Tmean, bs = "cr", k = 3) + s(Tseasonality, bs = "cr", k = 3) +
                                          s(PT, bs = "cr", k = 3),
                           data = envdata,
                           select = TRUE,
                           family = "gaussian",
                           method = "REML")  
summary(gam.MTL.best)

saveRDS(gam.MTL.best, "outputs/SI_IndividualMetricsAnalysis/gam_MTL_reduce_model.rds")
