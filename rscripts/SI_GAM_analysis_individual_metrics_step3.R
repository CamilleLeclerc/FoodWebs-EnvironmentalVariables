rm(list = ls()) #Removes all objects from the current workspace (R memory)
mypath <- rprojroot::find_package_root_file


##------------------------------------
##LOADING PACKAGES, FUNCTIONS AND DATA
##------------------------------------
library(ggplot2)
library(gratia)
library(cowplot)
library(mgcv)


##FUNCTIONS##
source("rfunctions/misc.R")
source("rfunctions/theme_niwot.R")


##DATA##
gam.dim1 <- readRDS("outputs/gam_dim1.rds")
gam.srtotal <- readRDS("outputs/SI_IndividualMetricsAnalysis/gam_srtotal.rds")
gam.srinvertebrate <- readRDS("outputs/SI_IndividualMetricsAnalysis/gam_srinvertebrate.rds")
gam.LD <- readRDS("outputs/SI_IndividualMetricsAnalysis/gam_LD.rds")

gam.dim2 <- readRDS("outputs/gam_dim2.rds")
gam.srvertebrate <- readRDS("outputs/SI_IndividualMetricsAnalysis/gam_srvertebrate.rds")
gam.MFCL <- readRDS("outputs/SI_IndividualMetricsAnalysis/gam_MFCL.rds")
gam.MTL <- readRDS("outputs/SI_IndividualMetricsAnalysis/gam_MTL.rds")




##-------------
##PLOTTING GAMs
##-------------
p.gam.dim1 <- draw(gam.dim1, scales = "fixed", residuals = TRUE) 
p.gam.srtotal <- draw(gam.srtotal, scales = "fixed", residuals = TRUE)
p.gam.srinvertebrate <- draw(gam.srinvertebrate, scales = "fixed", residuals = TRUE)
p.gam.LD <- draw(gam.LD, scales = "fixed", residuals = TRUE)

p.gam.dim2 <- draw(gam.dim2, scales = "fixed", residuals = TRUE)
p.gam.srvertebrate <- draw(gam.srvertebrate, scales = "fixed", residuals = TRUE)
p.gam.MFCL <- draw(gam.MFCL, scales = "fixed", residuals = TRUE)
p.gam.MTL <- draw(gam.MTL, scales = "fixed", residuals = TRUE)



cowplot::plot_grid(p.gam.dim1[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Lake area (km²)", title = ""),
                   p.gam.srtotal[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Taxa richness", x = "Lake area (km²)", title = ""),
                   p.gam.srinvertebrate[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Invertebrate richness", x = "Lake area (km²)", title = ""),
                   p.gam.LD[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Link density", x = "Lake area (km²)", title = ""),
                   
                   p.gam.dim1[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Shoreline development index", title = ""),
                   p.gam.srtotal[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Taxa richness", x = "Shoreline development index", title = ""),
                   p.gam.srinvertebrate[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Invertebrate richness", x = "Shoreline development index", title = ""),
                   p.gam.LD[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Link density", x = "Shoreline development index", title = ""),
                   
                   p.gam.dim1[[10]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Dissolved organic carbon (mg/L)", title = ""),
                   p.gam.srtotal[[10]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Taxa richness", x = "Dissolved organic carbon (mg/L)", title = ""),
                   p.gam.srinvertebrate[[10]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Invertebrate richness", x = "Dissolved organic carbon (mg/L)", title = ""),
                   p.gam.LD[[10]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Link density", x = "Dissolved organic carbon (mg/L)", title = ""),
                   
                   p.gam.dim1[[6]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Isothermality", title = ""),
                   p.gam.srtotal[[6]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Taxa richness", x = "Isothermality", title = ""),
                   p.gam.srinvertebrate[[6]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Invertebrate richness", x = "Isothermality", title = ""),
                   p.gam.LD[[6]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Link density", x = "Isothermality", title = ""),
                   
                   p.gam.dim1[[9]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Nitrate (mg/L)", title = ""),
                   p.gam.srtotal[[9]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Taxa richness", x = "Nitrate (mg/L)", title = ""),
                   p.gam.srinvertebrate[[9]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Invertebrate richness", x = "Nitrate (mg/L)", title = ""),
                   p.gam.LD[[9]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Link density", x = "Nitrate (mg/L)", title = ""),
                   
                   nrow = 5,
                   align = "hv")



                    
cowplot::plot_grid(p.gam.dim2[[7]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical diversity", x = "Temperature saisonality (°C)", title = ""),
                   p.gam.srvertebrate[[7]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Fish richness", x = "Temperature saisonality (°C)", title = ""),
                   p.gam.MFCL[[7]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Mean food chain length", x = "Temperature saisonality (°C)", title = ""),
                   p.gam.MTL[[7]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Mean trophic level", x = "Temperature saisonality (°C)", title = ""),
                   
                   p.gam.dim2[[5]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical diversity", x = "Annual mean temperature (°C)", title = ""),
                   p.gam.srvertebrate[[5]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Fish richness", x = "Annual mean temperature (°C)", title = ""),
                   p.gam.MFCL[[5]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Mean food chain length", x = "Annual mean temperature (°C)", title = ""),
                   p.gam.MTL[[5]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Mean trophic level", x = "Annual mean temperature (°C)", title = ""),
                   
                   p.gam.dim2[[3]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical diversity", x = "Maximal depth (m)", title = ""),
                   p.gam.srvertebrate[[3]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Fish richness", x = "Maximal depth (m)", title = ""),
                   p.gam.MFCL[[3]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Mean food chain length", x = "Maximal depth (m)", title = ""),
                   p.gam.MTL[[3]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Mean trophic level", x = "Maximal depth (m)", title = ""),
                   
                   p.gam.dim2[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical diversity", x = "Shoreline development index", title = ""),
                   p.gam.srvertebrate[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Fish richness", x = "Shoreline development index", title = ""),
                   p.gam.MFCL[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Mean food chain length", x = "Shoreline development index", title = ""),
                   p.gam.MTL[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Mean trophic level", x = "Shoreline development index", title = ""),
                   
                   p.gam.dim2[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical diversity", x = "Lake area (km²)", title = ""),
                   p.gam.srvertebrate[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Fish richness", x = "Lake area (km²)", title = ""),
                   p.gam.MFCL[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Mean food chain length", x = "Lake area (km²)", title = ""),
                   p.gam.MTL[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "Mean trophic level", x = "Lake area (km²)", title = ""),
                   
                   nrow = 5,
                   align="hv")
#20x16
