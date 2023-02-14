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
gam.dim2 <- readRDS("outputs/gam_dim2.rds")




##-------------
##PLOTTING GAMs
##-------------
p.gam.dim1 <- draw(gam.dim1, scales = "fixed", residuals = TRUE)
p.gam.dim2 <- draw(gam.dim2, scales = "fixed", residuals = TRUE)


#p.gam.dim1[[1]] -> lake area
#p.gam.dim1[[2]] -> COD
#p.gam.dim1[[3]] -> isothermality
#p.gam.dim1[[4]] -> max depth
#p.gam.dim1[[5]] -> NO3
#p.gam.dim1[[6]] -> PT
#p.gam.dim1[[7]] -> shoreline development index
#p.gam.dim1[[8]] -> sp_resid_dim
#p.gam.dim1[[9]] -> Tmean
#p.gam.dim1[[10]] -> Tseasonality


cowplot::plot_grid( p.gam.dim1[[1]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Lake area (km²)", title = ""),
                    p.gam.dim1[[7]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Shoreline development index", title = ""),
                    p.gam.dim1[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Dissolved organic carbon (mg/L)", title = ""),         
                    p.gam.dim1[[3]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Isothermality", title = ""),
                    p.gam.dim1[[5]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Nitrate (mg/L)", title = ""),
                    labels = "AUTO",
                    nrow = 3,
                    align = "hv")


cowplot::plot_grid( p.gam.dim2[[10]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical structure", x = "Temperature saisonality (°C)", title = ""),
                    p.gam.dim2[[9]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical structure", x = "Annual mean temperature (°C)", title = ""),
                    p.gam.dim2[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical structure", x = "Maximal depth (m)", title = ""),
                    p.gam.dim2[[7]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical structure", x = "Shoreline development index", title = ""),
                    p.gam.dim2[[1]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical structure", x = "Lake area (km²)", title = ""),
                    labels = "AUTO",
                    nrow = 3,
                    align = "hv")
#SAVE PDF 12 x 8 PORTRAIT


cowplot::plot_grid( p.gam.dim1[[10]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Temperature saisonality (°C)", title = ""),
                    p.gam.dim1[[9]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Annual mean temperature (°C)", title = ""),
                    p.gam.dim1[[6]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Phosphore total (mg/L)", title = ""),         
                    p.gam.dim1[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Maximal depth (m)", title = ""),
                    p.gam.dim1[[8]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Spatial residuals", title = ""),
                    labels = "AUTO",
                    nrow = 3,
                    align = "hv")


cowplot::plot_grid( p.gam.dim2[[3]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical structure", x = "Isothermalité", title = ""),
                    p.gam.dim2[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical structure", x = "Dissolved organic carbon (mg/L)", title = ""),
                    p.gam.dim2[[5]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical structure", x = "Nitrate (mg/L)", title = ""),
                    p.gam.dim2[[6]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical structure", x = "Total phosphorus (mg/L)", title = ""),
                    p.gam.dim2[[8]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical structure", x = "Spatial residuals", title = ""),
                    labels = "AUTO",
                    nrow = 3,
                    align = "hv")
