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


#p.gam.dim1[[1]] -> sp_resid_dim
#p.gam.dim1[[2]] -> lake area
#p.gam.dim1[[3]] -> max depth
#p.gam.dim1[[4]] -> shoreline development index
#p.gam.dim1[[5]] -> Tmean
#p.gam.dim1[[6]] -> isothermality
#p.gam.dim1[[7]] -> Tseasonality
#p.gam.dim1[[8]] -> PT
#p.gam.dim1[[9]] -> NO3
#p.gam.dim1[[10]] -> COD


cowplot::plot_grid( p.gam.dim1[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Lake area (km²)", title = ""),
                    p.gam.dim1[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Shoreline development index", title = ""),
                    p.gam.dim1[[10]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Dissolved organic carbon (mg/L)", title = ""),         
                    p.gam.dim1[[6]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Isothermality", title = ""),
                    p.gam.dim1[[9]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC1 - Trophic diversity", x = "Nitrate (mg/L)", title = ""),
                    labels = "AUTO",
                    nrow = 3,
                    align = "hv")


cowplot::plot_grid( p.gam.dim2[[7]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical diversity", x = "Temperature saisonality (°C)", title = ""),
                    p.gam.dim2[[5]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical diversity", x = "Annual mean temperature (°C)", title = ""),
                    p.gam.dim2[[3]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical diversity", x = "Maximal depth (m)", title = ""),
                    p.gam.dim2[[4]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical diversity", x = "Shoreline development index", title = ""),
                    p.gam.dim2[[2]] + geom_hline(yintercept = 0, linetype = 'dashed', color = "#8c8c8c") + theme_niwot_gam() + labs(y = "PC2 - Vertical diversity", x = "Lake area (km²)", title = ""),
                    labels = "AUTO",
                    nrow = 3,
                    align = "hv")
#SAVE PDF 12 x 8 PORTRAIT
