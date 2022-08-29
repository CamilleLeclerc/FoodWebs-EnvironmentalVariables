rm(list = ls()) #Removes all objects from the current workspace (R memory)
mypath <- rprojroot::find_package_root_file


##------------------------------
##LOADING PACKAGES AND FUNCTIONS
##------------------------------
##PACKAGES##
library(devtools)
library(factoextra)
library(FactoMineR)
library(ggbiplot)


##FUNCTIONS##
source("rfunctions/misc.R")
source("rfunctions/geom_flat_violin.R")
source("rfunctions/theme_niwot_ggplot2.R")


##DATA##
myload(topological_metrics, dir = mypath("outputs"))
rownames(topological_metrics) <- topological_metrics$cd.lac
topological_metrics <- topological_metrics %>% dplyr::select(-cd.lac)
colnames(topological_metrics)




##-----------------------------
##PRINCIPAL COMPONENTS ANALYSIS
##-----------------------------
res.pca <- PCA(topological_metrics, scale.unit = TRUE)

eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50), ggtheme = theme_classic()) +
  theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black")) +
  theme(axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"))

var <- get_pca_var(res.pca)
var
head(var$coord)
head(var$cos2)
head(var$contrib)

res.desc <- dimdesc(res.pca, axes = c(1:4), proba = 0.05)
res.desc$Dim.1
res.desc$Dim.2

res.pca$ind
lake_pca_coord <- res.pca$ind$coord[,1:4]
mysave(lake_pca_coord, dir = mypath("outputs"), overwrite = TRUE)




data.pca <- prcomp(topological_metrics, scale. = TRUE, retx = -1)
data.pca$rotation <- -1*data.pca$rotation
data.pca$center <- -1*data.pca$center
data.pca$scale <- -1*data.pca$scale
data.pca$x <- -1*data.pca$x
ggbiplot(data.pca, obs.scale = 1, var.scale = 1, circle = TRUE, pointsize = 3) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  theme_bw() + scale_y_continuous(breaks = seq(-5, 5, by = 2.5)) + scale_x_continuous(breaks = seq(-12, 8, by = 4)) +
  theme(axis.title.x = element_text(face = "bold", colour = "black", size = 16, vjust = 0), axis.text.x = element_text(colour = "black", size = 12)) +
  theme(axis.title.y = element_text(face = "bold", colour = "black", size = 16, vjust = 0.2), axis.text.y  = element_text(colour = "black", size = 12))
