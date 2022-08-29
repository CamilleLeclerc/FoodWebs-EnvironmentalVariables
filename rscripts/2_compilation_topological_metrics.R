rm(list = ls()) #Removes all objects from the current workspace (R memory)
mypath <- rprojroot::find_package_root_file


##------------------------------
##LOADING PACKAGES AND FUNCTIONS
##------------------------------
##PACKAGES##
library(plyr)
library(stringr)


##FUNCTIONS##
source("rfunctions/misc.R")




##----------------------------------------------
##COMPILING A UNIQUE FILE OF TOPOLOGICAL METRICS
##----------------------------------------------
mydir <- "outputs/TopologicalMetrics"
myfiles <- list.files(path = mydir, pattern = "*.txt", full.names = TRUE)
topological_metrics <- ldply(myfiles, read.table, sep = "", fill = TRUE, header = TRUE)
cd.lac <- as.data.frame(str_match(myfiles, "_\\s*(.*?)\\s*.txt")[,2])
colnames(cd.lac) <- "cd.lac"
topological_metrics <- cbind(cd.lac, topological_metrics)
mysave(topological_metrics, dir = mypath("outputs"), overwrite = TRUE)
