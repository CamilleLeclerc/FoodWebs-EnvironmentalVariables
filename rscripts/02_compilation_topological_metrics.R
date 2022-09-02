rm(list = ls()) #Removes all objects from the current workspace (R memory)
mypath <- rprojroot::find_package_root_file


##------------------------------
##LOADING PACKAGES AND FUNCTIONS
##------------------------------
##PACKAGES##
library(dplyr)
library(ggplot2)
library(PupillometryR)
library(plyr)
library(stringr)
library(tidyr)


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




##--------------------------------------------
##PLOTTING DITRISBUTION OF TOPOLOGICAL METRICS
##--------------------------------------------
colnames(topological_metrics)
rownames(topological_metrics) <- topological_metrics$cd.lac; topological_metrics$cd.lac <- NULL
colnames(topological_metrics) <- c("Taxa rich.", "Fish rich.", "Invertebr. rich.",
                                    "Body size", "Body mass ratio",
                                    "No. links", "Link dens.", "Connect.",
                                    "Gen.", "Vul.", "Gen. SD", "Vul. SD",
                                    "Frac. bas.", "Frac. int.", "Frac. top",
                                    "Max. Sim.", "MFCL",
                                    "Mean TL", "Max. TL", "Cluster. coef.")

colnames(topological_metrics)
sapply(topological_metrics, summary, na.rm = TRUE)
sapply(topological_metrics, sd, na.rm = TRUE)

data <- topological_metrics %>% gather(key = "text", value = "value")
summary(data)
data <- mutate_at(data, vars(text), as.factor)
data <- data %>% filter(!(text %in% c("Taxa rich.", "Fish rich.", "Invertebr. rich.")))


p <- data %>%
      ggplot(aes(x = text, y = value), color = "#bebebe", fill = "#000000") +
      geom_point(aes(y = value), position = position_jitter(width = 0.25, height = 0, seed = NULL), size = 1, alpha = 0.5) +
      geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.8) +
      PupillometryR::geom_flat_violin(trim = FALSE, position = position_nudge(x = 0.3, y = 0), alpha = 0.8) +
      facet_wrap(~text, ncol = 6, scales = "free") +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.position = "none",
                         axis.text = element_text(size = 16, colour = "#000000"), 
                         axis.line.x = element_line(color = "#000000"), 
                         axis.line.y = element_line(color = "#000000"),
                         strip.text.x = element_text(size = 20, face = "bold"),
                         axis.ticks.x = element_blank(),
                         axis.text.x = element_blank()) +
      labs(y = NULL, x = NULL)
p
