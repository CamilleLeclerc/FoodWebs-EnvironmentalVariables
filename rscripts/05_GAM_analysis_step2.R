rm(list = ls()) #Removes all objects from the current workspace (R memory)
mypath <- rprojroot::find_package_root_file


##------------------------------------
##LOADING PACKAGES, FUNCTIONS AND DATA
##------------------------------------
##PACKAGES##
library(mgcv)
library(stringi)
library(dplyr)
library(purrr)


##FUNCTIONS##
source("rfunctions/misc.R")
source("rfunctions/relative_contributions.R")


##DATA##
top_models <- lapply(c("/gam_dim1_reduce_model.rds", "gam_dim2_reduce_model.rds"), function(mods) { readRDS(paste0("outputs/", mods)) })




##-----------------
##GAM SUMMARY TABLE
##-----------------
#Script based on https://doi.org/10.5281/zenodo.596810

model_names <- c("Dim1", "Dim2")

model_tables <- map2(top_models, model_names, function(modd, model_name) {  summ <- summary(modd)
                                                                            summ$p.table
                                                                            summ$s.table
                                                                            rel_dev <- get_relative_contribs(modd)
                                                                            bind_rows(data_frame(Term = stri_extract_first_regex(rownames(summ$p.table), "(?<=\\()[^\\)]+(?=\\))"),
                                                                                                 Value = round(summ$p.table[,1], 3),
                                                                                                 `Z statistic` = round(summ$p.table[,3], 3),
                                                                                                 `Chi-sq statistic` = NA,
                                                                                                 `P-value` = ifelse(summ$p.table[,4] > 0.001, as.character(round(summ$p.table[,4], digits=3)), "<0.001"),
                                                                                                 `Effective Degrees of Freedom` = NA,
                                                                                                 `Total Dev. Explained` = as.character(NA),
                                                                                                 `Relative Dev. Explained` = as.character(NA),
                                                                                                 model = model_name),
                                                                                      data_frame(Term = stri_extract_first_regex(rownames(summ$s.table), "(?<=s\\()[^\\)]+(?=\\))"),
                                                                                                 Value = NA,
                                                                                                 `Z statistic` = NA,
                                                                                                 `Chi-sq statistic` = round(summ$s.table[,3], 3),
                                                                                                 `P-value` = ifelse(summ$s.table[,4] > 0.001, as.character(round(summ$s.table[,4], digits=3)), "<0.001"),
                                                                                                 `Effective Degrees of Freedom` = round(summ$s.table[,1], 3),
                                                                                                 `Total Dev. Explained` = paste0(round(summary(modd)$dev.expl*100, 1), "%"),
                                                                                                 `Relative Dev. Explained` = paste0(round(100*rel_dev$rel_deviance_explained, 1), "%"),
                                                                                                 model = model_name))
  
                                                                        })

model_rows <- map_int(model_tables, nrow)
model_tables2 <- model_tables %>%
  map(~ rbind(.[1,], .)) %>%
  map(function(x) {
    x$model <- c(x$model[1], rep(NA, nrow(x) -1))
    return(x)
  }) %>%
  bind_rows %>%
  mutate_each(funs(as.character), -Term, -model) %>%
  #arrange(model, Term !="Intercept") %>%
  dplyr::select(9, 1:8)

names(model_tables2)[1] <- ""
