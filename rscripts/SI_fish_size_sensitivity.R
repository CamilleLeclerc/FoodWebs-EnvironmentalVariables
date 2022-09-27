rm(list = ls()) #Removes all objects from the current workspace (R memory)
mypath <- rprojroot::find_package_root_file


##------------------------------
##LOADING PACKAGES AND FUNCTIONS
##------------------------------
##PACKAGES##
library(dplyr)
library(sjstats)
library(tidyr)


##FUNCTIONS##
source("rfunctions/misc.R")
source("rfunctions/perc_overlap.R")


##DATA##
myload(ind_size, list_lake, taxa_presence, dir = mypath("data")) 




##--------------------------------------------------------------
##"IND_SIZE" AND "TAXA_PRESENCE"/"LIST_LAKE": HARMONIZE DATASETS
##--------------------------------------------------------------
##TAXA
sort(unique(ind_size$species)) #62 species
unique(taxa_presence %>% dplyr::filter(biol.compart == "vertebrate") %>% dplyr::select(taxon) %>% as.data.frame(.)) #46 species
ind_size <- ind_size %>% filter(!species %in% c("Acipenser_ruthenus", "Alburnoides_bipunctatus", "Carassius_gibelio", 
                                                "Chelon_auratus", "Cyprinidae", "Gambusia_affinis", "Gasterosteus_aculeatus", 
                                                "Hypophthalmichthys_molitrix", "Mugilidae", "Neogobius_melanostomus",
                                                "Percidae", "Salvelinus_namaycush", "Hybride_br�me-gardon", "Hybrides_de_cyprinid�s")) 
ind_size$species[ind_size$species == "Abramis"] <- "Abramis_brama"
ind_size$species[ind_size$species == "Salmo_trutta_fario"] <- "Salmo_trutta"
ind_size$species[ind_size$species == "Salmo_trutta_lacustris"] <- "Salmo_trutta"

#SAMPLING YEARS
sort(unique(ind_size$camp_annee)) #2005-2020
ind_size <- ind_size %>% dplyr::filter(camp_annee < 2020)
sort(unique(ind_size$camp_annee)) #2005-2019

#SAMPLING LAKES
length(unique(ind_size$code_lac)) #285 lakes
ind_size <- ind_size %>% dplyr::filter(code_lac %in% list_lake$cd.lac)
length(unique(ind_size$code_lac)) #67 lakes




##----------------------------------------------
##COMPUTING REALIZED FISH BODY SIZE WITHIN LAKES
##----------------------------------------------
realized_ind_size <- ind_size %>%
                      drop_na(.) %>%
                      dplyr::select(species, fish) %>% #code_lac
                      group_by(species) %>% #code_lac
                      summarise_at(vars(fish), list('N' = ~length(.),
                                                      'Mean' = ~ mean(.) %>% round(., digits = 3),
                                                      'Std. Dev.' = ~ se(.) %>% round(., digits = 3),
                                                      'Median' = ~ median(.) %>% round(., digits = 3),
                                                      'Min' = ~ min(.),
                                                      'Max' = ~ max(.)))




##--------------------------------------------------------------------------------
##GETTING NICHE ATTRIBUTES FROM THE ALLOMETRIC NICHE MODEL OF VAGNON ET AL. (2021)
##--------------------------------------------------------------------------------
## see https://doi.org/10.1002/ecs2.3420
## see https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecs2.3420&file=ecs23420-sup-0004-AppendixS4.pdf
## see https://github.com/chloevagnon/aNM_method

source("rfunctions/get_niche_attributes.R") #the function “get_Niche_attributes” from Vagnon et al. (2021) is a modified function from Gravel et al. (2013) to infer niche attributes for vertebrate and invertebate consumers where the argument “consumer_category” must be either “vertebrate” or “invertebrate”.
load("data/Param_regvert.Rdata") # niche attributes for vertebrate consumers get from the allometric niche model of Vagnon et al. (2021)

realized_ind_size$log.mean.size <- log10((realized_ind_size$Mean)*1000) #to mm and to µm
realized_ind_size$log.max.size <- log10((realized_ind_size$Max)*1000) #to mm and to µm
realized_ind_size$log.min.size <- log10((realized_ind_size$Min)*1000) #to mm and to µm

#The niche attributes are: n = consumer size (niche position) / c = center of the consumer feeding range (r) / low and upp = minimal and maximal sizes of preys that determine the consumer feeding range (r)
#Mean body size used
niche_attributes_realized_meansize <- data.frame(matrix(NA, nrow = nrow(realized_ind_size), ncol = 4))
rownames(niche_attributes_realized_meansize) <- realized_ind_size$species #paste(realized_ind_meansize$species, realized_ind_meansize$code_lac, sep = "_")
colnames(niche_attributes_realized_meansize) <- c("n", "c", "low", "upp")
for (j in 1:nrow(realized_ind_size)){ niche_attributes_realized_meansize[j,] <- get_niche_attributes(consumer_size = realized_ind_size$log.mean.size[j], consumer_category = "vertebrate") }
rm(j)


#Max body size used
niche_attributes_realized_maxsize <- data.frame(matrix(NA, nrow = nrow(realized_ind_size), ncol = 4))
rownames(niche_attributes_realized_maxsize) <- realized_ind_size$species
colnames(niche_attributes_realized_maxsize) <- c("n", "c", "low", "upp")
for (j in 1:nrow(realized_ind_size)){ niche_attributes_realized_maxsize[j,] <- get_niche_attributes(consumer_size = realized_ind_size$log.max.size[j], consumer_category = "vertebrate") }
rm(j)


#Min body size used
niche_attributes_realized_minsize <- data.frame(matrix(NA, nrow = nrow(realized_ind_size), ncol = 4))
rownames(niche_attributes_realized_minsize) <- realized_ind_size$species
colnames(niche_attributes_realized_minsize) <- c("n", "c", "low", "upp")
for (j in 1:nrow(realized_ind_size)){ niche_attributes_realized_minsize[j,] <- get_niche_attributes(consumer_size = realized_ind_size$log.min.size[j], consumer_category = "vertebrate") }
rm(j)




##---------------------------------------------------------------------------------
##COMPARING NICHE ATTRIBUTES FROM LITERATURE BASED-BODY SIZE AND REALIZED BODY SIZE
##---------------------------------------------------------------------------------
myload(niche_attributes, taxa_size, dir = mypath("data"))
niche_attributes <- niche_attributes[1:46,]
taxa_size <- taxa_size[1:46,]


#Mean realized body size used
overlap_niche_meansize <- data.frame(matrix(NA, nrow = nrow(niche_attributes), ncol = 12))
colnames(overlap_niche_meansize) <- c("species", "logsize.lit", "n.lit", "c.lit", "low.lit", "upp.lit", "logsize.rea", "n.rea", "c.rea", "low.rea", "upp.rea", "perc_overlap")

for (i in 1:nrow(overlap_niche_meansize)){ 
overlap_niche_meansize[i, 1] <- taxa_size$taxon[i]
overlap_niche_meansize[i, 2] <- taxa_size$log.mean.size[i]
overlap_niche_meansize[i, 3] <- niche_attributes$n[i]
overlap_niche_meansize[i, 4] <- niche_attributes$c[i]
overlap_niche_meansize[i, 5] <- niche_attributes$low[i]
overlap_niche_meansize[i, 6] <- niche_attributes$upp[i]
overlap_niche_meansize[i, 7] <- realized_ind_size$log.mean.size[i]
overlap_niche_meansize[i, 8] <- niche_attributes_realized_meansize$n[i]
overlap_niche_meansize[i, 9] <- niche_attributes_realized_meansize$c[i]
overlap_niche_meansize[i, 10] <- niche_attributes_realized_meansize$low[i]
overlap_niche_meansize[i, 11] <- niche_attributes_realized_meansize$upp[i]
overlap_niche_meansize[i, 12] <- perc_overlap(niche_attributes$low[i], niche_attributes$upp[i], niche_attributes_realized_meansize$low[i], niche_attributes_realized_meansize$upp[i])
}
rm(i)


#Max realized body size used
overlap_niche_maxsize <- data.frame(matrix(NA, nrow = nrow(niche_attributes), ncol = 12))
colnames(overlap_niche_maxsize) <- c("species", "logsize.lit", "n.lit", "c.lit", "low.lit", "upp.lit", "logsize.rea", "n.rea", "c.rea", "low.rea", "upp.rea", "perc_overlap")

for (i in 1:nrow(overlap_niche_maxsize)){ 
  overlap_niche_maxsize[i, 1] <- taxa_size$taxon[i]
  overlap_niche_maxsize[i, 2] <- taxa_size$log.mean.size[i]
  overlap_niche_maxsize[i, 3] <- niche_attributes$n[i]
  overlap_niche_maxsize[i, 4] <- niche_attributes$c[i]
  overlap_niche_maxsize[i, 5] <- niche_attributes$low[i]
  overlap_niche_maxsize[i, 6] <- niche_attributes$upp[i]
  overlap_niche_maxsize[i, 7] <- realized_ind_size$log.max.size[i]
  overlap_niche_maxsize[i, 8] <- niche_attributes_realized_maxsize$n[i]
  overlap_niche_maxsize[i, 9] <- niche_attributes_realized_maxsize$c[i]
  overlap_niche_maxsize[i, 10] <- niche_attributes_realized_maxsize$low[i]
  overlap_niche_maxsize[i, 11] <- niche_attributes_realized_maxsize$upp[i]
  overlap_niche_maxsize[i, 12] <- perc_overlap(niche_attributes$low[i], niche_attributes$upp[i], niche_attributes_realized_maxsize$low[i], niche_attributes_realized_maxsize$upp[i])
}
rm(i)


#Min realized body size used
overlap_niche_minsize <- data.frame(matrix(NA, nrow = nrow(niche_attributes), ncol = 12))
colnames(overlap_niche_minsize) <- c("species", "logsize.lit", "n.lit", "c.lit", "low.lit", "upp.lit", "logsize.rea", "n.rea", "c.rea", "low.rea", "upp.rea", "perc_overlap")

for (i in 1:nrow(overlap_niche_minsize)){ 
  overlap_niche_minsize[i, 1] <- taxa_size$taxon[i]
  overlap_niche_minsize[i, 2] <- taxa_size$log.mean.size[i]
  overlap_niche_minsize[i, 3] <- niche_attributes$n[i]
  overlap_niche_minsize[i, 4] <- niche_attributes$c[i]
  overlap_niche_minsize[i, 5] <- niche_attributes$low[i]
  overlap_niche_minsize[i, 6] <- niche_attributes$upp[i]
  overlap_niche_minsize[i, 7] <- realized_ind_size$log.min.size[i]
  overlap_niche_minsize[i, 8] <- niche_attributes_realized_minsize$n[i]
  overlap_niche_minsize[i, 9] <- niche_attributes_realized_minsize$c[i]
  overlap_niche_minsize[i, 10] <- niche_attributes_realized_minsize$low[i]
  overlap_niche_minsize[i, 11] <- niche_attributes_realized_minsize$upp[i]
  overlap_niche_minsize[i, 12] <- perc_overlap(niche_attributes$low[i], niche_attributes$upp[i], niche_attributes_realized_minsize$low[i], niche_attributes_realized_minsize$upp[i])
}
rm(i)






##---------------------------------------------------------------------------------
##COMPARING FOOD WEB METRICS FROM LITERATURE BASED-BODY SIZE AND REALIZED BODY SIZE
##---------------------------------------------------------------------------------
myload(topological_metrics, topological_metrics_based_on_empirical_data, dir = mypath("outputs"))


overlap_metric <- data.frame(matrix(NA, nrow = ncol(topological_metrics)-1, ncol = 10))
colnames(overlap_metric) <- c("metric", "mean.lit", "sd.lit", "min.lit", "max.lit", "mean.rea", "sd.rea", "min.rea", "max.rea", "perc_overlap")

for (i in 2:ncol(topological_metrics)){ 
  overlap_metric[i-1, 1] <- colnames(topological_metrics)[i]
  overlap_metric[i-1, 2] <- mean(topological_metrics[,i]) %>% round(., digits = 2)
  overlap_metric[i-1, 3] <- sd(topological_metrics[,i]) %>% round(., digits = 2)
  overlap_metric[i-1, 4] <- min(topological_metrics[,i]) %>% round(., digits = 2)
  overlap_metric[i-1, 5] <- max(topological_metrics[,i]) %>% round(., digits = 2)
  overlap_metric[i-1, 6] <- mean(topological_metrics_based_on_empirical_data[,i]) %>% round(., digits = 2)
  overlap_metric[i-1, 7] <- sd(topological_metrics_based_on_empirical_data[,i]) %>% round(., digits = 2)
  overlap_metric[i-1, 8] <- min(topological_metrics_based_on_empirical_data[,i]) %>% round(., digits = 2)
  overlap_metric[i-1, 9] <- max(topological_metrics_based_on_empirical_data[,i]) %>% round(., digits = 2)
  overlap_metric[i-1, 10] <- perc_overlap(overlap_metric[i-1, 4], 
                                         overlap_metric[i-1, 5],
                                         overlap_metric[i-1, 8],
                                         overlap_metric[i-1, 9]) %>% round(., digits = 1)
}
rm(i)






##---------------------------------------------------------------------------
##COMPARING PREY ITEMS FROM LITERATURE BASED-BODY SIZE AND REALIZED BODY SIZE
##---------------------------------------------------------------------------
myload(metaweb, metaweb_based_on_empirical_data, dir = mypath("data"))
sort(colSums(metaweb))
sort(colSums(metaweb_based_on_empirical_data))

class(metaweb)
class(metaweb_based_on_empirical_data)
metaweb_based_on_empirical_data <- as.data.frame(metaweb_based_on_empirical_data)
colnames(metaweb) == colnames(metaweb_based_on_empirical_data)

common_prey <- data.frame(matrix(NA, nrow = 46, ncol = 3))
colnames(common_prey) <- c("Species", "common_prey", "total_prey")

for (i in 1:46){ 

sub_metaweb <- metaweb[,i]
sub_metaweb_based_on_empirical_data <- metaweb_based_on_empirical_data[,i]
intersect <- cbind(sub_metaweb, sub_metaweb_based_on_empirical_data)

common_prey[i, 1] <- colnames(metaweb)[i]
common_prey[i, 2] <- nrow(intersect %>% as.data.frame(.) %>% dplyr::filter(sub_metaweb == 1 & sub_metaweb_based_on_empirical_data == 1))
common_prey[i, 3] <- nrow(intersect %>% as.data.frame(.) %>% dplyr::filter(sub_metaweb == 1 | sub_metaweb_based_on_empirical_data == 1))

}

rm(i)



