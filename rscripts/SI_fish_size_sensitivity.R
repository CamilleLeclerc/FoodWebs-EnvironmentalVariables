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


##DATA##
myload(ind_size, list_lake, taxa_presence, dir = mypath("data")) 




##--------------------------------------------------------------
##"IND_SIZE" AND "TAXA_PRESENCE"/"LIST_LAKE": HARMONIZE DATASETS
##--------------------------------------------------------------
sort(unique(ind_size$species)) #60 species
unique(taxa_presence %>% dplyr::filter(biol.compart == "vertebrate") %>% dplyr::select(taxon) %>% as.data.frame(.)) #46 species
ind_size <- ind_size %>% filter(!species %in% c("Acipenser_ruthenus", "Alburnoides_bipunctatus", "Carassius_gibelio", 
                                                "Chelon_auratus", "Cyprinidae", "Gambusia_affinis", "Gasterosteus_aculeatus", 
                                                "Hypophthalmichthys_molitrix", "Mugilidae", "Neogobius_melanostomus",
                                                "Percidae", "Salvelinus_namaycush")) 
ind_size$species[ind_size$species == "Abramis"] <- "Abramis_brama"
ind_size$species[ind_size$species == "Salmo_trutta_fario"] <- "Salmo_trutta"
ind_size$species[ind_size$species == "Salmo_trutta_lacustris"] <- "Salmo_trutta"

sort(unique(ind_size$camp_annee)) #2005-2018
length(unique(ind_size$code_lac)) #284 lakes
ind_size <- ind_size %>% dplyr::filter(code_lac %in% list_lake$cd.lac)
length(unique(ind_size$code_lac)) #67 lakes
sort(unique(ind_size$camp_annee)) #2005-2018




##----------------------------------------------
##COMPUTING REALIZED FISH BODY SIZE WITHIN LAKES
##----------------------------------------------
realized_ind_size <- ind_size %>%
                      drop_na(.) %>%
                      dplyr::select(code_lac, species, fish) %>%
                      group_by(code_lac, species) %>%
                      summarise_at(vars(fish), list('N' = ~length(.),
                                                      'Mean' = ~ mean(.) %>% round(., digits = 3),
                                                      'Std. Dev.' = ~ se(.) %>% round(., digits = 3),
                                                      'Median' = ~ median(.) %>% round(., digits = 3),
                                                      'Min' = ~ min(.),
                                                      'Max' = ~ max(.)))




##-------------------------------------------------------------------------------
##GETING NICHE ATTRIBUTES FROM THE ALLOMETRIC NICHE MODEL OF VAGNON ET AL. (2021)
##-------------------------------------------------------------------------------
## see https://doi.org/10.1002/ecs2.3420
## see https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecs2.3420&file=ecs23420-sup-0004-AppendixS4.pdf
## see https://github.com/chloevagnon/aNM_method

source("rfunctions/get_niche_attributes.R") #the function “get_Niche_attributes” from Vagnon et al. (2021) is a modified function from Gravel et al. (2013) to infer niche attributes for vertebrate and invertebate consumers where the argument “consumer_category” must be either “vertebrate” or “invertebrate”.
load("data/Param_regvert.Rdata") # niche attributes for vertebrate consumers get from the allometric niche model of Vagnon et al. (2021)

realized_ind_size$log.mean.size <- log10((realized_ind_size$Mean)*1000) #to mm and to µm

#The niche attributes are: n = consumer size (niche position) / c = center of the consumer feeding range (r) / low and upp = minimal and maximal sizes of preys that determine the consumer feeding range (r)
niche_attributes_realized_size <- data.frame(matrix(NA, nrow = nrow(realized_ind_size), ncol = 4))
rownames(niche_attributes_realized_size) <- paste(realized_ind_size$species, realized_ind_size$code_lac, sep = "_")
colnames(niche_attributes_realized_size) <- c("n", "c", "low", "upp")

for (j in 1:nrow(realized_ind_size)){ niche_attributes_realized_size[j,] <- get_niche_attributes(consumer_size = realized_ind_size$log.mean.size[j], consumer_category = "vertebrate") }
rm(j)

niche_attributes_realized_size$taxon <- realized_ind_size$species
niche_attributes_realized_size$code_lac <- realized_ind_size$code_lac
rownames(niche_attributes_realized_size) <- NULL

