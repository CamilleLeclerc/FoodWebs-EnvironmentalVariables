rm(list = ls()) #Removes all objects from the current workspace (R memory)
mypath <- rprojroot::find_package_root_file


##------------------------------------
##LOADING PACKAGES, FUNCTIONS AND DATA
##------------------------------------
##PACKAGES##
library(cheddar)
library(dplyr)
library(fishualize)
library(ggplot2)
library(gridExtra)
library(igraph)
library(readr)
library(tibble)


##FUNCTIONS##
source("rfunctions/misc.R")
source("rfunctions/export_network3D2.R") #the function “export_network3D2” is modified from Hudson (2018) to export files to Network3D windows software (Yoon et al. 2004).
source("rfunctions/code_network_metrics.R") #the function to calculate different network metrics provided by Arnaud
                                                            #"AveragePPMR" ; "CalculateCommunityStats" ; "CalculateFoodWebStats" ; "CalculatePredatorOverlap" ; "Degree" ; "FindAllPaths" ; "FoodChainStats" ; "FoodChainStatsOld"       
                                                            #"FractionOfBasal" ; "FractionOfCannibalism" ; "FractionOfIntermediate" ; "FractionOfTop" ; "Generality" ; "InDegree" ; "MaximumSimilarity" ; "MeanFoodChainLength"     
                                                            #"MeanGenerality" ; "MeanVulnerability" ; "NormalisedGenerality" ; "NormalisedVulnerability" ; "NormalizeMatrix" ; "NumberOfBasal" ; "NumberOfCosumers" ; "NumberOfIntermediate"    
                                                            #"NumberOfResources" ; "NumberOfTop" ; "Omnivory" ; "OutDegree" ; "SDGenerality" ; "SDVulnerability" ; "TrophicGenerality" ; "TrophicLevels"           
                                                            #"TrophicPositions" ; "TrophicVulnerability" ; "Vulnerability"


##DATA##
myload(list_lake, dir = mypath("data")) ## List of studied lakes
myload(taxa_presence, dir = mypath("data")) ## List of taxa present within each lake
myload(taxa_size, dir = mypath("data")) ## Size information for each taxa
myload(taxa_weight, dir = mypath("data")) ## Weight information for each taxa
myload(niche_attributes, dir = mypath("data")) ## niche attributes for vertebrate and invertebrate consumers get from the allometric niche model of Vagnon et al. (2021)
myload(inferred_pred_prey_links, dir = mypath("data")) ## predator-prey links inferred from the allometric niche model of Vagnon et al. (2021)
myload(metaweb, dir = mypath("data")) ## metaweb gets from the allometric niche model of Vagnon et al. (2021)
#                                     # see https://doi.org/10.1002/ecs2.3420
#                                     # see https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecs2.3420&file=ecs23420-sup-0004-AppendixS4.pdf
#                                     # see https://github.com/chloevagnon/aNM_method




##-----------------------
##FOOD WEB RECONSTRUCTION
##-----------------------
for (i in 1:nrow(list_lake)){ 

  SP <- taxa_presence %>% filter(cd.lac == list_lake[i,])
  SS <- taxa_size %>% filter(taxon %in% SP$taxon)
  SS$length <- as.numeric(as.character(SS$length))
  SS$log.mean.size <- as.numeric(as.character(SS$log.mean.size))

  
  niche_attributes_lake <- niche_attributes %>% filter(rownames(niche_attributes) %in% SS$taxon)
  inferred_pred_prey_links_lake <- inferred_pred_prey_links %>% filter(predator %in% SS$taxon & prey %in% SS$taxon)
  fwmatrix_lake <- metaweb %>% filter(rownames(metaweb) %in% SS$taxon) %>% select(SS$taxon)


  ## Get topological metrics of the foodweb
    #We had to create a cheddar ‘community’ object that is requested for using several functions from the package “cheddar”.
  nodes <- data.frame(SS$taxon, SS$log.mean.size)
  colnames(nodes) <- c("node", "M")

  trophic_links <- inferred_pred_prey_links_lake %>% select(prey, predator)
  colnames(trophic_links) <- c("resource", "consumer")

  property <- list(title = c("foodweb"), M.units = "Log10(μm)")
  
  foodweb <- Community(nodes, properties = property, trophic.links = trophic_links)

  topological_metrics <- as.matrix(c(nrow(nodes), #Species richness
                                     nrow(SS %>% filter(biol.compart == "vertebrate")), #Species richness of vertebrates
                                     nrow(SS %>% filter(biol.compart == "invertebrate"))-2, #Species richness of invertebrates
                                     mean(SS$length, na.rm = TRUE), #Average body size (cm) of the community
                                     mean(inferred_pred_prey_links_lake$mass.ratio), #Average mass ratio of the community
                                     nrow(inferred_pred_prey_links_lake), #Number of links
                                     nrow(inferred_pred_prey_links_lake)/nrow(nodes), #Linkage density // ?LinkageDensity #returns the NumberOfTrophicLinks / NumberOfNodes, including cannibalistic links and isolated nodes ; The number of trophic links in Community
                                     DirectedConnectance(foodweb), #Connectance // DirectedConnectance #returns NumberOfTrophicLinks / NumberOfNodes^2, including cannibalistic links and isolated nodes ; The number of trophic links in Community
                                     Generality(fwmatrix_lake), #Generality	Representing the mean number of prey species per predator
                                     Vulnerability(fwmatrix_lake), #Vulnerability : Representing the mean number of consumer species per prey species
                                     SDGenerality(fwmatrix_lake),
                                     SDVulnerability(fwmatrix_lake),
                                     FractionOfBasal(fwmatrix_lake), #report the connectivity of nodes in a food web - No resources and one or more consumers
                                     FractionOfIntermediate(fwmatrix_lake), #report the connectivity of nodes in a food web - Nodes not fitting any of the above categories (i.e. isolated / basal / top-level)
                                     FractionOfTop(fwmatrix_lake), #report the connectivity of nodes in a food web - One or more resources and no consumers, other than possibly itself
                                     Maxsim(fwmatrix_lake), #Maximum similarity: The average maximum trophic similarity across species in the network
                                     MeanFoodChainLength(fwmatrix_lake), #Mean food chain length: Average length (i.e. number of links) of all the paths (food chains) running from each basal species to each top predator species in the food web
                                     mean(PreyAveragedTrophicLevel(RemoveCannibalisticLinks(foodweb, title='community'))), # Mean trophic level according to Braga et al. 2019 GEB
                                     max(PreyAveragedTrophicLevel(RemoveCannibalisticLinks(foodweb, title='community'))), # Maximum trophic level according to Braga et al. 2019 GEB
                                     transitivity(graph.adjacency(as.matrix(fwmatrix_lake), mode="directed")))) #weighted = "TRUE" ; #Clustering coefficient: Probability of linkage of two species, given that both are linked to a third species

  rownames(topological_metrics) <- c("SR.total", "SR.vertebrate", "SR.invertebrate", "Average.body.size.cm", "Average.mass.ratio.log10", "Number.of.Links", "Linkage.density", "Connectance", 
                                     "Generality", "Vulnerability", "SD.generality", "SD.vulnerability",
                                     "Fraction.basal.nodes", "Fraction.intermediate.nodes", "Fraction.top.nodes",
                                     "Maximum.similarity", "Mean.food.chain.length", "Mean.trophic.level", "Maximum.trophic.level", "Clustering.coefficient")
  colnames(topological_metrics) <- list_lake[i,]
  topological_metrics <- t(topological_metrics)


  ## Exportation the Cheddar community object and the associated species metrics to ‘.web’ and ‘.txt’ files that can be used by the Network 3D Windows software for additionnal investigations.
  links <- cbind(PredatorID = NodeNameIndices(foodweb, TLPS(foodweb)[,'consumer']), PreyID = NodeNameIndices(foodweb, TLPS(foodweb)[,'resource']))
  write.table(links, file.path("outputs/LakeFoodWebs/", paste('network_', list_lake[i,],'.web', sep = '')), row.names = FALSE, sep = ' ')

  write.table(topological_metrics, file.path("outputs/TopologicalMetrics/", paste('topologicalmetrics_', list_lake[i,],'.txt', sep = '')), row.names = TRUE, col.names = TRUE, sep = ' ')

  rm(SP, SS, niche_attributes_lake, fwmatrix_lake, inferred_pred_prey_links_lake, nodes, trophic_links, property, foodweb, topological_metrics, links)

print(i)

}

