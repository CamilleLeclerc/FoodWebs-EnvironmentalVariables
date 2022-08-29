library(hash)
library(igraph)
library(cheddar)

#Generality function from Arnaud's script
Generality <- function(M){ return(sum(colSums(M))/sum((colSums(M) != 0))); }

#Vulnerability function from Arnaud's script
Vulnerability <- function(M){ return(sum(rowSums(M))/sum((rowSums(M) != 0))); }

#SD of generality function from Arnaud's script
SDGenerality <- function(M){ return(sd(colSums(M))); }

#SD of vulnerality function from Arnaud's script
SDVulnerability <- function(M){ return(sd(rowSums(M))); }

InDegree <- Trophic.Generality <- NumberOfResources <- function(M){ return(colSums(M)); }
OutDegree <- Trophic.Vulnerability <- NumberOfCosumers <- function(M){ return(rowSums(M)); }
Degree <- function(M){ return(InDegree(M)+ OutDegree(M)); }

#Fraction of basal taxa function from Arnaud's script
FractionOfBasal <- function(M){ M_temp <- M;
                                diag(M_temp) <- 0;
                                b_sps <- sum(which(InDegree(M_temp) == 0) %in% which(OutDegree(M_temp) >= 1));
                                return(b_sps / dim(M)[1]);
                                }

#Number of basal taxa function from Arnaud's script
NumberOfBasal <- function(M){   M_temp <- M;
                                diag(M_temp) <- 0;
                                b_sps <- sum(which(InDegree(M_temp) == 0) %in% which(OutDegree(M_temp) >= 1));
                                return(b_sps);
                                }
                                
#Fraction of top-level taxa function from Arnaud's script
FractionOfTop <- function(M){   M_temp <- M;
                                diag(M_temp) <- 0;
                                t_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0));
                                return(t_sps / dim(M)[1]);
                                }

#Number of top-level taxa function from Arnaud's script
NumberOfTop <- function(M){     M_temp <- M;
                                diag(M_temp) <- 0;
                                t_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0));
                                return(t_sps);
                                }

#Faction of top-level taxa function from Arnaud's script
FractionOfIntermediate <- function(M){        M_temp <- M;
                                              diag(M_temp) <- 0;
                                              i_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1));
                                              return(i_sps / dim(M)[1]);
                                              }

#Number of intermediate taxa function from Arnaud's script
NumberOfIntermediate <- function(M){          M_temp <- M;
                                              diag(M_temp) <- 0;
                                              i_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1));
                                              return(i_sps);
                                              }

#Mean food chain length function from Arnaud's script
MeanFoodChainLength <- function(M){           # Use PredationMatrixToLinks() to create a Cheddar community from a predation matrix
                                              M <- t(M)
                                              node <- 1:dim(M)[1];
                                              for(n in 1:length(node)){ node[n] <- paste(node[n],'-'); }
                                              pm <- matrix(M, ncol = dim(M)[2], dimnames = list(node, node), byrow = TRUE);
                                              community <- Community(nodes = data.frame(node = node), trophic.links = PredationMatrixToLinks(pm), properties = list(title = 'Community'));
                                              community <- RemoveCannibalisticLinks(community, title = 'community');
                                              #community is a Cheddar community
                                              #community
                                              #NPS(community)
                                              #TLPS(community)
                                              #TrophicLevels(community)
                                              #You can add node properties such as category:
                                              #category <- c('producer', 'invertebrate', 'vert.endo')
                                              #community <- Community(nodes = data.frame(node = node, category = category),
                                              #trophic.links=PredationMatrixToLinks(pm),
                                              #properties=list(title = 'Test community'))
                                              #NPS(community)
                                              #chs <- TrophicChains(community);
                                              #ch_lens <- ChainLength(chs);
                                              chain.stats <- TrophicChainsStats(community)
                                              ch_lens <- (chain.stats$chain.lengths + 1)
                                              return(sum(ch_lens)/length(ch_lens));
                                              }

#Maximum trophic similiarity from Braga et al.'s script [https://github.com/JfvBraga/FoodwebSpace/blob/master/Main%20analysis/Functions.R]
#Function from https://github.com/opetchey/dumping_ground/blob/master/random_cascade_niche/FoodWebFunctions.r
Maxsim <- function(web){ sims <- matrix(0, length(web[,1]), length(web[,1]))
                         for(i in 1:length(web[,1]))
                               for(j in 1:length(web[,1]))
                                     sims[i,j] <- T.sim.ij(web, i, j)
                                     diag(sims) <- NA
                                     mean(apply(sims, 1, function(x) max(x[!is.na(x)])))
                                     }

T.sim.ij <- function(web, i, j){     same <- sum(web[i,] & web[j,]) + sum(web[,i] & web[,j])
                                     total <- sum(web[i,] | web[j,]) + sum(web[,i] | web[,j])
                                     same / total
                                     }
                                     