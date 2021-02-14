library(BIEN)
library(reshape2)
library(tidyr)

pca_trait_list <- function()
{
  return(c("EMPIRICAL_VOL", "NORMALISATION", "BETA", "GAMMA", "TIPS_MEAN", "TIPS_VARIANCE", "N_MAX", "BRANCH_ANGLE",
    "leaf.area", "leaf.area.per.leaf.dry.mass", "leaf.carbon.content.per.leaf.dry.mass", 
    "leaf.carbon.content.per.leaf.nitrogen.content",
    "leaf.dry.mass.per.leaf.fresh.mass", "leaf.nitrogen.content.per.leaf.area", "leaf.nitrogen.content.per.leaf.dry.mass", 
    "leaf.phosphorus.content.per.leaf.dry.mass", "leaf.photosynthetic.rate.per.leaf.area", "leaf.photosynthetic.rate.per.leaf.dry.mass",
    "leaf.stomatal.conductance.for.H2O.per.leaf.area", "maximum.whole.plant.height", "seed.mass", "stem.wood.density"))
}

#Takes the vector of species names and returns all traits in BIEN for the species
#Uses nested for-loops to retrieve traits because vectorized trait list not supported
get_tree_species_means <- function(name, trait_list)
{
  traits <- BIEN_trait_mean(name, trait_list[1])[,1:6]

  for(x in seq(2, length(trait_list)))
  {
    print(paste("Querying trait: ", trait_list[x]))
    try(traits <- rbind(traits, BIEN_trait_mean(name, trait_list[x])[,1:6]))
  }
  #Throwing away important information about taxonomic resolution and sample size here
  traits <- traits[,c(1:3)]
  traits <- spread(traits, key=trait, value=mean_value)
  nans <- which(trait_results == 'NaN')
  traits[nans] <- NA
  return(traits)
}

#Data is too sparse for a simple PCA - use a different method
compute_species_means <- function(names)
{
  all_species_traits <- BIEN_trait_species(names)
  all_species_traits$trait_value <- as.numeric(all_species_traits$trait_value)
  trait_means <- aggregate(trait_value ~ scrubbed_species_binomial + trait_name, all_species_traits, FUN = mean, na.action = na.omit)
  melted_means <- melt(trait_means)
  casted_means <- dcast(melted_means, formula=scrubbed_species_binomial~trait_name)
  return(casted_means)
}

get_trait <- function(names, trait)
{
  traits <- BIEN_trait_mean(names, trait)
  traits$mean_value <- as.numeric(as.character(traits$mean_value))
  return(traits)
}

get_leaf_area <- function(names)
{
  traits <- get_trait(names, "leaf area")
  return(traits)
}

get_wood_density <- function(names)
{
  traits <- get_trait(names, "stem wood density")
  return(traits)
}

get_dry_mass <- function(names)
{
  traits <- get_trait(names, "leaf dry mass")
  return(traits)
}

get_sla <- function(names)
{
  traits <- get_trait(names, "leaf area per leaf dry mass")
  return(traits)
}

get_raw_species_traits <- function(names)
{
  trait_list <- BIEN_trait_list()
}