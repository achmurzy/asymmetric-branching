library(BIEN)

#Takes the vector of species names and returns all traits in BIEN for the species
get_tree_species_traits <- function(names)
{
  trait_list <- BIEN_trait_list()
  traits <- BIEN_trait_mean(names, trait_list[1])
  for(x in seq(2, length(trait_list)))
  {
    traits <- rbind(traits, BIEN_trait_mean(names, trait_list[x]))
  }
  traits$mean_value <- as.numeric(as.character(traits$mean_value))
  return(traits)
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