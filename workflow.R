#Sample file names for exploratory analysis/debugging of single tree:
#cyl_data_Ery_01.pcd_t1_m1_D0.15_DA0.06_DI0.02_L3_F3.txt
#cyl_data_Pinus03.pcd_t1_m1_D0.15_DA0.06_DI0.02_L3_F3.txt
#
#test <- open_TreeQSM("data/results/cyl_data_Ery_01.pcd_t1_m1_D0.15_DA0.06_DI0.02_L3_F3.txt")
#test <- munge_TreeQSM(test)
#test_out <- branching_analysis(test)

#test_branch <- open_TreeQSM("data/results/branch_data_Ery_01.pcd_t1_m1_D0.15_DA0.06_DI0.02_L3_F3.txt")

source("munge.R")
source("analyse.R")
source("plot.R")

tree_metadata <- read.csv("data/tree_metadata.csv")
tree_metadata$FILENAME <- as.character(tree_metadata$FILENAME)
excluded_fits <- scan("data/excluded.txt", what=character())

###
# This workflow loops through a directory of fit QSM cylinders and combines all of the branch-level and whole-tree
# data with available metadata on the trees.
###

branch_data_names <- list.files(path="data/results/", pattern="^branch_data_")
stripped_names <- sub(pattern=".*branch_data_ *(.*?) *.pcd.*", "\\1", x=branch_data_names)
branch_data_names <- paste("data/results/", branch_data_names, sep="")
branch_data <- open_TreeQSM(branch_data_names[1])
branch_column_names = c("ORDER", "PARENT", "VOLUME", "LENGTH", "ANGLE", "HEIGHT", "AZIMUTH", "DIAMETER")

#For some reason we like going through the first loop iteration manually to initialize everything properly
cylinder_data_names <- list.files(path="data/results/", pattern="^cyl_data_")
stripped_names <- sub(pattern=".*cyl_data_ *(.*?) *.pcd.*", "\\1", x=cylinder_data_names)
cylinder_data_names <- paste("data/results/", cylinder_data_names, sep="")
cylinder_data <- open_TreeQSM(cylinder_data_names[1])
cylinder_data <- munge_TreeQSM(cylinder_data)

#guyana_names <- c("GUY01", "GUY02", "GUY04", "GUY05", "GUY06", "GUY07", "GUY09")
#wytham_names <- c("wytham_winter_10", "wytham_winter_45", "wytham_winter_70", "wytham_winter_81", "wytham_winter_195", 
#                  "wytham_winter_229", "wytham_winter_291", "wytham_winter_424", "wytham_winter_446", "wytham_winter_600",
#                  "wytham_winter_869", "wytham_winter_990", "wytham_winter_1551", "wytham_winter_1777", "wytham_winter_2225")

#The branching analysis outputs the external form of the tree data for whole-tree scaling analysis and trait reconciliation
out <- branching_analysis(cylinder_data)
cyl_data <- out$branches
scaling <- out$scaling
cylinder_data$FILENAME <- rep(stripped_names[1], nrow(cylinder_data))

num_trees = length(cylinder_data_names)

#We are excluding 5 of the 16 potential TreeQSM Tree variables, because these are not always present in the output
tree_data_names <- list.files(path="data/results/", pattern="^tree_data_")
tree_data_names <- paste("data/results/", tree_data_names, sep="")

### !!!
#EXTREMELY IMPORTANT NOTE: 
#IF YOU ADD A VARIABLE TO BE RETURNED BY THE SCALING ANALYSIS, IT HAS TO BE SPECIFIED HERE IN THE EXACT SAME ORDER !!!
### !!!
tree_data <- data.frame(FILENAME=character(num_trees), TREE_VOLUME=numeric(num_trees), TRUNK_VOLUME=numeric(num_trees),
                        BRANCH_VOLUME=numeric(num_trees), TREE_HEIGHT=numeric(num_trees), TRUNK_LENGTH=numeric(num_trees), 
                    BRANCH_LENGTH=numeric(num_trees), NUM_BRANCHES=integer(num_trees), MAX_BRANCH_ORDER=integer(num_trees), 
                    TREE_AREA=numeric(num_trees), DBH_QSM=numeric(num_trees), DBH_CYL=numeric(num_trees), BETA=numeric(num_trees),
                    GAMMA=numeric(num_trees), D_BETA=numeric(num_trees), D_GAMMA=numeric(num_trees), FIB_R=numeric(num_trees),
                    FIB_L=numeric(num_trees), WBE=numeric(num_trees), NODE_WBE=numeric(num_trees), THETA=numeric(num_trees), 
                    NODE_THETA=numeric(num_trees), EMPIRICAL=numeric(num_trees), NETWORK_N=integer(num_trees), 
                    CI_MIN=numeric(num_trees), CI_MAX=numeric(num_trees),
                    PATH_FRAC=numeric(num_trees), ASYM_FRAC=numeric(num_trees), stringsAsFactors = FALSE)
tree_data[1,] <- c(stripped_names[1], scan(tree_data_names[1])[1:11], as.numeric(unlist(scaling)))
tree_column_names <- c("TREE_VOLUME", "TRUNK_VOLUME", "BRANCH_VOLUME", "TREE_HEIGHT", "TRUNK_LENGTH", "BRANCH_LENGTH", "NUM_BRANCHES", 
                       "MAX_BRANCH_ORDER", "TREE_AREA", "DBH_QSM", "DBH_CYL")
###
# Begin Looping
###

for(i in seq(2, length(cylinder_data_names)))
{
  if(stripped_names[i] %in% excluded_fits)
  {
    print(paste("Skipping poorly fit tree:", stripped_names[i]))
    next
  }
  else
    print(paste("Working on tree:", stripped_names[i]))
  
  cyl_data <- open_TreeQSM(cylinder_data_names[i])
  cyl_data <- munge_TreeQSM(cyl_data)
  out <- branching_analysis(cyl_data, verbose_report = FALSE)
  cyl_data <- out$branches
  scaling <- out$scaling
  cyl_data$FILENAME <- rep(stripped_names[i], nrow(cyl_data))
  cylinder_data <- rbind(cylinder_data, cyl_data)
  
  #19 columns as specified above
  tree_data[i,] <- c(stripped_names[i], scan(tree_data_names[i])[1:11], as.numeric(unlist(scaling)))
  #tree_data[i,13:18] <- unlist(scaling)
  
  #plot_volume(cyl_data, scaling, stripped_names[i])
  #plot_exponents(cyl_data, scaling, stripped_names[i])
  #plot_asymmetry(cyl_data, stripped_names[i])
  #plot_symmetry(cyl_data, stripped_names[i])
  #plot_lengths(cyl_data, stripped_names[i])
  #plot_radii(cyl_data, stripped_names[i])
}

#cast everything back to its proper type because R is horrible
tree_data$NUM_BRANCHES <- as.integer(tree_data$NUM_BRANCHES)
tree_data$MAX_BRANCH_ORDER <- as.integer(tree_data$MAX_BRANCH_ORDER)
tree_data[,3:ncol(tree_data)] <- sapply(tree_data[,3:ncol(tree_data)], as.numeric)

#Manually add Guyanese trees to tree_data

#Join with metadata for whole-tree parameters
#With the all flag, filter trees not present in the metadata, and metadata not present in the trees
joined_data <- merge(tree_data, tree_metadata, by="FILENAME", all=FALSE)
joined_data <- joined_data[which(joined_data$FILENAME != ""),]

#Finite size correction is super tiny - within 95% of 3/4
smallest = which(min(joined_data$TREE_VOLUME, na.rm=TRUE) == joined_data$TREE_VOLUME)
small_tips = exp(joined_data[smallest,]$NETWORK_N * log(2))
largest = which(max(joined_data$TREE_VOLUME, na.rm=TRUE) == joined_data$TREE_VOLUME)

correction = 1 - (2^(-1/3) * small_tips^(-1/3) / log(joined_data[largest,]$TREE_VOLUME / joined_data[smallest,]$TREE_VOLUME))


source("traits.R")
species_names <- paste(joined_data$GENUS, joined_data$SPECIES)
unique_names <- unique(species_names)
#species_traits <- get_tree_species_traits(species_name)

#Be wary of the number of family-level trait means you're getting back from BIEN
leaf_area <- get_leaf_area(unique_names)
wood_density <- get_wood_density(unique_names)
joined_data$LEAF_AREA <- as.numeric(rep(0, nrow(joined_data)))
joined_data$WOOD_DENSITY <- as.numeric(rep(0, nrow(joined_data)))

for(x in seq(1, nrow(joined_data)))
{
  joined_data[x,]$LEAF_AREA <- leaf_area[which(leaf_area$species == species_names[x]),]$mean_value
  joined_data[x,]$WOOD_DENSITY <- wood_density[which(wood_density$species == species_names[x]),]$mean_value
}

#ggplot(joined_data, aes(x=WBE, y=LEAF_AREA)) + geom_point()
#ggplot(joined_data, aes(x=WBE, y=WOOD_DENSITY)) + geom_point()
