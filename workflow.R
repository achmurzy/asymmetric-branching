#Sample file names for exploratory analysis/debugging of single tree:
#cyl_data_Ery_01.pcd_t1_m1_D0.15_DA0.06_DI0.02_L3_F3.txt
#cyl_data_Pinus03.pcd_t1_m1_D0.15_DA0.06_DI0.02_L3_F3.txt

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

#For some reason we like going through the first loop iteration manually to initialize everything properly
cylinder_data_names <- list.files(path="data/results/", pattern="^cyl_data_")
stripped_names <- sub(pattern=".*cyl_data_ *(.*?) *.pcd.*", "\\1", x=cylinder_data_names)
cylinder_data_names <- paste("data/results/", cylinder_data_names, sep="")
cylinder_data <- open_TreeQSM(cylinder_data_names[1])
cylinder_data <- munge_TreeQSM(cylinder_data)
out <- branching_analysis(cylinder_data)
cyl_data <- out$branches
scaling <- out$scaling
cylinder_data$FILENAME <- rep(stripped_names[1], nrow(cylinder_data))

num_trees = length(cylinder_data_names)

#We are excluding 5 of the 16 potential TreeQSM Tree variables, because these are not always present in the output
tree_data_names <- list.files(path="data/results/", pattern="^tree_data_")
tree_data_names <- paste("data/results/", tree_data_names, sep="")
tree_data <- data.frame(FILENAME=character(num_trees), TREE_VOLUME=numeric(num_trees), TRUNK_VOLUME=numeric(num_trees), 
                        BRANCH_VOLUME=numeric(num_trees), TREE_HEIGHT=numeric(num_trees), TRUNK_LENGTH=numeric(num_trees), 
                    BRANCH_LENGTH=numeric(num_trees), NUM_BRANCHES=integer(num_trees), MAX_BRANCH_ORDER=integer(num_trees), 
                    TREE_AREA=numeric(num_trees), DBH_QSM=numeric(num_trees), DBH_CYL=numeric(num_trees), BETA=numeric(num_trees),
                    GAMMA=numeric(num_trees), D_BETA=numeric(num_trees), D_GAMMA=numeric(num_trees), WBE=numeric(num_trees),
                    THETA=numeric(num_trees), EMPIRICAL=numeric(num_trees), stringsAsFactors = FALSE)
tree_data[1,] <- c(stripped_names[1], scan(tree_data_names[1])[1:11], unlist(scaling))

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

#Join with metadata for whole-tree parameters
#With the all flag, filter trees not present in the metadata, and metadata not present in the trees
joined_data <- merge(tree_data, tree_metadata, by="FILENAME", all=FALSE)