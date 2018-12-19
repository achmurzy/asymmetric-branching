source("munge.R")
source("analyse.R")

tree_metadata <- read.csv("data/tree_metadata.csv", stringsAsFactors = FALSE)

cylinder_data_names <- list.files(path="data/results/", pattern="^cyl_data_")
stripped_names <- sub(pattern=".*cyl_data_ *(.*?) *.pcd.*", "\\1", x=cylinder_data_names)
cylinder_data_names <- paste("data/results/", cylinder_data_names, sep="")
cylinder_data <- open_TreeQSM(cylinder_data_names[1])
cylinder_data$FILENAME <- rep(stripped_names[1], nrow(cylinder_data))

num_trees = length(cylinder_data_names)

#We are excluding 5 of the 16 potential TreeQSM Tree variables, because these are not always present in the output
tree_data_names <- list.files(path="data/results/", pattern="^tree_data_")
tree_data_names <- paste("data/results/", tree_data_names, sep="")
tree_data <- data.frame(FILENAME=character(num_trees), TREE_VOLUME=numeric(num_trees), TRUNK_VOLUME=numeric(num_trees), 
                        BRANCH_VOLUME=numeric(num_trees), TREE_HEIGHT=numeric(num_trees), TRUNK_LENGTH=numeric(num_trees), 
                    BRANCH_LENGTH=numeric(num_trees), NUM_BRANCHES=integer(num_trees), MAX_BRANCH_ORDER=integer(num_trees), 
                    TREE_AREA=numeric(num_trees), DBH_QSM=numeric(num_trees), DBH_CYL=numeric(num_trees), stringsAsFactors = FALSE)
tree_data[1,] <- c(stripped_names[1], scan(tree_data_names[1])[1:11])

#We might try the entire analysis with different definitions of branch, rather than using
#cylinders which are fine-grained internodes. Ignores biology of bifurcations somewhat (a lot)
#branch_data_names <- list.files(path="data/results/", pattern="^branch_data_")

for(i in seq(2, length(cylinder_data_names)))
{
  print(paste("Working on tree:", stripped_names[i]))
  
  tree_data[i,] <- c(stripped_names[i], scan(tree_data_names[i])[1:11])
  
  #We would need to add code to take output from cylinder analysis and add it to the whole-tree data
  #cyl_data <- open_TreeQSM(cylinder_data_names[i])
  #cyl_data <- branching_analysis(cyl_data)
  #cyl_data$FILENAME <- rep(stripped_names[i], nrow(cyl_data))
  #cylinder_data <- rbind(cylinder_data, cyl_data)
}

#Interact with metadata for whole-tree parameters
joined_data <- merge(tree_data, tree_metadata, by="FILENAME", all=TRUE)
