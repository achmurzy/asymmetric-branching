library(ggnetwork)
library(network)
library(ggpmisc)

#Takes an internal representation and plots the tree as a network
#Takes the name of a value to color nodes in the tree
plot_tree_network <- function(tree, node_val)
{
  #tree[which(tree$CHILD_IDS == "NA"),]$CHILD_IDS <- "_" 
  v1 <- c()
  v2 <- c()
  node_value <- c()
  vol <- c()
  
  for(i in seq(1, nrow(tree)))
  {
    children <- unlist(strsplit(tree[i,]$CHILD_IDS, split="_"))
    if(length(children) > 0)
    {
      for(x in seq(1, length(children)))
      {
        child <- as.numeric(children[x])
        if(!is.na(child))
        {
          v1 <- c(v1, child)
          v2 <- c(v2, i)
          node_value <- c(node_value, tree[child,node_val])
          vol <- c(vol, tree[child,]$VOLUME)
        }
      } 
    }
  }

  edgelist <- matrix(1, nrow=length(v1), ncol=4)
  edgelist[,1] = as.numeric(v1)
  edgelist[,2] = as.numeric(v2) 
  edgelist[,3] = as.numeric(node_value)
  edgelist[,4] = as.numeric(vol)

  #adjacency <- matrix(0, nrow=length(v1)+1, ncol=length(v1)+1)
  #for(x in seq(1, length(v1)))
  #{
  #  adjacency[v1[x], v2[x]] = 1
  #}
  
  n <- network(edgelist[,1:2])
  set.edge.attribute(n, "node_value", edgelist[,3])
  set.edge.attribute(n, "volume", edgelist[,4])
  nn <- ggnetwork(n, layout = 'kamadakawai', niter=1000, weights = "length")
  ggplot(nn, aes(x = x, y = y, xend = xend, yend = yend)) +
         #geom_edges() +
         geom_nodes(aes(color=node_value, size=volume)) +
         theme_blank()
}

example_network <- function()
{
  n <- network(rgraph(100, tprob = 0.2), directed = FALSE)
  n %v% "family" <- sample(letters[1:3], 100, replace = TRUE)
  n %v% "importance" <- sample(1:3, 100, replace = TRUE)
  e <- network.edgecount(n)
  set.edge.attribute(n, "type", sample(letters[24:26], e, replace = TRUE))
  set.edge.attribute(n, "day", sample(1:3, e, replace = TRUE))
  ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(linetype = type), color = "grey50") +
    geom_nodes(aes(color = family, size = importance)) +
    theme_blank()
}

###
# Cylinder model plotting functions
###

#Takes a data frame of the tree cylinders, and optionally,
#computed whole-tree scaling ratios for comparison to the empirical slope
plot_volume <- function(tree, scaling, name)
{
	#form <- y ~ C*(1-exp(k*x))	
	form_line <- y~x
	fit = "tonybananza"
	normal_factor <- min(tree$V_TOT)

	v_tips<-ggplot(data=tree, 
		aes(x=log(V_TOT/normal_factor), y=log(TIPS))) +
	geom_point(alpha=0.5) + scale_y_continuous(limits= c(0, NA)) +
	scale_x_continuous(limits=c(0,NA)) +
	geom_smooth(method = 'lm', formula=form_line, se=FALSE) +
	stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"), outfit=fit<<-..eq.label..),
	       label.x.npc = "left", label.y.npc = 0.7, formula=form_line, parse = TRUE, size = 3)
  intercept = as.numeric(unlist(regmatches(fit,gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*", fit, perl=TRUE))))[1]

	#Comparing computed whole-tree scaling ratios
	if(!missing(scaling))
	{
		v_tips <- v_tips + 
		  geom_abline(slope=scaling$asym, intercept=intercept, colour='#E41A1C') +
		  geom_abline(slope=scaling$wbe, intercept=intercept, colour='#00ff00')
	}	
	#print(v_tips)

	#Test against regression line?
	ggsave(paste(name, "empirical_scaling.png",sep="_"), v_tips, path="figures/")
	return(v_tips)
}

plot_exponents <- function(tree, name)
{
	tree <- subset(tree, !is.na(THETA))
	exponents <- ggplot(data=tree, aes(x=BRANCH_ORDER, y=THETA)) +
	geom_point(alpha=0.5) +
	scale_y_continuous(limits=c(0, 1))
	print(exponents)
	ggsave(paste(name, "subtree_thetas.png",sep="_"), exponents, path="figures/")
}

plot_asymmetry <- function(tree, name)
{
	tree <- subset(tree, !is.na(THETA))
	asymmetry <- ggplot(data=tree, aes(x=D_BETA, y=D_GAMMA, color=THETA)) + scale_x_continuous(limits=c(-1, 1)) + scale_y_continuous(limits=c(-1, 1)) + scale_colour_gradient(low="purple", high="yellow") +
	geom_point(size=1)
	print(asymmetry)
	ggsave(paste(name, "subtree_asymmetry.png",sep="_"), asymmetry, path="figures/")
}
	

plot_symmetry <- function(tree, name)
{
	tree <- subset(tree, !is.na(THETA))
	symmetry <- ggplot(data=tree, aes(x=BETA, y=GAMMA, color=THETA)) + scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1)) + scale_colour_gradient(low="purple", high="yellow") +
	geom_point(size=1)
	print(symmetry)
	ggsave(paste(name, "subtree_symmetry.png",sep="_"), symmetry, path="figures/")
}

plot_lengths <- function(tree, name)
{
	tree <- subset(tree, !is.na(THETA))
	lengths <- ggplot(data=tree, aes(x=BRANCH_ORDER, y=LENGTH)) + 		geom_point(size=1)
	print(lengths)
	ggsave(paste(name, "subtree_lengths.png",sep="_"), lengths, path="figures/")
}

plot_radii <- function(tree, name)
{
	tree <- subset(tree, !is.na(THETA))
	radii <- ggplot(data=tree, aes(x=BRANCH_ORDER, y=RADIUS)) + 		geom_point(size=1)
	print(radii)
	ggsave(paste(name, "subtree_radii.png",sep="_"), radii, path="figures")
}

###
# 
###

###
# All-tree plotting functions
###

plot_empirical_wbe <- function(trees)
{
  ggplot(data = trees, aes(x = WBE, y = EMPIRICAL)) + geom_point(aes(color=SITE)) +
      #scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1))+
    geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 label.x.npc = "left", label.y.npc = 0.7, formula=y~x, parse = TRUE, size = 3)
  ggsave("symmetrical_scaling_all_trees.png", path="figures/")
}

plot_empirical_asym <- function(trees)
{
  ggplot(data = trees, aes(x = THETA, y = EMPIRICAL)) + geom_point(aes(color=SITE)) +
    scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1)) +
    geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 label.x.npc = "left", label.y.npc = 0.7, formula=y~x, parse = TRUE, size = 3)
  ggsave("asymmetrical_scaling_all_trees.png", path="figures/")
}

plot_symmetry_all <- function(trees)
{
  ggplot(data = trees, aes(x = BETA, y = GAMMA)) + geom_point() +
    scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1)) + 
    scale_colour_gradient(low="purple", high="yellow")
  ggsave("symmetry_ratios_all_trees.png", path="figures/")
}

plot_asymmetry_all <- function(trees)
{
  ggplot(data = trees, aes(x = D_BETA, y = D_GAMMA)) + geom_point() +
    scale_x_continuous(limits=c(-1, 1)) + scale_y_continuous(limits=c(-1, 1)) + 
    scale_colour_gradient(low="purple", high="yellow")
  ggsave("asymmetry_ratios_all_trees.png", path="figures/")
  
}

plot_length_scaling <- function(trees)
{
  ggplot(data = trees, aes(x = GAMMA, y = D_GAMMA)) + geom_point() +
    scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(-1, 1)) + 
    scale_colour_gradient(low="purple", high="yellow")
  ggsave("length_scaling_all_trees.png", path="figures/")
}

plot_radius_scaling <- function(trees)
{
  ggplot(data = trees, aes(x = BETA, y = D_BETA)) + geom_point() +
    scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(-1, 1)) + 
    scale_colour_gradient(low="purple", high="yellow")
  ggsave("radius_scaling_all_trees.png", path="figures/")
}

