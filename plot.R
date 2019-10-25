#library(ggnetwork)
library(ggpmisc)
library(plyr)
library(dplyr)

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
# Cylinder model plotting functions - single trees
###

#Takes a data frame of the tree cylinders
#Plots the original WBE derivation of tip count against network volume
plot_tips <- function(tree, slope)
{
  v_tips<-ggplot(data=tree, aes(x=log(V_TOT), y=log(TIPS))) +
  geom_point(alpha=0.5) +
  xlab("log Subtree Volume") + ylab("log Tip Count") +
  annotate(geom='text', x=0.5, y=2, color="blue", label=round(slope,2)) +
  geom_line(aes(x=log(V_TOT), y=PRED), color="blue")
  return(v_tips)
}

faceted_tips <- function(trees, cylinders)
{
  #Pull slopes from tree data to label them on the plot, and sort them
  cylinders$TIP_FACTOR <- NA
  for(tree in 1:nrow(trees))
  {
    cyls <- which(cylinders$FILENAME == trees[tree,]$FILENAME)
    cylinders[cyls,]$TIP_FACTOR <- trees[tree,]$EMPIRICAL
  }
  cylinders$TIP_FACTOR <- factor(cylinders$TIP_FACTOR, levels=sort(unique(cylinders$TIP_FACTOR)))
  tree_labeller <- function(variable)
  {return("")}
  v_tips<-ggplot(data=cylinders, aes(x=log(V_TOT), y=log(TIPS))) +
  geom_point(alpha=0.05) +
  geom_text(x=0.5, y=4.5, color="blue", hjust=0, vjust=1, aes(label=round(as.numeric(levels(TIP_FACTOR))[TIP_FACTOR],2))) +
  geom_line(aes(x=log(V_TOT), y=PRED), color="blue") +
  xlab("log Subtree Volume") + ylab("log Tip Count") + 
  facet_wrap(~TIP_FACTOR, labeller=tree_labeller)
  return(v_tips)
}

plot_tip_volume <- function(tree, slope)
{
  v_tips<-ggplot(data=tree, aes(x=log(V_TOT), y=log(TIPS_VOL))) +
  geom_point(alpha=0.5) +
  xlab("log Subtree Volume") + ylab("log Distal Tip Volume") +
  annotate(geom='text', x=0.5, y=7, color="green", label=round(slope,2)) +
  geom_line(aes(x=log(V_TOT), y=PRED_VOL), color="green")
  return(v_tips)
}

faceted_tip_volume <- function(trees, cylinders)
{
  #Pull slopes from tree data to label them on the plot, and sort them
  cylinders$VOL_FACTOR <- NA
  for(tree in 1:nrow(trees))
  {
    cyls <- which(cylinders$FILENAME == trees[tree,]$FILENAME)
    cylinders[cyls,]$VOL_FACTOR <- trees[tree,]$EMPIRICAL_VOL
  }
  cylinders$VOL_FACTOR <- factor(cylinders$VOL_FACTOR, levels=sort(unique(cylinders$VOL_FACTOR)))
  tree_labeller <- function(variable)
  {return("")}
  non_tips <- which(cylinders$TIPS > 1)
  v_tips<-ggplot(data=cylinders[non_tips,], aes(x=log(V_TOT), y=log(TIPS_VOL))) +
  geom_point(alpha=0.05) +
  geom_text(x=0.5, y=10, color="green", hjust=0, vjust=1, aes(label=round(as.numeric(levels(VOL_FACTOR))[VOL_FACTOR],2))) +
  geom_line(aes(x=log(V_TOT), y=PRED_VOL), color="green") +
  xlab("log Subtree Volume") + ylab("log Distal Tip Volume") + 
  facet_wrap(~VOL_FACTOR, labeller=tree_labeller)
  return(v_tips)
}

plot_average_volume <- function(tree, slope)
{
  v_tips<-ggplot(data=tree, aes(x=log(VOL_MEAN), y=log(TIPS))) +
  geom_point(alpha=0.5) +
  xlab("log (Subtree Volume / Average Tip Volume)") + ylab("log Tip Count") +
  annotate(geom='text', x=0.5, y=2, color="red", label=round(slope,2)) +
  geom_line(aes(x=log(VOL_MEAN), y=PRED_MA), color="red")
  return(v_tips)
}

faceted_average_volume <- function(trees, cylinders)
{
  #Pull slopes from tree data to label them on the plot, and sort them
  cylinders$AVG_FACTOR <- NA
  for(tree in 1:nrow(trees))
  {
    cyls <- which(cylinders$FILENAME == trees[tree,]$FILENAME)
    cylinders[cyls,]$AVG_FACTOR <- trees[tree,]$EMPIRICAL_VOL_MEAN
  }
  cylinders$AVG_FACTOR <- factor(cylinders$AVG_FACTOR, levels=sort(unique(cylinders$AVG_FACTOR)))
  tree_labeller <- function(variable)
  {return("")}
  non_tips <- which(cylinders$TIPS > 1)
  v_tips<-ggplot(data=cylinders[non_tips,], aes(x=log(VOL_MEAN), y=log(TIPS))) +
  geom_point(alpha=0.05) +
  geom_text(x=0.5, y=4.5, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(AVG_FACTOR))[AVG_FACTOR],2))) +
  geom_line(aes(x=log(VOL_MEAN), y=PRED_MA), color="red") +
  xlab("log (Subtree Volume / Average Tip Volume)") + ylab("log Distal Tips") + 
  facet_wrap(~AVG_FACTOR, labeller=tree_labeller)
  return(v_tips)
}

#Facet the volume plot below with simple regression lines
faceted_volume_plot <- function(trees, cylinders)
{
  cylinders$SLOPE_FACTOR <- NA
  cylinders$TIP_FACTOR <- NA
  cylinders$MEAN_FACTOR <- NA
  cylinders <- cylinders[which(cylinders$TIPS > 1),]
  for(tree in 1:nrow(trees))
  {
    cyls <- which(cylinders$FILENAME == trees[tree,]$FILENAME)
    cylinders[cyls,]$SLOPE_FACTOR <- trees[tree,]$EMPIRICAL_VOL
    cylinders[cyls,]$TIP_FACTOR <- trees[tree,]$EMPIRICAL
    cylinders[cyls,]$MEAN_FACTOR <- trees[tree,]$EMPIRICAL_VOL_MEAN
  }
  cylinders$SLOPE_FACTOR <- factor(cylinders$SLOPE_FACTOR, levels=sort(unique(cylinders$SLOPE_FACTOR)))
  cylinders$TIP_FACTOR <- factor(cylinders$TIP_FACTOR, levels=sort(unique(cylinders$TIP_FACTOR)))
  cylinders$MEAN_FACTOR <- factor(cylinders$MEAN_FACTOR, levels=sort(unique(cylinders$MEAN_FACTOR)))
  ggplot(data=cylinders, aes(x=log(V_TOT), y=log(TIPS_VOL))) + 
    geom_point(alpha=0.01) + 
    #geom_text(x=12, y=2.5, color="blue", hjust=0, vjust=1, aes(label=round(as.numeric(levels(TIP_FACTOR))[TIP_FACTOR],2))) +
    #geom_line(aes(x=log(V_TOT), y=PRED), color="blue") +
    geom_text(x=0, y=12, color="green", hjust=0, vjust=1, aes(label=round(as.numeric(levels(SLOPE_FACTOR))[SLOPE_FACTOR],2))) +
    geom_line(aes(x=log(V_TOT), y=PRED_VOL), color="green") +
    geom_text(x=12, y=2.5, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(MEAN_FACTOR))[MEAN_FACTOR],2))) +
    geom_line(aes(x=log(V_TOT), y=PRED_MA), color="red") + 
    geom_line(aes(x=log(VOL_MEAN), y=PRED_MEAN), color="purple") + 
    facet_wrap(~SLOPE_FACTOR)
}

#Idea for faceting: density plots for all nodes, 
#all nodes stacked by branch ratio (need to reshape to long format)
faceted_scaling_plot <- function(trees)
{
  valid <- which(!trees$INVALID)
  ggplot(trees[valid,]) + geom_density(aes(x=WBE_THETA), color="yellow") + 
  geom_density(aes(x=THETA), color="purple") + 
  geom_vline(xintercept=0.7937, linetype="solid", color="purple") +
  facet_wrap(~FILENAME)
}

faceted_asymmetry_plot <- function(trees)
{
  valid <- which(!trees$INVALID)
  ggplot(trees[valid,]) + geom_density(aes(x=D_BETA), color="orange") + 
  geom_density(aes(x=D_GAMMA), color="blue") + geom_vline(xintercept=0) +
  scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(limits=c(0,5)) + 
  facet_wrap(~FILENAME)
}
#

faceted_symmetry_plot <- function(trees)
{
  valid <- which(!trees$INVALID)
  ggplot(trees[valid,]) + geom_density(aes(x=BETA), color="red") + 
  geom_density(aes(x=GAMMA), color="green") + 
  geom_vline(xintercept = 0.7071, linetype="solid", color="red") + 
  geom_vline(xintercept=0.7937, linetype="solid", color="green") +
  scale_x_continuous(limits=c(0,1)) + facet_wrap(~FILENAME)
}

faceted_volume_scaling_plot <- function(trees)
{
  valid <- which(!trees$INVALID) 
  ggplot(trees[valid,]) + 
  geom_density(aes(x=4*BETA*D_BETA*D_GAMMA + 2*GAMMA*(D_BETA^2), color=POS)) + 
  geom_vline(xintercept=0) +
  scale_x_continuous(limits=c(-0.5,0.5)) + scale_y_continuous(limits=c(0,10)) + 
  facet_wrap(~FILENAME)
}

faceted_asymmetry_contribution_plot <- function(trees)
{
  valid <- which(!trees$INVALID & !is.na(trees$POS)) 
  ggplot(trees[valid,], aes(x=POS,y=4*BETA*D_BETA*D_GAMMA + 2*GAMMA*(D_BETA^2), color=POS)) + 
  geom_histogram(stat="sum") + 
  facet_wrap(~FILENAME)
}

faceted_symmetry_contribution_plot <- function(trees)
{
  valid <- which(!trees$INVALID & !is.na(trees$POS)) 
  ggplot(trees[valid,], aes(x=POS,y=2*(BETA^2)*GAMMA, color=POS)) + 
  geom_histogram(stat="sum") + 
  facet_wrap(~FILENAME)
}

faceted_invalid_nodes_plot <- function(trees)
{
  invalid <- which(trees$INVALID)
  ggplot(trees[valid,]) + geom_density(aes(x=D_BETA), color="orange") + 
  geom_density(aes(x=D_GAMMA), color="blue") + geom_vline(xintercept=0) +
  #scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(limits=c(0,5)) + 
  facet_wrap(~FILENAME)
}

#density lines stacked by branch order, also need to reshape to long format
faceted_branch_order_plot <- function(trees)
{
  ggplot(trees, aes(x=BRANCH_ORDER, y=THETA)) + geom_point(alpha=0.1) + facet_wrap(~FILENAME)
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
	

plot_symmetry <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	symmetry <- ggplot(data=tree, aes(x=BETA, y=GAMMA, color=THETA)) + 
		scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1)) + 
		scale_colour_gradient(low="purple", high="yellow") +
		geom_point(size=1)
	print(symmetry)
	width_i <- 2.13
	height_i <- 1.21
	ggsave("subtree_symmetry.png", plot=symmetry, path="figures/", device=png(), path="figures", width=width_i, height=height_i, units="in", dpi=600, scale=3.5)
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

#Takes a string argument to plot a cylinder variable by branch orderings
plot_branch_ordering <- function(tree, var, name)
{
  library(tidyr)
  gathered <- gather(tree, key=ORDERING, value=ORDER, HACK_ORDER, BRANCH_ORDER, STRAHLER_ORDER, N_GEN, BIN_ORDER)
  hack <- which(gathered$ORDERING == "HACK_ORDER")
  strahler <- which(gathered$ORDERING == "STRAHLER_ORDER")
  gathered[hack,"ORDER"] <- gathered[hack,"ORDER"] + 0.2
  gathered[strahler,"ORDER"] <- gathered[strahler,"ORDER"] + 0.4
  gathered[,var] <- log(gathered[,var])
  gg <- ggplot(gathered, aes_string(x="ORDER", y=var, color="ORDERING", group="ORDERING")) + geom_point(alpha=0.35)
  print(gg)
}

###
# All-node plotting functions
###
#Make sure to feed valid cylinders to these functions for interpretation
plot_asymmetry_all_nodes <- function(cylinders)
{
  gg <- ggplot(cylinders, aes(x=D_BETA, y=D_GAMMA)) + geom_point(aes(color=POS), alpha=0.1)
  print(gg)
  #ggsave('asymmetry_all_nodes.png')
}


###
# All-tree plotting functions
###

faceted_volume_scaling_contributions <- function(trees)
{
  ggplot(trees) + geom_boxplot(aes(y=ASYM_VOL)) + 
  geom_boxplot(aes(y=SYM_VOL)) + facet_wrap(~FILENAME)
}

plot_scaling_comparison <- function(trees)
{
  library(tidyr)
  gathered <- gather(trees, key=SCALING, value=EXPONENT, EMPIRICAL_VOL, WBE, THETA)
  gg <- ggplot(gathered, aes(x=SCALING, y=EXPONENT)) + geom_boxplot() +
  geom_hline(yintercept=0.75, linetype="longdash") +
	scale_x_discrete(labels=c("EMPIRICAL", "ASYMMETRIC", "SYMMETRIC"))
  gg <- gg + labs(title = "Average Scaling Values for All Trees")
  print(gg)
  width_i <- 2.13
  height_i <- 1.21
  #ggsave("scaling_all_trees.png", plot=gg, device=png(), path="figures", width=width_i, height=height_i, units="in", dpi=600, scale=3.5)
}

plot_N_comparison <- function(trees)
{
  library(tidyr)
  gathered <- gather(trees, key=DEPTH, value=N, NETWORK_N, N_ASYM, N_MEAN, N_MIN, N_MAX)
  gg <- ggplot(gathered, aes(x=DEPTH, y=N)) + geom_boxplot()
  #geom_hline(yintercept=0.75, linetype="longdash") +
  #scale_x_discrete(labels=c("SYMMETRICAL", "ASYMMETRICAL", "MEAN", "MIN", "MAX"))
  gg <- gg + labs(title = "Network Depth Algorithms")
  print(gg)
  width_i <- 2.13
  height_i <- 1.21
  #ggsave("scaling_all_trees.png", plot=gg, device=png(), path="figures", width=width_i, height=height_i, units="in", dpi=600, scale=3.5)
}

flattened_species_scaling_comparison <- function(trees)
{
  freq_dist <- plyr::count(trees$SPP)
  spp_counts <- which(freq_dist$freq > 5)
  maxed_spp <- freq_dist[spp_counts,]$x
  min_spp <- freq_dist[which(freq_dist$freq <= 5),]$x
  min_rows <- c()
  for(x in seq(1, length(min_spp)))
  {
    spp <- which(trees$SPP == min_spp[x])
    min_rows <- c(min_rows, spp)
  }
  flattened_trees <- trees[min_rows,]
  for(x in seq(1, length(maxed_spp)))
  {
    spp <- which(trees$SPP == maxed_spp[x])
    flattened_trees <- rbind(flattened_trees, trees[sample(spp, 5),])
  }
  #print(flattened_trees)
  #print(ggplot(flattened_trees, aes(x=SPP)) + geom_histogram(stat="count"))
  plot_scaling_comparison(flattened_trees)
}

plot_empirical_wbe <- function(trees)
{
  gg <- ggplot(data = trees, aes(x = WBE, y = EMPIRICAL)) + geom_point(aes(color=SITE)) +
      geom_errorbar(aes(ymin = CI_MIN, ymax = CI_MAX, color=SITE)) + 
      scale_x_continuous(limits=c(0.2, 0.8)) + scale_y_continuous(limits=c(0.2, 0.8)) +
    geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 label.x.npc = "left", label.y.npc = 0.9, formula=y~x, parse = TRUE, size = 5)
  print(gg)
  #ggsave("symmetrical_scaling_all_trees.png", path="figures/")
}

plot_empirical_node_wbe <- function(trees)
{
  gg <- ggplot(data = trees, aes(x = NODE_WBE, y = EMPIRICAL)) + geom_point(aes(color=SITE)) +
    geom_errorbar(aes(ymin = CI_MIN, ymax = CI_MAX)) + 
    #scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1))+
    geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 label.x.npc = "left", label.y.npc = 0.9, formula=y~x, parse = TRUE, size = 5)
  print(gg)
  #ggsave("symmetrical_node_all_trees.png", path="figures/")
}

plot_empirical_asym <- function(trees)
{
  gg <- ggplot(data = trees, aes(x = THETA, y = EMPIRICAL)) + geom_point(aes(color=SITE)) +
          geom_errorbar(aes(ymin = CI_MIN, ymax = CI_MAX, color=SITE)) +
    scale_x_continuous(limits=c(0.1, 0.8)) + scale_y_continuous(limits=c(0.1, 0.8)) +
    geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 label.x.npc = "left", label.y.npc = 0.9, formula=y~x, parse = TRUE, size = 5)
  gg <- gg + labs(title = "Asymmetry vs. Empirical Scaling", x = "ASYMMETRY")
  #print(gg)
  width_i <- 2.13
  height_i <- 1.21
  ggsave("asymmetrical_scaling_all_trees.png", plot=gg, device=png(), path="figures", width=width_i, height=height_i, units="in", dpi=600, scale=3.5)
}

plot_asymptotic_formula <- function (trees)
{
	nodes <- trees[,c("SYM_VOL", "ASYM_VOL", "NETWORK_N", "THETA")]
	nodes$REAL <- rep(1, nrow(nodes))
	nodes$NEAR <- rep(0.05, nrow(nodes))
	vol_scaling = seq(0.1, 1.5, 0.01)
	network_n = seq(2, 100)

	for(x in vol_scaling)
	{
		for(y in network_n)
		{
			res <- asymptotic_formula(x, y) # Returns NaNs for vol_scaling > 0.9
			if(is.na(res))
			{}
			else
			{
				nearness <- abs(0.75 - res)
				nodes <- rbind(nodes, c(0, x, y, res, 0.5, nearness))
			}
		}
	}
	nodes<- nodes[seq(dim(nodes)[1],1),]

	ggplot(nodes, aes(x=SYM_VOL+ASYM_VOL, y=NETWORK_N)) + 
    geom_point(aes(color=THETA, alpha=(1-NEAR)*REAL, size=1/(NEAR+0.01))) #+
		#scale_color_gradient(low="purple", high="yellow")
		#scale_x_continuous(breaks=seq(0.1, 1.5, by=0.1)) + scale_y_continuous(breaks=seq(2, 50, by=2))
}

plot_symmetrical_volume_scaling <- function (trees)
{
	nodes <- trees[,c("SYM_VOL", "BETA", "GAMMA")]
	nodes$REAL <- rep(1, nrow(nodes))
	nodes$NEAR <- rep(0.05, nrow(nodes))
	beta = seq(0.1, 1, 0.025)
	gamma = seq(0.1, 1, 0.025)

	for(x in beta)
	{
		for(y in gamma)
		{
			res <- 2*(x^2)*y
			nearness <- abs(0.7927 - res)
			nodes <- rbind(nodes, c(res, x, y, 0.5, nearness))
		}
	}
	nodes<- nodes[seq(dim(nodes)[1],1),]

	ggplot(nodes, aes(x=BETA, y=GAMMA)) + 
  geom_point(aes(color=SYM_VOL, alpha=(1-NEAR)*REAL, size=1/(NEAR+0.01))) +
    geom_vline(xintercept=0.707106781) + geom_hline(yintercept=0.793700526) +
		#scale_color_gradient(low="purple", high="yellow")
		scale_x_continuous(breaks=seq(0.1, 1, by=0.1)) + scale_y_continuous(breaks=seq(0.1, 1, by=0.1))
}

plot_asymmetrical_volume_scaling <- function (trees)
{
	nodes <- trees[,c("D_BETA", "D_GAMMA")]
	nodes$ASYM_VOL <- 4*0.7*nodes$D_BETA*nodes$D_GAMMA + 2*0.79*(nodes$D_BETA^2)	
	nodes$REAL <- rep(1, nrow(nodes))
	nodes$NEAR <- abs(nodes$ASYM_VOL)
	d_beta = seq(-0.3, 0.3, 0.005)
	d_gamma = seq(0.0, 0.3, 0.005)

	for(x in d_beta)
	{
		for(y in d_gamma)
		{
			res <- 4*0.7*x*y + 2*0.79*(x^2)
			nearness <- abs(res)
			nodes <- rbind(nodes, c(x, y, res, 0.55, nearness))
		}
	}
	nodes<- nodes[seq(dim(nodes)[1],1),]

	ggplot(nodes, aes(x=D_BETA, y=D_GAMMA)) + geom_point(aes(color=ASYM_VOL, alpha=(1-NEAR)*REAL, size=1/(NEAR+0.05))) +
		#scale_color_gradient(low="purple", high="yellow")
		scale_x_continuous(breaks=seq(-0.5, 0.5, by=0.1)) + scale_y_continuous(breaks=seq(-0.5, 0.5, by=0.1)) +
		scale_color_gradient(low="purple", high="yellow")
}


plot_normalized_asym <- function(trees)
{
  gg <- ggplot(data = trees, aes(x = THETA, y = EMPIRICAL_N_V)) + geom_point(aes(color=SITE)) +
          geom_errorbar(aes(ymin = CI_N_MIN, ymax = CI_N_MAX, color=SITE)) +
    #scale_x_continuous(limits=c(0.1, 0.8)) + scale_y_continuous(limits=c(0.1, 0.8)) +
    geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 label.x.npc = "left", label.y.npc = 0.9, formula=y~x, parse = TRUE, size = 5)
  gg <- gg + labs(title = "Asymmetry vs. Empirical Scaling", x = "ASYMMETRY")
  #print(gg)
  width_i <- 2.13
  height_i <- 1.21
  ggsave("normalized_scaling_all_trees.png", plot=gg, device=png(), path="figures", width=width_i, height=height_i, units="in", dpi=600, scale=3.5)
}

plot_empirical_node_asym <- function(trees)
{
  gg <- ggplot(data = trees, aes(x = NODE_THETA, y = EMPIRICAL)) + geom_point(aes(color=SITE)) +
    geom_errorbar(aes(ymin = CI_MIN, ymax = CI_MAX)) + 
    #scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1))+
    geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 label.x.npc = "left", label.y.npc = 0.9, formula=y~x, parse = TRUE, size = 5)
  print(gg)
  #ggsave("symmetrical_node_all_trees.png", path="figures/")
}

plot_wbe_asym <- function(trees)
{
  gg <- ggplot(joined_data, aes(x=WBE, y=THETA)) + geom_point(aes(color=SITE)) + 
    geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 label.x.npc = "left", label.y.npc = 0.7, formula=y~x, parse = TRUE, size = 3)
  print(gg)
}

plot_empirical_network_n <- function(trees)
{
  gg <- ggplot(joined_data, aes(x=NETWORK_N, y=EMPIRICAL, color=SITE)) + geom_point()
  print(gg)
  #ggsave('network_size_empirical_scaling.png')
}

#BETA:0.707106781
#GAMMA:0.793700526
plot_symmetry_all <- function(trees)
{
  gg <- ggplot(data = trees, aes(x = BETA, y = GAMMA)) + geom_point(aes(color=SITE)) +
    scale_x_continuous(limits=c(0.4, 0.9)) + scale_y_continuous(limits=c(0.4, 0.9)) 
  gg <- gg + labs(title = "Symmetrical Branching Traits across Trees")
  #print(gg)
  width_i <- 2.01
  height_i <- 1.14
  ggsave("symmetrical_ratios_all_trees.png", plot=gg, device=png(), path="figures", width=width_i, height=height_i, units="in", dpi=600, scale=3.5)
  print(gg)
  #ggsave("symmetry_ratios_all_trees.png", path="figures/")
}

plot_asymmetry_all <- function(trees)
{
  gg <- ggplot(data = trees, aes(x = D_BETA, y = D_GAMMA)) + geom_point() +
    scale_x_continuous(limits=c(-1, 1)) + scale_y_continuous(limits=c(-1, 1)) + 
    scale_colour_gradient(low="purple", high="yellow")
  print(gg)
  #ggsave("asymmetry_ratios_all_trees.png", path="figures/")
  
}

plot_asymmetry_proportion <- function(trees)
{
  gg <- ggplot(trees, aes(x=ASYM_FRAC, y=EMPIRICAL, color=SITE)) + geom_point()
  print(gg)
  #ggsave("asymmetry_proportion_all_trees.png")
}

plot_length_asymmetry_path_frac <- function(trees)
{
  gg <- ggplot(joined_data, aes(x=D_GAMMA, y=PATH_FRAC, color=SITE)) + geom_point()
  print(gg)
  #ggsave("length_asymmetry_path_fraction_all_trees.png")
}

plot_radius_asymmetry_path_frac <- function(trees)
{
  gg <- ggplot(joined_data, aes(x=D_BETA, y=PATH_FRAC, color=SITE)) + geom_point()
  print(gg)
  #ggsave("radius_asymmetry_path_fraction_all_trees.png")
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

