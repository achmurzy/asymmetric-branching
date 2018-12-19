library(ggnetwork)
library(network)

#Takes an internal representation and plots the tree as a network
#Takes the name of a value to color nodes in the tree
plot_tree_network <- function(tree, node_val)
{
  #tree[which(tree$CHILD_IDS == "NA"),]$CHILD_IDS <- "_" 
  v1 <- c()
  v2 <- c()
  node_value <- c()
  
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
        }
      } 
    }
  }

  edgelist <- matrix(1, nrow=length(v1), ncol=4)
  edgelist[,1] = as.numeric(v1)
  edgelist[,2] = as.numeric(v2) 
  edgelist[,3] = as.numeric(node_value)

  #adjacency <- matrix(0, nrow=length(v1)+1, ncol=length(v1)+1)
  #for(x in seq(1, length(v1)))
  #{
  #  adjacency[v1[x], v2[x]] = 1
  #}
  
  n <- network(edgelist[,1:2])
  set.edge.attribute(n, "node_value", edgelist[,3])
  nn <- ggnetwork(n, layout = 'kamadakawai', niter=1000, weights = "length")
  ggplot(nn, aes(x = x, y = y, xend = xend, yend = yend)) +
         #geom_edges() +
         geom_nodes(aes(color=node_value)) +
         theme_blank()
}

plot_volumes <- function(ll)
{
	count = 0
	names <- names(ll$species)
	lapply(ll$species, function(x) 
	{  
		count = count + 1
		plot_volume(x, TRUE, names[count]) 
	})
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

#Plots currently heavily weighted by network tips with high volume
plot_volume <- function(tree, single, name)
{
	form <- y ~ C*(1-exp(k*x))	
	form_line <- y~x
	
	normal_factor <- min(tree$V_TOT)
	print(normal_factor)
	#tree <- subset(tree, TIPS > 1 | V_TOT < normal_factor*100)
	v_tips<-ggplot(data=tree, 
		aes(x=log(V_TOT/normal_factor), y=log(TIPS))) +
	geom_point(alpha=0.5) + scale_y_continuous(limits= c(0, NA)) +
	scale_x_continuous(limits=c(0,NA)) +
	geom_smooth(method = 'lm', formula=form_line, se=FALSE) +
	stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., 			sep = "~~~")),label.x.npc = "left", label.y.npc = 0.3,
       		formula=form_line, parse = TRUE, size = 3)
	if(single)
	{
		ss <- mean_theta(summary_stats(tree))
		v_tips <- v_tips + 
		geom_abline(slope=ss$asym, intercept=0, 				colour='#E41A1C') +
		geom_abline(slope=ss$wbe, intercept=0, colour='#00ff00')
	}	
	#print(v_tips)

	#Test against regression line?
	ggsave(paste(name, "Vtot_tips_log.png",sep="_"), v_tips)
}

prediction_plot <- function(tree)
{
	form_line <- y~x
	
	normal_factor <- min(tree$V_TOT)
	tree["Volume_cm^3"]<- log(tree$V_TOT/normal_factor)
	tree["Num_Tips"] <- log(tree$TIPS)
	v_tips<-ggplot(data=tree,
	aes_string(x=quote(tree["Volume_cm^3"]), y=quote(tree["Num_Tips"]))) +
	geom_point(alpha=0.5) + scale_y_continuous(limits= c(0, NA)) +
	scale_x_continuous(limits=c(0,NA)) +
	geom_smooth(method = 'lm', formula=form_line, se=FALSE) + 		geom_abline(slope=0.35, intercept=-1.5, colour='#E41A1C') +
	geom_abline(slope=0.15, intercept=-1.5, colour='#00ff00') 
	
	
	print(v_tips)	
	ggsave(paste("Prediction_Vtot_tips_log.png",sep="_"), v_tips)
}

plot_exponents <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	exponents <- ggplot(data=tree, aes(x=BRANCH_ORDER, y=THETA)) +
	geom_point(alpha=0.5) +
	scale_y_continuous(limits=c(0, 1))
	print(exponents)
	ggsave("subtree_exponent.png", exponents)
}

plot_asymmetry <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	asymmetry <- ggplot(data=tree, aes(x=D_BETA, y=D_GAMMA, color=THETA)) + scale_x_continuous(limits=c(-1, 1)) + scale_y_continuous(limits=c(-1, 1)) + scale_colour_gradient(low="purple", high="yellow") +
	geom_point(size=1)
	print(asymmetry)
	ggsave("subtree_asymmetry.png", asymmetry)
}
	

plot_symmetry <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	symmetry <- ggplot(data=tree, aes(x=BETA, y=GAMMA, color=THETA)) + scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1)) + scale_colour_gradient(low="purple", high="yellow") +
	geom_point(size=1)
	print(symmetry)
	ggsave("subtree_symmetry.png", symmetry)
}

plot_lengths <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	lengths <- ggplot(data=tree, aes(x=BRANCH_ORDER, y=LENGTH)) + 		geom_point(size=1)
	print(lengths)
	ggsave("subtree_lengths.png", lengths)
}

plot_radii <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	radii <- ggplot(data=tree, aes(x=BRANCH_ORDER, y=RADIUS)) + 		geom_point(size=1)
	print(radii)
	ggsave("subtree_radii.png", radii)
}
