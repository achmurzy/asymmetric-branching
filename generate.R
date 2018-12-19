source("asymmetric_analysis.R")

length_ratio_def = 0.3
length_ratio_sim = seq(0.1, 1, by=0.1)

radius_ratio_def = 0.5
radius_ratio_sim = seq(0.25, 0.75, by=0.05)

length_asym_ratio_def = 0
length_asym_ratio_sim = seq(-0.5, 0.5, 0.1)

radius_asym_ratio_def = 0
radius_asym_ratio_sim = seq(-0.5, 0.5, 0.1)

root_length_def = 1000
root_length_sim = seq(100, 1000, by=100) 

root_width_def = 100
root_width_sim = seq(10, 100, by=10)

leaf_length_def = 10
leaf_length_sim = seq(0.1, 10, by=0.1)

leaf_width_def = 1
leaf_width_sim = seq(0.01, 1, by=0.01)

root <- list(L=root_length_def, R=root_width_def)
petiole <- list(L=leaf_length_def, R=leaf_width_def)

epiphyte_length_def = 1
epiphyte_lingth_sim = seq(1, 100, by=1)

epiphyte_radius_def = 0.1
epiphye_radius_sim = seq(0.1, 10, by = 0.1)

build_child_branches <- function(tree, p_id, counter, pos)
{
  #print(counter)
  #print(tree[counter,])
  #print(typeof(tree[counter,]$LENGTH))
  #print(typeof(tree[counter,]$GAMMA))
  if(pos)
  {
    len = tree[p_id,]$LENGTH * (tree[p_id,]$GAMMA + tree[p_id,]$D_GAMMA)
    rad = tree[p_id,]$RADIUS * (tree[p_id,]$BETA + tree[p_id,]$D_BETA)
  }else
  {
    len = tree[p_id,]$LENGTH * (tree[p_id,]$GAMMA - tree[p_id,]$D_GAMMA)
    rad = tree[p_id,]$RADIUS * (tree[p_id,]$BETA - tree[p_id,]$D_BETA)
  }
  
  gamma = tree[p_id,]$GAMMA
  beta = tree[p_id,]$BETA
  
  d_gamma = tree[p_id,]$D_GAMMA
  d_beta = tree[p_id,]$D_BETA
  
  order = tree[p_id,]$BRANCH_ORDER
  
  #Branch exceeds at least one leaf dimension
  if(len > tree[nrow(tree),]$LENGTH && rad > tree[nrow(tree),]$RADIUS)
  {
      counter = counter + 1
      tree[counter,] <- list(order+1, counter, "", paste(p_id), 0, 0, 0, tree[p_id,]$L_TOT+len, 
                             len, rad, beta, gamma, d_beta, d_gamma, 0, 0, TRUE)
      p_id = counter
         
      c_1 = counter + 1
      res <- build_child_branches(tree, p_id, counter, TRUE)
      tree <- res$data
      counter <- res$counter
      if(counter < c_1)
      { c_1 = "" }
      
      c_2 = counter + 1
      res <- build_child_branches(tree, p_id, counter, FALSE)
      tree <- res$data
      counter <- res$counter
      if(counter < c_2)
      { c_2 = "" }
      
      tree[p_id,]$CHILD_IDS <- paste(c_1, c_2, sep="_")
  }else
  {}
  
  return (list("data"=tree, "counter"=counter))
}

build_tree <- function(trunk=root, leaf=petiole, beta=radius_ratio_def, gamma=length_ratio_def, 
                                  d_beta=radius_asym_ratio_def, d_gamma=length_asym_ratio_def)
{
  counter = 1
  
  #Compute dimensions based on trunk and leaves
  #assumes asymmetry parameters positive and less than symmetry parameters
  if(d_beta < 0 || d_beta > beta || d_gamma < 0 || d_gamma > gamma)
  {
    print("insepct branching parameters for valididty")
    return()
  }
  len <- log(leaf$L/trunk$L)/log(gamma - d_gamma)
  rad <- log(leaf$R/trunk$R)/log(beta - d_beta)
  order <- ceiling(max(len, rad))
  
  #Estimate rows needed based on minimum branch order
  branches <- sum(2 ^ seq(0, order-1))
  tree <- make_internal_frame(branches+1)
  
  tree[1,] <- list(0, 1, "2_3", "0", 0, 0, 0, trunk$L, trunk$L, trunk$R, beta, gamma, d_beta, d_gamma, 0, 0, TRUE)
  tree[nrow(tree),] <- list(1, 1, "", "", 0, 0, 0, leaf$L, leaf$L, leaf$R, beta, gamma, d_beta, d_gamma, 0, 0, FALSE)
  
  c_1 = counter+1
  res <- build_child_branches(tree, 1, counter, TRUE)
  tree <- res$data
  counter <- res$counter
  if(counter < c_1)
  { c_1 = "" }
  
  tree <- res$data
  counter <- res$counter
  
  c_2 = counter+1
  res <- build_child_branches(tree, 1, counter, FALSE)
  tree <- res$data
  counter <- res$counter
  if(counter < c_2)
  { c_2 = "" }
  
  tree[1,]$CHILD_IDS <- paste(c_1, c_2, sep="_")
  tree <- tree[which(tree$POS),]
  return(tree)
}

#Perform centripetal ordering such that all network tips have order 0
centripetal_ordering_scheme <- function(tree)
{
  children <- as.integer(which(tree$CHILD_IDS == "_"))
  tree[children,]$BRANCH_ORDER = 0
  orders <- tree[children,]$BRANCH_ORDER
  zeros <- children != 0
  while(zeros[1])
  {
    children <- as.integer(split_ids(tree[children,]$PARENT_ID))
    zeros <- children != 0
    
    children <- children[zeros]
    orders <- orders[zeros]
    
    tree[children,]$BRANCH_ORDER <- orders + 1
    orders <- tree[children,]$BRANCH_ORDER
  }
  return(tree)
}

#Performs binomial experiment across tree hierarchy with continuous distribution
#Default to Poisson process - exponential distribution - density function only, or
tree_distribution <- function(tree, cd, ...)
{
  size <- range(tree$BRANCH_ORDER)
  tree_mask <- vector(mode="list", length=(size[2]+1))
  names(tree_mask) <- as.character(seq(size[1], size[2]))
  
  for(i in size[1]:size[2])
  {
    rows <- which(tree$BRANCH_ORDER == i)
    n = length(rows)
    pr = do.call(cd, list(i, ...))
    tree_mask[[as.character(i)]] <- rbinom(n, 1, prob=pr)
    #print(pr)
    #print(tree_mask[[as.character(i)]])
  }
  return(tree_mask)
}

#ellipsis operator allows arbitrary number of arguments
distribute_food <- function(tree, cd, ...)
{
  mask <- tree_distribution(tree, cd=cd, ...)
  size <- range(tree$BRANCH_ORDER)
  tree <- cbind(tree, FOOD=-1)
  for(i in size[1]:size[2])
  {
    rows <- which(tree$BRANCH_ORDER == i)
    tree[rows,]$FOOD <- mask[[as.character(i)]]
  }
  return(tree)
}

distribute_nests <- function(tree, cd, ...)
{
  mask <- tree_distribution(tree, cd=cd, ...)
  size <- range(tree$BRANCH_ORDER)
  tree <- cbind(tree, NEST=-1)
  for(i in size[1]:size[2])
  {
    rows <- which(tree$BRANCH_ORDER == i)
    tree[rows,]$NEST <- mask[[as.character(i)]]
  }
  return(tree)
}

#Can we store the "forest" in a separate structure that contains references to tree frames somehow?
#This is why people use pointers... need to have access to my names as numbers...
build_epiphytic_forest <- function(trees, cd, ...)
{
  for(i in nrow(trees))
  {
    mask <- tree_distribution(trees[i,])
  }
}

arboreal_dijkstra <- function(tree, start, end)
{
  if(start == end)
  {return(0)}
  library(hash)
  dist <- hash(seq(1, nrow(tree)), rep(Inf, nrow(tree)))
  current = start
  dist[as.integer(start)] = 0
  
  new_op <- paste(tree[start,]$CHILD_IDS, tree[start,]$PARENT_ID, sep="_")
  new_op <- unlist(strsplit(new_op, split="_"))
  new_op <- trimws(new_op[which(new_op != "" & new_op != "0")])
  options <- as.integer(new_op)
  
  while(!match(end, options, FALSE))
  {
    if(length(new_op) > 0)
    {.set(dist, new_op, dist[[as.character(current)]]+tree[as.integer(new_op),]$LENGTH)}
    current <- as.integer(names(which.min(values(dist[as.character(options)]))))
    
    new_op <- get_connected_branches(tree[current,])
    
    infs <- which(values(dist[new_op]) == Inf)
    new_op <- new_op[match(names(infs), new_op)]
    options <- options[which(options != current)]
    options <- c(as.integer(new_op), options)
  }

  return(dist[[as.character(current)]]+tree[end,]$LENGTH )
}

tree_distances <- function(tree)
{
  nests <- which(as.logical(tree$NEST))
  food <- which(as.logical(tree$FOOD))  
  
  distances <- matrix(nrow=(length(nests)), ncol=length(food))
  dimnames(distances) <- list(nests, food)
  distances <- data.frame(distances)
  colnames <- food
  for(i in nests)
  {
    n <- which(nests == i)
    for(j in food)
    {
      distances[n,which(food == j)] <- arboreal_dijkstra(tree, i, j)
    }
  }
  return(distances)
}

#Takes a list, vector whatever of tree distance matrices and collapses them to plot distribution of paths
plot_distances <- function(tree_distances)
{
  dist <- c(0)
  for(i in seq(1, length(tree_distances)))
  {
    dd <- unlist(tree_distances[[i]], use.names = FALSE)
    dist <- c(dist, dd)
  }
  
  dist <- dist[which(dist != 0)]
  food_nest <- data.frame(LENGTH=numeric(length(dist)))
  food_nest$LENGTH <- dist
  ggplot(food_nest, aes(x=LENGTH)) + geom_histogram(binwidth = 500)
}

plot_tree_stats <- function(tree)
{
  gj <- ggplot(data=tree[which(as.logical(tree$FOOD)),], aes(x=L_TOT)) + geom_histogram()
  print(gj)
}

### Simulation trials

symmetric_uniform_tree <- function(trunk=root, leaf=petiole, beta=radius_ratio_def, gamma=length_ratio_def, 
                                   d_beta=radius_asym_ratio_def, d_gamma=length_asym_ratio_def)
{
   sym <- build_tree(trunk, leaf, beta, gamma, d_beta, d_gamma)
   sym_tree <- centripetal_ordering_scheme(sym)
   maxmin = range(sym_tree$BRANCH_ORDER)
   sym_tree <- distribute_food(sym_tree, dunif, maxmin[1], maxmin[2])
   sym_tree <- distribute_nests(sym_tree, dunif, maxmin[1], maxmin[2])
   return(sym_tree)
}

asymmetric_uniform_tree <- function(trunk=root, leaf=petiole, beta=radius_ratio_def, gamma=length_ratio_def, 
                                    d_beta=0.25, d_gamma=0.25)
{
  asym <- build_tree(trunk, leaf, beta, gamma, d_beta, d_gamma)
  asym_tree <- centripetal_ordering_scheme(asym)
  maxmin = range(asym_tree$BRANCH_ORDER)
  asym_tree <- distribute_food(asym_tree, dunif, maxmin[1], maxmin[2])
  asym_tree <- distribute_nests(asym_tree, dunif, maxmin[1], maxmin[2])
  return(asym_tree)
}

#Extreme cases use geometric distributions to scale probabilities downward from root -> tip or vice versa
symmetric_marginal_tree <- function(trunk=root, leaf=petiole, beta=radius_ratio_def, gamma=length_ratio_def, 
                                   d_beta=radius_asym_ratio_def, d_gamma=length_asym_ratio_def)
{
  sym <- build_tree(trunk, leaf, beta, gamma, d_beta, d_gamma)
  uni_prob = (1/length(unique(sym$BRANCH_ORDER)))
  sym <- distribute_nests(sym, dgeom, uni_prob)
  sym[1,]$NEST <- 1
  sym_tree <- centripetal_ordering_scheme(sym)
  sym_tree <- distribute_food(sym_tree, dgeom, uni_prob)
  return(sym_tree)
}

asymmetric_marginal_tree <- function(trunk=root, leaf=petiole, beta=radius_ratio_def, gamma=length_ratio_def, 
                                    d_beta=0.25, d_gamma=0.25)
{
  asym <- build_tree(trunk, leaf, beta, gamma, d_beta, d_gamma)
  uni_prob = (1/length(unique(asym$BRANCH_ORDER)))
  asym <- distribute_nests(asym, dgeom, uni_prob)
  asym[1,]$NEST <- 1
  asym_tree <- centripetal_ordering_scheme(asym)
  asym_tree <- distribute_food(asym_tree, dgeom, uni_prob)
  return(asym_tree)
}

tree_trial <- function(trials = 10, tree_method, ...)
{
  trees <- vector(mode="list", length=trials)
  distances <- vector(mode="list", length = trials)
  for(i in seq(1, trials))
  {
    tree <- do.call(tree_method, list(...))
    trees [[i]] <- tree
    distances[[i]] <- tree_distances(trees[[i]])
  }
  return(list("trees"=trees, "dist"=distances))
}
