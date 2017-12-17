source("asymmetric_analysis.R")

test_trees <- function(trees, decision)
{
  process <- c(0)
  random <- c(0)
  for(i in seq(1, length(trees)))
  {
    process <- c(process, run_ant(trees[[i]], decision))
    #random <- c(random, run_ant(trees[[i]], choose_random))
  }
  
  return(list("process"=process, "random"=random))
}

#Takes list returned from testing above
plot_strategy <- function(result)
{
  result$process <- result$process[which(result$process != 0)]
  df <- data.frame(PROCESS=numeric(length(result$process)))
  df$PROCESS <- result$process
  #df$RANDOM <- result$random
  #df$RANDOM <- df[which(df$RANDOM != 0),]
  print(df)
  ggplot(df, aes(x=PROCESS)) + geom_histogram(binwidth = 500)
  #ggplot(df) + geom_histogram(aes(x=PROCESS), binwidth = 500)  
                      #geom_histogram(aes(x=RANDOM), color="black", binwidth = 1000)
}

# A strongly simplifying assumption would be to allow 
# our ants to have memory of all the branches they've visited
run_ant <- function(tree, decision)
{
  nests <- which(as.logical(tree$NEST))
  food <- which(as.logical(tree$FOOD))
  
  if(length(which(as.logical(tree$NEST))) == 0)
  { print("Warning: no nests on tree") 
    return(0)}
  if(length(which(as.logical(tree$FOOD))) == 0)
  { print("Warning: no food in tree") 
    return (0)}
  
  #print("Nesting sites")
  #print(nests) 
  
  search_lengths <- c(0)
  for(i in seq(1, length(nests)))
  {
    tree$TRAVERSED <- rep(FALSE, nrow(tree))
    index = nests[i]
    tree[index,]$TRAVERSED = TRUE
    #print("Started at nest: ")
    #print(index)
    path_length = 0
    counter = 0
    while(!as.logical(tree[index,]$FOOD) && counter < 1000)
    {
      counter = counter + 1
      index <- decision(tree, index)
      tree[index,]$TRAVERSED = TRUE
      #print("Step index: ")
      #print(index)
      path_length = path_length + tree[index,]$LENGTH
    }
    
    #print("Ended at food:")
    #print(index)
    
    #print("Steps:")
    #print(counter)
    if(counter < 1000)
      {search_lengths <- c(search_lengths, path_length)}
  }
  
  return(search_lengths)
}

#depth-first with memory on tree frame and p
choose_depth <- function(tree, branch)
{
  options <- get_connected_branches(tree[branch,])
  
  if(length(options) == 1) #Parent only
  {
    return(options) 
  }
  else if(length(options == 2)) #Single Child
  {
    if(!tree[options[1],]$TRAVERSED)
    {
      return(options[1])    
    }
    else
    {
      return(options[2]) 
    }
  } 
  else                  #Full structure
  {
    if(!tree[options[1],]$TRAVERSED && !tree[options[2],]$TRAVERSED)
    {
      return(options[which.min(tree[options,]$LENGTH)])  
    }
    else if(!tree[options[1],]$TRAVERSED)
    {
      return(options[1])
    }
    else if(!tree[options[2],]$TRAVERSED)
    {
      return(options[2])
    }
    else
    {
      return(options[3])
    }
  }
}

#These are both quite clearly depth-first search given
#the way branching parameters currently operate
choose_skinny <- function(tree, branch)
{
  options <- get_connected_branches(tree[branch,])
  print(options)
  if(length(options) == 1)
    return(options)
  else
    return(options[which.min(tree[options,]$RADIUS)])
}
choose_short <- function(tree, branch)
{
  options <- get_connected_branches(tree[branch,])
  if(length(options) == 1)
    return(options)
  else
    return(options[which.min(tree[options,]$LENGTH)])
}

choose_random <- function(tree, branch)
{
  options <- get_connected_branches(tree[branch,])
  if(length(options) == 1)
    return(options)
  else
    return(as.integer(sample(options, 1)))
}