open_TreeQSM <- function(file)
{
  return(read.csv(file, sep="\t", header=FALSE))
}

## Parses QSMs output from TreeQSM
## TreeQSM provides many outputs but this parses
## the cylinder data (structure-array)
##  0 row number is implicit ID
##  1 radius
##  2 length
##  3 4 5 start xyz
##  6 7 8 axis xyz
##  9 parent
##  10 extension (child, bifurcations implicit)
##  11 branch ID   
##  12 branch order                                 
##  13 branch position (running number in branch belonging)
##  14 added (after 'normal' cylinder fitting)
##  15 UnmodRadius (not present in some data sets from earlier versions, avoid)

#TreeQSM uses the classical ordering scheme (Hack's stream order), which is based
#on a geometrical definition assigning the lowest order to the longest segment, then
#iterating branch orders based on bifurcations afterward.  
#https://en.wikipedia.org/wiki/Stream_order
#When we parse this data, every child (bifurcation) iterates the order measure - Topological ordering
#We greatly increase the maximum branch ordering by applying this algorithm. 
#Certainly different orderings answer different questions and we may wish to accomodate multiple ones
munge_TreeQSM <- function(x)
{
	#print(paste("BRANCH_IDS_BEFORE: ", max(x[,11])))
	#print(paste("DEEP_ORDER_BEFORE: ", max(x[,12])))
	
	frame <- make_internal_frame(nrow(x))
	max_order = max(x[,12])
	#This loop initializes parent/child relationships and extracts relevant/cylinder.branch parameters
	for(i in 1:nrow(x)) 
	{
		child_rows = which(i == x[,9]) #Get all the cylinders listing this one as a parent
		c_ids=""	
		if(x[i,10] == 0) #Child ID of zero indicates a terminal branch or bifurcation
		{
			if(length(child_rows) > 0)
			{
			  c_ids <- paste(child_rows, collapse="_")
			}
		}
		else  #otherwise just add the extension child, to be collapsed later
		{
			if(length(child_rows) > 1) #if more than one child just add the last one?
			{
			  c_ids <- paste(child_rows[-1], collapse="_")
			}
			c_ids <- paste(x[i,10], c_ids, sep="_")
		}
	
		frame$RADIUS[i] = x[i,1]
		frame$LENGTH[i] = x[i,2]
		frame$PARENT_ID[i] = x[i,9]
		frame$CHILD_IDS[i] = c_ids
		frame$HACK_ORDER[i] = (-1*(x[i,12] - max_order))+1 #Make minimum branch valued at 1
	}
	
	result <- recurse_TreeQSM(frame, 1, 1, 1)
	x <- result$data
	max_ID = max(x$BRANCH_ID)
	
	frame <- make_internal_frame(max_ID)
	
	#Aggregate cylinders according to branch order for a morphologically meaningful model
	#Looks at raw fitted TreeQSM output to see which cylinders we're talking about here
	#Corresponds to the objects in the branch_data output of TreeQSM workflow
	for(i in 1:(max_ID))  
	{
		cylinders <- subset(x, x$BRANCH_ID == i)
		frame$LENGTH[i] = sum(cylinders$LENGTH)
		frame$RADIUS[i] = mean(cylinders$RADIUS)
		frame$VOLUME[i] = pi * frame$LENGTH[i] *(frame$RADIUS[i]^2) 
		frame$BRANCH_ORDER[i] = cylinders$BRANCH_ORDER[1]
		frame$HACK_ORDER[i] = cylinders$HACK_ORDER[1]
		frame$STRAHLER_ORDER[i] = cylinders$STRAHLER_ORDER[1]
		frame$CHILD_IDS[i]=cylinders$CHILD_IDS[nrow(cylinders)]
		children <- as.integer(unlist(strsplit(frame$CHILD_IDS[i], split="_")))
		if(length(children) > 0)
		  frame[children,]$PARENT_ID = i
	}
	return(frame)
}

#This function is what allows us to aggregate the smaller fitted cylinders.
#Uses a counter and a topological ordering assumption to fit together cylinders by a new branch id

#This function is also responsible for constructing alternative branching orders,
#particularly the main BRANCH_ORDER column which corresponds to a simple topological ordering
#iterated at every bifurcation
recurse_TreeQSM <- function(data, row_id, count, order)
{
	data$BRANCH_ID[row_id] = count
	data$BRANCH_ORDER[row_id] = order
	
	children <- unlist(strsplit(data$CHILD_IDS[row_id], split="_"))
	child_n = length(children)
	
	if(child_n == 0)
	{ 
	  data[row_id,]$STRAHLER_ORDER = 1 #Need to think about 1-indexing across branch orderings
	  return(list("data"=data, "count"=count)) 
	}
	else if(child_n == 1)
	{
	  ind = as.numeric(children[1])
	  res <- recurse_TreeQSM(data, ind, count, order)
	  data <- res$data
	  count = res$count
	  data$STRAHLER_ORDER[row_id] = data$STRAHLER_ORDER[ind]
	}		
	else
	{
		child_ids=c()
		for(i in 1:child_n)
		{
  		ind = as.numeric(children[i])
  		child_ids = c(child_ids, count+1)
  		res <- recurse_TreeQSM(data, ind, count+1, order+1)  
  		data <- res$data
  		count = res$count
		}
		data$CHILD_IDS[row_id] = paste(child_ids, collapse="_")
		child_strahler = data$STRAHLER_ORDER[as.numeric(children)]
		first = child_strahler[1]
		if(all(child_strahler == first))
		  data$STRAHLER_ORDER[row_id] = first+1
		else
		  data$STRAHLER_ORDER[row_id] = max(child_strahler)
	}
			
	return(list("data"=data, "count"=count))	
}

branch_TreeQSM <- function(x)
{
  branch_names <- c("ORDER", "PARENT", "VOLUME", "LENGTH", "ANGLE", "HEIGHT", "AZIMUTH", "DIAMETER")  
}

## Parses QSMs from SimpleTree
#Computree uses Horton ordering, we generally prefer topological ordering (see order.png) 
computree_munge <- function(x)
{	
  max_ID = max(x$segment_ID)
  frame <- make_internal_frame(max_ID+1)
  x$segment_children_id <- as.character(x$segment_children_id)
  
  for(i in 1:(max_ID+1))  #Aggregate cylinders according to branch order for a morphologically meaningful model
  {
    cylinders <- subset(x, x$segment_ID == (i-1))
    frame$LENGTH[i] = sum(cylinders$length)
    frame$RADIUS[i] = mean(cylinders$radius)
    frame$VOLUME[i] = pi * frame$LENGTH[i] *(frame$RADIUS[i]^2) 
    frame$BRANCH_ORDER[i] = cylinders$branch_order[1] 
    frame$BRANCH_ID[i] = i
    
    children <- unlist(strsplit(cylinders$segment_children_id[1], split="_"))
    children <- as.numeric(children) + 1
    children <- paste(children, collapse="_")
    frame$CHILD_IDS[i] = children
  }
  return(frame)
}
