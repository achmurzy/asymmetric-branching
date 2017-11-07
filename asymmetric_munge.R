#Takes a computree output csv and transforms to internal frame type
computree_munge <- function(x)
{	
	max_ID = max(x$segment_ID)
	frame <- make_internal_frame(max_ID+1)
	x$segment_children_id <- as.character(x$segment_children_id)
	
	for(i in 1:(max_ID+1))
	{
		cylinders <- subset(x, x$segment_ID == (i-1))
		frame$LENGTH[i] = sum(cylinders$length)
		frame$RADIUS[i] = mean(cylinders$radius)
		frame$VOLUME[i] = 
			pi * frame$LENGTH[i] *(frame$RADIUS[i]^2) 
		frame$BRANCH_ORDER[i] = cylinders$branch_order[1]
		frame$BRANCH_ID[i] = i
		children <- 
	unlist(strsplit(cylinders$segment_children_id[1], split="_"))
		children <- as.numeric(children) + 1
		children <- paste(children, collapse="_")
		frame$CHILD_IDS[i] = children
	}
	return(frame)
}

#Takes an insane-o oxford main-axis tree structure and transforms
#to morphological branch model. Lots of problems aggregating cylinders
#along the length dimension.
oxford_munge <- function(x)
{
	print(paste("BRANCH_IDS_BEFORE: ", max(x[,11])))
	print(paste("DEEP_ORDER_BEFORE: ", max(x[,12])))
	
	frame <- make_internal_frame(nrow(x))
	for(i in 1:nrow(x))
	{
		child_rows = which(i == x[,9])	
		c_ids=""	
		if(x[i,10] == 0)
		{
			if(length(child_rows) > 0)
			{
			c_ids <- paste(child_rows, collapse="_")
			}
		}
		else
		{
			if(length(child_rows) > 1)
			{
			c_ids <- paste(child_rows[-1], collapse="_")
			}
			c_ids <- paste(x[i,10], c_ids,sep="_")
		}
		
		frame$RADIUS[i] = x[i,1]
		frame$LENGTH[i] = x[i,2]
		frame$PARENT_ID[i] = x[i,9]
		frame$CHILD_IDS[i] = c_ids
	}
	result <- oxford_recurse(frame, 1, 1, 1)
	x <- result$data
	max_ID = max(x$BRANCH_ID)
	print(paste("BRANCH_IDS_AFTER: ", max(x$BRANCH_ID)))
	print(paste("DEEP_ORDER_AFTER: ", max(x$BRANCH_ORDER)))
	frame <- make_internal_frame(max_ID)
	for(i in 1:(max_ID))
	{
		cylinders <- subset(x, x$BRANCH_ID == i)
		frame$LENGTH[i] = sum(cylinders$LENGTH)
		frame$RADIUS[i] = mean(cylinders$RADIUS)
		frame$VOLUME[i] = 
			pi * frame$LENGTH[i] *(frame$RADIUS[i]^2) 
		frame$BRANCH_ORDER[i] = cylinders$BRANCH_ORDER[1]
		frame$CHILD_IDS[i]=cylinders$CHILD_IDS[nrow(cylinders)]
	}
	
	return(frame)
}

oxford_recurse <- function(data, row_id, count, order)
{
	#print(paste("Branch ID: ", count))
	data$BRANCH_ID[row_id] = count
	data$BRANCH_ORDER[row_id] = order
	children <- unlist(strsplit(data$CHILD_IDS[row_id], split="_"))
	child_n = length(children)
	if(child_n == 0)
	{ return(list("data"=data, "count"=count)) }
	else if(child_n == 1)
	{
	#print(paste("Branch Row: ", row_id))
	ind = as.numeric(children[1])
	res <- oxford_recurse(data, ind, count, order)
	data <- res$data
	count = res$count
	}		
	else
	{
		#print(paste("Node Row: ", row_id))
		child_ids=c()
		for(i in 1:child_n)
		{
		ind = as.numeric(children[i])
		#print(paste("Recurse child: ", ind))
		child_ids = c(child_ids, count+1)
		res <- oxford_recurse(data, ind, count+1, order+1)
		data <- res$data
		count = res$count
		#print(paste("Stack trace child: ", ind))
		}
		#print(paste("BRANCH_ID: ", data$BRANCH_ID[row_id])) 
		data$CHILD_IDS[row_id] = paste(child_ids, collapse="_")
		#print(paste("CHILD_IDS: ", data$CHILD_IDS[row_id]))
	}
			
	return(list("data"=data, "count"=count))	
}

