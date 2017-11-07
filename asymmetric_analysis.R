library(gdata)
library(ggplot2)
library(ggpmisc)
library(reshape)
library(smatr)


make_internal_frame <- function(rows)
{
frame <- data.frame(BRANCH_ORDER=integer(rows), BRANCH_ID=integer(rows), CHILD_IDS=character(rows), TIPS=integer(rows), V_TOT=numeric(rows), VOLUME = numeric(rows), LENGTH=numeric(rows), RADIUS=numeric(rows), BETA=numeric(rows), GAMMA=numeric(rows), D_BETA=numeric(rows), D_GAMMA=numeric(rows), THETA=numeric(rows), WBE_THETA=numeric(rows), POS=logical(rows), stringsAsFactors=FALSE)
return(frame)
}

#Aggregating cylinders by length is causing problems
#Need a more intelligent way to think about this
asymmetric_branching <- function(id, branches, c_v, c_t)
{
	branch <- branches[id,]
	children <- unlist(strsplit(branch$CHILD_IDS[1], split="_"))
	child_n = length(children) 
	if(child_n == 0)
	{
		branch$TIPS = 1
		branch$V_TOT = branch$VOLUME
		branch$BETA = NA	
		branch$GAMMA = NA
		branch$D_BETA = NA
		branch$D_GAMMA = NA
		branch$WBE_THETA = NA
		branch$THETA = NA
		branches[id,] <- branch
		return (list("data"=branches, "c_v"=c_v, "c_t"=c_t))
	}
	else if (child_n == 1)
	{
		ind = as.numeric(children[1])
		out <- asymmetric_branching(ind, branches, c_v, c_t)
		branches = out$data
		c_v = out$c_v
		c_t = out$c_t
		
		branch$BETA = branches[ind,]$RADIUS/branch$RADIUS
		branch$GAMMA = branches[ind,]$LENGTH/branch$LENGTH

		branch$TIPS = branches[ind,]$TIPS
		branch$V_TOT = branches[ind,]$V_TOT + branch$VOLUME

		branch$D_BETA = 0
		branch$D_GAMMA = 0
	}
	else
	{
		tips = 0
		vol = 0
		lengths = vector(length = child_n)
		radii = vector(length = child_n)
		volumes = vector(length = child_n)
		tips = vector(length = child_n)
		for(i in 1:child_n)
		{
			ind = as.numeric(children[i])
			out <- 
			asymmetric_branching(ind, branches, c_v, c_t)
			branches <- out$data
			c_v = out$c_v
			c_t = out$c_t
			tips[i] = branches[ind,]$TIPS
			volumes[i] = branches[ind,]$V_TOT
			lengths[i] = branches[ind,]$LENGTH
			radii[i] = branches[ind,]$RADIUS
		}

		#If trifurcating or worse, choose children
		#with most volume to contribute to network
		ind_1 = 1
		ind_2 = 2		
		if(child_n > 2)
		{
		ind_1 = which(max(volumes) == volumes)
		ind_2 = which(max(volumes[-ind_1]) == volumes)
		c_v = c_v + sum(volumes[-c(ind_1, ind_2)])
		}
		
		branch$TIPS = tips[ind_1] + tips[ind_2]
		branch$V_TOT = volumes[ind_1] + 
				volumes[ind_2] + 
				branch$VOLUME

		#Average scale factors
		r_c = (radii[ind_1] + radii[ind_2])/2
		l_c = (lengths[ind_1] + lengths[ind_2])/2
		
		#Difference scale factors
		d_r = (radii[ind_1] - radii[ind_2])/2
		d_l = (lengths[ind_1] - lengths[ind_2])/2

		branch$BETA = r_c/branch$RADIUS
		branch$GAMMA = l_c/branch$LENGTH

		branch$D_BETA = d_r/branch$RADIUS
		branch$D_GAMMA = d_l/branch$LENGTH
	
		if(radii[ind_1] > radii[ind_2])
		{
			if(lengths[ind_1] > lengths[ind_2])
			{branch$POS=TRUE}
			else
			{branch$POS=FALSE}
		}
		else
		{	
			if(lengths[ind_1] > lengths[ind_2])
			{branch$POS=FALSE}
			else
			{branch$POS=TRUE}
		}		
	}

	theta <- asym_scaling_exponent(branch)
	if(theta[1] <= 0 || theta [1] >= 1)
	{
		c_t = c_t + 1
		branch$WBE_THETA = NA
		branch$THETA = NA
	}
	else
	{
		branch$WBE_THETA = theta[1]
		branch$THETA = theta[2]
	} 
	
	
	branches[id,] <- branch
	return (list("data"=branches, "c_v"=c_v, "c_t"=c_t))
}

wbe_scaling_exponent <- function(branch)
{
	t_1 = (branch$BETA^2)*branch$GAMMA
	wbe = -1	
	if(t_1 > 0.5)
	{}
	else {wbe = log(t_1) / log(2)}
	
	return(wbe)	
}

asym_scaling_exponent <- function(branch)
{
	wbe = wbe_scaling_exponent(branch)
	t_1 = (2*branch$D_BETA*branch$D_GAMMA/branch$BETA*branch$GAMMA)
	t_2 = ((branch$D_BETA^2)/branch$BETA^2)
	log_factor = 1 + t_1 + t_2 
	asym = wbe + (log(log_factor)/log(2))
	return(c(-(wbe^-1), -(asym^-1)))
}

branching_analysis <- function(branches)
{
	out <- asymmetric_branching(1, branches, 0, 0)
	branches <- out$data
	print(paste("Volume pruned: ", out$c_v))
	print(paste("Pruned total volume: ", branches[1,]$V_TOT))

	print(paste("Delinquent thetas: ", out$c_t))
print(paste("Candidate thetas: ", nrow(subset(branches, TIPS>1))))
	
	#Compute whole-tree statistics	

	return (branches)
}


#data/dat_dat/Tectona_grandis_4.ply/detailed/Tectona4_2_detailed.csv
munge_oxford <- function(name, out)
{
	source("asymmetric_munge.R")
	comp_data <- read.csv(name, header=FALSE, sep="\t")
	
	data <- oxford_munge(comp_data)

	write.table(data, file=out, sep=",")
	return(data)
}

munge_computree <- function(name, out)
{
	source("asymmetric_munge.R")		
	comp_data <- read.csv(name, sep=";")
	data <- computree_munge(comp_data)
	write.table(data, file=out, sep=",")
	return(data)
}

analyse_one <- function(name)
{
	data <- read.csv(name, sep=",")
	data$CHILD_IDS <- as.character(data$CHILD_IDS)
	data <- branching_analysis(data)
	data$THETA <- pmax(pmin(data$THETA, 1), 0)	
	write.table(data, file=paste("result", name, sep="_"), sep=",")
	return(data)
}

analyse_set <- function(prefix, size)
{
	name <- paste("sample_cloud", 0, sep="")
	name <- paste(name, "_detailed.csv", sep="")
	name <- paste(prefix, name, sep="")
	print("Analysing: ")
	print(name)
	stand_data <- analyse_one(name)
	for(i in seq(1, size))
	{
		name <- paste("sample_cloud", i, sep="")
		name <- paste(name, "_detailed.csv", sep="")
		name <- paste(prefix, name, sep="")
		tree_data <- analyse_one(name)		
		stand_data <- rbind(stand_data, tree_data)
	}
	return (stand_data)
}

get_result <- function(name)
{
	val <- read.csv(paste("result_", name, sep=""))
	return (val)
}

get_results <- function(file)
{
	names <- parse(file)
	print(names)
	trees <- NULL
	ll <- list()
	for(i in seq(1:length(names)))
	{
		name <- paste(names[i])
		res <- get_result(name)
		trees <- rbind(trees, res)
		ll[[name]] <- res
	}
	return (list(frame=trees, species=ll))
}

summary_stats <- function(tree)
{
	eval <- aggregate(tree, list(order=tree$BRANCH_ORDER), function(x) mean(x,na.rm=TRUE))
	return(eval)
}

mean_val <- function(res, field)
{
	rows <- which(!is.nan(res[field]))	
	at <- mean(res[field][rows])
	return(at)
}

test_theory <- function(tree)
{
	formula <- TIPS~V_TOT
	theta <- summary_stats(tree)
	asym <- mean_val(theta, "THETA")
	wbe <- mean_val(theta, "WBE_THETA")

	d_beta <- mean_val(theta, "D_BETA")
	d_gamma <- mean_val(theta, "D_GAMMA")

	res1 <- sma(formula=formula, data=tree, na.action=na.omit,
			log="xy", slope.test=asym)
	res2 <- sma(formula=formula, data=tree, na.action=na.omit,
			log="xy", slope.test=wbe)
	return(list(asym=data.frame(res1$slopetest), 
			sym=data.frame(res2$slopetest)))
}

assemble_tests <- function(trees)
{
	names <- names(trees)
	frame <- NULL
	for(i in seq(1:length(names)))
	{
		print(names[i])
		res <- test_theory(trees[[names[i]]])
		val <- data.frame(SPECIES=names[i],
			REG_SLOPE=res$asym$b, 				ASYM_SLOPE=res$asym$test.value,
			ASYM_SIG=format(res$asym$p, digits=5),
			ASYM_DIFF=abs(res$asym$test.value-res$asym$b),
			SYM_SLOPE=res$sym$test.value, 				SYM_SIG=format(res$sym$p, digits=5),
			SYM_DIFF=abs(res$sym$test.value-res$sym$b))
		frame <- rbind(frame, val)
	}
	write.csv(frame, file="results.csv")
	return(frame)
}
