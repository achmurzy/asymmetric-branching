source("asymmetric_analysis.R")

library(hash)
library(plyr)
#We are representing trees as data frames
#Framework for D0L-Systems
#May need to implement recursive derivation functions to more easily compute whole-tree statistics
#Adapt from your previous recursive function?
#All the same functions for L-Systems, QSMs, etc using a common tree frame structure.
#Recursively write down derivation rules according to integer N, depth-first
#Perhaps less efficient than massively parallel application of rules? Also lose-whole organism structure in ontogeny (emphasis on final result)

#Parallel application can still be recursive?
#Only need to replace for-loops below
#Previous framework was analysing branches based on row-id vectors

#Creates an axiom from a given symbol
#Parallel re-writing requires 'do not derive' flag on new productions, so we make this default
create_axiom <- function(symbol)
{
	axiom <- make_internal_frame(1)
	axiom$BRANCH_ID = 0
	axiom$MODULE = symbol
	axiom$WRITE = FALSE	
	axiom$LENGTH = 1
	axiom$RADIUS = 0.1
	return(axiom)
} 

create_production <- function(LHS, RHS, prodHash)
{
	left <- as.character(LHS)
	RHS <- unlist(strsplit(RHS, split=""))
	right <- make_internal_frame(length(RHS))
	right$MODULE <- RHS
	right$RADIUS = rep(2^(-1/2), nrow(right))
	right$LENGTH = rep(2^(-1/3), nrow(right))
	
	prodHash[[left]] <- right
	return(prodHash)
}

create_system <- function(alphabet, ax, prod)
{
	axiom <- create_axiom(ax)
	axiom$WRITE = TRUE
	productions <- hash()
	for(i in seq(1, length(alphabet)))
	{create_production(alphabet[i], prod[i], productions)}
	return(list(AXIOM=axiom, PRODUCTIONS=productions))
}

#recursive workhorse for l-system derivation
apply_production <- function(LHS, sys, derivation)
{
	mods <- which(derivation$MODULE == LHS & derivation$WRITE)
	if(length(mods) > 0)
	{
		ancestors <- derivation[mods,]
		descendants <- adply(ancestors, 1, function(x){production_descent(x, sys$PRODUCTIONS)})
		derivation <- derivation[-mods,]
		derivation <- rbind(derivation, descendants)
	}
	return(derivation)
}

#Performs direct extraction of production rows and transfer of values from parent module
production_descent <- function(row, prods)
{
	children <- prods[[as.character(row$MODULE)]]
	children <- data.frame(children)
	#Generate child Ids
	#Add child Ids to parent
	#Add parent Id to children
	#Scale length and radius by production values
	children$BRANCH_ORDER = rep(row$BRANCH_ORDER+1, nrow(children))
	children$LENGTH = row$LENGTH * children$LENGTH
	children$RADIUS = row$RADIUS * children$RADIUS
	row$V_TOT = sum(pi * children$RADIUS * children$RADIUS * children$LENGTH)
	return(children)
}

#Perform a single parallel rewrite
perform_derivation <- function(sys, derivation)
{
	LHS <- keys(sys$PRODUCTIONS)
	for(i in seq(1, length(sys$PRODUCTIONS)))
	{derivation <- apply_production(LHS[i], sys, derivation)}
	
	#Only write non-terminals: A, B, C
	nonterminals <- which(derivation$MODULE == 'A')
	derivation[nonterminals,]$WRITE <- rep(TRUE, length(nonterminals))	
	return(derivation)
}

#Takes an integer gens and an L-System 2-list c(axiom, productions)
#Returns a derivation tree
derive_generations <- function(gens, sys, growth=FALSE)
{
	if(growth)
	{
		LHS <- keys(sys$PRODUCTIONS)
		growthMatrix = matrix(nrow=gens, ncol=length(LHS)+1)
		colnames(growthMatrix) <- c("Generation", LHS)
	}
	derive = sys$AXIOM
	for(gen in seq(1, gens))
	{
		derive <- perform_derivation(sys, derive)
		if(growth)
		{
			growthMatrix[gen,1] = gen
			for(i in seq(2, length(LHS)+1))
			{
				growthMatrix[gen, i] = length(which(derive$MODULE == LHS[i-1]))
			}
		}		
	}
	if(growth)
		return(list(GROWTH=growthMatrix, STRING=derive))
	else
		return(derive)
}

derive_recursively <- function(generations, system)
{
	derive <- system$AXIOM
	id <- system$AXIOM$BRANCH_ID
	out <- asymmetric_derivation(id, derive, system, generations, 0, 0, 0)
	branches <- out$data
	return (branches)
}

#Recursively derive an L-System to easily compute sub-tree statistics and aggregate across the whole structure
#Takes the branch id under analysis, the total set of tree branches, id counter, and cumulative geometric values (volume, petioles) on a recursive internode

#L-Systems are replacement grammars - must replace the left-hand side or nothing makes sense
#Use the internode 'I' symbol to compute geometry, non-terminals are meristems, model things this way.
#The identity of modules comes from the meristem that produced them. Therefore all branches can be the same
#geometric 'I' symbol, with scaling parameters produced by their parent meristem.

#Unfortunately, separating meristematic non-terminals from geometric terminals greatly complicates 
#recursive computation of tree geometry, maybe
asymmetric_derivation <- function(id, branches, system, gen, c_c, c_i)
{
	row_id <- which(branches$BRANCH_ID==id)
	branch <- branches[row_id,]

	c_c = c_c+1	#Iterate the child counter for a unique ID at this meristematic derivation
	
	children <- system$PRODUCTIONS[[branch$MODULE]] #Apply production from meristem immediately
	children$PARENT_ID = branch$BRANCH_ID #The meristems should have the same parent ID as the local internode
	children$BRANCH_ID = c_c	#meristems share ID with internode - establish parent-child relations		

	branches <- rbind(branches, children) #Add all the children, including future meristems
	child_ids <- which(branches$PARENT_ID == branch$BRANCH_ID) #Save their row names 
	internode <- which(branches[child_ids,]$MODULE == 'I')

	branches[child_ids,]$BRANCH_ORDER = branch$BRANCH_ORDER+1
	branches[internode,]$LENGTH = branch$LENGTH * children$LENGTH
	branches[internode,]$RADIUS = branch$RADIUS * children$RADIUS
	if(length(internode) > 1)
	{
		print("More than one child internode, check L-System")
	}
	
	child_ids <- child_ids[! child_ids %in% internode]
	child_n = length(child_ids)		

	if(gen == 0) #Our current internode is a terminal element and we'll stop the derivation here
	{
		branch <- branches[internode,]
		branch$TIPS = 1
		branch$V_TOT = branch$VOLUME
		branch$BETA = NA	
		branch$GAMMA = NA
		branch$D_BETA = NA
		branch$D_GAMMA = NA
		branch$WBE_THETA = NA
		branch$THETA = NA
		branches[internode,] <- branch
		#Also delete meristems here?
		return (list("data"=branches, "c_c"=c_c, "c_i"=internode))
	}
	else
	{
		lengths = vector(length = child_n)
		radii = vector(length = child_n)
		volumes = vector(length = child_n)
		tips = vector(length = child_n)
		ids = vector(length = child_n)
		for(i in 1:child_n)					#Recursively compute subtree statistics
		{
			child <- branches[child_ids[i],]
			out <- asymmetric_derivation(child$BRANCH_ID, branches, system, gen-1, c_c, c_v, c_t)
			branches <- out$data
			c_c = out$c_c
			c_i = out$c_i	
			
			#c_i contains the most recently derived internode
			if(c_i == NULL)
			{
				ids[i] = ""
				tips[i] = branches[c_i,]$TIPS
				volumes[i] = branches[c_i,]$V_TOT
				lengths[i] = branches[c_i,]$LENGTH
				radii[i] = branches[c_i,]$RADIUS	
			}
			else
			{
				ids[i] = branches[c_i,]$BRANCH_ID
				tips[i] = branches[c_i,]$TIPS
				volumes[i] = branches[c_i,]$V_TOT
				lengths[i] = branches[c_i,]$LENGTH
				radii[i] = branches[c_i,]$RADIUS
			}
		}

		#Remove meristems from derivation when??

		branches[internode,]$TIPS = sum(tips)
		branches[internode,]$V_TOT = sum(volumes) +  branch$VOLUME
		branches[internode,]$L_TOT = sum(lengths) + branch$LENGTH
		branches[internode,]$CHILD_IDS = paste(ids, collapse='_')

		#Average scale factors
		r_c = mean(radii)
		l_c = mean(lengths)

		#Difference scale factors - no formulae for n > 2
		#d_r = (radii[ind_1] - radii[ind_2])/2
		#d_l = (lengths[ind_1] - lengths[ind_2])/2

		branch$BETA = r_c/branch$RADIUS
		branch$GAMMA = l_c/branch$LENGTH

		#Tokunaga, Turcotte 1998
		#branch$D_BETA = d_r/branch$RADIUS
		#branch$D_GAMMA = d_l/branch$LENGTH
		#if(radii[ind_1] > radii[ind_2])
		#{
		#	if(lengths[ind_1] > lengths[ind_2])
		#	{branch$POS=TRUE}
		#	else
		#	{branch$POS=FALSE}
		#}
		#else
		#{	
		#	if(lengths[ind_1] > lengths[ind_2])
		#	{branch$POS=FALSE}
		#	else
		#	{branch$POS=TRUE}
		#}
		#theta <- asym_scaling_exponent(branch)
		#if(theta[1] <= 0 || theta [1] >= 1)
		#{
		#	c_t = c_t + 1
		#	branch$WBE_THETA = NA
		#	branch$THETA = NA
		#}
		#else
		#{
		#	branch$WBE_THETA = theta[1]
		#	branch$THETA = theta[2]
		#}
	}
	
	branches[internode,] <- branch
	return (list("data"=branches, "c_c"=c_c, "c_i"=internode))		
}