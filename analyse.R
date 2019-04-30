library(gdata)
library(reshape)
library(smatr)


make_internal_frame <- function(rows)
{
frame <- data.frame(BRANCH_ORDER=integer(rows), HACK_ORDER=integer(rows), STRAHLER_ORDER=integer(rows), N_GEN=numeric(rows),
                    BRANCH_ID=integer(rows), CHILD_IDS=character(rows), PARENT_ID=character(rows), 
                    TIPS=integer(rows), V_TOT=numeric(rows), VOLUME = numeric(rows), L_TOT=numeric(rows), LENGTH=numeric(rows), 
                    RADIUS=numeric(rows), FIB_R=numeric(rows), FIB_L=numeric(rows), BETA=numeric(rows), GAMMA=numeric(rows), D_BETA=numeric(rows), D_GAMMA=numeric(rows), 
                    THETA=numeric(rows), WBE_THETA=numeric(rows), POS=logical(rows), INVALID=logical(rows), stringsAsFactors=FALSE)
return(frame)
}

get_connected_branches <- function(branch)
{
  connect <- paste(branch$CHILD_IDS, branch$PARENT_ID, sep="_")
  connect <- trimws(unlist(strsplit(connect, split="_")))
  return(connect[which(connect != "" & connect != "0")])
}

#Recursively computes scaling relationships for every viable node
#Computed nodal values:
#BETA - average radius scale factor 
#GAMMA - average length scale factor
#D_BETA - difference radius scale factor
#D_GAMMA - difference length scale factor
#WBE_THETA - symmetric scaling exponent from average scale factor
#THETA - asymmetric scaling exponent from difference and average scale factors
asymmetric_branching <- function(id, branches, c_v)
{
	branch <- branches[id,]
	children <- unlist(strsplit(branch$CHILD_IDS[1], split="_"))
	child_n = length(children)
	
	#Get path length before recursions terminate
	p_id <- branch$PARENT_ID
	parent_l = 0
	if(p_id != "")
	  parent_l = branches[p_id,]$L_TOT  
	branch$L_TOT <- parent_l + branch$LENGTH
	
	if(child_n == 0)
	{
		branch$TIPS = 1
		branch$N_GEN = 0
		branch$V_TOT = branch$VOLUME
		branch$FIB_R = NA
		branch$FIB_L = NA
		branch$BETA = NA	
		branch$GAMMA = NA
		branch$D_BETA = NA
		branch$D_GAMMA = NA
		branch$WBE_THETA = NA
		branch$THETA = NA
		branch$POS = NA
		
		branches[id,] <- branch
		return (list("data"=branches, "c_v"=c_v))
	}
	else if (child_n == 1) #This should be considered an aberrant case as far as TreeQSM output is concerned
	{ #2/10/2019 our branch data demonstrates no cases like this in TreeQSM output
		ind = as.numeric(children[1])
		out <- asymmetric_branching(ind, branches, c_v)
		branches = out$data
		c_v = out$c_v

		branch$TIPS = branches[ind,]$TIPS
		branch$V_TOT = branches[ind,]$V_TOT + branch$VOLUME

		branch$BETA = branches[ind,]$RADIUS/branch$RADIUS
		branch$GAMMA = branches[ind,]$LENGTH/branch$LENGTH
		
		branch$D_BETA = 0
		branch$D_GAMMA = 0
		branch$FIB_R = 1
		branch$FIB_L = 1
		branch$POS = NA
	}
	else
	{
		tips = 0
		vol = 0
		lengths = vector(length = child_n)
		radii = vector(length = child_n)
		volumes = vector(length = child_n)
		tips = vector(length = child_n)
		
		#Recurse through every branch, even though we only include 2 in the scaling analysis
		for(i in 1:child_n)
		{
			ind = as.numeric(children[i])
			out <- asymmetric_branching(ind, branches, c_v)
			branches <- out$data
			c_v = out$c_v
			
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
		
		if(radii[ind_1] > radii[ind_2])
		{
		  branch$FIB_R = radii[ind_2] / radii[ind_1]
		  if(lengths[ind_1] > lengths[ind_2])
		  {
		    branch$FIB_L = lengths[ind_2] / lengths[ind_1]
		    branch$POS=TRUE
		  }
		  else #swap indices - length will be negative
		  {
		    tmp = ind_2
		    ind_2 = ind_1
		    ind_1 = tmp
		    branch$FIB_L = lengths[ind_2] / lengths[ind_1]
		    branch$POS=FALSE
		  }
		}
		else
		{	
		  branch$FIB_R = radii[ind_1] / radii[ind_2]
		  if(lengths[ind_1] > lengths[ind_2])
		  {
		    branch$FIB_L = lengths[ind_2] / lengths[ind_1]
		    branch$POS=FALSE
		  }
		  else #swap indices - bigger vessel is ind_2
		  {
		    tmp = ind_2
		    ind_2 = ind_1
		    ind_1 = tmp
		    branch$FIB_L = lengths[ind_2] / lengths[ind_1]
		    branch$POS=TRUE
		  }
		}		

		#Average scale factors
		r_c = (radii[ind_1] + radii[ind_2])/2
		l_c = (lengths[ind_1] + lengths[ind_2])/2
		
		#Difference scale factors
		#if the node shows negative asymmetry, you'll have a negative difference factor
		#Ensure length-difference factor is always positive by swapping indices above
		d_r = (radii[ind_1] - radii[ind_2])/2
		d_l = (lengths[ind_1] - lengths[ind_2])/2

		branch$BETA = r_c/branch$RADIUS
		branch$GAMMA = l_c/branch$LENGTH

		branch$D_BETA = d_r/branch$RADIUS
		branch$D_GAMMA = d_l/branch$LENGTH
	}
	
	branch$INVALID = branch$BETA >= 1 || branch$GAMMA >= 1
	if(branch$INVALID)
	{
	  branch$N_GEN <- log(branch$TIPS) / log(2)
	  branch$WBE_THETA <- NA
	  branch$THETA <- NA
	}
	else
	{
	  #Node-level finite scaling analyses
	  fsa <- finite_size_analysis(branch)
	  branch$WBE_THETA <- fsa$wbe
	  branch$THETA <- fsa$asym
	  branch$N_GEN <- fsa$N
	}
	
	#branch$WBE_THETA <- nodal_wbe(2, branch$BETA, branch$GAMMA)
	#branch$THETA <- nodal_asym(2, branch$BETA, branch$GAMMA, branch$D_BETA, branch$D_GAMMA)
	#branch$INVALID = branch$WBE_THETA <= 0 || branch$WBE_THETA >= 1 || branch$THETA <= 0 || branch$THETA >= 1
	branches[id,] <- branch
	return (list("data"=branches, "c_v"=c_v))
}

##These nodal equations are derived from infinite-network assumptions and so display problematic asymptotic behavior
##Use at your own risk

##Node-based theta from Bentley et al 2013. Robust to higher-order branching (trifurcations) if needed
##Assumes symmetry by using 'average scale factors', in Brummer's terminology
nodal_wbe <- function(n, beta, gamma)
{
  a = -log(beta) / log(n)
  b = -log(gamma) / log(n)
  wbe = 1 / (2*a + b)
  return(wbe)	
}

##Using Brummer's additive form to segment symmetric (average) and asymmetric (difference) contributions to scaling
##Only applicable to bifurcations
nodal_asym <- function(n, beta, gamma, d_beta, d_gamma)
{
  sym = log(gamma*(beta^2))/log(2)
  t_1 = (2*d_beta*d_gamma/(beta*gamma))
  t_2 = ((d_beta^2)/beta^2)
  asym = (1 + t_1 + t_2)/log(2) 
  theta = -(sym + asym)^(-1)
  return(theta)
}

finite_size_analysis <- function(subtree)
{
  N = log(subtree$TIPS)/log(2)  #Computing N this way... is it circular? Also assumes bifurcation at equal depths just because
  sym_vol = 2*(subtree$BETA^2)*subtree$GAMMA
  asym_vol = sym_vol + 4*subtree$BETA*subtree$D_BETA*subtree$D_GAMMA + 2*subtree$GAMMA*(subtree$D_BETA^2)
  wbe = asymptotic_formula(sym_vol, N)
  asym = asymptotic_formula(asym_vol, N)
  return(list("wbe"=wbe, "asym"=asym, "N"=N))
}

asymptotic_formula <- function(volumetric_scaling, N)
{
  if(abs(1-volumetric_scaling)>0.1) #Large network expression - heuristic threshold to protect against asymptotic behavior
  {
    res = log(2^N) / (log(2^N) + 
                        log(1 - volumetric_scaling^(N+1)) - #Generates NaN - volumetric^(N+1) > 1 when sym_vol > 1
                        log((1-volumetric_scaling)*volumetric_scaling^(N))) #Generates NaN - sym_vol > 1
  }
  else
  {
    res = log(2^N) / ((log(N+1)*(2^N)) - N*log(volumetric_scaling))
  }
  return(res)
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#Any tree-level parameter we need to compute from branches
whole_tree_scaling <- function(branches, verbose_report = TRUE)
{
  #Compute branching generations N (network depth) assuming bifurcation
  tips <- which(branches$TIPS == 1)
  
  #IMPORTANT !!!
  #The below two are NOT the same because you're pruning volume from the tree !!!
  #num_tips = length(tips)   <-- And this one really messes up our estimates of scaling !!! Because N is so wrong !!!
  num_tips = branches[1,]$TIPS
  N = log(num_tips)/log(2)  #Computing N this way... is it circular? Also assumes bifurcation at equal depths just because
  
  branches <- branches[which(!branches$INVALID),]
  
  beta = mean(branches$BETA, na.rm = TRUE)
  gamma = mean(branches$GAMMA, na.rm = TRUE)
  
  #Path fraction (avg L / max L) - Excluding nodes could be affecting this
  avg_L <- mean(branches[tips,]$L_TOT, na.rm = TRUE)
  max_L <- max(branches[tips,]$L_TOT, na.rm = TRUE)
  path_frac = avg_L/max_L
  
  #Asymmetry fraction (pos_asym / nodes) - since we include even volume-pruned sections, use length(tips)
  asym_frac = length(which(branches[-tips,]$POS)) / (nrow(branches) - length(tips))
  
  #Fibbonaci ratios - might tell us something else about asymmetry types (related to d_beta, d_gamma)?
  fib_r = mean(branches[-tips,]$FIB_R, na.rm = TRUE)
  fib_l = mean(branches[-tips,]$FIB_L, na.rm = TRUE)

  #No straightforward way to average D_BETA - but arithmetic is fine because radial asymmetry is usually small
  d_beta = mean(abs(branches$D_BETA), na.rm = TRUE)
  #Geometric mean because we ensure our length asymmetries are always positive by convention
  d_gamma = mean(branches$D_GAMMA, na.rm = TRUE)
  
  #Volume scaling expressions - ratio of total volume from sibling branches to parent branches
  sym_vol = 2*(beta^2)*gamma
  node_wbe = mean(branches$WBE_THETA, na.rm=TRUE)
  
  #Distribution of volume factors instead of dealing with 

  asym_vol = sym_vol + 4*beta*d_beta*d_gamma + 2*gamma*(d_beta^2)
  node_asym = mean(branches$THETA, na.rm=TRUE)
  
  wbe = asymptotic_formula(sym_vol, N)
  asym = asymptotic_formula(asym_vol, N)
  
  empirical = sma(formula=TIPS~V_TOT, data=branches, log='xy', method="MA")
  ci = empirical$slopetest[[1]]$ci
  empirical = empirical$coef[[1]][2,1]
  
  if(verbose_report)
  {
    print("Finite network correction assuming bifurcation...")
    print(paste("Inferred network generations:", N))
    print(paste("1-labeled terminal elements", length(tips)))
    print(paste("Trunk-distal tips", num_tips))
    print(paste("Tip discrepancy", length(tips) - branches[1,]$TIPS))
    
    print(paste("Only using valid nodes for whole-tree scaling:", nrow(branches)))
    
    print(paste("Average beta:", beta))
    print(paste("Average gamma:", gamma))
    print(paste("Average d_beta:", d_beta))
    print(paste("Average d_gamma:", d_gamma))
    
    print(paste("Symmetric volume scaling:", sym_vol))
    print(paste("Asymmetric volume scaling:", asym_vol))
  }
  #This should always be returning simple scalar summary stats
  return(list("beta"=beta, "gamma"=gamma, "d_beta"=d_beta, "d_gamma"=d_gamma, "fib_r" = fib_r, "fib_l"=fib_l,  
              "wbe"=wbe, "node_wbe"=node_wbe, "asym"=asym, "node_asym"=node_asym, "empirical"=empirical, "network_n"=N,
              "ci_min"=ci[1], "ci_max"=ci[2], "path_frac"=path_frac, "asym_frac"=asym_frac))
}

#Entry point for the recursive analysis where we can intialize what we need to, and report results
branching_analysis <- function(branches, verbose_report = TRUE)
{
  branches[1,]$L_TOT = branches[1,]$LENGTH  #Lets us additively compute path length for every node
  
	out <- asymmetric_branching(1, branches, 0)
	branches <- out$data
	
	#Compute whole-tree scaling ratios - compare to average scaling ratios from nodes	
	scaling = whole_tree_scaling(branches, verbose_report)
  
  if(verbose_report)
  {
    print(paste("Volume pruned: ", out$c_v))
    print(paste("Total volume after pruning: ", branches[1,]$V_TOT))
    
    print(paste("Delinquent thetas: ", length(which(branches$INVALID))))
    print(paste("Candidate thetas: ", nrow(subset(branches, TIPS>1)))) 
    
    print(scaling)
  }
  
	return (list("branches"=branches, "scaling"=scaling))
}
