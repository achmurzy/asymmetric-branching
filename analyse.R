library(gdata)
library(reshape)
library(smatr)
library(plyr)


make_internal_frame <- function(rows)
{
frame <- data.frame(BRANCH_ORDER=integer(rows), HACK_ORDER=integer(rows), STRAHLER_ORDER=integer(rows), BIN_ORDER=integer(rows),
					BRANCH_ID=integer(rows), CHILD_IDS=character(rows), PARENT_ID=character(rows), 
					N_SYM=numeric(rows), N_ASYM=numeric(rows), N_MEAN=numeric(rows), N_MIN=numeric(rows), N_MAX=numeric(rows), 
					TIPS=integer(rows), TIPS_VOL=numeric(rows), TIPS_MEAN=numeric(rows), VOL_MEAN=numeric(rows), A_BR = numeric(rows),
                    V_TOT=numeric(rows), VOLUME = numeric(rows), L_TOT=numeric(rows), LENGTH=numeric(rows), RADIUS=numeric(rows), 
                    FIB_R=numeric(rows), FIB_L=numeric(rows), BETA=numeric(rows), GAMMA=numeric(rows), D_BETA=numeric(rows), D_GAMMA=numeric(rows), 
                    THETA=numeric(rows), WBE_THETA=numeric(rows), POS=logical(rows), INVALID=logical(rows), stringsAsFactors=FALSE)
return(frame)
}

test_tree <- function(tree="cyl_data_Ery_01.pcd_t1_m1_D0.5_DA0.075_DI0.025_L3_F5.txt")
{
  test <- open_TreeQSM(paste("data/results/",tree,sep=""))
  test <- munge_TreeQSM(test)
  test_out <- branching_analysis(test)
  return(test_out)
}

get_connected_branches <- function(branch)
{
  connect <- paste(branch$CHILD_IDS, branch$PARENT_ID, sep="_")
  connect <- trimws(unlist(strsplit(connect, split="_")))
  return(connect[which(connect != "" & connect != "0")])
}

#Recursively computes scaling relationships for every viable node in a single tree
asymmetric_branching <- function(id, branches, c_v, children)
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
		branch$TIPS_VOL = branch$VOLUME / min_element
		branch$TIPS = 1
		branch$N_SYM = 0
		branch$N_ASYM = 0
		branch$N_MEAN = branch$BIN_ORDER#0
		branch$N_MIN = branch$BIN_ORDER#0
		branch$N_MAX = branch$BIN_ORDER#0
		
		branch$V_TOT = 0 #branch$VOLUME - Exclude terminal branch volume from V_TOT to ensure empirical regressions are independent
		
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
		return (list("data"=branches, "c_v"=c_v, "children"=children))
	}
	else if (child_n == 1) #This should be considered an aberrant case as far as TreeQSM output is concerned
	{ #2/10/2019 our branch data demonstrates no cases like this in TreeQSM output
		ind = as.numeric(children[1])
		out <- asymmetric_branching(ind, branches, c_v, children)
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
		tips_vol = vector(length = child_n)
		n_mean = vector(length = child_n)
		n_min = vector(length = child_n)
		n_max = vector(length = child_n)

		#Recurse through every branch, even though we only include 2 in the scaling analysis
		#This lets us keep track of how much of the network we are forced to discard
		for(i in 1:child_n)
		{
			ind = as.numeric(children[i])
			out <- asymmetric_branching(ind, branches, c_v, c())
			branches <- out$data
			c_v = out$c_v
			#Append subtree child_ids to immediate child_ids
			children <- c(children, out$children)
			
			tips[i] = branches[ind,]$TIPS
			tips_vol[i] = branches[ind,]$TIPS_VOL
			volumes[i] = branches[ind,]$V_TOT
			lengths[i] = branches[ind,]$LENGTH
			radii[i] = branches[ind,]$RADIUS
			
			n_mean[i] = branches[ind,]$N_MEAN
			n_min[i] = branches[ind,]$N_MEAN
			n_max[i] = branches[ind,]$N_MEAN
		}

		#If trifurcating or worse, choose children
		#with most volume to contribute to network
		ind_1 = 1
		ind_2 = 2		
		if(child_n > 2)
		{
		  volsort = sort(volumes)
		  ind_1 = child_n
		  ind_2 = child_n-1
		  c_v = c_v + sum(volumes[-c(ind_1, ind_2)])
		}
		
		branch$TIPS = tips[ind_1] + tips[ind_2]
		branch$TIPS_VOL = tips_vol[ind_1] + tips_vol[ind_2]

		#Reminder that these are still 'symmetrical' in that size differences are averaged over when we iterate
		#the N counter by 1 at each bifurcation. In a symmetrical tree these would all be equivalent
		branch$N_MEAN = mean(n_mean) + 1
	  	branch$N_MIN = min(n_min) + 1
	  	branch$N_MAX = max(n_max) + 1
		
		branch$V_TOT = volumes[ind_1] + 
				volumes[ind_2] + 
				(branch$VOLUME/min_element)

		child_tips = children[which(branches[children,]$TIPS == 1)]
		
		avg_tip_vol = gm_mean(branches[child_tips,]$TIPS_VOL, na.rm=TRUE)
		#avg_tip_vol = mean(branches[child_tips,]$VOLUME, na.rm=TRUE)
		branch$VOL_MEAN = branch$V_TOT / avg_tip_vol
		branch$TIPS_MEAN = avg_tip_vol
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
	
	#Can't asymptotic formulae handle this with volume scaling? 
	branch$INVALID = branch$BETA >= 1 || branch$GAMMA >= 1
	#if(branch$INVALID)
	if(FALSE)
	{
	  branch$N_GEN <- log(branch$TIPS) / log(2)
	  branch$WBE_THETA <- NA
	  branch$THETA <- NA
	}
	else
	{
	  #Node-level finite scaling analyses - Not using Bentley et al 2013 node level metrics
	  #This formula for handling low - N trees may not be really good for aggregating scaling ratios
	  fsa <- finite_size_analysis(branch)
	  branch$WBE_THETA <- fsa$wbe
	  branch$THETA <- fsa$asym
	  branch$N_SYM = fsa$N
	  
	  branches[id,] <- branch
	  #Calculate the asymmetric branching ratio on logarithmically binned branch ratios
	  if(length(children) > 1)
	  {
	  	#We might try doing this with other branch ordering schemes

		#This step is not very well behaved - slopes close to 1 or diverging misbehave	  	
	  	#Estimating branch ratios this way doesn't work because smallest branches are underrepresented
	  	#branch_dist <- count(branches[as.integer(c(id,children)),]$BIN_ORDER)
	  	#Its the same problem as using tips to estimate the MSE
	  	#slope = lm(log(freq)~x, branch_dist)
	  	#a_b_r = exp(abs(slope$coefficients[[2]]))

	  	#Use the same solution - use tip volumes instead 
	  	#Won't this just reverse the logarithmic binning and give approx lambda (2)?
		valid_subtree <- c(id,children)[which(branches[as.integer(c(id,children)),]$VOLUME > 0)]
		slope = lm(log(TIPS_VOL)~BIN_ORDER, branches[valid_subtree,])
	  	a_b_r = exp(abs(slope$coefficients[[2]]))
	  	#if(is.na(a_b_r)) Sometimes parents/children occupy the same bin
	  	#{
	  	#	print(branch)
	  	#	print(slope)
	  	#	print(branch_dist)
	  	#}
	  	branch$N_ASYM = log(branch$TIPS) / log(a_b_r)
	  	branch$A_BR = a_b_r
	  }
	  else
	  {
	  	branch$N_ASYM = NA	
	  }
	}
	
	#branch$WBE_THETA <- nodal_wbe(2, branch$BETA, branch$GAMMA)
	#branch$THETA <- nodal_asym(2, branch$BETA, branch$GAMMA, branch$D_BETA, branch$D_GAMMA)
	#branch$INVALID = branch$WBE_THETA <= 0 || branch$WBE_THETA >= 1 || branch$THETA <= 0 || branch$THETA >= 1
	branches[id,] <- branch
	return (list("data"=branches, "c_v"=c_v, "children"=children))
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
  N = log(subtree$TIPS)/log(2)  
  sym_vol = 2*(subtree$BETA^2)*subtree$GAMMA
  asym_vol = sym_vol + 4*subtree$BETA*subtree$D_BETA*subtree$D_GAMMA + 2*subtree$GAMMA*(subtree$D_BETA^2)
  wbe = asymptotic_formula(sym_vol, N)
  asym = asymptotic_formula(asym_vol, N)
  return(list("wbe"=wbe, "asym"=asym, "N"=N))
}

#Remember - every 2 here is a branching ratio that might be replaced
#Every 2^N is the number of tips in our networks - network size, and
#we're throwing away data on variation in volume. None of this theory
#is well-suited to empirical parameterization - that's your job.
asymptotic_formula <- function(volumetric_scaling, N)
{
  if(volumetric_scaling <= 0.99) #Large network expression - heuristic threshold to protect against asymptotic behavior
  {
    res = log(2^N) / (log(2^N) + 
                        log(1 - (volumetric_scaling^(N+1))) - 
                        log((1-volumetric_scaling)*(volumetric_scaling^(N))))
  }
  else #Functions under asymptotic behavior?
  {
    res = log(2^N) / ((log((N+1)*(2^N))) - N*log(volumetric_scaling))
  }
  return(res)
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#Any tree-level parameter we need to compute from branches
whole_tree_scaling <- function(branches, verbose_report = TRUE)
{
  #Make subsets of the branch network for various analyses...
  tips <- which(branches$TIPS == 1)
  non_tips <- which(branches$TIPS > 1)
  #Invalid nodes have scaling ratios outside of theoretical bounds (child dimensions greater than parents)
  valid_nodes <- which(!branches$INVALID)

  #This regression permits an estimate of branching ratio
  #if we assume that branch order is proportional to total subtree volume
  empirical_v = sma(formula=TIPS~VOLUME, data=branches[non_tips,], log='xy', method="SMA")
  empirical_v = empirical_v$coef[[1]][2,1]

  empirical_o = sma(formula=TIPS~BRANCH_ORDER, data=branches[non_tips,], log='xy', method="SMA")
  empirical_o = empirical_o$coef[[1]][2,1]

  #Using Horsfield 1976 to compute branching ratio
  #He uses different definitions from WBE, so might be crossing our wires
  #Have the feeling his is the correct one: 
  #"In an asymmetrical tree... the number of divisions down from the stem
  #is not the same to every end branch". His whole-tree branching ratio
  #based on regression is more in the spirit of what we're doing.
  #That is, branching ratio (n) and network depth (N) are convenient 
  #mathematical scaffolding for working with 'averaged networks', but with
  #datasets like these, we are losing a lot of information by working with
  #'averaged networks'. Surprisingly, asymmetry is damaging these concepts 
  #more than the actual geometry of the network elements
  #Completely reformulating theory on the basis of volume as THE touchstone
  #of size within the network might allow linking to data more easily.

  #Once we do this, might actually be able to compute delta using the
  #branching ratio to link volume and branch order.
  num_tips = branches[1,]$TIPS #Could replace this with total tip volume too
  #Doing this increases estimates of scaling, period.
  #branching_ratio = exp(empirical) 
  #Alternatives:
  #branching_ratio = exp(empirical_v)
  #branching_ratio = exp(empirical_o)
  branching_ratio = 2

  #The only reason we're doing this is because incorrectly calculating N
  #may be seriously affecting our estimates of scaling exponents
  N = log(num_tips)/log(branching_ratio)  
  N_SYM = branches[1,]$N_SYM
  N_ASYM = branches[1,]$N_ASYM
  N_MEAN = branches[1,]$N_MEAN
  N_MIN = branches[1,]$N_MIN
  N_MAX = branches[1,]$N_MAX
  
  #Path fraction (avg L / max L) - Excluding nodes could be affecting this
  avg_L <- mean(branches[tips,]$L_TOT, na.rm = TRUE)
  max_L <- max(branches[tips,]$L_TOT, na.rm = TRUE)
  path_frac = avg_L/max_L
  
  #Asymmetry fraction (pos_asym / nodes) - since we include even volume-pruned sections, use length(tips)
  asym_frac = length(which(branches[non_tips,]$POS)) / length(non_tips)
  
  #Fibbonaci ratios - might tell us something else about asymmetry types (related to d_beta, d_gamma)?
  fib_r = mean(branches[non_tips,]$FIB_R, na.rm = TRUE)
  fib_l = mean(branches[non_tips,]$FIB_L, na.rm = TRUE)

  beta = mean(branches[valid_nodes,]$BETA, na.rm = TRUE)
  gamma = mean(branches[valid_nodes,]$GAMMA, na.rm = TRUE)

  #No straightforward way to average D_BETA - but arithmetic is fine because radial asymmetry is usually small
  d_beta = mean(abs(branches[valid_nodes,]$D_BETA), na.rm = TRUE)
  #Geometric mean because we ensure our length asymmetries are always positive by convention
  d_gamma = mean(branches[valid_nodes,]$D_GAMMA, na.rm = TRUE)
  
  #Volume scaling expressions - ratio of total volume from sibling branches to parent branches
  sym_vol = 2*(beta^2)*gamma
  node_wbe = mean(branches[valid_nodes,]$WBE_THETA, na.rm=TRUE)
  
  #This is the mathematical reason asymmetric estimates increase scaling
  #But using agregated d_beta is definitely a source of error - biases our results upward
  asym_vol = 4*beta*d_beta*d_gamma + 2*gamma*(d_beta^2)
  vol_scaling = sym_vol + asym_vol
  node_asym = mean(branches[valid_nodes,]$THETA, na.rm=TRUE)
  
  wbe = asymptotic_formula(sym_vol, N)
  asym = asymptotic_formula(vol_scaling, N)
  
  #finite_size_correction = (1 / (1 + (1/4)*(branching_ratio*num_tips)^(-1/3)))  
  #Standardize data for regressions - normalize and scale so that we can compare everything
  branches[non_tips,]$V_TOT <- branches[non_tips,]$V_TOT / min(branches[non_tips,]$V_TOT)
  branches[non_tips,]$VOL_MEAN <- branches[non_tips,]$VOL_MEAN / min(branches[non_tips,]$VOL_MEAN)
  branches[non_tips,]$TIPS_VOL <- branches[non_tips,]$TIPS_VOL / min(branches[non_tips,]$TIPS_VOL)
  branches[non_tips,]$TIPS <- (branches[non_tips,]$TIPS / 2)
  #TIPS on the same scale as tip volume for comparng regressions - scaling the axes won't matter for log-log regressions - right?
  #tip_vol_factor = (max(branches[non_tips,]$TIPS_VOL)/max(branches[non_tips,]$TIPS))
  #branches[non_tips,]$TIPS <- branches[non_tips,]$TIPS * tip_vol_factor

  empirical = sma(formula=log(TIPS)~log(V_TOT), data=branches[non_tips,], method="SMA")
  ci = empirical$slopetest[[1]]$ci
  branches$PRED <- NA
  branches[non_tips,]$PRED = (empirical$coef[[1]][1][1,]) + (log(branches[non_tips,]$V_TOT) * empirical$coef[[1]][1][2,])
  empirical = empirical$coef[[1]][2,1]

  #Don't fit terminal tips for empirical regressions so that data is 'independent'
  empirical_vol = sma(formula=log(TIPS_VOL)~log(V_TOT), data=branches[non_tips,], method="SMA")
  ci_vol = empirical_vol$slopetest[[1]]$ci
  branches$PRED_VOL <- NA
  branches[non_tips,]$PRED_VOL = empirical_vol$coef[[1]][1][1,] + (log(branches[non_tips,]$V_TOT) * empirical_vol$coef[[1]][1][2,])
  empirical_vol = empirical_vol$coef[[1]][2,1]  

  empirical_vol_mean= sma(formula=log(TIPS)~log(VOL_MEAN), data=branches[non_tips,], method="SMA")
  ci_ma = empirical_vol_mean$slopetest[[1]]$ci
  branches$PRED_MA <- NA
  branches[non_tips,]$PRED_MA = empirical_vol_mean$coef[[1]][1][1,] + (log(branches[non_tips,]$VOL_MEAN) * empirical_vol_mean$coef[[1]][1][2,])
  empirical_vol_mean = empirical_vol_mean$coef[[1]][2,1]

  empirical_vol_tip = sma(formula=log(TIPS_MEAN*TIPS)~log(V_TOT), data=branches[non_tips,], method="SMA")
  branches$PRED_TIP <- NA
  branches[non_tips,]$PRED_TIP = empirical_vol_tip$coef[[1]][1][1,] + (log(branches[non_tips,]$V_TOT) * empirical_vol_tip$coef[[1]][1][2,])
  empirical_vol_tip = empirical_vol_tip$coef[[1]][2,1]

  empirical_tip_mean = sma(formula=log(TIPS_VOL/TIPS_MEAN)~log(V_TOT), data=branches[non_tips,], method="SMA")
  branches$PRED_MEAN <- NA
  branches[non_tips,]$PRED_MEAN = empirical_tip_mean$coef[[1]][1][1,] + (log(branches[non_tips,]$V_TOT) * empirical_tip_mean$coef[[1]][1][2,])
  empirical_tip_mean = empirical_tip_mean$coef[[1]][2,1]  

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
              "wbe"=wbe, "node_wbe"=node_wbe, "asym"=asym, "node_asym"=node_asym, 
			  "sym_vol" = sym_vol, "asym_vol" = asym_vol, "network_n"=N_SYM, "n_asym"=N_ASYM,
			  "n_mean"=N_MEAN, "n_min"=N_MIN, "n_max"=N_MAX,
			  "empirical"=empirical, "empirical_vol"=empirical_vol, "empirical_vol_mean"=empirical_vol_mean,
			  "empirical_vol_tip"=empirical_vol_tip, "empirical_tip_mean"=empirical_tip_mean,
			  "ci_min"=ci[1], "ci_max"=ci[2], "ci_vol_min"=ci_vol[1], "ci_vol_max"=ci_vol[2],
			  "ci_ma_min"=ci_ma[1], "ci_ma_max"=ci_ma[2], "path_frac"=path_frac, "asym_frac"=asym_frac, 
			  "preds_tips"=branches))
}

## ENTRY POINT ##
#Entry point for the recursive analysis where we can intialize what we need to, and report results
branching_analysis <- function(branches, verbose_report = TRUE)
{
  	branches[1,]$L_TOT = branches[1,]$LENGTH  #Lets us additively compute path length for every node
	
	non_zero <- which(branches$VOLUME > 0)
	min_element <<- min(branches[non_zero,]$VOLUME)
	print(min_element)

	#Initiate recursion to process the network
	out <- asymmetric_branching(1, branches, 0, c())
	branches <- out$data
	
	#Compute whole-tree scaling ratios - compare to average scaling ratios from nodes	
	scaling = whole_tree_scaling(branches, verbose_report)
	branches <- scaling$preds_tips
  	scaling <- scaling[-length(scaling)] #Get rid of last element there
  	
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
