library(gdata)
library(reshape)
library(smatr)
library(plyr)


make_internal_frame <- function(rows)
{
frame <- data.frame(BRANCH_ORDER=integer(rows), HACK_ORDER=integer(rows), STRAHLER_ORDER=integer(rows), BIN_ORDER=integer(rows),
					BRANCH_ID=integer(rows), CHILD_IDS=character(rows), PARENT_ID=character(rows), 
					N_SYM=numeric(rows), N_ASYM=numeric(rows), N_MEAN=numeric(rows), N_MIN=numeric(rows), N_MAX=numeric(rows), 
					TIPS=integer(rows), TIPS_VOL=numeric(rows), RAW_TIPS_VOL=numeric(rows), RAW_V_TOT=numeric(rows), 
					FULL_TIPS=numeric(rows), FULL_TIPS_VOL=numeric(rows), FULL_V_TOT=numeric(rows), 
					TIPS_MEAN=numeric(rows), TIPS_VAR=numeric(rows), TIPS_CV=numeric(rows), TIPS_GEO_CV=numeric(rows),
					RADS_MEAN=numeric(rows), VOL_MEAN=numeric(rows), A_BR = numeric(rows),
                    CHILD_LENGTHS=numeric(rows), CHILD_RADII=numeric(rows), CHILD_VOLS=numeric(rows), 
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
	#p_id <- branch$PARENT_ID
	#parent_l = 0
	#if(p_id != "")
	#  parent_l = branches[p_id,]$L_TOT  
	#branch$L_TOT <- parent_l + branch$LENGTH
	
	if(child_n == 0 | branch$INVALID)
	{
		branch$TIPS = 1
		branch$FULL_TIPS = 1

		branch$TIPS_VOL = branch$VOLUME / min_element
		branch$RAW_TIPS_VOL = branch$VOLUME
		branch$FULL_TIPS_VOL = branch$VOLUME

		branch$RAW_V_TOT = 0
		branch$V_TOT = 0 #branch$VOLUME - Exclude terminal branch volume from V_TOT to ensure empirical regressions are independent
		branch$FULL_V_TOT

		branch$N_SYM = 0
		branch$N_ASYM = 0
		branch$N_MEAN = branch$BIN_ORDER#0
		branch$N_MIN = branch$BIN_ORDER#0
		branch$N_MAX = branch$BIN_ORDER#0
		
		branch$TIPS_VAR = NA
		branch$TIPS_CV = NA
		branch$TIPS_GEO_CV = NA
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
	
		branch$V_TOT = branches[ind,]$V_TOT + (branch$VOLUME / min_element)

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

		vol_tots = vector(length = child_n)
		raw_volumes = vector(length = child_n)
		tips = vector(length = child_n)
		tips_vol = vector(length = child_n)
		raw_tips_vol = vector(length = child_n)
		
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
			raw_tips_vol[i] = branches[ind,]$RAW_TIPS_VOL
			vol_tots[i] = branches[ind,]$V_TOT
			raw_volumes[i] = branches[ind,]$RAW_V_TOT

			lengths[i] = branches[ind,]$LENGTH
			radii[i] = branches[ind,]$RADIUS
			volumes[i] = branches[ind,]$VOLUME 
			
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
		  volsort = sort(vol_tots)
		  ind_1 = child_n
		  ind_2 = child_n-1
		  pruned_volume = sum(vol_tots[-c(ind_1, ind_2)])
		  c_v = c_v + pruned_volume
		}
		
		branch$TIPS = tips[ind_1] + tips[ind_2]
		branch$FULL_TIPS = sum(tips)

		branch$TIPS_VOL = tips_vol[ind_1] + tips_vol[ind_2]	#+ tips_vol[i]
		branch$RAW_TIPS_VOL = raw_tips_vol[ind_1] + raw_tips_vol[ind_2]
		branch$FULL_TIPS_VOL = sum(raw_tips_vol)

		#These should be summed in the same way as the FULL_TIPS/FULL_TIP_VOL variables...
		#branch$CHILD_LENGTHS = lengths[ind_1] + lengths[ind_2]
		#branch$CHILD_RADII = radii[ind_1] + radii[ind_2]
		#branch$CHILD_VOLS = volumes[ind_1] + volumes[ind_2]

		#These should be summed in the same way as the FULL_TIPS/FULL_TIP_VOL variables...
		branch$CHILD_LENGTHS = sum(lengths)
		branch$CHILD_RADII = sum(radii)
		branch$CHILD_VOLS = sum(volumes)

		#Reminder that these are still 'symmetrical' in that size differences are averaged over when we iterate
		#the N counter by 1 at each bifurcation. In a symmetrical tree these would all be equivalent
		branch$N_MEAN = mean(n_mean) + 1
	  	branch$N_MIN = min(n_min) + 1
	  	branch$N_MAX = max(n_max) + 1
		
		branch$V_TOT = vol_tots[ind_1] + 
				vol_tots[ind_2] + 
				(branch$VOLUME/min_element)

		branch$RAW_V_TOT = raw_volumes[ind_1] + 
				raw_volumes[ind_2] + branch$VOLUME

		branch$FULL_V_TOT = sum(raw_volumes) + branch$VOLUME

		child_tips = children[which(branches[children,]$TIPS == 1)]
		
		branch$TIPS_MEAN = gm_mean(branches[child_tips,]$RAW_TIPS_VOL, na.rm=TRUE)
		branch$TIPS_VAR = var(branches[child_tips,]$RAW_TIPS_VOL, na.rm=TRUE)

		#Sample standard deviation divided by sample geometric mean - simplest form of CV
		branch$TIPS_CV = sqrt(branch$TIPS_VAR) / mean(branches[child_tips,]$RAW_TIPS_VOL, na.rm = TRUE)
		# 'Geometric coefficient of variation' for log-normal data
		#https://en.wikipedia.org/wiki/Coefficient_of_variation#Log-normal_data
		branch$TIPS_GEO_CV = sqrt(exp(var(log(branches[child_tips,]$RAW_TIPS_VOL), na.rm = TRUE)) - 1)

		#Expected number of terminal tips
		branch$VOL_MEAN = branch$V_TOT / branch$TIPS_MEAN
		
		branch$RADS_MEAN = gm_mean(branches[child_tips,]$RADIUS, na.rm=TRUE)

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

		#Average scale factors - holdover from asymmetry analyses. Better way to do this?
		#Different gamma for each child branch? Include more child branches in the average?
		#r_c = (radii[ind_1] + radii[ind_2])/2
		#l_c = (lengths[ind_1] + lengths[ind_2])/2
		r_c = mean(radii)
		l_c = mean(lengths, na.rm=TRUE)

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
		#valid_subtree <- c(id,children)[which(branches[as.integer(c(id,children)),]$VOLUME > 0)]
		#slope = lm(log(TIPS_VOL)~BIN_ORDER, branches[valid_subtree,])
	  	#a_b_r = exp(abs(slope$coefficients[[2]]))
	  	#if(is.na(a_b_r)) Sometimes parents/children occupy the same bin
	  	#{
	  	#	print(branch)
	  	#	print(slope)
	  	#	print(branch_dist)
	  	#}
	  	#branch$N_ASYM = log(branch$TIPS) / log(a_b_r)
	  	#branch$A_BR = a_b_r
	  	branch$N_ASYM = NA
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

basic_wbe <- function(beta, gamma)
{
	return(-log(2) / log(beta*beta*gamma))
}

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
whole_tree_scaling <- function(branches)
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
  
  #Tree-level variance in the size of terminal tips
  tips_mean = branches[1,]$TIPS_MEAN
  #tips_var = mean(branches[non_tips,]$TIPS_VAR, na.rm=TRUE)
  tips_var = branches[1,]$TIPS_VAR
  tips_cv = mean(branches[non_tips,]$TIPS_CV, na.rm=TRUE)
  tips_geo_cv = mean(branches[non_tips,]$TIPS_GEO_CV, na.rm=TRUE)
  
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
  beta_unfiltered = mean(branches$BETA, na.rm = TRUE)
  beta_ci = t.test(branches[valid_nodes,]$BETA, mu=0.7071)
  beta_ci_min = beta_ci$conf.int[1]
  beta_ci_max = beta_ci$conf.int[2]

  gamma = mean(branches[valid_nodes,]$GAMMA, na.rm = TRUE)
  inf = which(branches$GAMMA == Inf)
  gamma_unfiltered = -1
  if(length(inf) > 0)
  	gamma_unfiltered = mean(branches$GAMMA[-inf], na.rm = TRUE)
  else
  	gamma_unfiltered = mean(branches$GAMMA, na.rm = TRUE)
  gamma_ci = t.test(branches[valid_nodes,]$GAMMA, mu=0.7938)
  gamma_ci_min = gamma_ci$conf.int[1]
  gamma_ci_max = gamma_ci$conf.int[2]

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
  
  simple_theta = basic_wbe(beta, gamma)
  wbe = asymptotic_formula(sym_vol, N)
  asym = asymptotic_formula(vol_scaling, N)
  
  #Standardize data for regressions - normalize and scale so that we can compare everything
  branches[non_tips,]$V_TOT <- branches[non_tips,]$V_TOT / min(branches[non_tips,]$V_TOT)
  branches[non_tips,]$VOL_MEAN <- branches[non_tips,]$VOL_MEAN / min(branches[non_tips,]$VOL_MEAN)
  branches[non_tips,]$TIPS_VOL <- branches[non_tips,]$TIPS_VOL / min(branches[non_tips,]$TIPS_VOL)
  #TIPS on the same scale as tip volume for comparng regressions - scaling the axes won't matter for log-log regressions
  branches[non_tips,]$TIPS <- (branches[non_tips,]$TIPS / 2) 
  
  empirical_vol = sma(formula=log(TIPS_VOL)~log(V_TOT), data=branches[non_tips,], method="SMA")
  ci_vol = empirical_vol$slopetest[[1]]$ci
  branches$PRED_VOL <- NA
  branches[non_tips,]$PRED_VOL = empirical_vol$coef[[1]][1][1,] + (log(branches[non_tips,]$V_TOT) * empirical_vol$coef[[1]][1][2,])
  empirical_vol = empirical_vol$coef[[1]][2,1]  

  raw_empirical_vol = sma(formula=log(RAW_TIPS_VOL)~log(RAW_V_TOT), data=branches[non_tips,], method="SMA")
  branches$RAW_PRED_VOL <- NA
  branches[non_tips,]$RAW_PRED_VOL = raw_empirical_vol$coef[[1]][1][1,] + (log(branches[non_tips,]$RAW_V_TOT) * raw_empirical_vol$coef[[1]][1][2,])
  empirical_elevation = raw_empirical_vol$coef[[1]][1,1]

  ### Failed analyses I'm too afraid to delete

  tips_volume_mean = sma(formula=log(TIPS_VOL/TIPS)~log(TIPS_MEAN), data=branches[non_tips,], method="SMA")
  branches$PRED_TVM <- NA
  branches[non_tips,]$PRED_TVM = (tips_volume_mean$coef[[1]][1][1,]) + (log(branches[non_tips,]$TIPS_MEAN) * tips_volume_mean$coef[[1]][1][2,])
  tips_volume_mean = tips_volume_mean$coef[[1]][2,1]

  tips_branch_ratio = sma(formula=log(TIPS_VOL/TIPS_MEAN)~log(TIPS), data=branches[non_tips,], method="SMA")
  branches$PRED_BR <- NA
  branches[non_tips,]$PRED_BR = (tips_branch_ratio$coef[[1]][1][1,]) + (log(branches[non_tips,]$TIPS) * tips_branch_ratio$coef[[1]][1][2,])
  tips_branch_ratio = tips_branch_ratio$coef[[1]][2,1]

  tips_volume = sma(formula=log((TIPS_MEAN^(tips_volume_mean))*(TIPS^(tips_branch_ratio)))~log(TIPS_VOL), data=branches[non_tips,], method="SMA")
  branches$PRED_TV <- NA
  branches[non_tips,]$PRED_TV = (tips_volume$coef[[1]][1][1,]) + (log(branches[non_tips,]$TIPS_VOL) * tips_volume$coef[[1]][1][2,])
  tips_volume = tips_volume$coef[[1]][2,1]

  empirical = sma(formula=log(TIPS)~log(V_TOT), data=branches[non_tips,], method="SMA")
  ci = empirical$slopetest[[1]]$ci
  branches$PRED <- NA
  branches[non_tips,]$PRED = (empirical$coef[[1]][1][1,]) + (log(branches[non_tips,]$V_TOT) * empirical$coef[[1]][1][2,])
  empirical = empirical$coef[[1]][2,1]
	
  corrected_empirical = sma(formula=log((TIPS_MEAN^(tips_volume_mean))*(TIPS^(tips_branch_ratio)))~log(V_TOT), data=branches[non_tips,], method="SMA")
  branches$PRED_CE <- NA
  branches[non_tips,]$PRED_CE = corrected_empirical$coef[[1]][1][1,] + (log(branches[non_tips,]$V_TOT) * corrected_empirical$coef[[1]][1][2,])
  corrected_empirical = corrected_empirical$coef[[1]][2,1]	  

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

  normalisation = branches[1,]$RAW_TIPS_VOL / (branches[1,]$RAW_V_TOT ^ empirical_vol)

  #This should always be returning simple scalar summary stats
  return(list("beta"=beta, "beta_ci_min"=beta_ci_min, "beta_ci_max"=beta_ci_max, "gamma"=gamma, "gamma_ci_min"=gamma_ci_min, "gamma_ci_max"=gamma_ci_max, 
  			  "d_beta"=d_beta, "d_gamma"=d_gamma, "fib_r" = fib_r, "fib_l"=fib_l, 
  			  "simple_theta"=simple_theta, "wbe"=wbe, "node_wbe"=node_wbe, "asym"=asym, "node_asym"=node_asym, "sym_vol" = sym_vol, "asym_vol" = asym_vol, 
  			  "network_n"=N_SYM, "n_asym"=N_ASYM, "n_mean"=N_MEAN, "n_min"=N_MIN, "n_max"=N_MAX, 
  			  "tips_mean"=tips_mean, "tips_variance"=tips_var, "tips_cv"=tips_cv, "tips_geo_cv"=tips_geo_cv,
			  "tips_volume"=tips_volume, "tips_volume_mean"=tips_volume_mean, "tips_branch_ratio"=tips_branch_ratio,
			  "empirical"=empirical, "empirical_corrected"=corrected_empirical, "empirical_vol"=empirical_vol, 
			  "empirical_elevation"=empirical_elevation, "normalisation"=normalisation, "empirical_vol_mean"=empirical_vol_mean,
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

	#Initiate recursion to process the network
	out <- asymmetric_branching(1, branches, 0, c())
	branches <- out$data
	
	#Compute whole-tree scaling ratios - compare to average scaling ratios from nodes	
	scaling = whole_tree_scaling(branches)
	branches <- scaling$preds_tips
  	scaling <- scaling[-length(scaling)] #Get rid of last element there
  	
	  if(verbose_report)
	  {
	    print(paste("Volume pruned: ", out$c_v))
	    print(paste("Total volume after pruning: ", branches[1,]$V_TOT))
	    
	    print(paste("Percent invalid nodes: ", length(which(branches$INVALID)) / nrow(subset(branches, TIPS>1)) ))
	    
	    print(scaling)
	  }
  
	return (list("branches"=branches, "scaling"=scaling))
}
