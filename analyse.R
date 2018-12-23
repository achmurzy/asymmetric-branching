library(gdata)
library(reshape)
library(smatr)


make_internal_frame <- function(rows)
{
frame <- data.frame(BRANCH_ORDER=integer(rows), BRANCH_ID=integer(rows), CHILD_IDS=character(rows), PARENT_ID=character(rows), 
                    TIPS=integer(rows), V_TOT=numeric(rows), VOLUME = numeric(rows), L_TOT=numeric(rows), LENGTH=numeric(rows), 
                    RADIUS=numeric(rows), BETA=numeric(rows), GAMMA=numeric(rows), D_BETA=numeric(rows), D_GAMMA=numeric(rows), 
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
		return (list("data"=branches, "c_v"=c_v))
	}
	else if (child_n == 1)
	{
		ind = as.numeric(children[1])
		out <- asymmetric_branching(ind, branches, c_v)
		branches = out$data
		c_v = out$c_v
		
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
	
	  #Not sure about this
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

	branch$WBE_THETA <- nodal_wbe(2, branch$BETA, branch$GAMMA)
	branch$THETA <- nodal_asym(branch$BETA, branch$GAMMA, branch$D_BETA, branch$D_GAMMA)
	#branch$INVALID = branch$WBE_THETA <= 0 || branch$WBE_THETA >= 1 || branch$THETA <= 0 || branch$THETA >= 1
	branch$INVALID = branch$BETA >= 1 || branch$GAMMA >= 1
	branches[id,] <- branch
	return (list("data"=branches, "c_v"=c_v))
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
nodal_asym <- function(beta, gamma, d_beta, d_gamma)
{
  sym = log(gamma*(beta^2))/log(2)
	t_1 = (2*d_beta*d_gamma/(beta*gamma))
	t_2 = ((d_beta^2)/beta^2)
	asym = (1 + t_1 + t_2)/log(2) 
	theta = -(sym + asym)^(-1)
	return(theta)
}

#Whole-tree scaling exponents computed from average ratios across the network
#according to Brummer et al 2019 in prep
whole_tree_scaling <- function(branches, verbose_report = TRUE)
{
  #Compute branching generations N (network depth) assuming bifurcation
  num_tips = branches[1,]$TIPS
  N = log(num_tips)/log(2)
  
  branches <- branches[which(!branches$INVALID),]
  
  beta = mean(branches$BETA, na.rm = TRUE)
  gamma = mean(branches$GAMMA, na.rm = TRUE)
  
  #Positive/negative asymmetry may be washing out our aggregation of asymmetry
  d_beta = mean(branches$D_BETA, na.rm = TRUE)
  d_gamma = mean(branches$D_GAMMA, na.rm = TRUE)
  
  #Volume scaling expressions - ratio of total volume from sibling branches to parent branches
  sym_vol = 2*(beta^2)*gamma
  asym_vol = 2*(beta^2)*gamma + 4*gamma*d_beta*d_gamma + 2*gamma*(d_beta^2)
  
  #Scaling exponents
  if(abs(1-sym_vol)>0.1)
  {
    wbe = log(2^N) / (log(2^N) + 
                        log(1 - sym_vol^(N+1)) - #Generates NaN - sym_vol^(N+1) > 1 when sym_vol > 1
                          log((1-sym_vol)*sym_vol^(N))) #Generates NaN - sym_vol > 1
  }
  else
  {
    wbe = log(2^N) / ((log(N+1)*(2^N)) - N*log(sym_vol))
  }
  
  if(abs(1-asym_vol)>0.1)
  {
    asym = log(2^N)/ (log(2^N) + 
                        log(1 - asym_vol^(N+1)) - #Generates NaN
                          log((1-asym_vol)*asym_vol^(N))) #Generates NaN
  }
  else
  {
    asym = log(2^N) / ((log(N+1)*(2^N)) - N*log(asym_vol))
  }
  
  empirical = sma(formula=TIPS~V_TOT, data=branches, log='xy')
  empirical = empirical$coef[[1]][2,1]
  
  if(verbose_report)
  {
    print("Finite network correction assuming bifurcation...")
    print(paste("Total number of tips:", num_tips))
    print(paste("Inferred network generations:", N))
    
    print(paste("Only using valid nodes for whole-tree scaling:", nrow(branches)))
    
    print(paste("Average beta:", beta))
    print(paste("Average gamma:", gamma))
    print(paste("Average d_beta:", d_beta))
    print(paste("Average d_gamma:", d_gamma))
    
    print(paste("Symmetric volume scaling:", sym_vol))
    print(paste("Asymmetric volume scaling:", asym_vol))
  }
  
  return(list("beta"=beta, "gamma"=gamma, "d_beta"=d_beta, "d_gamma"=d_gamma, "wbe"=wbe, "asym"=asym, "empirical"=empirical))
}

branching_analysis <- function(branches, verbose_report = TRUE)
{
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
