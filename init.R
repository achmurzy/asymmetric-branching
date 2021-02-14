library(rmarkdown)
library(MuMIn)
library(smatr)

source('munge.R')
source('analyse.R')
source('plot.R')

#This script is meant to reproduce the main results/figures of the paper without re-running the entire recursive tree analysis in workflow.R

load('data/cylinder_data.RData')
#load('data/branch_data.RData')
load('data/tree_data.RData')

get_data <- function()
{
	gg <- merge_guyanese_trees()
	tree_data <- gg$trees
	cylinder_data <- gg$cyls
	tree_data <- faceted_full_volume_scaling(tree_data, cylinder_data) # These should all just get done in analyse.R
	tree_data <- faceted_full_tip_scaling(tree_data, cylinder_data)
	tree_data <- length_radius_cvs(tree_data)
	hand_trees <- analyse_destructive_sampling()
	return (list("trees"=tree_data, "cyls"=cylinder_data))
}

merge_guyanese_trees <- function()
{
	tree_data <- tree_data[,c(1:66)] #Exclude traits to rbind guyana trees
	gg <- analyse_guyanese_trees()
	guy_trees <- gg$trees
	tree_data <- rbind(tree_data, gg$trees)
	cylinder_data <- rbind(cylinder_data, gg$cyls)
	return (list("trees"=tree_data, "cyls"=cylinder_data))	
}

branching_ratios <- function(trees, ordering)
{
	filenames <- unique(cylinder_data$FILENAME)
	trees$BR <- rep(0, nrow(trees))
	for(x in filenames)
	{
		branches <- which(cylinder_data$FILENAME == x)
		order_count <- count(cylinder_data[branches,], eval(sym(ordering)))
		colnames(order_count) <- c(ordering, "count")
		mm <- lm(log(count)~eval(sym(ordering)), order_count)
		br = exp(abs(mm$coefficients[[2]]))
		
		tree <- which(trees$FILENAME == x)
		if(length(tree) > 0)
			trees[tree,]$BR <- br
	}
	return(trees)
}

#Required to run figure F - TODO: integrate CVs into regular analysis workflow
length_radius_cvs <- function(tree_data)
{
	unique_names <- unique(tree_data$FILENAME)
	tree_data$R_CV <- rep(0, nrow(tree_data))
	tree_data$L_CV <- rep(0, nrow(tree_data))
	tree_data$H_CV <- rep(0, nrow(tree_data))
	cylinder_data$H_CV <- rep(0, nrow(cylinder_data))
	cylinder_data$TIP_BIN <- rep(0, nrow(cylinder_data))
	for(x in unique_names)
	{
		cylinders <- cylinder_data[which(cylinder_data$FILENAME == x),]
		tips <- cylinders[which(cylinders$TIPS == 1),]
		
		lambda = 1.2
		tree_log_bin <- log_bin(cylinders$TIPS,lambda=lambda)	
		subtrees = log(tree_log_bin/min(cylinders$TIPS)) / log(lambda)
		cylinders$TIP_BIN <- floor(subtrees)
		bins <- unique(cylinders$TIP_BIN)
		h_cvs <- vector(length = length(bins))
		ind = 0
		for(y in bins)
		{
			ind = ind + 1
			subs <- cylinders[which(cylinders$TIP_BIN == y),]
			h_cvs[ind] = sd(subs$LENGTH, na.rm = TRUE) / gm_mean(subs$LENGTH, na.rm = TRUE)
			cylinders[which(cylinders$TIP_BIN == y),]$H_CV <- h_cvs[ind] 
		}
		cylinder_data[which(cylinder_data$FILENAME == x),] <- cylinders
		r_cv = sd(tips$RADIUS, na.rm = TRUE) / gm_mean(tips$RADIUS, na.rm = TRUE)
		l_cv = sd(tips$LENGTH, na.rm = TRUE) / gm_mean(tips$LENGTH, na.rm = TRUE)
		
		h_cv = mean(h_cvs, na.rm = TRUE)
		tree_data[which(tree_data$FILENAME == x),]$R_CV = r_cv
		tree_data[which(tree_data$FILENAME == x),]$L_CV = l_cv
		tree_data[which(tree_data$FILENAME == x),]$H_CV = h_cv
	}

	ggplot(cylinder_data[which(!cylinder_data$INVALID),], aes(x=H_CV, y=GAMMA)) + geom_point() + geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
 		stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
              label.x.npc = "left", label.y.npc = 0.7, formula=y~x, parse = TRUE, size = 5) +
 		facet_wrap(~FILENAME)
 	#ggsave('figures/H_CV', device='png', width=32, height=32)   	
	
	return(tree_data)
}

wideScreen <- function(howWide=Sys.getenv("COLUMNS")) {
  options(width=as.integer(howWide))
}

convert_to_wbe <- function(val)
{
	return(exp(-val * log(2)))
}

table_one_results <- function()
{
	#Apply the 10cm branch diameter cut-off from Lau et al
	gg <- analyse_guyanese_trees(TRUE)
	guy_trees_cutoff <- faceted_full_volume_scaling(gg$trees, gg$cyls)
	guy_trees_cutoff <- faceted_full_tip_scaling(guy_trees_cutoff, gg$cyls)

	#Whole-tree medians from table A.3 of Lau et al 2019 - Manual and LiDAR versus theoretical expectation and versus each other
	#Betas
	t.test(convert_to_wbe(c(0.49, 0.45, 0.43, 0.44, 0.47, 0.5, 0.5, 0.52, 0.33)), mu=0.7071)
	t.test(convert_to_wbe(c(0.66, 0.62, 0.69, 0.68, 0.62, 0.65, 0.66, 0.62, 0.53)), mu=0.7071)
	t.test(convert_to_wbe(c(0.49, 0.45, 0.43, 0.44, 0.47, 0.5, 0.5, 0.52, 0.33)), convert_to_wbe(c(0.66, 0.62, 0.69, 0.68, 0.62, 0.65, 0.66, 0.62, 0.53)))

	#Gammas
	t.test(convert_to_wbe(c(0.21, 0.59, 0.35, 0.17, 0.55, 0.62, 0.44, 0.57, 0.23)), mu=0.7938)
	t.test(convert_to_wbe(c(-0.12, 0.3, 0.13, -0.33, 0.42, -0.19, 0.28, 0.05, -0.11)), mu=0.7938)
	t.test(convert_to_wbe(c(0.21, 0.59, 0.35, 0.17, 0.55, 0.62, 0.44, 0.57, 0.23)), convert_to_wbe(c(-0.12, 0.3, 0.13, -0.33, 0.42, -0.19, 0.28, 0.05, -0.11)))
	
	#Thetas
	t.test(c(0.69, 0.56, 0.38, 0.24, 0.51, 0.75, 0.53, 0.59, 0.77), mu=0.75)
	t.test(c(0.52, 0.64, 0.47, 0.38, 0.38, 0.76, 0.27, 0.71, 0.45), mu=0.75)
	t.test(c(0.52, 0.64, 0.47, 0.38, 0.38, 0.76, 0.27, 0.71, 0.45), c(0.69, 0.56, 0.38, 0.24, 0.51, 0.75, 0.53, 0.59, 0.77))

	t.test(guy_trees_cutoff$BETA, mu=0.7071)
	t.test(guy_trees_cutoff$GAMMA, mu=0.7938)
	t.test(guy_trees_cutoff$FULL_TIPS, mu=0.75)
	t.test(guy_trees_cutoff$FULL_VOL, mu=0.75)

	t.test(unique(hand_trees$mean_beta), mu=0.7071)
	t.test(unique(hand_trees$mean_gamma), mu=0.7938)
	t.test(unique(hand_trees$tip_exponent), mu=0.75)
	t.test(unique(hand_trees$vol_exponent), mu=0.75)

	t.test(unique(hand_trees$mean_beta), guy_trees_cutoff$BETA)
	t.test(unique(hand_trees$mean_gamma), guy_trees_cutoff$GAMMA)
	t.test(unique(hand_trees$tip_exponent), guy_trees_cutoff$FULL_TIPS)
	t.test(unique(hand_trees$vol_exponent), guy_trees_cutoff$FULL_VOL)

	t.test(tree_data$BETA, mu=0.7071)
	t.test(tree_data$LENGTH_PRES, mu=0.7938) #Requires column from Figure C to run
	t.test(tree_data$FULL_TIPS, mu=0.75)
	t.test(tree_data$FULL_VOL, mu=0.75)
}

#Apply a radius cutoff to reproduce Lau et al 2019, or not
analyse_guyanese_trees <- function(cutoff = FALSE)
{
	#Manually add Guyanese trees to tree_data
	guyanese_names <- c("GUY01", "GUY02", "GUY03", "GUY04", "GUY05", "GUY06", "GUY07", "GUY08", "GUY09", "GUY10")
	guyana_files <- list.files(path="data/guyana/")
	num_trees <- length(guyanese_names)
	guy_trees <- data.frame(FILENAME=character(num_trees), NUM_BRANCHES=integer(num_trees), MAX_BRANCH_ORDER=integer(num_trees),
                        BRANCH_VOLUME=numeric(num_trees), TREE_HEIGHT=numeric(num_trees), TREE_VOLUME=numeric(num_trees),  
                    TRUNK_VOLUME=numeric(num_trees), TRUNK_LENGTH=numeric(num_trees), BRANCH_LENGTH=numeric(num_trees),  
                    TREE_AREA=numeric(num_trees), DBH_QSM=numeric(num_trees), DBH_CYL=numeric(num_trees), BETA=numeric(num_trees),
                    BETA_CI_MIN=numeric(num_trees), BETA_CI_MAX=numeric(num_trees), GAMMA=numeric(num_trees), GAMMA_CI_MIN=numeric(num_trees),
                    GAMMA_CI_MAX=numeric(num_trees), D_BETA=numeric(num_trees), D_GAMMA=numeric(num_trees), FIB_R=numeric(num_trees),
                    FIB_L=numeric(num_trees), SIMPLE_THETA=numeric(num_trees), WBE=numeric(num_trees), NODE_WBE=numeric(num_trees), THETA=numeric(num_trees), 
                    NODE_THETA=numeric(num_trees), SYM_VOL=numeric(num_trees), ASYM_VOL=numeric(num_trees), NETWORK_N=integer(num_trees),
                    N_ASYM=numeric(num_trees), N_MEAN=numeric(num_trees), N_MIN=numeric(num_trees), N_MAX=numeric(num_trees), 
                    TIPS_MEAN=numeric(num_trees), TIPS_VARIANCE=numeric(num_trees), TIPS_CV = numeric(num_trees), TIPS_GEO_CV=numeric(num_trees),
                    TIPS_VOLUME=numeric(num_trees), TIPS_VOLUME_MEAN=numeric(num_trees), TIPS_BRANCH_RATIO=numeric(num_trees),
			              EMPIRICAL=numeric(num_trees), EMPIRICAL_CORRECTED=numeric(num_trees), EMPIRICAL_VOL=numeric(num_trees), 
                    EMPIRICAL_INTERCEPT=numeric(num_trees), NORMALISATION=numeric(num_trees),
                    EMPIRICAL_VOL_MEAN=numeric(num_trees), EMPIRICAL_VOL_TIP=numeric(num_trees), EMPIRICAL_TIP_MEAN=numeric(num_trees),
		 	              CI_MIN=numeric(num_trees), CI_MAX=numeric(num_trees), CI_VOL_MIN=numeric(num_trees), CI_VOL_MAX=numeric(num_trees),		
		                CI_MA_MIN=numeric(num_trees), CI_MA_MAX=numeric(num_trees), PATH_FRAC=numeric(num_trees), ASYM_FRAC=numeric(num_trees), 
                    BRANCH_ANGLE=numeric(num_trees),
                    stringsAsFactors = FALSE)
	guy_cyls <- data.frame()
	for(i in seq(1, length(guyanese_names)))
	{
	  cyl_data <- open_TreeQSM(paste("data/guyana", guyana_files[i], sep="/"))
	  cyl_data <- munge_TreeQSM(cyl_data)
	  if(cutoff)
	  {cyl_data[which(cyl_data$RADIUS*2 < 0.10),]$INVALID <- TRUE}
	  out <- branching_analysis(cyl_data, verbose_report = FALSE)
	  cyl_data <- out$branches
	  scaling <- out$scaling
	  cyl_data$FILENAME <- rep(guyanese_names[i], nrow(cyl_data))
	  guy_cyls <- rbind(guy_cyls, cyl_data)

	  new_tree <- c(guyanese_names[i], rep(NA, 11), as.numeric(unlist(scaling)), NA)
  	  guy_trees[i,] <- new_tree
	}
	guy_trees[,3:ncol(guy_trees)] <- sapply(guy_trees[,3:ncol(guy_trees)], as.numeric)
	tree_metadata <- read.csv("data/tree_metadata.csv")
	tree_metadata$FILENAME <- as.character(tree_metadata$FILENAME)
	guy_trees <- merge(guy_trees, tree_metadata, by="FILENAME", all=FALSE)
	guy_trees$BINOMIAL <- paste(guy_trees$GENUS, guy_trees$SPECIES)
	return(list("trees"=guy_trees, "cyls"=guy_cyls))
	#Compare LiDAR to hand-measured trees for both regression types. 
}

analyse_destructive_sampling <- function()
{
	destructive_data = "data/destructive_sampling"
	files = list.files(path=destructive_data, pattern="*.csv")
	myfiles = do.call(rbind, lapply(files, FUN=function(x) read.csv(paste(destructive_data, x, sep="/"), stringsAsFactors = FALSE)))
	myfiles$r_avg <- (myfiles$Diameter_base + myfiles$Diameter_top) / 4	#Convert to radius
	myfiles$volume <- pi * myfiles$r_avg * myfiles$r_avg * myfiles$Length
	myfiles$total_volume <- rep(NA, nrow(myfiles))
	myfiles$child_volume <- rep(NA, nrow(myfiles))
	myfiles$tip_volume <- rep(NA, nrow(myfiles))
	myfiles$tips <- rep(NA, nrow(myfiles))
	myfiles$beta <- rep(NA, nrow(myfiles))
	myfiles$gamma <- rep(NA, nrow(myfiles))
	myfiles$mean_beta <- rep(NA, nrow(myfiles))
	myfiles$mean_gamma <- rep(NA, nrow(myfiles))
	myfiles$tip_exponent <- rep(NA, nrow(myfiles))
	myfiles$vol_exponent <- rep(NA, nrow(myfiles))
	myfiles$child_parent_scaling <- rep(NA, nrow(myfiles))

	filenames <- unique(myfiles$TreeID)
	for(x in filenames)
	{
		print(paste("Processing tree", x))
		nodes <- myfiles[which(myfiles$TreeID == x),]
		trunk <- which(nodes$branch_order == 0)
		new_nodes <- recurse_tree(trunk, nodes)
		tip_scaling = sma(formula=log(tips)~log(total_volume), data=new_nodes, method="SMA")
		vol_scaling = sma(formula=log(tip_volume)~log(total_volume), data=new_nodes, method="SMA")
		child_parent_scaling = sma(formula=log(child_volume)~log(volume), data=new_nodes, method="SMA")
		myfiles[which(myfiles$TreeID == x),] <- new_nodes
		myfiles[which(myfiles$TreeID == x),]$tip_exponent <- tip_scaling$coef[[1]][2,1]
		myfiles[which(myfiles$TreeID == x),]$vol_exponent <- vol_scaling$coef[[1]][2,1]
		myfiles[which(myfiles$TreeID == x),]$child_parent_scaling <- child_parent_scaling$coef[[1]][2,1]
		myfiles[which(myfiles$TreeID == x),]$mean_beta <- mean(new_nodes$beta, na.rm=TRUE)
		myfiles[which(myfiles$TreeID == x),]$mean_gamma <- mean(new_nodes$gamma, na.rm=TRUE)
	}

	return(myfiles)
}

recurse_tree <- function(n_id, nodes)
{
	children <- which(nodes$Parent_Node == nodes[n_id,]$NodeID)
	if(length(children) > 0)
	{
		for(child in children)
		{
			nodes <- recurse_tree(child, nodes)
		}
		nodes[n_id,]$tips <- sum(nodes[children,]$tips)
		nodes[n_id,]$tip_volume <- sum(nodes[children,]$tip_volume)
		nodes[n_id,]$child_volume <- sum(nodes[children,]$volume)
		nodes[n_id,]$total_volume <- sum(nodes[children,]$total_volume) + nodes[n_id,]$volume
		nodes[n_id,]$beta <- mean(nodes[children,]$r_avg) / nodes[n_id,]$r_avg
		nodes[n_id,]$gamma <- mean(nodes[children,]$Length) / nodes[n_id,]$Length
	}
	else
	{
		nodes[n_id,]$tips <- 1
		nodes[n_id,]$tip_volume <- nodes[n_id,]$volume
		nodes[n_id,]$total_volume <- nodes[n_id,]$volume	
	}
	return (nodes)
}

#June 10th - THETA + CV is the most competitive model, but remains absolute garbage
#There's a way to throw everything in and iterate every permutation
# find it and justify the inclusion of each parameter
perform_model_selection <- function()
{
	model_theta <- lm(EMPIRICAL_VOL~THETA, data=tree_data)
	model_n_max <- lm(EMPIRICAL_VOL~N_MAX, data=tree_data)
	model_n_theta <- lm(EMPIRICAL_VOL~N_MAX*THETA, data=tree_data)
	model_n_beta_gamma <- lm(EMPIRICAL_VOL~N_MAX*BETA*GAMMA, data=tree_data)
	model_theta_var <- lm(EMPIRICAL_VOL~THETA+TIPS_VARIANCE, data=tree_data)
	model_theta_mean <- lm(EMPIRICAL_VOL~THETA+TIPS_MEAN, data=tree_data)
	model_theta_mean_var <- lm(EMPIRICAL_VOL~THETA+TIPS_MEAN+TIPS_VARIANCE, data=tree_data)
	model_theta_cv <- lm(EMPIRICAL_VOL~THETA+TIPS_CV, data=tree_data)
	
	out <- model.sel(model_theta, model_theta_mean, model_theta_var, model_theta_mean_var, model_theta_cv)
	
	#Most successful model: BETA and N_MAX by AIC
	options(na.action = "na.fail")
	inclusions = tree_data[,c("EMPIRICAL_VOL", "N_MAX", "NETWORK_N", "N_MIN", "BETA","GAMMA", "TIPS_VARIANCE", "TIPS_MEAN")]
	all <- lm(EMPIRICAL_VOL ~ ., data=inclusions)
	fishing <- dredge(all)
	return(fishing)
}

tip_scaling_ratios <- function(tree_data, cylinder_data)
{
	unique_names <- unique(tree_data$FILENAME)
	tree_data$EMP_PRED <- rep(0, nrow(tree_data))
	cyls <- data.frame()
	tree_s <- data.frame()
	for(i in seq(1, length(unique_names)))
	{
		x <- unique_names[i]
		cylinders <- cylinder_data[which(cylinder_data$FILENAME == x),]
		max_order <- max(cylinders$BRANCH_ORDER)
		tree_lookup <- data.frame(orders=seq(1,max_order), avg_lengths=rep(0, max_order))
		for(y in seq(1, max_order))	#TODO: Calculate tip variability for each branch order to implement another definition of gamma star
		{
			tree_lookup[y,2] = mean(cylinders[which(cylinders$BRANCH_ORDER == y),]$LENGTH)
		}
		cylinders$GAMMA_STAR = cylinders$LENGTH / tree_lookup[cylinders$BRANCH_ORDER,]$avg_lengths 
		infs <- which(cylinders$GAMMA_STAR == -Inf)
		tree_sub <- which(tree_data$FILENAME == x)
		xs = -log(2) / log(cylinders$GAMMA_STAR)
		mean_val = mean(xs[-1])
		print(mean_val)
		tree_data[tree_sub,]$EMP_PRED = mean_val
		cyls <- rbind(cyls, cylinders)
		tree_s <- rbind(tree_s, tree_data[tree_sub,])
	}
	return(list("trees"=tree_s, "cyls"=cyls))
}

#Comparing branch level TreeQSM output to cylinder-level data.
#While there are deviations, the cylinders well-approximate branch-level volume data
#branch_vol = c()
#for(name in tree_data$FILENAME)
#{
#	branch_names <- which(branch_data$FILENAME == name)
#	branch_vol <- c(branch_vol, sum(branch_data[branch_names,]$VOLUME))
#}
#Branch level volume data is given in liters (1000cm^3)
#ggplot(tree_data, aes(x=raw_v_tot, branch_vol/1000)) + geom_point()

volume_scaling_vectorized <- Vectorize(asymptotic_formula)

#What is this telling us?
ggplot(tree_data, aes(x=BETA, y=TIPS_CV)) + geom_point() + geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
     stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                  label.x.npc = "left", label.y.npc = 0.7, formula=y~x, parse = TRUE, size = 5) 
