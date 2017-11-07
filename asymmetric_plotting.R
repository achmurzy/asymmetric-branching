plot_volumes <- function(ll)
{
	count = 0
	names <- names(ll$species)
	lapply(ll$species, function(x) 
	{  
		count = count + 1
		plot_volume(x, TRUE, names[count]) 
	})
}

#Plots currently heavily weighted by network tips with high volume
plot_volume <- function(tree, single, name)
{
	form <- y ~ C*(1-exp(k*x))	
	form_line <- y~x
	
	normal_factor <- min(tree$V_TOT)
	print(normal_factor)
	#tree <- subset(tree, TIPS > 1 | V_TOT < normal_factor*100)
	v_tips<-ggplot(data=tree, 
		aes(x=log(V_TOT/normal_factor), y=log(TIPS))) +
	geom_point(alpha=0.5) + scale_y_continuous(limits= c(0, NA)) +
	scale_x_continuous(limits=c(0,NA)) +
	geom_smooth(method = 'lm', formula=form_line, se=FALSE) +
	stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., 			sep = "~~~")),label.x.npc = "left", label.y.npc = 0.3,
       		formula=form_line, parse = TRUE, size = 3)
	if(single)
	{
		ss <- mean_theta(summary_stats(tree))
		v_tips <- v_tips + 
		geom_abline(slope=ss$asym, intercept=0, 				colour='#E41A1C') +
		geom_abline(slope=ss$wbe, intercept=0, colour='#00ff00')
	}	
	#print(v_tips)

	#Test against regression line?
	ggsave(paste(name, "Vtot_tips_log.png",sep="_"), v_tips)
}

prediction_plot <- function(tree)
{
	form_line <- y~x
	
	normal_factor <- min(tree$V_TOT)
	tree["Volume_cm^3"]<- log(tree$V_TOT/normal_factor)
	tree["Num_Tips"] <- log(tree$TIPS)
	v_tips<-ggplot(data=tree,
	aes_string(x=quote(tree["Volume_cm^3"]), y=quote(tree["Num_Tips"]))) +
	geom_point(alpha=0.5) + scale_y_continuous(limits= c(0, NA)) +
	scale_x_continuous(limits=c(0,NA)) +
	geom_smooth(method = 'lm', formula=form_line, se=FALSE) + 		geom_abline(slope=0.35, intercept=-1.5, colour='#E41A1C') +
	geom_abline(slope=0.15, intercept=-1.5, colour='#00ff00') 
	
	
	print(v_tips)	
	ggsave(paste("Prediction_Vtot_tips_log.png",sep="_"), v_tips)
}

plot_exponents <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	exponents <- ggplot(data=tree, aes(x=BRANCH_ORDER, y=THETA)) +
	geom_point(alpha=0.5) +
	scale_y_continuous(limits=c(0, 1))
	print(exponents)
	ggsave("subtree_exponent.png", exponents)
}

plot_asymmetry <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	asymmetry <- ggplot(data=tree, aes(x=D_BETA, y=D_GAMMA, color=THETA)) + scale_x_continuous(limits=c(-1, 1)) + scale_y_continuous(limits=c(-1, 1)) + scale_colour_gradient(low="purple", high="yellow") +
	geom_point(size=1)
	print(asymmetry)
	ggsave("subtree_asymmetry.png", asymmetry)
}
	

plot_symmetry <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	symmetry <- ggplot(data=tree, aes(x=BETA, y=GAMMA, color=THETA)) + scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1)) + scale_colour_gradient(low="purple", high="yellow") +
	geom_point(size=1)
	print(symmetry)
	ggsave("subtree_symmetry.png", symmetry)
}

plot_lengths <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	lengths <- ggplot(data=tree, aes(x=BRANCH_ORDER, y=LENGTH)) + 		geom_point(size=1)
	print(lengths)
	ggsave("subtree_lengths.png", lengths)
}

plot_radii <- function(tree)
{
	tree <- subset(tree, !is.na(THETA))
	radii <- ggplot(data=tree, aes(x=BRANCH_ORDER, y=RADIUS)) + 		geom_point(size=1)
	print(radii)
	ggsave("subtree_radii.png", radii)
}
