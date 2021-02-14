#library(ggnetwork)
library(ggpmisc)
library(plyr)
library(dplyr)
library(tidyr)
library(smatr)
library(raster)
###
# Manuscript Plots
###

manuscript_theme = theme(
  text = element_text(size=24),
  strip.text = element_text(face = "italic")
  )

figure_M <- function()
{
  comp_tree <- cylinder_data[which(cylinder_data$FILENAME == "Ery_11"),]
  tree_tips <- which(comp_tree$TIPS > 1)
  comp_tree <- comp_tree[tree_tips,]
  
  mm <- sma(data=comp_tree, formula=log(TIPS)~log(V_TOT))
  print(mm)
  tip_slope <- tree_data[which(tree_data$FILENAME == filename),]$EMPIRICAL
  
  figure_aa <- ggplot(data=comp_tree, aes(x=log(V_TOT), y=log(TIPS))) +
  geom_point(alpha=0.5) +
  #geom_line(aes(x=log(V_TOT), y=PRED), color="blue") + 
  geom_abline(slope=mm$coef[[1]][1][2,], intercept=mm$coef[[1]][1][1,], color="blue") +
  #scale_x_continuous(limits = c(-1, 16)) + scale_y_continuous(limits= c(-1, 12)) +
  labs(x="log Subtree Volume", y="log Tip Count", title="Distal terminal tips by subtree")  
  #annotate(geom='text', x=0.5, y=7, color="blue", size=12, label=round(tip_slope,2))

  figure_aa <- figure_aa + manuscript_theme

  ggsave("figures/Figure_M.png")
}

#c("Prunus avium", "Erythrophleum fordii", "Erythrophleum suaveolens", "Ulmus americana", "Acer pseudoplatanus")
figure_A <- function(filenames)
{
  full_tree <- data.frame()
  for(file in filenames)
  {
    cyls <- which(cylinder_data$FILENAME == file)
    comp_tree <- cylinder_data[cyls,]
    t_tree <- tree_data[which(tree_data$FILENAME == file),]
    comp_tree$GENUS <- t_tree$GENUS
    comp_tree$SPECIES <- t_tree$SPECIES
    
    tree_tips <- which(comp_tree$FULL_TIPS > 1)
    comp_tree <- comp_tree[tree_tips,]
    mm <- sma(data=comp_tree, formula=log(FULL_TIPS)~log(FULL_V_TOT))
    step = 17 / (nrow(comp_tree)-1)
    print(paste(file, mm$coef[[1]][1][2,]))
    comp_tree$XS <- seq(-15, 2, step)
    comp_tree$TIP_PREDS <- mm$coef[[1]][1][1,] + (mm$coef[[1]][1][2,] * comp_tree$XS)
    
    mma <- sma(data=comp_tree, formula=log(FULL_TIPS_VOL)~log(FULL_V_TOT))
    print(paste(file, mma$coef[[1]][1][2,]))
    comp_tree$VOL_PREDS <- mma$coef[[1]][1][1,] + (mma$coef[[1]][1][2,] * comp_tree$XS)

    full_tree <- rbind(full_tree, comp_tree)
  }

  full_tree$BINOMIAL <- factor(paste(full_tree$GENUS, full_tree$SPECIES), 
    levels = c("Prunus avium", "Erythrophleum fordii", "Erythrophleum suaveolens", "Ulmus americana", "Acer pseudoplatanus"))

  figure_aa <- ggplot(data=full_tree, aes(x=log(FULL_V_TOT), y=log(FULL_TIPS))) +
  geom_point(alpha=0.5) +
  geom_line(aes(x=XS, y=TIP_PREDS), color="blue") +
  scale_x_continuous(limits = c(-15, 2)) + scale_y_continuous(limits= c(0, 7)) +
  labs(x="log Subtree Volume", y="log Total Subtree Tip Count", title="Distal terminal tips\nby subtree") + 
  facet_wrap(~BINOMIAL, ncol=1)

  figure_aa <- figure_aa + manuscript_theme
  
  ggsave("figures/Figure_Aa.png", width=5, height=20)

  figure_ab <- ggplot(data=full_tree, aes(x=log(FULL_V_TOT), y=log(FULL_TIPS_VOL))) +
  geom_line(aes(x=XS, y=VOL_PREDS), color="green") +
  geom_line(aes(x=XS, y=TIP_PREDS-10), linetype="dashed", color="blue") +
  scale_x_continuous(limits = c(-15, 2)) + scale_y_continuous(limits= c(-15, 0)) +
  labs(x="log Subtree Volume", y="log Total Subtree Tip Volume", title="Distal terminal volume\nby subtree") + 
  facet_wrap(~BINOMIAL, ncol=1)
  
  figure_ab <- figure_ab + manuscript_theme
  
  ggsave("figures/Figure_Ab.png", width=5, height=20)
  return(full_tree)
}

#We are using tree averages of strahler-based space-filling with order 2 excluded to center the length scaling on the theoretical prediction
figure_B <- function()
{
  library(ggridges)
  
  trees <- gather(tree_data, key=DIM, value=EXPONENT, FULL_TIPS, FULL_VOL, LENGTH_PRES, RADIUS_PRES)
  trees$GRP <- NA
  trees$GRP[which(trees$DIM == "FULL_VOL" | trees$DIM == "FULL_TIPS")] <- "Metabolic scaling"
  trees$GRP[which(trees$DIM == "LENGTH_PRES")] <- "Length scaling"
  trees$GRP[which(trees$DIM == "RADIUS_PRES")] <- "Radius scaling"
  trees$GRP <- factor(trees$GRP, levels = c("Radius scaling", "Length scaling", "Metabolic scaling"))
  trees$DIM <- factor(trees$DIM, levels = c("FULL_VOL", "FULL_TIPS", "LENGTH_PRES", "RADIUS_PRES"), 
                                 labels=c("Volume scaling", "Tip scaling", "Length scaling", "Radius scaling"))
  
  #figure_br <- ggplot(trees, aes(x=EXPONENT, y=SCALING, group=SCALING, fill=..x..)) + 
  figure_br <- ggplot(trees, aes(x=EXPONENT, y=DIM,  fill=..x..)) +
  geom_density_ridges_gradient(jittered_points=TRUE, point_size=1, point_alpha = 0.25) + 
  #geom_segment(aes(x=0.75, xend=0.75, y=0, yend=3)) +
  #geom_segment(aes(x=0.70, xend=0.70, y=0, yend=3), linetype="dashed") +
  #geom_segment(aes(x=0.79, xend=0.79, y=0, yend=3), linetype="dashed") +
  scale_fill_gradient(low="purple", high="yellow") + 
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1.0), limits=c(0.25, 1.25)) +  
  facet_wrap(~GRP, scales="free_y", ncol=1) + 
  labs(x="Scaling Exponent", title="Observed scaling exponents across trees")

  figure_br <- figure_br + theme(text = element_text(size = 18), legend.position = "none", axis.title.y = element_blank()) #theme_ridges(center_axis_labels = TRUE)

  ggsave('figures/Figure_B.png')
}

figure_C <- function()
{
  mean_LENGTHs <- data.frame(STRAHLER=numeric(1), MEAN_LENGTH=numeric(1), MEAN_RADIUS=numeric(1), COUNT = numeric(1), FILENAME=character(1), stringsAsFactors = FALSE)
  for(tree in 1:nrow(tree_data))
  {
    current_set <- cylinder_data[which(cylinder_data$FILENAME == tree_data[tree,]$FILENAME),]
    unique_strahler <- unique(current_set$STRAHLER_ORDER)
    for(strahler in unique_strahler)
    {
      if(strahler > 1)
      {
        current_cyls <- which(current_set$STRAHLER_ORDER == strahler & current_set$VOLUME != 0 & !current_set$INVALID)
        next_cyls <- which(current_set$STRAHLER_ORDER == (strahler-1) & current_set$VOLUME != 0 & !current_set$INVALID)
        if(length(current_cyls) > 0 && length(next_cyls) > 0)
        {
            LENGTH_preservation_mean <- mean(current_set$LENGTH[next_cyls]) / mean(current_set$LENGTH[current_cyls])
            RADIUS_preservation_mean <- mean(current_set$RADIUS[next_cyls]) / mean(current_set$RADIUS[current_cyls])
            insertion <- c(strahler, LENGTH_preservation_mean, RADIUS_preservation_mean, length(next_cyls), tree_data[tree,]$FILENAME)
            mean_LENGTHs <- rbind(mean_LENGTHs, insertion)
        }
      }
    }

    file = which(mean_LENGTHs$FILENAME == tree_data[tree,]$FILENAME)
    mean_LENGTHs$STRAHLER <- as.numeric(mean_LENGTHs$STRAHLER)
    mean_LENGTHs$MEAN_LENGTH <- as.numeric(mean_LENGTHs$MEAN_LENGTH)
    mean_LENGTHs$MEAN_RADIUS <- as.numeric(mean_LENGTHs$MEAN_RADIUS)
    mean_LENGTHs$COUNT <- as.numeric(mean_LENGTHs$COUNT)

    tree_data$RADIUS_PRES[tree] <- weighted.mean(mean_LENGTHs$MEAN_RADIUS[file], mean_LENGTHs$COUNT[file])
    tree_data$LENGTH_PRES[tree] <- weighted.mean(mean_LENGTHs$MEAN_LENGTH[file], mean_LENGTHs$COUNT[file])
  }

  lml <- lm(formula = MEAN_LENGTH~STRAHLER, data=mean_LENGTHs, weights = mean_LENGTHs$COUNT)
  print(summary(lml))
  print(confint(lml))
  lmr <- lm(formula = MEAN_RADIUS~STRAHLER, data=mean_LENGTHs, weights = mean_LENGTHs$COUNT)
  print(summary(lmr))
  print(confint(lmr))

  filter <- which(mean_LENGTHs$MEAN_RADIUS < 2) #Filter two outliers for plotting, make sure to keep for analysis
  mean_LENGTHs$STRAHLER <- factor(mean_LENGTHs$STRAHLER, levels=sort(unique(mean_LENGTHs$STRAHLER)))
  gathered <- gather(mean_LENGTHs[filter,], key=DIM, value=MEAN, MEAN_LENGTH, MEAN_RADIUS)
  gathered$DIM <- factor(gathered$DIM, levels=c("MEAN_LENGTH", "MEAN_RADIUS"), 
                                        labels=c("Length", "Radius"))
  figure_gg <- ggplot(gathered, aes(x=STRAHLER, y=MEAN, fill=DIM)) + 
    geom_hline(yintercept=2^(-1/2), linetype="dashed", color="#C77CFF") +
    geom_hline(yintercept=2^(-1/3), linetype="dashed", color="#F8766D") +
    geom_abline(slope=lmr[[1]][2], intercept=lmr[[1]][1], color="#C77CFF") +
    geom_abline(slope=lml[[1]][2], intercept=lml[[1]][1], color="#F8766D") +
    geom_violin(scale="count", alpha=0.5) +   
    scale_fill_manual(values=c("#F8766D", "#C77CFF")) + 
    labs(x="Strahler order", y="Branch scaling ratio", title="Branch scaling by branch level")
  figure_ee <- figure_gg + manuscript_theme + theme(legend.title=element_blank())
  ggsave('figures/Figure_C.png', width=15, height=15)  
  
  return(tree_data)
}

figure_S1 <- function()
{
  #tree_data$STRAHLER_THETA <- -log(2) / log(tree_data$RADIUS_PRES*tree_data$RADIUS_PRES*tree_data$LENGTH_PRES)
  tree_data$STRAHLER_THETA <- -log(2) / log(tree_data$BETA*tree_data$BETA*tree_data$LENGTH_PRES)
  gathered <- gather(tree_data, key=MST, value=EXP, SIMPLE_THETA, STRAHLER_THETA)
  mm <- lm(FULL_VOL~SIMPLE_THETA, data=tree_data)
  print(summary(mm))
  ms <- lm(FULL_VOL~STRAHLER_THETA, data=tree_data)
  print(summary(ms))
  figure_c <- ggplot(gathered, aes(x = EXP, y = FULL_VOL, group=MST, color=MST)) + geom_point() +
          #geom_errorbar(aes(ymin = CI_VOL_MIN, ymax = CI_VOL_MAX)) +
          geom_abline(slope=mm[[1]][2], intercept=mm[[1]][1], color="red") +
          geom_abline(slope=ms[[1]][2], intercept=ms[[1]][1], color="blue") +
          geom_hline(yintercept=0.75, linetype="dashed") + 
          geom_vline(xintercept=0.75, linetype="dashed") +
    #scale_x_continuous(limits=c(0.35, 0.65)) + scale_y_continuous(limits=c(0.4, 1.2)) + 
    labs(x="Theoretical scaling exponent", y="Empirical scaling exponent", title="Theoretical and empirical \n metabolic scaling exponents") 
  figure_c <- figure_c + manuscript_theme + theme(legend.position = c(0.7, 0.2)) 

  ggsave('figures/Figure_S1.png')
}

figure_S2 <- function()
{
  nodes <- tree_data[,c("SYM_VOL", "BETA", "BETA_CI_MIN", "BETA_CI_MAX", "GAMMA", "GAMMA_CI_MIN", "GAMMA_CI_MAX")]
  nodes$REAL <- rep(1, nrow(nodes))
  nodes$NEAR <- rep(0.25, nrow(nodes))
  beta = seq(0.1, 1, 0.01)
  gamma = seq(0.1, 1, 0.01)

  for(x in beta)
  {
    for(y in gamma)
    {
      res <- 2*(x^2)*y
      nearness <- abs(0.7927 - res)
      nodes <- rbind(nodes, c(res, x, x, x, y, y, y, 0.25, nearness))
    }
  }
  nodes<- nodes[seq(dim(nodes)[1],1),]

  figure_d <- ggplot(nodes, aes(x = BETA, y = GAMMA)) + 
    geom_density(aes(x=BETA, y=clamp((..scaled..*0.5)-0.435, 0, 1)), adjust=0.5, fill="red", inherit.aes=FALSE) + 
    geom_density(aes(x=clamp((..scaled..*0.5)-0.445, 0, 1), y=GAMMA), adjust=0.5, fill="darkblue", inherit.aes=FALSE) + 
  geom_errorbar(aes(xmin = BETA_CI_MIN, xmax = BETA_CI_MAX)) +
  geom_errorbar(aes(ymin = GAMMA_CI_MIN, ymax = GAMMA_CI_MAX)) +
  geom_point(aes(color=2*(BETA^2)*GAMMA, alpha=(1-NEAR)*REAL, size=1/(NEAR+0.01)), show.legend = FALSE) +
  
    geom_vline(xintercept=0.707106781, linetype="longdash") + geom_hline(yintercept=0.793700526, linetype="longdash") +
    scale_x_continuous(breaks=seq(0.0, 1, by=0.1)) + scale_y_continuous(breaks=seq(0.0, 1, by=0.1)) +
    labs(x="Radius branching ratio", y="Length branching ratio", title="Average length and radius \n branching ratios")

  figure_d <- figure_d + manuscript_theme

  ggsave('figures/Figure_S1.png')
}

figure_S3 <- function()
{
  nodes <- tree_data[,c("SYM_VOL", "ASYM_VOL", "NETWORK_N", "THETA")]
  nodes$REAL <- rep(1, nrow(nodes))
  nodes$NEAR <- rep(0.05, nrow(nodes))
  vol_scaling = seq(0.1, 1.5, 0.02)
  network_n = seq(2, 100)

  for(x in vol_scaling)
  {
    for(y in network_n)
    {
      res <- asymptotic_formula(x, y) # Returns NaNs for vol_scaling > 0.9
      if(is.na(res))
      {}
      else
      {
        nearness <- abs(0.75 - res)
        nodes <- rbind(nodes, c(0, x, y, res, 0.5, nearness))
      }
    }
  }
  nodes<- nodes[seq(dim(nodes)[1],1),]

  figure_e <- ggplot(nodes, aes(x=SYM_VOL+ASYM_VOL, y=NETWORK_N)) + 
    geom_point(aes(color=THETA, alpha=(1-NEAR)*REAL, size=1/(NEAR+0.01)), show.legend = FALSE) + 
    labs(x="Branching ratio factor (chi)", y="Network depth", title="Network scaling from branching \n geometry and bifurcation depth")

  figure_e <- figure_e + manuscript_theme

  ggsave('figures/Figure_S2.png')
}

figure_S4 <- function()
{
  mm <- sma(GAMMA~L_CV, data=tree_data, method="SMA")
  summary(mm)
  tree_data$PRED_LCV <- NA
  tree_data$PRED_LCV <- mm$coef[[1]][1][1,] + (tree_data$L_CV * mm$coef[[1]][1][2,])
  
  figure_f <- ggplot(tree_data, aes(x=L_CV, y=GAMMA)) + geom_point() + 
    geom_line(aes(x=L_CV, y=PRED_LCV), color="blue") + 
    geom_abline(slope=mm$coef[[1]][1][2,], intercept=mm$coef[[1]][1][1,], linetype="dashed", color="blue") +
    geom_hline(yintercept=0.793700526, linetype="longdash") + 
    scale_y_continuous(limits= c(0, 1.0)) + scale_x_continuous(limits = c(0.0, 2.5)) + 
    labs(x="Terminal tip length coefficient of variation", y = "Mean length scaling ratio", title="Branch length variability and \n whole-tree length scaling")
     
  figure_f <- figure_f + manuscript_theme
  ggsave('figures/Figure_S3.png')
}

faceted_volume_scaling_by_tips <- function(trees, cylinders)
{
  for(tree in 1:nrow(trees))
  {
    current_set <- cylinders[which(cylinders$FILENAME == trees[tree,]$FILENAME),]
    unique_tips <- unique(current_set$FULL_TIPS)
    trees$TIP_COUNT[tree] <- length(unique_tips)
    #trees$TIP_COUNT[tree] <- max(current_set$BRANCH_ORDER, na.rm=TRUE)
  }
  
  gg <- ggplot(trees, aes(x=log(TIP_COUNT), y=FULL_VOL)) + geom_point() + 
    geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
              label.x.npc = "left", label.y.npc = 0.7, formula=y~x, parse = TRUE, size = 5)
  print(gg)
  return(trees)
}

faceted_slenderness_scaling_by_tips <- function(trees, cylinders)
{
  mean_lr <- data.frame(FULL_TIPS=numeric(1), MEAN_LENGTH=numeric(1), MEAN_RADIUS=numeric(1), MEAN_SLENDER=numeric(1), COUNT = numeric(1), SLOPE=numeric(1), PREDICTION=numeric(1), FILENAME=character(1), stringsAsFactors = FALSE)
  cylinders$VAR_FACTOR <- NA
  cylinders$PRED_VAR <- NA
  
  trees$LR_SLOPE <- NA
  for(tree in 1:nrow(trees))
  {
    current_set <- cylinders[which(cylinders$FILENAME == trees[tree,]$FILENAME),]
    unique_tips <- unique(current_set$FULL_TIPS)
    for(tips in unique_tips)
    {
      tip_cyls <- which(current_set$FULL_TIPS == tips & current_set$VOLUME != 0)
      mean_length <- mean(current_set$LENGTH[tip_cyls])
      mean_radius <- mean(current_set$RADIUS[tip_cyls])
      mean_slender <- mean(current_set$LENGTH[tip_cyls] / current_set$RADIUS[tip_cyls])
      insertion <- c(tips, mean_length, mean_radius, mean_slender, length(tip_cyls), NA, NA, trees[tree,]$FILENAME)
      mean_lr <- rbind(mean_lr, insertion)
    }

    file = which(mean_lr$FILENAME == trees[tree,]$FILENAME)
    mean_lr$FULL_TIPS <- as.numeric(mean_lr$FULL_TIPS)
    mean_lr$MEAN_LENGTH <- as.numeric(mean_lr$MEAN_LENGTH)
    mean_lr$MEAN_RADIUS <- as.numeric(mean_lr$MEAN_RADIUS)
    mean_lr$MEAN_SLENDER <- as.numeric(mean_lr$MEAN_SLENDER)
    mean_lr$COUNT <- as.numeric(mean_lr$COUNT)

    #reggie <- lm(formula = log(MEAN_LENGTH)~log(MEAN_RADIUS), data=mean_lr[file,], weights = mean_lr$COUNT[file])
    #mean_lr$PREDICTION[file] <- reggie$coefficients[[1]] + (log(mean_lr$MEAN_RADIUS[file]) * reggie$coefficients[[2]])
    #mean_lr$SLOPE[file] <- reggie$coefficients[[2]]

    reggie <- lm(formula=log(MEAN_SLENDER)~log(MEAN_RADIUS), data=mean_lr[file,], weights = mean_lr$COUNT[file])
    mean_lr$PREDICTION[file] <- reggie$coefficients[[1]] + (log(mean_lr$MEAN_RADIUS[file]) * reggie$coefficients[[2]])
    mean_lr$SLOPE[file] <- reggie$coefficients[[2]]

    trees$LR_SLOPE[tree] <- reggie$coefficients[[2]]
  }
  
  mean_lr$SLOPE <- factor(mean_lr$SLOPE, levels=sort(unique(mean_lr$SLOPE)))
  mean_lr$PREDICTION <- as.numeric(mean_lr$PREDICTION)
  #gg <- ggplot(mean_lr, aes(x=log(MEAN_RADIUS), y=log(MEAN_LENGTH))) + geom_point() + 
  gg <- ggplot(mean_lr, aes(x=log(MEAN_RADIUS), y=log(MEAN_SLENDER))) + geom_point() + 
    geom_line(aes(y=PREDICTION), color="blue") +
    geom_text(x=-2, y=2.5, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(SLOPE))[SLOPE],2))) +
    #geom_text(x=-6, y=2, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(SLOPE))[SLOPE],2))) +
    facet_wrap(~FILENAME)
  print(gg)
  return(trees)
}

faceted_slenderness_scaling <- function(trees, cylinders)
{
  cylinders$VAR_FACTOR <- NA
  cylinders$PRED_VAR <- NA
  cylinders$ELEVATION <- NA
  cylinders$TIP <- NA

  trees$BASE_SLOPE <- NA
  trees$TIP_SLOPE <- NA
  
  for(tree in 1:nrow(trees))
  {
    cyls <- which(cylinders$FILENAME == trees[tree,]$FILENAME)
    tips <- which(cylinders[cyls,]$TIPS == 1)
    negative_inf <- which(log(cylinders$LENGTH[cyls]) == -Inf)
    if(length(negative_inf) > 0)
      cyls <- cyls[-negative_inf]
    reggie <- sma(formula=log(LENGTH/RADIUS)~log(RADIUS), data=cylinders[cyls,][-tips,], method="SMA")
    tip_reggie <- sma(formula=log(LENGTH/RADIUS)~log(RADIUS), data=cylinders[cyls,][tips,], method="SMA")
    
    elev <- lm(formula=log(LENGTH/RADIUS)+(0.33*log(RADIUS))~1, data=cylinders[cyls,])

    cylinders[cyls,][-tips,]$PRED_VAR <- (reggie$coef[[1]][1][1,]) + (log(cylinders$RADIUS[cyls][-tips]) * reggie$coef[[1]][1][2,])
    cylinders[cyls,][tips,]$PRED_VAR <- (tip_reggie$coef[[1]][1][1,]) + (log(cylinders$RADIUS[cyls][tips]) * tip_reggie$coef[[1]][1][2,])
    cylinders[cyls,][-tips,]$TIP <- FALSE
    cylinders[cyls,][tips,]$TIP <- TRUE
    
    trees$BASE_SLOPE[tree] <- reggie$coef[[1]][1][2,]
    trees$TIP_SLOPE[tree] <- tip_reggie$coef[[1]][1][2,]
    #cylinders[cyls,]$VAR_FACTOR <- reggie$coef[[1]][2,1]
    cylinders[cyls,]$ELEVATION <- elev$coef[[1]][1] - (log(cylinders$RADIUS[cyls]) * 0.33)
  }
  #cylinders$VAR_FACTOR <- factor(cylinders$VAR_FACTOR, levels=sort(unique(cylinders$VAR_FACTOR)))
  tree_labeller <- function(variable)
  {return("")}
  v_tips<-ggplot(data=cylinders, aes(y=log(LENGTH/RADIUS), x=log(RADIUS), color=TIP)) +
  geom_point(alpha=0.025) +
  #geom_text(x=-2, y=-2.0, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(VAR_FACTOR))[VAR_FACTOR],2))) +
  #geom_text(x=0.5, y=35.5, color="blue", hjust=0, vjust=1, aes(label=FILENAME)) +
  geom_line(aes(x=log(RADIUS), y=PRED_VAR)) +
  #geom_line(aes(x=log(RADIUS), y=ELEVATION), color="green") +
  xlab("log Branch Radius") + ylab("log Slenderness ratio") + 
  facet_wrap(~FILENAME, labeller=tree_labeller)
  trees$SLOPE_DIFF <- trees$BASE_SLOPE-trees$TIP_SLOPE
  exclude <- which(trees$SLOPE_DIFF > 1)
  gg<-ggplot(trees[-exclude,], aes(x=SLOPE_DIFF, y=EMPIRICAL_VOL))+geom_point()+geom_smooth(method = 'lm', formula=y~x, se=FALSE) + 
     stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), position="dodge",
        formula=y~x, parse = TRUE, size = 8)
  print(gg)
  #ggsave("faceted_length_radius_scaling.png")
  return(trees)
}

faceted_parent_child_volume_scaling <- function(trees, cylinders)
{
  cylinders$VAR_FACTOR <- NA
  cylinders$PRED_VAR <- NA

  trees$VOL_SLOPE <- NA
  for(tree in 1:nrow(trees))
  {
    cyls <- which(cylinders$FILENAME == trees[tree,]$FILENAME & !cylinders$INVALID)
    non_tips <- which(cylinders$TIPS[cyls] > 1)
    negative_inf <- which(log(cylinders$CHILD_VOLS[cyls][non_tips]) == -Inf | log(cylinders$VOLUME[cyls][non_tips]) == -Inf)
    if(length(negative_inf) > 0)
      cyls <- cyls[-negative_inf]
      non_tips <- which(cylinders$TIPS[cyls] > 1)
    if(nrow(cylinders[cyls,][non_tips,]) > 2)
    {
    reggie <- sma(formula=log(CHILD_VOLS)~log(VOLUME), data=cylinders[cyls,][non_tips,], method="SMA")
    
    cylinders[cyls,][non_tips,]$PRED_VAR <- (reggie$coef[[1]][1][1,]) + (log(cylinders$VOLUME[cyls][non_tips]) * reggie$coef[[1]][1][2,])
    cylinders[cyls,]$VAR_FACTOR <- reggie$coef[[1]][1][2,]
    trees$VOL_SLOPE[tree] <- reggie$coef[[1]][1][2,]
    }
  }
  cylinders$VAR_FACTOR <- factor(cylinders$VAR_FACTOR, levels=sort(unique(cylinders$VAR_FACTOR)))
  non_tips <- which(cylinders$TIPS > 1)
  gg <- ggplot(cylinders[non_tips,], aes(x=log(VOLUME), y=log(CHILD_VOLS))) + geom_point(alpha=0.05) + 
    geom_line(aes(x=log(VOLUME), y=PRED_VAR)) +
    geom_text(x=-15, y=15, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(VAR_FACTOR))[VAR_FACTOR],2))) +
    facet_wrap(~FILENAME)
  print(gg)
  return(trees)
}

max_order_comparison <- function(trees, cylinders)
{
  for(tree in 1:nrow(trees))
  {
    cyls <- which(cylinder_data$FILENAME == trees[tree,]$FILENAME)
    trees$MAX_SUCTIONS[tree] <- max(cylinder_data$BRANCH_ORDER[cyls])
    trees$MAX_STRAHLER[tree] <- max(cylinder_data$STRAHLER_ORDER[cyls])  
  }
  ggplot(trees, aes(x=MAX_STRAHLER, y=MAX_SUCTIONS)) + geom_point()
}

faceted_tip_strahler_scaling <- function(trees, cylinders)
{
  cylinders$VAR_FACTOR <- NA
  cylinders$PRED_VAR <- NA

  trees$TIP_SLOPE <- NA
  for(tree in 1:nrow(trees))
  {
    cyls <- which(cylinders$FILENAME == trees[tree,]$FILENAME)
    
    reggie <- sma(formula=log(FULL_TIPS)~STRAHLER_ORDER, data=cylinders[cyls,], method="SMA")
    
    cylinders[cyls,]$PRED_VAR <- (reggie$coef[[1]][1][1,]) + (cylinders$STRAHLER_ORDER[cyls] * reggie$coef[[1]][1][2,])
    cylinders[cyls,]$VAR_FACTOR <- reggie$coef[[1]][1][2,]
    trees$TIP_SLOPE[tree] <- reggie$coef[[1]][1][2,]
  }
  cylinders$VAR_FACTOR <- factor(cylinders$VAR_FACTOR, levels=sort(unique(cylinders$VAR_FACTOR)))
  gg <- ggplot(cylinders, aes(x=STRAHLER_ORDER, y=log(FULL_TIPS))) + geom_point(alpha=0.05) + 
    geom_line(aes(x=STRAHLER_ORDER, y=PRED_VAR)) +
    geom_text(x=2, y=6, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(VAR_FACTOR))[VAR_FACTOR],2))) +
    facet_wrap(~FILENAME)
  print(gg)
  return(trees)
}

#Try the ratio of the sums, instead of the average of the ratios.
faceted_parent_child_volume_scaling_by_tips <- function(trees, cylinders)
{
  mean_vols <- data.frame(FULL_TIPS=numeric(1), MEAN_VOL=numeric(1), COUNT = numeric(1), SLOPE=numeric(1), PREDICTION=numeric(1), FILENAME=character(1), stringsAsFactors = FALSE)
  cylinders <- cylinders[which(cylinders$FULL_TIPS > 1),]

  trees$VOL_PRES <- NA
  for(tree in 1:nrow(trees))
  {
    current_set <- cylinders[which(cylinders$FILENAME == trees[tree,]$FILENAME),]
    unique_tips <- unique(current_set$FULL_TIPS)
    for(tips in unique_tips)
    {
      tip_cyls <- which(current_set$FULL_TIPS == tips & current_set$VOLUME != 0)
      #vol_preservation_mean <- mean(log(current_set$CHILD_VOLS[tip_cyls] / current_set$VOLUME[tip_cyls]))
      vol_preservation_mean <- log(sum(current_set$CHILD_VOLS[tip_cyls] / sum(current_set$VOLUME[tip_cyls])))
      insertion <- c(tips, vol_preservation_mean, length(tip_cyls), NA, NA, trees[tree,]$FILENAME)
      mean_vols <- rbind(mean_vols, insertion)
    }

    file = which(mean_vols$FILENAME == trees[tree,]$FILENAME)
    mean_vols$FULL_TIPS <- as.numeric(mean_vols$FULL_TIPS)
    mean_vols$MEAN_VOL <- as.numeric(mean_vols$MEAN_VOL)
    mean_vols$COUNT <- as.numeric(mean_vols$COUNT)
    
    #reggie <- sma(formula=MEAN_VOL~log(FULL_TIPS), data=mean_vols[file,], method="SMA")    
    #mean_vols$PREDICTION[file] <- reggie$coef[[1]][1][1,] + (log(mean_vols$FULL_TIPS[file]) * reggie$coef[[1]][1][2,])
    #mean_vols$SLOPE[file] <- reggie$coef[[1]][1][2,]
    #trees$VOL_PRES[tree] <- reggie$coef[[1]][1][2,]

    reggie <- lm(formula = MEAN_VOL~log(FULL_TIPS), data=mean_vols[file,], weights = mean_vols$COUNT[file])
    mean_vols$PREDICTION[file] <- reggie$coefficients[[1]] + (log(mean_vols$FULL_TIPS[file]) * reggie$coefficients[[2]])
    mean_vols$SLOPE[file] <- reggie$coefficients[[2]]
    #trees$VOL_PRES[tree] <- reggie$coefficients[[2]]
    #trees$VOL_PRES[tree] <- reggie$coefficients[[1]]
    trees$VOL_PRES[tree] <- mean(mean_vols$MEAN_VOL[file])
  }
  
  mean_vols$SLOPE <- factor(mean_vols$SLOPE, levels=sort(unique(mean_vols$SLOPE)))
  mean_vols$PREDICTION <- as.numeric(mean_vols$PREDICTION)
  gg <- ggplot(mean_vols, aes(x=log(FULL_TIPS), y=MEAN_VOL)) + geom_point(alpha=0.25) + 
    geom_line(aes(x=log(FULL_TIPS), y=PREDICTION), color="blue") +
    geom_hline(yintercept=0, linetype="dashed") + 
    geom_text(x=5, y=-5, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(SLOPE))[SLOPE],2))) +
    facet_wrap(~FILENAME)
  print(gg)
  return(trees)
}

faceted_parent_child_volume_scaling_by_strahler <- function(trees, cylinders)
{
  mean_vols <- data.frame(STRAHLER=numeric(1), MEAN_VOL=numeric(1), COUNT = numeric(1), SLOPE=numeric(1), PREDICTION=numeric(1), FILENAME=character(1), stringsAsFactors = FALSE)
  trees$VOL_PRES <- NA
  
  for(tree in 1:nrow(trees))
  {
    current_set <- cylinders[which(cylinders$FILENAME == trees[tree,]$FILENAME),]
    unique_strahler <- unique(current_set$STRAHLER_ORDER)
    for(strahler in unique_strahler)
    {
      if(strahler != 1)
      {
        current_cyls <- which(current_set$STRAHLER_ORDER == strahler & current_set$VOLUME != 0)
        next_cyls <- which(current_set$STRAHLER_ORDER == strahler-1 & current_set$VOLUME != 0)
        vol_preservation_mean <- sum(current_set$VOLUME[next_cyls]) / sum(current_set$VOLUME[current_cyls])
        insertion <- c(strahler, vol_preservation_mean, length(current_cyls), NA, NA, trees[tree,]$FILENAME)
        mean_vols <- rbind(mean_vols, insertion)
      }
    }

    file = which(mean_vols$FILENAME == trees[tree,]$FILENAME)
    mean_vols$STRAHLER <- as.numeric(mean_vols$STRAHLER)
    mean_vols$MEAN_VOL <- as.numeric(mean_vols$MEAN_VOL)
    mean_vols$COUNT <- as.numeric(mean_vols$COUNT)
    
    #reggie <- sma(formula=MEAN_VOL~log(FULL_TIPS), data=mean_vols[file,], method="SMA")    
    #mean_vols$PREDICTION[file] <- reggie$coef[[1]][1][1,] + (log(mean_vols$FULL_TIPS[file]) * reggie$coef[[1]][1][2,])
    #mean_vols$SLOPE[file] <- reggie$coef[[1]][1][2,]
    #trees$VOL_PRES[tree] <- reggie$coef[[1]][1][2,]

    reggie <- lm(formula = MEAN_VOL~STRAHLER, data=mean_vols[file,], weights = mean_vols$COUNT[file])
    mean_vols$PREDICTION[file] <- reggie$coefficients[[1]] + (mean_vols$STRAHLER[file] * reggie$coefficients[[2]])
    mean_vols$SLOPE[file] <- mean(mean_vols$MEAN_VOL[file])
    #trees$VOL_PRES[tree] <- reggie$coefficients[[2]]
    #trees$VOL_PRES[tree] <- reggie$coefficients[[1]]
    trees$VOL_PRES[tree] <- mean(mean_vols$MEAN_VOL[file])
  }
  
  mean_vols$SLOPE <- factor(mean_vols$SLOPE, levels=sort(unique(mean_vols$SLOPE)))
  mean_vols$PREDICTION <- as.numeric(mean_vols$PREDICTION)
  gg <- ggplot(mean_vols, aes(x=STRAHLER, y=MEAN_VOL)) + geom_point(alpha=0.25) + 
    geom_line(aes(x=STRAHLER, y=PREDICTION), color="blue") +
    geom_hline(yintercept=0.79, linetype="dashed") + 
    geom_text(x=6, y=2, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(SLOPE))[SLOPE],2))) +
    facet_wrap(~FILENAME)
  print(gg)
  return(trees)
}

faceted_parent_child_radius_scaling_by_tips <- function(trees, cylinders)
{
  mean_rads <- data.frame(FULL_TIPS=numeric(1), MEAN_RAD=numeric(1), COUNT = numeric(1), SLOPE=numeric(1), PREDICTION=numeric(1), FILENAME=character(1), stringsAsFactors = FALSE)
  cylinders <- cylinders[which(cylinders$FULL_TIPS > 1),]

  trees$RAD_SLOPE <- NA
  for(tree in 1:nrow(trees))
  {
    current_set <- cylinders[which(cylinders$FILENAME == trees[tree,]$FILENAME),]
    unique_tips <- unique(current_set$FULL_TIPS)
    for(tips in unique_tips)
    {
      tip_cyls <- which(current_set$FULL_TIPS == tips & current_set$VOLUME != 0)
      rad_preservation_mean <- mean(log(current_set$CHILD_RADII[tip_cyls] / current_set$RADIUS[tip_cyls]))
      insertion <- c(tips, rad_preservation_mean, length(tip_cyls), NA, NA, trees[tree,]$FILENAME)
      mean_rads <- rbind(mean_rads, insertion)
    }

    file = which(mean_rads$FILENAME == trees[tree,]$FILENAME)
    mean_rads$FULL_TIPS <- as.numeric(mean_rads$FULL_TIPS)
    mean_rads$MEAN_RAD <- as.numeric(mean_rads$MEAN_RAD)
    mean_rads$COUNT <- as.numeric(mean_rads$COUNT)
    
    reggie <- lm(formula = MEAN_RAD~log(FULL_TIPS), data=mean_rads[file,], weights = mean_rads$COUNT[file])
    mean_rads$PREDICTION[file] <- reggie$coefficients[[1]] + (log(mean_rads$FULL_TIPS[file]) * reggie$coefficients[[2]])
    mean_rads$SLOPE[file] <- reggie$coefficients[[2]]
    trees$RAD_SLOPE[tree] <- reggie$coefficients[[2]]
  }
  mean_rads$SLOPE <- factor(mean_rads$SLOPE, levels=sort(unique(mean_rads$SLOPE)))
  mean_rads$PREDICTION <- as.numeric(mean_rads$PREDICTION)
  gg <- ggplot(mean_rads, aes(x=log(FULL_TIPS), y=MEAN_RAD)) + geom_point(alpha=0.05) + 
    geom_line(aes(x=log(FULL_TIPS), y=PREDICTION), color="blue") +
    geom_text(x=5, y=-1, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(SLOPE))[SLOPE],2))) +
    facet_wrap(~FILENAME)
  print(gg)
  return(trees)
}


faceted_parent_child_length_scaling_by_tips <- function(trees, cylinders)
{
  mean_LENGTHs <- data.frame(FULL_TIPS=numeric(1), MEAN_LENGTH=numeric(1), COUNT = numeric(1), SLOPE=numeric(1), PREDICTION=numeric(1), FILENAME=character(1), stringsAsFactors = FALSE)
  cylinders <- cylinders[which(cylinders$FULL_TIPS > 1),]

  trees$LENGTH_PRES <- NA
  for(tree in 1:nrow(trees))
  {
    current_set <- cylinders[which(cylinders$FILENAME == trees[tree,]$FILENAME),]
    unique_tips <- unique(current_set$FULL_TIPS)
    for(tips in unique_tips)
    {
      tip_cyls <- which(current_set$FULL_TIPS == tips & current_set$VOLUME != 0)
      #LENGTH_preservation_mean <- mean(current_set$CHILD_LENGTHS[tip_cyls] / current_set$LENGTH[tip_cyls])
      LENGTH_preservation_mean <- sum(current_set$CHILD_LENGTHS[tip_cyls]) / sum(current_set$LENGTH[tip_cyls])
      insertion <- c(tips, LENGTH_preservation_mean, length(tip_cyls), NA, NA, trees[tree,]$FILENAME)
      mean_LENGTHs <- rbind(mean_LENGTHs, insertion)
    }

    file = which(mean_LENGTHs$FILENAME == trees[tree,]$FILENAME)
    mean_LENGTHs$FULL_TIPS <- as.numeric(mean_LENGTHs$FULL_TIPS)
    mean_LENGTHs$MEAN_LENGTH <- as.numeric(mean_LENGTHs$MEAN_LENGTH)
    mean_LENGTHs$COUNT <- as.numeric(mean_LENGTHs$COUNT)
    
    reggie <- lm(formula = MEAN_LENGTH~log(FULL_TIPS), data=mean_LENGTHs[file,], weights = mean_LENGTHs$COUNT[file])
    mean_LENGTHs$PREDICTION[file] <- reggie$coefficients[[1]] + (log(mean_LENGTHs$FULL_TIPS[file]) * reggie$coefficients[[2]])
    mean_LENGTHs$SLOPE[file] <- reggie$coefficients[[2]]
    
    #trees$LENGTH_SLOPE[tree] <- reggie$coefficients[[2]]
    #trees$LENGTH_PRES[tree] <- reggie$coefficients[[1]]
    trees$LENGTH_PRES[tree] <- mean(mean_LENGTHs$MEAN_LENGTH[file])
  }
  mean_LENGTHs$SLOPE <- factor(mean_LENGTHs$SLOPE, levels=sort(unique(mean_LENGTHs$SLOPE)))
  mean_LENGTHs$PREDICTION <- as.numeric(mean_LENGTHs$PREDICTION)
  gg <- ggplot(mean_LENGTHs, aes(x=log(FULL_TIPS), y=MEAN_LENGTH)) + geom_point(alpha=0.25) + 
    geom_line(aes(x=log(FULL_TIPS), y=PREDICTION), color="blue") +
    geom_hline(yintercept=0) +
    geom_text(x=5, y=-1, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(SLOPE))[SLOPE],2))) +
    facet_wrap(~FILENAME)
  print(gg)
  return(trees)
}

#Average ratio between levels is near, but less than 1. Doesn't disclude 1 or 0.79, wide CIs. Perhaps its time we unfiltered gammas.
#If we remove 3 trees with length preservation averages greater than 3, we get mean of 0.77, not signif diff from 0.79.
#I mean... this is straightforwardly a measure of gamma that fits predictions. Seems worth sharing.
##Bugged code - and you shared it with colleagues. shape up
faceted_parent_child_length_scaling_by_strahler <- function(trees, cylinders)
{
  mean_LENGTHs <- data.frame(STRAHLER=numeric(1), MEAN_LENGTH=numeric(1), COUNT = numeric(1), SLOPE=numeric(1), PREDICTION=numeric(1), FILENAME=character(1), stringsAsFactors = FALSE)
  #cylinders <- cylinders[which(cylinders$FULL_TIPS > 1),]

  trees$LENGTH_PRES <- NA
  for(tree in 1:nrow(trees))
  {
    current_set <- cylinders[which(cylinders$FILENAME == trees[tree,]$FILENAME),]
    unique_strahler <- unique(current_set$STRAHLER_ORDER)
    for(strahler in unique_strahler)
    {
      if(strahler != 1)
      {
        current_cyls <- which(current_set$STRAHLER_ORDER == strahler & current_set$VOLUME != 0)
        next_cyls <- which(current_set$STRAHLER_ORDER == (strahler-1) & current_set$VOLUME != 0)
        LENGTH_preservation_mean <- sum(current_set$LENGTH[current_cyls]) / sum(current_set$LENGTH[next_cyls])
        insertion <- c(strahler, LENGTH_preservation_mean, length(current_cyls), NA, NA, trees[tree,]$FILENAME)
        mean_LENGTHs <- rbind(mean_LENGTHs, insertion)
      }
    }

    file = which(mean_LENGTHs$FILENAME == trees[tree,]$FILENAME)
    mean_LENGTHs$STRAHLER <- as.numeric(mean_LENGTHs$STRAHLER)
    mean_LENGTHs$MEAN_LENGTH <- as.numeric(mean_LENGTHs$MEAN_LENGTH)
    mean_LENGTHs$COUNT <- as.numeric(mean_LENGTHs$COUNT)
    
    reggie <- lm(formula = MEAN_LENGTH~STRAHLER, data=mean_LENGTHs[file,], weights = mean_LENGTHs$COUNT[file])
    mean_LENGTHs$PREDICTION[file] <- reggie$coefficients[[1]] + (mean_LENGTHs$STRAHLER[file] * reggie$coefficients[[2]])
    mean_LENGTHs$SLOPE[file] <- reggie$coefficients[[2]]
    trees$LENGTH_PRES[tree] <- mean(mean_LENGTHs$MEAN_LENGTH[file])
  }
  mean_LENGTHs$SLOPE <- factor(mean_LENGTHs$SLOPE, levels=sort(unique(mean_LENGTHs$SLOPE)))
  mean_LENGTHs$PREDICTION <- as.numeric(mean_LENGTHs$PREDICTION)
  gg <- ggplot(mean_LENGTHs, aes(x=STRAHLER, y=MEAN_LENGTH)) + geom_point(alpha=0.25) + 
    geom_line(aes(x=STRAHLER, y=PREDICTION), color="blue") +
    geom_hline(yintercept=1) + 
    geom_text(x=6, y=3, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(SLOPE))[SLOPE],2))) +
    facet_wrap(~FILENAME)
  #print(gg)
  return(trees)
}

get_better_gammas <- function(trees, cylinders)
{
  trees$BETTER_GAMMA <- NA
  for(tree in 1:nrow(trees))
  {
    current_set <- cylinders[which(cylinders$FILENAME == trees[tree,]$FILENAME),]
    trees[tree,]$BETTER_GAMMA <- mean(current_set$GAMMA, na.rm = TRUE)
  }
  return(trees)
}

vol_rad_preservation <- function()
{
   ggplot(tt, aes(x=RAD_SLOPE, y=VOL_SLOPE)) + geom_point() + geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
              label.x.npc = "left", label.y.npc = 0.7, formula=y~x, parse = TRUE, size = 5)
}

vol_len_preservation <- function()
{
  ggplot(tt, aes(x=LENGTH_SLOPE, y=VOL_SLOPE)) + geom_point() + geom_smooth(method = 'lm', formula=y~x, se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
              label.x.npc = "left", label.y.npc = 0.7, formula=y~x, parse = TRUE, size = 5) 
}

faceted_full_volume_scaling <- function(trees, cylinders)
{
  cylinders$VAR_FACTOR <- NA
  cylinders$PRED_VAR <- NA

  trees$FULL_VOL <- NA
  for(tree in 1:nrow(trees))
  {
    cyls <- which(cylinders$FILENAME == trees[tree,]$FILENAME)
    non_tips <- which(cylinders$TIPS[cyls] > 1)
    
    reggie <- sma(formula=log(FULL_TIPS_VOL)~log(FULL_V_TOT), data=cylinders[cyls,][non_tips,], method="SMA")
    cylinders[cyls,][non_tips,]$PRED_VAR <- (reggie$coef[[1]][1][1,]) + (log(cylinders$FULL_V_TOT[cyls][non_tips]) * reggie$coef[[1]][1][2,])
    cylinders[cyls,]$VAR_FACTOR <- reggie$coef[[1]][1][2,]
    ci_vol = reggie$slopetest[[1]]$ci
    trees$CI_VOL_MIN[tree] <- ci_vol[1]
    trees$CI_VOL_MAX[tree] <- ci_vol[2]
    trees$FULL_VOL[tree] <- reggie$coef[[1]][1][2,]
  }
  cylinders$VAR_FACTOR <- factor(cylinders$VAR_FACTOR, levels=sort(unique(cylinders$VAR_FACTOR)))
  non_tips <- which(cylinders$TIPS > 1)
  gg <- ggplot(cylinders[non_tips,], aes(x=log(FULL_V_TOT), y=log(FULL_TIPS_VOL))) + geom_point(alpha=0.05) + 
    geom_line(aes(x=log(FULL_V_TOT), y=PRED_VAR)) +
    geom_text(x=-12, y=-2, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(VAR_FACTOR))[VAR_FACTOR],2))) +
    facet_wrap(~FILENAME)
  print(gg)
  return(trees)
}

faceted_full_volume_scaling_by_tips <- function(trees, cylinders)
{
  mean_vols <- data.frame(FULL_TIPS=numeric(1), MEAN_V_TOT=numeric(1), MEAN_TIPS_VOL=numeric(1), COUNT = numeric(1), SLOPE=numeric(1), PREDICTION=numeric(1), FILENAME=character(1), stringsAsFactors = FALSE)
  cylinders <- cylinders[which(cylinders$FULL_TIPS > 1),]

  trees$FULL_VOL_SLOPE <- NA
  for(tree in 1:nrow(trees))
  {
    current_set <- cylinders[which(cylinders$FILENAME == trees[tree,]$FILENAME),]
    unique_tips <- unique(current_set$FULL_TIPS)
    for(tips in unique_tips)
    {
      tip_cyls <- which(current_set$FULL_TIPS == tips & current_set$VOLUME != 0)
      if(length(tip_cyls) > 0)
      {
        full_v_tot <- sum(current_set$FULL_V_TOT[tip_cyls], na.rm=TRUE)  #Sum these instead
        full_tips_vol <- sum(current_set$FULL_TIPS_VOL[tip_cyls], na.rm=TRUE)  #
        insertion <- c(tips, full_v_tot, full_tips_vol, length(tip_cyls), NA, NA, trees[tree,]$FILENAME)
        mean_vols <- rbind(mean_vols, insertion)
      }
    }

    file = which(mean_vols$FILENAME == trees[tree,]$FILENAME)
    mean_vols$FULL_TIPS <- as.numeric(mean_vols$FULL_TIPS)
    mean_vols$MEAN_V_TOT <- as.numeric(mean_vols$MEAN_V_TOT)
    mean_vols$MEAN_TIPS_VOL <- as.numeric(mean_vols$MEAN_TIPS_VOL)
    mean_vols$COUNT <- as.numeric(mean_vols$COUNT)
    
    reggie <- sma(formula = log(MEAN_TIPS_VOL)~log(MEAN_V_TOT), data=mean_vols[file,], method="SMA")
    mean_vols$PREDICTION[file] <- reggie$coef[[1]][1][1,] + (log(mean_vols$MEAN_V_TOT[file]) * reggie$coef[[1]][1][2,])
    mean_vols$SLOPE[file] <- reggie$coef[[1]][1][2,]
    trees$FULL_VOL_SLOPE[tree] <- reggie$coef[[1]][1][2,]
  }
  
  mean_vols$SLOPE <- factor(mean_vols$SLOPE, levels=sort(unique(mean_vols$SLOPE)))
  mean_vols$PREDICTION <- as.numeric(mean_vols$PREDICTION)

  gg <- ggplot(mean_vols, aes(x=log(MEAN_V_TOT), y=log(MEAN_TIPS_VOL))) + geom_point(alpha=0.25) + 
    geom_line(aes(x=log(MEAN_V_TOT), y=PREDICTION), color="blue") +
    geom_text(x=-8, y=-2, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(SLOPE))[SLOPE],2))) +
    facet_wrap(~FILENAME)
  print(gg)
  return(trees)
}

faceted_full_tip_scaling <- function(trees, cylinders)
{
  cylinders$VAR_FACTOR <- NA
  cylinders$PRED_VAR <- NA

  trees$FULL_TIPS <- NA
  for(tree in 1:nrow(trees))
  {
    cyls <- which(cylinders$FILENAME == trees[tree,]$FILENAME)
    non_tips <- which(cylinders$TIPS[cyls] > 1)
    
    if(nrow(cylinders[cyls,][non_tips,]) > 2)
    {
      reggie <- sma(formula=log(FULL_TIPS)~log(FULL_V_TOT), data=cylinders[cyls,][non_tips,], method="SMA")
      cylinders[cyls,][non_tips,]$PRED_VAR <- (reggie$coef[[1]][1][1,]) + (log(cylinders$FULL_V_TOT[cyls][non_tips]) * reggie$coef[[1]][1][2,])
      cylinders[cyls,]$VAR_FACTOR <- reggie$coef[[1]][1][2,]
      trees$FULL_TIPS[tree] <- reggie$coef[[1]][1][2,]
    }
  }
  cylinders$VAR_FACTOR <- factor(cylinders$VAR_FACTOR, levels=sort(unique(cylinders$VAR_FACTOR)))
  non_tips <- which(cylinders$TIPS > 1)
  gg <- ggplot(cylinders[non_tips,], aes(x=log(FULL_V_TOT), y=log(FULL_TIPS))) + geom_point(alpha=0.05) + 
    geom_line(aes(x=log(FULL_V_TOT), y=PRED_VAR)) +
    geom_text(x=-12, y=5, color="red", hjust=0, vjust=1, aes(label=round(as.numeric(levels(VAR_FACTOR))[VAR_FACTOR],2))) +
    facet_wrap(~FILENAME)
  print(gg)
  return(trees)
}

###
# All-node plotting functions
###
#Make sure to feed valid cylinders to these functions for interpretation

plot_all_volume_scaling_slopes_raw <- function(trees, cylinders)
{
  ggplot(cylinders, aes(x=log(RAW_V_TOT), y=RAW_PRED_VOL)) + geom_line(aes(group=FILENAME, color=FILENAME), alpha=0.5) + 
    geom_abline(slope=0.75, intercept=-2.5, linetype="dashed")
}


plot_all_volume_scaling_slopes_trees <- function(trees, cylinders)
{
  ggplot(cylinders, aes(x=RAW_V_TOT, y=RAW_PRED_VOL)) + geom_blank() + geom_abline(data=trees, aes(slope=EMPIRICAL_VOL, intercept=EMPIRICAL_INTERCEPT-5))
}

plot_three_parameter_theta <- function(trees)
{
  nodes <- trees[,c("BETA", "GAMMA", "THETA")]
  nodes$B_R <- rep(2, nrow(nodes))
  nodes$REAL <- rep(1, nrow(nodes))
  nodes$NEAR <- rep(0.05, nrow(nodes))
  beta = seq(0.1, 1, 0.05)
  gamma = seq(0.1, 1, 0.05)
  b_r = seq(1, 3, 0.1)

  for(x in beta)
  {
    for(y in gamma)
    {
      for(z in b_r)
      {
        res <- -log(z) / log((x^2) * y)
        nearness <- abs(0.75 - res)
        nodes <- rbind(nodes, c(x, y, res, z, 0.01, nearness))
      }
    }
  }
  nodes<- nodes[seq(dim(nodes)[1],1),]
  theta = 225
  phi = 0

  ggplot(nodes, aes(x=BETA, y=GAMMA, z=B_R, color=THETA, alpha = (1-NEAR)*REAL, size = 2/(NEAR+0.01))) + 
    #geom_point(aes(color=THETA, alpha=(1-NEAR)*REAL, size=1/(NEAR+0.01))) +
    axes_3D(theta=theta, phi=phi) + stat_3D(theta=theta, phi=phi) + 
    labs_3D(theta=theta, phi=phi, labs=c("BETA", "GAMMA", "BRANCHING RATIO")) + 
    theme_void() + theme(legend.position = "none") 
    #scale_color_gradient(low="purple", high="yellow")
    #scale_x_continuous(breaks=seq(0.1, 1.5, by=0.1)) + scale_y_continuous(breaks=seq(2, 50, by=2))
}

plot_asymptotic_formula <- function (trees)
{
	nodes <- trees[,c("SYM_VOL", "ASYM_VOL", "NETWORK_N", "THETA")]
	nodes$REAL <- rep(1, nrow(nodes))
	nodes$NEAR <- rep(0.05, nrow(nodes))
	vol_scaling = seq(0.1, 1.5, 0.01)
	network_n = seq(2, 100)

	for(x in vol_scaling)
	{
		for(y in network_n)
		{
			res <- asymptotic_formula(x, y) # Returns NaNs for vol_scaling > 0.9
			if(is.na(res))
			{}
			else
			{
				nearness <- abs(0.75 - res)
				nodes <- rbind(nodes, c(0, x, y, res, 0.5, nearness))
			}
		}
	}
	nodes<- nodes[seq(dim(nodes)[1],1),]

	ggplot(nodes, aes(x=SYM_VOL+ASYM_VOL, y=NETWORK_N)) + 
    geom_point(aes(color=THETA, alpha=(1-NEAR)*REAL, size=1/(NEAR+0.01))) #+
		#scale_color_gradient(low="purple", high="yellow")
		#scale_x_continuous(breaks=seq(0.1, 1.5, by=0.1)) + scale_y_continuous(breaks=seq(2, 50, by=2))
}

plot_symmetrical_volume_scaling <- function (trees)
{
	nodes <- trees[,c("SYM_VOL", "BETA", "GAMMA")]
	nodes$REAL <- rep(1, nrow(nodes))
	nodes$NEAR <- rep(0.05, nrow(nodes))
	beta = seq(0.1, 1, 0.025)
	gamma = seq(0.1, 1, 0.025)

	for(x in beta)
	{
		for(y in gamma)
		{
			res <- 2*(x^2)*y
			nearness <- abs(0.7927 - res)
			nodes <- rbind(nodes, c(res, x, y, 0.5, nearness))
		}
	}
	nodes<- nodes[seq(dim(nodes)[1],1),]

	ggplot(nodes, aes(x=BETA, y=GAMMA)) + 
  geom_point(aes(color=SYM_VOL, alpha=(1-NEAR)*REAL, size=1/(NEAR+0.01))) +
    geom_vline(xintercept=0.707106781) + geom_hline(yintercept=0.793700526) +
		#scale_color_gradient(low="purple", high="yellow")
		scale_x_continuous(breaks=seq(0.1, 1, by=0.1)) + scale_y_continuous(breaks=seq(0.1, 1, by=0.1))
}

plot_asymmetrical_volume_scaling <- function (trees)
{
	nodes <- trees[,c("D_BETA", "D_GAMMA")]
	nodes$ASYM_VOL <- 4*0.7*nodes$D_BETA*nodes$D_GAMMA + 2*0.79*(nodes$D_BETA^2)	
	nodes$REAL <- rep(1, nrow(nodes))
	nodes$NEAR <- abs(nodes$ASYM_VOL)
	d_beta = seq(-0.3, 0.3, 0.005)
	d_gamma = seq(0.0, 0.3, 0.005)

	for(x in d_beta)
	{
		for(y in d_gamma)
		{
			res <- 4*0.7*x*y + 2*0.79*(x^2)
			nearness <- abs(res)
			nodes <- rbind(nodes, c(x, y, res, 0.55, nearness))
		}
	}
	nodes<- nodes[seq(dim(nodes)[1],1),]

	ggplot(nodes, aes(x=D_BETA, y=D_GAMMA)) + geom_point(aes(color=ASYM_VOL, alpha=(1-NEAR)*REAL, size=1/(NEAR+0.05))) +
		#scale_color_gradient(low="purple", high="yellow")
		scale_x_continuous(breaks=seq(-0.5, 0.5, by=0.1)) + scale_y_continuous(breaks=seq(-0.5, 0.5, by=0.1)) +
		scale_color_gradient(low="purple", high="yellow")
}

plot_trait_pca <- function(trees)
{
  library(ggfortify)
  architectural_traits <- c("BETA", "GAMMA", "EMPIRICAL_VOL")
  pca_traits <- c("leaf area", "leaf area per leaf dry mass", "leaf nitrogen content per leaf area", "leaf life span", "maximum whole plant height", "seed mass", "stem wood density")
  pca_data <- trees[,c(architectural_traits, pca_traits)]
  agg_tree_data <- aggregate(pca_data, by=list(trees$BINOMIAL), FUN=mean)
  pcares <- prcomp(agg_tree_data[,-1], scale. = TRUE)
  autoplot(pcares, data=agg_tree_data, color='Group.1', loadings = TRUE, loadings.label = TRUE)
}

plot_vegan_trait_significance <- function(trees)
{
  library(BiodiversityR)
  agg_tree_data <- aggregate(trees, by=list(trees$BINOMIAL), FUN=mean)
  pca_model <- rda(agg_tree_data[,pca_trait_list()], scale = TRUE)
  biodivr <- PCAsignificance(pca_model)

  #June 1st - the first six axes are signficant
  #Axis one appears to be an MST axis - architectural traits and leaf photosynthetic rate
  #2 - LES
  #3 - Seed mass and wood density
  #4 - Size 
  png(filename='figures/pc12.png', width=800, height=800)
  plot1 <- ordiplot(pca_model, display="species",choices = c(1,2), scaling = 1)
  ordiequilibriumcircle(pca_model,plot1)
  text(plot1, "species", col="blue", cex=0.75)
  dev.off()

  png(filename='figures/pc34.png', width=800, height=800)
  plot2 <- ordiplot(pca_model, display="species",choices = c(3,4), scaling = 1)
  ordiequilibriumcircle(pca_model,plot2)
  text(plot2, "species", col="blue", cex=0.9)
  dev.off()

  png(filename='figures/pc56.png', width=800, height=800)
  plot3 <- ordiplot(pca_model, display="species",choices = c(5,6), scaling = 1)
  ordiequilibriumcircle(pca_model,plot3)
  text(plot3, "species", col="blue", cex=0.9)
  dev.off()

  return(biodivr)
}

#Takes an internal representation and plots the tree as a network
#Takes the name of a value to color nodes in the tree
plot_tree_network <- function(tree, node_val)
{
  #tree[which(tree$CHILD_IDS == "NA"),]$CHILD_IDS <- "_" 
  v1 <- c()
  v2 <- c()
  node_value <- c()
  vol <- c()
  
  for(i in seq(1, nrow(tree)))
  {
    children <- unlist(strsplit(tree[i,]$CHILD_IDS, split="_"))
    if(length(children) > 0)
    {
      for(x in seq(1, length(children)))
      {
        child <- as.numeric(children[x])
        if(!is.na(child))
        {
          v1 <- c(v1, child)
          v2 <- c(v2, i)
          node_value <- c(node_value, tree[child,node_val])
          vol <- c(vol, tree[child,]$VOLUME)
        }
      } 
    }
  }

  edgelist <- matrix(1, nrow=length(v1), ncol=4)
  edgelist[,1] = as.numeric(v1)
  edgelist[,2] = as.numeric(v2) 
  edgelist[,3] = as.numeric(node_value)
  edgelist[,4] = as.numeric(vol)

  #adjacency <- matrix(0, nrow=length(v1)+1, ncol=length(v1)+1)
  #for(x in seq(1, length(v1)))
  #{
  #  adjacency[v1[x], v2[x]] = 1
  #}
  
  n <- network(edgelist[,1:2])
  set.edge.attribute(n, "node_value", edgelist[,3])
  set.edge.attribute(n, "volume", edgelist[,4])
  nn <- ggnetwork(n, layout = 'kamadakawai', niter=1000, weights = "length")
  ggplot(nn, aes(x = x, y = y, xend = xend, yend = yend)) +
         #geom_edges() +
         geom_nodes(aes(color=node_value, size=volume)) +
         theme_blank()
}

example_network <- function()
{
  n <- network(rgraph(100, tprob = 0.2), directed = FALSE)
  n %v% "family" <- sample(letters[1:3], 100, replace = TRUE)
  n %v% "importance" <- sample(1:3, 100, replace = TRUE)
  e <- network.edgecount(n)
  set.edge.attribute(n, "type", sample(letters[24:26], e, replace = TRUE))
  set.edge.attribute(n, "day", sample(1:3, e, replace = TRUE))
  ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(linetype = type), color = "grey50") +
    geom_nodes(aes(color = family, size = importance)) +
    theme_blank()
}