source("asymmetric_analysis.R")

#For all tree types we want to try to compute the number of rows
#(branches) before starting generation.

length_ratio_def = 0.3
length_ratio_sim = seq(0.1, 1, by=0.1)

radius ratio_def = 0.5
radius_ratio_sim = seq(0.25, 0.75, by=0.05)

length_asym_ratio_def = 0
length_asym_ratio_sim = seq(-0.5, 0.5, 0.1)

radius_asym_ratio_def = 0
radius_asym_ratio_sim = seq(-0.5, 0.5, 0.1)

root_length_def = 1000
root_length_sim = seq(100, 1000, by=100) 

root_width_def = 100
root_width_sim = seq(10, 100, by=10)

leaf_length_def = 10
leaf_length_sim = seq(0.1, 10, by=0.1)

leaf_width_def = 1
leaf_width_sim = seq(0.01, 1, by=0.01)

epiphyte_length_def = 1
epiphyte_lingth_sim = seq(1, 100, by=1)

epiphyte_radius_def = 0.1
epiphye_radius_sim = seq(0.1, 10, by = 0.1)

#Nest/food distributions? Dictionary of probabilities keyed by order?
#Epiphytic connection distributions?
#We have no clue which distribution to use.
#If we had one (say a log normal one) how would we apply it? Code that.

#...parameterize sample function, dbinom function?
#for i=1..n branch orders:
#	sample(c(0, 1), 2^i, replace=TRUE, prob=probs[i])
#	or
#	dbinom(...)
#	Or explicitly generating probabilistic elements isn't
#	necessary, and only sample probabilities simulating ants

#When asking, Does this branch have food or not? Or,
#Does this branch have nesting space, or an epiphyte?
#We can use a binomial distribution with probabilities
#defined by branch order. This is a decent first approximation.

#We can generate the probabilities themselves from another random
#distribution defined by the size of the tree and other parameters
#that skew things as needed. What is that underlying distribution?

#trunk, leaf are lists: list("L"=length, "R"=radius)
build_symmetric_tree(trunk, leaf)
{
	
}

build_asymmetric_tree(trunk, leaf, l_ratio, r_ratio, asym_l, asym_r)
{

}

build_epiphytic_forest()
{

}
