#The string-based approach doesn't mesh with out previous work, or the data structures
#R excels at - we need a dataframe or table-based implementation

#Growing frames by row will be very slow, but our simulated L-Systems will never be 
#that large, ideally. Simulations may be intractable around 12-15 branching generations

#One problem is that rows (observations) in QSMs are always geometric units, while
#symbols in L-Systems can represent more abstract developmental potential. R just
#seems like the wrong tool for this. Maybe L-Systems too.
#We need to plan more what this needs to be able to do.

source("l_system.R")
library(stringr)

#A - apex
#B - branch
#C - second-order branch
#K - inflorescence
#I - internode
#[ - branch open
#] - branch close
alphabet <- c('A', 'B', 'C', 'K', 'I', '[', ']')

#Models taken from (Prusinkiewicz and Remphrey 2000)
#Worth double-checking and ground-truthing with (Halle et al 1978)
#There appear to be a few errors

corner_key <- c('I[K]A', 'B', 'C', 'K', 'I', '[', ']')
corner <- create_system(alphabet, 'A', corner_key)

chamberlain_key <- c('I[A]K', 'B', 'C', 'K', 'I', '[', ']')
chamberlain <- create_system(alphabet, 'A', chamberlain_key)

roux_cook_key <- c('I[B]A', 'I[K]B', 'C', 'K', 'I', '[', ']')
roux_cook <- create_system(alphabet, 'A', roux_cook_key)

build_massart <- function(n)
{
	branch_whorl <- '[B]'
	apical <- c(paste(paste('I', str_dup(branch_whorl, n)), 'A'))
	lateral<- c('I[C]B')
	second <- c('I[K]C')
	return (c(apical, lateral, second, 'K', 'I', '[', ']'))
}
massart_key <- build_massart(5)
massart <- create_system(alphabet, 'A', massart_key)

#Probabilistic systems can model variation in differentiation, but also leads
#to unifying systems. Also differences in branching (i.e. corner/chamberlain)
#may not affect modular scaling in the end.
rauh_attim_key <- c('I[A]A', 'B', 'C', 'K', 'I', '[', ']')
rauh_attim <- create_system(alphabet, 'A', rauh_attim_key)

shoute_key <- c('I[A][A]', 'B', 'C', 'K', 'I', '[', ']')
shoute <- create_system(alphabet, 'A', shoute_key)

leeuwenberg_key <- c('I[A][A]K', 'B', 'C', 'K', 'I', '[', ']')
leeuwenberg <- create_system(alphabet, 'A', leeuwenberg_key)

#Stone and Scarrone differ in plagiotropic orientation
#We cannot analyze this kind of variation with L-Systems
petit_key <- c('I[B]K', 'I[B][B]K', 'C', 'K', 'I', '[', ']')
petit <- create_system(alphabet, 'A', petit_key)

#Again, differences in branching (here compared to massart) won't change observed
#modular scaling, but there is much biology cloaked here without data appended
#to our modules (parametric systems)
build_fagerlind <- function(n)
{
	branch_whorl <- '[B]'
	apical <- c(paste(paste('I', str_dup(branch_whorl, n)), 'A'))
	lateral<- c('I[B][B]C')
	second <- 'IK'
	return (c(apical, lateral, second, 'K', 'I', '[', ']'))
}
fagerlind_key <- build_fagerlind(5)
fagerlind <- create_system(alphabet, 'A', fagerlind_key)

build_aubreville <- function(n)
{
	branch_whorl <- '[B]'
	apical <- c(paste(paste('I', str_dup(branch_whorl, n)), 'A'))
	lateral<- c('I[B][B]C')
	second <- 'I[C]K'
	return (c(apical, lateral, second, 'K', 'I', '[', ']'))
}
aubreville_key <- build_aubreville(5)
aubreville <- create_system(alphabet, 'A', aubreville_key)

build_koriba <- function(n)
{
	branch_whorl <- '[B]'
	apical <- c(paste(paste('I[A]', str_dup(branch_whorl, n)), 'K'))
	lateral<- 'B'
	second <- 'C'
	return (c(apical, lateral, second, 'K', 'I', '[', ']'))
}
koriba_key <- build_koriba(5)
koriba <- create_system(alphabet, 'A', koriba_key)

build_prevost <- function(n)
{
	branch_whorl <- '[B]'
	apical <- c(paste(paste('I[A]I', str_dup(branch_whorl, n)), 'K'))
	lateral<- 'I[B][B]K'
	second <- 'C'
	return (c(apical, lateral, second, 'K', 'I', '[', ']'))
}
prevost_key <- build_prevost(5)
prevost <- create_system(alphabet, 'A', prevost_key)

build_nozeran <- function(n)
{       
	branch_whorl <- '[B]'
	apical <- c(paste(paste('I[A]I', str_dup(branch_whorl, n)), 'K'))
	lateral<- 'I[C][C]'
	second <- 'I[K]B'
	return (c(apical, lateral, second, 'K', 'I', '[', ']'))
}
nozeran_key <- build_nozeran(5)
nozeran <- create_system(alphabet, 'A', nozeran_key)

#Here we aren't exercising the full power of L-Systems, as orthotropy and 
#plagiotropy are probably valid functional types of branches, and the symbolic
#mode can model this transition (As this model is intended to do)
#In this attempt, we just want to count 'biomass' modules
mangenot_key <- c('I[A]B', 'I[KC]B', 'I[KB]C', 'K', 'I', '[', ']')
mangenot <- create_system(alphabet, 'A', mangenot_key)

#Should be identical to mangenot given we aren't recognizing branch differentiation
#Authors do not provide information on flowering
champagnat_key <- c('I[A]B', 'I[C]B', 'IB', 'K', 'I', '[', ']')
champagnat <- create_system(alphabet, 'A', champagnat_key)

#Same with plagiotropy
troll_key <- c('I[A]B', 'I[C]B', 'IB', 'K', 'I', '[', ']')
troll <- create_system(alphabet, 'A', troll_key)

#Rhythmic and continuous growth are un-examined, and ultimately related to the way
#we differentially apply productions to model growth and differentiation across time

build_wbe <- function(n)
{
	branch_whorl <- '[IA]'
	apical <- c(str_dup(branch_whorl, n))
	return(c(apical, 'B', 'C', 'K', 'I', '[', ']'))
}
wbe_key <- build_wbe(2)
wbe <- create_system(alphabet, 'A', wbe_key)

simple_key <- c('IA', 'B', 'C', 'K', 'I', '[', ']')
simple <- create_system(alphabet, 'A', simple_key)
 
