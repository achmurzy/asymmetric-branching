Three projects here:
-QSM analysis
-Tree foraging simulation - arboreal communities
-L-Systems

Instructions:
-Run 'workflow.R' to pull all of the fit QSMs and perform a complete scaling analysis. Filenames are unique tree identifiers
for the entire workflow - data from other sources is never modified, filenames are used exactly.
-Testing/debugging - use the single test tree 'Ery_01' with guide at the top of workflow.R

-Report: use library(rmarkdown) then: 
render('scaling_report.Rmd')

TO-DO:	
-Add stats to every analysis included in your Tuesday seminar talk. Confidence intervals on regression data points and
slopes you're drawing conclusions from.
          
-Excluding nodes/branches based on scaling consideratinos might be affecting non-scaling related computations. 
Make sure this isn't the case, especially for path fraction - 
mainly check how your INVALID column is affecting calculations, reducing tips, etc.

-Think about how to modularize nodes and collections of nodes. 
  -Do a within-tree clustering analysis on raw geometry and scaling ratios to classify.
  	-Better, some kind of network-based (hierarchical?) clustering to take account of connectivity between branches/scaling ratios to detect reiterated units.

-Port the analysis to treestruct. Provide an automated way to feed columns (and desired operations) to the  tree recursion.
-Integrate with L-Systems library (applying transformation rules to branch/cylinder data rows)

DONE/COMMENTS:
-Compare symmetric and asymmetric values for volume scaling. Relative magnitudes?
	SOLUTION: Symmetric values greatly outweigh asymmetric contributions by 5-6 times. Plot the asymptotic equation to see how this 'nudges' scaling values for a given network N.

-Create a faceting function for comparing distributions of traits across nodes, for large numbers of trees. 'Small multiples'
	-Need to visualize regression lines not point clouds in small multiples, and sort by the magnitude of the slopes to see influences on global averages of scaling values

Negative asymmetry emerges as pervasive in trees. This might be a more important result than we realize.

-Make sure 0-scaling trees arent fixing your regression in some way. - Yep, this was happening and Brian could tell. Someone's seen too many scaling plots in his life >.<

-Get confidence intervals from individual SMA regressions. Plot them on the scaling comparisons. Theory over-predicts scaling exponents in the symmetrical case - accords with out observations about negative asymmetry

-Fix L_TOT to compute path fraction, or get it from branch-level TreeQSM output like you supposedly intended?
  -Properly add Parent IDs 
  
-While you're doing that, add a Fibonacci ratio to the length and radius of child branches to test for the golden ratio. 

-Double-check that asymmetry is properly aggregated (reconciliing ranges of values in positive/negative asymmetry)
  -Make sure that averaging asymmetry ratios does not cancel out positive/negative values from a theoretical POV
    -Do node-level asymmetry computations then aggregate those - there is no closed form solution for a network that
    mixes asymmetry types (but check with Alex).
      -Advantage of node-level computations: accomodate ternary and quaternary branching - only symmetric case!
        -Factor asymptotic equations for scaling and compute N at every node, apply finite-sized forms. Also apply
          finite-scaling approximation for each node
          RESULT: eh
