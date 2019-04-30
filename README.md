Three projects here:
-QSM analysis
-Tree foraging simulation - arboreal communities
-L-Systems

Instructions:
-Run 'workflow.R' to pull all of the fit QSMs and perform a complete scaling analysis. Filenames are unique tree identifiers
for the entire workflow - data from other sources is never modified, filenames are used exactly.
-Testing/debugging - use the single test tree 'Ery_01' with guide at the top of workflow.R

TO-DO:
-Make column names readable and meaningful so that plotted figures are immediately readable
          
-Add columns for every conceivable branch ordering scheme - how does TreeQSM order their branches?
  -We MUST use the geometry of the tree to order these branches in some way that isn't biased.
    -Use maximum network size (N) estimates from branch ordering schemes, evaluate changes to scaling predictions
      -Find some way to correct the terminal element size of the network, maybe normalize by the size of the smallest element,
      or compare to a theoretical prediction of the capillary size.
      -Branch ordering schemes define terminal elements, then work backwards to identify the main axis. We'd like to do the       opposite, potentially doing something like using the smallest tip, or the longest path length tip, as the "true" tip, and defining other tips from there based on a similarity tolerance to the "true" tip (some kind of normalization).
      -Problem: WBE is built on the foundation of a topological ordering assumption - and that it precisely describes size scaling
      -Additional Problem: Both topo and Strahler ordering seriously misbehave under asymmetric conditions
        -Moving towards the biological meaning of (N), follow (Horn 2000) and graph ln(#branches) against branch order for diff
          branch ordering schemes to evaluate how branching intensity changes across the tree network
          
-Add stats to every analysis included in your Tuesday seminar talk. Confidence intervals on regression data points and
slopes you're drawing conclusions from.          
          
-Excluding nodes/branches based on scaling considerations might be affecting non-scaling related computations. Make sure this isn't the case, especially for path fraction - mainly check how your INVALID column is affecting calculations, reducing tips, etc.

-Think about how to modularize nodes and collections of nodes. 
  -Do a within-tree clustering analysis on raw geometry and scaling ratios to classify

-Integrate with L-Systems library (applying transformation rules to branch/cylinder data rows)

DONE/COMMENTS:
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
