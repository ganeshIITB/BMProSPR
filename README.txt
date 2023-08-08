DESCRIPTION
-----------

This is a set of MATLAB codes for implementing BM-ProSPR algorithm detailed in paper titled "Modulating the dynamics of NFkB and PI3K enhances the ensemble-level TNFR1 signaling mediated apoptotic response"  

INPUT FILES

1. nodes_t_cell.mat- includes fixed nodes as zero, input nodes as 1 and update nodes as dynamically varying nodes


FILES

1. ROABooleanSSPrediction.m- Main files to run BM-ProSPR algorithm

2. GenerateStateSpace.m- generates a statespace of dynamically varying nodes

3. order_map_generate.m- Generates random unique permutations

4. binary_permutations.m- generates binary digits for constructing states

5. transition_states.m- calculates transition probabilities between all pair of states

6. t_cell.m- Boolean functions for tLGL network (taken from Saadatpour et al, 2011)

7. stg_index.m- saves the state transitions as string format in database

8. unique_indexing.m- Enumerates the number of state transitions between unique state pairs

9. transition_probability.m- calculates transition probabilities between all pair of states

10. mat_comparison.m- Compares the state transition matrix of two successive permutations

11. temporality_measure_individual.m- Temporality measure is taken from Aming Li et al, 2020

12. pgrnk_estimation.m- taken from David F. Gleich et al, 2010 (an inner-outer iteration for computing pagerank). One should download the publicly available script to calculate PageRank (https://www.cs.purdue.edu/homes/dgleich/codes/innout/)

13. tau_measure.m- measures kendall's tau based on PageRank attained by states

14. find_perm_cut_off.m- Find permutation cut-off using error threshold

15. ss_probability.m- Generates absorption probability to reach fixed points from every transient states

16. abs_probability.m- Generates absorption probability to reach fixed points from every states