%% calculates transition probabilities between all pair of states

function trans_prob= transition_probability(adjacency_mat_cmpl,iter)

trans_prob= adjacency_mat_cmpl/(iter-1);    %calculating transition probabilities
