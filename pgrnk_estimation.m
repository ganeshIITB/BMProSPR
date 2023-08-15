%% taken from David F. Gleich et al, 2010 (an inner-outer iteration for computing pagerank)

function pr= pgrnk_estimation(trans_prob)

trans_prob= trans_prob-diag(diag(trans_prob));  % removes diagonal elements to ensure correct pagerank
n= length(trans_prob);
a= 0.85;
b= a/2;
maxiter= 1000;
tol= 1e-6;
itol= 1e-2;
v= ones(n,1)/n;
pr= inoutpr(trans_prob,a,v,tol,maxiter,b,itol);     % compute pagerank

