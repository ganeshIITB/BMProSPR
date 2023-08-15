%% Compares the state transition matrix of two successive permutations

function [n_d_u,n_d_f,n_d_i,n_m_u,n_m_f,n_m_i,n_a_d_u,n_a_d_f,n_a_d_i,n_a_m_u,n_a_m_f,n_a_m_i,adj_curr,T_curr]= mat_comparison(adj_curr,adj_prev,T_curr,T_prev)

comp= T_curr-T_prev;        % difference in elements of transition probabilites

n_d_u= sum(sum(abs(comp))); % unity norm of difference in prob.
n_d_f= sqrt(sum(sum(comp.*comp)));    % Frobenius norm of diff. in prob.
n_d_i= max(max(comp));        %infinity norm of diff. in prob.

n_m_u= sum(sum(abs(T_curr))); % unity, Frob, inifinity norm of current prob.
n_m_f= sqrt(sum(sum(T_curr.*T_curr)));
n_m_i= max(max(T_curr));

comp_adj= adj_curr-adj_prev;    % difference in connectivity

n_a_d_u= sum(sum(abs(comp_adj)));    %unity, Frobenius and infinity norm of diff in adj.
n_a_d_f= sqrt(sum(sum(comp_adj.*comp_adj)));
n_a_d_i= max(max(comp_adj));

n_a_m_u= sum(sum(abs(adj_curr)));    %unity, Frobenius and infinity norm of adj. 
n_a_m_f= sqrt(sum(sum(adj_curr.*adj_curr)));
n_a_m_i= max(max(adj_curr));
