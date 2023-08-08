%% Temporality measure is taken from Aming Li et al, 2020

function T= temporality_measure_individual(n_a_d_u,n_a_m_u)

T_temp= cumsum(n_a_d_u./n_a_m_u);   % cumulative relative difference of connectivity

T= T_temp(end)/length(T_temp);