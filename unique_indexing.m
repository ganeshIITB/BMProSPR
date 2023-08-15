%% Enumerates the number of state transitions between unique state pairs

function [st_order,weights,n_index_w]= unique_indexing(n_index1,n_index2)

st_order= n_index1; % st order same as source index
n1= st_order;   
[~,n2]= cellfun(@size,n_index2);    % size of targets from sourcestate

j1=1;               %initializing for the loop

n_index_c= (repelem(n1,n2))';   % pairs of all source target state

n3= cell2mat(n_index2);     % convert to mat file
n_index_c(:,2)= n3';    %avoiding dummy values in lookup table

[n_index_w,~,freq_w]= unique(n_index_c,'rows');    %discarding redundant transitions

w_1= unique(freq_w);
weights=histc(freq_w(:),w_1);      %  weights of all unique pairs
