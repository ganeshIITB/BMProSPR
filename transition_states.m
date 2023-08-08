%% Captures transitions from each states based on random unique permutations 

function [Map_states,Map_source_target] = transition_states(StateSpace,fh,Map_states,Map_source_target,order_map,order_ind,input_nodes)

str_index=size(Map_states,1)+1;

for i=1:size(StateSpace,1)  %loop for all state-space
   
    stg= StateSpace(i,:)';      %the source state
    
    order_ind1= order_ind(i); 
    order_str= order_map(order_ind1);
    order= sscanf(order_str,'%d');  %  permutation

    Node1= fh(stg, order,input_nodes);  % update step
    str_index= stg_index(stg,Node1,str_index,Map_states,Map_source_target); % keeping index variables in maps
end
