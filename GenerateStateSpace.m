%% This generates a statespace of dynamically varying nodes

function StateSpace= GenerateStateSpace(ndn,up_nd,in_nd,fx_nd)

nn= binary_permutations(length(up_nd)); % all permutations of update nodes

parfor i=1:length(nn)
    Node = zeros(ndn,1);
    Node(in_nd')= ones(length(in_nd),1);  %creating vectors with index of nodes and give the value as 1 to them
    Node(fx_nd)= zeros(length(fx_nd),1);    % nodes having fixed state remained as 0 binary state
    Node(up_nd)= nn(i,:);       % all other nodes are given random binary states
    StateSpace(i,:)= Node;      % Save every random binary states to create network state
end