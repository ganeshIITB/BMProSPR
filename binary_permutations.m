%% generates binary digits for constructing states

function nn= binary_permutations(controller_nodes)

x = controller_nodes; %number of nodes
nn = zeros(1,x);
parfor i = 1:2^x
    a = dec2bin(i-1);
    n = [];
    k=x;
    j= length(a);
    while j>=1
     n(k)= str2num(a(j));
     k=k-1;
     j=j-1;
    end 
    nn(i,:)=n;
end
% save('StateSpace.mat','nn');        