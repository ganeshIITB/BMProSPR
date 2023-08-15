%% Boolean functions for tLGL network (taken from Saadatpour et al, 2011)

function Node1= t_cell(Node, order,input_nodes)

%Node(1)= S1P
%Node(2)= FLIP
%Node(3)= FAS
%Node(4)= Ceramide
%Node(5)= DISC
%Node(6)= Apoptosis

for i=1:length(order)
if any(input_nodes==order(i))
 else
switch order(i)
    case 1
        out= ~(Node(4) | Node(6));
    case 2
        out= ~(Node(5) | Node(6));
    case 3
        out= ~(Node(1) | Node(6));
    case 4
        out= Node(3) & ~(Node(1) | Node(6));
    case 5
        out= (Node(4) | (Node(3) & ~Node(2))) & ~Node(6);
    case 6
        out= Node(5) | Node(6);
end
Node(order(i))=out;
end
end
Node1= Node;