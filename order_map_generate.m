%% Generates random unique permutations 

function order_map_generate(r,seq,iterations,state_space,replica,seq1)
rng shuffle
mkdir order_files   % folder to save random order map and random index of maps

for it=replica         % the number of files for map and index
    ord_map= containers.Map('0 0 0 0',0); %temporary map to make template of orders
    i1=1;   %initialize the index of maps
    
    for i=1:seq     % loop to generate orders
        order(i,:)= randperm(r);    % random order sequence generated
        str_order= char(join(string(order(i,:))));  % convert into string or char format
    
        if isKey(ord_map,str_order)==1  % to check if this order sequence is already taken or not
        else    
            ord_map(str_order)= i1;    % if not, order sequence is added in mapped
            i1=i1+1;
        end
    end
    
    v= values(ord_map);     % collect index from map
    v(1)=[];                        % first one is template, needs to be removed
    k= keys(ord_map);       % collect the order sequence from map
    k(1)=[];
    
    clear ord_map
    order_map= containers.Map(v,k);     % map with index as key and order sequence as values
    
    order_ind= zeros(iterations,state_space);   % initialize of size state-space
    rng shuffle
    parfor i=1:state_space         % loop for random index of maps
        order_ind(:,i)= randperm(seq1,iterations);
    end
    order_ind= single(order_ind);

    save(sprintf('order_files/order_map_%d',it),'order_map')    %saving maps and index
    save(sprintf('order_files/order_ind_%d',it),'order_ind')
    clearvars -except r seq iterations state_space seq1
end