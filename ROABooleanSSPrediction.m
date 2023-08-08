% Codes for article titled "Modulating the dynamics of NFkB and PI3K enhances the ensemble-level TNFR1 signaling mediated apoptotic response" 
% Authors: Shubhank Sherekar, Chaitra S Todankar, Ganesh A Viswanathan
% Department of Chemical Engineering, Indian Institute of Technology Bombay, Mumbai, India.

% Codes for running BM-ProSPR algorithm on 6-node tLGL network (Saadatpour et al, 2011)
% This code runs on R2018b or later version
% Minimum 8 GB RAM is required to generate STG of minimal length. It uses Inner-outer iterations approach to calculate PageRank (Gliech et al, 2010). 
% One should download the publicly available script to calculate PageRank

% Note that if system permits enough RAM, uncomment line no. 113 and 114 and comment 116 to construct DiGraph network of respective permutations


%% generate state-space from biological network
load('nodes_t_cell.mat')

r=6;
StateSpace= GenerateStateSpace(r,update_nodes,input_nodes,fixed_nodes);
save('StateSpace','StateSpace')

%% generate permutations for state-space
rng shuffle                     % to shucdffle the random numbers from generator

seq=10000;
seq1=720;
iterations=720;
state_space_size=2^length(update_nodes);
rep=1:3;
order_map_generate(r,seq,iterations,state_space_size,rep,seq1)

%% Generate transitions, STG and 1-STM from...
... different no. of random permutations ...
... and different relative norm distance from ... 
.... analytical 1-STM and from previous no. of random permutations

mkdir('temp_pr')     % make directory for saving temp and pgrnk
mkdir('abs_probs')   % make directory for saving absorption probabilities
for it=rep
    tic     % time evaluation starts from here

    load(sprintf('order_files/order_map_%d',it))           % order list in container map
    order_obj= matfile(sprintf('order_files/order_ind_%d',it));             % random indexes of order map
    
    mkdir(sprintf('Output_replica_%d',it))      % make directories for saving files
    tic     % time evaluation starts from here
    
    fh= @t_cell;        % function handle (input)
    init_state_size= state_space_size;    % size of network
    net_nodes= size(order_map,1);           % no. of permutations

    Map_source_target= containers.Map('KeyType','single','ValueType','any');
    % initialize map with keys as index of source state and values as index of target state
    
    for i=1:init_state_size     % initializing map with source state containes zero target states
        Map_source_target(i)=[];
    end
    
    load('StateSpace')
    
    parfor i=1:init_state_size                      % keep the states in form of strings
        str_states{i}= sprintf('%d',StateSpace(i,:));
    end
    
    Map_states= containers.Map(str_states,1:init_state_size);   % map state-space with index
    clear str_states StateSpace

    transition_st{1}= sparse(init_state_size,init_state_size);      % initialize transition matrix
    adjacency_un{1}= sparse(init_state_size,init_state_size);   %initialize steady-state transition probability
    tau= zeros(iterations,1);


    for iter=2:iterations                    % loop for different permutations
        splits=5;        % number of subsets  in which state-space splits
        state_size_subset= ceil(init_state_size/splits); % subset size 
        itr=1;
        st_order=1;  
        weight=1;
        n_index_w= zeros(1,2);      % initialize pairs of source and target index
        order_ind1= order_obj.order_ind;
        order_ind= order_ind1(iter,:);
        clear order_ind1
        load('StateSpace')
        for itrt=1:splits   % loop for subsets
            incr=1;             % initialize increments
            
            if itrt<splits     % do the transitions till the condition is reached
                [Map_states,Map_source_target]= transition_states(StateSpace(itr:itr+state_size_subset,:),fh,Map_states,Map_source_target,order_map,order_ind(1,itr:itr+state_size_subset),input_nodes); % evaluate function for transition of states      
                key_m= keys(Map_source_target);     % extract index of source state
                key_map= key_m{itr}:key_m{itr+state_size_subset}; %extract  index of source state of only subsets  
                value_map= values(Map_source_target,key_m(itr:itr+state_size_subset));  % extract index of target state of only subsets
                
                [st_order1,weight1,n_index_w1]= unique_indexing(key_map,value_map);     % unique pairs of index of source and target states
            else     % Repeat the transition step for last subset (unknown size)
                [Map_states,Map_source_target]= transition_states(StateSpace(itr:end,:),fh,Map_states,Map_source_target,order_map,order_ind(1,itr:end),input_nodes);     
                key_m= keys(Map_source_target);
           
                key_map= key_m{itr}:key_m{length(Map_source_target)};
                value_map= values(Map_source_target,key_m(itr:length(Map_source_target)));
           
                [st_order1,weight1,n_index_w1]= unique_indexing(key_map,value_map);
            end
            
            st_order= [st_order, st_order1];    % updating source-target orders
            weight= [weight; weight1];      % weights of source-target pairs
            n_index_w= [n_index_w; n_index_w1];     % index of source-target pairs
            
            itr= itr+state_size_subset+1;       % updating index for next iterations
            st_order1=[];  weight1=[]; n_index_w1=[];
            clear key_map value_map
        end
        clear key_m StateSpace
        n_index_w(1,:)=[]; st_order(1)=[]; weight(1)=[]; % ignore dummy variables       
        
%         H= digraph(n_index_w(:,1),n_index_w(:,2),weight); % graph construction
%         adjacency_mat= adjacency(H,'weighted');

        adjacency_mat=sparse(double(n_index_w(:,1)),double(n_index_w(:,2)),weight,init_state_size,init_state_size);     % adjacency matrix  (adj)  with weights
        clear n_index_w weight
        adjacency_mat1= adjacency_mat(st_order,:);  % adj with consistent index
        adjacency_mat=[];
        adjacency_mat_cmpl= adjacency_mat1(:,st_order); % rearranging according to st orders
        adjacency_mat1=[];
        
        trans_prob= transition_probability(adjacency_mat_cmpl,iter); % transition probability from adjacency matrix
        save(sprintf('Output_replica_%d/graph_plot_%d',it,iter-1),'adjacency_mat_cmpl')
        adjacency_mat_cmpl=[];
        
        fp= find(diag(trans_prob)==1);
        adjacency_mat_un_cmpl=logical(trans_prob);
        
        [n_d_u(iter-1),n_d_f(iter-1),n_d_i(iter-1),n_m_u(iter-1),n_m_f(iter-1),n_m_i(iter-1),n_a_d_u(iter-1),n_a_d_f(iter-1),n_a_d_i(iter-1),n_a_m_u(iter-1),n_a_m_f(iter-1),n_a_m_i(iter-1),adjacency_un{iter},transition_st{iter}]= mat_comparison(adjacency_mat_un_cmpl,adjacency_un{iter-1},trans_prob,transition_st{iter-1});
        adjacency_mat_un_cmpl=[];
        if iter>2
            adjacency_un{iter-2}=[];
            transition_st{iter-2}=[];               % clear residual data for more space
        end
        
        pr1= pgrnk_estimation(trans_prob);
        pr(:,iter)=pr1';
        T_measure(iter)= temporality_measure_individual(n_a_d_u,n_a_m_u);
        
        [perm_cut_off,tau]= find_perm_cut_off(T_measure,pr(:,iter-1:iter),tau,iter);
        
            
        save(sprintf('Output_replica_%d/norms_all_t_%d',it,iter-1),'n_d_u','n_d_f','n_d_i','n_m_u','n_m_f','n_m_i','n_a_d_u','n_a_d_f','n_a_d_i','n_a_m_u','n_a_m_f','n_a_m_i')
        
        if isempty(perm_cut_off)==0
            cut_off_ss_prob(it)=perm_cut_off;
            break;
        end
        
        clear adjacency_mat_cmpl adjacency_mat_un_cmpl st_order weight n_index_w H
    end
    disc_pairs= (1-tau)/2;
    ss_prob= ss_probability(trans_prob);
    abs_prob= abs_probability(ss_prob,trans_prob);

    save(sprintf('temp_pr/T_measure_%d',it),'T_measure')
    save(sprintf('temp_pr/pgrnk_%d',it),'pr')
    save(sprintf('temp_pr/tau_%d',it),'tau')
    save(sprintf('temp_pr/disc_pairs_%d',it),'disc_pairs')
    save(sprintf('abs_probs/abs_prob_%d',it),'abs_prob')
    toc
    clearvars -except r seq state_space_size iterations input_nodes rep StateSpace cut_off_ss_prob
end

save('cut_off_ss_prob','cut_off_ss_prob')
