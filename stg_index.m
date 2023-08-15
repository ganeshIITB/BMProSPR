%% saves the state transitions as string format in database 

function str_index= stg_index(stg,Node1,str_index,Map_states,Map_source_target)

str_s= sprintf('%d',stg);
str_t= sprintf('%d',Node1);

index_s= values(Map_states,{str_s});
index_t= values(Map_states,{str_t});

v= values(Map_source_target,index_s);
index_t_ex= [v{1},index_t{1}];
Map_source_target(index_s{1})= index_t_ex;

