%% Generating absorption probability to reach fixed points from every states

function abs_prob= abs_probability(ss_prob,trans_prob)

fp= find(diag(trans_prob)==1);
range= 1:length(trans_prob);

range= setdiff(range,fp);
abs_prob= zeros(length(trans_prob),length(fp));


abs_prob(range,:)=ss_prob;

for i=1:length(fp)
    abs_prob(fp(i),i)=1;
end
