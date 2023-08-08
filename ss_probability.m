%% Generates absorption probability to reach fixed points from every transient states

function ss_prob= ss_probability(trans_prob)

fp= find(diag(trans_prob)==1);
non_fp= setdiff(1:length(trans_prob),fp);

new_ind= [fp; non_fp'];
new_c= trans_prob(new_ind,:);
new_r= new_c(:,new_ind);

trans_prob_new= new_r;

I= trans_prob_new(1:length(fp),1:length(fp));
Z= trans_prob_new(length(fp)+1,1:length(fp));
R= trans_prob_new(length(fp)+1:end,1:length(fp));
Q= trans_prob_new(length(fp)+1:end,length(fp)+1:end);
Y= speye(length(Q),length(Q));
M= Y-Q;

b= ones(length(trans_prob),1)/length(trans_prob);
tol=1e-4;
maxit=1000;

for i=1:length(fp)
    ss_prob(:,i)= bicgstab(M,R(:,i),tol,maxit);
end