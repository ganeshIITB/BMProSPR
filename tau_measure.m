%% measures kendall's tau based on PageRank attained by states

function tau= tau_measure(pr)

for i1=2:size(pr,2)
    t= py.scipy.stats.kendalltau(pr(:,i1-1),pr(:,i1));
    % py used for only predefined modules; if there is...
    ...need to use python script, use system('python3 <py file name>')
    tau(i1-1)=t{1};
end
