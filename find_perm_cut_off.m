%% Find permutation cut-off using error threshold

function [perm_cut_off,tau]= find_perm_cut_off(T_measure,pr,tau,iter)


t= py.scipy.stats.kendalltau(pr(:,1),pr(:,2));
tau(iter-1)=t{1};


tm_diff= abs(diff(T_measure));
tau_cut_off= 0.9999;
D_cut_off=1e-4;
temp_cut_off= 1e-4;
tau_critical1= find(tau(2:end)>tau_cut_off);
tau_critical2= diff(tau_critical1)==1;
tau_critical3= bwlabel(tau_critical2);
props= regionprops(tau_critical3,'Area','PixelList');
tau_critical4= [props.Area];
tau_critical5= find(tau_critical4>=10,1);

if isempty(tau_critical5)==1
    tau_critical5=0;
end

tau_critical= tau_critical1(find(tau_critical3==tau_critical5,1));

perm_cut_off= max(find(tm_diff<temp_cut_off,1),tau_critical);