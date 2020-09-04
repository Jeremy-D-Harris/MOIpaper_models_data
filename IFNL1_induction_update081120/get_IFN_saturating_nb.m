function IFN_foldchange = get_IFN_saturating_nb(moi,params_saturating,alpha,this_hpi,k_dispersion)

matlab_r = k_dispersion;

% cell death rate: time-dependent, input-independent: k,b,mu
k = alpha(1); % shape parameter
b = alpha(2); % scale parameters
mu = alpha(3); % constant death rate

dt = 0.01;
time = 0:dt:this_hpi;

% time-dependent, input-independent cell death
int_psi = b*time.^k + mu*time; % integrating Weibull hazard function

prob_cellsremain_t = exp(-int_psi); % function of time

N_W = 1000;
num_WT_vector = 0:N_W;

m = params_saturating.r;
% m = params_MM.m;
K = params_saturating.K;

response_i_vector = 1+m*num_WT_vector./(K+num_WT_vector);

for i =1:length(moi)
    
    this_moi = moi(1,i);
    
    matlab_p = 1-this_moi/(matlab_r + this_moi); % actually probability of failure, I think, in matlab speak
    
    % set up WT vector
    prob_psiW_vector = nbinpdf(num_WT_vector,matlab_r,matlab_p);%poisspdf(num_WT_vector,psiW);
    prob_above_N_W = 1 - nbincdf(N_W,matlab_r,matlab_p);%poisscdf(N_W,psiW);
    prob_psiW_vector(1,length(num_WT_vector)) = prob_psiW_vector(end)+prob_above_N_W;
    
    % probability of no DI and at least one WT
    IFN_foldchange(i) = prob_psiW_vector(1)+sum(response_i_vector(2:end).*prob_psiW_vector(2:end))*sum(prob_cellsremain_t)*dt;
    
    
end