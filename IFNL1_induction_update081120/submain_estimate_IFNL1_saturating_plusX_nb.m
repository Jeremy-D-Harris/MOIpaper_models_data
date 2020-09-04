function results = submain_estimate_IFNL1_saturating_plusX_nb(data,deathrate_pars,k_dispersion)

z_score = 1.96;

% this is the submain file that fits three models to 8 and 18 hpi IFNL1
% data with model combinations:
% 0-8hpi:  (1) input-independent
% 8-18hpi: (1) input-independent, (2) linear, (3) saturating
% ---------------------------------- % ---------------------------------- %

% saturating (0-8 hpi) + model X (8-18 hpi)
disp('Fitting saturating (0-8 hpi) + model X (8-18 hpi) to data...'); fprintf('\n');

% MOI array: three replicates (cols), increaseing MOI (rows)
actual_moi = data.MOI;
size_inoculum = size(actual_moi);
actual_moi_collect = reshape(actual_moi',[1 size_inoculum(1)*size_inoculum(2)]);

% 8 hpi data
data_IFN_lambda1_8hpi = max(data.IFN_lambda1_8hpi,1);
data_IFN_lambda1_8hpi(1,2) = NaN;
size_data_8hpi = size(data_IFN_lambda1_8hpi);
data_IFN_lambda1_8hpi_vector = reshape(data_IFN_lambda1_8hpi',[1 size_data_8hpi(1)*size_data_8hpi(2)]);

% 18 hpi data
data_IFN_lambda1_18hpi = data.IFN_lambda1_18hpi;
size_data_18hpi = size(data_IFN_lambda1_18hpi);
data_IFN_lambda1_18hpi_vector = reshape(data_IFN_lambda1_18hpi',[1 size_data_18hpi(1)*size_data_18hpi(2)]);

% collect 8 and 18 hpi
data_IFNL1_collect = [data_IFN_lambda1_8hpi_vector;data_IFN_lambda1_18hpi_vector];
log2_data_IFNL1_collect = log2(data_IFNL1_collect);

% ---------------------------------- % ---------------------------------- %

% count number of NaN entries
n_pts_NaN = sum(sum(isnan(data_IFNL1_collect)));

% number of data points used
n_pts = size_data_8hpi(1)*size_data_8hpi(2)+ size_data_18hpi(1)*size_data_18hpi(2) - n_pts_NaN;

% ---------------------------------- % ---------------------------------- %

% MOI array to plot models IFNL1 fold change
actual_moi_finer = 10.^linspace(-2,2,101);

% plot on single-cell level
infection_events = 0:10;
infection_events_finer = linspace(0,10,101);

% rgb values - yellow, violet
default_color_rgb = [0.9290    0.6940    0.1250; 0.4940    0.1840    0.5560];
Marksize = [30 30 30];
datapoint=['.','.','.'];

% ----------  saturating  response (0-8hpi) + X (8-18 hpi)  ------------ %

for model_X = 1:3
    
    params_inputindependent.r = [];
    params_linear.r = [];
    params_saturating.r = [];
    params_saturating.K = [];
    
    m_saturating_init = 1115;
    K_saturating_init = 830;
    
    % specify which model (from 8-18 hpi)
    
    if model_X==1
        
        
        % 3-1: fit saturating (0-8 hpi) + input-independent (8-18 hpi)
        inputindependent_display = 'saturating (0-8 hpi) +  input-independent (8-18 hpi):';
        disp(inputindependent_display); fprintf('\n');
        
        l_init = 616; % input-independent rate (8-18hpi)
        y0 = log([m_saturating_init,K_saturating_init,l_init]);
        k_31 = length(y0);
        
        [log_lambda_inputindependent, LSE_IFNL1_31, exitflag, output] = fminsearch(@(y)getfit_IFN_saturating_plusX_nb(y,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion),y0);
        
        lambda_inputindependent = exp(log_lambda_inputindependent');
        
        var_log_lambda_m31 = getvar_IFN_saturating_plusX_nb(log_lambda_inputindependent,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion);
         
        CI_m31 = get_confidinterval(log_lambda_inputindependent',var_log_lambda_m31,z_score);

        disp(['m = ',num2str(lambda_inputindependent(1),'%1.2e'),' (0-8 hpi)', ', CI[', num2str(CI_m31(1,1),'%1.2e'),', ',num2str(CI_m31(1,2),'%1.2e'), '] ']);
        disp(['K = ',num2str(lambda_inputindependent(2),'%1.2e'),' (0-8 hpi)', ', CI[', num2str(CI_m31(2,1),'%1.2e'),', ',num2str(CI_m31(2,2),'%1.2e'), '] ']);
        disp(['l = ',num2str(lambda_inputindependent(3),'%2.0f'),' (8-18 hpi)', ', CI[', num2str(CI_m31(3,1),'%2.0f'),', ',num2str(CI_m31(3,2),'%2.0f'), '] ']);
        disp(['Least squares error = ',num2str(LSE_IFNL1_31,'%2.1f')]);  fprintf('\n\n');
        
        m_31_0to8hpi = lambda_inputindependent(1);
        K_31_0to8hpi = lambda_inputindependent(2);
        l_31_8to18hpi = lambda_inputindependent(3);
        
        % AIC is calculated from the sum LSE of 8 and 18 hpi
        AIC_inputindependent = n_pts*log(LSE_IFNL1_31)+2*k_31;
        
        % get 8 hpi
        params_saturating.r = m_31_0to8hpi;
        params_saturating.K = K_31_0to8hpi;
        
        % saturating model fold change: 0-8 hpi
        IFNL1_31_8hpi= get_IFN_saturating_nb(actual_moi_finer,params_saturating,deathrate_pars,data.hpi(1),k_dispersion);
        % get 18 hpi
        params_inputindependent.r = l_31_8to18hpi;
        % input-independent model fold change: 0-8 hpi
        IFNL1_31_0to8hpi = get_IFN_inputindependent_nb(actual_moi_finer,params_inputindependent,deathrate_pars,data.hpi(1),k_dispersion);
        % input-independent model fold change: 0-18 hpi
        IFNL1_31_0to18hpi = get_IFN_inputindependent_nb(actual_moi_finer,params_inputindependent,deathrate_pars,data.hpi(2),k_dispersion);
        % input-independent model fold change: 8-18 hpi
        IFNL1_31_8to18hpi= IFNL1_31_0to18hpi - IFNL1_31_0to8hpi;
        % saturating (0-8 hpi) + input-independent (8-18 hpi)
        IFNL1_31_18hpi = IFNL1_31_8hpi + IFNL1_31_8to18hpi;
        
        
%         pause;
        
    elseif model_X==2
        
        % 3-2: fit saturating (0-8 hpi) + linear dependent (8-18 hpi)
        linear_legend = 'saturating (0-8 hpi) +  linear (8-18 hpi):';
        disp(linear_legend); fprintf('\n');
        
        r_init = 161; % linear rate (8-18hpi)
        y0 = log([m_saturating_init,K_saturating_init,r_init]);
        k_32 = length(y0);
        
        [log_lambda_linear, LSE_IFNL1_32, exitflag, output] = fminsearch(@(y)getfit_IFN_saturating_plusX_nb(y,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion),y0);

        lambda_linear = exp(log_lambda_linear');
        
        var_log_lambda_m32 = getvar_IFN_saturating_plusX_nb(log_lambda_linear,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion);
         
        CI_m32 = get_confidinterval(log_lambda_linear',var_log_lambda_m32,z_score);

        disp(['m = ',num2str(lambda_linear(1),'%1.2e'),' (0-8 hpi)', ', CI[', num2str(CI_m32(1,1),'%1.2e'),', ',num2str(CI_m32(1,2),'%1.2e'), '] ']);
        disp(['K = ',num2str(lambda_linear(2),'%1.2e'),' (0-8 hpi)', ', CI[', num2str(CI_m32(2,1),'%1.2e'),', ',num2str(CI_m32(2,2),'%1.2e'), '] ']);
        disp(['r = ',num2str(lambda_linear(3),'%2.0f'),' (8-18 hpi)', ', CI[', num2str(CI_m32(3,1),'%2.0f'),', ',num2str(CI_m32(3,2),'%2.0f'), '] ']);
        disp(['Least squares error = ',num2str(LSE_IFNL1_32,'%2.1f')]);  fprintf('\n\n');
        
        m_32_0to8hpi = lambda_linear(1);
        K_32_0to8hpi = lambda_linear(2);
        r_32_8to18hpi = lambda_linear(3);
        
        % AIC is calculated from the sum LSE of 8 and 18 hpi
        AIC_linear = n_pts*log(LSE_IFNL1_32)+2*k_32;
        
        % get 8 hpi
        params_saturating.r = m_32_0to8hpi;
        params_saturating.K = K_32_0to8hpi;
        
        % saturating model fold change: 0-8 hpi
        IFNL1_32_8hpi= get_IFN_saturating_nb(actual_moi_finer,params_saturating,deathrate_pars,data.hpi(1),k_dispersion);
        % get 18 hpi
        params_linear.r = r_32_8to18hpi;
        % linear model fold change: 0-8 hpi
        IFNL1_32_0to8hpi = get_IFN_linear_nb(actual_moi_finer,params_linear,deathrate_pars,data.hpi(1),k_dispersion);
        % linear model fold change: 0-18 hpi
        IFNL1_32_0to18hpi= get_IFN_linear_nb(actual_moi_finer,params_linear,deathrate_pars,data.hpi(2),k_dispersion);
        % linear model fold change: 8-18 hpi
        IFNL1_32_8to18hpi= IFNL1_32_0to18hpi - IFNL1_32_0to8hpi;
        % saturating (0-8 hpi) + linear (8-18 hpi)
        IFNL1_32_18hpi = IFNL1_32_8hpi + IFNL1_32_8to18hpi;
        
%         pause;
        
    else
        
        % 3-3: fit input-independent (0-8 hpi) + saturating (8-18 hpi)
        saturating_legend = 'saturating (0-8 hpi) +  saturating (8-18 hpi):';
        disp(saturating_legend); fprintf('\n');
        
        m_init = 1190; % max saturating rate (8-18hpi)
        K_init = 1.94;   % input at half max saturating rate (8-18hpi)
        y0 = log([m_saturating_init,K_saturating_init,m_init,K_init]);
        k_33 = length(y0);

        [log_lambda_saturating, LSE_IFNL1_33, exitflag, output] = fminsearch(@(y)getfit_IFN_saturating_plusX_nb(y,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion),y0);

        lambda_saturating = exp(log_lambda_saturating');
        
        var_log_lambda_m33 = getvar_IFN_saturating_plusX_nb(log_lambda_saturating,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion);
         
        CI_m33 = get_confidinterval(log_lambda_saturating',var_log_lambda_m33,z_score);

        disp(['m = ',num2str(lambda_saturating(1),'%1.2e'),' (0-8 hpi)', ', CI[', num2str(CI_m33(1,1),'%1.2e'),', ',num2str(CI_m33(1,2),'%1.2e'), '] ']);
        disp(['K = ',num2str(lambda_saturating(2),'%1.2e'),' (0-8 hpi)', ', CI[', num2str(CI_m33(2,1),'%1.2e'),', ',num2str(CI_m33(2,2),'%1.2e'), '] ']);
        disp(['m = ',num2str(lambda_saturating(3),'%2.0f'),' (8-18 hpi)', ', CI[', num2str(CI_m33(3,1),'%2.0f'),', ',num2str(CI_m33(3,2),'%2.0f'), '] ']);
        disp(['K = ',num2str(lambda_saturating(4),'%2.2f'),' (8-18 hpi)', ', CI[', num2str(CI_m33(4,1),'%2.3f'),', ',num2str(CI_m33(4,2),'%2.2f'), '] ']);
        disp(['Least squares error = ',num2str(LSE_IFNL1_33,'%2.1f')]);  fprintf('\n\n');
        
        m_33_0to8hpi = lambda_saturating(1);
        K_33_0to8hpi = lambda_saturating(2);
        m_33_8to18hpi = lambda_saturating(3);
        K_33_8to18hpi = lambda_saturating(4);
        
        % AIC is calculated from the sum LSE of 8 and 18 hpi
        AIC_saturating = n_pts*log(LSE_IFNL1_33)+2*k_33;
        
        % get 8 hpi
        params_saturating.r = m_33_0to8hpi;
        params_saturating.K = K_33_0to8hpi;
        
        % saturating model fold change: 0-8 hpi
        IFNL1_33_8hpi= get_IFN_saturating_nb(actual_moi_finer,params_saturating,deathrate_pars,data.hpi(1),k_dispersion);
        % get 18 hpi
        params_saturating.r = m_33_8to18hpi;
        params_saturating.K = K_33_8to18hpi;
        % saturating model fold change: 0-8 hpi
        IFNL1_33_0to8hpi= get_IFN_saturating_nb(actual_moi_finer,params_saturating,deathrate_pars,data.hpi(1),k_dispersion);
        % saturating model fold change: 0-18 hpi
        IFNL1_33_0to18hpi= get_IFN_saturating_nb(actual_moi_finer,params_saturating,deathrate_pars,data.hpi(2),k_dispersion);
        % saturating model fold change: 8-18 hpi
        IFNL1_33_8to18hpi= IFNL1_33_0to18hpi - IFNL1_33_0to8hpi;
        % saturating (0-8 hpi) + saturating (8-18 hpi)
        IFNL1_33_18hpi = IFNL1_33_8hpi + IFNL1_33_8to18hpi;
        
        
        
    end
    
end


% two plots will appear:
% figure 1 - the model fits of INFL1 fold change vs. MOI (at 8 and 18 hpi)
% figure 2 - model fitted IFNL1 induction rates vs. hpi
% right panel:  model 2 + models 1-3


% figure 1: plot IFN fold change vs. increasing MOI
figure(1);
subplot(1,3,3);
for X=1:length(actual_moi)
    for i=1:3
        if i == 1
            
            if X ==1
                
                % plot the 8 hpi data
                h(1) = loglog(actual_moi(X,i), data_IFN_lambda1_8hpi(X,i), datapoint(i),'Color',default_color_rgb(1,:),'MarkerSize',Marksize(i),'LineWidth',1.5); hold on;
                
            end
            
        else
            
            loglog(actual_moi(X,i), data_IFN_lambda1_8hpi(X,i), datapoint(i),'Color',default_color_rgb(1,:),'MarkerSize',Marksize(i),'LineWidth',1.5); hold on;
            
        end
        
    end
end

% plot the input-independent model (0-8 hpi)
h(3) = loglog(actual_moi_finer, IFNL1_31_8hpi,'k', 'MarkerSize',12,'LineWidth',2); hold on;
h(4) = loglog(actual_moi_finer, IFNL1_32_8hpi,'Color',[0.5 0.5 0.5],'MarkerSize',12,'LineWidth',2); hold on;
h(5) = loglog(actual_moi_finer, IFNL1_33_8hpi,'--','Color',[0 0 0], 'MarkerSize',12,'LineWidth',2); hold on;


% plot the 18 hpi data
figure(1); subplot(1,3,3);
for X=1:length(actual_moi)
    for i=1:3
        if i == 1
            
            if X ==1
                
                h(2) = loglog(actual_moi(X,i), data_IFN_lambda1_18hpi(X,i), datapoint(i),'Color',default_color_rgb(2,:),'MarkerSize',Marksize(i),'LineWidth',1.5); hold on;
                
            end
            
        else
            
            loglog(actual_moi(X,i), data_IFN_lambda1_18hpi(X,i), datapoint(i),'Color',default_color_rgb(2,:),'MarkerSize',Marksize(i),'LineWidth',1.5); hold on;
        end
        
    end
end


loglog(actual_moi_finer, IFNL1_31_18hpi,'k', 'MarkerSize',12,'LineWidth',2); hold on;
loglog(actual_moi_finer, IFNL1_32_18hpi,'Color',[0.5 0.5 0.5],'MarkerSize',12,'LineWidth',2); hold on;
loglog(actual_moi_finer, IFNL1_33_18hpi,'--','Color',[0 0 0], 'MarkerSize',12,'LineWidth',2); hold on;
xlabel('actual MOI'); ylabel('IFN Fold Change');
title('saturating model (0-8 hpi)');
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
legend(h,{'8 hpi data','18 hpi data', 'input-independent (8-18 hpi)','linear (8-18 hpi)','saturating (8-18 hpi)'},'Location','Northwest');


% ---------------------------------- % ---------------------------------- %

% figure 2: plot rate models over time


time_hpi = linspace(0,18,101);
locs = find(time_hpi>8);

% time from 0-8hpi
time_0to8hpi = time_hpi(1:(locs(1)-1));
% time from 8-18hpi
time_8to18hpi = time_hpi(locs(1):end);

i_vary = [1 10 20];
for n=1:length(i_vary)
    
    i=i_vary(n);
    
    % saturating+input-independent
    IFNresponse_31(n,1:length(time_hpi)) = 1+m_31_0to8hpi*i/(K_31_0to8hpi+i);
    IFNresponse_31(n,1) = 0;
    IFNresponse_31(n,locs(1):end) = 1+l_31_8to18hpi;
    
    % saturating+linear
    IFNresponse_32(n,1:length(time_hpi)) = 1+m_32_0to8hpi*i/(K_32_0to8hpi+i);
    IFNresponse_32(n,1) = 0;
    IFNresponse_32(n,locs(1):end) = 1+r_32_8to18hpi*i;
    
    % saturating+saturating
    IFNresponse_33(n,1:length(time_hpi)) = 1+m_33_0to8hpi*i/(K_33_0to8hpi+i);
    IFNresponse_33(n,1) = 0;
    IFNresponse_33(n,locs(1):end) = 1+m_33_8to18hpi*i/(K_33_8to18hpi+i);
    
end


figure(2);
subplot(1,3,3);
g(1) = semilogy(time_hpi, IFNresponse_31(1,:),'k', 'MarkerSize',12,'LineWidth',2); hold on;
semilogy(time_hpi, IFNresponse_31(2:3,:),'k','MarkerSize',12,'LineWidth',2); hold on;
g(2) = semilogy(time_hpi, IFNresponse_32(1,:),'Color',[0.5 0.5 0.5], 'MarkerSize',12,'LineWidth',2); hold on;
semilogy(time_hpi, IFNresponse_32(2:3,:),'Color',[0.5 0.5 0.5], 'MarkerSize',12,'LineWidth',2); hold on;
g(3) = semilogy(time_hpi, IFNresponse_33(1,:),'--','Color',[0 0 0], 'MarkerSize',12,'LineWidth',2);
semilogy(time_hpi, IFNresponse_33(2:3,:),'--','Color',[0 0 0], 'MarkerSize',12,'LineWidth',2);
xlabel('hpi'); ylabel('IFN induction rate');
title('saturating rate (0-8 hpi)');
axis([0 18 10^0 1e5])
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
legend(g,{'input-independent (8-18 hpi)','linear (8-18 hpi)','saturating (8-18 hpi)'},'Location','Northwest');


% ---------------------------------- % ---------------------------------- %

% return results to main m-file

AIC_values = [AIC_inputindependent,AIC_linear,AIC_saturating];
results.AIC_values = AIC_values;

results.estimate_m31 = lambda_inputindependent;
results.estimate_m32 = lambda_linear;
results.estimate_m33 = lambda_saturating;

results.CI_m31 = CI_m31;
results.CI_m32 = CI_m32;
results.CI_m33 = CI_m33;

% model results: concatenate
results.actual_moi_finer = actual_moi_finer';
% 8 hpi
results.IFNL1_3X_8hpi = [IFNL1_31_8hpi',IFNL1_32_8hpi',IFNL1_33_8hpi'];
% 18 hpi
results.IFNL1_3X_18hpi = [IFNL1_31_18hpi',IFNL1_32_18hpi',IFNL1_33_18hpi'];

results.time_hpi = time_hpi';
results.IFNresponse_3X = [IFNresponse_31',IFNresponse_32',IFNresponse_33'];

