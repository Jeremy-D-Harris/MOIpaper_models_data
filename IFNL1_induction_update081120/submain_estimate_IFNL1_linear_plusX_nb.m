function results = submain_estimate_IFNL1_linear_plusX_nb(data,deathrate_pars,k_dispersion)

z_score = 1.96;

% this is the submain file that fits three models to 8 and 18 hpi IFNL1
% data with model combinations:
% 0-8hpi:  (2) linear
% 8-18hpi: (1) input-independent, (2) linear, (3) saturating
% ---------------------------------- % ---------------------------------- %

fprintf('\n')
disp('% ---------------------------------- % ---------------------------------- %');
fprintf('\n')

% linear (0-8 hpi) + model X (8-18 hpi)
disp('Fitting linear (0-8 hpi) + model X (8-18 hpi) to data...'); fprintf('\n');

fprintf('\n')
disp('% ---------------------------------- % ---------------------------------- %');
fprintf('\n')

% ---------------------------------- % ---------------------------------- %

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

% -------------  linear response (0-8hpi) + X (8-18 hpi)  --------------- %

for model_X = 1:3
    
    params_inputindependent.r = [];
    params_linear.r = [];
    params_saturating.r = [];
    params_saturating.K = [];
    
    r_linear_init = 1.34;
    
    % specify which model (from 8-18 hpi)
    
    if model_X==1
        
        
        % 2-1: fit linear (0-8 hpi) + input-independent (8-18 hpi)
        
        r_init = 616; % input-independent rate (8-18hpi)
        y0 = log([r_linear_init,r_init]);
        k_21 = length(y0);

        inputindependent_display = 'linear (0-8 hpi) +  input-independent (8-18 hpi):';
        disp(inputindependent_display); fprintf('\n');
        
        [log_lambda_inputindependent, LSE_IFNL1_21, exitflag, output] = fminsearch(@(y)getfit_IFN_linear_plusX_nb(y,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion),y0);

        lambda_inputindependent = exp(log_lambda_inputindependent');
        
        var_log_lambda_m21 = getvar_IFN_linear_plusX_nb(log_lambda_inputindependent,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion);
        
        CI_m21 = get_confidinterval(log_lambda_inputindependent',var_log_lambda_m21,z_score);

        disp(['r = ',num2str(lambda_inputindependent(1),'%2.2f'),' (0-8 hpi)', ', CI[', num2str(CI_m21(1,1),'%2.2f'),', ',num2str(CI_m21(1,2),'%2.2f'), '] ']);
        disp(['l = ',num2str(lambda_inputindependent(2),'%2.0f'),' (8-18 hpi)',', CI[', num2str(CI_m21(2,1),'%2.2f'),', ',num2str(CI_m21(2,2),'%2.2f'), '] ']);
        disp(['Least squares error = ',num2str(LSE_IFNL1_21,'%2.1f')]);  fprintf('\n\n');
        
        r_21_0to8hpi = lambda_inputindependent(1);
        l_21_8to18hpi = lambda_inputindependent(2);
        
        % AIC is calculated from the sum LSE of 8 and 18 hpi
        AIC_inputindependent = n_pts*log(LSE_IFNL1_21)+2*k_21;
        
        % get 8 hpi
        params_linear.r = r_21_0to8hpi;
        
        % linear model fold change: 0-8 hpi
        IFNL1_21_8hpi= get_IFN_linear_nb(actual_moi_finer,params_linear,deathrate_pars,data.hpi(1),k_dispersion);
        % get 18 hpi
        params_inputindependent.r = l_21_8to18hpi;
        % input-independent model fold change: 0-8 hpi
        IFNL1_21_0to8hpi = get_IFN_inputindependent_nb(actual_moi_finer,params_inputindependent,deathrate_pars,data.hpi(1),k_dispersion);
        % input-independent model fold change: 0-18 hpi
        IFNL1_21_0to18hpi = get_IFN_inputindependent_nb(actual_moi_finer,params_inputindependent,deathrate_pars,data.hpi(2),k_dispersion);
        % input-independent model fold change: 8-18 hpi
        IFNL1_21_8to18hpi= IFNL1_21_0to18hpi - IFNL1_21_0to8hpi;
        % linear (0-8 hpi) + input-independent (8-18 hpi)
        IFNL1_21_18hpi = IFNL1_21_8hpi + IFNL1_21_8to18hpi;
        
%         pause;
        
    elseif model_X==2
        
        % 2-2: fit linear (0-8 hpi) + linear dependent (8-18 hpi)
        linear_display = 'linear (0-8 hpi) +  linear (8-18 hpi):';
        disp(linear_display); fprintf('\n');
        
        r_init = 161; % linear rate (8-18hpi)
        y0 = log([r_linear_init,r_init]);
        k_22 = length(y0);
        
        [log_lambda_linear, LSE_IFNL1_22, exitflag, output] = fminsearch(@(y)getfit_IFN_linear_plusX_nb(y,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion),y0);

        lambda_linear = exp(log_lambda_linear');
        
        var_log_lambda_m22 = getvar_IFN_linear_plusX_nb(log_lambda_linear,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion);
        
        CI_m22 = get_confidinterval(log_lambda_linear',var_log_lambda_m22,z_score);

        disp(['r = ',num2str(lambda_linear(1),'%2.2f'),' (0-8 hpi)', ', CI[', num2str(CI_m22(1,1),'%2.3f'),', ',num2str(CI_m22(1,2),'%2.2f'), '] ']);
        disp(['r = ',num2str(lambda_linear(2),'%2.0f'),' (8-18 hpi)', ', CI[', num2str(CI_m22(2,1),'%2.0f'),', ',num2str(CI_m22(2,2),'%2.0f'), '] ']);
        disp(['Least squares error = ',num2str(LSE_IFNL1_22,'%2.1f')]);  fprintf('\n\n');
        
        r_22_0to8hpi = lambda_linear(1);
        r_22_8to18hpi = lambda_linear(2);
        
        % AIC is calculated from the sum LSE of 8 and 18 hpi
        AIC_linear = n_pts*log(LSE_IFNL1_22)+2*k_22;
        
        % get 8 hpi
        params_linear.r = r_22_0to8hpi;
        % linear model fold change: 0-8 hpi
        IFNL1_22_8hpi= get_IFN_linear_nb(actual_moi_finer,params_linear,deathrate_pars,data.hpi(1),k_dispersion);
        
        % get 18 hpi
        params_linear.r = r_22_8to18hpi;
        % linear model fold change: 0-8 hpi
        IFNL1_22_0to8hpi = get_IFN_linear_nb(actual_moi_finer,params_linear,deathrate_pars,data.hpi(1),k_dispersion);
        % linear model fold change: 0-18 hpi
        IFNL1_22_0to18hpi= get_IFN_linear_nb(actual_moi_finer,params_linear,deathrate_pars,data.hpi(2),k_dispersion);
        % linear model fold change: 8-18 hpi
        IFNL1_22_8to18hpi= IFNL1_22_0to18hpi - IFNL1_22_0to8hpi;
        % linear (0-8 hpi) + linear (8-18 hpi)
        IFNL1_22_18hpi = IFNL1_22_8hpi + IFNL1_22_8to18hpi;
        
%         pause;
        
    else
        
        % 2-3: fit linear (0-8 hpi) + saturating (8-18 hpi)
        saturating_display = 'linear (0-8 hpi) +  saturating (8-18 hpi):';
        disp(saturating_display); fprintf('\n');
        
        m_init = 1188; % max saturating rate (8-18hpi)
        K_init = 1.94; % input at half max saturating rate (8-18hpi)
        y0 = log([r_linear_init,m_init,K_init]);
        k_23 = length(y0);
        
        [log_lambda_saturating, LSE_IFNL1_23, exitflag, output] = fminsearch(@(y)getfit_IFN_linear_plusX_nb(y,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion),y0);
        
        lambda_saturating = exp(log_lambda_saturating');
        
        var_log_lambda_m23 = getvar_IFN_linear_plusX_nb(log_lambda_saturating,actual_moi_collect,log2_data_IFNL1_collect,params_inputindependent,params_linear,params_saturating,deathrate_pars,data.hpi,model_X,k_dispersion);
        
        CI_m23 = get_confidinterval(log_lambda_saturating',var_log_lambda_m23,z_score);
        
        disp(['r = ',num2str(lambda_saturating(1),'%2.2f'),' (0-8 hpi)', ', CI[', num2str(CI_m23(1,1),'%2.2f'),', ',num2str(CI_m23(1,2),'%2.2f'), '] ']);
        disp(['m = ',num2str(lambda_saturating(2),'%2.0f'),' (8-18 hpi)', ', CI[', num2str(CI_m23(2,1),'%2.0f'),', ',num2str(CI_m23(2,2),'%2.0f'), '] ']);
        disp(['K = ',num2str(lambda_saturating(3),'%2.2f'),' (8-18 hpi)', ', CI[', num2str(CI_m23(3,1),'%2.2f'),', ',num2str(CI_m23(3,2),'%2.2f'), '] ']);
        disp(['Least squares error = ',num2str(LSE_IFNL1_23,'%2.1f')]);  fprintf('\n\n');
        
        r_23_0to8hpi = lambda_saturating(1);
        m_23_8to18hpi = lambda_saturating(2);
        K_23_8to18hpi = lambda_saturating(3);
        
        % AIC is calculated from the sum LSE of 8 and 18 hpi
        AIC_saturating = n_pts*log(LSE_IFNL1_23)+2*k_23;
        
        % get 8 hpi
        params_linear.r = r_23_0to8hpi;
        % linear model fold change: 0-8 hpi
        IFNL1_23_8hpi= get_IFN_linear_nb(actual_moi_finer,params_linear,deathrate_pars,data.hpi(1),k_dispersion);
        
        
        % get 18 hpi
        params_saturating.r = m_23_8to18hpi;
        params_saturating.K = K_23_8to18hpi;
        % saturating model fold change: 0-8 hpi
        IFNL1_23_0to8hpi= get_IFN_saturating_nb(actual_moi_finer,params_saturating,deathrate_pars,data.hpi(1),k_dispersion);
        % saturating model fold change: 0-18 hpi
        IFNL1_23_0to18hpi= get_IFN_saturating_nb(actual_moi_finer,params_saturating,deathrate_pars,data.hpi(2),k_dispersion);
        % saturating model fold change: 8-18 hpi
        IFNL1_23_8to18hpi= IFNL1_23_0to18hpi - IFNL1_23_0to8hpi;
        % linear (0-8 hpi) + saturating (8-18 hpi)
        IFNL1_23_18hpi = IFNL1_23_8hpi + IFNL1_23_8to18hpi;
        
        
        
    end
    
end

% ---------------------------------- % ---------------------------------- %


% two plots will appear:
% figure 1 - the model fits of INFL1 fold change vs. MOI (at 8 and 18 hpi)
% figure 2 - model fitted IFNL1 induction rates vs. hpi
% middle panel:  model 2 + models 1-3


% figure 1: plot IFN fold change vs. increasing MOI
figure(1);
subplot(1,3,2);
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

% plot the linear model (0-8 hpi)
h(3) = loglog(actual_moi_finer, IFNL1_21_8hpi,'k', 'MarkerSize',12,'LineWidth',2); hold on;
h(4) = loglog(actual_moi_finer, IFNL1_22_8hpi,'Color',[0.5 0.5 0.5],'MarkerSize',12,'LineWidth',2); hold on;
h(5) = loglog(actual_moi_finer, IFNL1_23_8hpi,'--','Color',[0 0 0], 'MarkerSize',12,'LineWidth',2); hold on;


% plot the 18 hpi data
figure(1); subplot(1,3,2);
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


loglog(actual_moi_finer, IFNL1_21_18hpi,'k', 'MarkerSize',12,'LineWidth',2); hold on;
loglog(actual_moi_finer, IFNL1_22_18hpi,'Color',[0.5 0.5 0.5],'MarkerSize',12,'LineWidth',2); hold on;
loglog(actual_moi_finer, IFNL1_23_18hpi,'--','Color',[0 0 0], 'MarkerSize',12,'LineWidth',2); hold on;
xlabel('actual MOI'); ylabel('IFN Fold Change');
title('linear model (0-8 hpi)');
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
    
    % linear+input-independent
    IFNresponse_21(n,1:length(time_hpi)) = 1+r_21_0to8hpi*i;
    IFNresponse_21(n,1) = 0;
    IFNresponse_21(n,locs(1):end) = 1+l_21_8to18hpi;

    % linear+linear
    IFNresponse_22(n,1:length(time_hpi)) = 1+r_22_0to8hpi*i;
    IFNresponse_22(n,1) = 0;
    IFNresponse_22(n,locs(1):end) = 1+r_22_8to18hpi*i;
    
    % linear+saturating
    IFNresponse_23(n,1:length(time_hpi)) = 1+r_23_0to8hpi*i;
    IFNresponse_23(n,1) = 0;
    IFNresponse_23(n,locs(1):end) = 1+m_23_8to18hpi*i/(K_23_8to18hpi+i);
    
end


figure(2);
subplot(1,3,2);
g(1) = semilogy(time_hpi, IFNresponse_21(1,:),'k', 'MarkerSize',12,'LineWidth',2); hold on;
semilogy(time_hpi, IFNresponse_21(2:3,:),'k', 'MarkerSize',12,'LineWidth',2); hold on;
g(2) = semilogy(time_hpi, IFNresponse_22(1,:),'Color',[0.5 0.5 0.5], 'MarkerSize',12,'LineWidth',2); hold on;
semilogy(time_hpi, IFNresponse_22(2:3,:),'Color',[0.5 0.5 0.5], 'MarkerSize',12,'LineWidth',2); hold on;
g(3) = semilogy(time_hpi, IFNresponse_23(1,:),'--','Color',[0 0 0], 'MarkerSize',12,'LineWidth',2);
semilogy(time_hpi, IFNresponse_23(2:3,:),'--','Color',[0 0 0], 'MarkerSize',12,'LineWidth',2);
xlabel('hpi'); ylabel('IFN induction rate');
title('linear rate (0-8 hpi)');
axis([0 18 10^0 1e5])
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
legend(g,{'input-independent (8-18 hpi)','linear (8-18 hpi)','saturating (8-18 hpi)'},'Location','Northwest');


% ---------------------------------- % ---------------------------------- %

% return results to main m-file

AIC_values = [AIC_inputindependent,AIC_linear,AIC_saturating];
results.AIC_values = AIC_values;

results.estimate_m21 = lambda_inputindependent;
results.estimate_m22 = lambda_linear;
results.estimate_m23 = lambda_saturating;

results.CI_m21 = CI_m21;
results.CI_m22 = CI_m22;
results.CI_m23 = CI_m23;

% model results: concatenate
results.actual_moi_finer = actual_moi_finer';
% 8 hpi
results.IFNL1_2X_8hpi = [IFNL1_21_8hpi',IFNL1_22_8hpi',IFNL1_23_8hpi'];
% 18 hpi
results.IFNL1_2X_18hpi = [IFNL1_21_18hpi',IFNL1_22_18hpi',IFNL1_23_18hpi'];

results.time_hpi = time_hpi';
results.IFNresponse_2X = [IFNresponse_21',IFNresponse_22',IFNresponse_23'];