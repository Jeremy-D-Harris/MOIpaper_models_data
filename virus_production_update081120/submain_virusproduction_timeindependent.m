function results_timeindependent = submain_virusproduction_timeindependent(data,deathrate_pars,k_dispersion,nu_init_timeindependent,cell_line)

z_score = 1.96;

which_time_model = 0; % constant production rate in time

% update 02/20/20

% script to fit three time-independent models of virus production with 
% input-dependent assumptions, models 1-3:
% 'model 1': input-independent
% 'model 2': linear input-dependence
% 'model 3': saturating input-dependence

% ---------------------------------- % ---------------------------------- %

fprintf('\n')
disp('% ---------------------------------- % ---------------------------------- %');
fprintf('\n')

disp('Fitting time-independent models of virus production to data...'); fprintf('\n');

fprintf('\n')
disp('% ---------------------------------- % ---------------------------------- %');
fprintf('\n')

% ---------------------------------- % ---------------------------------- %

% previously calculated Bulk MOI values
actual_moi = data.actual_moi;

% reshape treatments into a 1D vector
moi_vector = reshape(actual_moi',[1 data.n_datapoints]);

moi = 10.^linspace(log10(0.001),log10(10),101); % finer moi array for plotting models

% ---------------------------------- % ---------------------------------- %

% include virus production in row vector
data_6hpi = reshape(data.hpi_6hrs',[1, data.n_datapoints]);
data_12hpi = reshape(data.hpi_12hrs',[1, data.n_datapoints]);
data_18hpi = reshape(data.hpi_18hrs',[1, data.n_datapoints]);

% concatenate 6,12,18 hpi by row
data_allhpi = [data_6hpi; data_12hpi; data_18hpi];
size_data_allhpi = size(data_allhpi);
n_pts_allhpi = size_data_allhpi(1)*size_data_allhpi(2);

% count the number of NaNs in the data
number_NaN = sum(sum(isnan(data_allhpi)));

% total number of data points that go into the fitting
n_pts_allhpi = n_pts_allhpi - number_NaN;

Ncells = data.Ncells; % cell count ~ 2 million cells

% ---------------------------------- % ---------------------------------- %


% 'model 1': time-independent, input-independent rate of virus production
model_num = 1;

inputindependent_legend = 'time-independent, input-independent rate of virus production:';
disp(inputindependent_legend); fprintf('\n');


x0_inputindependent = log(nu_init_timeindependent(1)); % nu(t,i) = l
k_inputindependent = length(x0_inputindependent);

options = optimset('MaxFunEvals',500,'MaxIter',500,'TolFun',1e-8,'TolX',1e-8);
% minimize residual sum squared error
[log_nu_inputindependent, LSE_GE_inputindependent] = fminsearch(@(x)getfit_GEmodel(x, moi_vector, data_allhpi, data.hpi, deathrate_pars, k_dispersion, Ncells, cell_line, model_num, which_time_model),x0_inputindependent,options);

AIC_model1 = n_pts_allhpi*log(LSE_GE_inputindependent)+2*k_inputindependent;

nu_m1 = exp(log_nu_inputindependent');
var_log_nu_inputindependent = get_var_timeindependent(log_nu_inputindependent, moi_vector, data_allhpi, data.hpi, deathrate_pars, k_dispersion, Ncells, cell_line, model_num);

CI_m1 = get_confidinterval(log_nu_inputindependent',var_log_nu_inputindependent,z_score);


disp(['l = ',num2str(nu_m1(1),'%2.3f'), ', CI[', num2str(CI_m1(1,1),'%2.3f'),', ',num2str(CI_m1(1,2),'%2.3f'), '] ']);
disp(['Least squares error = ',num2str(LSE_GE_inputindependent,'%2.1f')]);  fprintf('\n\n');


% pause;
% ---------------------------------- % ---------------------------------- %


% 'model 2': time-independent, linearly input-dependent rate of virus production
model_num = 2;

linear_legend = 'time-independent, linear input-dependent rate of virus production:';
disp(linear_legend); fprintf('\n');


x0_linear = log(nu_init_timeindependent(2));% r*i
k_linear = length(x0_linear);
% minimize residual sum squared error
[log_nu_linear, LSE_GE_linear] = fminsearch(@(x)getfit_GEmodel(x, moi_vector, data_allhpi, data.hpi, deathrate_pars, k_dispersion, Ncells, cell_line, model_num,which_time_model),x0_linear,options);

AIC_model2 = n_pts_allhpi*log(LSE_GE_linear)+2*k_linear;

nu_m2 = exp(log_nu_linear');
var_log_nu_linear = get_var_timeindependent(log_nu_linear, moi_vector, data_allhpi, data.hpi, deathrate_pars, k_dispersion, Ncells, cell_line, model_num);

CI_m2 = get_confidinterval(log_nu_linear',var_log_nu_linear,z_score);

disp(['r = ',num2str(nu_m2(1),'%2.3f'), ', CI[', num2str(CI_m2(1,1),'%2.3f'),', ',num2str(CI_m2(1,2),'%2.3f'), '] ']);
disp(['Least squares error = ',num2str(LSE_GE_linear,'%2.1f')]);  fprintf('\n\n');


% pause;
% ---------------------------------- % ---------------------------------- %


% 'model 3': time-independent, saturating input-dependent rate of virus production
model_num = 3;

saturating_legend = 'time-independent, saturating input-dependent rate of virus production:';
disp(saturating_legend); fprintf('\n');


x0_saturating = log(nu_init_timeindependent(3:4));% m*i/(i+K)
k_saturating = length(x0_saturating);
% minimize residual sum squared error
[log_nu_saturating, LSE_GE_saturating] = fminsearch(@(x)getfit_GEmodel(x, moi_vector, data_allhpi, data.hpi, deathrate_pars, k_dispersion, Ncells, cell_line, model_num,which_time_model),x0_saturating,options);

AIC_model3 = n_pts_allhpi*log(LSE_GE_saturating)+2*k_saturating;

nu_m3 = exp(log_nu_saturating');
var_log_nu_saturating = get_var_timeindependent(log_nu_saturating, moi_vector, data_allhpi, data.hpi, deathrate_pars, k_dispersion, Ncells, cell_line, model_num);

CI_m3 = get_confidinterval(log_nu_saturating',var_log_nu_saturating,z_score);

if cell_line ==1
    disp(['m = ',num2str(nu_m3(1),'%2.3e'), ', CI[', num2str(CI_m3(1,1),'%2.3e'),', ',num2str(CI_m3(1,2),'%2.3e'), '] ']);
    disp(['K = ',num2str(nu_m3(2),'%2.3e'), ', CI[', num2str(CI_m3(2,1),'%2.3e'),', ',num2str(CI_m3(2,2),'%2.3e'), '] ']);
else
    disp(['m = ',num2str(nu_m3(1),'%2.3f'), ', CI[', num2str(CI_m3(1,1),'%2.3f'),', ',num2str(CI_m3(1,2),'%2.3f'), '] ']);
    disp(['K = ',num2str(nu_m3(2),'%2.3f'), ', CI[', num2str(CI_m3(2,1),'%2.3f'),', ',num2str(CI_m3(2,2),'%2.3f'), '] ']);
end
disp(['Least squares error = ',num2str(LSE_GE_saturating,'%2.1f')]);  fprintf('\n\n');



% pause;
% ---------------------------------- % ---------------------------------- %

% return AIC values to main m-file
AIC_values = [AIC_model1,AIC_model2,AIC_model3];

% ---------------------------------- % ---------------------------------- %
% two plots will appear:
% figure 1 - the model fits vs. MOI (at 6, 12, 18 hpi)
% time-independent models plotted in the bottom panels
% figure 2 - the model fitted production rates, left panel:
% rate vs. input (for 6, 12, 18 hpi); right panel: rate vs. hpi (for i = 1, 4, 9)
% time-independent rates plotted in bottom panels

% order of rgb values: blue, orange, yellow, green
rgb_values = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4660, 0.6740, 0.1880];


figure(1);

for n = 1:length(data.hpi)
    
    subplot(2,3,n);
    semilogy(moi_vector,data_allhpi(n,:),'.','Color',rgb_values(n+1,:),'MarkerSize',30); hold on;
end


for n = 1:length(data.hpi)
    
    this_hpi = data.hpi(n);
    
    % models under best-fit parameter values: virus production at t = 6, 12, 18 hpi
    GE_model_inputindependent(n,:) = get_GEmodel(nu_m1, moi,Ncells,deathrate_pars,k_dispersion,this_hpi,cell_line,1,which_time_model);
    GE_model_linear(n,:) = get_GEmodel(nu_m2, moi,Ncells,deathrate_pars,k_dispersion,this_hpi,cell_line,2,which_time_model);
    GE_model_saturating(n,:) = get_GEmodel(nu_m3, moi,Ncells,deathrate_pars,k_dispersion,this_hpi,cell_line,3,which_time_model);
    
    
    figure(1); subplot(2,3,n);
    if n == 1
        
        p(1) = semilogy(moi,GE_model_inputindependent(n,:),'Color',rgb_values(n+1,:)); hold on;
        p(2) = semilogy(moi,GE_model_linear(n,:),'--','Color',rgb_values(n+1,:)); hold on;
        p(3) = semilogy(moi,GE_model_saturating(n,:),'Color',rgb_values(n+1,:),'LineWidth',3); hold on;
        
    elseif n == 2
        
        q(1) = semilogy(moi,GE_model_inputindependent(n,:),'Color',rgb_values(n+1,:)); hold on;
        q(2) = semilogy(moi,GE_model_linear(n,:),'--','Color',rgb_values(n+1,:)); hold on;
        q(3) = semilogy(moi,GE_model_saturating(n,:),'Color',rgb_values(n+1,:),'LineWidth',3); hold on;
        
    else
        
        r(1) = semilogy(moi,GE_model_inputindependent(n,:),'Color',rgb_values(n+1,:)); hold on;
        r(2) = semilogy(moi,GE_model_linear(n,:),'--','Color',rgb_values(n+1,:)); hold on;
        r(3) = semilogy(moi,GE_model_saturating(n,:),'Color',rgb_values(n+1,:),'LineWidth',3); hold on;
        
    end
    
    title([num2str(data.hpi(n),'%2.0f'),' hpi']);
    axis([0 10 1e3 1e9]);
    xlabel('Bulk MOI'); ylabel('GE/mL');
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    ax.FontName = 'Arial';
    
    
end


% subtitle: time-independent output models
text(-2,1.2,'time-independent virus production models','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'FontWeight','bold', 'FontName', 'Arial')

legend(p,{'input-independent','linear input-dependent','saturating input-dependent'},'Location','SouthEast');
legend(q,{'input-independent','linear input-dependent','saturating input-dependent'},'Location','SouthEast');
legend(r,{'input-independent','linear input-dependent','saturating input-dependent'},'Location','SouthEast');


% ---------------------------------- % ---------------------------------- %

% models: virus production rates vs. input
ilist = 0:10;

t_vary = [6 12 18];
for n=1:length(t_vary)
    t=t_vary(n);
    celloutput_constant(n,:) = nu_m1(1)*ones(size(ilist)); % r*(t-s)
    celloutput_constant(n,1) = 0;
    celloutput_linear(n,:) = nu_m2(1)*ilist; % r*(t-s)*i
    celloutput_MM(n,:) = nu_m3(1)*ilist./(nu_m3(2)+ilist); % m*(t-s)*i/(i+K)
end



% ---------------------------------- % ---------------------------------- %

figure(2);


time = linspace(0,10,101);

% models: virus production rates vs. hpi
viraloutput_inputindependent_rate = nu_m1(1)*ones(1,length(time));
viraloutput_inputindependent_rate(1) = 0;

i_vary = [1 4 8];
for i = 1:length(i_vary)
    this_i = i_vary(i);
    viraloutput_linear_rate(i,:) = nu_m2(1)*this_i*ones(1,length(time));
    viraloutput_saturating_rate(i,:) = nu_m3(1)*this_i/(nu_m3(2)+this_i)*ones(1,length(time));
end

viraloutput_linear_rate(:,1) = 0;
viraloutput_saturating_rate(:,1) = 0;

subplot(1,2,1);
t(1) = plot(time, viraloutput_inputindependent_rate,'Color',rgb_values(1,:)); hold on;
for i = 1:length(i_vary)
    if i==1
        t(2) = plot(time, viraloutput_linear_rate(i,:),'--','Color',rgb_values(1,:)); hold on;
        t(3) = plot(time, viraloutput_saturating_rate(i,:),'Color',rgb_values(1,:),'LineWidth',3); hold on;
    else
        plot(time, viraloutput_linear_rate(i,:),'--','Color',rgb_values(1,:)); hold on;
        plot(time, viraloutput_saturating_rate(i,:),'Color',rgb_values(1,:),'LineWidth',3); hold on;
    end
end
xlabel('time (hpi)'); ylabel('virus production rate');
title('time-independent virus production rates');
if cell_line == 1
axis([0 20 0 5]);
else
axis([0 20 0 0.5]);
end
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
ax.FontName = 'Arial';

legend(t,{'input-independent','linear input-dependent','saturating input-dependent'},'Location','NorthWest');


% % subtitle: time-independent rates
% text(-1,1.2,'time-independent virus production rates','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',18,'FontWeight','bold', 'FontName', 'Arial')
% 
% % label linear model: i = 1, 4, 8
% text(1.01,0.135,'i = 1','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'FontWeight','bold', 'FontName', 'Arial')
% text(1.01,0.45,'i = 4','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'FontWeight','bold', 'FontName', 'Arial')
% text(1.01,0.9,'i = 8','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'FontWeight','bold', 'FontName', 'Arial')
% 
% % label saturating model: i = 1, 4, 8
% text(0.8,0.24,'i = 1','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'FontWeight','bold', 'FontName', 'Arial')
% text(0.8,0.455,'i = 4','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'FontWeight','bold', 'FontName', 'Arial')
% text(0.8,0.59,'i = 8','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'FontWeight','bold', 'FontName', 'Arial')

% ---------------------------------- % ---------------------------------- %


results_timeindependent.AIC_values = AIC_values;

results_timeindependent.bulk_MOI = moi';

results_timeindependent.virus_prod_6hpi = [GE_model_inputindependent(1,:)',GE_model_linear(1,:)',GE_model_saturating(1,:)'];
results_timeindependent.virus_prod_12hpi = [GE_model_inputindependent(2,:)',GE_model_linear(2,:)',GE_model_saturating(2,:)'];
results_timeindependent.virus_prod_18hpi = [GE_model_inputindependent(3,:)',GE_model_linear(3,:)',GE_model_saturating(3,:)'];

results_timeindependent.hpi = time';

results_timeindependent.virusproduction_rates = [viraloutput_inputindependent_rate',viraloutput_linear_rate', viraloutput_saturating_rate'];

results_timeindependent.estime_m1 = nu_m1;
results_timeindependent.estime_m2 = nu_m2;
results_timeindependent.estime_m3 = nu_m3;

results_timeindependent.CI_m1 = CI_m1;
results_timeindependent.CI_m2 = CI_m2;
results_timeindependent.CI_m3 = CI_m3;
