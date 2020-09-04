function results = submain_celldeath_zeroinflatedpoisson(data,FACs_data_18hpi_vector,alpha_init,Ncells_init,r_dispersion_init,cell_line)

% zero-inflated Poisson
which_distribution = 3;

z_score = 1.96; % for 95% confidence intervals

% submain script: models 1-4
% 'model 1': constant rate of cell death, input-independent
% 'model 2': Weibull hazard function cell death rate, input-independent
% 'model 3': constant rate of cell death (time-independent), input-dependent
% 'model 4': time-dependent + input-dependent

% ---------------------------------- % ---------------------------------- %

% previously calculated Bulk MOI values
actual_moi = data.actual_moi; % (acutal adsorbed virus)/2e6: 6,12,18 hpi
actual_moi_3hpi = data.actual_moi_3hpi; % (acutal adsorbed virus)/2e6: 3 hpi

% reshape treatments into a 1D vector
moi_vector = reshape(actual_moi',[1 data.n_datapoints]);
moi_vector_3hpi = reshape(actual_moi_3hpi',[1 data.n_datapoints]);

moi = linspace(0,10,101); % finer moi array for plotting models



% ---------------------------------- % ---------------------------------- %


if cell_line == 1
    
    fitto_FACs_data_18hpi_vector = FACs_data_18hpi_vector(1:(end-3)); % only use first four replicates
    moi_vector_FACs_data = moi_vector(4:(end-3));
    fig_title = 'MDCK cells';
    
else
    
    fitto_FACs_data_18hpi_vector = FACs_data_18hpi_vector(1:(end-3)); % only use first four replicates
    moi_vector_FACs_data = moi_vector(4:(end-3));
    fig_title = 'A549 cells';
    
end


% ---------------------------------- % ---------------------------------- %


% include fraction cell loss in row vector
data_3hpi = reshape(data.hpi_3hrs',[1, data.n_datapoints]);
data_6hpi = reshape(data.hpi_6hrs',[1, data.n_datapoints]);
data_12hpi = reshape(data.hpi_12hrs',[1, data.n_datapoints]);
data_18hpi = reshape(data.hpi_18hrs',[1, data.n_datapoints]);

% concatenate percent cell death: 3,6,12,18 hpi by row
counts_cellsremaining_data_allhpi = [data_3hpi; data_6hpi; data_12hpi; data_18hpi];
size_data_allhpi = size(counts_cellsremaining_data_allhpi);

% check: total number of data points across all treatments and replicates
n_pts_allhpi_FACs = size_data_allhpi(1)*size_data_allhpi(2)+length(fitto_FACs_data_18hpi_vector);


% ---------------------------------- % ---------------------------------- %

fprintf('\n')
disp('% ---------------------------------- % ---------------------------------- %');
fprintf('\n')

disp('Fitting cell death models under a zero-inflated Poisson (p, prob of extra zeros) to data...')

fprintf('\n')
disp('% ---------------------------------- % ---------------------------------- %');
fprintf('\n\n')

% ---------------------------------- % ---------------------------------- %



% 'model 1': constant rate of cell death, input-independent
model_num = 1;

disp(['time-independent, input-independent cell death rate: model ',num2str(model_num),' (zero-inflated Poisson)']); fprintf('\n');

% constant death rate:
mu_m1 = 0.02;

% model 1: dispersion parameter:
p_extrazeros_m1 = 0.25%r_dispersion_init(1);

% initialize minimization - cell death rate
x0_timeindependent_inputindependent = log([alpha_init(1),mu_m1,Ncells_init,p_extrazeros_m1]); % 'a': alpha(t,i) = a

k_timeindependent_inputindependent = length(x0_timeindependent_inputindependent);

% options = optimset('Display','iter');
% minimize sum squared error to get estimates of log(a)
% [log_alpha_dispersion_timeindependent_inputindependent, LSE_percent_celldeath_timeindependent_inputindependent] = fminsearch(@(x)getfit_deathrate_dispersion(x,moi_vector_3hpi,moi_vector,moi_vector_FACs_data,counts_cellsremaining_data_allhpi,fitto_FACs_data_18hpi_vector,data.hpi,model_num,which_distribution),x0_timeindependent_inputindependent,options);
[log_alpha_dispersion_timeindependent_inputindependent, LSE_percent_celldeath_timeindependent_inputindependent] = fminsearch(@(x)getfit_deathrate_dispersion(x,moi_vector_3hpi,moi_vector,moi_vector_FACs_data,counts_cellsremaining_data_allhpi,fitto_FACs_data_18hpi_vector,data.hpi,model_num,which_distribution),x0_timeindependent_inputindependent);

AIC_model1 = n_pts_allhpi_FACs*log(LSE_percent_celldeath_timeindependent_inputindependent)+2*k_timeindependent_inputindependent;

var_log_alpha_m1 = get_var_zeroinflatedpoisson(log_alpha_dispersion_timeindependent_inputindependent, moi_vector_3hpi, moi_vector, moi_vector_FACs_data,data.hpi, counts_cellsremaining_data_allhpi, fitto_FACs_data_18hpi_vector, model_num, cell_line);
alpha_m1 = exp(log_alpha_dispersion_timeindependent_inputindependent');
p_extrazeros_m1_estimate = alpha_m1(end);

CI_m1 = get_confidinterval(log_alpha_dispersion_timeindependent_inputindependent',var_log_alpha_m1,z_score);

disp(['a = ',num2str(alpha_m1(1),'%2.4f'), ', CI[', num2str(CI_m1(1,1),'%2.4f'),', ',num2str(CI_m1(1,2),'%2.4f'), '] ']);
disp(['mu = ',num2str(alpha_m1(2),'%2.4f'), ', CI[', num2str(CI_m1(2,1),'%2.4f'),', ',num2str(CI_m1(2,2),'%2.4f'), '] ']);
disp(['N cells init = ',num2str(alpha_m1(3),'%1.3e'), ', CI[', num2str(CI_m1(3,1),'%1.3e'),', ',num2str(CI_m1(3,2),'%1.3e'), '] ']);
disp(['probability extra zeros = ',num2str(p_extrazeros_m1_estimate,'%2.4f'), ', CI[', num2str(CI_m1(4,1),'%2.4f'),', ',num2str(CI_m1(4,2),'%2.4f'), '] ']);
disp(['Least squares error = ',num2str(LSE_percent_celldeath_timeindependent_inputindependent,'%2.0f')]);  fprintf('\n\n');


% pause;
% ---------------------------------- % ---------------------------------- %

% 'model 2': Weibull hazard function cell death rate, input-independent
model_num = 2;

disp(['time-dependent, input-independent cell death rate: model ',num2str(model_num),' (zero-inflated Poisson)']); fprintf('\n');

mu_m2 = 0.02;

% model 2: dispersion parameter:
p_extrazeros_m2 = 0.3%r_dispersion_init(2);

% initialize minimization - shape and scale parameters of the Weibull hazard function
x0_timedependent_inputindependent = log([alpha_init(2:3),mu_m2,Ncells_init,p_extrazeros_m2]); % 'k','b': alpha(t,i) = b*k*t^(k-1)

k_timedependent_inputindependent = length(x0_timedependent_inputindependent);

% options = optimset('Display','iter');
% minimize sum squared error to get estimates of log(k) and log(b)
% [log_alpha_dispersion_timedependent_inputindependent, LSE_percent_celldeath_timedependent_inputindependent] = fminsearch(@(x)getfit_deathrate_dispersion(x,moi_vector_3hpi,moi_vector,moi_vector_FACs_data,counts_cellsremaining_data_allhpi,fitto_FACs_data_18hpi_vector,data.hpi,model_num,which_distribution),x0_timedependent_inputindependent,options);
[log_alpha_dispersion_timedependent_inputindependent, LSE_percent_celldeath_timedependent_inputindependent] = fminsearch(@(x)getfit_deathrate_dispersion(x,moi_vector_3hpi,moi_vector,moi_vector_FACs_data,counts_cellsremaining_data_allhpi,fitto_FACs_data_18hpi_vector,data.hpi,model_num,which_distribution),x0_timedependent_inputindependent);

AIC_model2 = n_pts_allhpi_FACs*log(LSE_percent_celldeath_timedependent_inputindependent)+2*k_timedependent_inputindependent;

var_log_alpha_m2 = get_var_zeroinflatedpoisson(log_alpha_dispersion_timedependent_inputindependent, moi_vector_3hpi, moi_vector, moi_vector_FACs_data,data.hpi, counts_cellsremaining_data_allhpi, fitto_FACs_data_18hpi_vector, model_num, cell_line);
alpha_m2 = exp(log_alpha_dispersion_timedependent_inputindependent');
p_extrazeros_m2_estimate = alpha_m2(end);

CI_m2 = get_confidinterval(log_alpha_dispersion_timedependent_inputindependent',var_log_alpha_m2,z_score);

disp(['k = ',num2str(alpha_m2(1),'%2.4f'), ', CI[', num2str(CI_m2(1,1),'%2.4f'),', ',num2str(CI_m2(1,2),'%2.4f'), '] ']);
disp(['b = ',num2str(alpha_m2(2),'%2.4f'), ', CI[', num2str(CI_m2(2,1),'%2.4f'),', ',num2str(CI_m2(2,2),'%2.4f'), '] ']);
disp(['mu = ',num2str(alpha_m2(3),'%2.4f'), ', CI[', num2str(CI_m2(3,1),'%2.4f'),', ',num2str(CI_m2(3,2),'%2.4f'), '] ']);
disp(['N cells init = ',num2str(alpha_m2(4),'%1.3e'), ', CI[', num2str(CI_m2(4,1),'%1.3e'),', ',num2str(CI_m2(4,2),'%1.3e'), '] ']);
disp(['probability extra zeros  = ',num2str(p_extrazeros_m2_estimate,'%2.4f'), ', CI[', num2str(CI_m2(5,1),'%2.4f'),', ',num2str(CI_m2(5,2),'%2.4f'), '] ']);
disp(['Least squares error = ',num2str(LSE_percent_celldeath_timedependent_inputindependent,'%2.0f')]);  fprintf('\n\n');


% pause;
% ---------------------------------- % ---------------------------------- %


% 'model 3': constant rate of cell death (time-independent), input-dependent
model_num = 3;

disp(['time-independent, input-dependent cell death rate: model ',num2str(model_num),' (zero-inflated Poisson)']); fprintf('\n');

mu_m3 = 0.02;

% model 3: dispersion parameter:
p_extrazeros_m3 = 0.3%r_dispersion_init(3);

% initialize minimization - per capita death rate
x0_timeindependent_inputdependent = log([alpha_init(4:5),mu_m3, Ncells_init,p_extrazeros_m3]); % 'c': alpha(t,i) = c*i^eps

k_timeindependent_inputdependent = length(x0_timeindependent_inputdependent);

% minimize sum squared error to get estimates of log(c)
% [log_alpha_dispersion_timeindependent_inputdependent, LSE_percent_celldeath_timeindependent_inputdependent] = fminsearch(@(x)getfit_deathrate_dispersion(x,moi_vector_3hpi,moi_vector,moi_vector_FACs_data,counts_cellsremaining_data_allhpi,fitto_FACs_data_18hpi_vector,data.hpi,model_num,which_distribution),x0_timeindependent_inputdependent,options);
[log_alpha_dispersion_timeindependent_inputdependent, LSE_percent_celldeath_timeindependent_inputdependent] = fminsearch(@(x)getfit_deathrate_dispersion(x,moi_vector_3hpi,moi_vector,moi_vector_FACs_data,counts_cellsremaining_data_allhpi,fitto_FACs_data_18hpi_vector,data.hpi,model_num,which_distribution),x0_timeindependent_inputdependent);

AIC_model3 = n_pts_allhpi_FACs*log(LSE_percent_celldeath_timeindependent_inputdependent)+2*k_timeindependent_inputdependent;

var_log_alpha_m3 = get_var_zeroinflatedpoisson(log_alpha_dispersion_timeindependent_inputdependent, moi_vector_3hpi, moi_vector, moi_vector_FACs_data,data.hpi, counts_cellsremaining_data_allhpi, fitto_FACs_data_18hpi_vector, model_num, cell_line);
alpha_m3 = exp(log_alpha_dispersion_timeindependent_inputdependent');
p_extrazeros_m3_estimate = alpha_m3(end);

CI_m3 = get_confidinterval(log_alpha_dispersion_timeindependent_inputdependent',var_log_alpha_m3,z_score);

disp(['c = ',num2str(alpha_m3(1),'%2.4f'), ', CI[', num2str(CI_m3(1,1),'%2.4f'),', ',num2str(CI_m3(1,2),'%2.4f'), '] ']);
disp(['epsilon = ',num2str(alpha_m3(2),'%2.3e'), ', CI[', num2str(CI_m3(2,1),'%2.3e'),', ',num2str(CI_m3(2,2),'%2.3e'), '] ']);
disp(['mu = ',num2str(alpha_m3(3),'%2.4f'), ', CI[', num2str(CI_m3(3,1),'%2.4f'),', ',num2str(CI_m3(3,2),'%2.4f'), '] ']);
disp(['N cells init = ',num2str(alpha_m3(4),'%1.3e'), ', CI[', num2str(CI_m3(4,1),'%1.3e'),', ',num2str(CI_m3(4,2),'%1.3e'), '] ']);
disp(['probability extra zeros = ',num2str(p_extrazeros_m3_estimate,'%2.4f'), ', CI[', num2str(CI_m3(5,1),'%2.4f'),', ',num2str(CI_m3(5,2),'%2.4f'), '] ']);
disp(['Least squares error = ',num2str(LSE_percent_celldeath_timeindependent_inputdependent,'%2.0f')]);  fprintf('\n\n');


% pause;
% ---------------------------------- % ---------------------------------- %


% 'model 4': time-dependent + input-dependent
model_num = 4;

disp(['time-dependent, input-dependent cell death rate: model ',num2str(model_num),' (zero-inflated Poisson)']); fprintf('\n');

mu_m4 = 0.02;

% model 4: dispersion parameter:
p_extrazeros_m4 = 0.3 %r_dispersion_init(4);

% minimize sum squared error to get estimates of 'k' and 'd'
% initialize minimization - per capita death rate
x0_timedependent_inputdependent = log([alpha_init(6:8),mu_m4,Ncells_init,p_extrazeros_m4]); % 'k','d': alpha(t,i) = b4*k4*t^(k4-1)*i^eps

k_timedependent_inputdependent = length(x0_timedependent_inputdependent);

% options = optimset('Display','iter','MaxFunEvals', 1000);
% minimize sum squared error to get estimate of per capita cell death rate
% [log_alpha_dispersion_timedependent_inputdependent, LSE_percent_celldeath_timedependent_inputdependent] = fminsearch(@(x)getfit_deathrate_dispersion(x,moi_vector_3hpi,moi_vector,moi_vector_FACs_data,counts_cellsremaining_data_allhpi,fitto_FACs_data_18hpi_vector,data.hpi,model_num,which_distribution),x0_timedependent_inputdependent,options);
[log_alpha_dispersion_timedependent_inputdependent, LSE_percent_celldeath_timedependent_inputdependent] = fminsearch(@(x)getfit_deathrate_dispersion(x,moi_vector_3hpi,moi_vector,moi_vector_FACs_data,counts_cellsremaining_data_allhpi,fitto_FACs_data_18hpi_vector,data.hpi,model_num,which_distribution),x0_timedependent_inputdependent);

AIC_model4 = n_pts_allhpi_FACs*log(LSE_percent_celldeath_timedependent_inputdependent)+2*k_timedependent_inputdependent;

var_log_alpha_m4 = get_var_zeroinflatedpoisson(log_alpha_dispersion_timedependent_inputdependent, moi_vector_3hpi, moi_vector, moi_vector_FACs_data,data.hpi, counts_cellsremaining_data_allhpi, fitto_FACs_data_18hpi_vector, model_num, cell_line);
alpha_m4 = exp(log_alpha_dispersion_timedependent_inputdependent');
p_extrazeros_m4_estimate = alpha_m4(end);

CI_m4 = get_confidinterval(log_alpha_dispersion_timedependent_inputdependent',var_log_alpha_m4,z_score);


disp(['k = ',num2str(alpha_m4(1),'%2.4f'), ', CI[', num2str(CI_m4(1,1),'%2.4f'),', ',num2str(CI_m4(1,2),'%2.4f'), '] ']);
disp(['d = ',num2str(alpha_m4(2),'%2.4f'), ', CI[', num2str(CI_m4(2,1),'%2.4f'),', ',num2str(CI_m4(2,2),'%2.4f'), '] ']);
disp(['eps = ',num2str(alpha_m4(3),'%2.3e'), ', CI[', num2str(CI_m4(3,1),'%2.3e'),', ',num2str(CI_m4(3,2),'%2.3e'), '] ']);
disp(['mu = ',num2str(alpha_m4(4),'%2.4f'), ', CI[', num2str(CI_m4(4,1),'%2.4f'),', ',num2str(CI_m4(4,2),'%2.4f'), '] ']);
disp(['N cells init = ',num2str(alpha_m4(5),'%1.3e'), ', CI[', num2str(CI_m4(5,1),'%1.3e'),', ',num2str(CI_m4(5,2),'%1.3e'), '] ']);
disp(['probability extra zeros = ',num2str(p_extrazeros_m4_estimate,'%2.4f'), ', CI[', num2str(CI_m4(6,1),'%2.4f'),', ',num2str(CI_m4(6,2),'%2.4f'), '] ']);
disp(['Least squares error = ',num2str(LSE_percent_celldeath_timedependent_inputdependent,'%2.0f')]);  fprintf('\n\n');



% pause;
% ---------------------------------- % ---------------------------------- %

% three plots will appear:
% figure 3 - model fits to percent cells remaining vs. Bulk MOI
% figure 4 - cell death rate vs. hpi
% figure 5 (bottom panel) - model fits to percent of alive cells that are infected vs. Bulk MOI

% color order of RGB values: blue, orange, yellow, green (darker)
rgb_values = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4660, 0.6740, 0.1880];


figure(5);
% plot percent cells remaining vs. MOI (3 6 12 18 hpi)
for k = 1:length(data.hpi)
    
    if k == 1
        % plot percent cells remaining data at all hours post infection
        subplot(1,4,k); plot(moi_vector_3hpi,counts_cellsremaining_data_allhpi(k,:),'.','Color',rgb_values(k,:),'MarkerSize',30); hold on;
    else
        % plot percent cells remaining data at all hours post infection
        subplot(1,4,k); plot(moi_vector,counts_cellsremaining_data_allhpi(k,:),'.','Color',rgb_values(k,:),'MarkerSize',30); hold on;
    end
end

Ncells_init_m1 = alpha_m1(end-1);
Ncells_init_m2 = alpha_m2(end-1);
Ncells_init_m3 = alpha_m3(end-1);
Ncells_init_m4 = alpha_m4(end-1);

for k = 1:length(data.hpi)
    
    this_hpi = data.hpi(k);
    
    % models: percent cells remaining at t = 3, 6, 12, 18 hpi
    m1_cellsremaining = (Ncells_init_m1/100)*get_percentalive(alpha_m1,data.hpi(k),moi,p_extrazeros_m1_estimate,1,which_distribution); % last two: which model, which distribution
    m2_cellsremaining = (Ncells_init_m2/100)*get_percentalive(alpha_m2,data.hpi(k),moi,p_extrazeros_m2_estimate,2,which_distribution);
    m3_cellsremaining = (Ncells_init_m3/100)*get_percentalive(alpha_m3,data.hpi(k),moi,p_extrazeros_m3_estimate,3,which_distribution);
    m4_cellsremaining = (Ncells_init_m4/100)*get_percentalive(alpha_m4,data.hpi(k),moi,p_extrazeros_m4_estimate,4,which_distribution);
    
    
    % fill in percent remaining (increasing MOI going down row) for each hpi (columns)
    m1_allhpi(:,k) = m1_cellsremaining';
    m2_allhpi(:,k) = m2_cellsremaining';
    m3_allhpi(:,k) = m3_cellsremaining';
    m4_allhpi(:,k) = m4_cellsremaining';
    
    figure(5); subplot(1,4,k);
    if k ==1
        
        r(1) = plot(moi,m1_cellsremaining,'--','Color',rgb_values(k,:),'LineWidth',1); hold on;
        r(2) = plot(moi,m2_cellsremaining,'Color',rgb_values(k,:),'LineWidth',1); hold on;
        r(3) = plot(moi,m3_cellsremaining,'--','Color',rgb_values(k,:),'LineWidth',3); hold on;
        r(4) = plot(moi,m4_cellsremaining,'Color',rgb_values(k,:),'LineWidth',3); hold on;
        
    else
        
        plot(moi,m1_cellsremaining,'--','Color',rgb_values(k,:),'LineWidth',1); hold on;
        plot(moi,m2_cellsremaining,'Color',rgb_values(k,:),'LineWidth',1); hold on;
        plot(moi,m3_cellsremaining,'--','Color',rgb_values(k,:),'LineWidth',3); hold on;
        plot(moi,m4_cellsremaining,'Color',rgb_values(k,:),'LineWidth',3); hold on;
    end
    title([num2str(this_hpi,'%2.0f'),' hpi']);
    axis([0 10 0 2.5e5]);
    xlabel('Bulk MOI'); ylabel('Number of cells remaining');
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    ax.FontName = 'Arial';
    
    
    % include figure 3 to compare the FACs data with model predicted from the
    % 18hpi time point data -- MDCK or A549
    % easier to put inside the loop
    if k == 4 % 18 hpi time point
        
        % now, can we try to compare the FACs data with a negative-binomial
        % model instead of a Poisson model of infection distribution across
        
        
        % dispersion parameter for negative-binomial
        m1_percentremaining_infected_zip = get_percentaliveinfected(alpha_m1,data.hpi(k),moi,p_extrazeros_m3_estimate,1,which_distribution); % which model, which distirbution
        m2_percentremaining_infected_zip = get_percentaliveinfected(alpha_m2,data.hpi(k),moi,p_extrazeros_m3_estimate,2,which_distribution);
        m3_percentremaining_infected_zip = get_percentaliveinfected(alpha_m3,data.hpi(k),moi,p_extrazeros_m3_estimate,3,which_distribution);
        m4_percentremaining_infected_zip = get_percentaliveinfected(alpha_m4,data.hpi(k),moi,p_extrazeros_m3_estimate,4,which_distribution);
        
        figure(7); subplot(3,1,3);
        
%         fitto_FACs_data_18hpi_vector = FACs_data_18hpi_vector(1:(end-3)); % only use first four replicates
%         moi_vector_FACs_data = moi_vector(1:(end-3));
        
        
        q(1) = plot(moi_vector_FACs_data,fitto_FACs_data_18hpi_vector,'ko','MarkerSize',12); hold on;
        q(2) = plot(moi_vector((end-2):end),FACs_data_18hpi_vector((end-2):end),'ro','MarkerSize',12); hold on;
        q(3) = plot(moi,m1_percentremaining_infected_zip,'k--','LineWidth',1); hold on;
        q(4) = plot(moi,m2_percentremaining_infected_zip,'k','LineWidth',1); hold on;
        q(5) = plot(moi,m3_percentremaining_infected_zip,'k--','LineWidth',3); hold on;
        q(6) = plot(moi,m4_percentremaining_infected_zip,'k','LineWidth',3); hold on;
        axis([0 10 0 100]);
        xlabel('Bulk MOI'); ylabel('Percent of alive cells that are infected');
        title([fig_title,': zero-inflated Poisson']);
        ax = gca;
        ax.FontSize = 10;
        ax.FontWeight = 'bold';
        ax.FontName = 'Arial';
    end
    
    
end

figure(5);
sp = suptitle([fig_title,': zero-inflated Poisson']);
sp.FontSize = 12;
sp.FontWeight = 'Bold';

legend(r,{'time-independent, input-independent','time-dependent, input-independent','time-independent, input-dependent','time-dependent, input-dependent'},'Location','SouthWest'); %[0.25 0.65 .1 .05]
% legend(q,{'FACs 18 hpi (included)','FACs 18 hpi (excluded)','time-independent, input-independent','time-dependent, input-independent','time-independent, input-dependent','time-dependent, input-dependent'},'Location','SouthEast');

% create a finer time array
time_hpi = linspace(0,18,1000);

% model 1: time-independent, input-independent death rate = a
cell_deathrate_m1 = get_celldeathrate_model(alpha_m1,time_hpi,1,1);
uninfectedcell_deathrate_m1 = alpha_m1(2)*ones(1,length(time_hpi));

% model 2: time-dependent, input-independent death rate = b*k*t^(k-1)
cell_deathrate_m2 = get_celldeathrate_model(alpha_m2,time_hpi,2,1);
uninfectedcell_deathrate_m2 = alpha_m2(3)*ones(1,length(time_hpi));

% model 3: time-independent, input-dependent death rate = c*i^eps
uninfectedcell_deathrate_m3 = alpha_m3(3)*ones(1,length(time_hpi));

% model 4: time-dependent, input-dependent death rate = b*k*t^(k-1)*i^eps
uninfectedcell_deathrate_m4 = alpha_m4(4)*ones(1,length(time_hpi));

i_vary = [1 4 8];
for n=1:length(i_vary)
    this_i = i_vary(n);
    cell_deathrate_m3(n,:) = get_celldeathrate_model(alpha_m3,time_hpi,3,this_i);
    cell_deathrate_m4(n,:) = get_celldeathrate_model(alpha_m4,time_hpi,4,this_i);
end
cell_deathrate_m3(:,1) = 0;
cell_deathrate_m4(:,1) = 0;

figure(6); subplot(1,2,1);
% plot cell death rate vs. time (hpi)
p(1) = plot(time_hpi,cell_deathrate_m1,'--','Color',rgb_values(1,:),'LineWidth',1); hold on;
p(2) = plot(time_hpi,cell_deathrate_m2,'Color',rgb_values(1,:),'LineWidth',1); hold on;
p(3) = plot(time_hpi,cell_deathrate_m3(1,:),'--','Color',rgb_values(1,:),'LineWidth',3); hold on;
plot(time_hpi,cell_deathrate_m3(2,:),'--','Color',rgb_values(1,:),'LineWidth',3); hold on;
plot(time_hpi,cell_deathrate_m3(3,:),'--','Color',rgb_values(1,:),'LineWidth',3); hold on;
p(4) = plot(time_hpi,cell_deathrate_m4(1,:),'Color',rgb_values(1,:),'LineWidth',3); hold on;
plot(time_hpi,cell_deathrate_m4(2,:),'Color',rgb_values(1,:),'LineWidth',3); hold on;
plot(time_hpi,cell_deathrate_m4(3,:),'Color',rgb_values(1,:),'LineWidth',3); hold on;
xlabel('time (hpi)'); ylabel('Cell death rate (per hour)');
title([fig_title, ': zero-inflated Poisson']);
if cell_line == 1
axis([0 18 0 0.5]);
else
axis([0 18 0 0.5]);
end
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
ax.FontName = 'Arial';
% legend(p,{'time-independent, input-independent','time-dependent, input-independent','time-independent, input-dependent'},'Location','NorthEast');
legend(p,{'time-independent, input-independent','time-dependent, input-independent','time-independent, input-dependent','time-dependent, input-dependent'},'Location','NorthEast');
% labels (i = 1, 4, 8) for models 3,4
% text(0.8,0.15,['i = ',num2str(i_vary(1))],'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',18,'FontWeight','bold', 'FontName', 'Arial')
% text(0.8,0.35,['i = ',num2str(i_vary(2))],'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',18,'FontWeight','bold', 'FontName', 'Arial')
% text(0.8,0.65,['i = ',num2str(i_vary(3))],'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',18,'FontWeight','bold', 'FontName', 'Arial')

figure(6); subplot(1,2,2);
% plot cell death rate vs. time (hpi)
q(1) = plot(time_hpi,uninfectedcell_deathrate_m1,'--','Color',rgb_values(3,:),'LineWidth',1); hold on;
q(2) = plot(time_hpi,uninfectedcell_deathrate_m2,'Color',rgb_values(3,:),'LineWidth',1); hold on;
q(3) = plot(time_hpi,uninfectedcell_deathrate_m3,'--','Color',rgb_values(3,:),'LineWidth',3); hold on;
q(4) = plot(time_hpi,uninfectedcell_deathrate_m4,'Color',rgb_values(3,:),'LineWidth',3); hold on;
xlabel('time (hpi)'); ylabel('Uninfected cell death rate (per hour)');
title([fig_title, ': zero-inflated Poisson']);
if cell_line == 1
axis([0 18 0 0.02]);
else
axis([0 18 0 0.02]);
end
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
ax.FontName = 'Arial';


% best parameter estimates
results.estimates_m1 = alpha_m1;
results.estimates_m2 = alpha_m2;
results.estimates_m3 = alpha_m3;
results.estimates_m4 = alpha_m4;

% corresponding confidence intervals
results.CIs_m1 = CI_m1;
results.CIs_m2 = CI_m2;
results.CIs_m3 = CI_m3;
results.CIs_m4 = CI_m4;

% collect results and create structure
results.hpi = time_hpi';
results.death_rate = [cell_deathrate_m1',cell_deathrate_m2',cell_deathrate_m3',cell_deathrate_m4'];

% percent remaining - each hpi
results.bulk_MOI = moi';
results.percentcells_remaining_3hpi = [m1_allhpi(:,1),m2_allhpi(:,1),m3_allhpi(:,1),m4_allhpi(:,1)];
results.percentcells_remaining_6hpi = [m1_allhpi(:,2),m2_allhpi(:,2),m3_allhpi(:,2),m4_allhpi(:,2)];
results.percentcells_remaining_12hpi = [m1_allhpi(:,3),m2_allhpi(:,3),m3_allhpi(:,3),m4_allhpi(:,3)];
results.percentcells_remaining_18hpi = [m1_allhpi(:,4),m2_allhpi(:,4),m3_allhpi(:,4),m4_allhpi(:,4)];

% percent remaining and alive
results.percentcells_alive_infected = [m1_percentremaining_infected_zip',m2_percentremaining_infected_zip',m3_percentremaining_infected_zip',m4_percentremaining_infected_zip'];

% collect and return AIC values
results.AIC_values = [AIC_model1,AIC_model2,AIC_model3,AIC_model4];

