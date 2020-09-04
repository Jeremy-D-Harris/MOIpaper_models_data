% function main_celldeath_viraldistribution(void)

% main file for fitting cell death and viral distribution 
% models to data -- updated 02/07/20

% the data include:
% percent cells surviving vs. actual MOI (3, 6, 12, 18 hpi);
% percent of alive cells that are infected vs. actual MOI (18 hpi)

% for each cell line, MDCK or A549 cells, fit four models of cell death
% rates to cell death and FACs data (18 hpi), simultaneously, by minimizing
% the residual sum of squares between the model and data

% four cell death models considered:
% 'model 1': constant rate of cell death (time-independent), input-independent
% 'model 2': Weibull hazard function cell death rate (time-dependent), input-independent
% 'model 3': constant rate of cell death (time-independent), input-dependent
% 'model 4': time-dependent, input-dependent

% this script: loads data and obtains AIC values from fitting results;
% deltaAIC values are calculated and cleanly displayed for models 1-4 under
% Poisson and negative binomial distributions

% ---------------------------------- % ---------------------------------- %

clear all; close all; clc;

% 1: MDCK cells; 2: A549 cells
cell_line = 1;

% ---------------------------------- % ---------------------------------- %

% switch over cell line
switch cell_line
    
    case 1 % MDCK cells
        
        
        disp('Fitting cell death models to MDCK cell data...');
        
        % load MDCK data: actual MOI and percent cell death measurements
        load('celldeath_MDCK.mat');
        
        % updated FACs data from - 'Martin2019_sourcedata_A549_10202019'
        FACs_data_18hpi = data.FACs_data_18hpi; 
        
        FACs_data_18hpi_vector = reshape(FACs_data_18hpi',1,data.n_datapoints);
        
        % model 1:  'a': alpha(t,i) = a
        alpha0_timeindependent_inputindependent = 0.07; % 'a'
        
        % model 2: 'k','b': alpha(t,i) = b*k*t^(k-1)
        alpha0_timedependent_inputindependent = [0.5, 0.3]; % 'k','b'
        
        % model 3: 'c', 'eps': alpha(t,i) = c*i^eps
        alpha0_timeindependent_inputdependent = [0.07, 1e-3]; % 'c', 'eps'
        
        % model 4: 'k','d','eps': alpha(t,i) = i^eps*d*k*t^(k-1)
        alpha0_timedependent_inputdependent = [0.5, 0.3, 1e-3]; % 'k','d','eps'
        
        alpha_init = [alpha0_timeindependent_inputindependent,alpha0_timedependent_inputindependent,alpha0_timeindependent_inputdependent,alpha0_timedependent_inputdependent];
        
        % approximate dispersion parameter of negative binomial:
        r_dispersion_init = [0.74,0.63,0.74,0.63]; % models 1,2,3,4
        
        
    case 2 % A549 cells
        
        
        disp('Fitting cell death models to A549 cell data...');
        
        % load MDCK data: actual MOI and percent cell death measurements
        load('celldeath_A549.mat');
        
        % updated FACs data from - 'Martin2019_sourcedata_A549_10202019'
        FACs_data_18hpi = data.FACs_data_18hpi;
        FACs_data_18hpi_vector = reshape(FACs_data_18hpi',1,data.n_datapoints);
        
        % model 1:  'a': alpha(t,i) = a
        alpha0_timeindependent_inputindependent = 0.02; % 'a'
        
        % model 2: 'k','b': alpha(t,i) = b*k*t^(k-1)
        alpha0_timedependent_inputindependent = [1.2,  0.02]; % 'k','b'
        
        % model 3: 'c', 'eps': alpha(t,i) = c*i^eps
        alpha0_timeindependent_inputdependent = [0.03, 1e-3]; % 'c', 'eps'
        
        % model 4: 'k','d','eps': alpha(t,i) = i^eps*d*k*t^(k-1)
        alpha0_timedependent_inputdependent = [1.2, 0.02, 1e-3]; % 'k','d','eps'
        
        alpha_init = [alpha0_timeindependent_inputindependent,alpha0_timedependent_inputindependent,alpha0_timeindependent_inputdependent,alpha0_timedependent_inputdependent];
        
        % approximate dispersion parameter of negative binomial:
        r_dispersion_init = [0.35,0.35,0.35,0.35]; % models 1,2,3,4
        
        
end


% ---------------------------------- % ---------------------------------- %

% now fit the models with these initial parameter values under Poisson and
% negative binomial distributions to cell death data
for k = 1:2 % can change to run Poisson or negative binomial only
    
    if k ==1
        
        % fit Poisson distributed models
        results_Poisson = submain_celldeath_poisson(data,FACs_data_18hpi_vector,alpha_init,cell_line);
        
    else
        
        pause;
        % fit negative binomial distributed models
        results_nb = submain_celldeath_negbinomial(data,FACs_data_18hpi_vector,alpha_init,r_dispersion_init,cell_line);
        
    end
    
end

% ---------------------------------- % ---------------------------------- %

% collect AIC values from Poisson and negative binomial distribution models
AIC_values = [results_Poisson.AIC_values,results_nb.AIC_values];

% best AIC value = lowest value
min_AIC_value = min(AIC_values);

% deltaAIC values with respect to the best model
deltaAIC_values = AIC_values - min_AIC_value;

% ---------------------------------- % ---------------------------------- %

% to cleanly display deltaAIC values with respect to the best model
chardelta = char(916);
this_string = [chardelta, 'AIC (Poisson: time-independent, input-independent) = ',num2str(deltaAIC_values(1),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (Poisson: time-dependent, input-independent) = ',num2str(deltaAIC_values(2),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (Poisson: time-independent, input-dependent) = ',num2str(deltaAIC_values(3),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (Poisson: time-dependent, input-dependent) = ',num2str(deltaAIC_values(4),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (Negative binomial: time-independent, input-independent) = ',num2str(deltaAIC_values(5),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (Negative binomial: time-dependent, input-independent) = ',num2str(deltaAIC_values(6),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (Negative binomial: time-independent, input-dependent) = ',num2str(deltaAIC_values(7),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (Negative binomial: time-dependent, input-dependent) = ',num2str(deltaAIC_values(8),'%2.1f')];
disp(this_string)

