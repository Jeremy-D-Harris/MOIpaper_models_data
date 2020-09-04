% function void = main_estimate_virusproduction(void)

% main file for fitting virus production models to data -- updated 02/07/20

% the data include:
% percent virus production (GE/mL) vs. actual MOI (6, 12, 18 hpi)

% for each cell line, MDCK or A549 cells, we fit three models of virus
% production with consideration of time-independent (constant rate over time)
% and time-dependent (piece-wise linear production rate over time)
% rates to virus production data at 6, 12, 18 hpi, simultaneously, by 
% minimizing the residual sum of squares between the model and data

% note that the best model fits of cell death and viral distribution have 
% been brought forward for each cell line -- time-dependent cell death rate
% for MDCK cells and time-independent cell death rate for A549 cells

% six total cell death rate models considered: time-independent (0) and 
% time-dependent (1) models with the following input-depencies:
% 'model 1': input-independent
% 'model 2': linear input-dependence
% 'model 3': saturating input-dependence

% this script: loads data and obtains AIC values from fitting results;
% deltaAIC values are calculated and cleanly displayed for models 1-3
% for time-independent and time-dependent models


% ---------------------------------- % ---------------------------------- %

clear all; close all; clc;

cell_line = 1; % 1: MDCK cells; 2: A549 cells

% which_time_model = 0;

if cell_line == 1 % MDCK cells
    
    disp('Fitting virus production models to MDCK cell data...');
    
    % load the data structure that includes MOI and virus production data
    load('virusproduction_MDCK.mat');
    
    % dispersion parameter of negative binomial
    r_dispersion = 0.597422438658409; % as of 08/11/20
    
    % time-dependent death rate best-fit parameter values
    % Weibull model parameters: [k,b,mu] 
    deathrate_pars = [0.307563294635010,0.447205463276431,0.0145642125211200]; % as of 08/11/20
    
    for which_time_model = 0:1
        
        switch which_time_model
            
            case 0 % time-independent
                
                % constant model: 'l'
                nu_constant = 1.15;
                
                % linear model: 'r'
                nu_linear = 0.37;
                
                % MM model: 'm' 'K'
                nu_saturating = [5.6e14, 15e14]; 
                
                nu_init_timeindependent = [nu_constant,nu_linear,nu_saturating];
                
            case 1  % time-dependent
                
                % constant model: 'l', 's'
                nu_constant = [2, 5];
                
                % linear model: 'r', 's'
                nu_linear = [0.67, 5];
                
                % MM model: 'm', 's', 'K'
                nu_saturating = [1.7, 5.2, 2.6];
                
                nu_init_timedependent = [nu_constant,nu_linear,nu_saturating];
        end
        
    end
    
    
else % A549 cells
    
    disp('Fitting virus production models to A549 cell data...');
    
    % load the data structure that includes MOI and virus production data
    load('virusproduction_A549.mat');
    
    
    % k = 0.6005, CI[0.3620, 0.9960] 
    % b = 0.1180, CI[0.0529, 0.2632] 
    % mu = 0.0118, CI[0.0083, 0.0169] 
    % r (dispersion) = 0.3383, CI[0.2697, 0.4245] 

    % dispersion parameter of negative binomial
    r_dispersion = 0.338330683504483; % as of 08/11/20
    
    % constant death rate best-fit parameter value: [k,d,mu]
    deathrate_pars = [0.600487133857002,0.117960454130866,0.0118317417830537]; % as of 08/11/20
    
    
    for which_time_model = 0:1
        
        
        switch which_time_model
            
            case 0 % time-independent
                
                % constant model: 'l'
                nu_constant = 0.13;
                
                % linear model: 'r'
                nu_linear = 0.044;
                
                % MM model: 'm', 'K'
                nu_saturating = [0.27, 2.06];
                
                nu_init_timeindependent = [nu_constant,nu_linear,nu_saturating];
                
            case 1 % time-dependent
                
                % constant model: 'l', 's'
                nu_constant = [0.17, 5];
                
                % linear model: 'r', 's'
                nu_linear = [0.056, 5];
                
                % MM model: 'm', 's', 'K'
                nu_saturating = [0.3, 5, 1.8];
                
                nu_init_timedependent = [nu_constant,nu_linear,nu_saturating];
                
        end
        
    end
    
end


% ---------------------------------- % ---------------------------------- %


for k = 2 % can change to run time-independent or time-dependent only
    
    if k ==1
        
        % fit time-independent models
        results_timeindependent = submain_virusproduction_timeindependent(data,deathrate_pars,r_dispersion,nu_init_timeindependent,cell_line);
        % AIC_values = results_timeindependent.AIC_values;
        
    else
        
        pause;
        % fit time-dependent models
        results_timedependent = submain_virusproduction_timedependent(data,deathrate_pars,r_dispersion,nu_init_timedependent,cell_line);
        % AIC_values = results_timedependent.AIC_values;
        
    end
    
end

% for best models under NB
if 1
    
    if cell_line == 1
        
    save('results_MDCK_timedependentvirusproduction','results_timedependent');
    
    else
        
    save('results_A549_timedependentvirusproduction','results_timedependent');
        
    end
    
    
end


if 0
% ---------------------------------- % ---------------------------------- %

% concatenate all AIC scores
AIC_values = [results_timeindependent.AIC_values,results_timedependent.AIC_values];

% get model with the best AIC (lowest value)
min_AIC_value = min(AIC_values);

% deltaAIC values with respect to the best model
deltaAIC_values = AIC_values - min_AIC_value;

% ---------------------------------- % ---------------------------------- %


% to cleanly display deltaAIC values with respect to the best model
chardelta = char(916);
this_string = [chardelta, 'AIC (time-independent, input-independent) = ',num2str(deltaAIC_values(1),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (time-independent, linear input-dependent) = ',num2str(deltaAIC_values(2),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (time-independent, saturating input-dependent) = ',num2str(deltaAIC_values(3),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (time-dependent, input-independent) = ',num2str(deltaAIC_values(4),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (time-dependent, linear input-dependent) = ',num2str(deltaAIC_values(5),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (time-dependent, saturating input-dependent) = ',num2str(deltaAIC_values(6),'%2.1f')];
disp(this_string)

end