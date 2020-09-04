% function main_celldeath_viraldistribution(void)

% only meant for negative binomial with fixed dispersion parameter


% updated to also fit N cells init - the initial cell population

% updated 08/03/20
% copied from 'cell_death_update072720' -- does not include high MOI of
% FACs when fitting to percent of alive cells that are infected

% update 07/27/20 to include 0-inflated Poisson (suggested by reviewer)

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
% -- want to add 0-inflated Poisson distribution - jdh, 07/27/20

% ---------------------------------- % ---------------------------------- %

clear all; close all; clc;

% 1: MDCK cells; 2: A549 cells
cell_line = 2;

% ---------------------------------- % ---------------------------------- %

% switch over cell line
switch cell_line
    
    case 1 % MDCK cells
        
        
        disp('Fitting cell death models to MDCK cell data...');
        
        % load MDCK data: actual MOI and percent cell death measurements
        %         load('celldeath_MDCK.mat');
        %         load('celldeath_MDCK_cellcount_no0data.mat');
        load('celldeath_MDCK_cellcount_update080420.mat');
        
        % updated FACs data from - 'Martin2019_sourcedata_A549_10202019'
        FACs_data_18hpi = data.FACs_data_18hpi;
        
        FACs_data_18hpi_vector = reshape(FACs_data_18hpi',1,(data.n_datapoints-3));
        
        % model 1:  'a': alpha(t,i) = a
        alpha0_timeindependent_inputindependent = 0.07; % 'a'
        
        % model 2: 'k','b': alpha(t,i) = b*k*t^(k-1)
        alpha0_timedependent_inputindependent = [0.5, 0.3]; % 'k','b'
        
        % model 3: 'c', 'eps': alpha(t,i) = c*i^eps
        alpha0_timeindependent_inputdependent = [0.07, 1e-3]; % 'c', 'eps'
        
        % model 4: 'k','d','eps': alpha(t,i) = i^eps*d*k*t^(k-1)
        alpha0_timedependent_inputdependent = [0.5, 0.3, 1e-3]; % 'k','d','eps'
        
        alpha_init = [alpha0_timeindependent_inputindependent,alpha0_timedependent_inputindependent,alpha0_timeindependent_inputdependent,alpha0_timedependent_inputdependent];
        
        Ncells_init = 250000;
        
        % approximate dispersion parameter of negative binomial:
        %         r_dispersion_init = [0.74,0.63,0.74,0.63]; % models 1,2,3,4
        r_dispersion_fixed = 2
        %         r_dispersion_init = 2
        
    case 2 % A549 cells
        
        
        disp('Fitting cell death models to A549 cell data...');
        
        % load MDCK data: actual MOI and percent cell death measurements
        %         load('celldeath_A549.mat');
        load('celldeath_A549_cellcount_update080420.mat');
        
        % updated FACs data from - 'Martin2019_sourcedata_A549_10202019'
        FACs_data_18hpi = data.FACs_data_18hpi;
        FACs_data_18hpi_vector = reshape(FACs_data_18hpi',1,(data.n_datapoints-3));
        
        % model 1:  'a': alpha(t,i) = a
        alpha0_timeindependent_inputindependent = 0.02; % 'a'
        
        % model 2: 'k','b': alpha(t,i) = b*k*t^(k-1)
        alpha0_timedependent_inputindependent = [1.2,  0.02]; % 'k','b'
        
        % model 3: 'c', 'eps': alpha(t,i) = c*i^eps
        alpha0_timeindependent_inputdependent = [0.03, 1e-3]; % 'c', 'eps'
        
        % model 4: 'k','d','eps': alpha(t,i) = i^eps*d*k*t^(k-1)
        alpha0_timedependent_inputdependent = [1.2, 0.02, 1e-3]; % 'k','d','eps'
        
        alpha_init = [alpha0_timeindependent_inputindependent,alpha0_timedependent_inputindependent,alpha0_timeindependent_inputdependent,alpha0_timedependent_inputdependent];
        
        Ncells_init = 250000;
        
        % approximate dispersion parameter of negative binomial:
        %         r_dispersion_init = [0.35,0.35,0.35,0.35]; % models 1,2,3,4
%                 r_dispersion_fixed = 0.05
%         r_dispersion_fixed = 1
        r_dispersion_fixed = 0.1
        
end


% ---------------------------------- % ---------------------------------- %

% now fit the models with these initial parameter values under Poisson and
% negative binomial distributions to cell death data
% for k = 2 % can change to run Poisson, negative binomial, or zip only
%
%     if k ==1
%
%         % fit Poisson distributed models
%         results_Poisson = submain_celldeath_poisson(data,FACs_data_18hpi_vector,alpha_init,Ncells_init,cell_line);
%
%         if cell_line == 1
%             outfile = 'results_MDCK_Poisson';
%         else
%             outfile = 'results_A549_Poisson';
%         end
%
%         save(outfile,'results_Poisson');
%
%     elseif k == 2

%         pause;
% fit negative binomial distributed models
results_nb = submain_celldeath_negbinomial_fixed(data,FACs_data_18hpi_vector,alpha_init,Ncells_init,r_dispersion_fixed,cell_line);

results_nb.r_dispersion_fixed = r_dispersion_fixed;

if cell_line == 1
    %             outfile = 'results_MDCK_negbinomial';
    %     outfile = 'results_MDCK_negbinomial_rpt2';
    outfile = 'results_MDCK_negbinomial_r2';
else
    %     outfile = 'results_A549_negbinomial';
    %                 outfile = 'results_A549_negbinomial_rpt05';
%     outfile = 'results_A549_negbinomial_r1';
    outfile = 'results_A549_negbinomial_rpt1';
end

save(outfile,'results_nb');

%     else
%
%         %         pause;
%         % fit negative binomial distributed models
%         results_zip = submain_celldeath_zeroinflatedpoisson(data,FACs_data_18hpi_vector,alpha_init,Ncells_init,r_dispersion_init,cell_line);
%
%         if cell_line == 1
%             outfile = 'results_MDCK_zip';
%         else
%             outfile = 'results_A549_zip';
%         end
%
%         save(outfile,'results_zip');
%     end
%
%
%
% end

% pause;

% % for best models under NB
% if 1
%
%     if cell_line == 1
%
%         %     save('results_MDCK_Poisson','results_Poisson');
%         %     save('results_MDCK_negbinomial','results_nb');
%
%
%     else
%
%         %     save('results_A549_Poisson','results_Poisson');
%         %     save('results_A549_negbinomial','results_nb');
%         save('results_A549_zip','results_zip');
%
%     end
%
% end



if 0
    
    % ---------------------------------- % ---------------------------------- %
    
    % collect AIC values from Poisson and negative binomial distribution models
    AIC_values = [results_Poisson.AIC_values,results_nb.AIC_values,results_zip.AIC_values];
    
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
    this_string = [chardelta, 'AIC (zero-inflated Poisson: time-independent, input-independent) = ',num2str(deltaAIC_values(9),'%2.1f')];
    disp(this_string)
    this_string = [chardelta, 'AIC (zero-inflated Poisson: time-dependent, input-independent) = ',num2str(deltaAIC_values(10),'%2.1f')];
    disp(this_string)
    this_string = [chardelta, 'AIC (zero-inflated Poisson: time-independent, input-dependent) = ',num2str(deltaAIC_values(11),'%2.1f')];
    disp(this_string)
    this_string = [chardelta, 'AIC (zero-inflated Poisson: time-dependent, input-dependent) = ',num2str(deltaAIC_values(12),'%2.1f')];
    disp(this_string)
    
end


