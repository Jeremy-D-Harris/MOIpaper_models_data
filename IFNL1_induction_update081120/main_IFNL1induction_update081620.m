% function void = main_IFNL1induction_update081620(void)

% updated 05/11/20

% update on 11/13/19

% main file to estimate IFNL1 induction rates of
% time-independent models, simultaneously fitting to 8 and 18 hpi data
% with the following input-depent models:
% 'model 1': input-independent
% 'model 2': linear input-dependence
% 'model 3': saturating input-dependence
% from these three models, we applied model combination to that data
% from 0-8 hpi and from 8-18 hpi. Altogether, this gave us nine model
% combinations to consider.

% this 'main' file calls the three 'submain' files that do most of the work
% and collects the AIC values are returned to display deltaAIC values
% between the nine model combinations :
% 'submain_IFNL1_inputindependent_plusX':  model 1 (from 0-8 hpi) + models 1-3 (from 8-18 hpi)
% 'submain_IFNL1_linear_plusX': model 2 (from 0-8 hpi) + models 1-3 (from 8-18 hpi)
% 'submain_IFNL1_saturating_plusX': model 3 (from 0-8 hpi) + models 1-3 (from 8-18 hpi)


% ---------------------------------- % ---------------------------------- %

clear all; close all; clc;

% k = 0.6005, CI[0.3620, 0.9960] 
% b = 0.1180, CI[0.0529, 0.2632] 
% mu = 0.0118, CI[0.0083, 0.0169] 
% r (dispersion) = 0.3383, CI[0.2697, 0.4245] 
    
% dispersion parameter of negative binomial
r_dispersion = 0.338330683504483; % as of 08/11/20

% constant death rate best-fit parameter values: [k,d,mu]
deathrate_pars = [0.600487133857002,0.117960454130866,0.0118317417830537]; % as of 08/11/20

% load the data structure that includes MOI and virus titer data
load('IFNL1_8and18hpi.mat');
% updated IFN data from Brigitte - 11/21/19; fewer points below 1!!!!

% model 1 (from 0-8 hpi) + models 1-3 (from 8-18 hpi)
results_inputindependent_plusX = submain_estimate_IFNL1_inputindependent_plusX_nb(data,deathrate_pars,r_dispersion);

AIC_values_inputindependent_plusX = results_inputindependent_plusX.AIC_values;

pause;

% % model 2 (from 0-8 hpi) + models 1-3 (from 8-18 hpi)
results_linear_plusX = submain_estimate_IFNL1_linear_plusX_nb(data,deathrate_pars,r_dispersion);

AIC_values_linear_plusX = results_linear_plusX.AIC_values;

pause;

% model 3 (from 0-8 hpi) + models 1-3 (from 8-18 hpi)
results_saturating_plusX = submain_estimate_IFNL1_saturating_plusX_nb(data,deathrate_pars,r_dispersion);

AIC_values_saturating_plusX = results_saturating_plusX.AIC_values;

% pause;


% concatenate all AIC scores
AIC_values = [AIC_values_inputindependent_plusX, AIC_values_linear_plusX,AIC_values_saturating_plusX];

% get model with the best AIC (lowest value)
min_AIC_value = min(AIC_values);

% deltaAIC values with respect to the best model
deltaAIC_values = AIC_values - min_AIC_value;

% ----------------  model selection - 8 and 18 hpi data  ---------------- %

% to cleanly display deltaAIC values with respect to the best model
chardelta = char(916);
this_string = [chardelta, 'AIC (input-independent+input-independent) = ',num2str(deltaAIC_values(1),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (input-independent+linear input-dependent) = ',num2str(deltaAIC_values(2),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (input-independent+saturating input-dependent) = ',num2str(deltaAIC_values(3),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (linear input-dependent+input-independent) = ',num2str(deltaAIC_values(4),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (linear input-dependent+linear input-dependent) = ',num2str(deltaAIC_values(5),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (linear input-dependent+ saturating input-dependent) = ',num2str(deltaAIC_values(6),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (saturating input-dependent+ input-independent) = ',num2str(deltaAIC_values(7),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (saturating input-dependent+linear input-dependent) = ',num2str(deltaAIC_values(8),'%2.1f')];
disp(this_string)
this_string = [chardelta, 'AIC (saturating input-dependent+saturating input-dependent) = ',num2str(deltaAIC_values(9),'%2.1f')];
disp(this_string)
