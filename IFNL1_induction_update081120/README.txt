updated 04/11/20 -- IFN induction model fits to data

We fit IFNL1 induction rates of time-independent models to 8 and 18 hpi
data simultaneously with the following input-depent models: 
'model 1': input-independent 
'model 2': linear input-dependence 
'model 3': saturating input-dependence 

From these three models, we applied model combinations to that data from
0-8 hpi and from 8-18 hpi. Altogether, this gave us 9 model combinations
to consider. 

The data for IFNL1 induction in A549 cells are contained in the .mat
file "IFNL1_8and18hpi." The actual MOI values and the IFNL1 fold change
values can be found in the source data file under the tab "Figure 6C-F
(A549)"

This folder contains the code that fits 9 model combinations of
IFNL1 induction in A549 cells: the 'main' file calls the three 'submain'
files that do most of the work. To run the code, open the m-file
"main_IFNL1induction" (main) in Matlab, and press 'Run.' This script
loads the data and passes it to the three submain m-files:

'submain_IFNL1_inputindependent_plusX':  
model 1 (from 0-8 hpi) + models 1-3 (from 8-18 hpi) 

'submain_IFNL1_linear_plusX': 
model 2 (from 0-8 hpi) + models 1-3 (from 8-18 hpi)

'submain_IFNL1_saturating_plusX': 
model 3 (from 0-8 hpi) + models 1-3 (from 8-18 hpi) 


Within each of the submain files, the file "getfit_IFN_ " is called using
the Matlab function "fminsearch" function to minimize the residual sum
of squares between the models and the data with parameters
log-transformed. The variance in the error is calculated to obtain 95%
confidence intervals around the parameter estimates that minimize the
error.

Command window: the parameter estimates and least squares error
should display for each model, and the code should pause after each
model fit. So the user should press any key when prompted (making sure
the cursor is in the command window), and the code will start fitting
the next model. Once these are all run, the AIC values are returned to
the main file and deltaAIC values are calculated between the 9 model
combinations. 

What figures display:

Figure 1: IFN fold change models fit to data vs. bulk MOI at 8 and 18 hpi
Figure 2: estimated IFN induction rates vs. hpi
left: model 1 (from 0-8 hpi) + models 1-3 (from 8-18 hpi) 
middle: model 2 (from 0-8 hpi) + models 1-3 (from 8-18 hpi)
right: model 3 (from 0-8 hpi) + models 1-3 (from 8-18 hpi) 

