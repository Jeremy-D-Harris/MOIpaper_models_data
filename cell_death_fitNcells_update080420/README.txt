updated 04/09/20 -- cell death and virus distribution model fits to data

This folder contains the code that fits 8 model of cell death in total:
4 models under the Poisson distribution assumption and 4 models under
the negative binomial distribution assumption. The data for MDCK cells
and A549 cells are contained in the .mat files "celldeath_MDCK.mat" and
"celldeath_A549.mat." The actual MOI values are the same as those
contained in the source data file tab "Fig1D,E." The cell death data of
MDCK cells is from tab "Figure 2A-D (MDCK)" and in A549 cells, is from
tab "Figure 3B-E (A549)." The flow data give the percent of alive that
are infected at 18 hpi: MDCK cell data is from tab "Figure 2F (MDCK)"
and A549 cell data is in "Figure 3F (A549)"

To run the code, open the m-file "main_estimate_virusproduction" (main)
in Matlab and press 'Run'. This script loads the data and passes it to
the two submain m-files, "submain_celldeath_Poisson" and
"submain_celldeath_negbinomial," each of which fit the cell death models
but under Poisson and negative binomial distribution assumptions.

Within each of the submain files, the file "getfit_deathrate_dispersion"
is called using the Matlab function "fminsearch" function to minimize
the residual sum of squares between the models and the data with
parameters log-transformed. The variance in the error is calculated to
obtain 95% confidence intervals around the parameter estimates that
minimize the error.

Command window: the cell line data to which the models are fitting
should first appear. The parameter estimates and least squares error
should displayed for each model, and the code should pause after each
model fit. So the user should press any key when prompted (making sure
the cursor is in command window), and the code will start fitting the
next model. Lastly, the deltaAIC values between the 8 models are
displayed.

What figures display: 

Figure 1: percent cells remaining model (Poisson) fits to data vs. 
bulk MOI at 3, 6, 12, and 18 hpi 

Figure 2: estimated cell death rates (Poisson) vs. hpi

Figure 3: percent cells remaining model (negative binomial) fits to 
data vs. bulk MOI at 3, 6, 12, and 18 hpi 

Figure 4: estimated cell death rates (negative binomial) vs. hpi

Figure 5: percent of alive cells that are infected at 18 hpi vs. bulk MOI
top: assuming virions are Poisson-distributed across cells
bottom: assuming virions are negative binomial distributed across cells

Some of the models take more time than others to run, for example, the 
parameter estimates of epsilon are close to zero.

