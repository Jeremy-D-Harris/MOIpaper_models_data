updated 04/09/20 -- virus production model fits to data

This folder contains the code that fits 6 models of virus production to
data. The data for MDCK cells and A549 cells are contained in the .mat
files "virusproduction_MDCK.mat" and "virusproduction_A549.mat." The
actual MOI values are the same as those contained in the source data
file tab "Fig1D,E." For MDCK cells, the viral output data are from tab
"Figure 4A,C-E (MDCK)," and for A549 cells, the data are in tab "Figure
4B,G-I (A549)"

To run the code, open the m-file "main_estimate_virusproduction" (main)
in Matlab and press 'Run'. This script loads the data and passes it to
the two submain m-files "submain_virusproduction_ ," which separately
fit the time-independent and time-dependent models of virus production.

Within each of the submain files, the files "getfit_ " are called using
the Matlab function "fminsearch" function to minimize the residual sum
of squares between the models and the data with parameters
log-transformed. The variance in the error is calculated to obtain 95%
confidence intervals around minimizing parameter estimates.

Command window: the cell line data to which the models are fitting
should first appear. The parameter estimates and least squares error is
displayed for each model. The code should pause after each model fit, so
the user should press any key when prompted (making sure the cursor is
in command window), and the code will start fitting the next model.
Lastly, the deltaAIC values between the 6 models are displayed.

What figures display:

Figure 1: GE model fits to data vs. bulk MOI at 6, 12, and 18 hpi; top
panels: show time-independent model fits bottom panels: show
time-dependent model fits

Figure 2: estimated virus production rates vs. hpi left panel:
time-independent models parametrized with best-fit values right panel:
time-dependent models parametrized with best-fit values

Some of the models take more time than others to run, for example,
fitting the time-dependent, saturating model to MDCK cell data

