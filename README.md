# Joint modelling of COVID-19 antibody decay and risk of infection
This code accompanies the paper "Improved estimates of COVID-19 correlates of protection, antibody decay and vaccine efficacy waning: a joint modelling approach". \
We run a two-stage joint model for post-vaccination COVID-19 antibody decay and the risk of subsequent COVID-19 infection. We fit the model to data from the COV002 trial of the ChAdOx1 nCoV-19 vaccine. The two-stage model approximates a Bayesian joint model and accounts for uncertainty in the latent antibody trajectories by multiple imputation.
## Project files
### Main files
 - **Longitudinal_hierarchical_model_normalised_ri_rs_exp_slope_t.stan** \
   Stan code for the longitudinal model for log antibody levels. The model assumes log antibody levels decay linearly over time after the peak, and the slope is constrained to be negative (antibodies cannot increase over time). Individual random intercepts and slopes are included. A t-distributed random error is assumed, with Gaussian individual random effects. The data is normalised prior to fitting the model to improve mixing, and then the parameters are transformed back to those for the unnormalised data to improve interpretability. a_0 denotes the random intercept and a_1 denotes the random slope of the log antibody levels. Parameters are defined in the Supplementary Information for the paper.
 - **Longitudinal_antibody_model.R** \
   This samples from the posterior for the antibody trajectories by sampling from the longitudinal model for log antibody levels (Longitudinal_hierarchical_model_normalised_ri_rs_exp_slope_t.stan).
 - **Cox_infection_model.R** \
   For each sample from Longitudinal_antiboy_model.R, this runs a Cox model for risk of COVID-19 infection, where the risk depends on the sample of antibody trajectories. Runs in parallel using the "snowfall"/"parallel" packages in R. We run on a SLURM server. Outputs an array of the MLEs and covariance matrices for the Cox model parameters for each sample of the antibody trajectories, as well as the model design matrix.
 - **Posterior_estimates.R** \
   This pools the samples from the longitudinal antibody model and the Cox infection model to calculate estimates of posterior quantities of interest. Specific outputs are described in the section on this file below. Also runs in parallel using the "snowfall"/"parallel" packages in R, which we ran on a SLURM server.
 - **Antibody_plots.R** \
This file produces plots, numbers and tables for the paper, using the raw data and the output from Longitudinal_antibody_model.R, including the random intercept and slope. The plots, numbers and tables mainly relate to antibodies. The random intercept and slope array is large, hence the plots files were split in two (to avoid issues with memory). Locations of code for specific plots and tables are detailed in the section on this file below.
 - **VE_plots.R** \
This file produces plots, numbers and tables for the paper, using the raw data and the output from Posterior_estimates.R and Cox_infection_model.R. The plots, numbers and tables mainly relate to vaccine efficacy, though some relate to antibody levels. Locations of code for specific plots and tables are detailed in the section on this file below.
### Minor files
 - **functions.R** \
   Contains some user-written functions used in some of the code. Some functions are old and not used any more, only worth looking at if there is a specific function you want to understand.
 - **read_data_set_factors.R** \
   Reads the joint_correlates dataset which contains a row for each individual, and the long_correlates dataset which contains a row for each antibody observation. Sets factors to be correct and calculates a few new simple variables from the data. Note the main data cleaning file is not included in this repository.
### Files not present
 - Read_directories_antibody_plots.R, Read_directories_VE_plots.R \
   These files are referred to in Antibody_plots.R and VE_plots.R, but not provided in this repository. They contain code defining the directories where files were read from and saved to, to run the code in Antibody_plots.R and VE_plots.R. The user should set their own directories.
 - The main data cleaning file is not provided in this repository, nor is the data itself.

## Locations of Tables and Figures

