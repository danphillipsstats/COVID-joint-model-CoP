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
The best way to find a Table/Figure is searching for the "Phrase in code" in the File. The "Phrase in code" gives a phrase written in the name of the file which will find the code used to make that file. We aim to keep the line in code correct but mistakes may be present - using "Phrase in code" is more likely to be reliable. For example, searching for the string "Antibody_vs_time_obs_pred_mean" in the file "Antibody_plots.R" will find the code for Figure 1. \
Tables 1-3 and Supplementary Table 8 are not directly produced in the code, so omitted. Supplementary Table 1 and Supplementary Figure 1 is produced in the data cleaning code, which is not provided.
### Main text
| Figure | File  | Phrase in code | Line in code |
| ------------- | ------------- | ---- | ---- |
| Fig. 1  | Antibody_plots  | Antibody_vs_time_obs_pred_mean | 371 |
| Fig. 2  | VE_plots | Mean_VE_vs_antibody | 246 |
| Fig. 3  | VE_plots | Mean_VE_vs_time | 364 |
| Fig. 4  | Antibody_plots | Covariate_effects_PB28_halflife | 898 |
| Fig. 5  | VE_plots | Covariate_effects_mean_VE_vs_time_ | 788 |


### Supplementary Information
| Table | File  | Phrase in code | Line in code |
| ------------- | ------------- | ---- | ---- |
| Supp. Table 1  | Data cleaning | Not provided |  |
| Supp. Table 2  | Antibody_plots  | Antibody_observations | 315 |
| Supp. Table 3  | Antibody_plots | Antibody_outliers | 323 |
| Supp. Table 4  | Antibody_plots | Longitudinal_model_parameters | 930 |
| Supp. Table 5  | Antibody_plots | Covariate_effects_PB28_half_life | 919 |
| Supp. Table 6  | VE_plots | Covariate_effects_hazard.csv | 237 |
| Supp. Table 7  | VE_plots | Antibody_at_VE.csv | 312 |
| Supp. Table 8  | Not in code |  |  |

| Figure | File  | Phrase in code | Line in code |
| ------------- | ------------- | ---- | ---- |
| Supp. Fig. 1  | Data cleaning | Not provided |  |
| Supp. Fig. 2  | Antibody_plots  | Observed_antibody_vs_time_incl_control&outliers | 184 |
| Supp. Fig. 3  | VE_plots | Baseline_hazard | 460 |
| Supp. Fig. 4  | VE_plots | Numbers_at_risk_over_time | 503 |
| Supp. Fig. 5  | VE_plots | Antibody_time_hist | 476 |
| Supp. Fig. 6  | VE_plots | Cases_hist | 523 |
| Supp. Fig. 7  | VE_plots | Covariate_effects_hazard.png | 210 |
| Supp. Fig. 8  | VE_plots | VE_at_antibody_quantiles_vs_time | 595 |
| Supp. Fig. 9  | VE_plots | Antibody_quantiles_vs_time | 629 |
| Supp. Fig. 10  | VE_plots | Covariate_effects_mean_VE_vs_time_ | 788 |
| Supp. Fig. 11 | VE_plots | Omicron_mean_VE_vs_time | 857 |
| Supp. Fig. 12 | VE_plots | Individuals_random_0123_ | 490 |
