##############
# This file uses the output from the Cox model in Cox_infection_model.R
# to calculate the plots, figures and results needed for the paper.

# Create a list to store the time taken to run different parts of the file
runtimes <- list()
t_init <- Sys.time()
print(t_init)
runtimes$start <- t_init

# Set what the outcome for this analysis will be 
# (pos = all COVID-19 infections, prim = primary symptomatic COVID-19 infections)
event_outcome <- c("pos","prim")[2] # To be run with both infection outcomes separately

# Set which type of antibodies we are using
antibody_type <- c("S","neuts")[1]
long_type <- c("7_day_incubator","0_day_incubator")[1]
# The name of the longitudinal model
long_file_name <- "ri_rs_7inc_pos_prim_exp_slope_t_all_covariates_t0PB28"
# The name for this file
model_name <- "_submission"
if(!(substring(long_file_name,1,5)=="neuts") == (antibody_type=="neuts")){warning("Your long_file_name and antibody_type do not match. You may be running a neuts analysis on S data or vice versa. Please double check.")}

# Define the directories in which files are saved - these may need adjusting for other users.
directory_correlates <- getwd()
directory_correlates_project <- paste0(directory_correlates,"/3_Programs")
output_directory <- paste0(directory_correlates,"/4_Output")
main_plot_directory <- paste0(output_directory,"/Plots")
# Create a plot folder for the specific model run in Longitudinal_antibody_model.R
if (!dir.exists(main_plot_directory)){dir.create(main_plot_directory)}
plot_directory <- paste0(main_plot_directory,"/",long_file_name)
if (!dir.exists(plot_directory)){
  dir.create(plot_directory)
  dir.create(paste0(plot_directory,"/Antibody_plots"))
}

#load packages
library(rstan)
library(survival)
library(parallel)

#######
# Load user-written functions
source(paste0(directory_correlates_project,"/functions.R"))
################################################################################
# Two-stage model, second stage
cox_model_name <- paste0(event_outcome,"_site_parallel_simple_7inc")
cox_model_namelong <- paste0("_",cox_model_name)
cox_stage_two <- readRDS(paste0(output_directory,"/Cox_infection_model_long_",long_file_name,cox_model_namelong,".RDS"))

# Longitudinal model
# For an output which is a list, first element being (possibly thinned) extracted samples
# the second element being the data
# Take a subsample of size
nsamples <- dim(cox_stage_two$cox_model_pred)[1]
nchains <- dim(cox_stage_two$cox_model_pred)[2]
long_out_rs <- readRDS(paste0(output_directory,"/Correlates_long_",long_file_name,"_reffects.RDS"))
# File is a list, first entry is the Stan output, second entry the inputted data.
data_long <- long_out_rs[[2]]
long_out_rs <- long_out_rs[[1]]
totalsample <- dim(long_out_rs)[1]
nsamples <- min(nsamples,totalsample)
if (nsamples<totalsample){
  sample_choices <- seq(from=1, to =totalsample,by=totalsample/nsamples)
  long_out_rs <- long_out_rs[sample_choices,,]
}
# Find the population parameters for intercept and slope
a_0_array <- long_out_rs[,,which(substring(dimnames(long_out_rs)$parameters,1,4)=="a_0[")]
a_1_array <- long_out_rs[,,which(substring(dimnames(long_out_rs)$parameters,1,4)=="a_1[")]
n <- dim(a_0_array)[3]
rm(long_out_rs)
runtimes$a_array <- Sys.time();

# Generate correlated gamma and zeta
set.seed(1234)
Z0_norm <- matrix(rnorm(n=nsamples*nchains),nrow=nsamples,ncol=nchains) # Standard normal random variables
Z1_norm <- matrix(rnorm(n=nsamples*nchains),nrow=nsamples,ncol=nchains)
# Calculate the correlation
cor_gamma_zeta <- cox_stage_two$cox_model_var[,,"antibody.As_vaccinated_arm_2ChAdOx1"]/sqrt(cox_stage_two$cox_model_var[,,"antibody.antibody"]*cox_stage_two$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"])
# Transform the standard normal random variables to normal random variables with correct mean, and covariances
gammas <- cox_stage_two$cox_model_pred[,,"antibody"] + Z0_norm*sqrt(cox_stage_two$cox_model_var[,,"antibody.antibody"])
zetas <- cox_stage_two$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"] + (Z0_norm*cor_gamma_zeta + Z1_norm*sqrt(1-cor_gamma_zeta^2))*sqrt(cox_stage_two$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"])

###############
# Read joint_correlates and long_correlates
# Set their factors to be correct
source(paste0(directory_correlates_project,"/read_data_set_factors.R"))
if(antibody_type=="S"){ # For anti-spike IgG
  # Create long_correlates to be the long_correlates dataset with the appropriate antibodies and log antibodies
  long_correlates <- long_correlates_S
  long_correlates$antibody <- long_correlates$S.PPD
  long_correlates$log_antibody <- log(long_correlates$antibody)
  # Limit of detection
  LOD <- 33
  # Conversion factor from arbitrary units AU/mL in our data to international unit BAU/mL
  conv_factor <- 0.00645
  print(LOD*conv_factor)
  antibody_name <- "S IgG"
  antibody_units <- "BAU/mL"
  # Marking outliers
  long_correlates$outlier <- F
  # Exclude the exceptionally small results as they appear to be outliers and this appears to be preventing good mixing in the model (rho = -1)
  long_correlates$outlier[which(long_correlates$log_antibody<7 & long_correlates$As_vaccinated_arm_2=="ChAdOx1")] <- T
  # Exclude the exceptionally large result at PB90
  long_correlates$outlier[which((long_correlates$antibody_time>50 &
                                   long_correlates$antibody_time<130 &
                                   long_correlates$log_antibody>13) &
                                  long_correlates$As_vaccinated_arm_2=="ChAdOx1")] <- T
} else if(antibody_type=="neuts"){ # If using neutralising antibodies
  long_correlates <- long_correlates_neuts
  long_correlates$antibody <- long_correlates$neuts
  long_correlates$log_antibody <- log(long_correlates$antibody)
  # Limit of detection
  LOD <- 40
  # Conversion factor from ND50 in our data to international unit IU/mL
  conv_factor <- 0.1458
  print(LOD*conv_factor)
  antibody_name <- "nAb"
  antibody_units <- "IU/mL"
  # Marking outliers
  long_correlates$outlier <- F # No outliers
}
## Left-censoring below the limit of detection
# Check how many antibody observations are less than, or equal to, the limit of detection, and treat them as left-censored
nantibody_LOD <- sum(long_correlates$log_antibody==log(LOD))
nantibody_leq_LOD <- sum(long_correlates$log_antibody<=log(LOD))
if(nantibody_LOD<nantibody_leq_LOD & nantibody_LOD>0){suppressWarnings(paste0("There are ",nantibody_LOD," antibody observations equal to the limit of detection. We will treat these as though they are not censored, but observed exactly."))}
if(nantibody_LOD==nantibody_leq_LOD & nantibody_LOD>0){warning(paste0("There are ",nantibody_LOD," antibody observations equal to the limit of detection. We will treat these as though they are not censored, but observed exactly."))}
long_correlates$left_censored <- FALSE
long_correlates$left_censored[which(long_correlates$antibody<LOD)] <- TRUE
long_correlates$antibody[which(long_correlates$left_censored)] <- LOD
long_correlates$log_antibody[which(long_correlates$left_censored)] <- log(LOD)
long_correlates <- long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),] 

# Defining variables
long_correlates <- long_correlates[,c("redcap_data_access_group","As_vaccinated_arm_2","sc_repeat_pid","record_id","sc_age","age_group","cor2dose_positive","cor2dose_primary","cor2dose_outcome_pos_prim","antibody","log_antibody","antibody_time","ll","timepoint","visit2","left_censored","outlier")]

# times
# Final antibody time we will report
amax_time <- 209-7*(long_type=="7_day_incubator") # Minus 7 days if we're using the 7_day_incubator method
# Earliest antibody time we will report
amin_time <- 28
adsb <- amin_time:amax_time # Days since boost (ds_t0+21)
if(long_type=="7_day_incubator") {
  t0 <- 21 + data_long$t_center
  # VE time is 7 days after antibody time due to assuming a 7 day gap between exposure and reported infection
  VEdsb <- adsb + 7 # Days since boost (ds_t0+21)
} else if(long_type=="0_day_incubator") {t0 <- 14 + data_long$t_center; VEdsb <- adsb} # Days since boost (ds_t0+21)
# Number of days from boost to day 0 in the model (i.e. 28)
ds_t0 <- adsb-t0 # Days since day 0 in the model for the antibodies
nts <- length(ds_t0)
time_units <- c("Days","Weeks","Months")[1]
# div is the Amount to divide days by to get the time units we're using (1 if days, 7 if weeks, 365/12 if months)
if (time_units=="Days"){div <- 1}
if (time_units=="Weeks"){div <- 7}
if (time_units=="Months"){div <- 365/12}
atsb <- adsb/div # time since boost in the appropriate units
VEtsb <- VEdsb/div # time since boost in the appropriate units

# The data for vaccinated individuals only
joint_correlates_ChAd <- joint_correlates[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1"),]
# Set the column names for a_0_array and a_1_array
dimnames(a_0_array)$parameters <- joint_correlates_ChAd$sc_repeat_pid
dimnames(a_1_array)$parameters <- joint_correlates_ChAd$sc_repeat_pid

# Create histogram of observed antibody values and predicted antibody values at end of study day 209
antibody_breaks <- seq(-5,6,by=0.2)
obs_antibodies <- long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),]$antibody*conv_factor
pred_antibody_values_28 <- exp(a_0_array+(28-t0)*a_1_array)*conv_factor
pred_antibody_values_28_hist <- hist(log10(pred_antibody_values_28), plot=F, breaks = antibody_breaks)
# Calculate proportion of predicted antibodies less than 10, 100, 1000 BAU/mL
round(mean(pred_antibody_values_28<10)*100,1)
round(mean(pred_antibody_values_28<100)*100,1)
round(mean(pred_antibody_values_28<1000)*100,1)
rm(pred_antibody_values_28)
pred_antibody_values_182 <- exp(a_0_array+(182-t0)*a_1_array)*conv_factor
pred_antibody_values_182_hist <- hist(log10(pred_antibody_values_182), plot=F,breaks=antibody_breaks)
# Calculate proportion of predicted antibodies less than 10, 100, 1000 BAU/mL
round(mean(pred_antibody_values_182<10)*100,1)
round(mean(pred_antibody_values_182<100)*100,1)
round(mean(pred_antibody_values_182<1000)*100,1)
rm(pred_antibody_values_182)
obs_antibody_values_hist <- hist(log10(obs_antibodies), plot=F,breaks=antibody_breaks)
obs_antibody_values_hist$counts <- obs_antibody_values_hist$density/max(obs_antibody_values_hist$density)*40
# Change histogram breaks to be in terms of antibodies not log antibodies
obs_antibody_values_hist$breaks <- 10^(obs_antibody_values_hist$breaks)
pred_antibody_values_28_hist$counts <- pred_antibody_values_28_hist$density/max(obs_antibody_values_hist$density)*40
pred_antibody_values_28_hist$breaks <- 10^(pred_antibody_values_28_hist$breaks)
pred_antibody_values_182_hist$counts <- pred_antibody_values_182_hist$density/max(obs_antibody_values_hist$density)*40
pred_antibody_values_182_hist$breaks <- 10^(pred_antibody_values_182_hist$breaks)
rm(obs_antibodies)

# Read data of curve from Wei et al. for extrapolation to Omicron
# Central estimate of antibody relative efficacy relationship
omicron_data <- read.csv(paste0(directory_correlates,"/2_Data_for_Analysis/Omicron_correlates_Wei_et_al_estimate.csv"),col.names = c("VE","antibody"))
omicron_data <- as.data.frame(omicron_data)
omicron_data <- omicron_data[order(omicron_data$antibody),]
omicron_data$yy <- log(1-omicron_data$VE/100)
loess_estimate <- loess(yy~antibody,data=omicron_data,span=0.3)
plot(VE~antibody,data=omicron_data, ylim=c(0,100))
predictions <- 100*(1-exp(predict(loess_estimate,omicron_data$antibody)))
lines(predictions~omicron_data$antibody,lty=1)
100*(1-exp(predict(loess_estimate,0)))

# Lower 95% estimate of antibody relative efficacy relationship
omicron_data <- read.csv(paste0(directory_correlates,"/2_Data_for_Analysis/Omicron_correlates_Wei_et_al_lower.csv"),col.names = c("VE","antibody"))
omicron_data <- as.data.frame(omicron_data)
omicron_data <- omicron_data[order(omicron_data$antibody),]
omicron_data$yy <- log(1-omicron_data$VE/100)
loess_lower_estimate <- loess(yy~antibody,data=omicron_data,span=0.3)
points(VE~antibody,data=omicron_data, col="blue",lty=2)
predictions <- 100*(1-exp(predict(loess_lower_estimate,omicron_data$antibody)))
lines(predictions~omicron_data$antibody, col="blue",lty=2)
100*(1-exp(predict(loess_lower_estimate,0)))

# Upper 95% estimate of antibody relative efficacy relationship
omicron_data <- read.csv(paste0(directory_correlates,"/2_Data_for_Analysis/Omicron_correlates_Wei_et_al_upper.csv"),col.names = c("VE","antibody"))
omicron_data <- as.data.frame(omicron_data)
omicron_data <- omicron_data[order(omicron_data$antibody),]
omicron_data$yy <- log(1-omicron_data$VE/100)
loess_upper_estimate <- loess(yy~antibody,data=omicron_data,span=0.3)
points(VE~antibody,data=omicron_data, col="blue",lty=2)
predictions <- 100*(1-exp(predict(loess_upper_estimate,omicron_data$antibody)))
lines(predictions~omicron_data$antibody, col="blue",lty=2)
100*(1-exp(predict(loess_upper_estimate,0)))

# Generate the antibody values needed to measure VE vs antibody
quants <- c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)
nquants <- length(quants)
# Convert BAU/ml antibody levels of 1,10,...,10^5 to log AU/ml units
tens_BAU <- 10^seq(-1,7,1)
tens <- log(10^seq(-1,7,1)/conv_factor)
# Exclude those less than the smallest observed value or greater than the largest value
which_tens <- which(tens>min(long_correlates$log_antibody) & tens<max(long_correlates$log_antibody))
tens <- tens[which_tens]
tens_BAU <- tens_BAU[which_tens]

# Set the antibody values at which to test VE vs antibody
# Observe log antibody levels at an even sequence from min antibodies to max antibodies,
# and at the quantiles of the observed antibody levels, and at the powers of 10
log_antibody_values <- c(seq(from=min(long_correlates$log_antibody), to=max(long_correlates$log_antibody), by = 0.01),quantile(long_correlates$log_antibody,c(0,quants,1)),quantile(a_0_array+max(ds_t0)*a_1_array,seq(0.01,1-0.01,0.01)))
log_antibody_values <- log_antibody_values[!duplicated(log_antibody_values)]
antibody_values <- exp(log_antibody_values)
BAU_antibodies <- antibody_values*conv_factor
BAU_antibodies <- c(BAU_antibodies,tens_BAU)
log_antibody_values <- c(log_antibody_values,tens)
log_antibody_values <- log_antibody_values[!duplicated(log_antibody_values)] # Remove duplicated values
log_antibody_values <- log_antibody_values[order(log_antibody_values)]
antibody_values <- exp(log_antibody_values)
BAU_antibodies <- BAU_antibodies[!duplicated(BAU_antibodies)]
BAU_antibodies <- BAU_antibodies[order(BAU_antibodies)]
nant_val <- length(log_antibody_values)
# Record the time
runtimes$data_setup <- Sys.time();
nwork <- length(parallelly::availableWorkers())

# Calculate VE against antibody
VE_summary_chains <- array(dim=c(nchains,nant_val,4))
for (chain in seq_len(nchains)){ # Calculate VE from that chain
  VE_samples <- 100*(1 - exp(outer(gammas[,chain],antibody_values) + zetas[,chain])) # Given as a percentage
  VE_summary_chains[chain,,] <- cbind(colMeans(VE_samples),t(apply(VE_samples,2,quantile,c(0.5,0.025,0.975))))
  rm(VE_samples)
}
dimnames(VE_summary_chains)[[2]] <- BAU_antibodies
dimnames(VE_summary_chains)[[3]] <- c("mean","median","lower95CI","upper95CI")
# Plot the summary for the four chains on the same graph to check they look the same.

# Overall summary
VE_samples <- 100*(1 - exp(outer(c(gammas),antibody_values) + c(zetas))) # Given as a percentage
VE_summary <- cbind(colMeans(VE_samples),t(apply(VE_samples,2,quantile,c(0.5,0.025,0.975))))
rm(VE_samples)
rownames(VE_summary) <- BAU_antibodies
colnames(VE_summary) <- c("mean","median","lower95CI","upper95CI")

runtimes$VE_summary <- Sys.time();

# Calculate the antibody associated with x% VE
VE <- c(seq(0,1-0.05,0.05),0.975,0.99,0.67)
VE <- VE[order(VE)]
quants <- c(0.5,0.025,0.975)
# For the different chains
antibody_at_VE_summary_chains <- array(dim=c(nchains,length(VE),length(quants)))
# Name the dimensions of the 
dimnames(antibody_at_VE_summary_chains)[[1]] <- paste0("chain::",1:4)
dimnames(antibody_at_VE_summary_chains)[[2]] <- VE
dimnames(antibody_at_VE_summary_chains)[[3]] <- quants
for (chain in seq_len(nchains)){ # Calculate the antibody level required for a certain VE
  antibody_at_VE_samples <- sweep(-outer(zetas[,chain],log(1-VE),FUN="-"),1,gammas[,chain],"/")
  colnames(antibody_at_VE_samples) <- VE
  antibody_at_VE_summary_chains[chain,,] <- t(apply(antibody_at_VE_samples,2,quantile,quants))
  rm(antibody_at_VE_samples)
}
antibody_at_VE_summary_chains[which(antibody_at_VE_summary_chains<0)] <- NaN
antibody_at_VE_summary_chains <- antibody_at_VE_summary_chains*conv_factor
signif(antibody_at_VE_summary_chains[,as.character(c(0.5,0.7,0.9)),],3)

# Overall
antibody_at_VE_samples <- sweep(-outer(c(zetas),log(1-VE),FUN="-"),1,c(gammas),"/")
colnames(antibody_at_VE_samples) <- VE
antibody_at_VE_summary <- t(apply(antibody_at_VE_samples,2,quantile,quants))
rm(antibody_at_VE_samples)
antibody_at_VE_summary[which(antibody_at_VE_summary<0)] <- NaN
colnames(antibody_at_VE_summary) <- c("median","lower95CI","upper95CI")
# Note negative antibody is sometimes produced but antibody cannot be negative.
antibody_at_VE_summary <- antibody_at_VE_summary*conv_factor
signif(antibody_at_VE_summary[as.character(c(0.5,0.7,0.9)),],2)

runtimes$antibody_at_VE_summary <- Sys.time();

# Value at 0 antibodies/direct effect
direct_effects <- 100*(1-exp(zetas))
direct_effect_CI_chains <- cbind(colMeans(direct_effects),t(apply(direct_effects,2,quantile,c(0.5,0.025,0.975))))
colnames(direct_effect_CI_chains) <- c("mean","median","lower95CI","upper95CI")
direct_effect_CI <- c(mean(direct_effects),quantile(direct_effects,c(0.5,0.025,0.975)))
names(direct_effect_CI) <- c("mean","median","lower95CI","upper95CI")
rm(direct_effects)

# Antibody indirect effect
VE_indirect_samples <- 100*(1 - exp(outer(gammas,antibody_values)))
library(abind)
VE_indirect_summary_chains <- abind(apply(VE_indirect_samples,c(2,3),mean),aperm(apply(VE_indirect_samples,c(2,3),quantile,c(0.5,0.025,0.975)),c(2,3,1)), along=3)
rm(VE_indirect_samples)
dimnames(VE_indirect_summary_chains)[[2]] <- BAU_antibodies
dimnames(VE_indirect_summary_chains)[[3]] <- c("mean","median","lower95CI","upper95CI")

# Overall
VE_indirect_samples <- 100*(1 - exp(outer(c(gammas),antibody_values)))
VE_indirect_summary <- cbind(apply(VE_indirect_samples,2,mean),t(apply(VE_indirect_samples,2,quantile,c(0.5,0.025,0.975))))
rm(VE_indirect_samples)
rownames(VE_indirect_summary) <- BAU_antibodies
colnames(VE_indirect_summary) <- c("mean","median","lower95CI","upper95CI")

runtimes$direct_indirect_VE <- Sys.time();


##############################################
# Calculate parameter effects
# Find the variances (diagonal entries of covariance matrix)
npred <- dim(cox_stage_two$cox_model_pred)[3]
# Draw independently sampled effects, in an array of the same size.
# Note they don't need to be correlated as this is for forest plot
draw_effect_samples_chain <- function(chain,effect_name,effect_means,effect_covariances){
  effect_var_name <- paste(effect_name,effect_name,sep=".")
  effect_mean <- effect_means[,chain,effect_name]
  effect_sd <- sqrt(effect_covariances[,chain,effect_var_name])
  return(rnorm(n=length(effect_mean),mean=effect_mean,sd=effect_sd))
}
draw_effect_samples <- function(effect_name,effect_means,effect_covariances){
  return(sapply(seq_len(nchains),draw_effect_samples_chain,effect_name,effect_means,effect_covariances))
}
# vapply allows you to determine the shape of the object as an array
set.seed(1234)
effect_samples <- vapply(dimnames(cox_stage_two$cox_model_pred)[[3]],FUN=draw_effect_samples,effect_means = cox_stage_two$cox_model_pred, effect_covariances = cox_stage_two$cox_model_var, FUN.VALUE = array(0,dim=c(nsamples,nchains)) )
# Check the number of samples is sufficient for the hazard parameters
apply(effect_samples,3,Rhat)
apply(effect_samples,3,ess_tail)
apply(effect_samples,3,ess_bulk)
# Check the answers give the same. 
exp(apply(effect_samples,c(2,3),mean))
exp(aperm(apply(effect_samples,c(2,3),quantile,c(0.5,0.025,0.975)),c(2,1,3)))
# Slight differences at two significant figures for the 2.5 and 97.5% quantiles
effect_summary <- cbind(apply(effect_samples,c(3),mean),t(apply(effect_samples,c(3),quantile,c(0.5,0.025,0.975))))
colnames(effect_summary) <- c("mean","median","lower95CI","upper95CI")
# Adjust the antibody effect
if(antibody_name=="S IgG"){
  antibody_multiplier_WHO <- 100
} else if(antibody_name=="nAb"){antibody_multiplier_WHO <- 10}

antibody_multiplier <- antibody_multiplier_WHO/conv_factor
effect_summary["antibody",] <- effect_summary["antibody",]*antibody_multiplier
effect_summary <- cbind(effect_summary,exp(effect_summary))
colnames(effect_summary) <- paste0(rep(c("log",""),each=4),rep(c("mean","median","lower95CI","upper95CI"),2))
effect_summary <- round(as.data.frame(effect_summary),5)
# Give a list of all possible rownames from all models which have been considered at some point
{rownames_full <- c(paste0("Antibody level\n (Effect due to increase of ",antibody_multiplier_WHO," ",antibody_units,")"),
                    paste0("Square root of antibody level\n (Effect due to increase of sqrt(",(antibody_multiplier_WHO),") sqrt(BAU/mL))"),
                    "Vaccination direct effect\n (ChAdOx1 nCoV-19 vs control)",
                    "Dose schedule\n (Two low doses vs two standard doses)",
                    "Dose schedule\n (Low dose followed by low dose vs two standard doses)",
                    "Age\n (56-69 vs 18-55 years)", "Age\n (\u226570 vs 18-55 years)", 
                    "Sex\n (Female vs Male)", "Ethnicity\n (Non-white vs White)",
                    "Comorbidity\n (Comorbidity vs None)", "BMI\n (\u226530 vs <30)",
                    "Interval between first and second dose\n (9-11 weeks vs \u226512 weeks)",
                    "Interval between first and second dose\n (6-8 weeks vs \u226512 weeks)",
                    "Interval between first and second dose\n (<6 weeks vs \u226512 weeks)", 
                    "Healthcare worker (HCW) status\n (HCW facing <1 COVID-19 patient per day vs Not a HCW)",
                    "Healthcare worker (HCW) status\n (HCW worker facing \u22651 COVID-19 patient per day vs Not a HCW)",
                    "First dose\n (Low dose vs Standard dose)",
                    paste0("Antibody LDLD interaction\n (Additional effect due to increase of ",antibody_multiplier_WHO," BAU/mL for those LDLD)"),
                    paste0("Antibody LDSD interaction\n (Additional effect due to increase of ",antibody_multiplier_WHO," BAU/mL for those LDSD)"),
                    paste0("Antibody interval interaction\n (Additional effect due to increase of ",antibody_multiplier_WHO," BAU/mL for those 9-11 weeks)"),
                    paste0("Antibody interval interaction\n (Additional effect due to increase of ",antibody_multiplier_WHO," BAU/mL for those 6-8 weeks)"),
                    paste0("Antibody interval interaction\n (Additional effect due to increase of ",antibody_multiplier_WHO," BAU/mL for those <6 weeks)"),
                    "Vaccination LDLD interaction\n (Additional effect due to ChAdOx1 vaccination for those LDLD)",
                    "Vaccination LDSD interaction\n (Additional effect due to ChAdOx1 vaccination for those LDSD)",
                    "Vaccination interval interaction\n (Additional effect due to ChAdOx1 vaccination for those 9-11 weeks)",
                    "Vaccination interval interaction\n (Additional effect due to ChAdOx1 vaccination for those 6-8 weeks)",
                    "Vaccination interval interaction\n (Additional effect due to ChAdOx1 vaccination for those <6 weeks)")
  rowlabels_full <- c("antibody","sqrt(antibody)","As_vaccinated_arm_2ChAdOx1",
                      "cor2dose_scheduleLDLD","cor2dose_scheduleLDSD",
                      "age_group56-69", "age_group>=70", 
                      "sc_genderFemale", "cor2dose_non_whiteOther",
                      "cor2dose_comorbiditiesComorbidity", "cor2dose_bmi_geq_30BMI_geq_30",
                      "cor2dose_interval9-11",
                      "cor2dose_interval6-8",
                      "cor2dose_interval<6", 
                      "cor2dose_hcw_statusLess_than_1_covid",
                      "cor2dose_hcw_statusMore_than_1_covid",
                      "cor2dose_primescheduleLD",
                      "antibody:cor2dose_scheduleLDLD",
                      "antibody:cor2dose_scheduleLDSD",           
                      "antibody:cor2dose_interval9-11",
                      "antibody:cor2dose_interval6-8",                
                      "antibody:cor2dose_interval<6",
                      "As_vaccinated_arm_2ChAdOx1:cor2dose_scheduleLDLD",
                      "As_vaccinated_arm_2ChAdOx1:cor2dose_scheduleLDSD",
                      "As_vaccinated_arm_2ChAdOx1:cor2dose_interval9-11",
                      "As_vaccinated_arm_2ChAdOx1:cor2dose_interval6-8",
                      "As_vaccinated_arm_2ChAdOx1:cor2dose_interval<6")
  names(rownames_full) <- rowlabels_full} # Match the rowname used to the rowlabel
effect_summary$rownames <- rownames_full[match(rownames(effect_summary),names(rownames_full))]
effect_summary$rownames <- factor(effect_summary$rownames,levels=effect_summary$rownames[nrow(effect_summary):1])

runtimes$effect_summary <- Sys.time();
rm(cox_stage_two)

##############################################
# Set up the analysis of covariate effects on VE over time
# That is, predicting VE for new individuals
long_out <- readRDS(paste0(output_directory,"/Correlates_long_",long_file_name,"_noreffects.RDS"))
# File is a list, first entry is the Stan output, second entry the inputted data.
data_long <- long_out[[2]]
long_out <- long_out[[1]]
long_out <- long_out[,,grepl( "_tf", dimnames(long_out)$parameters, fixed = TRUE)] # Keep only the transformed parameters _tf

long_covariates <- c("age_group","sc_gender","cor2dose_non_white","cor2dose_comorbidities","cor2dose_bmi_geq_30","cor2dose_interval","cor2dose_hcw_status","cor2dose_primeschedule","cor2dose_outcome_pos_prim","NelsonAalen")
n_new_ind <- 13
new_data <- joint_correlates[1:n_new_ind,long_covariates]

for (i in 1:n_new_ind){ # Set them all to be baseline in the covariates
  new_data[i,which(!is.element(long_covariates,c("NelsonAalen","cor2dose_outcome_pos_prim")))] <- c("18-55","Female","White","None","BMI_less_than_30",">=12","Not_hcw","SD")
}
# Temporarily make this row different to avoid problems with model.matrix
new_data[n_new_ind,which(long_covariates!="NelsonAalen")] <- c("18-55","Female","White","None","BMI_less_than_30",">=12","Not_hcw","SD","Positive")
# Create rows which differ from the baseline covariates
new_data_names <- c("baseline","56-69",">=70","Male","Other","Comorbidity","BMI_geq_30","9-11","6-8","<6","Less_than_1_covid","More_than_1_covid","LD")
new_data[2,1] <- new_data_names[2]
new_data[3,1] <- new_data_names[3]
new_data[4,2] <- new_data_names[4]
new_data[5,3] <- new_data_names[5]
new_data[6,4] <- new_data_names[6]
new_data[7,5] <- new_data_names[7]
new_data[8,6] <- new_data_names[8]
new_data[9,6] <- new_data_names[9]
new_data[10,6] <- new_data_names[10]
new_data[11,7] <- new_data_names[11]
new_data[12,7] <- new_data_names[12]
new_data[13,8] <- new_data_names[13]
rownames(new_data) <- new_data_names
new_rows <- as.data.frame(model.matrix(~.,data=new_data))[,-1]
new_rows_tf <- sweep(sweep(new_rows,2, colMeans(data_long$X_l),FUN="-"),2,apply(data_long$X_l,2,sd),FUN="/")
# Create dataset with indicator columns for positive instead of a factor
new_data2 <- new_data[,!is.element(colnames(new_data),"cor2dose_outcome_pos_prim")]
new_data2[,"cor2dose_positive"] <- as.character(new_data$cor2dose_outcome_pos_prim)
new_data2[which(new_data2[,"cor2dose_positive"]=="Primary"),"cor2dose_positive"] <- "Positive"
new_data2[,"cor2dose_positive"] <- factor(new_data2[,"cor2dose_positive"],levels=c("Negative","Positive"))
# Create dataset with indicator columns for primary instead of a factor
new_data3 <- new_data[,!is.element(colnames(new_data),"cor2dose_outcome_pos_prim")]
new_data3[,"cor2dose_primary"] <- as.character(new_data$cor2dose_outcome_pos_prim)
new_data3[which(new_data3[,"cor2dose_primary"]=="Positive"),"cor2dose_primary"] <- "Negative"
new_data3[,"cor2dose_primary"] <- factor(new_data3[,"cor2dose_primary"],levels=c("Negative","Positive"))
levels(new_data3[,"cor2dose_primary"]) <- c("Non-case","Case")

# Impute the missing outcome and Nelson-Aalen estimates
long_covariates2 <- c("age_group","sc_gender","cor2dose_non_white","cor2dose_comorbidities","cor2dose_bmi_geq_30","cor2dose_interval","cor2dose_hcw_status","cor2dose_primeschedule","cor2dose_positive")
long_covariates3 <- c("age_group","sc_gender","cor2dose_non_white","cor2dose_comorbidities","cor2dose_bmi_geq_30","cor2dose_interval","cor2dose_hcw_status","cor2dose_primeschedule","cor2dose_primary")
joint_correlates_ChAd <- joint_correlates[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1"),]
# Logistic regression for probability of testing positive on all other covariates except NelsonAalen
logmI <- glm(cor2dose_positive~., data = joint_correlates_ChAd[,long_covariates2], family = binomial)
chol_coef_logmI <- chol(vcov(logmI))
# Logistic regression for probability of symptoms (given test positive) on all other covariates except NelsonAalen
logmS <- glm(cor2dose_primary~., data = joint_correlates_ChAd[which(joint_correlates_ChAd$cor2dose_positive=="Positive"),long_covariates3], family = binomial)
chol_coef_logmS <- chol(vcov(logmS))
# Linear regression for NelsonAalen estimate given all other covariates (including testing positive and symptoms)
linmNA <- lm(NelsonAalen~., data = joint_correlates_ChAd[,long_covariates])
chol_coef_linmNA <- chol(vcov(linmNA))
plot(density(linmNA$fitted.values))

mmlogm <- model.matrix(cor2dose_positive~.,new_data2[,long_covariates2])
mmlinm <- model.matrix(NelsonAalen~.,new_data[,long_covariates])

# Columns which contain the transformed beta variables
beta_l0s_tf_cols <- which(is.element(substr(dimnames(long_out)$parameters,start=1,stop=8),c("beta_l0_")))
beta_l1s_tf_cols <- which(is.element(substr(dimnames(long_out)$parameters,start=1,stop=8),c("beta_l1_")))

# Number of samples of the random error to draw per iteration when calculating covariate effects
n_error_samples <- n

print(object.size(a_0_array),units="GB")
parameters_list <- lapply(seq_len(nchains), function(chain) lapply(seq_len(nsamples),function(i) {
  ls <- list(a_0_array[i,chain,],a_1_array[i,chain,],gammas[i,chain],zetas[i,chain],long_out[i,chain,])
  names(ls) <- c("a_0_vec","a_1_vec","gamma","zeta","long_out")
  return(ls)
}))
print(object.size(parameters_list),units="GB")
sort( sapply(ls(),function(x){object.size(get(x))})) 
rm(a_0_array,a_1_array,Z0_norm,Z1_norm,cor_gamma_zeta,gammas,zetas,long_out)

quant_diff <- 0.05
effect_vs_time_quants <- seq(0+quant_diff,1-quant_diff,quant_diff)
#
dimnames_meanVE <- c("mean_VE",effect_vs_time_quants)
dimnames_omicron <- paste("Omicron_mean_VE",c("mean_effect","upper95_effect","lower95_effect"),sep="_")
dimnames_new_ind <- new_data_names
dimnames_antibody_quantiles <- paste0("antibody_quantiles_",effect_vs_time_quants)


##############################################
# Analysis of effects over time
# We run a parallelised apply function to run the whole analysis for each parameter draw
# then later pool them together
outputs_dim <- 1+2*length(effect_vs_time_quants)+3+n_new_ind
sfLapply_VE_vs_time_analysis <- function(chain, my_seed = 1234){
  parameters_sublist <- parameters_list[[chain]]
  # Begin a cluster in snowfall package
  cluster <- snowfall::sfInit(nwork,type="SOCK", parallel=T)   
  # Load the relevant objects in the parallel workspaces - including 
  snowfall::sfExport("loess_estimate","loess_lower_estimate","loess_upper_estimate",
                     "antibody_values","ds_t0","atsb","effect_vs_time_quants","conv_factor","nts",
                     "VEtsb",
                     "new_data","new_rows_tf","n_new_ind","new_data_names",
                     "mmlinm","mmlogm",
                     "chol_coef_logmI","chol_coef_logmS","chol_coef_linmNA",
                     "logmI","logmS","linmNA",
                     "beta_l0s_tf_cols","beta_l1s_tf_cols","data_long",
                     "n_error_samples",
                     "dimnames_meanVE","dimnames_omicron","dimnames_new_ind")
  # Load libraries
  snowfall::sfLibrary(survival)
  snowfall::sfLibrary(rstan)
  # Set up the random number generation
  snowfall::sfClusterSetupRNG(seed=my_seed)
  print(paste("Number of workers:",nwork))
  
  # Runs the function on each element of the list in parallel
  # That is, runs a cox model for each imputation of antibody levels.
  out_list <- snowfall::sfLapply(x = parameters_sublist, fun = VE_vs_time_analysis)
  snowfall::sfStop()
  out_chain <- vapply(out_list,function(x) matrix(x,nrow=nts,dimnames=list(VEtsb,c(dimnames_meanVE,dimnames_antibody_quantiles,dimnames_omicron,dimnames_new_ind))),FUN.VALUE=matrix(0,nts,outputs_dim))
  out_chain <- aperm(out_chain,c(3,1,2))
  return(out_chain)
}

VE_vs_time_analysis <- function(parameter_draw_list){
  attach(parameter_draw_list)
  antibody_mat <- exp(a_0_vec+outer(a_1_vec,ds_t0))
  # Antibody quantiles
  antibody_quantiles <- t(apply(antibody_mat,2,quantile,effect_vs_time_quants)*conv_factor)
  
  ################
  # Mean VE vs time
  # Vaccination effect at t (overall and means for each individual)
  ind_VEs <- 100*(1-exp(gamma*antibody_mat+zeta))
  VE_vs_time <- cbind(colMeans(ind_VEs),t(apply(ind_VEs,2,quantile,effect_vs_time_quants)))
  
  ############
  # Extrapolation to Omicron
  antibody_BAU_mat <- antibody_mat*conv_factor
  # Truncate values greater than the largest value of antibodies from the Wei study
  antibody_BAU_mat <- pmin(antibody_BAU_mat,max(loess_estimate$x))
  
  ind_VEs_estimate <- 100*(1-exp(apply(antibody_BAU_mat,2,function(x) predict(loess_estimate,x))))
  loess_estimate_means <- colMeans(ind_VEs)
  ind_VEs_lower_estimate <- 100*(1-exp(apply(antibody_BAU_mat,2,function(x) predict(loess_lower_estimate,x))))
  loess_lower_estimate_means <- colMeans(ind_VEs)
  ind_VEs_upper_estimate <- 100*(1-exp(apply(antibody_BAU_mat,2,function(x) predict(loess_upper_estimate,x))))
  loess_upper_estimate_means <- colMeans(ind_VEs)
  omicron_means <- cbind(colMeans(ind_VEs_estimate),colMeans(ind_VEs_lower_estimate),colMeans(ind_VEs_upper_estimate))
  rm(ind_VEs_estimate,ind_VEs_lower_estimate,ind_VEs_upper_estimate)
  
  #######################
  # Predicted risk for new individual outside of trial with 0,1,2,3 observations
  
  #########
  # Simulate antibody data for a new individual from the posterior
  # Draw a random indicator that they test positive
  coef_draw_logmI <- coef(logmI) + chol_coef_logmI%*%rnorm(length(coef(logmI)))
  # plogis is the logistic function - so this is a prediction from the imputation model including uncertainty about the regression coefficients
  infection_draw <- rbinom(n=n_new_ind,size=1,prob=plogis(mmlogm%*%coef_draw_logmI))
  # Draw a random indicator that they are symptomatic *given* they tested positive
  coef_draw_logmS <- coef(logmS) + chol_coef_logmS%*%rnorm(length(coef(logmS)))
  symp_given_inf_draw <- rbinom(n=n_new_ind,size=1,prob=plogis(mmlogm%*%coef_draw_logmS))
  # Draw a NelsonAalen value by predicting from the NA linear model *given* values of prim and pos
  mmlinm[,"cor2dose_outcome_pos_primPositive"] <- infection_draw
  mmlinm[,"cor2dose_outcome_pos_primPrimary"] <- infection_draw*symp_given_inf_draw
  # Draw from asymptotic distribution of parameter MLEs
  coef_draw_linmNA <- coef(linmNA) + chol_coef_linmNA%*%rnorm(length(coef(linmNA)))
  # Calculate residuals associated with draw
  resid_linmNA <- linmNA$model$NelsonAalen - model.matrix(linmNA)%*%coef_draw_linmNA
  # And predicted values for new individuals (no error)
  predict_linmNA <- mmlinm%*%coef_draw_linmNA
  # Add a randomly drawn residual to predicted value to estimate with error. Set to 0 if negative as NA estimate is non-negative
  NA_draw <- pmax(predict_linmNA+sample(resid_linmNA,size=n_new_ind,replace=T),0)
  # Generate their outcome
  new_rows_tf$cor2dose_outcome_pos_primPositive[1:n_new_ind] <- (infection_draw-mean(data_long$X_l[,"cor2dose_outcome_pos_primPositive"]))/sd(data_long$X_l[,"cor2dose_outcome_pos_primPositive"])
  new_rows_tf$cor2dose_outcome_pos_primPrimary[1:n_new_ind] <- (infection_draw*symp_given_inf_draw-mean(data_long$X_l[,"cor2dose_outcome_pos_primPrimary"]))/sd(data_long$X_l[,"cor2dose_outcome_pos_primPrimary"])
  new_rows_tf$NelsonAalen <- (NA_draw-mean(data_long$X_l[,"NelsonAalen"]))/sd(data_long$X_l[,"NelsonAalen"])
  
  # Simulate the random intercepts and slopes with uncertainty
  stdnorms <- array(rnorm(2*n_new_ind*n_error_samples),dim=c(n_error_samples,n_new_ind,2))
  # Sample values of each parameter from N(MLE,Var(MLE))
  eta_0_tf_sim <- stdnorms[,,1]
  eta_1_tf_sim <- stdnorms[,,1]*long_out["rho_tf"] + stdnorms[,,2]*sqrt(1-(long_out["rho_tf"])^2)
  rm(stdnorms)
  # a_0 = x*beta_0+alpha_0 + eta_0*tau_0
  a_0_tf_sim <- sweep(eta_0_tf_sim*long_out["tau_0_tf"],2,long_out[beta_l0s_tf_cols]%*%t(new_rows_tf) + long_out["alpha_0_tf"],FUN="+")
  a_0_tf_sim
  # a_1 = x*beta_1+alpha_1 + eta_1*tau_1
  a_1_tf_sim <- -exp(sweep(eta_1_tf_sim*long_out["tau_1_tf"],2,long_out[beta_l1s_tf_cols]%*%t(new_rows_tf) + long_out["alpha_1_tf"],FUN="+"))
  rm(eta_1_tf_sim,eta_0_tf_sim)
  a_0_sim <- a_0_tf_sim*sd(data_long$log_y)+mean(data_long$log_y)
  a_1_sim <- a_1_tf_sim*sd(data_long$log_y)/sd(data_long$t_l)
  rm(a_0_tf_sim,a_1_tf_sim)
  VE_vs_time_new_inds <- apply(100*(1-exp(gamma*exp(sweep(outer(a_1_sim,ds_t0),c(1,2),a_0_sim,FUN="+"))+zeta)),c(3,2),mean)
  dim(VE_vs_time_new_inds)
  
  out <- cbind(VE_vs_time,antibody_quantiles,omicron_means,VE_vs_time_new_inds)
  return(c(out))
}
VE_vs_time_all <- vapply(seq_len(nchains),sfLapply_VE_vs_time_analysis,FUN.VALUE = array(0,c(length(parameters_list[[1]]),nts,outputs_dim)))
VE_vs_time_all <- aperm(VE_vs_time_all,c(1,4,2,3))
dimnames(VE_vs_time_all)[[2]] <- paste0("chains:",1:nchains)
runtimes$finished_loop <- Sys.time();
rm(parameters_list)

######################################################
# Post snowfall analysis
# Results split by chain
VE_vs_time_chains <- abind(apply(VE_vs_time_all[,,,dimnames_meanVE],c(2,3,4),mean),
                           aperm(apply(VE_vs_time_all[,,,dimnames_meanVE],c(2,3,4),quantile,c(0.5,0.025,0.975)),c(2,3,4,1)),along=4)

antibody_quantiles_chains <- abind(apply(VE_vs_time_all[,,,dimnames_antibody_quantiles],c(2,3,4),mean),
                           aperm(apply(VE_vs_time_all[,,,dimnames_antibody_quantiles],c(2,3,4),quantile,c(0.5,0.025,0.975)),c(2,3,4,1)),along=4)
# Change the dimension names - antibody times were at ds_t0 + 7, not the same as VE_times later
dimnames(antibody_quantiles_chains)[[2]] <- atsb
VE_vs_time_omicron_chains <- abind(apply(VE_vs_time_all[,,,dimnames_omicron],c(2,3,4),mean),
                           aperm(apply(VE_vs_time_all[,,,dimnames_omicron],c(2,3,4),quantile,c(0.5,0.025,0.975)),c(2,3,4,1)),along=4)

VE_vs_time_new_ind_chains <- abind(apply(VE_vs_time_all[,,,dimnames_new_ind],c(2,3,4),mean),
                                   aperm(apply(VE_vs_time_all[,,,dimnames_new_ind],c(2,3,4),quantile,c(0.5,0.025,0.975)),c(2,3,4,1)),along=4)

# Aggregated results
VE_vs_time <- abind(apply(VE_vs_time_all[,,,dimnames_meanVE],c(3,4),mean),
                           aperm(apply(VE_vs_time_all[,,,dimnames_meanVE],c(3,4),quantile,c(0.5,0.025,0.975)),c(2,3,1)),along=3)

antibody_quantiles <- abind(apply(VE_vs_time_all[,,,dimnames_antibody_quantiles],c(3,4),mean),
                    aperm(apply(VE_vs_time_all[,,,dimnames_antibody_quantiles],c(3,4),quantile,c(0.5,0.025,0.975)),c(2,3,1)),along=3)
dimnames(antibody_quantiles)[[1]] <- atsb

VE_vs_time_omicron <- abind(apply(VE_vs_time_all[,,,dimnames_omicron],c(3,4),mean),
                                   aperm(apply(VE_vs_time_all[,,,dimnames_omicron],c(3,4),quantile,c(0.5,0.025,0.975)),c(2,3,1)),along=3)

VE_vs_time_new_ind <- abind(apply(VE_vs_time_all[,,,dimnames_new_ind],c(3,4),mean),
                                   aperm(apply(VE_vs_time_all[,,,dimnames_new_ind],c(3,4),quantile,c(0.5,0.025,0.975)),c(2,3,1)),along=3)


runtimes$post_snowfall_analysis <- Sys.time();

######################################################
# Find the difference between each time when we "checked the clock" in the file, to know how long each step took
runtimes_diff <- mapply(difftime,runtimes[2:length(runtimes)],runtimes[1:(length(runtimes)-1)])
names(runtimes_diff) <- names(runtimes)[-1]
output <- list("effect_summary" = effect_summary,
               "VE_summary_chains" = VE_summary_chains,
               "VE_summary" = VE_summary,
               "antibody_at_VE_summary" = antibody_at_VE_summary,
               "antibody_at_VE_summary_chains" = antibody_at_VE_summary_chains,
               "direct_effect_CI_chains" = direct_effect_CI_chains,
               "direct_effect_CI" = direct_effect_CI,
               "VE_indirect_summary" = VE_indirect_summary,
               "VE_indirect_summary_chains" = VE_indirect_summary_chains,
               "antibody_quantiles_chains" = antibody_quantiles_chains,
               "antibody_quantiles" = antibody_quantiles,
               "VE_vs_time_chains" = VE_vs_time_chains,
               "VE_vs_time" = VE_vs_time,
               "VE_vs_time_omicron_chains" = VE_vs_time_omicron_chains,
               "VE_vs_time_omicron" = VE_vs_time_omicron,
               "new_data" = new_data,
               "logmI" = logmI, "logmS" = logmS, "linmNA" = linmNA,
               "VE_vs_time_new_ind_chains" = VE_vs_time_new_ind_chains,
               "VE_vs_time_new_ind" = VE_vs_time_new_ind,
               "data_long" = data_long,
               "pred_antibody_values_28_hist" = pred_antibody_values_28_hist,
               "pred_antibody_values_182_hist" = pred_antibody_values_182_hist,
               "obs_antibody_values_hist" = obs_antibody_values_hist,
               "tens_BAU" = tens_BAU, "ds_t0" = ds_t0, "atsb" = atsb, "VEtsb" = VEtsb,
               "runtimes" = runtimes,
               "runtimes_diff" = runtimes_diff)
# Save the output
saveRDS(output,paste0(output_directory,"/Posterior_estimates_cox_infection_model_long_",long_file_name,cox_model_namelong,model_name,".RDS"))
