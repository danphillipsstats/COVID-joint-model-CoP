# This file runs the longitudinal Bayesian mixed effects log-linear antibody model.
# This is the first stage in the two-stage joint model.

t1 <- Sys.time()
print(t1)
# To run code on the server.

#############################################################################
# Define the directories in which files are saved - these may need adjusting for other users.
directory_correlates <- getwd()
directory_correlates_project <- paste0(directory_correlates,"/3_Programs")
output_directory <- paste0(directory_correlates,"/4_Output")

#load packages
library(rstan)
library(survival)

#######
# Load user-written functions
source(paste0(directory_correlates_project,"/functions.R"))
###############
# Read the objects joint_correlates and long_correlates_S
# Set their factors to be correct
source(paste0(directory_correlates_project,"/read_data_set_factors.R"))

# Defining variables

# Select the necessary comments
joint_correlates <- joint_correlates[,c("As_vaccinated_arm_2","sc_repeat_pid","sc_age","age_group","cor2dose_schedule","cor2dose_interval","cor2dose_non_white","cor2dose_bmi_geq_30","cor2dose_comorbidities","cor2dose_hcw_status","cor2dose_primary","cor2dose_positive","start_time","end_time","sc_gender","cor2dose_primeschedule","cor2dose_outcome_pos_prim")] 

################
# Removing outlying anti-spike IgG observations
# Exclude the exceptionally small results as they appear to be outliers and this appears to be preventing good mixing in the model (rho = -1)
long_correlates_S <- long_correlates_S[which(long_correlates_S$log_S>7),]
# Exclude the exceptionally large result at PB90
long_correlates_S <- long_correlates_S[which(!(long_correlates_S$antibody_time>50 & long_correlates_S$antibody_time<130 & long_correlates_S$log_S>13)),]
# Include only the ChAdOx1 group' in the antibody model's antibody data
long_correlates_S <- long_correlates_S[which(long_correlates_S$As_vaccinated_arm_2=="ChAdOx1"),]

###############
# Create the data list for stan
# Check there are sufficient cases for each longitudinal covariate to include them in the model
long_covariates <- c("age_group","sc_gender","cor2dose_non_white","cor2dose_comorbidities","cor2dose_bmi_geq_30","cor2dose_interval","cor2dose_primeschedule","cor2dose_positive")
for(i in 1:length(long_covariates)){
  print(table(joint_correlates[which(joint_correlates$cor2dose_positive=="Positive"),c(long_covariates[i],"As_vaccinated_arm_2")]))
}
# Check there are sufficient longitudinal observations for each covariate
match_long_joint <- match(long_correlates_S$sc_repeat_pid,joint_correlates$sc_repeat_pid)
long_correlates_S[,long_covariates] <- joint_correlates[match_long_joint,long_covariates]
for(i in 1:length(long_covariates)){
  print(table(long_correlates_S[,c(long_covariates[i],"As_vaccinated_arm_2")]))
} # These all look like they have plenty of observations
# We can't include schedule on the other hand - only 2 observations for LDLD
print(table(long_correlates_S[,c("cor2dose_schedule","As_vaccinated_arm_2")]))

# Create a variable for the end_time in terms of calendar time
joint_correlates$end_cal_time <- joint_correlates$start_time + joint_correlates$end_time
# A simple model to estimate the cumulative hazard using the Nelson-Aalen estimate
basic_basehaz <- survfit(Surv(start_time,end_cal_time,cor2dose_positive)~As_vaccinated_arm_2,id=sc_repeat_pid,data=joint_correlates)
haz <- basic_basehaz$n.event[,2]/basic_basehaz$n.risk[,1]
haz <- as.data.frame(haz)
colnames(haz) <- c("cumhaz")
# Set the arm variable
haz$As_vaccinated_arm_2 <- NA
haz$As_vaccinated_arm_2[1:basic_basehaz$strata[1]] <- levels(joint_correlates$As_vaccinated_arm_2)[1]
haz$As_vaccinated_arm_2[basic_basehaz$strata[1]+(1:basic_basehaz$strata[2])] <- levels(joint_correlates$As_vaccinated_arm_2)[2]
haz$time <- basic_basehaz$time
# Exclude unnecessary rows
discrete_haz <- haz[which(haz[,1]>0),]
# Create a NelsonAalen estimator variable. This is the same as the NelsonAalen estimator applied to the two arms separately.
joint_correlates$NelsonAalen <- NA
for (i in 1:nrow(joint_correlates)){
  joint_correlates$NelsonAalen[i] <- sum(discrete_haz$cumhaz[which(discrete_haz$time>joint_correlates$start_time[i] &
                                                                     discrete_haz$time<=joint_correlates$end_cal_time[i] &
                                                                     discrete_haz$As_vaccinated_arm_2==joint_correlates$As_vaccinated_arm_2[i])])
}
# t_center = 7 means centered at day 28 post second dose as antibody_time=x means a reading at 21+x days post boost.
data_long <- stan_long_model(long_data = long_correlates_S[which(long_correlates_S$As_vaccinated_arm_2=="ChAdOx1"),c("sc_repeat_pid","log_S","antibody_time")],
                             log_y_var = "log_S", obs_time_var = "antibody_time",
                             individual_data = joint_correlates[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1"),],
                             longitudinal_covariates = c("age_group","sc_gender","cor2dose_non_white","cor2dose_comorbidities","cor2dose_bmi_geq_30","cor2dose_interval","cor2dose_hcw_status","cor2dose_primeschedule","cor2dose_outcome_pos_prim","NelsonAalen"),
                             ind_var = "sc_repeat_pid",
                             t_center = 7)
print(paste0("mean(log_y) = ",mean(data_long$log_y)))
# Parse the stan model
long_model = rstan::stan_model(file=paste0(directory_correlates_project,'/Longitudinal_hierarchical_model_normalised_ri_rs_exp_slope_t.stan'))

warmup_iterations = 5e4
sampling_iterations = 15e3 + warmup_iterations #best to use 1e3 or higher

# Stan sampling step - run in parallel with 4 cores (can be done on a local computer but will take a very long time with this many iterations)
long_out = rstan::sampling(
  object = long_model,
  data = data_long,
  chains = 4,
  cores = 4,
  iter = sampling_iterations,
  warmup = warmup_iterations,
  refresh = sampling_iterations/10, #show an update @ each %10
  seed = 1000,
  include = FALSE, # Exclude the following parameters from the model output
  pars = c("a_0_tf","a_1_tf","eta_0_tf","eta_1_tf"))

file_name <- "ri_rs_7inc_pos_prim_exp_slope_t_all_covariates_t0PB28"
name <- paste0("Correlates_long_",file_name)

# Summaries of stan model output
summary_out <- summary(long_out)$summary
print(max(summary_out[,"Rhat"]))
print(min(summary_out[,"n_eff"]))

long_out <- extract(long_out,permuted=F)
# Save the random intercepts and slopes
long_out_a <- long_out[,,which(is.element(substring(dimnames(long_out)$parameters,1,4),c("a_0[","a_1[")))]
saveRDS(list(long_out_a,data_long,summary_out),paste0(output_directory,"/",name,"_reffects.RDS"))
# Save the other non-random effect parameters in a separate file
long_out <- long_out[,,which(!is.element(substring(dimnames(long_out)$parameters,1,2),c("a_","et")))]
saveRDS(list(long_out,data_long,summary_out),paste0(output_directory,"/",name,"_noreffects.RDS"))
t2 <- Sys.time()
print("Stan model and saving complete:")
print(difftime(t2,t1))