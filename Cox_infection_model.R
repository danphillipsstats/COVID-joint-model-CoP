# This file runs a Cox model based on the output from the Longitudinal_antibody_model.R file
# The second stage in the two-stage joint model.
# Models the time until testing positive for COVID-19, based on multiply imputed predicted antibody levels.

t_init <- Sys.time()
print(t_init)
# To run code on the server.

##############################################################################
# Cox infection model
##############################################################################
# Define the directories in which files are saved - these may need adjusting for other users.
directory_correlates <- getwd() 
directory_correlates_project <- paste0(directory_correlates,"/3_Programs")
output_directory <- paste0(directory_correlates,"/4_Output")

#load packages
library(rstan)
library(survival)
library(tidyr)
library(parallel)
library(snowfall)

#######
# Load user-written functions
source(paste0(directory_correlates_project,"/functions.R"))
##########################
# Set what the outcome for this analysis will be 
# (pos = all COVID-19 infections, prim = primary symptomatic COVID-19 infections)
event_outcome <- c("pos","prim")[2] # To be run with both infection outcomes separately

# Get the intercepts and slopes
set.seed(1)
# Name of file from Longitudinal_antibody_model.R
file_name <- "ri_rs_7inc_pos_prim_exp_slope_t_all_covariates_t0PB28_resid"
# Read the output from Longitudinal_antibody_model.R
long_out_rs <- readRDS(paste0(output_directory,"/Correlates_long_",file_name,"_reffects.RDS"))
# The list of objects put into the stan model
data_long <- long_out_rs[[2]]
# The random intercepts and slopes outputted from the stan model
long_out_rs <- long_out_rs[[1]]
dim(long_out_rs)
# Find the population parameters
a_0_array <- long_out_rs[,,which(substring(dimnames(long_out_rs)$parameters,1,4)=="a_0[")]
a_1_array <- long_out_rs[,,which(substring(dimnames(long_out_rs)$parameters,1,4)=="a_1[")]
rm(long_out_rs)
# Undo the time shift by t_center + 7 (so that t0 is now PB14, 7 days before t0 PB21 in the survival data)
# This means the hazard at time t will relate to the antibody level 7 days prior.
a_0_array <- a_0_array - (data_long$t_center + 7) * a_1_array
n <- dim(a_0_array)[3]
nsamples <- dim(a_0_array)[1]
nchains <- dim(a_0_array)[2]
print("Created a_0_array and a_1_array")

# Do the above when the new longitudinal model has been run.
###############
# Read joint_correlates
# Set their factors to be correct
source(paste0(directory_correlates_project,"/read_data_set_factors.R"))
joint_correlates$vaccine_allocation <- as.numeric(joint_correlates$As_vaccinated_arm_2=="ChAdOx1")

#######
# Prepare the times
surv_fix_columns <- c("sc_repeat_pid","As_vaccinated_arm_2","age_group","sc_age","cor2dose_hcw_status","cor2dose_non_white","cor2dose_bmi_geq_30","cor2dose_comorbidities","cor2dose_primeschedule","cor2dose_schedule","cor2dose_interval","sc_gender","geographic_big_region","geographic_region","site")
if (event_outcome == "prim"){ # if event is primary symptomatic COVID-19 infection
  times <- joint_correlates[,c(surv_fix_columns,"cor2dose_primary_ind","start_time","end_time")]
} else if (event_outcome == "pos"){ # if event is any COVID-19 infection
  times <- joint_correlates[,c(surv_fix_columns,"cor2dose_positive_ind","start_time","end_time")]
}
# Change end_time to be calendar time of ending the at-risk period, instead of time since the start of the at-risk period
times$end_time <- times$start_time+times$end_time
# Use tmerge
# Create survival dataset
if (event_outcome == "prim"){
  times$event <- as.numeric(times$cor2dose_primary_ind==1)
} else if (event_outcome == "pos"){
  times$event <- as.numeric(times$cor2dose_positive_ind==1)
}
times2 <- tmerge(times,times,id=sc_repeat_pid,endpt=event(end_time,event),tstart = start_time, tstop = end_time)
# Create antibody data
if (event_outcome == "prim"){ # if event is primary symptomatic COVID-19 infection
  event_times <- times$end_time[which(times$cor2dose_primary_ind==1)]
} else if (event_outcome == "pos"){ # if event is any COVID-19 infection
  event_times <- times$end_time[which(times$cor2dose_positive_ind==1)]
}
event_times <- c(0,unique(event_times[order(event_times)]))

##############################################################################
# Create sample dataset
a_0_sample <- a_0_array[1,1,]
a_1_sample <- a_1_array[1,1,]
# Prepare wide antibody data
antibody_data <- exp(a_0_sample%*%t(rep(1,length(event_times))) + a_1_sample%*%t(event_times))
colnames(antibody_data) <- event_times
antibody_data <- as.data.frame(antibody_data)
surv_time_id_columns <- c("sc_repeat_pid","start_time","end_time")
antibody_data[,surv_time_id_columns] <- joint_correlates[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1"),surv_time_id_columns]
# Change end_time to be calendar time of ending the at-risk period, instead of time since the start of the at-risk period
antibody_data$end_time <- antibody_data$end_time + antibody_data$start_time
# Transform antibody data from wide to long
last_time <- as.character(event_times[length(event_times)])
long_antibody_data <- gather(antibody_data,key=time,value=antibody,"0":last_time,factor_key=F)
long_antibody_data$time <- as.numeric(long_antibody_data$time)
long_antibody_data <- long_antibody_data[order(long_antibody_data$sc_repeat_pid),]
# Merge time-varying antibody data with survival dataset
surv_antibody_data <- tmerge(times2,long_antibody_data,id=sc_repeat_pid,antibody=tdc(time,antibody))
surv_antibody_data <- surv_antibody_data[,c(surv_fix_columns,"tstart","tstop","endpt","antibody")]
names(surv_antibody_data)[names(surv_antibody_data)=="tstart"] <- "start_time"
names(surv_antibody_data)[names(surv_antibody_data)=="tstop"] <- "end_time"
names(surv_antibody_data)[names(surv_antibody_data)=="endpt"] <- "event"

print("Created surv_antibody_data")

################
# Create an array of intercepts and slopes to be imputed for each iteration
length(joint_correlates$sc_repeat_pid[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1")])
length(joint_correlates$sc_repeat_pid)
dimnames(a_0_array)$parameters <- joint_correlates$sc_repeat_pid[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1")]
dimnames(a_1_array)$parameters <- joint_correlates$sc_repeat_pid[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1")]
# Create a list of 4 sublists for each chain
# In the sublists entry i corresponds to the ith draw from posterior for a_0, a_1
# Each entry contains a matrix where row 1 is a_0 (intercept) for each individual and row 2 is a_1 (gradient)
a_01_mat_list <- lapply(seq_len(nchains), function(chain) lapply(seq_len(nsamples),function(i) rbind(a_0_array[i,chain,],a_1_array[i,chain,])))
print(object.size(a_01_mat_list),units="GB")
rm(a_0_array,a_1_array)

print("Created a_01_mat_list")

# Set antibody levels to 0 for the control individuals
surv_antibody_data$antibody[which(surv_antibody_data$As_vaccinated_arm_2=="Control")] <- 0
vacc_group_ind <- which(surv_antibody_data$As_vaccinated_arm_2=="ChAdOx1")

# Create initial templates of outputs to fill in during the loop
a_0_sample_times <- a_01_mat_list[[1]][[1]][1,as.character(surv_antibody_data$sc_repeat_pid[vacc_group_ind])]
a_1_sample_times <- a_01_mat_list[[1]][[1]][2,as.character(surv_antibody_data$sc_repeat_pid[vacc_group_ind])]
surv_antibody_data$antibody[vacc_group_ind] <- exp(a_0_sample_times + surv_antibody_data$end_time[vacc_group_ind]*a_1_sample_times)
# Includes an effect due to antibodies, as well as a direct effect due to vaccination (As_vaccinated_arm_2).
# Includes covariates age, sex, ethnicity, comorbidity, BMI, healthcare worker
# all of which may affect the risk of infection independently of vaccination (i.e. for both vaccinated and control individuals)
cox_model_formula <- Surv(start_time,end_time,event)~antibody+age_group+sc_gender+cor2dose_non_white+cor2dose_comorbidities+cor2dose_bmi_geq_30+cor2dose_hcw_status+strata(site)
#####
# Create a model.matrix for use in later Output analysis
joint_correlates_mat <- joint_correlates
if (event_outcome == "prim"){ # if event is primary symptomatic COVID-19 infection
  joint_correlates_mat$event <- joint_correlates_mat$cor2dose_primary_ind
} else if (event_outcome == "pos"){ # if event is any COVID-19 infection
  joint_correlates_mat$event <- joint_correlates_mat$cor2dose_positive_ind
}
joint_correlates_mat$antibody <- 0
joint_correlates_mat$antibody[which(joint_correlates_mat$As_vaccinated_arm_2=="ChAdOx1")] <- 1 # Just set to one for ease of understanding interaction terms
joint_correlates_mat$end_time <- joint_correlates_mat$start_time + joint_correlates_mat$end_time # As end_time should be end time not length of time at risk
cox_model_mat <- model.matrix(cox_model_formula,data=joint_correlates_mat)
# Run a cox model, to create placeholders
cox_model <- try(coxph(cox_model_formula,
                     data=surv_antibody_data, id=sc_repeat_pid))
npred <- length(c(cox_model$coefficients,cox_model$var))
print("Created cox_model")

# Set up the cluster using snowfall package
nwork <- length(parallelly::availableWorkers()) # number of parallel workers

# Check the object size in memory
obj_size <- format(object.size(a_01_mat_list),units="GB",standard="SI")
obj_size
paste0("Object size if copied to all nodes:")
paste0(as.numeric(substring(obj_size,first=1,last=gregexpr(pattern="G",obj_size)[[1]]-1))*nwork," GB")

t3 <- Sys.time()
print("Starting apply...")
# Function runs the Cox model using antibody data from matrix a_01_mat
# Where row 1 of a_01_mat is a_0 (intercept) and row 2 is a_1 (gradient)
# Columns are all individuals.
# Returns a vector containing first the MLEs and then the asymptotic covariance matrix of Cox parameters
cox_VE_antibody <- function(a_01_mat){
  surv_antibody_data$antibody[vacc_group_ind] <- exp(a_01_mat[1,as.character(surv_antibody_data$sc_repeat_pid[vacc_group_ind])] + surv_antibody_data$start_time[vacc_group_ind]*a_01_mat[2,as.character(surv_antibody_data$sc_repeat_pid[vacc_group_ind])])
  cox_model <- try(coxph(cox_model_formula,
                       data=surv_antibody_data, id=sc_repeat_pid))
  return(try(c(cox_model$coefficients,cox_model$var)))
}

sfSapply_a_01_mat_list <- function(chain, my_seed = 1234){
  a_01_mat_sub_list <- a_01_mat_list[[chain]]
  # Begin a cluster in snowfall package
  cluster <- snowfall::sfInit(nwork,type="SOCK", parallel=T)   
  # Load the relevant objects in the parallel workspaces
  snowfall::sfExport("surv_antibody_data","vacc_group_ind","cox_model_formula")
  # Load the survival library
  snowfall::sfLibrary(survival)
  # Set up the random number generation
  snowfall::sfClusterSetupRNG(seed=my_seed)
  
  print(paste("Number of workers:",nwork))
  # Runs the function on each element of the list in parallel
  # That is, runs a cox model for each imputation of antibody levels.
  cox_model_pred_var <- t(snowfall::sfSapply(x = a_01_mat_sub_list, fun = cox_VE_antibody))
  snowfall::sfStop()
  return(cox_model_pred_var)
} # vapply ensures the outputted array is of the right dimensions
cox_model_pred_var <- vapply(seq_len(nchains),sfSapply_a_01_mat_list,FUN.VALUE = array(0,c(nsamples,npred)))

# Define the column names
variance_names <- paste(rep(names(cox_model$coefficients), each = length(names(cox_model$coefficients))),names(cox_model$coefficients), sep = ".")
dim(cox_model_pred_var)
dimnames(cox_model_pred_var)[[2]] <- c(names(cox_model$coefficients),variance_names)
# aperm reorders the dimensions in the array
output <- list("cox_model_pred" = aperm(cox_model_pred_var[,names(cox_model$coefficients),],c(1,3,2)),
               "cox_model_var" = aperm(cox_model_pred_var[,variance_names,],c(1,3,2)),
               "cox_model_mat" = cox_model_mat)

cox_model_name <- paste0(event_outcome,"_site_parallel_simple_7inc_nodirect")
saveRDS(output,paste0(output_directory,"/Cox_infection_model_long_",file_name,"_",cox_model_name,".RDS"))
pryr::mem_used()
t_final <- Sys.time()
print(difftime(t_final,t_init))
