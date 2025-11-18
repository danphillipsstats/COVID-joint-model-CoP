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
directory_correlates <- dirname(dirname(getwd()))
directory_correlates_project <- paste0(directory_correlates,"/3_Programs/COVID-joint-model-CoP")
output_directory <- paste0(directory_correlates,"/4_Output")
conv_factor <- 0.00645

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
file_name <- "ri_rs_7inc_pos_prim_exp_slope_t_all_covariates_t0PB28"
# Read the output from Longitudinal_antibody_model.R
long_out_rs <- readRDS(paste0("/data/dragon110/shug5415/Documents/Research/Correlates_joint_modelling/4_Output/Correlates_long_",file_name,"_reffects.RDS"))
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

thin <- seq(from=0,to=15000,length.out=6)[-1]
a_0_array <- a_0_array[thin,,]
a_1_array <- a_1_array[thin,,]
dim(a_1_array)

a_0_mat <- apply(a_0_array,3,c)
a_1_mat <- apply(a_1_array,3,c)
m <- nrow(a_0_mat)
inds_mat <- matrix(1:n,nrow=m,ncol=n,byrow=T)

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
a_0_sample <- c(a_0_mat) # Includes values for each individual from each imputed dataset stacked on top of each other
a_1_sample <- c(a_1_mat) # So m values for individual 1, then m values for individual 2 etc.
# Prepare wide antibody data
antibody_data <- exp(a_0_sample%*%t(rep(1,length(event_times))) + a_1_sample%*%t(event_times))
colnames(antibody_data) <- event_times
antibody_data <- as.data.frame(antibody_data)
surv_time_id_columns <- c("sc_repeat_pid","start_time","end_time")
antibody_data[,surv_time_id_columns] <- joint_correlates[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1"),surv_time_id_columns][c(inds_mat),]
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

# Set antibody levels to 0 for the control individuals
surv_antibody_data$antibody[which(surv_antibody_data$As_vaccinated_arm_2=="Control")] <- 0
vacc_group_ind <- which(surv_antibody_data$As_vaccinated_arm_2=="ChAdOx1")

# Make the directory
# if (!dir.exists(paste0(plot_directory,"/VE_plots/residuals2"))){dir.create(paste0(plot_directory,"/VE_plots/residuals2"))}

#####
# Includes an effect due to antibodies, as well as a direct effect due to vaccination (As_vaccinated_arm_2).
# Includes covariates age, sex, ethnicity, comorbidity, BMI, healthcare worker
# all of which may affect the risk of infection independently of vaccination (i.e. for both vaccinated and control individuals)
cox_model_direct_formula <- Surv(start_time,end_time,event)~antibody+As_vaccinated_arm_2+age_group+sc_gender+cor2dose_non_white+cor2dose_comorbidities+cor2dose_bmi_geq_30+cor2dose_hcw_status+strata(site)
cox_model_direct <- try(coxph(cox_model_direct_formula,
                     data=surv_antibody_data, id=sc_repeat_pid, weights = rep(1/m,nrow(surv_antibody_data))))

# Cox-Snell residuals
cox_snell_direct_ant <- surv_antibody_data$event - resid(cox_model_direct)
cox_snell_direct_ant_NA <- survfit(Surv(cox_snell_direct_ant,surv_antibody_data$event)~1)
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Cox-Snell_antibody_direct.png"),width = 1400,height=1000, pointsize=23)
plot(cox_snell_direct_ant_NA, fun="cumhaz", mark.time=F, main = "Cox-Snell residual plot for antibody model with direct effect")
abline(0,1,col=2)
dev.off()

# Model with no direct effect but antibody effect only
cox_model_formula <- Surv(start_time,end_time,event)~antibody+age_group+sc_gender+cor2dose_non_white+cor2dose_comorbidities+cor2dose_bmi_geq_30+cor2dose_hcw_status+strata(site)
cox_model <- try(coxph(cox_model_formula,
                                data=surv_antibody_data, id=sc_repeat_pid, weights = rep(1/m,nrow(surv_antibody_data))))

cox_model_interactions_formula <- Surv(start_time,end_time,event)~antibody + antibody:age_group+age_group+sc_gender+cor2dose_non_white+cor2dose_comorbidities+cor2dose_bmi_geq_30+cor2dose_hcw_status+strata(site)
cox_model_interactions <- try(coxph(cox_model_interactions_formula,
                                    data=surv_antibody_data, id=sc_repeat_pid, weights = rep(1/m,nrow(surv_antibody_data))))

# Cox-Snell residuals
cox_snell_ant <- surv_antibody_data$event - resid(cox_model)
cox_snell_ant_NA <- survfit(Surv(cox_snell_ant,surv_antibody_data$event)~1)
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Cox-Snell_antibody_only.png"),width = 1400,height=1000, pointsize=23)
plot(cox_snell_ant_NA, fun="cumhaz", mark.time=F, main = "Cox-Snell residual plot for antibody model without direct effect")
abline(0,1,col=2)
dev.off()

# Model using log(1 + antibodies)
surv_logantibody_data <- surv_antibody_data
surv_logantibody_data$logantibody <- 0
surv_logantibody_data$logantibody[which(surv_logantibody_data$antibody!=0)] <- log(1+surv_logantibody_data$antibody*conv_factor)[which(surv_logantibody_data$antibody!=0)]

cox_model_formula_logant <- Surv(start_time,end_time,event)~logantibody+As_vaccinated_arm_2+age_group+sc_gender+cor2dose_non_white+cor2dose_comorbidities+cor2dose_bmi_geq_30+cor2dose_hcw_status+strata(site)
cox_logantibody_model <- try(coxph(cox_model_formula_logant,
                                   data=surv_logantibody_data, id=sc_repeat_pid, weights = rep(1/m,nrow(surv_logantibody_data))))
cox_model_formula_logant_nodirect <- Surv(start_time,end_time,event)~logantibody+age_group+sc_gender+cor2dose_non_white+cor2dose_comorbidities+cor2dose_bmi_geq_30+cor2dose_hcw_status+strata(site)
cox_logantibody_model_nodirect <- try(coxph(cox_model_formula_logant_nodirect,
                                   data=surv_logantibody_data, id=sc_repeat_pid, weights = rep(1/m,nrow(surv_logantibody_data))))

# Cox-Snell residuals for log antibody model
cox_snell_logant <- surv_antibody_data$event - resid(cox_logantibody_model)
cox_snell_logant_NA <- survfit(Surv(cox_snell_logant,surv_antibody_data$event)~1)
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Cox-Snell_logantibody_direct.png"),width = 1400,height=1000, pointsize=23)
plot(cox_snell_logant_NA, fun="cumhaz", mark.time=F, main = "Cox-Snell residual plot for log antibody model")
abline(0,1,col=2)
dev.off()
write.csv(round(cbind(AIC(cox_model_direct,cox_model,cox_logantibody_model,cox_logantibody_model_nodirect),BIC(cox_model_direct,cox_model,cox_logantibody_model,cox_logantibody_model_nodirect))[,c(1,2,4)],1),
          paste0(output_directory,"/IC_compare_COVID_jm.csv"))
# Both prefer antibody with no direct effect over log(1+A) with or without direct effect, or antibody with direct effect.


# Martingale residuals for the effect of antibodies
# Fit cox model to vaccine group only (as these are the individuals we have antibody observations for)
cox_model_formula_noant <- Surv(start_time,end_time,event)~age_group+sc_gender+cor2dose_non_white+cor2dose_comorbidities+cor2dose_bmi_geq_30+cor2dose_hcw_status+strata(site)
#####
# Run a cox model, to create placeholders
surv_antibody_data_ChAdOx1 <- surv_antibody_data[which(surv_antibody_data$As_vaccinated_arm_2=="ChAdOx1"),]
cox_model_noant <- try(coxph(cox_model_formula_noant,
                       data=surv_antibody_data_ChAdOx1, id=sc_repeat_pid, weights = rep(1/m,nrow(surv_antibody_data_ChAdOx1))))
mart_resid_ant <- resid(cox_model_noant,type="martingale")
# Plot all
ord <- order(surv_antibody_data_ChAdOx1$antibody)
mart_resid <- mart_resid_ant[ord]
antibody_for_mart <- surv_antibody_data_ChAdOx1$antibody[ord]
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Martingale_residuals.png"),width = 1400,height=1000, pointsize=23)
plot(mart_resid~antibody_for_mart, xlab = "Antibody", ylab = "Martingale residuals", cex=0.1, log="x")
lines(lowess(mart_resid~antibody_for_mart),lwd=2,col=2)
dev.off()
# restrict to antibodies <= 1e5
ord <- order(surv_antibody_data_ChAdOx1$antibody)
mart_resid <- mart_resid_ant[ord]
antibody_for_mart <- surv_antibody_data_ChAdOx1$antibody[ord]
ord2 <- which(antibody_for_mart<1e5)
antibody_for_mart <- antibody_for_mart[ord2]
mart_resid <- mart_resid[ord2]
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Martingale_residuals_zoomed.png"),width = 1400,height=1000, pointsize=23)
plot(mart_resid~antibody_for_mart, xlab = "Antibody", ylab = "Martingale residuals", cex=0.1)
lowess.out <- lowess(mart_resid~antibody_for_mart)
lines(lowess(mart_resid~antibody_for_mart),lwd=2,col=2)
dev.off()
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Martingale_residuals_log.png"),width = 1400,height=1000, pointsize=23)
plot(mart_resid~log(antibody_for_mart), xlab = "Antibody", ylab = "Martingale residuals with antibody on log scale", cex=0.1)
lowess.out <- lowess(mart_resid~log(antibody_for_mart))
lines(lowess(mart_resid~log(antibody_for_mart)),lwd=2,col=2)
dev.off()

# Schoenfeld residuals - antibody and direct effect model
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Schoenfeld_residuals_antibody_direct_ant.png"),width = 1400,height=1000, pointsize=23)
z <- cox.zph(cox_model_direct)
plot(z[1], main = "Schoenfeld residuals for antibody in antibody with direct effect model")
abline(h=0)
dev.off()
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Schoenfeld_residuals_antibody_direct_vacc.png"),width = 1400,height=1000, pointsize=23)
plot(z[2], main = "Schoenfeld residuals for vaccination in antibody with direct effect model")
abline(h=0)
dev.off()
# Schoenfeld residuals - antibody effect only model
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Schoenfeld_residuals_antibody_only_ant.png"),width = 1400,height=1000, pointsize=23)
z_noant <- cox.zph(cox_model)
plot(z_noant[1], main = "Schoenfeld residuals for antibody in antibody without direct effect model")
abline(h=0)
dev.off()
# Schoenfeld residuals
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Schoenfeld_residuals_logantibody_direct_ant.png"),width = 1400,height=1000, pointsize=23)
z_logant <- cox.zph(cox_logantibody_model)
plot(z_logant[1], main = "Schoenfeld residuals for log antibody in log antibody model")
abline(h=0)
dev.off()
jpeg(paste0(plot_directory,"/VE_plots/residuals2/Schoenfeld_residuals_logantibody_direct_vacc.png"),width = 1400,height=1000, pointsize=23)
plot(z_logant[2], main = "Schoenfeld residuals for vaccination in log antibody model")
abline(h=0)
dev.off()


cox_model_name <- paste0(event_outcome,"_site_parallel_simple_7inc")
saveRDS(output,paste0(output_directory,"/Cox_infection_model_long_",file_name,"_",cox_model_name,".RDS"))
pryr::mem_used()
t_final <- Sys.time()
print(difftime(t_final,t_init))
