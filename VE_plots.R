##############
# This file uses the output from the Cox model in Cox_infection_model.R
# and the posterior estimates in Posterior_estimates.R
# to calculate the plots, figures and results needed for the paper.

t_init <- Sys.time()
print(t_init)

antibody_type <- c("S","neuts")[1]
long_type <- c("7_day_incubator","0_day_incubator")[1]
long_file_name <- "ri_rs_7inc_pos_prim_exp_slope_t_all_covariates_t0PB28"
cox_model_name <- c("prim_site_parallel_simple_7inc","pos_site_parallel_simple_7inc")
{n.outcomes <- length(cox_model_name)
  outcome_type_short <- substring(cox_model_name,1,3)
  outcome_type <- rep("unknown",n.outcomes); outcome_title <- rep("Unknown",n.outcomes)
  for (i in 1:n.outcomes){
    if (outcome_type_short[i]=="pos"){outcome_type[i] <- "allpositive"; outcome_title[i] <- "All positive cases"}
    else if (outcome_type_short[i]=="pri"){outcome_type[i] <- "symptomatic"; outcome_title[i] <- "Primary symptomatic cases"}
  }
}
model_name <- "submission"
if(!(substring(long_file_name,1,5)=="neuts") == (antibody_type=="neuts")){warning("Your long_file_name and antibody_type do not match. You may be running a neuts analysis on S data or vice versa. Please double check.")}
# Directories
SLURM <- FALSE # Is this being run on SLURM?
if(SLURM){
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
  if (!dir.exists(paste0(plot_directory,"/VE_plots"))){
    dir.create(paste0(plot_directory,"/VE_plots"))
  }
  if (!dir.exists(paste0(plot_directory,"/VE_plots/",model_name))){
    dir.create(paste0(plot_directory,"/VE_plots/",model_name))
  }
} else{ # This file has code defining where the directories are saved. Not included in the GitHub
  read_dir_directory <- getwd() # Replace with where Read_directories_VE_plots.R is saved
  source(paste0(read_dir_directory,"/Read_directories_VE_plots.R"))
}

#load packages
library(rstan)
library(survival)
library(ggplot2)
require(RColorBrewer)

#######
# Load user-written functions
source(paste0(directory_correlates_project,"/functions.R"))
################################################################################
# Output from Posterior_estimates
out <- list()
for (i in 1:n.outcomes){
  assign(paste0("cox_model_name_",outcome_type[i]),paste0("_",cox_model_name[i],"_",model_name))
  assign(paste0("cox_stage_two_",outcome_type[i]),readRDS(paste0(output_directory,"/Cox_infection_model_long_",long_file_name,"_",cox_model_name[i],".RDS")))
  # Read the output from Posteior_estimates.R
  out[[i]] <- readRDS(paste0(output_directory,"/Posterior_estimates_cox_infection_model_long_",long_file_name,"_",cox_model_name[i],"_",model_name,".RDS"))
  # Save the entries which don't depend on outcome separately
  data_long <- out[[i]]$data_long
  pred_antibody_values_28_hist <- out[[i]]$pred_antibody_values_28_hist
  pred_antibody_values_182_hist <- out[[i]]$pred_antibody_values_182_hist
  tens_BAU <- out[[i]]$tens_BAU
  tens_BAU2 <- tens_BAU[which(tens_BAU>1)]
  antibody_quantiles <- out[[i]]$antibody_quantiles
  VE_vs_time_omicron <- out[[i]]$VE_vs_time_omicron
}
names(out) <- outcome_type

###############
# Read joint_correlates and long_correlates
# Set their factors to be correct
source(paste0(directory_correlates_project,"/read_data_set_factors.R"))
if(antibody_type=="S"){
  # Create long_correlates to be the long_correlates dataset with the appropriate antibodies and log antibodies
  long_correlates <- long_correlates_S
  long_correlates$antibody <- long_correlates$S.PPD
  long_correlates$log_antibody <- log(long_correlates$antibody)
  # Limit of detection
  LOD <- 33
  # Conversion factor from arbitrary units AU/mL in our data to international unit BAU/mL
  conv_factor <- 0.00645
  print(LOD*conv_factor)
  antibody_name <- "anti-spike IgG"
  Antibody_name <- "Anti-spike IgG"
  antibody_units <- "BAU/mL"
  # Axis limits for plots of log antibodies
  # Hard-coding this in so can more easily compare models.
  log_antibody_lims <- c(-0.092,14.8)
  # log_antibody_lims <- c(min(log_A_mat_extreme,long_correlates$log_antibody[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1")]),max(log_A_mat_extreme,long_correlates$log_antibody[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1")]))
  # Marking outliers
  long_correlates$outlier <- F
  # Exclude the exceptionally small results as they appear to be outliers and this appears to be preventing good mixing in the model (rho = -1)
  long_correlates$outlier[which(long_correlates$log_antibody<7 & long_correlates$As_vaccinated_arm_2=="ChAdOx1")] <- T
  # Exclude the exceptionally large result at PB90
  long_correlates$outlier[which((long_correlates$antibody_time>50 &
                                   long_correlates$antibody_time<130 &
                                   long_correlates$log_antibody>13) &
                                  long_correlates$As_vaccinated_arm_2=="ChAdOx1")] <- T
} else if(antibody_type=="neuts"){
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
  # Axis limits for plots of log antibodies
  # Hard-coding this in so can more easily compare models.
  log_antibody_lims <- c(-1,9.2)
  # log_antibody_lims <- c(min(log_A_mat_extreme,long_correlates$log_antibody[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1")]),max(log_A_mat_extreme,long_correlates$log_antibody[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1")]))
  # Marking outliers
  long_correlates$outlier <- F # No outliers
}
## Left-censoring below the limit of detection
nantibody_LOD <- sum(long_correlates$log_antibody==log(LOD))
nantibody_leq_LOD <- sum(long_correlates$log_antibody<=log(LOD))
if(nantibody_LOD<nantibody_leq_LOD & nantibody_LOD>0){suppressWarnings(paste0("There are ",nantibody_LOD," antibody observations equal to the limit of detection. We will treat these as though they are not censored, but observed exactly."))}
# Individuals who were below the limit of detection were given a value equal to the LOD, so these should be treated as censored
if(nantibody_LOD==nantibody_leq_LOD & nantibody_LOD>0){warning(paste0("There are ",nantibody_LOD," antibody observations equal to the limit of detection. We will treat these as though they are not censored, but observed exactly."))}
long_correlates$left_censored <- FALSE
long_correlates$left_censored[which(long_correlates$antibody<LOD)] <- TRUE
long_correlates$antibody[which(long_correlates$left_censored)] <- LOD
long_correlates$log_antibody[which(long_correlates$left_censored)] <- log(LOD)
long_correlates <- long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),] 

# Defining variables
long_correlates <- long_correlates[,c("redcap_data_access_group","As_vaccinated_arm_2","sc_repeat_pid","record_id","sc_age","age_group","cor2dose_primary","cor2dose_positive","antibody","log_antibody","antibody_time","ll","timepoint","visit2","left_censored","outlier")]
# Log antibody levels (untransformed units)
# joint_correlates <- joint_correlates[,c("redcap_data_access_group","As_vaccinated_arm_2","sc_repeat_pid","sc_age","age_group","cor2dose_schedule","cor2dose_interval","cor2dose_non_white","cor2dose_bmi_geq_30","cor2dose_comorbidities","cor2dose_hcw_status","cor2dose_primary","cor2dose_positive","start_time","end_time")]

# The data for vaccinated individuals only
joint_correlates_ChAd <- joint_correlates[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1"),]

# times
amax_time <- 209-7*(long_type=="7_day_incubator") # Minus 7 days if we're using the 7_day_incubator method
amin_time <- 28
adsb <- amin_time:amax_time # Days since boost (ds_t0+21)
if(long_type=="7_day_incubator") {
  raw_data_t0 <- 21 # start_time in long_correlates data
  VEdsb <- adsb + 7 # Days since boost (ds_t0+21)
} else if(long_type=="0_day_incubator") {raw_data_t0 <- 14; VEdsb <- adsb} # Days since boost (ds_t0+21)
t0 <- raw_data_t0 + data_long$t_center;
# Number of days from boost to day 0 in the model (i.e. 28)
ds_t0 <- adsb-t0 # Days since day 0 in the model for the antibodies
nts <- length(ds_t0)
time_units <- c("Days","Weeks","Months")[1]
# div is the Amount to divide days by to get the time units we're using (1 if days, 7 if weeks, 365/12 if months)
if (time_units=="Days"){div <- 1}
if (time_units=="Weeks"){div <- 7}
if (time_units=="Months"){div <- 365/12}
atsb <- adsb/div # antibody time since boost in the appropriate units
VEtsb <- VEdsb/div # vaccine efficacy time since boost in the appropriate units (7 days later than antibody time)

#####
# antibody plots
# log_antibody_lims <- c(min(log_A_mat_extreme,long_correlates$log_antibody[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1")]),max(log_A_mat_extreme,long_correlates$log_antibody[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1")]))
long_correlates$antibody_time_since_boost <- (long_correlates$antibody_time+t0-data_long$t_center)/div

time_axis_lab <- paste0(time_units," since second dose vaccination")
antibody_lab <- paste0(Antibody_name," levels (",antibody_units,")")
time_axis_marks <- c(0,30,60,90,120,150,180,210) # Choose the marks for the x axis in weeks
long_time_axis_marks <- c(time_axis_marks,240,270,300,330,360)
# Place the marks for the y-axis in powers of 10 BAU/mL converted to the log AU/mL antibody scale
antibody_ax_mark_loc <- log(10^seq(floor(log10(exp(log_antibody_lims)*conv_factor))[1],ceiling(log10(exp(log_antibody_lims)*conv_factor))[2],by=1)/conv_factor)
# Write marks for the y-axis in powers of 10 BAU/mL
antibody_ax_mark <- 10^seq(floor(log10(exp(log_antibody_lims)*conv_factor))[1],ceiling(log10(exp(log_antibody_lims)*conv_factor))[2],by=1)
VE_lab <- "Vaccine efficacy (%)"
VE_lims <- c(-50,100)
VE_axis_marks <- c(seq(-60,100,20))
VE_axis_ticks <- c(seq(-60,100,10))
# Define colours
blue_for_histogram <- rgb(0.583,0.830,0.956,alpha=0.6)
green_for_histogram <- rgb(0.875,0.933,0.600,alpha=0.6)

################################################################################
# Calculate parameter effects - effect_summary
# Adjust the antibody effect
if(antibody_name=="anti-spike IgG"){
  antibody_multiplier_WHO <- 100
} else if(antibody_name=="nAb"){antibody_multiplier_WHO <- 10}


###############
# Plot effects of covariates on hazard of infection

effect_summary <- rbind(data.frame(as.data.frame(out$symptomatic$effect_summary),"Outcome" = "symptomatic"),
                        data.frame(as.data.frame(out$allpositive$effect_summary),"Outcome" = "allpositive"))
effect_summary$Outcome <- factor(effect_summary$Outcome, levels = c("symptomatic","allpositive"))

dim(effect_summary)
xmin <- 1/20
xmax <- 20

outcome_cols <- palette("Okabe-Ito")[c(2,7)]
outcome_cols <- colorRampPalette(brewer.pal(5,"Spectral"))(5)
outcome_cols <- outcome_cols[c(1,length(outcome_cols))]
names(outcome_cols) <- names(out)


xlabs <- 1:xmax
# Supplementary Figure 7
png(paste0(plot_directory,"/VE_plots/",model_name,"/","Covariate_effects_hazard.png"),width = 1400,height=1000, pointsize=23)
par(mar = c(6, 4, 6, 2) + 0.1,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
ggplot(effect_summary, aes(x=median, y=rownames, col=Outcome)) +
  geom_point(position=position_dodge(width = 0.6), cex=3) +
  geom_errorbarh(aes(xmin=lower95CI, xmax=upper95CI), height=.4, linewidth=0.7,  position=position_dodge(width = 0.6)) +
  geom_vline(xintercept=1, linetype=2, lwd=1) +
  scale_x_log10(limits = c(xmin,xmax),n.breaks=8) +
  scale_color_manual(values=outcome_cols)+
  ylab(NULL) + xlab(paste0("Multiplicative effect on hazard for COVID-19")) +
  theme_bw() + guides(color=F) +
  theme(text=element_text(size=23),axis.text.x=element_text(size=23, margin = margin(t=0,r=0,b=15,l=0)),
        axis.text.y=element_text(size=20),
        axis.ticks.length.x = unit(0.3,"cm"))
dev.off()

# Create table of the effects of covariates on hazard of infection
ncov <- nrow(effect_summary[which(effect_summary$Outcome=="symptomatic"),])
effect_summary[,c("rownames")] <- as.character(effect_summary[,c("rownames")])
effect_summary_tab <- rbind(c("Multiplicative effect on hazard for symptomatic COVID-19 infection",0,0,0,0),
                            effect_summary[which(effect_summary$Outcome=="symptomatic"),c("rownames","median","mean","lower95CI", "upper95CI")],
                            c("Multiplicative effect on hazard for any COVID-19 infection",0,0,0,0),
                            effect_summary[which(effect_summary$Outcome=="allpositive"),c("rownames","median","mean","lower95CI", "upper95CI")])
effect_summary_tab[c(1,ncov+2),c("median","mean","lower95CI", "upper95CI")] <- rep(c("Median","Mean","0.025 quantile","0.975 quantile"),each=2)
colnames(effect_summary_tab) <- 
  c("Covariate multiplicative effect",
    c("Posterior estimates","","",""))
# Supplementary Table 6
write.table(effect_summary_tab,paste0(directory_correlates,"/4_Output/Tables/Covariate_effects_hazard.csv"),row.names=F,sep=",")


###########################
# VE vs antibody plot

# Figure 2
# VE_summary
BAU_antibodies <- as.numeric(rownames(out[[1]]$VE_summary))

png(paste0(plot_directory,"/VE_plots/",model_name,"/","Mean_VE_vs_antibody.png"),width = 1400,height=1000, pointsize = 23)
par(mar=c(5.1,4.1,1.1,1.1))
plot(out[[1]]$VE_summary[,"lower95CI"]~BAU_antibodies, lty=2, log="x", col="NA",
     type="l", xlab = antibody_lab, yaxt="n", ylim=VE_lims, xaxt="n",
     ylab = VE_lab)
abline(v = 10^seq(-1,4,0.5),col="grey95")
abline(h = seq(-60,100,10),col="grey95")
abline(v = 10^seq(-1,4,1),col="grey90")
abline(h = seq(-60,100,20),col="grey90")
axis(2,at=VE_axis_marks, lab = paste0(VE_axis_marks,"%"), las=1)
log10_antibodies <- log10(BAU_antibodies)
min_anti <- floor(min(log10_antibodies))
max_anti <- ceiling(max(log10_antibodies))
axis(1,at=10^seq(0,max_anti,1),lab=gsub(" ","",format(10^seq(0,max_anti,1),scientific=F),fixed=T))

# Add a histogram of predicted antibody values at end of study day 209 and of observed antibody values
plot(pred_antibody_values_182_hist,add=TRUE,border=green_for_histogram,col=green_for_histogram)
plot(pred_antibody_values_28_hist,add=TRUE,border=green_for_histogram,col=blue_for_histogram)
for (i in 1:n.outcomes){
  polygon(c(BAU_antibodies,rev(BAU_antibodies)), c(out[[i]]$VE_summary[,"upper95CI"], rev(out[[i]]$VE_summary[,"lower95CI"])),
          col = alpha(outcome_cols[i],0.3),lty=0)
  lines(out[[i]]$VE_summary[,"median"]~BAU_antibodies, col = outcome_cols[i], lwd=2)
  lines(out[[i]]$VE_summary[,"lower95CI"]~BAU_antibodies,lty=2, col = outcome_cols[i], lwd=1.5)
  lines(out[[i]]$VE_summary[,"upper95CI"]~BAU_antibodies,lty=2, col = outcome_cols[i], lwd=1.5)
}
abline(v=LOD*conv_factor,lty=2,lwd=2)
abline(h=0,lty=2)
dev.off()

# Create csv of VE_vs_antibodies plot
VE_vs_antibodies <- cbind(out$symptomatic$VE_summary,out$allpositive$VE_summary)
colnames(VE_vs_antibodies)[1:4] <- paste("symptomatic","posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
colnames(VE_vs_antibodies)[5:8] <- paste("anyinfection","posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
write.csv(cbind(BAU_antibodies,VE_vs_antibodies),paste0(directory_correlates,"/4_Output/Plot_data/VE_vs_antibodies.csv"),row.names = F)

# Values
for (i in 1:n.outcomes){
  print(outcome_type[i])
  print(tens_BAU2)
  print(print_CIs(x=round(out[[i]]$VE_summary[as.character(tens_BAU2),],1),add_char="%"))
}

# Calculate the antibody associated with x% VE
VE <- rownames(out$allpositive$antibody_at_VE_summary)
VEs_want <- seq(0.5,0.9,0.2)
# Across the chains
for (i in 1:n.outcomes){
  print(outcome_type[i])
  print(round(aperm(out[[i]]$antibody_at_VE_summary_chains[,as.character(VEs_want),],c(1,3,2)),0))
  print(round(out[[i]]$antibody_at_VE_summary[as.character(VEs_want),"median"],0))
  print(round(out[[i]]$antibody_at_VE_summary[as.character(VEs_want),c("lower95CI","upper95CI")],0))
}

for (i in 1:n.outcomes){
  print(outcome_type[i])
  print(print_CIs(round(out[[i]]$antibody_at_VE_summary[as.character(VEs_want),])))
}
# Create table of antibody levels required for different vaccine efficacies
estimate_interval <- function(x){paste0(x[1]," (",x[2],", ",x[3],")")}
antibody_at_VE_tab <- cbind(round(as.numeric(VE)*100,0),
                            apply(round(out$symptomatic$antibody_at_VE_summary,0),1,estimate_interval),
                            apply(round(out$allpositive$antibody_at_VE_summary,0),1,estimate_interval))
antibody_at_VE_tab <- gsub("NaN","<0",antibody_at_VE_tab)
antibody_at_VE_tab <- rbind(c("", "Required anti-spike IgG antibody level",""),
                            c("Vaccine efficacy (%)", "Symptomatic COVID-19","Any COVID-19 infection"),
                            antibody_at_VE_tab)
# Supplementary Table 7
write.table(antibody_at_VE_tab,paste0(directory_correlates,"/4_Output/Tables/Antibody_at_VE.csv"),row.names=F,col.names=F,sep=",")

# Value at 0 antibodies/direct effect
for (i in 1:n.outcomes){
  print(outcome_type[i])
  print(paste0(round(out[[i]]$direct_effect_CI["median"],1),"% ",paste0("(",paste(round(out[[i]]$direct_effect_CI[c("lower95CI","upper95CI")],1),collapse=",~"),")")))
}

# Antibody indirect effect
for (i in 1:n.outcomes){
  print(outcome_type[i])
  print(print_CIs(x=round(out[[i]]$VE_indirect_summary[as.character(tens_BAU),],0),add_char="%"))
}

# # Plot
# png(paste0(plot_directory,"/VE_plots/",model_name,"/","Mean_VE_indirect_antibody_vs_antibody.png"),width = 1400,height=1000,pointsize = 23)
# par(mar=c(5.1,4.1,1.1,1.1))
# plot(out[[1]]$VE_indirect_summary[,"lower95CI"]~BAU_antibodies, lty=2, log="x",col="NA",
#      type="l", xlab = antibody_lab, yaxt="n", ylim=c(0,100), xaxt="n", yaxs="i",
#      ylab = paste0("Vaccine efficacy due to indirect ",antibody_name," effect (%)"))
# abline(v = 10^seq(-1,4,0.5),col="grey95")
# abline(h = seq(-60,100,10),col="grey95")
# abline(v = 10^seq(-1,4,1),col="grey90")
# abline(h = seq(-60,100,20),col="grey90")
# axis(2,at=VE_axis_marks, lab = paste0(VE_axis_marks,"%"), las=1)
# log10_antibodies <- log10(BAU_antibodies)
# min_anti <- floor(min(log10_antibodies))
# max_anti <- ceiling(max(log10_antibodies))
# axis(1,at=10^seq(0,max_anti,1),lab=gsub(" ","",format(10^seq(0,max_anti,1),scientific=F),fixed=T))
# 
# # Add a histogram of predicted antibody values at end of study day 209 and of observed antibody values
# plot(pred_antibody_values_182_hist,add=TRUE,border=green_for_histogram,col=green_for_histogram)
# plot(pred_antibody_values_28_hist,add=TRUE,border=green_for_histogram,col=blue_for_histogram)
# for (i in 1:n.outcomes){
#   polygon(c(BAU_antibodies,rev(BAU_antibodies)), c(out[[i]]$VE_indirect_summary[,"upper95CI"], rev(out[[i]]$VE_indirect_summary[,"lower95CI"])),
#           col = alpha(outcome_cols[i],0.3),lty=0)
#   lines(out[[i]]$VE_indirect_summary[,"median"]~BAU_antibodies, col = outcome_cols[i], lwd=2)
#   lines(out[[i]]$VE_indirect_summary[,"lower95CI"]~BAU_antibodies,lty=2, col = outcome_cols[i], lwd=1.5)
#   lines(out[[i]]$VE_indirect_summary[,"upper95CI"]~BAU_antibodies,lty=2, col = outcome_cols[i], lwd=1.5)
# }
# abline(h=0,lty=2)
# abline(h=0,lty=2)
# dev.off()

# Create csv of VE_vs_antibodies plot
VE_vs_antibodies_indirect <- cbind(out$symptomatic$VE_indirect_summary,out$allpositive$VE_indirect_summary)
colnames(VE_vs_antibodies_indirect)[1:4] <- paste("symptomatic","posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
colnames(VE_vs_antibodies_indirect)[5:8] <- paste("anyinfection","posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
write.csv(cbind(BAU_antibodies,VE_vs_antibodies_indirect),paste0(directory_correlates,"/4_Output/Plot_data/VE_vs_antibodies_indirect.csv"),row.names = F)

#####
# Figure 3
# Plot the mean VE over time
png(paste0(plot_directory,"/VE_plots/",model_name,"/","Mean_VE_vs_time.png"),width = 1400, height=1000, pointsize = 23)
par(mar=c(5.1,4.1,1.1,1.1))
plot(NA,type="l",ylim=c(0,100),
     ylab= "Mean vaccine efficacy (%)", xlab= time_axis_lab,col=NA,yaxt="n",xlim=c(0,max(VEtsb)), xaxt="n")
abline(v = seq(0,210,15),col="grey95")
abline(h = seq(-60,100,10),col="grey95")
abline(v = seq(0,210,30),col="grey90")
abline(h = seq(-60,100,20),col="grey90")
for (i in 1:n.outcomes){
  polygon(c(VEtsb,rev(VEtsb)), c(out[[i]]$VE_vs_time[,"mean_VE","97.5%"], rev(out[[i]]$VE_vs_time[,"mean_VE","2.5%"])),
          col = alpha(outcome_cols[i],0.3),lty=0)
  lines(out[[i]]$VE_vs_time[,"mean_VE","50%"]~VEtsb, col = outcome_cols[i], lwd=2)
  lines(out[[i]]$VE_vs_time[,"mean_VE","2.5%"]~VEtsb,lty=2, col = outcome_cols[i], lwd=1.5)
  lines(out[[i]]$VE_vs_time[,"mean_VE","97.5%"]~VEtsb,lty=2, col = outcome_cols[i], lwd=1.5)
}
axis(2,at=VE_axis_marks, lab = paste0(VE_axis_marks,"%"), las=1)
axis(2,at=VE_axis_ticks,labels=F)
axis(1,at=time_axis_marks)
abline(h=0,lty=2)
dev.off()

# # And with histograms
# png(paste0(plot_directory,"/VE_plots/",model_name,"/","Mean_VE_vs_time_with_hists.png"), width = 1400, height=1000, pointsize = 23)
# par(mar=c(5.1,4.1,1.1,1.1))
# plot(NA,type="l",ylim=c(0,100), yaxs="i",
#      ylab= "Mean vaccine efficacy (%)", xlab= time_axis_lab,col=NA,yaxt="n",xlim=c(0,max(VEtsb)), xaxt="n")
# abline(v = seq(0,210,15),col="grey95")
# abline(h = seq(-60,100,10),col="grey95")
# abline(v = seq(0,210,30),col="grey90")
# abline(h = seq(-60,100,20),col="grey90")
# axis(2,at=VE_axis_marks, lab = paste0(VE_axis_marks,"%"), las=1)
# axis(2,at=VE_axis_ticks,labels=F)
# axis(1,at=time_axis_marks)
# event_timing_hist <- hist(joint_correlates[which(joint_correlates$cor2dose_positive=="Positive"),c("end_time")]+raw_data_t0, breaks = seq(0,365,by=7), plot=F)
# maximum_height <- 40
# event_timing_hist$counts <- event_timing_hist$density/max(event_timing_hist$density)*maximum_height
# # event_timing_hist$counts <- event_timing_hist$counts*height_factor
# plot(event_timing_hist,add=TRUE,border=blue_for_histogram,col=blue_for_histogram)
# antibody_timing_hist <- hist(long_correlates[,c("antibody_time")]+raw_data_t0, breaks = seq(0,365,by=7), plot=F)
# antibody_timing_hist$counts <- antibody_timing_hist$density/max(antibody_timing_hist$density)*maximum_height*2
# # antibody_timing_hist$counts <- antibody_timing_hist$counts*height_factor
# plot(antibody_timing_hist,add=TRUE,border=blue_for_histogram,col=green_for_histogram)
# for (i in 1:n.outcomes){
#   polygon(c(VEtsb,rev(VEtsb)), c(out[[i]]$VE_vs_time[,"mean_VE","97.5%"], rev(out[[i]]$VE_vs_time[,"mean_VE","2.5%"])),
#           col = alpha(outcome_cols[i],0.3),lty=0)
#   lines(out[[i]]$VE_vs_time[,"mean_VE","50%"]~VEtsb, col = outcome_cols[i], lwd=2)
#   lines(out[[i]]$VE_vs_time[,"mean_VE","2.5%"]~VEtsb,lty=2, col = outcome_cols[i], lwd=1.5)
#   lines(out[[i]]$VE_vs_time[,"mean_VE","97.5%"]~VEtsb,lty=2, col = outcome_cols[i], lwd=1.5)
# }
# axis(2,at=VE_axis_marks, lab = paste0(VE_axis_marks,"%"), las=1)
# axis(2,at=VE_axis_ticks,labels=F)
# axis(1,at=time_axis_marks)
# dev.off()

# Plot of estimates of baseline hazard
require(survPen)
joint_correlates_Control <- joint_correlates[which(joint_correlates$As_vaccinated_arm_2=="Control"),]
joint_correlates_Control$cor2dose_primary_ind <- as.numeric(levels(joint_correlates_Control$cor2dose_primary_ind))[joint_correlates_Control$cor2dose_primary_ind]
joint_correlates_Control$cor2dose_positive_ind <- as.numeric(levels(joint_correlates_Control$cor2dose_positive_ind))[joint_correlates_Control$cor2dose_positive_ind]
df <- 20
# Primary symptomatic infections
# Estimate baseline hazard observed in the control group using penalized restricted cubic splines
f.prim <- ~ smf(end_cal_time,df=df) # knots placed at quantiles
mod.prim <- survPen(f.prim,data=joint_correlates_Control,t0=start_time,t1=end_cal_time,event=cor2dose_primary_ind)

# Use the output to estimate the cumulative hazard
tdiff <- 0.1
tt <- seq(0,360,tdiff)
pred.prim <- predict(mod.prim,data.frame(end_cal_time=tt))
min(joint_correlates_Control$end_cal_time[which(joint_correlates_Control$cor2dose_primary_ind==1)])
# Nelson-Aalen estimate for cumulative hazard
nelsonaalen_prim <- survfit(Surv(start_time,end_cal_time,cor2dose_primary_ind)~1,data=joint_correlates_Control)
# Plot both to check they give the same result
plot(nelsonaalen_prim, cumhaz=T,conf.int=T)
lines(cumsum(pred.prim$haz[which(tt>=60)]*tdiff) ~ tt[which(tt>=60)], col="blue")

# Any COVID-19 infection
# Estimate baseline hazard observed in the control group using penalized restricted cubic splines
f.pos <- ~ smf(end_cal_time,df=df) # knots placed at quantiles
mod.pos <- survPen(f.pos,data=joint_correlates_Control,t0=start_time,t1=end_cal_time,event=cor2dose_positive_ind)

# Use the output to estimate the cumulative hazard
tdiff <- 0.1
tt <- seq(0,360,tdiff)
pred.pos <- predict(mod.pos,data.frame(end_cal_time=tt), conf.int = 0.95)
min(joint_correlates_Control$end_cal_time[which(joint_correlates_Control$cor2dose_positive_ind==1)])
# Nelson-Aalen estimate for cumulative hazard
nelsonaalen_pos <- survfit(Surv(start_time,end_cal_time,cor2dose_positive_ind)~1,data=joint_correlates_Control)
# Plot both to check they give the same result
plot(nelsonaalen_pos, cumhaz=T,conf.int=T)
lines(cumsum(pred.pos$haz[which(tt>=60)]*tdiff) ~ tt[which(tt>=60)], col="blue")
lines(cumsum(pred.pos$haz.inf[which(tt>=60)]*tdiff) ~ tt[which(tt>=60)], col="blue",lty=2)
lines(cumsum(pred.pos$haz.sup[which(tt>=60)]*tdiff) ~ tt[which(tt>=60)], col="blue",lty=2)
difftime(as.Date(joint_correlates$cor2dose_end_date[1]),as.Date("2020-07-18"))

# Supplementary Figure 3
# Plot the hazard over time
png(paste0(plot_directory,"/VE_plots/",model_name,"/","Baseline_hazard.png"), width = 1400, height=1000, pointsize = 23)
par(mar=c(5.1,4.1,1.1,1.1), mfrow=c(1,1))
plot(NA,xlim=c(0,250),ylim=c(0,0.002), ylab = "Baseline hazard", xlab = "Calendar days since first participant at risk on 18 July 2020")
abline(v = seq(0,250,25),col="grey95")
abline(h = seq(0,0.002,0.00025),col="grey95")
abline(v = seq(0,250,50),col="grey90")
abline(h = seq(0,0.002,0.0005),col="grey90")
polygon(c(tt,rev(tt)), c(pred.pos$haz.sup, rev(pred.pos$haz.inf)),
        col = alpha(outcome_cols[which(names(outcome_cols)=="allpositive")],0.3),lty=0)
polygon(c(tt,rev(tt)), c(pred.prim$haz.sup, rev(pred.prim$haz.inf)),
        col = alpha(outcome_cols[which(names(outcome_cols)=="symptomatic")],0.3),lty=0)
lines(pred.pos$haz ~ tt, col = outcome_cols[which(names(outcome_cols)=="allpositive")],lwd=3)
lines(pred.prim$haz ~ tt, col = outcome_cols[which(names(outcome_cols)=="symptomatic")],lwd=3)
dev.off()

# Supplementary Figure 5
# Plot histogram of timings of antibody measurements over time since second dose
png(paste0(plot_directory,"/VE_plots/",model_name,"/","Antibody_time_hist.png"), width = 1400, height=1000, pointsize = 23)
par(mar=c(5.1,4.1,1.1,1.1))
antibody_time_hist <- hist(long_correlates[,c("antibody_time")]+raw_data_t0, breaks = seq(0,365,by=7),plot=F)
plot(NA, xlab = time_axis_lab, yaxs="i",
     ylab = "Number of antibody measurements", main = NA, xaxt="n", xlim = c(0,360), ylim = c(0,850))
abline(v = seq(0,360,15),col="grey95")
abline(h = seq(0,800,100),col="grey95")
abline(v = seq(0,360,30),col="grey90")
abline(h = seq(0,800,200),col="grey90")
plot(antibody_time_hist,add=TRUE, col = blue_for_histogram)
axis(1,at=long_time_axis_marks)
dev.off()

# Plot histogram of number of people at risk over time since second dose
times <- seq(0,max(long_time_axis_marks))
# Calculates the number at risk at time "time" when people exit study at "end_times"
# and enter at "start_times" - taken to be everyone enters at 0 if left null.
n_at_risk_fun <- function(time,end_times,start_times=NULL){
  if (!is.null(start_times)){
    return(sum(end_times>=time && start_times < time))
  }
  else{
    return(sum(end_times>=time))
  }
}
n_at_risk <- sapply(times,n_at_risk_fun, end_times = joint_correlates[,c("end_time")]+raw_data_t0)
n_at_risk_ChAd <- sapply(times,n_at_risk_fun, end_times = joint_correlates_ChAd[,c("end_time")]+raw_data_t0)
# Supplementary Figure 4
png(paste0(plot_directory,"/VE_plots/",model_name,"/","Numbers_at_risk_over_time.png"), width = 1400, height=1000, pointsize = 23)
par(mar=c(5.1,4.1,2.1,1.1))
plot(NA, type = "l",xlab = time_axis_lab, ylab = "Number of participants at risk", main = NA, xaxt="n", yaxs="i",
     ylim= c(0,10000), xlim= c(0,max(long_time_axis_marks)))
abline(v = seq(0,360,15),col="grey95")
abline(h = seq(0,10000,1000),col="grey95")
abline(v = seq(0,360,30),col="grey90")
abline(h = seq(0,10000,2000),col="grey90")
polygon(c(times,rev(times)), c(n_at_risk, rev(n_at_risk_ChAd)),
        col = green_for_histogram,lty=0)
polygon(c(times,rev(times)), c(n_at_risk_ChAd, rev(rep(0,length(times)))),
        col = blue_for_histogram,lty=0)
lines(n_at_risk)
lines(n_at_risk_ChAd)
axis(1,at=long_time_axis_marks, xaxs="i")
legend(x=270,y=4000,legend = c("Control", "ChAdOx1 nCoV-19"),
       col = c(green_for_histogram,blue_for_histogram), lwd=15)
dev.off()

# Supplementary Figure 6
# Plot histogram of timings of cases over time since second dose
png(paste0(plot_directory,"/VE_plots/",model_name,"/","Cases_hist.png"), width = 1400, height=800, pointsize = 23)
par(mar=c(5.1,4.1,3.1,1.1),mfrow=c(1,2))
# Positive cases
cases_positive_hist <- hist(joint_correlates[which(joint_correlates$cor2dose_positive=="Positive"),c("end_time")]+raw_data_t0, breaks = seq(0,365,by=7),plot=F)
plot(NA, xlab = time_axis_lab, yaxs="i", main = "(a) All COVID-19 cases",
     ylab = "Number of cases", xaxt="n", xlim = c(0,360), ylim = c(0,50))
abline(v = seq(0,360,15),col="grey95")
abline(h = seq(0,50,5),col="grey95")
abline(v = seq(0,360,30),col="grey90")
abline(h = seq(0,50,10),col="grey90")
plot(cases_positive_hist,add=TRUE, col = alpha(outcome_cols[which(names(outcome_cols)=="allpositive")],0.5))
axis(1,at=long_time_axis_marks)
# Primary cases
cases_primary_hist <- hist(joint_correlates[which(joint_correlates$cor2dose_primary=="Case"),c("end_time")]+raw_data_t0, breaks = seq(0,365,by=7),plot=F)
plot(NA, xlab = time_axis_lab, yaxs="i", main = "(b) Primary symptomatic COVID-19 cases",
     ylab = "Number of primary symptomatic cases", xaxt="n", xlim = c(0,360), ylim = c(0,50))
abline(v = seq(0,360,15),col="grey95")
abline(h = seq(0,50,5),col="grey95")
abline(v = seq(0,360,30),col="grey90")
abline(h = seq(0,50,10),col="grey90")
plot(cases_primary_hist,add=TRUE, col = alpha(outcome_cols[which(names(outcome_cols)=="symptomatic")],0.5))
axis(1,at=long_time_axis_marks)
dev.off()

# Create csv of VE_vs_time plot
{VE_vs_time <- cbind(out$symptomatic$VE_vs_time[,"mean_VE",],out$allpositive$VE_vs_time[,"mean_VE",])
colnames(VE_vs_time)[1:4] <- paste("symptomatic","posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
colnames(VE_vs_time)[5:8] <- paste("anyinfection","posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
VE_vs_time <- cbind(VEtsb,VE_vs_time); colnames(VE_vs_time)[1] <- "Days since second dose"
write.csv(VE_vs_time,paste0(directory_correlates,"/4_Output/Plot_data/VE_vs_time.csv"), row.names = F)}

# Print mean VE at times
times <- as.character(c(28,90,182)+7)
for (i in 1:n.outcomes){
  print(outcome_type[i])
  VE_vs_time_summary <- round(out[[i]]$VE_vs_time[times,"mean_VE",],1)
  print(times)
  print(paste0(paste0(paste0(VE_vs_time_summary[,"50%"],"% "),paste0("(",apply(VE_vs_time_summary[,c("2.5%","97.5%")],1,paste,collapse=",~"),")")),collapse = ", "))
}

# Print VE quantiles at times
quants_want <- c(0.25,0.5,0.75)
nquants_want <- length(quants_want)
print("Antibody quantiles")
antibody_quant <- paste("antibody_quantiles",quants_want,sep="_")
times_ant <- as.character(as.numeric(times)-7)
for (s in 1:length(times)){
  antibody_quants_times <- round(antibody_quantiles[times_ant[s],antibody_quant,],0)
  print(paste0("Day ",times[s],": ",paste0(paste0(paste0(antibody_quants_times[,"50%"]," "),paste0("(",apply(antibody_quants_times[,c("2.5%","97.5%")],1,paste,collapse=",~"),")")),collapse = ", ")))
}
for (i in 1:n.outcomes){
  print(outcome_type[i])
  print((quants_want))
  for (s in 1:length(times)){
    VE_quants_times <- round(out[[i]]$VE_vs_time[times[s],as.character(quants_want),],1)
    print(paste0("Day ",times[s],": ",paste0(paste0(paste0(VE_quants_times[,"50%"],"% "),paste0("(",apply(VE_quants_times[,c("2.5%","97.5%")],1,paste,collapse=",~"),")")),collapse = ", ")))
  }
}

# Use more quantiles for plots
quants_want <- c(0.05,0.25,0.5,0.75,0.95)
nquants_want <- length(quants_want)

# Quantiles of VE_vs_time
nhalf_quants <- (nquants_want-1)/2
quant_cols <- brewer.pal(9,"Spectral")[c(rev(1:nhalf_quants),nhalf_quants+1,10-(1:nhalf_quants))]
quant_cols <- quant_cols[(length(quant_cols):1)]
quant_cols[nhalf_quants+1] <- "black"

ab <- c("(a)","(b)")

# Supplementary Figure 8
png(paste0(plot_directory,"/VE_plots/",model_name,"/","VE_at_antibody_quantiles_vs_time.png"),width = 1700,height=1000, pointsize = 23)
par(mfrow=c(1,2), mar = c(4.5,3.5,3.1,0.5),oma=c(0,1.5,0,0))
for (i in 1:n.outcomes){
  outcome <- outcome_type[i]
  plot(NA,ylim=VE_lims, xlim=c(0,max(VEtsb)), xaxt="n", yaxt="n",ylab=NA,
       xlab = time_axis_lab, main = paste(ab[i],outcome_title[i]))
  abline(v = seq(0,210,15),col="grey95")
  abline(h = seq(-60,100,10),col="grey95")
  abline(v = seq(0,210,30),col="grey90")
  abline(h = seq(-60,100,20),col="grey90")
  for (k in 1:nquants_want){
    quant <- as.character(quants_want[k])
    polygon(c(VEtsb,rev(VEtsb)), c(out[[i]]$VE_vs_time[,quant,"97.5%"], rev(out[[i]]$VE_vs_time[,quant,"2.5%"])),
            col = alpha(quant_cols[k],0.3),lty=0)
  }
  for (k in 1:nquants_want){
    quant <- as.character(quants_want[k])
    lines(out[[i]]$VE_vs_time[,quant,"50%"]~VEtsb, col = quant_cols[k],lwd=3)
    lines(out[[i]]$VE_vs_time[,quant,"2.5%"]~VEtsb, col=quant_cols[k],lty=2,lwd=1.5)
    lines(out[[i]]$VE_vs_time[,quant,"97.5%"]~VEtsb, col=quant_cols[k],lty=2,lwd=1.5)
  }
  abline(h=0,lty=2,lwd=1.5)
  axis(2,at=VE_axis_marks, lab = paste0(VE_axis_marks,"%"), las=1)
  axis(1,at=time_axis_marks)
  axis(2,at=VE_axis_ticks,labels=F)
  legend(x=0,y=-1,legend=paste0(quants_want*100,"%")[nquants_want:1],col=quant_cols[nquants_want:1],lwd=3, title = "Quantiles")
}
mtext(VE_lab,side=2, outer=T, line=0)
dev.off()

{
# Quantiles of antibody levels over time
# Plot antibody quantiles over time
ab_marks <- seq(0,1000,100)
# Supplementary Figure 9
png(paste0(plot_directory,"/Antibody_plots/","/","Antibody_quantiles_vs_time.png"),width = 1400,height=1000, pointsize = 23)
par(mar=c(5.1,4.1,2.1,2.1))
plot(NA,ylim=c(0,min(ab_marks[which(ab_marks>max(antibody_quantiles))])), xlim=c(0,max(atsb)),
     ylab = antibody_lab, xaxt="n", yaxs="i", yaxt="n",
     xlab = time_axis_lab)
abline(v = seq(0,210,15),col="grey95")
abline(h = seq(0,1000,50),col="grey95")
abline(v = seq(0,210,30),col="grey90")
abline(h = seq(0,1000,100),col="grey90")
for (k in 1:nquants_want){
  antibody_quant <- paste("antibody_quantiles",quants_want[k],sep="_")
  polygon(c(atsb,rev(atsb)), c(antibody_quantiles[,antibody_quant,"97.5%"], rev(antibody_quantiles[,antibody_quant,"2.5%"])),
          col = alpha(quant_cols[k],0.3),lty=0)
  lines(antibody_quantiles[,antibody_quant,"50%"]~atsb, col=quant_cols[k], lwd=3)
  lines(antibody_quantiles[,antibody_quant,"2.5%"]~atsb, col=quant_cols[k],lty=2, lwd=1.5)
  lines(antibody_quantiles[,antibody_quant,"97.5%"]~atsb, col=quant_cols[k],lty=2, lwd=1.5)
}
abline(h=0,lty=1)
axis(1,at=time_axis_marks)
axis(2,seq(0,1000,100))
legend(x=1.5,y=290,legend=paste0(quants_want*100,"%")[nquants_want:1],col=quant_cols[nquants_want:1],lwd=3, title = "Quantiles")
dev.off()

# Plot antibody quantiles over time on log scale
png(paste0(plot_directory,"/Antibody_plots/","/","Antibody_quantiles_log_vs_time.png"),width = 1400,height=1000, pointsize = 23)
par(mar=c(5.1,4.1,2.1,2.1))
plot(NA,ylim=c(min(antibody_quantiles),max(antibody_quantiles)), xlim=c(0,max(atsb)),
     ylab = antibody_lab, xaxt="n", log="y", yaxt="n",
     xlab = time_axis_lab)
abline(v = seq(0,210,15),col="grey95")
abline(h = 10^seq(-3,5,0.5),col="grey95")
abline(v = seq(0,210,30),col="grey90")
abline(h = 10^seq(-3,5),col="grey90")
for (k in 1:nquants_want){
  antibody_quant <- paste("antibody_quantiles",quants_want[k],sep="_")
  polygon(c(atsb,rev(atsb)), c(antibody_quantiles[,antibody_quant,"97.5%"], rev(antibody_quantiles[,antibody_quant,"2.5%"])),
          col = alpha(quant_cols[k],0.3),lty=0)
  lines(antibody_quantiles[,antibody_quant,"50%"]~atsb, col=quant_cols[k], lwd=3)
  lines(antibody_quantiles[,antibody_quant,"2.5%"]~atsb, col=quant_cols[k],lty=2, lwd=1.5)
  lines(antibody_quantiles[,antibody_quant,"97.5%"]~atsb, col=quant_cols[k],lty=2, lwd=1.5)
}
abline(h=0,lty=1)
axis(1,at=time_axis_marks)
axis(2,antibody_ax_mark,antibody_ax_mark)
legend(x=-3,y=60,legend=paste0(quants_want*100,"%")[nquants_want:1],col=quant_cols[nquants_want:1],lwd=3, title = "Quantiles")
dev.off()
}

# Create csv of VE_vs_time quantiles plot
nquants_VEant_vs_time <- dim(out$symptomatic$VE_vs_time[,dimnames(out$symptomatic$VE_vs_time)[[2]]!="mean_VE",])[[2]]
primary_VE_vs_time_quantiles <- out$symptomatic$VE_vs_time[,dimnames(out$symptomatic$VE_vs_time)[[2]]!="mean_VE",]
primary_VE_vs_time_quantiles_mat <- t(apply(aperm(primary_VE_vs_time_quantiles,c(1,3,2)),1,c, use.names=T))
colnames(primary_VE_vs_time_quantiles_mat) <- paste("symptomatic",rep(dimnames(out$symptomatic$VE_vs_time)[[2]][-1],each=4),"posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
positive_VE_vs_time_quantiles <- out$allpositive$VE_vs_time[,dimnames(out$allpositive$VE_vs_time)[[2]]!="mean_VE",]
positive_VE_vs_time_quantiles_mat <- t(apply(aperm(positive_VE_vs_time_quantiles,c(1,3,2)),1,c, use.names=T))
colnames(positive_VE_vs_time_quantiles_mat) <- paste("anyinfection",rep(dimnames(out$allpositive$VE_vs_time)[[2]][-1],each=4),"posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
VE_quantiles_vs_time <- cbind(VEtsb,primary_VE_vs_time_quantiles_mat,positive_VE_vs_time_quantiles_mat)
colnames(VE_quantiles_vs_time)[1] <- "Days since second dose"
VE_quantiles_vs_time <- rbind(VE_quantiles_vs_time,c(""),
                              c("Outcome",rep(c("Symptomatic COVID-19","Any COVID-19 infection"),each=nquants_VEant_vs_time*4)),
                              c("Quantile of vaccine efficacy in study population",rep(rep(dimnames(out$symptomatic$VE_vs_time)[[2]][-1],each=4),2)),
                              c("Posterior summary",rep(c("mean","median","0.025_quantile","0.975_quantile"),nquants_VEant_vs_time*2)))
write.table(VE_quantiles_vs_time,paste0(directory_correlates,"/4_Output/Plot_data/VE_quantiles_vs_time.csv"), row.names = F,sep=",")

# Create csv of Antibody_vs_time quantiles plot
antibody_vs_time_quantiles <- out$symptomatic$antibody_quantiles[,dimnames(out$symptomatic$antibody_quantiles)[[2]]!="mean_VE",]
antibody_vs_time_quantiles_mat <- t(apply(aperm(antibody_vs_time_quantiles,c(1,3,2)),1,c, use.names=T))
colnames(antibody_vs_time_quantiles_mat) <- paste(rep(dimnames(out$symptomatic$antibody_quantiles)[[2]],each=4),"posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
antibody_vs_time_quantiles_mat <- cbind(rownames(antibody_vs_time_quantiles_mat),antibody_vs_time_quantiles_mat)
colnames(antibody_vs_time_quantiles_mat)[1] <- "Days since second dose"
quant_names <- dimnames(out$symptomatic$antibody_quantiles)[[2]]
quant_names <- gsub("antibody_quantiles_","",quant_names)
antibody_vs_time_quantiles_mat <- rbind(antibody_vs_time_quantiles_mat,c(""),
                                        c("Quantile of antibody level in study population",rep(quant_names,each=4)),
                              c("Posterior summary",rep(c("mean","median","0.025_quantile","0.975_quantile"),nquants_VEant_vs_time)))
write.table(antibody_vs_time_quantiles_mat,paste0(directory_correlates,"/4_Output/Plot_data/antibody_quantiles_vs_time.csv"), row.names = F,sep=",")

# Covariate effects on VE
rowlabels_list <- list(c("baseline","56-69",">=70"), c("baseline","Male"),
                       c("baseline","Other"), c("baseline","Comorbidity"),
                       c("baseline","BMI_geq_30"),
                       c("baseline","Less_than_1_covid","More_than_1_covid"),
                       c("baseline","9-11","6-8","<6"),
                       c("baseline","LD"))
do.call(c,rowlabels_list)
rownames_list <- list(c("18-55","56-69",">=70"), c("Female","Male"),
                      c("White","Other"), c("None","Comorbidity"),
                      c("BMI < 30","BMI \u2265 30"),
                      c("Non-HCW","HCW < 1 COVID","HCW \u2265 1 COVID"),
                      c("\u226512","9-11","6-8","<6"),
                      c("Standard dose","Low dose"))
variable_names <- c("Age (years)","Sex","Ethnicity","Comorbidity","BMI (kg/m^2)",
                    "Healthcare worker (HCW) status",
                    "Interval between\n first and second dose (weeks)",
                    "First dose \n(Low dose vs standard dose)")
names(rowlabels_list) <- variable_names
names(rownames_list) <- variable_names
nvar <- length(rowlabels_list)
grid_size <- 4

# Figure 5 and Supplementary Figure 10
cov_cols <- palette("Okabe-Ito")
for (k in 1:n.outcomes){
  png(paste0(plot_directory,"/VE_plots/",model_name,"/","Covariate_effects_mean_VE_vs_time_",outcome_type[k],".png"),width = 1800,height=1000,
      pointsize = 23)
  par(mfrow=c(2,grid_size), mar = c(2.5,2.5,2.5,0.1),oma=c(4.1,4.1,0.1,1))
  count <- 0
  for (i in 1:nvar){
    count <- count + 1
    plot(NA,type="l",ylim=c(0,100),
         main = variable_names[i],
         ylab= "Vaccine efficacy (%)", xlab= time_axis_lab,col=NA,yaxt="n",xlim=c(0,max(VEtsb)), xaxt="n")
    abline(v = seq(0,210,15),col="grey95")
    abline(h = seq(-60,100,10),col="grey95")
    abline(v = seq(0,210,30),col="grey90")
    abline(h = seq(-60,100,20),col="grey90")
    for (j in 1:length(rownames_list[[i]])){
      lines(out[[k]]$VE_vs_time_new_ind[,rowlabels_list[[i]][[j]],"2.5%"]~VEtsb,lty=2, col = cov_cols[j], lwd=2)
      lines(out[[k]]$VE_vs_time_new_ind[,rowlabels_list[[i]][[j]],"50%"]~VEtsb,lty=1, col = cov_cols[j], lwd=2)
      lines(out[[k]]$VE_vs_time_new_ind[,rowlabels_list[[i]][[j]],"97.5%"]~VEtsb,lty=2, col = cov_cols[j], lwd=2)
    }
    axis(1,time_axis_marks)
    legend(x=0,y=30,legend = rownames_list[[i]], col = cov_cols, lty=1, lwd=2)
    if ((((count-1)%%grid_size==0))){
      axis(2,at=VE_axis_marks, lab = paste0(VE_axis_marks,"%"), las=1)
      axis(2,at=VE_axis_ticks,labels=F)
    } else{axis(2,at=VE_axis_ticks,labels=F)}
    abline(h=0,lty=2, lwd=1.5)
  }
  mtext(time_axis_lab,side=1,outer=TRUE,cex=1.3,line=2)
  mtext("Estimated vaccine efficacy (%)",side=2,line=2,outer=TRUE,cex=1.3,las=0)
  dev.off()
}

# Print VE quantiles at times
print("Antibody quantiles")
covar <- dimnames(out$symptomatic$VE_vs_time_new_ind)[[2]][c(1,3,10)]
n.covar <- length(covar)
times_new_ind <- as.character(c(28,90,182)+7)
for (s in 1:length(times_new_ind)){
  print(times_new_ind[s])
  print("\n")
  for (i in 1:n.covar){
  print((covar[i]))
    VE_vs_time_new_ind_times <- round(out$symptomatic$VE_vs_time_new_ind[times_new_ind[s],covar[i],],1)
    print(paste0("Day ",times_new_ind[s],": ",paste0(paste0(paste0(VE_vs_time_new_ind_times["50%"],"% "),paste0("(",paste(VE_vs_time_new_ind_times[c("2.5%","97.5%")],collapse=",~"),")")),collapse = ", ")))
  }
}

# Create csv of VE_vs_time for new individuals (covariate effects) plot
n_new_inds_VE_vs_time <- dim(out$symptomatic$VE_vs_time_new_ind[,dimnames(out$symptomatic$VE_vs_time_new_ind)[[2]]!="mean_VE",])[[2]]
primary_VE_vs_time_new_inds <- out$symptomatic$VE_vs_time_new_ind[,dimnames(out$symptomatic$VE_vs_time_new_ind)[[2]]!="mean_VE",]
primary_VE_vs_time_new_inds_mat <- t(apply(aperm(primary_VE_vs_time_new_inds,c(1,3,2)),1,c, use.names=T))
new_ind_names <- paste0(c("",rep("age_",2),"","ethnicity_","","",rep("interval_",3),rep("healthcare_worker_",2),"first_dose_"),dimnames(out$symptomatic$VE_vs_time_new_ind)[[2]])
colnames(primary_VE_vs_time_new_inds_mat) <- paste("symptomatic",rep(new_ind_names,each=4),"posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
positive_VE_vs_time_new_inds <- out$allpositive$VE_vs_time_new_ind[,dimnames(out$allpositive$VE_vs_time_new_ind)[[2]]!="mean_VE",]
positive_VE_vs_time_new_inds_mat <- t(apply(aperm(positive_VE_vs_time_new_inds,c(1,3,2)),1,c, use.names=T))
colnames(positive_VE_vs_time_new_inds_mat) <- paste("allpositive",rep(new_ind_names,each=4),"posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
VE_new_inds_vs_time <- cbind(VEtsb,primary_VE_vs_time_new_inds_mat,positive_VE_vs_time_new_inds_mat)
colnames(VE_new_inds_vs_time)[1] <- "Days since second dose"
VE_new_inds_vs_time <- rbind(VE_new_inds_vs_time,c(""),
                             c("Outcome",rep(c("Symptomatic COVID-19","Any COVID-19 infection"),each=n_new_inds_VE_vs_time*4)),
                              c("Covariate for new individual",rep(rep(new_ind_names,each=4),2)),
                              c("Posterior summary",rep(c("mean","median","0.025_quantile","0.975_quantile"),n_new_inds_VE_vs_time*2)))
write.table(VE_new_inds_vs_time,paste0(directory_correlates,"/4_Output/Plot_data/VE_covariate_effects_vs_time.csv"), row.names = F,sep=",")


quant_cols_omicron <- colorRampPalette(brewer.pal(5,"Spectral"))(nquants_want)
quant_cols_omicron <- c("black",quant_cols_omicron[c(1,length(quant_cols_omicron))])
effects_omicron <- dimnames(VE_vs_time_omicron)[[2]]

# Supplementary Figure 11
# Extrapolation of VE vs time to Omicron variant
png(paste0(plot_directory,"/VE_plots/",model_name,"/","Omicron_mean_VE_vs_time.png"),width = 1400,height=1000, pointsize = 23)
par(mar=c(5.1,4.1,2.1,2.1))
plot(VE_vs_time_omicron[,"Omicron_mean_VE_mean_effect","50%"]~VEtsb,type="l",ylim=c(0,100),
     ylab= "Predicted mean relative efficacy against Omicron BA.4/5 (%)", xlab= time_axis_lab,col=NA,yaxt="n",xlim=c(0,max(VEtsb)), xaxt="n")
abline(v = seq(0,210,15),col="grey95")
abline(h = seq(-60,100,10),col="grey95")
abline(v = seq(0,210,30),col="grey90")
abline(h = seq(-60,100,20),col="grey90")
# Mean VE effect estimate
for (k in 1:length(effects_omicron)){
  polygon(c(VEtsb,rev(VEtsb)), c(VE_vs_time_omicron[,effects_omicron[k],"97.5%"], rev(VE_vs_time_omicron[,effects_omicron[k],"2.5%"])),
          col = alpha(quant_cols_omicron[k],0.3),lty=0, angle=45, lend=2)
  lines(VE_vs_time_omicron[,effects_omicron[k],"50%"]~VEtsb,lty=1, lwd=3, col = quant_cols_omicron[k])
  lines(VE_vs_time_omicron[,effects_omicron[k],"2.5%"]~VEtsb,lty=2, lwd=1, col = quant_cols_omicron[k])
  lines(VE_vs_time_omicron[,effects_omicron[k],"97.5%"]~VEtsb,lty=2, lwd=1, col = quant_cols_omicron[k])
}
axis(2,at=VE_axis_marks, lab = paste0(VE_axis_marks,"%"), las=1)
axis(2,at=VE_axis_ticks,labels=F)
axis(1,at=time_axis_marks)
abline(h=0,lty=2,lwd=2)
dev.off()

dimnames(VE_vs_time_omicron)
effect <- dimnames(VE_vs_time_omicron)[[2]][c(1)]
n.effect <- length(effect)
times_new_ind <- as.character(c(28,90,182)+7)
for (s in 1:length(times_new_ind)){
  print(times_new_ind[s])
  print("\n")
  for (i in 1:n.effect){
    print((effect[i]))
    VE_vs_time_new_ind_times <- round(VE_vs_time_omicron[times_new_ind[s],effect[i],],1)
    print(paste0("Day ",times_new_ind[s],": ",paste0(paste0(paste0(VE_vs_time_new_ind_times["50%"],"% "),paste0("(",paste(VE_vs_time_new_ind_times[c("2.5%","97.5%")],collapse=",~"),")")),collapse = ", ")))
  }
}
dimnames(VE_vs_time_omicron)
effect <- dimnames(VE_vs_time_omicron)[[2]][c(2,3)]
n.effect <- length(effect)
times_new_ind <- as.character(c(28,90,182)+7)
for (s in 1:length(times_new_ind)){
  print(times_new_ind[s])
  print("\n")
  for (i in 1:n.effect){
    print((effect[i]))
    VE_vs_time_new_ind_times <- round(VE_vs_time_omicron[times_new_ind[s],effect[i],],1)
    print(paste0("Day ",times_new_ind[s],": ",paste0(paste0(paste0(VE_vs_time_new_ind_times["50%"],"% "),paste0("(",paste(VE_vs_time_new_ind_times[c("2.5%","97.5%")],collapse=",~"),")")),collapse = ", ")))
  }
}

# Create csv of VE_vs_time for omicron variant plot
n_scenarios_VE_vs_time <- dim(out$symptomatic$VE_vs_time_omicron[,dimnames(out$symptomatic$VE_vs_time_omicron)[[2]]!="mean_VE",])[[2]]
VE_vs_time_omicron_tab <- out$symptomatic$VE_vs_time_omicron[,dimnames(out$symptomatic$VE_vs_time_omicron)[[2]]!="mean_VE",]
VE_vs_time_omicron_mat <- t(apply(aperm(VE_vs_time_omicron_tab,c(1,3,2)),1,c, use.names=T))
scenario_names <- dimnames(out$symptomatic$VE_vs_time_omicron)[[2]]
colnames(VE_vs_time_omicron_mat) <- paste("anyinfection",rep(scenario_names,each=4),"posterior",c("mean","median","0.025_quantile","0.975_quantile"),sep="_")
VE_vs_time_omicron_tab <- cbind(VEtsb,VE_vs_time_omicron_mat)
colnames(VE_vs_time_omicron_tab)[1] <- "Days since second dose"
scenario_names2 <- c("estimate","upper 95% confidence bound","lower 95% confidence bound")
VE_vs_time_omicron_tab <- rbind(VE_vs_time_omicron_tab,c(""),
                             c("Assumed VE vs antibody curve from Wei et al. (2023)",rep(scenario_names2,each=4)),
                             c("Posterior summary",rep(c("mean","median","0.025_quantile","0.975_quantile"),n_scenarios_VE_vs_time)))
write.table(VE_vs_time_omicron_tab,paste0(directory_correlates,"/4_Output/Plot_data/VE_vs_time_omicron.csv"), row.names = F,sep=",")

################################################
# Discussion

######
# Correlates of protection
antibody_levels <- as.numeric(dimnames(out$symptomatic$VE_summary)[[1]])
x <- 21 # which.min(abs(antibody_levels-x)) gives the index of antibody levels whose value is closest to x
antibody_levels[which.min(abs(antibody_levels-x))]
round(out$symptomatic$VE_summary[which.min(abs(antibody_levels-x)),c("median","lower95CI","upper95CI")],1)
x <- 25 # which.min(abs(antibody_levels-x)) gives the index of antibody levels whose value is closest to x
antibody_levels[which.min(abs(antibody_levels-x))]
round(out$symptomatic$VE_summary[which.min(abs(antibody_levels-x)),c("median","lower95CI","upper95CI")],1)
x <- 100 # which.min(abs(antibody_levels-x)) gives the index of antibody levels whose value is closest to x
antibody_levels[which.min(abs(antibody_levels-x))]
round(out$symptomatic$VE_summary[which.min(abs(antibody_levels-x)),c("median","lower95CI","upper95CI")],1)
x <- 400 # which.min(abs(antibody_levels-x)) gives the index of antibody levels whose value is closest to x
antibody_levels[which.min(abs(antibody_levels-x))]
round(out$symptomatic$VE_summary[which.min(abs(antibody_levels-x)),c("median","lower95CI","upper95CI")],1)
x <- 1000 # which.min(abs(antibody_levels-x)) gives the index of antibody levels whose value is closest to x
antibody_levels[which.min(abs(antibody_levels-x))]
round(out$symptomatic$VE_summary[which.min(abs(antibody_levels-x)),c("median","lower95CI","upper95CI")],1)

# Feng et al. comparison
round(out$symptomatic$antibody_at_VE_summary["0.5",],0)
round(out$symptomatic$antibody_at_VE_summary["0.6",],0)
round(out$symptomatic$antibody_at_VE_summary["0.7",],0)
round(out$symptomatic$antibody_at_VE_summary["0.8",],0)
round(out$symptomatic$antibody_at_VE_summary["0.9",],0)

# Generate correlated gamma and zeta
cox_stage_two <- cox_stage_two_symptomatic
nchains <- dim(cox_stage_two$cox_model_var)[2]
nsamples <- dim(cox_stage_two$cox_model_var)[1]
set.seed(1234) 
Z0_norm <- matrix(rnorm(n=nsamples*nchains),nrow=nsamples,ncol=nchains)
Z1_norm <- matrix(rnorm(n=nsamples*nchains),nrow=nsamples,ncol=nchains)
# Calculate the correlation
cor_gamma_zeta <- cox_stage_two$cox_model_var[,,"antibody.As_vaccinated_arm_2ChAdOx1"]/sqrt(cox_stage_two$cox_model_var[,,"antibody.antibody"]*cox_stage_two$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"])
gammas <- cox_stage_two$cox_model_pred[,,"antibody"] + Z0_norm*sqrt(cox_stage_two$cox_model_var[,,"antibody.antibody"])
zetas <- cox_stage_two$cox_model_pred[,,"As_vaccinated_arm_2ChAdOx1"] + (Z0_norm*cor_gamma_zeta + Z1_norm*sqrt(1-cor_gamma_zeta^2))*sqrt(cox_stage_two$cox_model_var[,,"As_vaccinated_arm_2ChAdOx1.As_vaccinated_arm_2ChAdOx1"])

gamma_quantiles <- quantile(gammas, c(0.5,0.025,0.975))
HR_log1increase <- function(A_vec,gamma_quants = gamma_quantiles){
  return(exp(outer(A_vec,gamma_quants)*(exp(1)-1)))
}
round(HR_log1increase(c(10000,1000,500,383,219,119,65,30,5,1)/conv_factor),2)



######
# Vaccine efficacy over time

# Comparison with Andrews et al.
times <- as.character(c(5,10,20)*7)
VE_vs_time_summary <- round(out$symptomatic$VE_vs_time[times,"mean_VE",],1)
print(paste0(paste0(paste0(VE_vs_time_summary[,"50%"],"% "),paste0("(",apply(VE_vs_time_summary[,c("2.5%","97.5%")],1,paste,collapse=",~"),")")),collapse = ", "))

VE_vs_time_summary <- round(out$symptomatic$VE_vs_time_new_ind[times,"baseline",],1)
print(paste0(paste0(paste0(VE_vs_time_summary[,"50%"],"% "),paste0("(",apply(VE_vs_time_summary[,c("2.5%","97.5%")],1,paste,collapse=",~"),")")),collapse = ", "))

times <- as.character(c(35,63))
VE_vs_time_summary <- round(out$symptomatic$VE_vs_time_new_ind[times,"56-69",],1)
print(paste0(paste0(paste0(VE_vs_time_summary[,"50%"],"% "),paste0("(",apply(VE_vs_time_summary[,c("2.5%","97.5%")],1,paste,collapse=",~"),")")),collapse = ", "))

VE_vs_time_summary <- round(out$symptomatic$VE_vs_time_new_ind[times,">=70",],1)
print(paste0(paste0(paste0(VE_vs_time_summary[,"50%"],"% "),paste0("(",apply(VE_vs_time_summary[,c("2.5%","97.5%")],1,paste,collapse=",~"),")")),collapse = ", "))

# Omicron
# Comparison with Hogan
dimnames(VE_vs_time_omicron)
effect <- dimnames(VE_vs_time_omicron)[[2]]
n.effect <- length(effect)
times_new_ind <- as.character(c(90,180))
for (s in 1:length(times_new_ind)){
  print(times_new_ind[s])
  print("\n")
  for (i in 1:n.effect){
    print((effect[i]))
    VE_vs_time_omicron_times <- round(VE_vs_time_omicron[times_new_ind[s],effect[i],],1)
    print(paste0("Day ",times_new_ind[s],": ",paste0(paste0(paste0(VE_vs_time_omicron_times["50%"],"% "),paste0("(",paste(VE_vs_time_omicron_times[c("2.5%","97.5%")],collapse=",~"),")")),collapse = ", ")))
  }
}

pryr::mem_used()
t_final <- Sys.time()
print(difftime(t_final,t_init))
