#############################################################################
# Antibody plots
# This file produces antibody plots from the longitudinal Bayesian mixed effects 
# log-linear antibody model in Longitudinal_antibody_model.R

t1 <- Sys.time()
print(t1)
# To run code on the server.

#############################################################################
antibody_type <- c("S","neuts")[1]
long_type <- c("7_day_incubator","0_day_incubator")[1]
# The name of the longitudinal model
file_name <- "ri_rs_7inc_pos_prim_exp_slope_t_all_covariates_hcw_cumhaz_t0PB28_6e4"
if(!(substring(file_name,1,5)=="neuts") == (antibody_type=="neuts")){warning("Your file_name and antibody_type do not match. You may be running a neuts analysis on S data or vice versa. Please double check.")}
name <- paste0("Correlates_long_",file_name)

# Directories
SLURM <- FALSE # Is this being run on SLURM? We need different directories depending on yes/no
if(SLURM){
  directory_correlates <- getwd()
  directory_correlates_project <- paste0(directory_correlates,"/3_Programs")
  output_directory <- paste0(directory_correlates,"/4_Output")
  output_directory2 <- output_directory
  main_plot_directory <- paste0(output_directory,"/Plots")
  # Create a plot folder for the specific model run in Longitudinal_antibody_model.R
  if (!dir.exists(main_plot_directory)){dir.create(main_plot_directory)}
  plot_directory <- paste0(main_plot_directory,"/",file_name)
  if (!dir.exists(plot_directory)){
    dir.create(plot_directory)
    dir.create(paste0(plot_directory,"/Antibody_plots"))
  }
} else{ # This file has code defining where the directories are saved. Not included in the GitHub.
  source(paste0(getwd(),"Read_directories_antibody_plots.R"))
}

#load packages
library(rstan)
library(Hmisc)
library(ggplot2)
library(gridExtra)

#######
# Load user-written functions
source(paste0(directory_correlates_project,"/functions.R"))
################################################################################
# Longitudinal model
# Read the longitudinal random effects
long_out_rs <- readRDS(paste0(output_directory,"/Correlates_long_",file_name,"_reffects.RDS"))
# File is a list, first entry is the Stan output, second entry the inputted data.
data_long <- long_out_rs[[2]]
long_out_rs <- long_out_rs[[1]]
chainsamples <- dim(long_out_rs)[1]
# Merge the chains into one matrix
# Find the population parameters for intercept and slope
a_0_array <- long_out_rs[,,which(substring(dimnames(long_out_rs)$parameters,1,2)=="a_0[")]
a_1_array <- long_out_rs[,,which(substring(dimnames(long_out_rs)$parameters,1,2)=="a_1[")]
# a_0_array gives the value at time 14 + t_0 days after the second dose.
n <- dim(a_0_array)[3]
nsamples <- dim(a_0_array)[1]
nchains <- dim(a_0_array)[2]
rm(long_out_rs)
###############
# Read joint_correlates and long_correlates
# Set their factors to be correct
source(paste0(directory_correlates_project,"/read_data_set_factors.R"))
if(antibody_type=="S"){ # If using anti-spike IgG
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
  log_antibody_lims <- c(3,14.8)
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
if(nantibody_LOD==nantibody_leq_LOD & nantibody_LOD>0){warning(paste0("There are ",nantibody_LOD," antibody observations equal to the limit of detection. We will treat these as though they are not censored, but observed exactly."))}
long_correlates$left_censored <- FALSE
long_correlates$left_censored[which(long_correlates$antibody<LOD)] <- TRUE
long_correlates$antibody[which(long_correlates$left_censored)] <- LOD
long_correlates$log_antibody[which(long_correlates$left_censored)] <- log(LOD)

# Defining variables
# long_correlates <- long_correlates[,c("redcap_data_access_group","As_vaccinated_arm_2","sc_repeat_pid","record_id","sc_age","age_group","cor2dose_primary","cor2dose_positive","antibody","log_antibody","antibody_time","ll","timepoint","visit2","left_censored","outlier","cor2dose_outcome")]
# Log antibody levels (untransformed units)
# joint_correlates <- joint_correlates[,c("redcap_data_access_group","As_vaccinated_arm_2","sc_repeat_pid","sc_age","age_group","cor2dose_schedule","cor2dose_interval","cor2dose_non_white","cor2dose_bmi_geq_30","cor2dose_comorbidities","cor2dose_hcw_status","cor2dose_primary","cor2dose_positive","start_time","end_time")]

# The data for vaccinated individuals only
joint_correlates_ChAd <- joint_correlates[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1"),]
# Set the column names for a_0_array and a_1_array
dimnames(a_0_array)[3]$parameters <- joint_correlates_ChAd$sc_repeat_pid
dimnames(a_1_array)[3]$parameters <- joint_correlates_ChAd$sc_repeat_pid

# times
max_time <- 209
min_time <- 28
adsb <- min_time:max_time # Days since boost for antibodies (ds_t0+14)
# Number of days from boost to day 0 in the model (i.e. day 0 in data + t_center the day 0 adjustment in longitudinal model)
if(long_type=="7_day_incubator") {
  raw_data_t0 <- 21 # start_time in long_correlates data
  VEdsb <- adsb + 7 # Days since boost (ds_t0+21)
} else if(long_type=="0_day_incubator") {raw_data_t0 <- 14; VEdsb <- adsb} # Days since boost (ds_t0+21)
t0 <- raw_data_t0 + data_long$t_center;
ds_t0 <- adsb-t0 # Days since day 0 in the model
nts <- length(ds_t0)
time_units <- c("Days","Weeks","Months")[1]
# div is the Amount to divide days by to get the time units we're using (1 if days, 7 if weeks, 365/12 if months)
if (time_units=="Days"){div <- 1}
if (time_units=="Weeks"){div <- 7}
if (time_units=="Months"){div <- 365/12}
atsb <- adsb/div # time since boost for antibodies in the appropriate units
long_correlates$antibody_time_since_boost <- (long_correlates$antibody_time+raw_data_t0)/div

# Raw data testing positive
table(joint_correlates[,c("cor2dose_outcome_pos_prim","As_vaccinated_arm_2")])
round(prop.table(table(joint_correlates[,c("cor2dose_outcome_pos_prim","As_vaccinated_arm_2")]),margin=2)*100,1)
table(joint_correlates[,c("cor2dose_positive","As_vaccinated_arm_2")])
round(prop.table(table(joint_correlates[,c("cor2dose_positive","As_vaccinated_arm_2")]),margin=2)*100,1)

#####
# antibody plots
# Raw data
# log_A_mat_extreme <- outer(a_0_array,rep(1,2)) + outer(a_1_array,ds_t0[c(1,nts)]) # Get all predicted antibody levels at the first and last times, to calculate the necessary limits for plots y-axis
log_A_mat <- outer(apply(a_0_array,c(3),mean),rep(1,nts)) + outer(apply(a_1_array,c(3),mean),ds_t0) # Mean antibody levels for each individual at each timepoint
# rm(log_A_mat_extreme)

time_axis_lab <- paste0(time_units," since second dose vaccination")
Antibody_lab <- paste0(Antibody_name," levels (",antibody_units,")")
antibody_lab <- paste0(antibody_name," levels (",antibody_units,")")
time_axis_marks <- c(0,30,60,90,120,150,180,210) # Choose the marks for the x axis in weeks
# Place the marks for the y-axis in powers of 10 BAU/mL converted to the log AU/mL antibody scale
antibody_ax_mark_loc <- log(10^seq(floor(log10(exp(log_antibody_lims)*conv_factor))[1],ceiling(log10(exp(log_antibody_lims)*conv_factor))[2],by=1)/conv_factor)
antibody_ax_mark_loc_subgrid <- log(10^seq(floor(log10(exp(log_antibody_lims)*conv_factor))[1],ceiling(log10(exp(log_antibody_lims)*conv_factor))[2],by=0.5)/conv_factor)
# Write marks for the y-axis in powers of 10 BAU/mL
antibody_ax_mark <- 10^seq(floor(log10(exp(log_antibody_lims)*conv_factor))[1],ceiling(log10(exp(log_antibody_lims)*conv_factor))[2],by=1)

# Marking outliers in the vaccinated group
long_correlates_Control <- long_correlates[which(long_correlates$As_vaccinated_arm_2=="Control"),]

# Plot observed antibody over time for control group amd vaccine group
# Observations below the limit of detection (LOD) have been given a value of LOD
cols <- data.frame("ChAdOx1"=c("black","grey"),
                   "Outlier"=c("red3","firebrick"),
                   "Control"=c("blue3","steelblue"))
rownames(cols) <- c("points","lines")
# Supplementary Figure 2
png(paste0(plot_directory,"/Antibody_plots/Observed_antibody_vs_time_incl_control&outliers.png"),width=1400,height=1000,pointsize=23)
par(mar=c(5.1,5.1,0.5,2.1))
plot(NA, cex = 0.2,ylim=c(min(antibody_ax_mark_loc_subgrid),log(10000/conv_factor)), yaxt="n", xlab=time_axis_lab,
     xlim=c(0,max(atsb)),col="NA", ylab = paste0("Observed ",antibody_lab),xaxt="n", cex.lab=1.3)
abline(v = seq(0,210,15),col="grey95")
abline(h = antibody_ax_mark_loc_subgrid,col="grey95")
abline(v = seq(0,210,30),col="grey90")
abline(h = antibody_ax_mark_loc,col="grey90")
for(ii in long_correlates[which(long_correlates$As_vaccinated_arm_2=="Control"),]$ll){ # Draw lines for the control individuals
  lines(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$ll==ii),], lwd=0.2,col=cols["lines","Control"])
} # Draw points for the Control individuals
points(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$As_vaccinated_arm_2=="Control"),], cex = 0.2,col=cols["points","Control"],pch=19)
for(ii in long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1" & !long_correlates$outlier),]$ll){ # Draw lines for the ChAdOx1 individuals who aren't outliers
  lines(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$ll==ii),], lwd=0.2, col = cols["lines","ChAdOx1"])
} # Draw points for the ChAdOx1 individuals who aren't outliers
points(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1" & !long_correlates$outlier),], cex = 0.2, col=cols["points","ChAdOx1"],pch=19)
for(ii in long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1" & long_correlates$outlier),]$ll){ # Draw lines for the ChAdOx1 outlying observations
  lines(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$ll==ii),], lwd=0.2,col=cols["lines","Outlier"])
}  # Draw points for the outliers
points(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1" & long_correlates$outlier),], cex = 0.5,col=cols["points","Outlier"],pch=19)
abline(h=log(LOD),lty=2,lwd=2)
axis(1,time_axis_marks, cex=1.5)
axis(2,antibody_ax_mark_loc,antibody_ax_mark)
if (sum(long_correlates$As_vaccinated_arm_2=="ChAdOx1" & long_correlates$outlier)>0){
  legend(x=150,y=log(10000/conv_factor)-0.5,lty=0,legend=c("ChAdOx1","ChAdOx1 outlier","Control"),col=as.character(cols["points",]),pch=19)
} else { # Legend with outliers included if there are any, else just ChAdOx1/Control
  legend(x=150,y=log(10000/conv_factor)-0.5,lty=0,legend=c("ChAdOx1","Control"),col=as.character(cols["points",c("ChAdOx1","Control")]),pch=19)
}
dev.off()

################
# Remove outliers from the dataset
outliers <- long_correlates[which(long_correlates$outlier),]
long_correlates <- long_correlates[which(!long_correlates$outlier),]
# Set those less than the LLOQ to be -1
long_correlates$antibody[which(long_correlates$antibody<=LOD)] <- -1/conv_factor
# Exclude the control group
long_correlates$visit2 <- factor(long_correlates$visit2, levels = c("PB28","PB90","PB182"))
levels(long_correlates$visit2) <- c("PB28","PB90","PB182")
long_correlates_Control <- long_correlates[which(long_correlates$As_vaccinated_arm_2=="Control"),]
long_correlates <- long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),]

#######################
# Create table of antibody observations at timepoints
visits.rownames <- rep(c(paste0("PB",c(28,90,182)),"Total"),each=3)
observations.rownames <- rep(c("n.available","median_IQR","n<LLOQ"),3)
table_observed_antibody <- matrix(nrow=15,ncol=7,dimnames = list(c("arm","case","n",paste(visits.rownames,observations.rownames,sep="_")),
                                                                 c("visit","observation","Non-case","Case","Symptomatic","ChAdOx1", "Control")))
table_observed_antibody[,] <- "~"
table_observed_antibody["arm",][c(3,7)] <- c("ChAdOx1","Control")
table_observed_antibody["case",][c(3:7)] <- c("Non-case","Case","Symptomatic case","Total ChAdOx1","Total Control")
table_observed_antibody[,"visit"][c(4,7,10,13)] <- c(paste0("PB",c(28,90,182)),"Total")
table_observed_antibody[,"observation"][4:15] <- rep(c("Number of antibody measurements available","Median (IQR) anti-spike IgG value BAU/mL","Number of antibody measurements < LLOQ"),4)
table_observed_antibody["n",][c(3,4)] <- paste0("n=",table(joint_correlates_ChAd$cor2dose_positive),
                                                " (",round(prop.table(table(joint_correlates_ChAd$cor2dose_positive))*100,1),"%)")
table_observed_antibody["n","Symptomatic"] <- paste0("n=",table(joint_correlates_ChAd$cor2dose_primary)["Case"],
                                                     " (",round(table(joint_correlates_ChAd$cor2dose_primary)["Case"]/n*100,1),"%)")
table_observed_antibody["n",][c(6,7)] <- paste0("n=",table(joint_correlates$As_vaccinated_arm_2)[c("ChAdOx1", "Control")])

n.table_observed_antibody <- c(table(joint_correlates_ChAd$cor2dose_positive),
                               table(joint_correlates_ChAd$cor2dose_primary)["Case"],
                               table(joint_correlates$As_vaccinated_arm_2)[c("ChAdOx1", "Control")])

table_observed_antibody[paste(paste0("PB",c(28,90,182)),"n.available",sep="_"),"ChAdOx1"] <- table(long_correlates$visit2)
table_observed_antibody[paste(paste0("PB",c(28,90,182)),"n.available",sep="_"),"Control"] <- table(long_correlates_Control$visit2)

for (visit in paste0("PB",c(28,90,182))){
  # Number of antibody measurements available
  table_observed_antibody[paste(visit,"n.available",sep="_"),][c(3,4)] <- table(long_correlates$cor2dose_positive[which(long_correlates$visit2==visit)])
  table_observed_antibody[paste(visit,"n.available",sep="_"),"Symptomatic"] <- table(long_correlates$cor2dose_primary[which(long_correlates$visit2==visit)])["Case"]
  visit.n.available <- as.numeric(table_observed_antibody[paste(visit,"n.available",sep="_"),3:7])
  table_observed_antibody[paste(visit,"n.available",sep="_"),3:7] <- 
    paste0(table_observed_antibody[paste(visit,"n.available",sep="_"),3:7]," (",
           round(as.numeric(table_observed_antibody[paste(visit,"n.available",sep="_"),3:7])/n.table_observed_antibody*100,1),"%)")
  # median and IQR of antibody observations
  visit_non_case_antibody <- long_correlates$antibody[which(long_correlates$visit2==visit & long_correlates$cor2dose_positive=="Negative")]*conv_factor
  table_observed_antibody[paste(visit,"median_IQR",sep="_"),]["Non-case"] <- paste0(round(median(visit_non_case_antibody),0)," (",paste0(round(quantile(visit_non_case_antibody,c(0.25,0.75)),0),collapse=", "),")")
  visit_case_antibody <- long_correlates$antibody[which(long_correlates$visit2==visit & long_correlates$cor2dose_positive=="Positive")]*conv_factor
  table_observed_antibody[paste(visit,"median_IQR",sep="_"),]["Case"] <- paste0(round(median(visit_case_antibody),0)," (",paste0(round(quantile(visit_case_antibody,c(0.25,0.75)),0),collapse=", "),")")
  visit_symp_antibody <- long_correlates$antibody[which(long_correlates$visit2==visit & long_correlates$cor2dose_primary=="Case")]*conv_factor
  table_observed_antibody[paste(visit,"median_IQR",sep="_"),]["Symptomatic"] <- paste0(round(median(visit_symp_antibody),0)," (",paste0(round(quantile(visit_symp_antibody,c(0.25,0.75)),0),collapse=", "),")")
  visit_ChAdOx1_antibody <- long_correlates$antibody[which(long_correlates$visit2==visit)]*conv_factor
  table_observed_antibody[paste(visit,"median_IQR",sep="_"),]["ChAdOx1"] <- paste0(round(median(visit_ChAdOx1_antibody),0)," (",paste0(round(quantile(visit_ChAdOx1_antibody,c(0.25,0.75)),0),collapse=", "),")")
  visit_Control_antibody <- long_correlates_Control$antibody[which(long_correlates_Control$visit2==visit)]*conv_factor
  table_observed_antibody[paste(visit,"median_IQR",sep="_"),]["Control"] <- paste0(round(median(visit_Control_antibody),1)," (",paste0(round(quantile(visit_Control_antibody,c(0.25,0.75)),1),collapse=", "),")")
  # Numbers <LLOQ
  table_observed_antibody[paste(visit,"n<LLOQ",sep="_"),]["Non-case"] <- sum(visit_non_case_antibody<0)
  table_observed_antibody[paste(visit,"n<LLOQ",sep="_"),]["Case"] <- sum(visit_case_antibody<0)
  table_observed_antibody[paste(visit,"n<LLOQ",sep="_"),]["Symptomatic"] <- sum(visit_symp_antibody<0)
  table_observed_antibody[paste(visit,"n<LLOQ",sep="_"),]["ChAdOx1"] <- sum(visit_ChAdOx1_antibody<0)
  table_observed_antibody[paste(visit,"n<LLOQ",sep="_"),]["Control"] <- sum(visit_Control_antibody<0)
  table_observed_antibody[paste(visit,"n<LLOQ",sep="_"),3:7] <- 
    paste0(table_observed_antibody[paste(visit,"n<LLOQ",sep="_"),3:7]," (",
           round(as.numeric(table_observed_antibody[paste(visit,"n<LLOQ",sep="_"),3:7])/visit.n.available*100,1),"%)")
}

table_observed_antibody[paste("Total","n.available",sep="_"),"ChAdOx1"] <- nrow(long_correlates)
table_observed_antibody[paste("Total","n.available",sep="_"),"Control"] <- nrow(long_correlates_Control)
# Number of antibody measurements available
table_observed_antibody[paste("Total","n.available",sep="_"),][c(3,4)] <- table(long_correlates$cor2dose_positive)
table_observed_antibody[paste("Total","n.available",sep="_"),"Symptomatic"] <- table(long_correlates$cor2dose_primary)["Case"]
visit.n.available <- as.numeric(table_observed_antibody[paste("Total","n.available",sep="_"),3:7])
table_observed_antibody[paste("Total","n.available",sep="_"),3:7] <- 
  paste0(table_observed_antibody[paste("Total","n.available",sep="_"),3:7]," (",
         round(as.numeric(table_observed_antibody[paste("Total","n.available",sep="_"),3:7])/n.table_observed_antibody*100,1),"%)")
# median and IQR of antibody observations
visit_non_case_antibody <- long_correlates$antibody[which(long_correlates$cor2dose_positive=="Negative")]*conv_factor
table_observed_antibody[paste("Total","median_IQR",sep="_"),]["Non-case"] <- paste0(round(median(visit_non_case_antibody),0)," (",paste0(round(quantile(visit_non_case_antibody,c(0.25,0.75)),0),collapse=", "),")")
visit_case_antibody <- long_correlates$antibody[which(long_correlates$cor2dose_positive=="Positive")]*conv_factor
table_observed_antibody[paste("Total","median_IQR",sep="_"),]["Case"] <- paste0(round(median(visit_case_antibody),0)," (",paste0(round(quantile(visit_case_antibody,c(0.25,0.75)),0),collapse=", "),")")
visit_symp_antibody <- long_correlates$antibody[which(long_correlates$cor2dose_primary=="Case")]*conv_factor
table_observed_antibody[paste("Total","median_IQR",sep="_"),]["Symptomatic"] <- paste0(round(median(visit_symp_antibody),0)," (",paste0(round(quantile(visit_symp_antibody,c(0.25,0.75)),0),collapse=", "),")")
visit_ChAdOx1_antibody <- long_correlates$antibody*conv_factor
table_observed_antibody[paste("Total","median_IQR",sep="_"),]["ChAdOx1"] <- paste0(round(median(visit_ChAdOx1_antibody),0)," (",paste0(round(quantile(visit_ChAdOx1_antibody,c(0.25,0.75)),0),collapse=", "),")")
visit_Control_antibody <- long_correlates_Control$antibody*conv_factor
table_observed_antibody[paste("Total","median_IQR",sep="_"),]["Control"] <- paste0(round(median(visit_Control_antibody),1)," (",paste0(round(quantile(visit_Control_antibody,c(0.25,0.75)),1),collapse=", "),")")
# Numbers <LLOQ
table_observed_antibody[paste("Total","n<LLOQ",sep="_"),]["Non-case"] <- sum(visit_non_case_antibody<0)
table_observed_antibody[paste("Total","n<LLOQ",sep="_"),]["Case"] <- sum(visit_case_antibody<0)
table_observed_antibody[paste("Total","n<LLOQ",sep="_"),]["Symptomatic"] <- sum(visit_symp_antibody<0)
table_observed_antibody[paste("Total","n<LLOQ",sep="_"),]["ChAdOx1"] <- sum(visit_ChAdOx1_antibody<0)
table_observed_antibody[paste("Total","n<LLOQ",sep="_"),]["Control"] <- sum(visit_Control_antibody<0)
table_observed_antibody[paste("Total","n<LLOQ",sep="_"),3:7] <- 
  paste0(table_observed_antibody[paste("Total","n<LLOQ",sep="_"),3:7]," (",
         round(as.numeric(table_observed_antibody[paste("Total","n<LLOQ",sep="_"),3:7])/visit.n.available*100,1),"%)")

# Set the -1 values to <LLOQ as these are those below the detection limit
table_observed_antibody <- gsub("-1","<LLOQ",table_observed_antibody)
# Set the IQR where there is only one observation to be NA
table_observed_antibody <- gsub("30, 30","NA",table_observed_antibody)
# Supplementary Table 2
write.table(table_observed_antibody,paste0(directory_correlates,"/4_Output/Tables/Antibody_observations.csv"),row.names=F,col.names = F, sep=",")

outliers_table <- outliers[,c("cor2dose_outcome_pos_prim","antibody_time_since_boost","visit2","antibody")]
outliers_table$antibody <- round(outliers_table$antibody*conv_factor,1)
outliers_table$antibody[which(outliers$antibody<=LOD)] <- "<LLOQ"
levels(outliers_table$cor2dose_outcome_pos_prim) <- c("Non-case","Asymptomatic case","Symptomatic case")
colnames(outliers_table) <- c("Case","Days since second dose","Visit","Anti-spike IgG value (BAU/mL)")
# Supplementary Table 3
write.csv(outliers_table,paste0(directory_correlates,"/4_Output/Tables/Antibody_outliers.csv"),row.names=F)

############
# # Plot observed antibody over time
# long_correlates$visit_col <- long_correlates$visit2
# visit_cols <- c("black","blue3","red3")
# levels(long_correlates$visit_col) <- visit_cols
# long_correlates$visit_col <- as.character(long_correlates$visit_col)
# 
# png(paste0(plot_directory,"/Antibody_plots/Observed_antibody_vs_time.png"),width=784,height=560)
# plot(NA, cex = 0.1,
#      main = paste0("Observed ",antibody_name," levels over time"), ylim=log_antibody_lims, yaxt="n", xlab=time_axis_lab,
#      xlim=c(0,max(atsb)),col="NA", ylab = paste0("Observed ",antibody_lab),xaxt="n")
# abline(v = seq(0,210,15),col="grey95")
# abline(h = antibody_ax_mark_loc_subgrid,col="grey95")
# abline(v = seq(0,210,30),col="grey90")
# abline(h = antibody_ax_mark_loc,col="grey90")
# for(ii in long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),]$ll){
#   lines(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$ll==ii),], lwd=0.2,col="grey")
# }
# points(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),], cex = 0.1, col = as.character(long_correlates$visit_col),pch=1)
# axis(1,time_axis_marks)
# axis(2,antibody_ax_mark_loc,antibody_ax_mark)
# abline(h=log(LOD),lty=2)
# # Hmisc::putKeyEmpty(x=long_correlates$antibody_time_since_boost,
# #             y=long_correlates$log_antibody,
# #             type="l",lty=0,labels=levels(long_correlates$visit2),col=visit_cols,pch=19, empty.method="area",
# #             xlim=c(0,max(atsb)),ylim=log_antibody_lims)
# legend(x=15,y=log(LOD)-0.1,col=visit_cols,legend=levels(long_correlates$visit2),pch=19)
# dev.off()

# # Plot predicted antibody over time
# png(paste0(plot_directory,"/Antibody_plots/Median_",file_name,".png"),width=784,height=560)
# plot(NA,xlim=c(0,max(atsb)),ylim=log_antibody_lims,xlab=time_axis_lab, ylab = paste0("Predicted ",antibody_lab),
#      main = paste0("Predicted ",antibody_name," levels over time"),xaxt="n",yaxt="n")
# abline(v = seq(0,210,15),col="grey95")
# abline(h = antibody_ax_mark_loc_subgrid,col="grey95")
# abline(v = seq(0,210,30),col="grey90")
# abline(h = antibody_ax_mark_loc,col="grey90")
# for (i in 1:n){
#   lines(log_A_mat[i,]~ atsb,col=c("grey"),lwd=0.2,lty=1)
# }
# axis(1,time_axis_marks)
# axis(2,antibody_ax_mark_loc,antibody_ax_mark)
# dev.off()

# Plot predicted and observed antibody over time next to each other

# png(paste0(plot_directory,"/Antibody_plots/Antibody_vs_time_obs_pred_mean.png"),width = 1800,height=1000, pointsize = 24)
# par(mfrow=c(1,2))
# plot(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),], cex = 0.1,
#      main = paste0("(a) Observed ",antibody_name," levels over time"), ylim=log_antibody_lims, yaxt="n", xlab=time_axis_lab,
#      xlim=c(0,max(atsb)),col="NA", ylab = paste0("Observed ",antibody_lab),xaxt="n")
# abline(v = seq(0,210,15),col="grey95")
# abline(h = antibody_ax_mark_loc_subgrid,col="grey95")
# abline(v = seq(0,210,30),col="grey90")
# abline(h = antibody_ax_mark_loc,col="grey90")
# for(ii in long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),]$ll){
#   lines(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$ll==ii),], lwd=0.2,col="grey")
# }
# points(log_antibody~antibody_time_since_boost,data = long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),], cex = 0.2, col = as.character(long_correlates$visit_col),pch=19)
# axis(1,time_axis_marks)
# axis(2,antibody_ax_mark_loc,antibody_ax_mark)
# abline(h=log(LOD),lty=2)
# legend(x=15,y=log(LOD)+2.5,col=visit_cols,legend=levels(long_correlates$visit2),pch=19)
# plot(NA,xlim=c(0,max(atsb)),ylim=log_antibody_lims,xlab=time_axis_lab, ylab = paste0("Predicted ",antibody_lab),
#      main = paste0("(b) Predicted ",antibody_name," levels over time"),xaxt="n",yaxt="n")
# abline(v = seq(0,210,15),col="grey95")
# abline(h = antibody_ax_mark_loc_subgrid,col="grey95")
# abline(v = seq(0,210,30),col="grey90")
# abline(h = antibody_ax_mark_loc,col="grey90")
# for (i in 1:n){
#   lines(log_A_mat[i,]~ atsb,col=c("grey"),lwd=0.2,lty=1)
# }
# axis(1,time_axis_marks)
# axis(2,antibody_ax_mark_loc,antibody_ax_mark)
# dev.off()

# Set characters to be factor variables
long_correlates[sapply(long_correlates, is.character)] <- lapply(long_correlates[sapply(long_correlates, is.character)], 
                                                                 as.factor)

# # Plot for different parameter samples
# set.seed(1)
# nsubsample <- 30
# sample_choices <- ceiling(seq(from=1, to =nsamples,by=nsamples/nsubsample*nchains)) # The second sample choice 335 contains the sample with the smallest value at day209 of them all, the very steep line which looks a bit off. Perhaps gives the impression the model is doing worse than it actually is.
# nsubsample <- length(sample_choices)
# 
# require(RColorBrewer)
# cols <- brewer.pal(n = 8, name = 'Set1')
# if (!dir.exists(paste0(plot_directory,"/Antibody_plots/Parameter_samples"))){
#   dir.create(paste0(plot_directory,"/Antibody_plots/Parameter_samples"))
# }
# set.seed(1)
# rand_order <- sample(n,n,replace=T)
# for (chain in 1:nchains){
#   for (s in 1:nsubsample){
#     log_A_mat_s <- outer(a_0_array[sample_choices[s],chain,],rep(1,nts)) + outer(a_1_array[sample_choices[s],chain,],ds_t0)
#     png(paste0(plot_directory,"/Antibody_plots/Parameter_samples/Parameter_samples_",file_name,"_chain_",chain,"_sample_",sample_choices[s],".png"), width=784, height=560)
#     plot(NA,xlim=c(0,max(atsb)),ylim=log_antibody_lims,xlab=time_axis_lab, ylab = paste0("Predicted ",antibody_lab),
#          main = paste0("Predicted ",antibody_name," levels over time"), xaxt="n",yaxt="n")
#     for (i in 8:n){
#       lines(log_A_mat_s[rand_order[i],]~atsb,col="grey",lwd=0.2,lty=1)
#     }
#     for (i in 1:7){
#       lines(log_A_mat_s[rand_order[i],]~atsb,col=cols[i],lwd=1,lty=1)
#     }
#     axis(1,time_axis_marks)
#     axis(2,antibody_ax_mark_loc,antibody_ax_mark)
#     dev.off()
#   }
#   
# }

# # Plot for individuals
# # Choose individuals corresponding to samples with largest and smallest 
# # intercepts, slopes, and values at d209
# a_0_min_ind <- arrayInd(which.min(a_0_array), dim(a_0_array))[3]
# a_0_max_ind <- arrayInd(which.max(a_0_array), dim(a_0_array))[3]
# a_1_min_ind <- arrayInd(which.min(a_1_array), dim(a_1_array))[3]
# a_1_max_ind <- arrayInd(which.max(a_1_array), dim(a_1_array))[3]
# end_min_ind <- arrayInd(which.min(a_0_array+a_1_array*ds_t0[nts]), dim(a_0_array+a_1_array*ds_t0[nts]))[3]
# end_max_ind <- arrayInd(which.max(a_0_array+a_1_array*ds_t0[nts]), dim(a_0_array+a_1_array*ds_t0[nts]))[3]
# minmax_names <- c("a_0_min","a_0_max","a_1_min","a_1_max","end_min","end_max")
# 
# png(paste0(plot_directory,"/Antibody_plots/Extreme_individuals_",file_name,".png"), width=784, height=560)
# par(mfrow=c(2,3), mar = c(1.5,2.5,2.5,0.1),oma=c(4.1,4.1,0,0))
# layout(mat = matrix(1:6,nrow=2,ncol=3))
# for (count in 1:length(minmax_names)){
#   j <- sapply(paste0(minmax_names,"_ind"),get)[count]
#   log_A_mat_s <- outer(a_0_array[,,j],rep(1,nts)) + outer(a_1_array[,,j],ds_t0)
#   log_A_mat_s_quants <- apply(log_A_mat_s,3,quantile,c(0.025,0.5,0.975))
#   plot(NA,xlim=c(0,max(atsb)),ylim=log_antibody_lims,xaxt="n",yaxt="n",xlab="",ylab="", main = paste0(minmax_names[count],"\nIndividual ",j))#,xlab=time_axis_lab, ylab = paste0("Predicted ",antibody_lab),
#   abline(v = seq(0,210,15),col="grey95")
#   abline(h = antibody_ax_mark_loc_subgrid,col="grey95")
#   abline(v = seq(0,210,30),col="grey90")
#   abline(h = antibody_ax_mark_loc,col="grey90")
#   #main = paste0("Predicted ",antibody_name," levels over time","\n individual ",j))
#   lines(log_A_mat_s_quants[2,]~atsb,lty=1)
#   lines(log_A_mat_s_quants[1,]~atsb,lty=2)
#   lines(log_A_mat_s_quants[3,]~atsb,lty=2)
#   long_j <- long_correlates[which(is.element(long_correlates$sc_repeat_pid,joint_correlates_ChAd$sc_repeat_pid[j])),]
#   long_j$atsb <- (long_j$antibody_time+14)/div
#   if (nrow(long_j)!=0) {
#     points(log_antibody~atsb, data=long_j, col = "blue3",pch=19)
#   }
#   if ((count%%2)==0){
#     axis(1,time_axis_marks)
#   } else{axis(1,time_axis_marks,rep("",length(time_axis_marks)))}
#   if ((count<3)){
#     axis(2,antibody_ax_mark_loc,antibody_ax_mark,las=1)
#   } else{axis(2,antibody_ax_mark_loc,rep("",length(antibody_ax_mark)))}
# }
# mtext(time_axis_lab,side=1,outer=TRUE,cex=1.3,line=2)
# mtext(paste0("Predicted ",antibody_lab),side=2,line=2,outer=TRUE,cex=1.3,las=0)
# dev.off()

# Choose 3 random invididuals with 0, 1, 2, 3 observations
set.seed(1234)
no_obs_ind <- sample(which(!is.element(joint_correlates_ChAd$sc_repeat_pid,long_correlates$sc_repeat_pid)),3)
one_obs_all_ind <- names(which(table(long_correlates$sc_repeat_pid)==1)) # Find the pids for individuals with one observation
one_obs_ind <- sample(which(is.element(joint_correlates_ChAd$sc_repeat_pid,one_obs_all_ind)),3) # Sample 3 individuals with one observation. Similar below.
two_obs_all_ind <- names(which(table(long_correlates$sc_repeat_pid)==2))
two_obs_ind <- sample(which(is.element(joint_correlates_ChAd$sc_repeat_pid,two_obs_all_ind)),3)
three_obs_all_ind <- names(which(table(long_correlates$sc_repeat_pid)==3))
three_obs_ind <- sample(which(is.element(joint_correlates_ChAd$sc_repeat_pid,three_obs_all_ind)),3)
inds_plot <- c(no_obs_ind,one_obs_ind,two_obs_ind,three_obs_ind)
# Supplementary Figure 12
png(paste0(plot_directory,"/Antibody_plots/Individuals_random_0123_",file_name,".png"), width=1400, height=1000, pointsize=23)
par(mfrow=c(3,4), mar = c(1.5,2.5,1.5,0.1),oma=c(4.1,4.1,0.1,0.5))
layout(mat = matrix(1:12,nrow=3,ncol=4))
for (count in 1:length(inds_plot)){
  j <- inds_plot[count]
  log_A_mat_s <- outer(a_0_array[,,j],rep(1,nts)) + outer(a_1_array[,,j],ds_t0)
  log_A_mat_s_quants <- apply(log_A_mat_s,3,quantile,c(0.025,0.5,0.975))
  plot(NA,xlim=c(0,max(atsb)),ylim=log_antibody_lims,xaxt="n",yaxt="n",xlab="",ylab="", main = paste0("Individual ",j))#,xlab=time_axis_lab, ylab = paste0("Predicted ",antibody_lab),
  abline(v = seq(0,210,15),col="grey95")
  abline(h = antibody_ax_mark_loc_subgrid,col="grey95")
  abline(v = seq(0,210,30),col="grey90")
  abline(h = antibody_ax_mark_loc,col="grey90")
  #main = paste0("Predicted ",antibody_name," levels over time","\n individual ",j))
  lines(log_A_mat_s_quants[2,]~atsb,lty=1,lwd=1.7)
  lines(log_A_mat_s_quants[1,]~atsb,lty=2,lwd=1.7)
  lines(log_A_mat_s_quants[3,]~atsb,lty=2,lwd=1.7)
  long_j <- long_correlates[which(is.element(long_correlates$sc_repeat_pid,joint_correlates_ChAd$sc_repeat_pid[j])),]
  long_j$atsb <- (long_j$antibody_time+14)/div
  if (nrow(long_j)!=0) {
    points(log_antibody~atsb, data=long_j, col = "blue3",pch=19)
  }
  if ((count%%3)==0){
    axis(1,time_axis_marks)
  } else{axis(1,time_axis_marks,rep("",length(time_axis_marks)))}
  if ((count<4)){
    axis(2,antibody_ax_mark_loc,antibody_ax_mark,las=1)
  } else{axis(2,antibody_ax_mark_loc,rep("",length(antibody_ax_mark)))}
}
mtext(time_axis_lab,side=1,outer=TRUE,cex=1.3,line=2)
mtext(paste0(Antibody_lab),side=2,line=2,outer=TRUE,cex=1.3,las=0)
dev.off()

# Data for these individuals
joint_correlates_ChAd[inds_plot,]
long_correlates[which(is.element(long_correlates$sc_repeat_pid,joint_correlates_ChAd$sc_repeat_pid[inds_plot])),]

# # Longitudinal residuals
# # Get predicted values at each observation
# a_0_pred <- a_0_array[,,as.character(long_correlates$sc_repeat_pid[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1")])]
# a_1_pred <- a_1_array[,,as.character(long_correlates$sc_repeat_pid[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1")])]
# 
# obs_times <- long_correlates$antibody_time - data_long$t_center
# pred_all <- a_0_pred + sweep(a_1_pred,3,obs_times,FUN="*")
# pred_means <- apply(pred_all,3,median)
# pred_sd <- apply(pred_all,3,FUN=sd)
# obs_atsb <- long_correlates$antibody_time_since_boost
# resid_post_pred <- (long_correlates$log_antibody-pred_means) # I think we don't want to divide by pred_sd
# resid_df <- data.frame("atsb" = obs_atsb,
#                        "residuals" = resid_post_pred,
#                        "ll" = long_correlates$ll)
# long_correlates$log_antibody[which(long_correlates$left_censored)] <- log(LOD)
# png(paste0(plot_directory,"/Antibody_plots/Longitudinal_residuals.png"),width=1400,height=1000,pointsize = 23)
# par(mar=c(5.1,4.1,0.5,0.5))
# plot(NA, ylim = max(abs(resid_post_pred))*c(-1,1),
#      ylab = paste0("Difference between predicted and observed log ",antibody_name),
#      xlab = time_axis_lab,xaxt="n",xlim=c(0,max(atsb)))
# abline(v = seq(0,210,15),col="grey95")
# abline(h = seq(-2,2,0.25),col="grey95")
# abline(v = seq(0,210,30),col="grey90")
# abline(h = seq(-2,2,0.5),col="grey90")
# for(ii in resid_df$ll){
#   lines(residuals~atsb,data = resid_df[which(resid_df$ll==ii),], lwd=0.5, col = alpha("grey",0.4))
# }
# points(resid_post_pred[which(!long_correlates$left_censored)]~obs_atsb[which(!long_correlates$left_censored)])
# points(resid_post_pred[which(long_correlates$left_censored)]~obs_atsb[which(long_correlates$left_censored)], col="blue")
# abline(h=0,lty=2)
# axis(1,time_axis_marks)
# if (sum(long_correlates$left_censored)>0){
#   Hmisc::putKeyEmpty(y=resid_post_pred,x=obs_atsb,type="p",pch=1,lty=0,labels=c("Observed exactly","Left-censored"),col=c("black","blue"))
# } # Include legend if left-censored observations are present.
# dev.off()
# visit <- c("PB28","PB90","PB182")[1]
# x <- resid_post_pred[which(long_correlates$visit2==visit)]
# left_censored <- long_correlates$left_censored[which(long_correlates$visit2==visit)]
# 
# if (sum(left_censored)>0){
#   require(EnvStats)
#   qqPlotCensored(x,censored=left_censored,add.line=T,
#                  distribution = "norm", main = paste0("Normal Q-Q plot for log antibody levels at visit ",visit))
#   t_df <- 4
#   qqPlotCensored(x,censored=left_censored,add.line=T,
#                  distribution = "t", param.list=list(df=t_df), main = paste0("Student-t Q-Q plot with ",t_df," dof for log antibody levels at visit ",visit))
# } else{
#   qqnorm(x)
#   qqline(x)
# }

###################
# Write observed antibody proportions and levels
# Check all are ChAdOx1
# No individuals have multiple observations at the same timepoint.
sum(duplicated(long_correlates[,c("sc_repeat_pid","visit2")]))
pos_table <- table(long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),c("visit2","cor2dose_positive")])
prim_table <- table(long_correlates[which(long_correlates$As_vaccinated_arm_2=="ChAdOx1"),c("visit2","cor2dose_primary")])
aprop_table <- cbind(pos_table,prim_table[,"Case"])
colnames(aprop_table)[3] <- "Primary"
aperc_table <- round(sweep(aprop_table,2,c(table(joint_correlates_ChAd$cor2dose_positive),table(joint_correlates_ChAd$cor2dose_primary)[2]),"/")*100,1)
aprop_table_print <- aprop_table
aprop_table_print[1:length(aprop_table_print)] <-  paste0(aprop_table, " (",aperc_table,"%)")
aprop_table_print <- t(aprop_table_print)
apply(aprop_table_print,1,function(x){paste(x, collapse=", ")})

table(long_correlates$As_vaccinated_arm_2)
antibody_numbers <- array(dim=c(3,3,2))
dimnames(antibody_numbers) <- list(c("Non-case","Positive","Primary"),c("PB28","PB90","PB182"),c("main","brackets"))
# Numbers of observations
long_correlates$antibody[which(long_correlates$left_censored)] <- -1/conv_factor # Set the censored values to be -1 international units
# Non-cases
# Median antibodies
antibody_numbers[1,,1] <- c(round(median(long_correlates$antibody[which(long_correlates$visit2=="PB28" & long_correlates$cor2dose_positive=="Negative")])*conv_factor,0),
                            round(median(long_correlates$antibody[which(long_correlates$visit2=="PB90" & long_correlates$cor2dose_positive=="Negative")])*conv_factor,0),
                            round(median(long_correlates$antibody[which(long_correlates$visit2=="PB182" & long_correlates$cor2dose_positive=="Negative")])*conv_factor,0))
# IQR antibodies
antibody_numbers[1,,2] <- c(paste(round(quantile(long_correlates$antibody[which(long_correlates$visit2=="PB28" & long_correlates$cor2dose_positive=="Negative")],c(0.25,0.75))*conv_factor,0), collapse=",~"),
                            paste(round(quantile(long_correlates$antibody[which(long_correlates$visit2=="PB90" & long_correlates$cor2dose_positive=="Negative")],c(0.25,0.75))*conv_factor,0), collapse=",~"),
                            paste(round(quantile(long_correlates$antibody[which(long_correlates$visit2=="PB182" & long_correlates$cor2dose_positive=="Negative")],c(0.25,0.75))*conv_factor,0), collapse=",~"))
# Positive cases
# Median antibodies
antibody_numbers[2,,1] <- c(round(median(long_correlates$antibody[which(long_correlates$visit2=="PB28" & long_correlates$cor2dose_positive=="Positive")])*conv_factor,0),
                            round(median(long_correlates$antibody[which(long_correlates$visit2=="PB90" & long_correlates$cor2dose_positive=="Positive")])*conv_factor,0),
                            round(median(long_correlates$antibody[which(long_correlates$visit2=="PB182" & long_correlates$cor2dose_positive=="Positive")])*conv_factor,0))
# IQR antibodies
antibody_numbers[2,,2] <- c(paste(round(quantile(long_correlates$antibody[which(long_correlates$visit2=="PB28" & long_correlates$cor2dose_positive=="Positive")],c(0.25,0.75))*conv_factor,0), collapse=",~"),
                            paste(round(quantile(long_correlates$antibody[which(long_correlates$visit2=="PB90" & long_correlates$cor2dose_positive=="Positive")],c(0.25,0.75))*conv_factor,0), collapse=",~"),
                            paste(round(quantile(long_correlates$antibody[which(long_correlates$visit2=="PB182" & long_correlates$cor2dose_positive=="Positive")],c(0.25,0.75))*conv_factor,0), collapse=",~"))
# Primary cases
# Median antibodies
antibody_numbers[3,,1] <- c(round(median(long_correlates$antibody[which(long_correlates$visit2=="PB28" & long_correlates$cor2dose_primary=="Case")])*conv_factor,0),
                            round(median(long_correlates$antibody[which(long_correlates$visit2=="PB90" & long_correlates$cor2dose_primary=="Case")])*conv_factor,0),
                            round(median(long_correlates$antibody[which(long_correlates$visit2=="PB182" & long_correlates$cor2dose_primary=="Case")])*conv_factor,0))
# IQR antibodies
antibody_numbers[3,,2] <- c(paste(round(quantile(long_correlates$antibody[which(long_correlates$visit2=="PB28" & long_correlates$cor2dose_primary=="Case")],c(0.25,0.75))*conv_factor,0), collapse=",~"),
                            paste(round(quantile(long_correlates$antibody[which(long_correlates$visit2=="PB90" & long_correlates$cor2dose_primary=="Case")],c(0.25,0.75))*conv_factor,0), collapse=",~"),
                            paste(round(quantile(long_correlates$antibody[which(long_correlates$visit2=="PB182" & long_correlates$cor2dose_primary=="Case")],c(0.25,0.75))*conv_factor,0), collapse=",~"))
antibody_numbers2 <- apply(antibody_numbers,c(1,2),paste,collapse = ", (")
antibody_numbers2[,] <- paste0(antibody_numbers2,")")
antibody_numbers2 <- gsub("-1","LOD",antibody_numbers2) # Replace the -1s with LOD so the censored values come out as being at the LOD
antibody_numbers3 <- apply(antibody_numbers2, 1, paste0, collapse=", ")
print(antibody_numbers3)

# Antibodies in control individuals
long_correlates_Control$antibody[which(long_correlates_Control$left_censored)] <- -1/conv_factor # Set the censored values to be -1 international units
nrow(long_correlates_Control)
length(unique(long_correlates_Control$sc_repeat_pid))
length(unique(long_correlates_Control$sc_repeat_pid))/nrow(joint_correlates_ChAd)*100
sum(long_correlates_Control$left_censored)
round(sum(long_correlates_Control$left_censored)/nrow(long_correlates_Control)*100,1)
median(long_correlates_Control$antibody)*conv_factor
quantile(long_correlates_Control$antibody,c(0.25,0.75))*conv_factor

# Half-lives
median_half_lives <- apply(-log(2)/a_1_array,c(1,2),median)
t(round(apply(median_half_lives,2,quantile,c(0.5,0.025,0.975)),0))
round(quantile(median_half_lives,c(0.5,0.025,0.975)),0)
Rhat(median_half_lives)
ess_tail(median_half_lives)
rm(median_half_lives)
# Predicted antibody values
pred_antibody_numbers <- array(dim=c(3,2))
dimnames(pred_antibody_numbers) <- list(c("PB28","PB90","PB182"),c("median","IQR"))
d28_median <- exp(apply(a_0_array + (28-t0)*a_1_array,c(1,2),median))*conv_factor
d90_median <- exp(apply(a_0_array + (90-t0)*a_1_array,c(1,2),median))*conv_factor
d182_median <- exp(apply(a_0_array + (182-t0)*a_1_array,c(1,2),median))*conv_factor
# Results from each chain
t(round(apply(d28_median,2,quantile,c(0.5,0.025,0.975)),1))
t(round(apply(d90_median,2,quantile,c(0.5,0.025,0.975)),1))
t(round(apply(d182_median,2,quantile,c(0.5,0.025,0.975)),1))
# Median and associated 95% credible interval
pred_antibody_numbers[,1] <- c(round(median(d28_median),0),
                               round(median(d90_median),0),
                               round(median(d182_median),0))
pred_antibody_numbers[,2] <- c(paste(round(quantile(d28_median,c(0.025,0.975)),0),collapse=",~"),
                               paste(round(quantile(d90_median,c(0.025,0.975)),0),collapse=",~"),
                               paste(round(quantile(d182_median,c(0.025,0.975)),0),collapse=",~"))
apply(log_A_mat[,which(is.element(atsb,c(28,90,182)))],2,quantile,c(0.25,0.5,0.75))
pred_antibody_numbers2 <- apply(pred_antibody_numbers,1,paste,collapse = " (")
(pred_antibody_numbers3 <- paste0(paste(pred_antibody_numbers2,collapse = "), "),")"))

# Antibody 25% and 75% quantiles at day 28, and associated credible interval
d28_quants <- exp(apply(a_0_array + (28-t0)*a_1_array,c(1,2),quantile,c(0.25,0.75)))*conv_factor
# Results from each chain
t(round(apply(d28_quants["25%",,],2,quantile,c(0.5,0.025,0.975)),1))
t(round(apply(d28_quants["75%",,],2,quantile,c(0.5,0.025,0.975)),1))
print("Larger deviation for the credible interval for the upper quantile (on the decimal places scale)")
antibody_quantiles_summary <- array(dim=c(2,2))
antibody_quantiles_summary[,1] <- c(round(median(d28_quants["25%",,]),0),
                                    round(median(d28_quants["75%",,]),0))
antibody_quantiles_summary[,2] <- c(paste(round(quantile(d28_quants["25%",,],c(0.025,0.975)),0), collapse=",~"),
                                    paste(round(quantile(d28_quants["75%",,],c(0.025,0.975)),0), collapse=",~"))
antibody_quantiles2 <- apply(antibody_quantiles_summary,1,paste,collapse = " (")
(antibody_quantiles3 <- paste0(paste(antibody_quantiles2,collapse = "), "),")"))
# Quantile of half-lives
half_lives <- -log(2)/a_1_array # Given in days
hl25 <- apply(half_lives,c(1,2),quantile,0.25)
hl75 <- apply(half_lives,c(1,2),quantile,0.75)
Rhat(hl25); ess_tail(hl25); Rhat(hl75); ess_tail(hl75)
t(round(apply(hl25,2,quantile,c(0.5,0.025,0.975)),1))
t(round(apply(hl75,2,quantile,c(0.5,0.025,0.975)),1))
round(median(hl25),0)
paste(round(quantile(hl25,c(0.025,0.975)),0),collapse = ",~")
round(median(hl75),0)
paste(round(quantile(hl75,c(0.025,0.975)),0),collapse = ",~")
rm(half_lives)

# Median PB28 for cases and non-cases
# Non-case
medianPB28noncase <- exp(apply((a_0_array + (28-t0)*a_1_array)[,,which(joint_correlates_ChAd$cor2dose_positive!="Positive")],c(1,2),median))*conv_factor
Rhat(medianPB28noncase)
ess_tail(medianPB28noncase)
t(round(apply(medianPB28noncase,2,quantile,c(0.5,0.025,0.975)),1))
round(median(medianPB28noncase),0)
paste(round(quantile(medianPB28noncase,c(0.025,0.975)),0),collapse = ",~")
# Positive
medianPB28case <- exp(apply((a_0_array + (28-t0)*a_1_array)[,,which(joint_correlates_ChAd$cor2dose_positive=="Positive")],c(1,2),median))*conv_factor
Rhat(medianPB28case)
ess_tail(medianPB28case)
t(round(apply(medianPB28case,2,quantile,c(0.5,0.025,0.975)),1))
print("Some variation here")
round(median(medianPB28case),0)
paste(round(quantile(medianPB28case,c(0.025,0.975)),0),collapse = ",~")
# Primary
medianPB28case <- exp(apply((a_0_array + (28-t0)*a_1_array)[,,which(joint_correlates_ChAd$cor2dose_primary=="Case")],c(1,2),median))*conv_factor
Rhat(medianPB28case)
ess_tail(medianPB28case)
t(round(apply(medianPB28case,2,quantile,c(0.5,0.025,0.975)),1))
print("Some variation here")
round(median(medianPB28case),0)
paste(round(quantile(medianPB28case,c(0.025,0.975)),0),collapse = ",~")


# Half-lives
# Non-case
half_lives_noncase <- apply((-log(2)/a_1_array)[,,which(joint_correlates_ChAd$cor2dose_positive!="Positive")],c(1,2),median)
Rhat(half_lives_noncase)
ess_tail(half_lives_noncase)
t(round(apply(half_lives_noncase,2,quantile,c(0.5,0.025,0.975)),0))
round(median(half_lives_noncase),0)
paste(round(quantile(half_lives_noncase,c(0.025,0.975)),0),collapse = ",~")
# Positive
half_lives_case <- apply((-log(2)/a_1_array)[,,which(joint_correlates_ChAd$cor2dose_positive=="Positive")],c(1,2),median)
Rhat(half_lives_case)
ess_tail(half_lives_case)
t(round(apply(half_lives_case,2,quantile,c(0.5,0.025,0.975)),0))
round(median(half_lives_case),0)
paste(round(quantile(half_lives_case,c(0.025,0.975)),0),collapse = ",~")
# Primary
half_lives_case <- apply((-log(2)/a_1_array)[,,which(joint_correlates_ChAd$cor2dose_primary=="Case")],c(1,2),median)
Rhat(half_lives_case)
ess_tail(half_lives_case)
t(round(apply(half_lives_case,2,quantile,c(0.5,0.025,0.975)),0))
round(median(half_lives_case),0)
paste(round(quantile(half_lives_case,c(0.025,0.975)),0),collapse = ",~")
#####
# Covariate effects on PB28 level and half-life
long_out <- readRDS(paste0(output_directory,"/Correlates_long_",file_name,"_noreffects.RDS"))
# File is a list, first entry is the Stan output, second entry the inputted data.
data_long <- long_out[[2]]
long_out <- long_out[[1]]
long_out_save <- long_out
dim(long_out)

rownames_full <- c("Age\n (56-69 vs 18-55 years)", "Age\n (\u226570 vs 18-55 years)", 
                   "Sex\n (Female vs Male)", "Ethnicity\n (Non-white vs White)",
                   "Comorbidity\n (Comorbidity vs None)", "BMI\n (\u226530 vs <30 kg/m^2)",
                   "Interval between first and second dose\n (9-11 weeks vs \u226512 weeks)",
                   "Interval between first and second dose\n (6-8 weeks vs \u226512 weeks)",
                   "Interval between first and second dose\n (<6 weeks vs \u226512 weeks)", 
                   "Healthcare worker (HCW) status\n (HCW facing < 1 COVID-19\n patient per day vs Not a HCW)",
                   "Healthcare worker (HCW) status\n (HCW worker facing \u2265 1\n COVID-19 patient per day vs Not a HCW)",
                   "First dose\n (Low dose vs Standard dose)",
                   "COVID-19 outcome\n (Positive vs Negative)",
                   "COVID-19 outcome\n (Primary symptomatic vs Negative)",
                   "COVID-19 outcome\n (Positive vs Negative)", # here twice as it now appears as a factor but used to appear as two separate indicators.
                   "COVID-19 outcome\n (Primary symptomatic vs Negative)",
                   "Age\n (18-69 years vs >=70)", # Come back and change this age - use 18-69 as baseline
                   "Age\n (\u226560 years vs 18-55 years)",
                   "Cumulative hazard\n (Nelson-Aalen estimator\n for cumulative hazard)") 
rowlabels_full <- c("age_group56-69","age_group>=70", "sc_genderFemale", "cor2dose_non_whiteOther",
                    "cor2dose_comorbiditiesComorbidity", "cor2dose_bmi_geq_30BMI_geq_30",
                    "cor2dose_interval9-11", "cor2dose_interval6-8", "cor2dose_interval<6",
                    "cor2dose_hcw_statusLess_than_1_covid", "cor2dose_hcw_statusMore_than_1_covid",
                    "cor2dose_primescheduleLD", "cor2dose_outcome_pos_primPositive", "cor2dose_outcome_pos_primPrimary",
                    "cor2dose_outcome_pos_primPositive", "cor2dose_outcome_pos_primPrimary",
                    "age_group_7018-69","age_group_56>=56","NelsonAalen")
row_reorder <- c("age_group56-69","age_group>=70", "sc_genderFemale", "cor2dose_non_whiteOther",
                 "cor2dose_comorbiditiesComorbidity", "cor2dose_bmi_geq_30BMI_geq_30",
                 "cor2dose_hcw_statusLess_than_1_covid", "cor2dose_hcw_statusMore_than_1_covid",
                 "cor2dose_interval9-11", "cor2dose_interval6-8", "cor2dose_interval<6",
                 "cor2dose_primescheduleLD", "cor2dose_outcome_pos_primPositive", "cor2dose_outcome_pos_primPrimary",
                 "NelsonAalen")
######
# # Plot the standardised parameters (same scale as priors)
# beta_0_tf_samples <- long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l0_")] # The effect on PB28
# beta_1_tf_samples <- long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l1_")] # Negative as larger parameter means smaller half life. The effect on half-life
# 
# beta_0_tf <- data.frame("median" = apply(beta_0_tf_samples,3,median),
#                         "lower" = apply(beta_0_tf_samples,3,quantile,0.025),
#                         "upper" = apply(beta_0_tf_samples,3,quantile,0.975))
# rownames(beta_0_tf) <- colnames(data_long$X_l)
# names(rownames_full) <- rowlabels_full
# beta_0_tf$rownames <- rownames_full[match(rownames(beta_0_tf),names(rownames_full))]
# beta_0_tf$rownames <- factor(beta_0_tf$rownames,levels=beta_0_tf$rownames)
# beta_0_tf <- beta_0_tf[row_reorder,]
# 
# # beta_0_tf <- exp(beta_0_tf_log) # This gives the multiplicative factors on the PB28 value
# 
# beta_1_tf <- data.frame("median" = apply(beta_1_tf_samples,3,median),
#                         "lower" = apply(beta_1_tf_samples,3,quantile,0.025),
#                         "upper" = apply(beta_1_tf_samples,3,quantile,0.975))
# rownames(beta_1_tf) <- colnames(data_long$X_l)
# # beta_1_tf <- exp(beta_1_tf_log) # This gives the multiplicative factors on the half life
# beta_1_tf$rownames <- rownames_full[match(rownames(beta_1_tf),names(rownames_full))]
# beta_1_tf$rownames <- factor(beta_1_tf$rownames,levels=beta_1_tf$rownames)
# beta_1_tf <- beta_1_tf[row_reorder,]
# # Create the plots
# beta_0_tf$rownames <- factor(beta_0_tf$rownames,levels=beta_0_tf$rownames[nrow(beta_1_tf):1])
# beta_1_tf$rownames <- factor(beta_1_tf$rownames,levels=beta_1_tf$rownames[nrow(beta_1_tf):1])
# plot1 <- ggplot(beta_0_tf, aes(x=median, y=rownames)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin=lower, xmax=upper), height=.3) +
#   geom_vline(xintercept=0, linetype=2) +
#   ggtitle("(a) Day 28 effect beta_0'") +
#   ylab(NULL) + xlab(paste0("Parameter for effect on PB28 ",antibody_name," level (",antibody_units,")")) +
#   theme_bw() + 
#   theme(axis.text=element_text(size=16),axis.title=element_text(size=16))
# plot2 <- ggplot(beta_1_tf, aes(x=median, y=rownames)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin=lower, xmax=upper), height=.3) +
#   geom_vline(xintercept=0, linetype=2) +
#   ggtitle("(b) Slope effect beta_1'") +
#   ylab(NULL) + xlab(paste0("Parameter for effect on ",antibody_name," slope (+ve means faster decay)")) +
#   theme_bw() + 
#   theme(axis.text=element_text(size=16),axis.title=element_text(size=16))
# png(paste0(plot_directory,"/Antibody_plots/Covariate_effects_beta_01_tf.png"),width=1700,height=1000)
# grid.arrange(plot1, plot2, ncol = 2)
# dev.off()

#######
long_out <- long_out_save
# Plot the multiplicative PB28 and half-life parameters
if (is.element("cor2dose_primaryCase",colnames(data_long$X_l))) {
  # Set the primary effect to be due to being both primary and positive (as it is impossible to be primary without being positive)
  long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l0[")][,,which(colnames(data_long$X_l)=="cor2dose_primaryCase")] <- rowSums(long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l0[")][,,which(is.element(colnames(data_long$X_l),c("cor2dose_primaryCase","cor2dose_positivePositive")))], dims=2)
  long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l1[")][,,which(colnames(data_long$X_l)=="cor2dose_primaryCase")] <- rowSums(long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l1[")][,,which(is.element(colnames(data_long$X_l),c("cor2dose_primaryCase","cor2dose_positivePositive")))], dims=2)
}
# Set the NelsonAalen effect to be due to an increase of 1 standard deviation
if (sum(colnames(data_long$X_l)=="NelsonAalen")!=0){
  long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l0[")][,,which(colnames(data_long$X_l)=="NelsonAalen")] <- long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l0[")][,,which(colnames(data_long$X_l)=="NelsonAalen")]*sd(data_long$X_l[,"NelsonAalen"])
  long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l1[")][,,which(colnames(data_long$X_l)=="NelsonAalen")] <- long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l1[")][,,which(colnames(data_long$X_l)=="NelsonAalen")]*sd(data_long$X_l[,"NelsonAalen"])
}
# Create data
beta_PB28_samples <- exp(long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l0[")]) # The effect on PB28
beta_halflife_samples <- exp(-long_out[,,which(substr(dimnames(long_out)$parameters,start=1,stop=8)=="beta_l1[")]) # Negative as larger parameter means smaller half life. The effect on half-life
round(apply(beta_PB28_samples,c(2,3),median),3)
round(apply(beta_PB28_samples,c(2,3),quantile,0.025),3)
round(apply(beta_PB28_samples,c(2,3),quantile,0.975),3)
round(apply(beta_halflife_samples,c(2,3),median),3)
round(apply(beta_halflife_samples,c(2,3),quantile,0.025),3)
round(apply(beta_halflife_samples,c(2,3),quantile,0.975),3)
beta_PB28 <- data.frame("median" = apply(beta_PB28_samples,3,median),
                        "mean" = apply(beta_PB28_samples,3,mean),
                        "lower" = apply(beta_PB28_samples,3,quantile,0.025),
                        "upper" = apply(beta_PB28_samples,3,quantile,0.975))
rownames(beta_PB28) <- colnames(data_long$X_l)
names(rownames_full) <- rowlabels_full
beta_PB28$rownames <- rownames_full[match(rownames(beta_PB28),names(rownames_full))]
beta_PB28$rownames <- factor(beta_PB28$rownames,levels=beta_PB28$rownames)
beta_PB28 <- beta_PB28[row_reorder,]

# beta_PB28 <- exp(beta_PB28_log) # This gives the multiplicative factors on the PB28 value

beta_halflife <- data.frame("median" = apply(beta_halflife_samples,3,median),
                            "mean" = apply(beta_halflife_samples,3,mean),
                            "lower" = apply(beta_halflife_samples,3,quantile,0.025),
                            "upper" = apply(beta_halflife_samples,3,quantile,0.975))
rownames(beta_halflife) <- colnames(data_long$X_l)
# beta_halflife <- exp(beta_halflife_log) # This gives the multiplicative factors on the half life
beta_halflife$rownames <- rownames_full[match(rownames(beta_halflife),names(rownames_full))]
beta_halflife$rownames <- factor(beta_halflife$rownames,levels=beta_halflife$rownames)
beta_halflife <- beta_halflife[row_reorder,]
# Create the plots
beta_PB28$rownames <- factor(beta_PB28$rownames,levels=beta_PB28$rownames[nrow(beta_halflife):1])
beta_halflife$rownames <- factor(beta_halflife$rownames,levels=beta_halflife$rownames[nrow(beta_halflife):1])
xmin <- floor(min(beta_PB28$lower)*10)/10
xmax <- ceiling((max(beta_PB28$upper)*10)+1)/10
plot1 <- ggplot(beta_PB28, aes(x=median, y=rownames)) +
  geom_point(cex=3) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.4, linewidth=0.7) +
  geom_vline(xintercept=1, linetype=2) +
  scale_x_log10(breaks =seq(xmin,xmax,0.1),limits = c(xmin,xmax)) +
  ggtitle("(a) Day 28 effect") +
  ylab(NULL) + xlab(paste0("Multiplicative effect on PB28 ",antibody_name," level (",antibody_units,")")) +
  theme_bw() +
  theme(text=element_text(size=16),axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.ticks.length.x = unit(0.3,"cm"))
xmin <- floor(min(beta_halflife$lower)*10)/10
xmax <- max(ceiling(max(beta_halflife$upper)*10)/10)
plot2 <- ggplot(beta_halflife, aes(x=median, y=rownames)) +
  geom_point(cex=3) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.4, linewidth=0.7) +
  geom_vline(xintercept=1, linetype=2) +
  scale_x_log10(breaks = seq(xmin,xmax,0.1),limits = c(xmin,xmax)) +
  ggtitle("(b) Half-life effect") +
  ylab(NULL) + xlab(paste0("Multiplicative effect on ",antibody_name," half-life (days)")) +
  theme_bw() +
  theme(text=element_text(size=16),axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16))
# Figure 4
png(paste0(plot_directory,"/Antibody_plots/Covariate_effects_PB28_halflife.png"),width=1700,height=1000)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()

#####
# Create table of PB28 and half-life effects
ncov <- nrow(beta_PB28)
beta_PB28[,"rownames"] <- as.character(beta_PB28[,"rownames"])
beta_halflife[,"rownames"] <- as.character(beta_halflife[,"rownames"])
beta_PB28_half_life <- rbind(c("Multiplicative effect on anti-spike IgG level 28 days post second dose",0,0,0,0),
                             beta_PB28[,c("rownames","median","mean","lower", "upper")],
                             c("Multiplicative effect on anti-spike IgG level half-life",0,0,0,0),
                             beta_halflife[,c("rownames","median","mean","lower", "upper")])
beta_PB28_half_life[,"rownames"] <- as.character(beta_PB28_half_life[,"rownames"])
beta_PB28_half_life[,c("median","mean","lower", "upper")] <- round(as.numeric(as.matrix(beta_PB28_half_life[,c("median","mean","lower", "upper")])),2)

beta_PB28_half_life[c(1,ncov+2),c("median","mean","lower", "upper")] <- rep(c("Median","Mean","0.025 quantile","0.975 quantile"),each=2)
colnames(beta_PB28_half_life) <- 
  c("Covariate multiplicative effect",
    c("Posterior estimates","","",""))
# Supplementary Table 5
write.table(beta_PB28_half_life,paste0(directory_correlates,"/4_Output/Tables/Covariate_effects_PB28_half_life.csv"),row.names=F,sep=",")

# Create table of parameter outputs from longitudinal model
long_out_parameters <- long_out_save[,,c("alpha_0","alpha_1","sigma_e","tau_0","tau_1","rho")]
# Convert to BAU/mL
long_out_parameters[,,"alpha_0"] <- long_out_parameters[,,"alpha_0"] + log(conv_factor)

parameter_output <- cbind(apply(long_out_parameters,3,mean),apply(long_out_parameters,3,sd)^2,t(apply(long_out_parameters,3,quantile,c(0.5,0.025,0.975))))
colnames(parameter_output) <- c("Mean","Variance","Median","0.025 quantile","0.975 quantile")
parameter_output <- round(parameter_output,3)
# Supplementary Table 4
write.csv(parameter_output,paste0(directory_correlates,"/4_Output/Tables/Longitudinal_model_parameters.csv"),row.names=T)

# Create table of parameter outputs from longitudinal model including covariate effects
long_out_parameters <- long_out_save[,,c("alpha_0","alpha_1","sigma_e","tau_0","tau_1","rho",
                                         dimnames(long_out_save)$parameters[which(substr(dimnames(long_out_save)$parameters,start=1,stop=8)=="beta_l0[")],
                                         dimnames(long_out_save)$parameters[which(substr(dimnames(long_out_save)$parameters,start=1,stop=8)=="beta_l1[")])]
# Convert to BAU/mL
long_out_parameters[,,"alpha_0"] <- long_out_parameters[,,"alpha_0"] + log(conv_factor)

parameter_output <- cbind(apply(long_out_parameters,3,mean),apply(long_out_parameters,3,sd)^2,t(apply(long_out_parameters,3,quantile,c(0.5,0.025,0.975))))
colnames(parameter_output) <- c("Mean","Variance","Median","0.025 quantile","0.975 quantile")
parameter_output <- round(parameter_output,5)

write.csv(parameter_output,paste0(directory_correlates,"/4_Output/Tables/Longitudinal_model_parameters_including_covariate.csv"),row.names=T)

#########

round(beta_PB28[c("cor2dose_interval9-11","cor2dose_interval6-8","cor2dose_interval<6"),c("median","lower","upper")]*100)
paste(paste0(round(beta_PB28[c("cor2dose_interval9-11","cor2dose_interval6-8","cor2dose_interval<6"),c("median")]*100),"% ",paste0("(",apply(round(beta_PB28[c("cor2dose_interval9-11","cor2dose_interval6-8","cor2dose_interval<6"),c("lower","upper")]*100),1,paste0,collapse=",~"),")")),collapse=", ")

pryr::mem_used()
print("File complete:")
t3 <- Sys.time()
print(difftime(t3,t1))
