# Read joint_correlates and long_correlates_S
# Set their factors to be correct
require(survival)
data_cleaning_commit_name <- "0c126098abccf39b82d0d8285b8659cd10c5b355"

# Baseline data (one row per participant)
joint_correlates <- read.csv(paste0(directory_correlates,"/2_Data_for_Analysis/joint_correlates_",data_cleaning_commit_name,".csv"))
# Longitudinal antibody observations (one row per observation)
long_correlates_S <- read.csv(paste0(directory_correlates,"/2_Data_for_Analysis/long_correlates_S_",data_cleaning_commit_name,".csv"))
# For pseudoneutralising antibodies
long_correlates_neuts <- read.csv(paste0(directory_correlates,"/2_Data_for_Analysis/long_correlates_neuts_",data_cleaning_commit_name,".csv"))

# Set factors
joint_correlates$site <- factor(joint_correlates$redcap_data_access_group)
joint_correlates$site <- relevel(joint_correlates$site,"oxford__ccvtm")
joint_correlates$As_vaccinated_arm_2 <- factor(joint_correlates$As_vaccinated_arm_2,levels = c("Control","ChAdOx1"))
joint_correlates$age_group <- factor(joint_correlates$age_group, levels = c("18-55","56-69",">=70"))
joint_correlates$cor2dose_schedule <- factor(joint_correlates$cor2dose_schedule, levels = c("SDSD","LDLD","LDSD"))
joint_correlates$cor2dose_primeschedule <- factor(joint_correlates$cor2dose_primeschedule,levels=c("SD","LD"))
joint_correlates$cor2dose_interval <- factor(joint_correlates$cor2dose_interval, levels = c(">=12","9-11","6-8","<6"))# Set baseline to be >=12 as this is what has mainly been used in practise.
joint_correlates$cor2dose_interval_comb <- factor(joint_correlates$cor2dose_interval_comb, levels = c(">=9","<8"))
joint_correlates$sc_gender <- factor(joint_correlates$sc_gender, levels = c("Male","Female"))
joint_correlates$cor2dose_non_white <- factor(joint_correlates$cor2dose_non_white, levels = c("White","Other"))
joint_correlates$cor2dose_comorbidities <- factor(joint_correlates$cor2dose_comorbidities, levels = c("None","Comorbidity"))
joint_correlates$cor2dose_bmi_geq_30 <- factor(joint_correlates$cor2dose_bmi_geq_30, levels = c("BMI_less_than_30","BMI_geq_30"))
joint_correlates$cor2dose_outcome <- factor(joint_correlates$cor2dose_outcome, levels = c("Negative","Primary","Non-Primary","Asymptomatic","Unknown"))
joint_correlates$cor2dose_outcome_pos_prim <- factor(joint_correlates$cor2dose_positive_ind+joint_correlates$cor2dose_primary_ind, levels = c(0,1,2))
levels(joint_correlates$cor2dose_outcome_pos_prim) <- c("Negative","Positive","Primary")
joint_correlates$cor2dose_positive_ind <- factor(joint_correlates$cor2dose_positive_ind, levels = c(0,1))
joint_correlates$cor2dose_positive <- factor(joint_correlates$cor2dose_positive, levels = c("Negative","Positive"))
joint_correlates$cor2dose_primary_ind <- factor(joint_correlates$cor2dose_primary_ind, levels = c(0,1))
joint_correlates$cor2dose_primary <- factor(joint_correlates$cor2dose_primary, levels = c("Non-case","Case"))
joint_correlates$cor2dose_hcw_status <- factor(joint_correlates$cor2dose_hcw_status, levels = c("Not_hcw","Less_than_1_covid","More_than_1_covid"))
joint_correlates$geographic_region <- factor(joint_correlates$geographic_region,levels=c("London","SouthEast","SouthWest","East","EastMidlands","WestMidlands","Wales","Yorkshire","NorthWest","NorthEast","Scotland"))
joint_correlates$geographic_big_region <- factor(joint_correlates$geographic_big_region,levels=c("London","SouthEnglandandWales","MidlandsandYorkshireEngland","NorthEnglandandScotland"))

# Create a variable for the end_time in terms of calendar time
joint_correlates$end_cal_time <- joint_correlates$start_time + joint_correlates$end_time
# A simple model to estimate the cumulative hazard
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
# Create a NelsonAalen estimator variable
joint_correlates$NelsonAalen <- NA
for (i in 1:nrow(joint_correlates)){
  joint_correlates$NelsonAalen[i] <- sum(discrete_haz$cumhaz[which(discrete_haz$time>joint_correlates$start_time[i] &
                                                                     discrete_haz$time<=joint_correlates$end_cal_time[i] &
                                                                     discrete_haz$As_vaccinated_arm_2==joint_correlates$As_vaccinated_arm_2[i])])
}

long_correlates_S$age_group <- factor(long_correlates_S$age_group, levels = c("18-55","56-69",">=70"))
long_correlates_S$cor2dose_schedule <- factor(long_correlates_S$cor2dose_schedule, levels = c("SDSD","LDLD","LDSD"))
long_correlates_S$cor2dose_primeschedule <- factor(long_correlates_S$cor2dose_primeschedule,levels=c("SD","LD"))
long_correlates_S$cor2dose_interval <- factor(long_correlates_S$cor2dose_interval, levels = c(">=12","9-11","6-8","<6"))# Set baseline to be >=12 as this is what has mainly been used in practise.
long_correlates_S$cor2dose_interval_comb <- factor(long_correlates_S$cor2dose_interval_comb, levels = c(">=9","<8"))
long_correlates_S$sc_gender <- factor(long_correlates_S$sc_gender, levels = c("Male","Female"))
long_correlates_S$cor2dose_non_white <- factor(long_correlates_S$cor2dose_non_white, levels = c("White","Other"))
long_correlates_S$cor2dose_comorbidities <- factor(long_correlates_S$cor2dose_comorbidities, levels = c("None","Comorbidity"))
long_correlates_S$cor2dose_bmi_geq_30 <- factor(long_correlates_S$cor2dose_bmi_geq_30, levels = c("BMI_less_than_30","BMI_geq_30"))
long_correlates_S$cor2dose_outcome <- factor(long_correlates_S$cor2dose_outcome, levels = c("Negative","Primary","Non-Primary","Asymptomatic","Unknown"))
long_correlates_S$cor2dose_outcome_pos_prim <- factor(long_correlates_S$cor2dose_positive_ind+long_correlates_S$cor2dose_primary_ind, levels = c(0,1,2))
levels(long_correlates_S$cor2dose_outcome_pos_prim) <- c("Negative","Positive","Primary")
long_correlates_S$cor2dose_positive_ind <- factor(long_correlates_S$cor2dose_positive_ind, levels = c(0,1))
long_correlates_S$cor2dose_positive <- factor(long_correlates_S$cor2dose_positive, levels = c("Negative","Positive"))
long_correlates_S$cor2dose_primary_ind <- factor(long_correlates_S$cor2dose_primary_ind, levels = c(0,1))
long_correlates_S$cor2dose_primary <- factor(long_correlates_S$cor2dose_primary, levels = c("Non-case","Case"))
long_correlates_S$cor2dose_hcw_status <- factor(long_correlates_S$cor2dose_hcw_status, levels = c("Not_hcw","Less_than_1_covid","More_than_1_covid"))
long_correlates_S$geographic_region <- factor(long_correlates_S$geographic_region,levels=c("London","SouthEast","SouthWest","East","EastMidlands","WestMidlands","Wales","Yorkshire","NorthWest","NorthEast","Scotland"))
long_correlates_S$geographic_big_region <- factor(long_correlates_S$geographic_big_region,levels=c("London","SouthEnglandandWales","MidlandsandYorkshireEngland","NorthEnglandandScotland"))
long_correlates_S$visit2 <- factor(long_correlates_S$visit2, levels = c("PB+28","PB+90","PB+182"))
levels(long_correlates_S$visit2) <- c("PB28","PB90","PB182")

long_correlates_neuts$age_group <- factor(long_correlates_neuts$age_group, levels = c("18-55","56-69",">=70"))
long_correlates_neuts$cor2dose_schedule <- factor(long_correlates_neuts$cor2dose_schedule, levels = c("SDSD","LDLD","LDSD"))
long_correlates_neuts$cor2dose_primeschedule <- factor(long_correlates_neuts$cor2dose_primeschedule,levels=c("SD","LD"))
long_correlates_neuts$cor2dose_interval <- factor(long_correlates_neuts$cor2dose_interval, levels = c(">=12","9-11","6-8","<6"))# Set baseline to be >=12 as this is what has mainly been used in practise.
long_correlates_neuts$cor2dose_interval_comb <- factor(long_correlates_neuts$cor2dose_interval_comb, levels = c(">=9","<8"))
long_correlates_neuts$sc_gender <- factor(long_correlates_neuts$sc_gender, levels = c("Male","Female"))
long_correlates_neuts$cor2dose_non_white <- factor(long_correlates_neuts$cor2dose_non_white, levels = c("White","Other"))
long_correlates_neuts$cor2dose_comorbidities <- factor(long_correlates_neuts$cor2dose_comorbidities, levels = c("None","Comorbidity"))
long_correlates_neuts$cor2dose_bmi_geq_30 <- factor(long_correlates_neuts$cor2dose_bmi_geq_30, levels = c("BMI_less_than_30","BMI_geq_30"))
long_correlates_neuts$cor2dose_outcome <- factor(long_correlates_neuts$cor2dose_outcome, levels = c("Negative","Primary","Non-Primary","Asymptomatic","Unknown"))
long_correlates_neuts$cor2dose_outcome_pos_prim <- factor(long_correlates_neuts$cor2dose_positive_ind+long_correlates_neuts$cor2dose_primary_ind, levels = c(0,1,2))
levels(long_correlates_neuts$cor2dose_outcome_pos_prim) <- c("Negative","Positive","Primary")
long_correlates_neuts$cor2dose_positive_ind <- factor(long_correlates_neuts$cor2dose_positive_ind, levels = c(0,1))
long_correlates_neuts$cor2dose_positive <- factor(long_correlates_neuts$cor2dose_positive, levels = c("Negative","Positive"))
long_correlates_neuts$cor2dose_primary_ind <- factor(long_correlates_neuts$cor2dose_primary_ind, levels = c(0,1))
long_correlates_neuts$cor2dose_primary <- factor(long_correlates_neuts$cor2dose_primary, levels = c("Non-case","Case"))
long_correlates_neuts$cor2dose_hcw_status <- factor(long_correlates_neuts$cor2dose_hcw_status, levels = c("Not_hcw","Less_than_1_covid","More_than_1_covid"))
long_correlates_neuts$geographic_region <- factor(long_correlates_neuts$geographic_region,levels=c("London","SouthEast","SouthWest","East","EastMidlands","WestMidlands","Wales","Yorkshire","NorthWest","NorthEast","Scotland"))
long_correlates_neuts$geographic_big_region <- factor(long_correlates_neuts$geographic_big_region,levels=c("London","SouthEnglandandWales","MidlandsandYorkshireEngland","NorthEnglandandScotland"))
long_correlates_neuts$visit2 <- factor(long_correlates_neuts$visit2, levels = c("PB+28","PB+90","PB+182"))
levels(long_correlates_neuts$visit2) <- c("PB28","PB90","PB182")
