#####
# These create a list of data to be inputted into the stan models

# Creates a list of data to input into Longitudinal_hierarchical_model_normalised_ri_rs_exp_slope_t.stan
stan_long_model <- function(long_data, log_y_var, obs_time_var, individual_data, 
                            longitudinal_covariates, ind_var, t_center = "mean"){
  n <- nrow(individual_data)
  
  X_l_df <- as.data.frame(individual_data[,longitudinal_covariates])
  X_l <- as.matrix(model.matrix(~.,data=X_l_df)[,-1],ncol=length(longitudinal_covariates))
  ll <- match(long_data[,c(ind_var)],
              individual_data[,c(ind_var)],
              incomparables = NA) # ll[i] is the individual corresponding to the i th antibody observation
  n_l <- nrow(long_data)
  p_l <- ncol(X_l)
  if (t_center=="mean") {t_center <- mean(long_data[,obs_time_var])}
  
  data <- list("n" =n, # number of individuals
               "n_l" = n_l, # number of antibody observations
               "p_l" = p_l, # number of covariates in X_l
               "ll" = ll,
               "X_l" = X_l, # covariate matrix
               "log_y" = long_data[,log_y_var], # the observed log antibodies
               "t_l" = long_data[,obs_time_var], # the timings of the observations
               "t_center" = t_center) # the transformed time 0 in the running of the stan model to improve mixing
  return(data)
}

stan_joint_model <- function(long_data, y_var, obs_time_var, individual_data, 
                             survival_covariates, longitudinal_covariates,
                             start_time_var, end_time_var, event_var, ind_var,
                             knots, h_0_prior_mean, h_0_prior_sd, t_center = "mean",
                             fixed_basehaz_fn = NULL){
  require(splines)
  n <- nrow(individual_data)
  start_time <- min(individual_data[,start_time_var])
  end_time <- max(individual_data[,end_time_var] + individual_data[,start_time_var])
  #####
  q_nodes.unadjusted <- c(
    -0.991455371120813,
    -0.949107912342758,
    -0.864864423359769,
    -0.741531185599394,
    -0.586087235467691,
    -0.405845151377397,
    -0.207784955007898,
    0,
    0.207784955007898,
    0.405845151377397,
    0.586087235467691,
    0.741531185599394,
    0.864864423359769,
    0.949107912342758,
    0.991455371120813
  )
  q_weights.unadjusted <- c(
    0.022935322010529,
    0.063092092629979,
    0.10479001032225,
    0.140653259715526,
    0.169004726639268,
    0.190350578064785,
    0.204432940075299,
    0.209482141084728,
    0.204432940075299,
    0.190350578064785,
    0.169004726639268,
    0.140653259715526,
    0.10479001032225,
    0.063092092629979,
    0.022935322010529
  )
  q <- length(q_nodes.unadjusted)
  #####
  vacc_times <- individual_data[,start_time_var]
  t <- individual_data[,end_time_var]
  ind_e <- which(individual_data[,event_var]==1)
  t_event <- t[ind_e]
  n_event <- length(ind_e)
  
  X_s_df <- as.data.frame(individual_data[,survival_covariates])
  X_s <- as.matrix(model.matrix(~.,data=X_s_df)[,-1],ncol=length(survival_covariates))
  X_l_df <- as.data.frame(individual_data[,longitudinal_covariates])
  X_l <- as.matrix(model.matrix(~.,data=X_l_df)[,-1],ncol=length(longitudinal_covariates))
  ll <- match(long_data[,c(ind_var)],
              individual_data[,c(ind_var)],
              incomparables = NA)
  if(sum(is.na(ll))>0) error("Not all longitudinal data rows have a corresponding row in the individual data")
  n_l <- nrow(long_data)
  p_l <- ncol(X_l)
  p_s <- ncol(X_s)
  if (t_center=="mean") {t_center <- mean(long_data[,obs_time_var])}
  
  q_nodes <- (matrix(t,n,1)/2)%*%(matrix(q_nodes.unadjusted,1,q)+1) + matrix(vacc_times,n,q)
  q_nodes <- c(t(q_nodes))
  q_weights <- (matrix(t,n,1)/2)%*%matrix(q_weights.unadjusted,1,q)
  q_weights <- c(t(q_weights))
  # q_nodes_scaled <- (b-a)/2*(q_nodes+1) + a
  # q_weights_scaled <- (b-a)/2*q_weights
  eval_points <- c(t_event + vacc_times[ind_e],q_nodes)
  B_q <- bs(eval_points,knots = knots, intercept = F, Boundary.knots = c(start_time,end_time))
  m <- dim(B_q)[2]
  
  data <- list("n" =n,
               "n_event" = n_event,
               "n_l" = n_l,
               "p_l" = p_l,
               "ll" = ll,
               "X_l" = X_l,
               "y" = long_data[,y_var],
               "t_l" = long_data[,obs_time_var],
               "p_s" = p_s,
               "m" = m,
               "X_s" = X_s,
               "t" = t,
               "t_event" = t_event,
               "event" = individual_data[,event_var],
               "q" = q,
               "n_times_q" = n*q,
               "n_eval" = length(eval_points),
               "t_q" = q_nodes,
               "q_weights" = q_weights,
               "B_q" = B_q,
               "ind_q" = rep(1:n,each=q),
               "ind_e" = ind_e,
               "h_0_prior_mean" = h_0_prior_mean,
               "h_0_prior_sd" = h_0_prior_sd,
               "t_center" = t_center)
  if (!is.null(fixed_basehaz_fn)){
    data$log_basehaz_etimes <- log(predict(fixed_basehaz_fn,t_event)$y)-mean(log(predict(fixed_basehaz_fn,q_nodes)$y))
    data$log_basehaz_qtimes <- log(predict(fixed_basehaz_fn,q_nodes)$y)-mean(log(predict(fixed_basehaz_fn,q_nodes)$y))
  }
  return(data)
}

# Includes vacc_var an indicator that they received the vaccine/treatment and so have longitudinal data available.
# 0 if received control, 1 if received vaccine.
# Separates the survival outcome into vaccine and control component.
stan_joint_model_vaccine <- function(long_data, y_var, obs_time_var, individual_data, 
                                     survival_covariates, longitudinal_covariates,
                                     start_time_var, end_time_var, event_var, ind_var, vacc_var,
                                     knots, h_0_prior_mean, h_0_prior_sd, t_center = "mean", fixed_basehaz_fn = NULL){
  require(splines)
  n <- nrow(individual_data)
  start_time <- min(individual_data[,start_time_var])
  end_time <- max(individual_data[,end_time_var] + individual_data[,start_time_var])
  #####
  q_nodes.unadjusted <- c(
    -0.991455371120813,
    -0.949107912342758,
    -0.864864423359769,
    -0.741531185599394,
    -0.586087235467691,
    -0.405845151377397,
    -0.207784955007898,
    0,
    0.207784955007898,
    0.405845151377397,
    0.586087235467691,
    0.741531185599394,
    0.864864423359769,
    0.949107912342758,
    0.991455371120813
  )
  q_weights.unadjusted <- c(
    0.022935322010529,
    0.063092092629979,
    0.10479001032225,
    0.140653259715526,
    0.169004726639268,
    0.190350578064785,
    0.204432940075299,
    0.209482141084728,
    0.204432940075299,
    0.190350578064785,
    0.169004726639268,
    0.140653259715526,
    0.10479001032225,
    0.063092092629979,
    0.022935322010529
  )
  q <- length(q_nodes.unadjusted)
  #####
  vacc_times <- individual_data[,start_time_var]
  t <- individual_data[,end_time_var]
  ind_e <- which(individual_data[,event_var]==1)
  t_event <- t[ind_e]
  n_event <- length(ind_e)
  ind_v <- which(individual_data[,vacc_var]==1)
  n_v <- sum(individual_data[,vacc_var])
  ind_ve <- which(individual_data[which(individual_data[,vacc_var]==1),event_var]==1) # Among those who had were vaccinated, who had an event 
  n_event_v <- length(ind_ve)
  ind_ev <- which(individual_data[which(individual_data[,event_var]==1),vacc_var]==1) # Among those who had an event, who was vaccinated
  ind_q <- rep(1:n,each=q)
  ind_vq <- rep(1:n_v,each=q)
  ind_qv <- which(is.element(ind_q,ind_vq))
  
  X_s_df <- as.data.frame(individual_data[,survival_covariates])
  X_s <- as.matrix(model.matrix(~.,data=X_s_df)[,-1],ncol=length(survival_covariates))
  X_l_df <- as.data.frame(individual_data[which(as.logical(ind_v)),longitudinal_covariates])
  X_l <- as.matrix(model.matrix(~.,data=X_l_df)[,-1],ncol=length(longitudinal_covariates))
  ll <- match(long_data[,c(ind_var)],
              individual_data[ind_v,c(ind_var)],
              incomparables = NA)
  n_l <- nrow(long_data)
  p_l <- ncol(X_l)
  p_s <- ncol(X_s)
  if (t_center=="mean") {t_center <- mean(long_data[,obs_time_var])}
  
  q_nodes <- (matrix(t,n,1)/2)%*%(matrix(q_nodes.unadjusted,1,q)+1) + matrix(vacc_times,n,q)
  q_nodes <- c(t(q_nodes))
  q_weights <- (matrix(t,n,1)/2)%*%matrix(q_weights.unadjusted,1,q)
  q_weights <- c(t(q_weights))
  # q_nodes_scaled <- (b-a)/2*(q_nodes+1) + a
  # q_weights_scaled <- (b-a)/2*q_weights
  eval_points <- c(t_event + vacc_times[ind_e],q_nodes)
  B_q <- bs(eval_points,knots = knots, intercept = F, Boundary.knots = c(start_time,end_time))
  m <- dim(B_q)[2]
  
  data <- list("n" =n,
               "n_v" = n_v,
               "n_event" = n_event,
               "n_event_v" = n_event_v,
               "n_l" = n_l,
               "p_l" = p_l,
               "ll" = ll,
               "X_l" = X_l,
               "log_y" = long_data[,y_var],
               "t_l" = long_data[,obs_time_var],
               "p_s" = p_s,
               "m" = m,
               "X_s" = X_s,
               "t" = t,
               "t_event" = t_event,
               "event" = individual_data[,event_var],
               "q" = q,
               "n_times_q" = n*q,
               "n_v_times_q" = n_v*q,
               "n_eval" = length(eval_points),
               "t_q" = q_nodes,
               "q_weights" = q_weights,
               "B_q" = B_q,
               "ind_q" = ind_q,
               "ind_e" = ind_e,
               "ind_ve" = ind_ve,
               "ind_ev" = ind_ev,
               "ind_vq" = ind_vq,
               "ind_qv" = ind_qv,
               "h_0_prior_mean" = h_0_prior_mean,
               "h_0_prior_sd" = h_0_prior_sd,
               "t_center" = t_center)
  if (!is.null(fixed_basehaz_fn)){
    data$log_basehaz_etimes <- log(predict(fixed_basehaz_fn,t_event)$y) - mean(log(predict(fixed_basehaz_fn,q_nodes)$y))
    data$log_basehaz_qtimes <- log(predict(fixed_basehaz_fn,q_nodes)$y) - mean(log(predict(fixed_basehaz_fn,q_nodes)$y))
  }
  return(data)
}

# Includes vacc_var an indicator that they received the vaccine/treatment and so have longitudinal data available.
# 0 if received control, 1 if received vaccine.
# Separates the survival outcome into vaccine and control component.
stan_joint_model_vaccine_separate <- function(long_data, y_var, obs_time_var, individual_data, 
                                              survival_covariates, longitudinal_covariates,
                                              start_time_var, end_time_var, event_var, ind_var, vacc_var,
                                              knots, h_0_prior_mean, h_0_prior_sd, t_center = "mean", fixed_basehaz_fn = NULL){
  require(splines)
  n <- nrow(individual_data)
  start_time <- min(individual_data[,start_time_var])
  end_time <- max(individual_data[,end_time_var] + individual_data[,start_time_var])
  #####
  q_nodes.unadjusted <- c(
    -0.991455371120813,
    -0.949107912342758,
    -0.864864423359769,
    -0.741531185599394,
    -0.586087235467691,
    -0.405845151377397,
    -0.207784955007898,
    0,
    0.207784955007898,
    0.405845151377397,
    0.586087235467691,
    0.741531185599394,
    0.864864423359769,
    0.949107912342758,
    0.991455371120813
  )
  q_weights.unadjusted <- c(
    0.022935322010529,
    0.063092092629979,
    0.10479001032225,
    0.140653259715526,
    0.169004726639268,
    0.190350578064785,
    0.204432940075299,
    0.209482141084728,
    0.204432940075299,
    0.190350578064785,
    0.169004726639268,
    0.140653259715526,
    0.10479001032225,
    0.063092092629979,
    0.022935322010529
  )
  q <- length(q_nodes.unadjusted)
  #####
  vacc_times <- individual_data[,start_time_var]
  t <- individual_data[,end_time_var]
  ind_e <- which(individual_data[,event_var]==1)
  t_event <- t[ind_e]
  n_event <- length(ind_e)
  ind_v <- which(individual_data[,vacc_var]==1)
  ind_c <- which(individual_data[,vacc_var]==0)
  n_v <- sum(individual_data[,vacc_var])
  n_c <- n-n_v
  ind_e_v <- which(individual_data[which(individual_data[,vacc_var]==1),event_var]==1) # Among those who were vaccinated, who had an event 
  ind_e_c <- which(individual_data[which(individual_data[,vacc_var]==0),event_var]==1) # Among those who were not vaccinated, who had an event 
  n_v_event <- length(ind_e_v)
  ind_q <- rep(1:n,each=q)
  ind_q_v <- rep(1:n_v,each=q)
  ind_q_c <- rep(1:n_c,each=q)
  
  X_s_df <- as.data.frame(individual_data[,survival_covariates])
  X_s <- as.matrix(model.matrix(~.,data=X_s_df)[,-1],ncol=length(survival_covariates))
  X_l_df <- as.data.frame(individual_data[which(as.logical(ind_v)),longitudinal_covariates])
  X_l <- as.matrix(model.matrix(~.,data=X_l_df)[,-1],ncol=length(longitudinal_covariates))
  ll <- match(long_data[,c(ind_var)],
              individual_data[ind_v,c(ind_var)],
              incomparables = NA)
  n_l <- nrow(long_data)
  p_l <- ncol(X_l)
  p_s <- ncol(X_s)
  if (t_center=="mean") {t_center <- mean(long_data[,obs_time_var])}
  
  q_nodes <- (matrix(t,n,1)/2)%*%(matrix(q_nodes.unadjusted,1,q)+1) + matrix(vacc_times,n,q)
  q_nodes <- c(t(q_nodes))
  q_weights <- (matrix(t,n,1)/2)%*%matrix(q_weights.unadjusted,1,q)
  q_weights <- c(t(q_weights))
  # q_nodes_scaled <- (b-a)/2*(q_nodes+1) + a
  # q_weights_scaled <- (b-a)/2*q_weights
  eval_points_v <- c(t_event + vacc_times[ind_e],q_nodes)
  B_q <- bs(eval_points,knots = knots, intercept = F, Boundary.knots = c(start_time,end_time))
  
  m <- dim(B_q)[2]
  
  data <- list("n" =n,
               "n_v" = n_v,
               "n_event" = n_event,
               "n_v_event" = n_v_event,
               "n_l" = n_l,
               "p_l" = p_l,
               "ll" = ll,
               "X_l" = X_l,
               "log_y" = long_data[,y_var],
               "t_l" = long_data[,obs_time_var],
               "p_s" = p_s,
               "m" = m,
               "X_s" = X_s,
               "t" = t,
               "t_event" = t_event,
               "event" = individual_data[,event_var],
               "q" = q,
               "n_times_q" = n*q,
               "n_v_times_q" = n_v*q,
               "n_eval" = length(eval_points),
               "t_q" = q_nodes,
               "q_weights" = q_weights,
               "B_q" = B_q,
               "ind_q" = ind_q,
               "ind_e" = ind_e,
               "ind_e_v" = ind_e_v,
               "ind_vq" = ind_vq,
               "ind_qv" = ind_qv,
               "h_0_prior_mean" = h_0_prior_mean,
               "h_0_prior_sd" = h_0_prior_sd,
               "t_center" = t_center)
  if (!is.null(fixed_basehaz_fn)){
    data$log_basehaz_etimes <- log(predict(fixed_basehaz_fn,t_event)$y) - mean(log(predict(fixed_basehaz_fn,q_nodes)$y))
    data$log_basehaz_qtimes <- log(predict(fixed_basehaz_fn,q_nodes)$y) - mean(log(predict(fixed_basehaz_fn,q_nodes)$y))
  }
  return(data)
}


# Sets up the data for the longitudinal stan model with left censored observations present
# Datapoints less than or equal to lower_limit are treated as left-censored at lower_limit (observed to be less than or equal to lower limit)
# All other data is assumed to be observed exactly. Other aspects of the function are like stan_long_model.
stan_long_model_left_censor <- function(long_data, y_var, obs_time_var, individual_data, 
                                        longitudinal_covariates, ind_var, t_center = "mean", lower_limit){
  n <- nrow(individual_data)
  
  X_l_df <- as.data.frame(individual_data[,longitudinal_covariates])
  X_l <- as.matrix(model.matrix(~.,data=X_l_df)[,-1],ncol=length(longitudinal_covariates))
  ll <- match(long_data[,c(ind_var)],
              individual_data[,c(ind_var)],
              incomparables = NA)
  n_l <- nrow(long_data)
  p_l <- ncol(X_l)
  if (t_center=="mean") {t_center <- mean(long_data[,obs_time_var])}
  which_censored_y <- which(long_data[,y_var]<=lower_limit) # The indices of the censored values of y
  which_exact_y <- which(long_data[,y_var]>lower_limit) # The indices of the non-censored values of y (observed exactly)
  
  data <- list("n" =n,
               "n_l" = n_l,
               "p_l" = p_l,
               "ll" = ll,
               "X_l" = X_l,
               "y" = long_data[,y_var],
               "t_l" = long_data[,obs_time_var],
               "t_center" = t_center,
               "n_cens" = length(which_censored_y),
               "which_censored_y" = which_censored_y,
               "which_exact_y" = which_exact_y,
               "lower_limit" = lower_limit)
  return(data)
}

stan_surv_model <- function(long_data, individual_data, 
                            survival_covariates,
                            start_time_var, end_time_var, event_var,
                            knots, h_0_prior_mean, h_0_prior_sd, t_center = "mean",
                            fixed_basehaz_fn = NULL){
  require(splines)
  n <- nrow(individual_data)
  start_time <- min(individual_data[,start_time_var])
  end_time <- max(individual_data[,end_time_var] + individual_data[,start_time_var])
  #####
  q_nodes.unadjusted <- c(
    -0.991455371120813,
    -0.949107912342758,
    -0.864864423359769,
    -0.741531185599394,
    -0.586087235467691,
    -0.405845151377397,
    -0.207784955007898,
    0,
    0.207784955007898,
    0.405845151377397,
    0.586087235467691,
    0.741531185599394,
    0.864864423359769,
    0.949107912342758,
    0.991455371120813
  )
  q_weights.unadjusted <- c(
    0.022935322010529,
    0.063092092629979,
    0.10479001032225,
    0.140653259715526,
    0.169004726639268,
    0.190350578064785,
    0.204432940075299,
    0.209482141084728,
    0.204432940075299,
    0.190350578064785,
    0.169004726639268,
    0.140653259715526,
    0.10479001032225,
    0.063092092629979,
    0.022935322010529
  )
  q <- length(q_nodes.unadjusted)
  #####
  vacc_times <- individual_data[,start_time_var]
  t <- individual_data[,end_time_var]
  ind_e <- which(individual_data[,event_var]==1)
  t_event <- t[ind_e]
  n_event <- length(ind_e)
  
  X_s_df <- as.data.frame(individual_data[,survival_covariates])
  X_s <- as.matrix(model.matrix(~.,data=X_s_df)[,-1],ncol=length(survival_covariates))
  p_s <- ncol(X_s)
  
  q_nodes <- (matrix(t,n,1)/2)%*%(matrix(q_nodes.unadjusted,1,q)+1) + matrix(vacc_times,n,q)
  q_nodes <- c(t(q_nodes))
  q_weights <- (matrix(t,n,1)/2)%*%matrix(q_weights.unadjusted,1,q)
  q_weights <- c(t(q_weights))
  # q_nodes_scaled <- (b-a)/2*(q_nodes+1) + a
  # q_weights_scaled <- (b-a)/2*q_weights
  eval_points <- c(t_event + vacc_times[ind_e],q_nodes)
  B_q <- bs(eval_points,knots = knots, intercept = F, Boundary.knots = c(start_time,end_time))
  m <- dim(B_q)[2]
  
  data <- list("n" =n,
               "n_event" = n_event,
               "p_s" = p_s,
               "m" = m,
               "X_s" = X_s,
               "t" = t,
               "t_event" = t_event,
               "event" = individual_data[,event_var],
               "q" = q,
               "n_times_q" = n*q,
               "n_eval" = length(eval_points),
               "t_q" = q_nodes,
               "q_weights" = q_weights,
               "B_q" = B_q,
               "ind_q" = rep(1:n,each=q),
               "ind_e" = ind_e,
               "h_0_prior_mean" = h_0_prior_mean,
               "h_0_prior_sd" = h_0_prior_sd,
               "t_center" = t_center)
  if (!is.null(fixed_basehaz_fn)){
    data$log_basehaz_etimes <- log(predict(fixed_basehaz_fn,t_event)$y)-mean(log(predict(fixed_basehaz_fn,q_nodes)$y))
    data$log_basehaz_qtimes <- log(predict(fixed_basehaz_fn,q_nodes)$y)-mean(log(predict(fixed_basehaz_fn,q_nodes)$y))
  }
  return(data)
}

# Print a point estimate and associated CIs from summary matrix x
# x must have columns containing the point_estimates and CI lower and upper bound
# The _col variables give either the name or the number for the associated column
print_CIs <- function(x,point_estimate_col = c("median","mean")[1],lowerCI_col = "lower95CI",upperCI_col = "upper95CI",
                      add_char="",initial_explanation=T,interval_type=c("credible","confidence")[1],alpha=0.05){
  if (interval_type=="credible"){
    interval_name <- "CrI"
  }
  else if (interval_type=="confidence"){
    interval_name <- "CI"
  }
  intervals <- apply(x[,c(lowerCI_col,upperCI_col)],1,paste,collapse=",~")
  point_estimates <- paste0(x[,point_estimate_col],add_char," ")
  if (initial_explanation){intervals[1] <- paste0(100*(1-alpha),"% ",interval_name," ",intervals[1])}
  paste0(paste0(point_estimates,paste0("(",intervals,")")),collapse = ", ")
}

# Gauss-Kronrod nodes and weights
q_nodes.unadjusted <- c(
  -0.991455371120813,
  -0.949107912342758,
  -0.864864423359769,
  -0.741531185599394,
  -0.586087235467691,
  -0.405845151377397,
  -0.207784955007898,
  0,
  0.207784955007898,
  0.405845151377397,
  0.586087235467691,
  0.741531185599394,
  0.864864423359769,
  0.949107912342758,
  0.991455371120813
)
q_weights.unadjusted <- c(
  0.022935322010529,
  0.063092092629979,
  0.10479001032225,
  0.140653259715526,
  0.169004726639268,
  0.190350578064785,
  0.204432940075299,
  0.209482141084728,
  0.204432940075299,
  0.190350578064785,
  0.169004726639268,
  0.140653259715526,
  0.10479001032225,
  0.063092092629979,
  0.022935322010529
)
q <- length(q_nodes.unadjusted)

# trace plots
# samples - extract(out, inc_warmup=T,permuted=F), iter the total number of iterations
# warmup - the no of warm up iterations, plot_grid - par(mfrow=plot_grid)
# discard - any chains to discard in final plot
# m - the number of spline parameters
# p_l - the number of longitudinal (fixed-effect) covariates
# p_s - the number of survival (time-fixed) covariates
# Plots three plots for each parameter - trace plots for all iterations and all chains on left,
# trace plots after burnin for all chains in the middle, trace plots after burnin for
# all chains except discard on the right.
my_trace_plots <- function(samples, m=NULL, p_l=NULL, p_s=NULL,
                           parameters = NULL, tf = F,
                           iter=NULL, warmup=NULL, discard=NULL, incl_warmup=F, plot_row = NULL, nchains=4,
                           order = NULL, model = "jm", xi=T,rs=T,ri=T){
  if (is.null(parameters)){
    if (model == "jm"){
      parameters <- get_jm_par_names(m,p_l,p_s,tf, xi_incl=xi, rs_incl=rs,ri_incl=ri)
    }
    else if (model == "surv"){
      parameters <- get_surv_par_names(m,p_s,tf)
    }
    else if (model == "long"){
      parameters <- get_long_par_names(p_l,tf)
    }
    else stop('model must be equal to"jm", "surv" or "long", for a joint model, survival (proportional hazards) model, and longitudinal hierarchical model respectively' )
  }
  if (is.null(iter)){
    iter = dim(samples)[1]
  }
  if (is.null(warmup)){
    warmup = iter/2
  }
  chains <- 1:nchains
  keep_chains <- chains[which(!is.element(chains,discard))]
  if (is.null(plot_row)){
    plot_grid <- c(2,1)
  }
  else {
    plot_grid <- c(plot_row,1)
  }
  if (is.null(order)) order <- 1:nchains
  
  par(mfrow=plot_grid)
  if (!is.null(discard)){
    if (incl_warmup){
      for (param in parameters){
        plot(NA,xlim = c(1,iter), ylim = c(min(samples[,,param]),max(samples[,,param])), main = paste0(param," w/ warmup"), ylab="")
        for (i in order){
          lines(samples[,i,param],col=i)
        }
      }
    }
    for (param in parameters){
      plot(NA,xlim = c(warmup+1,iter), ylim = c(min(samples[(warmup+1):(iter),,param]),max(samples[(warmup+1):(iter),,param])), main = paste0(param," all chains"), ylab="")
      for (i in order){
        lines(samples[(warmup+1):(iter),i,param]~c((warmup+1):(iter)),col=i)
      }
    }
    for (param in parameters){
      plot(NA,xlim = c(warmup+1,iter), ylim = c(min(samples[(warmup+1):(iter),keep_chains,param]),max(samples[(warmup+1):(iter),keep_chains,param])), main = paste0(param," kept chains"), ylab="")
      for (i in (order)[keep_chains]){
        lines(samples[(warmup+1):(iter),i,param]~c((warmup+1):(iter)),col=i)
      }
    }
  }
  else{
    if (incl_warmup){
      for (param in parameters){
        plot(NA,xlim = c(1,iter), ylim = c(min(samples[,,param]),max(samples[,,param])), main = paste0(param," w/ warmup"), ylab="")
        for (i in order){
          lines(samples[,i,param],col=i)
        }
      }
    }
    for (param in parameters){
      plot(NA,xlim = c(warmup+1,iter), ylim = c(min(samples[(warmup+1):(iter),,param]),max(samples[(warmup+1):(iter),,param])), main = paste0(param), ylab="")
      for (i in order){
        lines(samples[(warmup+1):(iter),i,param]~c((warmup+1):(iter)),col=i)
      }
    }
  }
  par(mfrow=c(1,1))
}


#######
# These are mostly redundant and not used 
joint.model.log.lik <- function(X_l,X_s,y,t_l,t_event,t_q,beta_l,
                                alpha_0,tau_0,eta_0,alpha_1,tau_1,rho,eta_1,
                                sigma_e,h_0,xi,B_q,beta_s,gamma,q_weights,
                                center, t_center=NULL){
  if(isTRUE(center)){ # Note this might well produce a different normalising constant. Should compare like with like where possible.
    X_l <- sweep(X_l, 2, colMeans(X_l), FUN = "-")
    X_l <- sweep(X_l, 2, apply(X_l,2,sd), FUN = "/")
    X_s <- sweep(X_s, 2, colMeans(X_s), FUN = "-")
    X_s <- sweep(X_s, 2, apply(X_s,2,sd), FUN = "/")
    y <- (y - mean(y))/sd(y)
    t_center <- mean(t_l)
    t_sd <- sd(t_l)
    if (is.null(t_center)) {t_center <- mean(t_l)}
    t_l <- (t_l - t_center)/t_sd
    t_event <- (t_event - t_center)/t_sd
    t_q <- (t_q - t_center)/t_sd
  }
  else if (!(isFALSE(center))){ stop("Center must be true of false")}
  log.lik <- 0
  a <- X_l%*%beta_l + alpha_0 + tau_0*eta_0
  b <- alpha_1 + tau_1*(rho*eta_0 + sqrt(1-rho^2)*eta_1)
  mu <- a[ll] + t_l * b[ll]
  log.lik <- log.lik + sum(log(dnorm(y, mean = mu, sd = sigma_e)))
  
  eta_s <- X_s%*%beta_s
  log_basehaz = B_q%*%xi;
  log_basehaz_etimes = head(log_basehaz,n_event);
  log_basehaz_qtimes = tail(log_basehaz,n_times_q);
  log_haz_etimes = log_basehaz_etimes + h_0 + eta_s[ind_e] + gamma * (a[ind_e] + b[ind_e] * t_event);
  log_haz_qtimes = log_basehaz_qtimes + h_0 + eta_s[ind_q] + gamma * (a[ind_q] + b[ind_q] * t_q);
  log.lik = log.lik + sum(log_haz_etimes) - q_weights%*%exp(log_haz_qtimes);
  return(log.lik)
}
get_jm_par_names <- function(m, p_l, p_s, tf=F, xi_incl,rs_incl,ri_incl){
  if (tf) {
    parameters <- "h_0_tf"
    if (xi_incl) {parameters <- c(parameters,paste0("xi[",1:m,"]"))}
    parameters <- c(parameters,paste0("beta_s_tf[",1:p_s,"]"),paste0("beta_l_tf[",1:p_l,"]"),"sigma_e_tf","alpha_0_tf","alpha_1_tf")
    if (ri_incl) {parameters <- c(parameters,"tau_0_tf")}
    if (rs_incl) {
      parameters <- c(parameters,"tau_1_tf")
      if (ri_incl) {parameters <- c(parameters,"rho_tf")}
    }
    parameters <- c(parameters,"gamma_tf","lp__")
  }
  else {
    parameters <- "h_0"
    if (xi_incl) {parameters <- c(parameters,paste0("xi[",1:m,"]"))}
    parameters <- c(parameters,paste0("beta_s[",1:p_s,"]"),paste0("beta_l[",1:p_l,"]"),"sigma_e","alpha_0","alpha_1")
    if (ri_incl) {parameters <- c(parameters,"tau_0")}
    if (rs_incl) {
      parameters <- c(parameters,"tau_1")
      if (ri_incl) {parameters <- c(parameters,"rho")}
    }
    parameters <- c(parameters,"gamma","lp__")
  }
  return(parameters)
}
get_surv_par_names <- function(m, p_s, tf=F){
  parameters <- c(paste0("xi[",1:m,"]"),"h_0",paste0("beta_s[",1:p_s,"]"),"lp__")
  if (tf){
    parameters <- c(paste0("xi[",1:m,"]"),"h_0_tf",paste0("beta_s_tf[",1:p_s,"]"),"lp__")
  }
  return(parameters)
}
get_long_par_names <- function(p_l, tf=F){
  parameters <- c(paste0("beta_l[",1:p_l,"]"),"sigma_e","tau_0","tau_1","alpha_0","alpha_1","rho","lp__")
  if (tf){
    parameters <- c(paste0("beta_l_tf[",1:p_l,"]"),"sigma_e_tf","tau_0_tf","tau_1_tf","alpha_0_tf","alpha_1_tf","rho_tf","lp__")
  }
  return(parameters)
}
get_jm_par_names_new <- function(m, p_l, p_s, tf=F, xi_incl,rs_incl,ri_incl){
  if (tf) {
    parameters <- "h_0_tf"
    if (xi_incl) {parameters <- c(parameters,paste0("xi[",1:m,"]"))}
    parameters <- c(parameters,paste0("beta_s_tf[",1:p_s,"]"),paste0("beta_l0_tf[",1:p_l,"]"),paste0("beta_l1_tf[",1:p_l,"]"),"sigma_e_tf","alpha_0_tf","alpha_1_tf")
    if (ri_incl) {parameters <- c(parameters,"tau_0_tf")}
    if (rs_incl) {
      parameters <- c(parameters,"tau_1_tf")
      if (ri_incl) {parameters <- c(parameters,"rho_tf")}
    }
    parameters <- c(parameters,"gamma_tf","lp__")
  }
  else {
    parameters <- "h_0"
    if (xi_incl) {parameters <- c(parameters,paste0("xi[",1:m,"]"))}
    parameters <- c(parameters,paste0("beta_s[",1:p_s,"]"),paste0("beta_l0[",1:p_l,"]"),paste0("beta_l1[",1:p_l,"]"),"sigma_e","alpha_0","alpha_1")
    if (ri_incl) {parameters <- c(parameters,"tau_0")}
    if (rs_incl) {
      parameters <- c(parameters,"tau_1")
      if (ri_incl) {parameters <- c(parameters,"rho")}
    }
    parameters <- c(parameters,"gamma","lp__")
  }
  return(parameters)
}
get_surv_par_names_new <- function(m, p_s, tf=F){
  parameters <- c(paste0("xi[",1:m,"]"),"h_0",paste0("beta_s[",1:p_s,"]"),"lp__")
  if (tf){
    parameters <- c(paste0("xi[",1:m,"]"),"h_0_tf",paste0("beta_s_tf[",1:p_s,"]"),"lp__")
  }
  return(parameters)
}
get_long_par_names_new <- function(p_l, tf=F){
  parameters <- c(paste0("beta_l0[",1:p_l,"]"),paste0("beta_l1[",1:p_l,"]"),"sigma_e","tau_0","tau_1","alpha_0","alpha_1","rho","lp__")
  if (tf){
    parameters <- c(paste0("beta_l0_tf[",1:p_l,"]"),paste0("beta_l1_tf[",1:p_l,"]"),"sigma_e_tf","tau_0_tf","tau_1_tf","alpha_0_tf","alpha_1_tf","rho_tf","lp__")
  }
  return(parameters)
}

# density plots
# samples - extract(out, inc_warmup=T,permuted=F),
# plots - a vector containing "chains" if the density of each chain is to be plotted separately
#         and containing "total" if the total density excluding the discarded chain is to be plotted.
# if tv = T, then
# true_values of the parameters will be added as a line on the plot
# iter the total number of iterations
# warmup - the no of warm up iterations, plot_grid - par(mfrow=plot_grid)
# discard - any chains to discard in final plot
# Plots three plots for each parameter - trace plots for all iterations and all chains on left,
# trace plots after burnin for all chains in the middle, trace plots after burnin for
# all chains except discard on the right.
my_density_plots <- function(samples,plots=c("chains","total"),
                             m=NULL, p_l=NULL, p_s=NULL,
                             parameters = NULL, tf=F, tv=F,
                             true_values = NULL,
                             iter=NULL, warmup=NULL, discard=0, plot_grid = c(3,3), nchains = 4,
                             model = "jm"){
  if (is.null(parameters)){
    if (model == "jm"){
      parameters <- get_jm_par_names(m,p_l,p_s,tf)
    }
    else if (model == "surv"){
      parameters <- get_surv_par_names(m,p_s,tf)
    }
    else if (model == "long"){
      parameters <- get_long_par_names(p_l,tf)
    }
    else stop('model must be equal to"jm", "surv" or "long", for a joint model, survival (proportional hazards) model, and longitudinal hierarchical model respectively' )
  }
  if (tv){
    if (is.null(true_values)){
      true_values <- apply(array(parameters), 1, FUN=get)
    }
    if (is.null(iter)){
      iter = dim(samples)[1]
    }
    if (is.null(warmup)){
      warmup = iter/2
    }
    chains <- 1:nchains
    keep_chains <- chains[which(!is.element(chains,discard))]
    
    par(mfrow=plot_grid)
    if (is.element("chains",plots)){
      param_count <- 0
      for (param in parameters){
        param_count <- param_count+1
        tv <- true_values[param_count]
        d <- numeric(nchains)
        xmin <- numeric(nchains)
        xmax <- numeric(nchains)
        for(i in 1:nchains){
          dens <- density(samples[(warmup+1):(iter),i,param])
          d[i] <- max(dens$y)
          xmin[i] <- min(dens$x,tv, na.rm=T)
          xmax[i] <- max(dens$x,tv, na.rm=T)
        }
        plot(NA,ylim = c(0,max(d)), xlim = c(min(xmin),max(xmax)), main = param)
        for (i in 1:nchains){
          lines(density(samples[(warmup+1):(iter),i,param]),col=i)
        }
        abline(v=tv,lty=2)
      }
    }
    par(mfrow=plot_grid)
    if (is.element("total",plots)){
      param_count <- 0
      for (param in parameters){
        param_count <- param_count+1
        tv <- true_values[param_count]
        dens <- density(samples[(warmup+1):(iter),keep_chains,param])
        d <- max(dens$y)
        xmin <- min(dens$x,tv, na.rm=T)
        xmax <- max(dens$x,tv, na.rm=T)
        plot(NA,ylim = c(0,max(d)), xlim = c(min(xmin),max(xmax)), main = param)
        lines(density(as.numeric(samples[(warmup+1):(iter),keep_chains,param])))
        abline(v=tv,lty=2)
      }
    }
  }
  else if (!tv){
    if (is.null(iter)){
      iter = dim(samples)[1]
    }
    if (is.null(warmup)){
      warmup = iter/2
    }
    chains <- 1:nchains
    keep_chains <- chains[which(!is.element(chains,discard))]
    
    par(mfrow=plot_grid)
    if (is.element("chains",plots)){
      param_count <- 0
      for (param in parameters){
        param_count <- param_count+1
        d <- numeric(nchains)
        xmin <- numeric(nchains)
        xmax <- numeric(nchains)
        for(i in 1:nchains){
          dens <- density(samples[(warmup+1):(iter),i,param])
          d[i] <- max(dens$y)
          xmin[i] <- min(dens$x, na.rm=T)
          xmax[i] <- max(dens$x, na.rm=T)
        }
        plot(NA,ylim = c(0,max(d)), xlim = c(min(xmin),max(xmax)), main = param)
        for (i in 1:nchains){
          lines(density(samples[(warmup+1):(iter),i,param]),col=i)
        }
      }
    }
    par(mfrow=plot_grid)
    if (is.element("total",plots)){
      param_count <- 0
      for (param in parameters){
        param_count <- param_count+1
        dens <- density(samples[(warmup+1):(iter),keep_chains,param])
        d <- max(dens$y)
        xmin <- min(dens$x, na.rm=T)
        xmax <- max(dens$x, na.rm=T)
        plot(NA,ylim = c(0,max(d)), xlim = c(min(xmin),max(xmax)), main = param)
        lines(density(as.numeric(samples[(warmup+1):(iter),keep_chains,param])))
      }
    }
  }
  par(mfrow=c(1,1))
}

my_spline_plots <- function(stan_output,n_sub_sample,n_spline = m,zoom.out = 3,
                            B=NULL, xi=NULL, h_0=NULL, points, real = T){
  par(mfrow=c(1,1))
  sub_samples <- as.matrix(stan_output)
  sub_samples <- sub_samples[sample(1:nrow(sub_samples),n_sub_sample),,drop=FALSE]
  sub_samples <- sub_samples[,c(paste0("xi[",1:n_spline,"]"),"h_0"),drop=FALSE]
  
  if (!real){
    true_values_spline <- B%*%xi+h_0
    min_tv <- min(true_values_spline)
    max_tv <- max(true_values_spline)
    
    plot(NA,xlim=c(min(points),max(points)), ylim = zoom.out*((max_tv - min_tv)*c(-1,1)/2)+(max_tv + min_tv)/2, type = "l", ylab="log baseline hazard",xlab="time")
    for (i in 1:n_sub_sample){
      lines((B%*%as.numeric(sub_samples[i,1:n_spline]) + sub_samples[i,n_spline+1]) ~ points, col = "grey", lty=1,lwd=0.5)
    }
    lines((B%*%xi+h_0)~points, type = "l", col="red",lwd=1)
  }
  else if (real){
    bh_values <- matrix(NA,nrow=n_sub_sample,ncol=length(points))
    for (i in 1:n_sub_sample){
      bh_values[i,] <- (B%*%as.numeric(sub_samples[i,1:n_spline]) + sub_samples[i,n_spline+1])
    }
    min_tv <- min(bh_values)
    max_tv <- max(bh_values)
    zoom.out <- 1.1
    plot(NA,xlim=c(min(points),max(points)), ylim = zoom.out*((max_tv - min_tv)*c(-1,1)/2)+(max_tv + min_tv)/2, type = "l", ylab="log baseline hazard",xlab="time")
    for (i in 1:n_sub_sample){
      lines(bh_values[i,] ~ points, col = "grey", lty=1,lwd=0.5)
    }
  }
  
}

my_summary <- function(samples_mat){
  quants <- c(0.025,0.1,0.25,0.5,0.75,0.9,0.975)
  p <- ncol(samples_mat)
  sum_mat <- matrix(NA,p,2+length(quants))
  rownames(sum_mat) <- colnames(samples_mat)
  sum_mat[,1] <- apply(samples_mat,2,mean)
  sum_mat[,2] <- apply(samples_mat,2,sd)
  sum_mat[,2+(1:length(quants))] <- t(apply(samples_mat,2,quantile, quants))
  colnames(sum_mat) <- c("mean","sd",quants)
  sum_mat <- round(sum_mat,3)
  return(sum_mat)
}

