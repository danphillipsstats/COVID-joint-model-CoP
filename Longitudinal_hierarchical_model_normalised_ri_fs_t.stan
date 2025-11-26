// A longitudinal hierarchical model with two mixed effects
// Specific to where there is one grouping, but multiple mixed effects parameters
data {
	int <lower=0> 		n_l; 			// number of observations
	int <lower=0> 		p_l;		 	// number of fixed effect parameters
	int <lower=0> 		n; 			// number of individuals
	int<lower=1, upper=n> 	ll[n_l];		// Group indicators
	matrix[n, p_l] 		X_l;	 		// Group- level predictor matrix for fixed effects
	vector[n_l] 		log_y;			// outcome vector
	vector[n_l] 		t_l;			// times of observations
	real  t_center;  // Shift for t under new parametrisation.
}

transformed data{
	matrix[n,p_l] X_l_tf;
	vector[p_l] X_l_col_means;
	vector[p_l] X_l_col_sd;
	vector[n_l] log_y_tf;
	real log_y_mean;
	real log_y_sd;
	real t_sd;
	vector[n_l] t_l_tf;
	
	for (j in 1:p_l){
		X_l_col_means[j] = mean(X_l[,j]);
		X_l_col_sd[j] = sd(X_l[,j]);
		X_l_tf[,j] = (X_l[,j] - X_l_col_means[j])/X_l_col_sd[j];
	}
	log_y_mean = mean(log_y);
	log_y_sd = sd(log_y);
	log_y_tf = (log_y - log_y_mean)/log_y_sd;
	
	t_sd = sd(t_l);
	t_l_tf = (t_l - t_center)/t_sd;
}

parameters {
	vector[p_l] 		beta_l_tf;			// coefficients for fixed effects
	real<lower=0> 		sigma_e_tf;		// error scale
	real<lower=0> 		tau_0_tf;			// random intercept scale
	real 			alpha_0_tf;		// Mean random intercept
	real 			alpha_1_tf;		// Mean random slope
	vector[n] 		eta_0_tf;			// Random intercepts
}
transformed parameters{
	vector[n] 		a_tf;		// predicted longitudinal outcome value at t'=0
	
	a_tf = X_l_tf*beta_l_tf + alpha_0_tf + tau_0_tf*eta_0_tf;
}

model {
	vector[n_l] 		log_mu_tf;			// Mean for data point
	
	// Priors
	tau_0_tf ~ student_t(2,0,1);
	sigma_e_tf ~ student_t(2,0,1);
	beta_l_tf ~ normal(0,2);
	
	// Random effects part
	alpha_0_tf ~ normal(0,2);
	alpha_1_tf ~ normal(0,2);
	eta_0_tf ~ std_normal();
	
	// Model
	log_mu_tf = a_tf[ll] + t_l_tf * alpha_1_tf;
	log_y_tf ~ normal(log_mu_tf,sigma_e_tf);
}

generated quantities {
	
	// Longitudinal component
	vector[p_l] beta_l;
	real<lower=0> sigma_e;
	real<lower=0> tau_0;
	real alpha_0;
	real alpha_1;
	vector[n] a;
	//vector[n] eta_0; We don't care about the eta values I argue, so I won't bother transforming them

	// Longitudinal component
	beta_l = log_y_sd * beta_l_tf ./ X_l_col_sd;
	sigma_e = log_y_sd * sigma_e_tf;
	tau_0 = log_y_sd*tau_0_tf;
	alpha_0 = log_y_sd * (alpha_0_tf - dot_product(beta_l_tf, X_l_col_means ./ X_l_col_sd) - alpha_1_tf * t_center / t_sd) + log_y_mean;
	alpha_1 = alpha_1_tf * log_y_sd / t_sd;
	a = log_y_sd * (a_tf - alpha_1_tf * t_center / t_sd) + log_y_mean;
}
