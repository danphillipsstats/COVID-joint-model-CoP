// A longitudinal hierarchical model for log antibody decay.
// Random intercept and random slope. Random slope is constrained to be negative.
// The model allows for covariate effects.
// Observed log antibodies is referred to as log_y
// Data is normalised before the model is fit to improve mixing
// then parameters are transformed back to the unnormalised versions to improve interpretability.
// t-distributed random error assumed, and Gaussian random effects. Weakly informative priors.
data {
	int <lower=0> 		n_l; 			// number of observations
	int <lower=0> 		p_l;		 	// number of fixed effect parameters
	int <lower=0> 		n; 			// number of individuals
	int<lower=1, upper=n> 	ll[n_l];		// Group indicators
	matrix[n, p_l] 		X_l;	 		// Predictor matrix for fixed effects
	vector[n_l] 		log_y;			// Observed log antibodies vector
	vector[n_l] 		t_l;			// times of observations
	real                    t_center;               // Shift for t so t_center is the new time 0
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
	t_l_tf = (t_l-t_center)/t_sd;
}

parameters {
	vector[p_l] 		beta_l0_tf;			// covariate effects on intercept
	vector[p_l] 		beta_l1_tf;			// covariate effects on slope
	real<lower=0> 		sigma_e_tf;		// error scale
	real<lower=0> 		tau_0_tf;			// random intercept scale
	real<lower=0> 		tau_1_tf;			// random slope scale
	real 			alpha_0_tf;		// Population intercept
	real 			alpha_1_tf;		// Population slope
	real<lower=-1,upper=1>  rho_tf;			// Correlation between random intercept and random slopes
	vector[n] 		eta_0_tf;			// Individual random intercept effect
	vector[n] 		eta_1_tf;			// Individual random slope effect
}

transformed parameters {
	vector[n] 		a_0_tf;		// predicted longitudinal outcome value at t=0
	vector[n] 		a_1_tf;		// predicted longitudinal slope
	vector[n_l] 		log_mu_tf;			// Mean for data point
	
	a_0_tf = X_l_tf*beta_l0_tf + alpha_0_tf + tau_0_tf*eta_0_tf;
	a_1_tf = -exp(X_l_tf*beta_l1_tf + alpha_1_tf + tau_1_tf*(rho_tf*eta_0_tf + sqrt(1-rho_tf^2)*eta_1_tf)); // Using Cholesky decomposition
	
	log_mu_tf = a_0_tf[ll] + t_l_tf .* a_1_tf[ll];
}

model {
	
	// Priors
	tau_0_tf ~ student_t(2,0,1);
	tau_1_tf ~ student_t(2,0,1);
	sigma_e_tf ~ student_t(2,0,1);
	beta_l0_tf ~ normal(0,2);
	beta_l1_tf ~ normal(0,2);
	
	// Random effects part
	alpha_0_tf ~ normal(0,2);
	alpha_1_tf ~ normal(0,2);
	eta_0_tf ~ std_normal();
	eta_1_tf ~ std_normal();
	rho_tf ~ uniform(-1,1);
	
	// Model
	log_y_tf ~ student_t(4,log_mu_tf,sigma_e_tf);
}

generated quantities {
	
	// Random intercept and slope
	vector[n] a_0;
	vector[n] a_1;
	// Fixed effect parameters
	real 			alpha_0;		// Population intercept
	real 			alpha_1;		// Population slope
	vector[p_l] 		beta_l0;			// coefficients for fixed effects
	vector[p_l] 		beta_l1;			// coefficients for fixed effects
	real<lower=0> 		sigma_e;		// error scale
	real<lower=0> 		tau_0;			// covariate effects on intercept
	real<lower=0> 		tau_1;			// covariate effects on slope
	real<lower=-1,upper=1>  rho;			// Correlation between random intercept and random slopes
	vector[n_l]     resid;        // Subject-specific longitudinal residuals
	
	// Random intercept and slope
	a_0 = log_y_sd * a_0_tf + log_y_mean; 
	a_1 = a_1_tf * log_y_sd / t_sd;
	// Fixed effect parameters
	alpha_0 = log_y_sd * (alpha_0_tf - dot_product(beta_l0_tf, X_l_col_means ./ X_l_col_sd)) + log_y_mean; // the intercept at time t_center
	alpha_1 = alpha_1_tf + log(log_y_sd/t_sd) - dot_product(beta_l1_tf, X_l_col_means ./ X_l_col_sd);
	beta_l0 = log_y_sd * beta_l0_tf ./ X_l_col_sd; // the covariate effect on value at time t_center
	beta_l1 = beta_l1_tf ./ X_l_col_sd;
	tau_0 = log_y_sd * tau_0_tf;
	tau_1 = tau_1_tf;
	rho = rho_tf;
	sigma_e = log_y_sd * sigma_e_tf;
	resid = (log_y_tf-log_mu_tf)/sigma_e_tf;
}
