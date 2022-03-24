library('Rcpp')
library('tidyr')
library('dplyr')
library('loo')
library('bridgesampling')
library('rstan')
library('rstanarm')
library('bayestestR')
library('insight')
library('ggplot2')
library('ggmcmc')
library('gridExtra')

# ---------------------------------------------------------

PATHTODATA <- 

# ---------------------------------------------------------

options( mc.cores = 4 ) 
rstan_options( auto_write = TRUE)

# ---------------------------------------------------------

stan_power_law                        <- stan_model( file = "/power_law.stan" )

# ---------------------------------------------------------

raw_data                              <- read.csv( PATHTODATA, header=FALSE)
names(raw_data)                       <- c( 'spk','RT')

raw_data$RT                           <- raw_data$RT / 1000

stan_data                             <- list( 'RT' = raw_data$RT, 
                                               'spk' = raw_data$spk, 
                                               'stim' = rep(1,nrow(raw_data)),
                                               'n' = nrow(raw_data),
                                               'n_stim' = 1,
                                               'PRIOR_ONLY' = 0 )

# ---------------------------------------------------------

model_power_law                       <- sampling( stan_power_law, 
                                                   data = stan_data, 
                                                   iter = 5000, 
                                                   warmup = 2000,
                                                   thin = 2,
                                                   seed = 2022,
                                                   cores = 4 )

# ---------------------------------------------------------

matrix_power_law                      <- as.matrix( model_power_law )

# ---------------------------------------------------------

ypreds_power_law                      <- matrix_power_law[,(5:(nrow(raw_data) + 4))]
lowers_power_law                      <- apply(ypreds_power_law, 2, 
                                               function(x) quantile(x, probs = 0.055))
uppers_power_law                      <- apply(ypreds_power_law, 2, 
                                               function(x) quantile(x, probs = 0.945))
medians_power_law                     <- apply(ypreds_power_law, 2, 
                                               function(x) quantile(x, probs = 0.5))

# ---------------------------------------------------------

posterior_power_law                   <- extract(model_power_law)
alphas_power_law                      <- posterior_power_law$alpha
betas_power_law                       <- posterior_power_law$beta
rs_power_law                          <- posterior_power_law$r

estimates_power_law                   <- matrix( nrow = nrow(raw_data), 
                                                 ncol = length(alphas_power_law) )

for(i in 1:length( alphas_power_law )){
  estimates_power_law[,i]             <- ( alphas_power_law[i] + 
                                             betas_power_law[i] * 
                                             raw_data$spk^(-rs_power_law[i]) )
}

ypreds_power_law_estimate             <- apply( estimates_power_law, 1, 
                                                median )
lowers_power_law_estimate             <- apply( estimates_power_law, 1, 
                                                function(x) quantile(x, probs = 0.055) )
uppers_power_law_estimate             <- apply( estimates_power_law, 1, 
                                                function(x) quantile(x, probs = 0.945) )

# ---------------------------------------------------------

bridge_power_law                      <- bridge_sampler( model_linear_regression )

# ---------------------------------------------------------

loo_power_law                      <- loo( model_linear_regression )

# ---------------------------------------------------------

plot_data_power_law                             <- cbind(raw_data, 
                                               lowers_power_law, 
                                               uppers_power_law,
                                               ypreds_power_law,
                                               medians_power_law,
                                               lowers_power_law_estimate, 
                                               uppers_power_law_estimate,
                                               ypreds_power_law_estimate)
