library(tidyverse)
library(brms)
library(cmdstanr)
# install_cmdstan()  # This downloads and builds CmdStan locally
# Set backend explicitly
options(brms.backend = "cmdstanr")

# read in data and set variables
model_df <- readr::read_csv('data/site_model_df.csv')

combos1_df <- readr::read_csv('data/site-model-combos_1var.csv')
combos2_df <- readr::read_csv('data/site-model-combos_2var.csv')

wq_max_2model_logistic <- model_df
response = 'FIB_gt35_01'

# SITE LEVEL MODEL

# Run 1 variable models -----------------------------------------------

run_fib_model_brm <- function(model_id, response = 'FIB_gt35_01'){
  # Response + predictors
  response <- response
  model_id_row <- which(combos1_df$model_id==model_id)
  fixed_vars <- combos1_df[['Var1']][model_id_row] 
  
  # Reformulate
  form <- stats::reformulate(fixed_vars, response)
  
  print(glue::glue('Running model {model_id}...'))
  # Fit the model
  brms_model <- brm(
    formula = form,
    data = wq_max_2model_logistic,
    family = bernoulli(),  # binomial with 1 trial = Bernoulli
    prior =   set_prior("normal(0, 1000)", class = "b"),  # flat for fixed effects
    chains = 4,            # number of MCMC chains
    cores = 4,             # for parallel computation
    iter = 2000,           # number of iterations per chain
    backend = "cmdstanr",  # optional if you've set the global option
    save_pars = save_pars(all = TRUE),
    seed = 123             # for reproducibility
  )
  loo_fit <- loo(brms_model, moment_match = TRUE, reloo = TRUE)
  model_r2 <- bayes_R2(brms_model)
  
  # Effect size 
  brms_model_summary <- summary(brms_model)
  intercept_row_id <- which(row.names(brms_model_summary$fixed)=="Intercept")
  var1_est <- brms_model_summary$fixed[-intercept_row_id,'Estimate']
  var1_lo <- brms_model_summary$fixed[-intercept_row_id,'l-95% CI']
  var1_hi <- brms_model_summary$fixed[-intercept_row_id,'u-95% CI']
  includes0 <- var1_lo * var1_hi < 0
  model_id_var1 <- combos1_df$Var1[which(combos1_df$model_id == model_id)]
  
  # save model itself
  fs::dir_create('output/01_SiteLevelModels/01_1varModels_v1/models')
  model_file <- glue::glue('output/01_SiteLevelModels/01_1varModels_v1/models/{model_id}.rds')
  saveRDS(brms_model, file = model_file)
  
  model_stats_df <- data.frame(model_id = model_id,
                               model_var1 = fixed_vars[1],
                               model_elpd = loo_fit$estimates[which(row.names(loo_fit$estimates)=='elpd_loo'),'Estimate'],
                               model_elpdSE = loo_fit$estimates[which(row.names(loo_fit$estimates)=='elpd_loo'),'SE'],
                               model_looic = loo_fit$estimates[which(row.names(loo_fit$estimates)=='looic'),'Estimate'],
                               model_looicSE = loo_fit$estimates[which(row.names(loo_fit$estimates)=='looic'),'SE'],
                               model_ploo = loo_fit$estimates[which(row.names(loo_fit$estimates)=='p_loo'),'Estimate'],
                               model_plooSE = loo_fit$estimates[which(row.names(loo_fit$estimates)=='p_loo'),'SE'],
                               model_r2 = model_r2[1,'Estimate'], model_r2stdev = model_r2[1,'Est.Error'],
                               model_r2q02 = model_r2[1,'Q2.5'], model_r2q97 = model_r2[1,'Q97.5'],
                               model_id_var1 = model_id_var1,
                               var1 = rownames(fixef(brms_model))[-intercept_row_id],
                               var1_est = var1_est,
                               var1_lo = var1_lo,
                               var1_hi = var1_hi,
                               var1_includes0 = includes0
  )
  # save csv of model stats
  # Use LOO-CV (elpd_loo) for model selection, and report Bayes R² and posterior estimates for interpretability.
  fs::dir_create('output/01_SiteLevelModels/01_1varModels_v1/model_stats')
  stats_file <- glue::glue('output/01_SiteLevelModels/01_1varModels_v1/model_stats/{model_id}_stats.csv')
  readr::write_csv(model_stats_df, file = stats_file)
}

# Run the model over all combinations of variables

combos1_df$model_id %>% purrr::map(~run_fib_model_brm(.x))

# check progress for how many stats files there are
# length(fs::dir_ls('output/01_SiteLevelModels/01_1varModels_v1/model_stats'))

# 2 variable models


# Run 2 variable models -----------------------------------------------


run_fib_model_brm2var <- function(model_id, response = 'FIB_gt35_01'){
  # Response + predictors
  response <- response
  model_id_row <- which(combos2_df$model_id==model_id)
  fixed_vars <- c(combos2_df[['var1']][model_id_row],
                  combos2_df[['var2']][model_id_row])
  
  # Reformulate
  form <- stats::reformulate(fixed_vars, response)
  
  print(glue::glue('Running model {model_id}...'))
  # Fit the model
  brms_model <- brm(
    formula = form,
    data = wq_max_2model_logistic,
    family = bernoulli(),  # binomial with 1 trial = Bernoulli
    prior =   set_prior("normal(0, 1000)", class = "b"),  # flat for fixed effects
    chains = 4,            # number of MCMC chains
    cores = 4,             # for parallel computation
    iter = 2000,           # number of iterations per chain
    backend = "cmdstanr",  # optional if you've set the global option
    save_pars = save_pars(all = TRUE),
    seed = 123             # for reproducibility
  )
  loo_fit <- loo(brms_model, moment_match = TRUE, reloo = TRUE)
  model_r2 <- bayes_R2(brms_model)
  
  # Effect size 
  brms_model_summary <- summary(brms_model)
  intercept_row_id <- which(row.names(brms_model_summary$fixed)=="Intercept")
  var1_row_id <- which(row.names(brms_model_summary$fixed)==fixed_vars[1])
  var2_row_id <- which(row.names(brms_model_summary$fixed)==fixed_vars[2])
  var1_est <- brms_model_summary$fixed[var1_row_id,'Estimate']
  var1_lo <- brms_model_summary$fixed[var1_row_id,'l-95% CI']
  var1_hi <- brms_model_summary$fixed[var1_row_id,'u-95% CI']
  var1_includes0 <- var1_lo * var1_hi < 0
  var2_est <- brms_model_summary$fixed[var2_row_id,'Estimate']
  var2_lo <- brms_model_summary$fixed[var2_row_id,'l-95% CI']
  var2_hi <- brms_model_summary$fixed[var2_row_id,'u-95% CI']
  var2_includes0 <- var2_lo * var2_hi < 0
  # save model itself
  fs::dir_create('output/01_SiteLevelModels/02_2varModels_v1/models/')
  model_file <- glue::glue('output/01_SiteLevelModels/02_2varModels_v1/models/{model_id}.rds')
  saveRDS(brms_model, file = model_file)
  
  model_stats_df <- data.frame(model_id = model_id,
                               model_var1 = fixed_vars[1],
                               model_var2 = fixed_vars[2],
                               model_elpd = loo_fit$estimates[which(row.names(loo_fit$estimates)=='elpd_loo'),'Estimate'],
                               model_elpdSE = loo_fit$estimates[which(row.names(loo_fit$estimates)=='elpd_loo'),'SE'],
                               model_looic = loo_fit$estimates[which(row.names(loo_fit$estimates)=='looic'),'Estimate'],
                               model_looicSE = loo_fit$estimates[which(row.names(loo_fit$estimates)=='looic'),'SE'],
                               model_ploo = loo_fit$estimates[which(row.names(loo_fit$estimates)=='p_loo'),'Estimate'],
                               model_plooSE = loo_fit$estimates[which(row.names(loo_fit$estimates)=='p_loo'),'SE'],
                               model_r2 = model_r2[1,'Estimate'], 
                               model_r2stdev = model_r2[1,'Est.Error'],
                               model_r2q02 = model_r2[1,'Q2.5'], 
                               model_r2q97 = model_r2[1,'Q97.5'],
                               var1_est = var1_est,
                               var1_lo = var1_lo,
                               var1_hi = var1_hi,
                               var1_includes0 = var1_includes0,
                               var2_est = var2_est,
                               var2_lo = var2_lo,
                               var2_hi = var2_hi,
                               var2_includes0 = var2_includes0
  )
  # save csv of model stats
  # Use LOO-CV (elpd_loo) for model selection, and report Bayes R² and posterior estimates for interpretability.
  fs::dir_create('output/01_SiteLevelModels/02_2varModels_v1/model_stats/')
  stats_file <- glue::glue('output/01_SiteLevelModels/02_2varModels_v1/model_stats/{model_id}_stats.csv')
  readr::write_csv(model_stats_df, file = stats_file)
}



combos2_df$model_id %>% purrr::map(~run_fib_model_brm2var(.x))

# check progress for how many stats files there are
length(fs::dir_ls('output/01_SiteLevelModels/02_2varModels_v1/model_stats'))



