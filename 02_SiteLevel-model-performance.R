# read in and combine output from 1 variable models
library(data.table)
library(brms)
library(tidyverse)

# # read in and combine output from 1 variable models

models1var_files <- fs::dir_ls('output/01_SiteLevelModels/01_1varModels_v1/model_stats/')
combined <- rbindlist(lapply(models1var_files, fread))

get_looR2 <- function(model_id){
  model_file <- glue::glue('output/01_SiteLevelModels/01_1varModels_v1/models/{model_id}.rds')
  brms_model <- readRDS(model_file)
  loor2_df <- loo_R2(brms_model)
  
  loo_r2 <- loor2_df[,'Estimate']
  loo_r2se <- loor2_df[,'Est.Error']
  loo_r2lo <- loor2_df[,'Q2.5']
  loo_r2hi <- loor2_df[,'Q97.5']
  
  return_df <- data.frame(model_id = model_id,
                          loo_r2 = loo_r2,
                          loo_r2se = loo_r2se,
                          loo_r2lo = loo_r2lo,
                          loo_r2hi = loo_r2hi)
  return(return_df)
}

models1var_loor2_df <- combos1_df_subset$model_id %>% purrr::map(~get_looR2(.x))
models1var_loor2_df <- models1var_loor2_df %>% bind_rows()

combined_join_1var <- combined %>% dplyr::left_join(models1var_loor2_df, by = 'model_id')

# # read in and combine output from 2 variable models

models2var_files <- fs::dir_ls('output/01_SiteLevelModels/02_2varModels_v1/model_stats/')
combined <- rbindlist(lapply(models2var_files, fread))

get_looR2 <- function(model_id){
  model_file <- glue::glue('output/01_SiteLevelModels/02_2varModels_v1/models/{model_id}.rds')
  brms_model <- readRDS(model_file)
  loor2_df <- loo_R2(brms_model)
  
  loo_r2 <- loor2_df[,'Estimate']
  loo_r2se <- loor2_df[,'Est.Error']
  loo_r2lo <- loor2_df[,'Q2.5']
  loo_r2hi <- loor2_df[,'Q97.5']
  
  return_df <- data.frame(model_id = model_id,
                          loo_r2 = loo_r2,
                          loo_r2se = loo_r2se,
                          loo_r2lo = loo_r2lo,
                          loo_r2hi = loo_r2hi)
  return(return_df)
}

models2var_loor2_df <- combos2_df$model_id %>% purrr::map_df(~get_looR2(.x))

combined_join <- combined %>% dplyr::left_join(models2var_loor2_df, by = 'model_id')
