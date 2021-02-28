args <- commandArgs(TRUE)
eval(parse(text = args))
t1 <- Sys.time()

sampsize <- 450
# USE current_case_ids AND SAMPSIZE AS ARGUMENT NAME
## Analyses simulation data settings

library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(EFAtools)
library(bayestestR)
library(rstanarm)

# setwd("/Users/msteiner/Documents/Doktorat/efacomp_OSF/r/scicore/simulation_settings/")


model_control <- expand.grid(
  case = names(population_models$loadings)[1:27],
  cors = names(population_models$phis_3),
  N = c(sampsize), # repeat with 180
  stringsAsFactors = FALSE)
model_control$case_ids <- 1:nrow(model_control)

settings <- expand.grid(
  comm_meth = c("unity", "mac", "smc"),  # initial communality estimation methods
  criterion_type = c("max_individual", "sums"), # citerion types
  abs_eigen = c(TRUE, FALSE), # absolute eigenvalues yes or no
  conv_crit = c(1e-3, 1e-6), # convergence criterion
  var_type = c("svd", "kaiser"), # varimax type
  p_type = c("unnorm", "norm"),
  k = c(3, 4),
  stringsAsFactors = FALSE)

# define k according to earlier simulation studies
settings$k[settings$k == 3 & settings$p_type == "norm"] <- 2

settings$setting_id <- 1:nrow(settings)


# read in data
recovery <- readRDS(paste0(path, "/recovery_", sampsize, ".RDS"))

# to decrease size of recovery object
recovery <- recovery %>%
  filter(case_ids %in% current_case_ids)


# cohen's d function
cohens_d <- function(x, y) {
  
  n1 <- length(x)
  n2 <- length(y)
  m1 <- mean(x)
  m2 <- mean(y)
  s1 <- var(x)
  s2 <- var(y)
  
  # calculate pooled standard deviation with pooled variances
  s <- sqrt(((n1 - 1)*s1 + (n2 - 1)*s2) / (n1 + n2 -2))
  
  # cohen's d
  d <- (m1 - m2) / s
  
  # return cohen's d
  d
}

### deviance ===================================================================

iter_max <- nrow(settings)
try_steps <- seq(10, 190, 10)

results <- list()
for(i in current_case_ids){
  dat_ids_excl <- recovery %>%
    filter(case == model_control$case[i] & cors == model_control$cors[i] &
             N == model_control$N[i] & !is.na(error)) %>%
    pull(dat_id) %>%
    unique()

  N_datasets <- 1000 - length(dat_ids_excl)

  if (N_datasets > 10) {


    # only keep relevant data
    temp_dat <- recovery %>%
      filter(case == model_control$case[i] & cors == model_control$cors[i] &
               N == model_control$N[i] & !(dat_id %in% dat_ids_excl)) %>%
      select(-starts_with("error")) %>%
      drop_na()

    # get ordered setting ids
    sett_ids <- temp_dat %>%
      group_by(setting_id) %>%
      summarise(
        g = mean(g, na.rm = TRUE)
      ) %>%
      arrange(g) %>%
      pull(setting_id)


    for (step_i in try_steps) {

      mod <- stan_glm(g ~ setting_id,
                      data = temp_dat %>%
                        filter(setting_id %in% c(sett_ids[c(1, step_i)])) %>%
                        mutate(setting_id = factor(setting_id,
                                                   levels = sett_ids[c(1, step_i)])),
                      chains = 4, iter = 10000, family = gaussian(),
                      prior_intercept = normal(location = 0, scale = 10),
                      prior = normal(location = 0, scale = 2.5))

      post_delta_m <- describe_posterior(mod, ci = .95, rope_ci = 1)[2,]
      evidence <- if (sign(post_delta_m$CI_low) == sign(post_delta_m$CI_high) &
                      post_delta_m$ROPE_Percentage < .025) {
        break
      }

    }

    if (step_i == 10) {
      iter <- 2
    } else {
      iter <- step_i - 9
    }


    while(TRUE) {

      mod <- stan_glm(g ~ setting_id,
                      data = temp_dat %>%
                        select(g, setting_id) %>%
                        filter(setting_id %in% c(sett_ids[c(1, iter)])) %>%
                        mutate(setting_id = factor(setting_id,
                                                   levels = sett_ids[c(1, iter)])),
                      chains = 4, iter = 10000, family = gaussian(),
                      prior_intercept = normal(location = 0, scale = 10),
                      prior = normal(location = 0, scale = 2.5))

      post_delta_m <- describe_posterior(mod, ci = .95, rope_ci = 1)[2,]
      evidence <- if (sign(post_delta_m$CI_low) == sign(post_delta_m$CI_high) &
                      post_delta_m$ROPE_Percentage < .025) {
        "Different"
      } else if (post_delta_m$ROPE_Percentage > .95) {
        "Equal"
      } else {
        "Inconclusive"
      }
      es <- cohens_d(temp_dat$g[temp_dat$setting_id == sett_ids[iter]],
                     temp_dat$g[temp_dat$setting_id == sett_ids[1]])

      if (iter == iter_max || evidence == "Different") {

        break

      }

      iter <- iter + 1

    }

    best_ids <- paste(sett_ids[1:(iter - 1)], collapse = ";")

    # save output
    results[[i]] <- data.frame(
      case = model_control$case[i],
      cors = model_control$cors[i],
      N = model_control$N[i],
      case_id = model_control$case_ids[i],
      m_best = mean(temp_dat$g[temp_dat$setting_id == sett_ids[1]]),
      m_diff = mean(temp_dat$g[temp_dat$setting_id == sett_ids[iter]]),
      post_delta_m = post_delta_m$Median,
      post_delta_m_low = post_delta_m$CI_low,
      post_delta_m_high = post_delta_m$CI_high,
      evidence = evidence,
      es = es,
      best_ids = best_ids,
      iter = iter,
      N_datasets = N_datasets
    )

  } else {

    # save output
    results[[i]] <- data.frame(
      case = model_control$case[i],
      cors = model_control$cors[i],
      N = model_control$N[i],
      case_id = model_control$case_ids[i],
      m_best = NA,
      m_diff = NA,
      post_delta_m = NA,
      post_delta_m_low = NA,
      post_delta_m_high = NA,
      evidence = NA,
      es = NA,
      best_ids = NA,
      iter = NA,
      N_datasets = N_datasets
    )

  }
}

res_deviance <- do.call(rbind, results)
saveRDS(res_deviance, paste0("output/ad/results_deviance_settings_",
        sampsize, "_", current_case_ids, ".RDS"))

### deviance datasets with negative eigenvalues ================================
# 
iter_max <- nrow(settings)

results <- list()
for(i in current_case_ids){
  dat_ids_excl <- recovery %>%
    filter(case == model_control$case[i] & cors == model_control$cors[i] &
             N == model_control$N[i] & !is.na(error)) %>%
    pull(dat_id) %>%
    unique()

  N_datasets <- length(dat_ids_excl)

  if (N_datasets > 10) {


    # only keep relevant data
    temp_dat <- recovery %>%
      filter(case == model_control$case[i] & cors == model_control$cors[i] &
               N == model_control$N[i] & (dat_id %in% dat_ids_excl)) %>%
      select(-starts_with("error")) %>%
      drop_na()

    # get ordered setting ids
    sett_ids <- temp_dat %>%
      group_by(setting_id) %>%
      summarise(
        g = mean(g, na.rm = TRUE)
      ) %>%
      arrange(g) %>%
      pull(setting_id)
    iter_max = length(sett_ids)
    try_steps <- seq(10, iter_max, 10)

    for (step_i in try_steps) {

      mod <- stan_glm(g ~ setting_id,
                      data = temp_dat %>%
                        filter(setting_id %in% c(sett_ids[c(1, step_i)])) %>%
                        mutate(setting_id = factor(setting_id,
                                                   levels = sett_ids[c(1, step_i)])),
                      chains = 4, iter = 10000, family = gaussian(),
                      prior_intercept = normal(location = 0, scale = 10),
                      prior = normal(location = 0, scale = 2.5))

      post_delta_m <- describe_posterior(mod, ci = .95, rope_ci = 1)[2,]
      evidence <- if (sign(post_delta_m$CI_low) == sign(post_delta_m$CI_high) &
                      post_delta_m$ROPE_Percentage < .025) {
        break
      }

    }

    if (step_i == 10) {
      iter <- 2
    } else {
      iter <- step_i - 9
    }


    while(TRUE) {

      mod <- stan_glm(g ~ setting_id,
                      data = temp_dat %>%
                        select(g, setting_id) %>%
                        filter(setting_id %in% c(sett_ids[c(1, iter)])) %>%
                        mutate(setting_id = factor(setting_id,
                                                   levels = sett_ids[c(1, iter)])),
                      chains = 4, iter = 10000, family = gaussian(),
                      prior_intercept = normal(location = 0, scale = 10),
                      prior = normal(location = 0, scale = 2.5))

      post_delta_m <- describe_posterior(mod, ci = .95, rope_ci = 1)[2,]
      evidence <- if (sign(post_delta_m$CI_low) == sign(post_delta_m$CI_high) &
                      post_delta_m$ROPE_Percentage < .025) {
        "Different"
      } else if (post_delta_m$ROPE_Percentage > .95) {
        "Equal"
      } else {
        "Inconclusive"
      }
      es <- cohens_d(temp_dat$g[temp_dat$setting_id == sett_ids[iter]],
                     temp_dat$g[temp_dat$setting_id == sett_ids[1]])

      if (iter == iter_max || evidence == "Different") {

        break

      }

      iter <- iter + 1

    }

    best_ids <- paste(sett_ids[1:(iter - 1)], collapse = ";")

    # save output
    results[[i]] <- data.frame(
      case = model_control$case[i],
      cors = model_control$cors[i],
      N = model_control$N[i],
      case_id = model_control$case_ids[i],
      m_best = mean(temp_dat$g[temp_dat$setting_id == sett_ids[1]]),
      m_diff = mean(temp_dat$g[temp_dat$setting_id == sett_ids[iter]]),
      post_delta_m = post_delta_m$Median,
      post_delta_m_low = post_delta_m$CI_low,
      post_delta_m_high = post_delta_m$CI_high,
      evidence = evidence,
      es = es,
      best_ids = best_ids,
      iter = iter,
      N_datasets = N_datasets
    )

  } else {

    # save output
    results[[i]] <- data.frame(
      case = model_control$case[i],
      cors = model_control$cors[i],
      N = model_control$N[i],
      case_id = model_control$case_ids[i],
      m_best = NA,
      m_diff = NA,
      post_delta_m = NA,
      post_delta_m_low = NA,
      post_delta_m_high = NA,
      evidence = NA,
      es = NA,
      best_ids = NA,
      iter = NA,
      N_datasets = N_datasets
    )

  }
}

res_deviance <- do.call(rbind, results)
saveRDS(res_deviance, paste0("output/ad/results_deviance_NE_settings_",
        sampsize, "_", current_case_ids, ".RDS"))


ov_t2 <- Sys.time()
ov_t2 - t1
