args <- commandArgs(TRUE)
eval(parse(text = args))
t1 <- Sys.time()
# Model recovery analysis

library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(EFAtools)
library(bayestestR)
library(rstanarm)

# read in data
recovery <-  readRDS(paste0(path, "/model_recovery_results.RDS"))

sim_control <- expand.grid(
  case = names(population_models$loadings)[1:27],
  cors = names(population_models$phis_3),
  N = c(180, 450),
  stringsAsFactors = FALSE)

# cohen's d function
cohens_d <- function(x, y) {
  
  if (any(is.na(x))) {
    xc <- x[stats::complete.cases(x), ]
    y <- y[stats::complete.cases(x), ]
    x <- xc
  }
  if (any(is.na(y))) {
    yc <- y[stats::complete.cases(y), ]
    x <- x[stats::complete.cases(y), ]
    y <- yc
  }
  
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


### ============================================================================
### CASE BY CASE ===============================================================
### ============================================================================

### deviance ===================================================================

res_list <- list()
for (i in 1:nrow(sim_control)) {

  # only keep relevant data
  temp_dat <- recovery %>%
    filter(case == sim_control$case[i] &cors == sim_control$cors[i] &
           N == sim_control$N[i]) %>%
    select(-starts_with("error")) %>%
    drop_na()

  # run test and calculate effect size
  mod <- stan_glm(g ~ implementation,
                  data = temp_dat %>%
                    select(g_psych, g_spss) %>%
                    pivot_longer(everything(), values_to = "g",
                                 names_to = "implementation", names_prefix = "g_"),
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
  es <- cohens_d(temp_dat$g_spss, temp_dat$g_psych)

  # save output
  res_list[[i]] <- data.frame(
    case = sim_control$case[i],
    cors = sim_control$cors[i],
    N = sim_control$N[i],
    m_psych = mean(temp_dat$g_psych, na.rm = TRUE),
    m_spss = mean(temp_dat$g_spss, na.rm = TRUE),
    post_delta_m = post_delta_m$Median,
    post_delta_m_low = post_delta_m$CI_low,
    post_delta_m_high = post_delta_m$CI_high,
    evidence = evidence,
    es = es
  )
}

res_deviance <- do.call(rbind, res_list)
saveRDS(res_deviance, "output/analyses/results_deviance_psych_spss.RDS")

### heywood cases ==============================================================

res_list <- list()
for (i in 1:nrow(sim_control)) {

  # only keep relevant data
  temp_dat <- recovery %>%
    filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
           N == sim_control$N[i]) %>%
    select(-starts_with("error")) %>%
    drop_na() %>% 
    select(heywood_psych, heywood_spss) %>%
    pivot_longer(everything(), names_to = "implementation", values_to = "heywood",
                 names_prefix = "heywood_")

  prop_heywood_psych <- mean(temp_dat$heywood[temp_dat$implementation == "psych"], na.rm = TRUE)
  prop_heywood_spss <- mean(temp_dat$heywood[temp_dat$implementation == "spss"], na.rm = TRUE)

  if (mean(temp_dat$heywood) > .05) {
    # run logistic regression
    mod <- stan_glm(heywood ~ implementation, data = temp_dat,
                    chains = 4, iter = 10000, family = binomial(),
                    prior_intercept = normal(location = 0, scale = 10),
                    prior = normal(location = 0, scale = 2.5))

    post_mod <- describe_posterior(mod, ci = .95, rope_ci = 1)

    evidence <- if (sign(post_mod$CI_low[2]) == sign(post_mod$CI_high[2]) &
                    post_mod$ROPE_Percentage[2] < .025) {
      "Different"
    } else if (post_mod$ROPE_Percentage[2] > .95) {
      "Equal"
    } else {
      "Inconclusive"
    }

    OR <- exp(post_mod$Median[2])
    OR_low <- exp(post_mod$CI_low[2])
    OR_high <- exp(post_mod$CI_high[2])

  } else {

    if (prop_heywood_psych == 0 & prop_heywood_spss == 0) {

      evidence <- "Equal"
      OR <- 1
      OR_low <- 1
      OR_high <- 1

    } else {

      OR <- (prop_heywood_spss / (1-prop_heywood_psych)) /
        (prop_heywood_psych / (1- prop_heywood_psych))
      OR_low <- NA
      OR_high <- NA

      evidence <- if (OR == 1) {
        "Equal"
      } else {
        "Inconclusive"
      }

    }

  }

  # save output
  res_list[[i]] <- data.frame(
    case = sim_control$case[i],
    cors = sim_control$cors[i],
    N = sim_control$N[i],
    OR = OR,
    OR_low = OR_low,
    OR_high = OR_high,
    evidence = evidence
  )
}

res_heywood <- do.call(rbind, res_list)

saveRDS(res_heywood, "output/analyses/results_heywood_psych_spss.RDS")

### variable to factor correspondences =========================================

res_list <- list()
for (i in 1:nrow(sim_control)) {

  temp_dat <- recovery %>%
    filter(case == sim_control$case[i] &cors == sim_control$cors[i] &
             N == sim_control$N[i]) %>%
    select(-starts_with("error")) %>%
    drop_na()

  if (length(unique(temp_dat$diff_factor_corres_psych)) >= 2 &
      length(unique(temp_dat$diff_factor_corres_spss)) >= 2) {

    # run test and calculate effect size
    mod <- stan_glm(factor_corres ~ implementation,
                    data = temp_dat %>%
                      select(diff_factor_corres_psych, diff_factor_corres_spss) %>%
                      pivot_longer(everything(), values_to = "factor_corres",
                                   names_to = "implementation",
                                   names_prefix = "diff_factor_corres_"),
                    chains = 4, iter = 10000, family = gaussian(),
                    prior_intercept = normal(location = 0, scale = 10),
                    prior = normal(location = 0, scale = 2.5))

    post_mod <- describe_posterior(mod, ci = .95, rope_ci = 1)[2,]

    evidence <- if (sign(post_mod$CI_low) == sign(post_mod$CI_high) &
                    post_mod$ROPE_Percentage < .025) {
      "Different"
    } else if (post_mod$ROPE_Percentage > .95) {
      "Equal"
    } else {
      "Inconclusive"
    }


    post_delta_m <- post_mod$Median
    post_delta_m_low <- post_mod$CI_low
    post_delta_m_high <- post_mod$CI_high

    es <- cohens_d(temp_dat$diff_factor_corres_spss,
                   temp_dat$diff_factor_corres_psych)

  } else {

    post_delta_m <- mean(temp_dat$diff_factor_corres_psych, na.rm = TRUE) -
      mean(temp_dat$diff_factor_corres_spss, na.rm = TRUE)
    post_delta_m_low <- NA
    post_delta_m_high <- NA
    evidence <- if (mean(temp_dat$diff_factor_corres_psych, na.rm = TRUE) ==
                    mean(temp_dat$diff_factor_corres_spss, na.rm = TRUE)) {
      "Equal"
    } else {
      "Inconclusive"
    }

    es <- if (mean(temp_dat$diff_factor_corres_psych, na.rm = TRUE) ==
                    mean(temp_dat$diff_factor_corres_spss, na.rm = TRUE)) {
      0
    } else {
      cohens_d(temp_dat$diff_factor_corres_spss,
               temp_dat$diff_factor_corres_psych)
    }
  }

  # save output
  res_list[[i]] <- data.frame(
    case = sim_control$case[i],
    cors = sim_control$cors[i],
    N = sim_control$N[i],
    m_psych = mean(temp_dat$diff_factor_corres_psych, na.rm = TRUE),
    m_spss = mean(temp_dat$diff_factor_corres_spss, na.rm = TRUE),
    post_delta_m = post_delta_m,
    post_delta_m_low = post_delta_m_low,
    post_delta_m_high = post_delta_m_high,
    evidence = evidence,
    es = es
  )
}

res_factor_corres <- do.call(rbind, res_list)

saveRDS(res_factor_corres, "output/analyses/results_factor_corres_psych_spss.RDS")

ov_t2 <- Sys.time()
ov_t2 - t1

