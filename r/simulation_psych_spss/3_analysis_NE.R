# Model recovery analysis

if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if(!require(EFAtools)) install.packages("EFAtools"); library(EFAtools)
if(!require(bayestestR)) install.packages("bayestestR"); library(bayestestR)
if(!require(rstanarm)) install.packages("rstanarm"); library(rstanarm)
options(mc.cores = parallel::detectCores())
if(!require(viridis)) install.packages("viridis"); library(viridis)


# read in data
recovery <- readRDS("output/simulation_psych_spss/model_recovery_results.RDS")

sim_control <- expand.grid(
  case = names(population_models$loadings)[1:27],
  cors = names(population_models$phis_3),
  N = c(180, 450),
  stringsAsFactors = FALSE)

recovery$dat_id <- rep(1:1000, nrow(sim_control))

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
  
  dat_ids_excl <- recovery %>%
    filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
             N == sim_control$N[i] & init_psych == "unity") %>%
    pull(dat_id)
  
  N_datasets <- 1000 - length(dat_ids_excl)
  
  if (N_datasets > 5) {
    # only keep relevant data
    temp_dat <- recovery %>%
      filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
               N == sim_control$N[i] & !(dat_id %in% dat_ids_excl)) %>%
      select(-starts_with("error")) %>%
      drop_na()
    
    # run test and calculate effect size
    mod <- stan_glm(g ~ implementation,
                    data = temp_dat %>%
                      select(g_psych, g_spss) %>%
                      pivot_longer(everything(), values_to = "g",
                                   names_to = "implementation", names_prefix = "g_"),
                    chains = 4, iter = 15000, family = gaussian(),
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
      m_psych = mean(temp_dat$g_psych),
      m_spss = mean(temp_dat$g_spss),
      post_delta_m = post_delta_m$Median,
      post_delta_m_low = post_delta_m$CI_low,
      post_delta_m_high = post_delta_m$CI_high,
      evidence = evidence,
      es = es,
      N_datasets = N_datasets
    )
  } else {
    
    # save output
    res_list[[i]] <- data.frame(
      case = sim_control$case[i],
      cors = sim_control$cors[i],
      N = sim_control$N[i],
      m_psych = NA,
      m_spss = NA,
      post_delta_m = NA,
      post_delta_m_low = NA,
      post_delta_m_high = NA,
      evidence = NA,
      es = NA,
      N_datasets = N_datasets
    )
    
  }
  
}

res_deviance_no_NE <- do.call(rbind, res_list)
saveRDS(res_deviance_no_NE, "output/simulation_psych_spss/results_deviance_psych_spss_no_NE.RDS")

### heywood cases ==============================================================

res_list <- list()
for (i in 1:nrow(sim_control)) {
  
  dat_ids_excl <- recovery %>%
    filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
             N == sim_control$N[i] & init_psych == "unity") %>%
    pull(dat_id)
  
  N_datasets <- 1000 - length(dat_ids_excl)
  
  if (N_datasets > 5) {
    
    # only keep relevant data
    temp_dat <- recovery %>%
      filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
               N == sim_control$N[i] & !(dat_id %in% dat_ids_excl)) %>%
      select(-starts_with("error")) %>%
      select(heywood_psych, heywood_spss) %>%
      pivot_longer(everything(), names_to = "implementation", values_to = "heywood",
                   names_prefix = "heywood_")
    
    prop_heywood_psych <- mean(temp_dat$heywood[temp_dat$implementation == "psych"], na.rm = TRUE)
    prop_heywood_spss <- mean(temp_dat$heywood[temp_dat$implementation == "spss"], na.rm = TRUE)
    
    if (mean(temp_dat$heywood) > .05) {
      # run logistic regression
      mod <- stan_glm(heywood ~ implementation, data = temp_dat,
                      chains = 4, iter = 15000, family = binomial(),
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
      evidence = evidence,
      N_datasets = N_datasets
    )
    
  } else {
    # save output
    res_list[[i]] <- data.frame(
      case = sim_control$case[i],
      cors = sim_control$cors[i],
      N = sim_control$N[i],
      OR = NA,
      OR_low = NA,
      OR_high = NA,
      evidence = NA,
      N_datasets = N_datasets
    )
  }
  
}

res_heywood_no_NE <- do.call(rbind, res_list)

saveRDS(res_heywood_no_NE, "output/simulation_psych_spss/results_heywood_psych_spss_no_NE.RDS")

### variable to factor correspondences =========================================

res_list <- list()
for (i in 1:nrow(sim_control)) {
  
  dat_ids_excl <- recovery %>%
    filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
             N == sim_control$N[i] & init_psych == "unity") %>%
    pull(dat_id)
  
  N_datasets <- 1000 - length(dat_ids_excl)
  
  if (N_datasets > 5) {
    
    temp_dat <- recovery %>%
      filter(case == sim_control$case[i] &cors == sim_control$cors[i] &
               N == sim_control$N[i] & !(dat_id %in% dat_ids_excl)) %>%
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
                      chains = 4, iter = 15000, family = neg_binomial_2,
                      prior_intercept = normal(location = 0, scale = 10),
                      prior = normal(location = 0, scale = 2.5), mean_PPD = FALSE)
      
      post_mod <- hdi(mod, ci = .95)
      
      evidence <- if (sign(post_mod$CI_low[2]) == sign(post_mod$CI_high[2])) {
        "Different"
      } else {
        "Inconclusive"
      }
      
      
      post_delta_m <- median(as.data.frame(mod)[,2])
      post_delta_m_low <- post_mod$CI_low[2]
      post_delta_m_high <- post_mod$CI_high[2]
      
    } else {
      
      post_delta_m <- NA
      post_delta_m_low <- NA
      post_delta_m_high <- NA
      evidence <- "Inconclusive"
      
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
      N_datasets = N_datasets
    )
    
  } else {
    
    # save output
    res_list[[i]] <- data.frame(
      case = sim_control$case[i],
      cors = sim_control$cors[i],
      N = sim_control$N[i],
      m_psych = NA,
      m_spss = NA,
      post_delta_m = NA,
      post_delta_m_low = NA,
      post_delta_m_high = NA,
      evidence = NA,
      N_datasets = N_datasets
    )
    
  }
  
}

res_factor_corres_no_NE <- do.call(rbind, res_list)

saveRDS(res_factor_corres_no_NE, "output/simulation_psych_spss/results_factor_corres_psych_spss_no_NE.RDS")


### Only cases with negative initial Eigenvalues ===============================

### deviance ===================================================================

res_list <- list()
for (i in 1:nrow(sim_control)) {
  
  dat_ids_excl <- recovery %>%
    filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
             N == sim_control$N[i] & init_psych == "smc") %>%
    pull(dat_id)
  
  N_datasets <- 1000 - length(dat_ids_excl)
  
  if (N_datasets > 5) {
    # only keep relevant data
    temp_dat <- recovery %>%
      filter(case == sim_control$case[i] &cors == sim_control$cors[i] &
               N == sim_control$N[i] & !(dat_id %in% dat_ids_excl)) %>%
      select(-starts_with("error")) %>%
      drop_na()
    
    # run test and calculate effect size
    mod <- stan_glm(g ~ implementation,
                    data = temp_dat %>%
                      select(g_psych, g_spss) %>%
                      pivot_longer(everything(), values_to = "g",
                                   names_to = "implementation", names_prefix = "g_"),
                    chains = 4, iter = 15000, family = gaussian(),
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
      es = es,
      N_datasets = N_datasets
    )
  } else {
    
    # save output
    res_list[[i]] <- data.frame(
      case = sim_control$case[i],
      cors = sim_control$cors[i],
      N = sim_control$N[i],
      m_psych = NA,
      m_spss = NA,
      post_delta_m = NA,
      post_delta_m_low = NA,
      post_delta_m_high = NA,
      evidence = NA,
      es = NA,
      N_datasets = N_datasets
    )
    
  }
  
}

res_deviance_NE <- do.call(rbind, res_list)
saveRDS(res_deviance_NE, "output/simulation_psych_spss/results_deviance_psych_spss_NE.RDS")

### heywood cases ==============================================================

res_list <- list()
for (i in 1:nrow(sim_control)) {
  
  dat_ids_excl <- recovery %>%
    filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
             N == sim_control$N[i] & init_psych == "smc") %>%
    pull(dat_id)
  
  N_datasets <- 1000 - length(dat_ids_excl)
  
  if (N_datasets > 5) {
    
    # only keep relevant data
    temp_dat <- recovery %>%
      filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
               N == sim_control$N[i] & !(dat_id %in% dat_ids_excl)) %>%
      select(-starts_with("error")) %>%
      select(heywood_psych, heywood_spss) %>%
      pivot_longer(everything(), names_to = "implementation", values_to = "heywood",
                   names_prefix = "heywood_")
    
    prop_heywood_psych <- mean(temp_dat$heywood[temp_dat$implementation == "psych"], na.rm = TRUE)
    prop_heywood_spss <- mean(temp_dat$heywood[temp_dat$implementation == "spss"], na.rm = TRUE)
    
    if (mean(temp_dat$heywood) > .05) {
      # run logistic regression
      mod <- stan_glm(heywood ~ implementation, data = temp_dat,
                      chains = 4, iter = 15000, family = binomial(),
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
      evidence = evidence,
      N_datasets = N_datasets
    )
    
  } else {
    # save output
    res_list[[i]] <- data.frame(
      case = sim_control$case[i],
      cors = sim_control$cors[i],
      N = sim_control$N[i],
      OR = NA,
      OR_low = NA,
      OR_high = NA,
      evidence = NA,
      N_datasets = N_datasets
    )
  }
  
}

res_heywood_NE <- do.call(rbind, res_list)

saveRDS(res_heywood_NE, "output/simulation_psych_spss/results_heywood_psych_spss_NE.RDS")

### variable to factor correspondences =========================================

res_list <- list()
for (i in 1:nrow(sim_control)) {
  
  dat_ids_excl <- recovery %>%
    filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
             N == sim_control$N[i] & init_psych == "smc") %>%
    pull(dat_id)
  
  N_datasets <- 1000 - length(dat_ids_excl)
  
  if (N_datasets > 5) {
    
    temp_dat <- recovery %>%
      filter(case == sim_control$case[i] & cors == sim_control$cors[i] &
               N == sim_control$N[i] & !(dat_id %in% dat_ids_excl)) %>%
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
                      chains = 4, iter = 15000, family = neg_binomial_2,
                      prior_intercept = normal(location = 0, scale = 10),
                      prior = normal(location = 0, scale = 2.5), mean_PPD = FALSE)
      
      post_mod <- hdi(mod, ci = .95)
      
      evidence <- if (sign(post_mod$CI_low[2]) == sign(post_mod$CI_high[2])) {
        "Different"
      } else {
        "Inconclusive"
      }
      
      
      post_delta_m <- median(as.data.frame(mod)[,2])
      post_delta_m_low <- post_mod$CI_low[2]
      post_delta_m_high <- post_mod$CI_high[2]
      
    } else {
      
      post_delta_m <- NA
      post_delta_m_low <- NA
      post_delta_m_high <- NA
      evidence <- "Inconclusive"
      
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
      N_datasets = N_datasets
    )
    
  } else {
    
    # save output
    res_list[[i]] <- data.frame(
      case = sim_control$case[i],
      cors = sim_control$cors[i],
      N = sim_control$N[i],
      m_psych = NA,
      m_spss = NA,
      post_delta_m = NA,
      post_delta_m_low = NA,
      post_delta_m_high = NA,
      evidence = NA,
      N_datasets = N_datasets
    )
    
  }
  
  
  
}

res_factor_corres_NE <- do.call(rbind, res_list)

saveRDS(res_factor_corres_NE, "output/simulation_psych_spss/results_factor_corres_psych_spss_NE.RDS")


### Analyses ===================================================================

res_deviance_ne <- readRDS("output/simulation_psych_spss/results_deviance_psych_spss_NE.RDS")
res_heywood_ne <- readRDS("output/simulation_psych_spss/results_heywood_psych_spss_NE.RDS")
res_factor_corres_ne <- readRDS("output/simulation_psych_spss/results_factor_corres_psych_spss_NE.RDS")

res_deviance_no_ne <- readRDS("output/simulation_psych_spss/results_deviance_psych_spss_no_NE.RDS")
res_heywood_no_ne <- readRDS("output/simulation_psych_spss/results_heywood_psych_spss_no_NE.RDS")
res_factor_corres_no_ne <- readRDS("output/simulation_psych_spss/results_factor_corres_psych_spss_no_NE.RDS")

# admissible solutions
six_factors <- names(population_models$loadings)[18:27]
recovery %>% 
  mutate(nfactors = case_when(case %in% six_factors ~ 6,
                              TRUE ~ 3)) %>% 
  group_by(nfactors) %>% 
  summarise(
    m_psych = mean(nfac_admiss_psych, na.rm = TRUE),
    m_spss = mean(nfac_admiss_spss, na.rm = TRUE)
  )

# distribution of the evidence

# RMSE no NE
table(res_deviance_no_ne$evidence)

res_deviance_no_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(post_delta_m > 0 ~ "psych",
                            post_delta_m < 0 ~ "SPSS")) %>%
  select(better) %>% 
  table()

res_deviance_no_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(post_delta_m > 0 ~ "psych",
                            post_delta_m < 0 ~ "SPSS")) %>%
  select(better, case, cors) %>%
  table()

# Heywood no NE
table(res_heywood_no_ne$evidence)

res_heywood_no_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(OR > 1 ~ "psych",
                            OR < 1 ~ "SPSS")) %>%
  select(better) %>% 
  table()

res_heywood_no_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(OR > 1 ~ "psych",
                            OR < 1 ~ "SPSS")) %>%
  select(better, case, cors) %>% 
  table()

# Factor correspondence no NE
table(res_factor_corres_no_ne$evidence)

res_factor_corres_no_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(post_delta_m > 0 ~ "psych",
                            post_delta_m < 0 ~ "SPSS")) %>%
  select(better) %>% 
  table()

res_factor_corres_no_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(post_delta_m > 0 ~ "psych",
                            post_delta_m < 0 ~ "SPSS")) %>%
  select(better, case, cors) %>% 
  table()

# RMSE NE
table(res_deviance_ne$evidence)

res_deviance_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(post_delta_m > 0 ~ "psych",
                            post_delta_m < 0 ~ "SPSS")) %>%
  select(better) %>% 
  table()

res_deviance_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(post_delta_m > 0 ~ "psych",
                            post_delta_m < 0 ~ "SPSS")) %>%
  select(better, case, cors) %>%
  table()

# Heywood NE
table(res_heywood_ne$evidence)

res_heywood_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(OR > 1 ~ "psych",
                            OR < 1 ~ "SPSS")) %>%
  select(better) %>% 
  table()

res_heywood_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(OR > 1 ~ "psych",
                            OR < 1 ~ "SPSS")) %>%
  select(better, case, cors) %>% 
  table()

# Factor correspondence NE
table(res_factor_corres_ne$evidence)

res_factor_corres_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(post_delta_m > 0 ~ "psych",
                            post_delta_m < 0 ~ "SPSS")) %>%
  select(better) %>% 
  table()

res_factor_corres_ne %>%
  filter(evidence == "Different") %>%
  mutate(better = case_when(post_delta_m > 0 ~ "psych",
                            post_delta_m < 0 ~ "SPSS")) %>%
  select(better, case, cors) %>% 
  table()


# summary statistics

temp <- recovery %>% 
  #  filter(init_psych == "smc") %>% 
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong"))) %>%
  group_by(case, cors) %>% 
  summarise(
    prop_psych = mean(heywood_psych, na.rm = TRUE),
    prop_spss = mean(heywood_spss, na.rm = TRUE),
    fit_psych = mean(g_psych, na.rm = TRUE),
    fit_spss = mean(g_spss, na.rm = TRUE),
    prop_fc_psych = mean(diff_factor_corres_psych > 0, na.rm = TRUE),
    prop_fc_spss = mean(diff_factor_corres_spss > 0, na.rm = TRUE),
    diff_fc_psych_spss = mean(diff_factor_corres_spss_psych, na.rm = TRUE),
    m_load_psych_spss = mean(m_delta_load_spss_psych, na.rm = TRUE),
    max_load_psych_spss = mean(max_delta_load_spss_psych, na.rm = TRUE),
    admiss_psych = mean(nfac_admiss_psych, na.rm = TRUE),
    admiss_spss = mean(nfac_admiss_spss, na.rm = TRUE)
  )

summary(temp)


### Printout tables ============================================================


nvec <- stringr::str_to_title(gsub("case_", " ", names(population_models$loadings)[1:27])) 
cases <- names(population_models$loadings)

get_col <- function(x, ind) {
  if (is.na(x$evidence[ind]) || x$evidence[ind] == "Inconclusive") {
    whichcol <- "incon"
  } else if (x$evidence[ind] == "Equal") {
    whichcol <- "equal"
  } else {
    if(x$post_delta_m[ind] > 0) {
      whichcol <- "psych"
    } else {
      whichcol <- "spss"
    }
  }
  
  paste0("\\cellcolor{", whichcol, "} ")
}

get_col_h <- function(x, ind) {
  if (is.na(x$evidence[ind]) || x$evidence[ind] == "Inconclusive") {
    whichcol <- "incon"
  } else if (x$evidence[ind] == "Equal") {
    whichcol <- "equal"
  } else {
    if(x$OR[ind] > 1) {
      whichcol <- "psych"
    } else {
      whichcol <- "spss"
    }
  }
  
  paste0("\\cellcolor{", whichcol, "} ")
}


cell <- function(x, ind, es) {
  
  if (is.na(x$post_delta_m[ind])) {
    "-"
  } else {
    paste0(EFAtools:::.numformat(x$post_delta_m[ind], 2), " [",
           EFAtools:::.numformat(x$post_delta_m_low[ind], 2), ", ",
           EFAtools:::.numformat(x$post_delta_m_high[ind], 2), "], $", es, "$ = ",
           round(x$es[ind], 2))
  }
  
}


cell_h <- function(x, ind) {
  
  if (is.na(x$OR[ind])) {
    "-"
  } else {
    
    ll <- ifelse(is.na(x$OR_low[ind]), "-",
                 EFAtools:::.numformat(x$OR_low[ind], 2))
    hh <- ifelse(is.na(x$OR_high[ind]), "-",
                 EFAtools:::.numformat(x$OR_high[ind], 2))
    paste0(EFAtools:::.numformat(x$OR[ind], 2), " [",
           ll , ", ",
           hh , "]")
  }
  
}

cell_f <- function(x, ind) {
  
  if (is.na(x$post_delta_m[ind])) {
    "-"
  } else {
    ll <- ifelse(is.na(x$post_delta_m_low[ind]), "-",
                 EFAtools:::.numformat(x$post_delta_m_low[ind], 2))
    hh <- ifelse(is.na(x$post_delta_m_high[ind]), "-",
                 EFAtools:::.numformat(x$post_delta_m_high[ind], 2))
    
    paste0(EFAtools:::.numformat(x$post_delta_m[ind], 2), " [",
           ll, ", ",
           hh, "]")
  }
  
}

# table deviance

N <- 450

rows <- c()
for (i in 1:length(nvec)) {
  ind_l <- res_deviance_no_ne$case == cases[i] & res_deviance_no_ne$cors == "zero" &
    res_deviance_no_ne$N == N
  ind_m <- res_deviance_no_ne$case == cases[i] & res_deviance_no_ne$cors == "moderate" &
    res_deviance_no_ne$N == N
  ind_mi <-res_deviance_no_ne$case == cases[i] & res_deviance_no_ne$cors == "mixed" &
    res_deviance_no_ne$N == N
  ind_h <- res_deviance_no_ne$case == cases[i] & res_deviance_no_ne$cors == "strong" &
    res_deviance_no_ne$N == N
  rows[i] <- paste0(nvec[i], " & ", get_col(res_deviance_no_ne, ind_l),
                    cell(res_deviance_no_ne, ind_l, "d"), " & ",
                    get_col(res_deviance_no_ne, ind_m),
                    cell(res_deviance_no_ne, ind_m, "d"),
                    " & ", get_col(res_deviance_no_ne, ind_mi),
                    cell(res_deviance_no_ne, ind_mi, "d"),
                    " & ", get_col(res_deviance_no_ne, ind_h),
                    cell(res_deviance_no_ne, ind_h, "d"), " \\\\")
  
  
}

cat(paste(rows, collapse = "\n"))


# table deviance NE

N <- 450

rows <- c()
for (i in 1:length(nvec)) {
  ind_l <- res_deviance_ne$case == cases[i] & res_deviance_ne$cors == "zero" &
    res_deviance_ne$N == N
  ind_m <- res_deviance_ne$case == cases[i] & res_deviance_ne$cors == "moderate" &
    res_deviance_ne$N == N
  ind_mi <-res_deviance_ne$case == cases[i] & res_deviance_ne$cors == "mixed" &
    res_deviance_ne$N == N
  ind_h <- res_deviance_ne$case == cases[i] & res_deviance_ne$cors == "strong" &
    res_deviance_ne$N == N
  rows[i] <- paste0(nvec[i], " & ", get_col(res_deviance_ne, ind_l),
                    cell(res_deviance_ne, ind_l, "d"), " & ",
                    get_col(res_deviance_ne, ind_m),
                    cell(res_deviance_ne, ind_m, "d"),
                    " & ", get_col(res_deviance_ne, ind_mi),
                    cell(res_deviance_ne, ind_mi, "d"),
                    " & ", get_col(res_deviance_ne, ind_h),
                    cell(res_deviance_ne, ind_h, "d"), " \\\\")
  
  
}

cat(paste(rows, collapse = "\n"))




# table factor correspondence

# no negative eigenvalues
N <- 450

rows <- c()
for (i in 1:length(nvec)) {
  ind_l <- res_factor_corres_no_ne$case == cases[i] & res_factor_corres_no_ne$cors == "zero" &
    res_factor_corres_no_ne$N == N
  ind_m <- res_factor_corres_no_ne$case == cases[i] & res_factor_corres_no_ne$cors == "moderate" &
    res_factor_corres_no_ne$N == N
  ind_mi <-res_factor_corres_no_ne$case == cases[i] & res_factor_corres_no_ne$cors == "mixed" &
    res_factor_corres_no_ne$N == N
  ind_h <- res_factor_corres_no_ne$case == cases[i] & res_factor_corres_no_ne$cors == "strong" &
    res_factor_corres_no_ne$N == N
  rows[i] <- paste0(nvec[i], " & ", get_col(res_factor_corres_no_ne, ind_l),
                    cell_f(res_factor_corres_no_ne, ind_l), " & ",
                    get_col(res_factor_corres_no_ne, ind_m),
                    cell_f(res_factor_corres_no_ne, ind_m),
                    " & ", get_col(res_factor_corres_no_ne, ind_mi),
                    cell_f(res_factor_corres_no_ne, ind_mi),
                    " & ", get_col(res_factor_corres_no_ne, ind_h),
                    cell_f(res_factor_corres_no_ne, ind_h), " \\\\")
  
  
}

cat(paste(rows, collapse = "\n"))


# negative eigenvalues
N <- 450

rows <- c()
for (i in 1:length(nvec)) {
  ind_l <- res_factor_corres_ne$case == cases[i] & res_factor_corres_ne$cors == "zero" &
    res_factor_corres_ne$N == N
  ind_m <- res_factor_corres_ne$case == cases[i] & res_factor_corres_ne$cors == "moderate" &
    res_factor_corres_ne$N == N
  ind_mi <-res_factor_corres_ne$case == cases[i] & res_factor_corres_ne$cors == "mixed" &
    res_factor_corres_ne$N == N
  ind_h <- res_factor_corres_ne$case == cases[i] & res_factor_corres_ne$cors == "strong" &
    res_factor_corres_ne$N == N
  rows[i] <- paste0(nvec[i], " & ", get_col(res_factor_corres_ne, ind_l),
                    cell_f(res_factor_corres_ne, ind_l), " & ",
                    get_col(res_factor_corres_ne, ind_m),
                    cell_f(res_factor_corres_ne, ind_m),
                    " & ", get_col(res_factor_corres_ne, ind_mi),
                    cell_f(res_factor_corres_ne, ind_mi),
                    " & ", get_col(res_factor_corres_ne, ind_h),
                    cell_f(res_factor_corres_ne, ind_h), " \\\\")
  
  
}

cat(paste(rows, collapse = "\n"))


# table heywood

# no negative eigenvalues
N <- 450

rows <- c()
for (i in 1:length(nvec)) {
  ind_l <- res_heywood_no_ne$case == cases[i] & res_heywood_no_ne$cors == "zero" &
    res_heywood_no_ne$N == N
  ind_m <- res_heywood_no_ne$case == cases[i] & res_heywood_no_ne$cors == "moderate" &
    res_heywood_no_ne$N == N
  ind_mi <-res_heywood_no_ne$case == cases[i] & res_heywood_no_ne$cors == "mixed" &
    res_heywood_no_ne$N == N
  ind_h <- res_heywood_no_ne$case == cases[i] & res_heywood_no_ne$cors == "strong" &
    res_heywood_no_ne$N == N
  rows[i] <- paste0(nvec[i], " & ", get_col_h(res_heywood_no_ne, ind_l),
                    cell_h(res_heywood_no_ne, ind_l), " & ",
                    get_col_h(res_heywood_no_ne, ind_m),
                    cell_h(res_heywood_no_ne, ind_m),
                    " & ", get_col_h(res_heywood_no_ne, ind_mi),
                    cell_h(res_heywood_no_ne, ind_mi),
                    " & ", get_col_h(res_heywood_no_ne, ind_h),
                    cell_h(res_heywood_no_ne, ind_h), " \\\\")
  
  
}

cat(paste(rows, collapse = "\n"))


# negative eigenvalues
N <- 450

rows <- c()
for (i in 1:length(nvec)) {
  ind_l <- res_heywood_ne$case == cases[i] & res_heywood_ne$cors == "zero" &
    res_heywood_ne$N == N
  ind_m <- res_heywood_ne$case == cases[i] & res_heywood_ne$cors == "moderate" &
    res_heywood_ne$N == N
  ind_mi <-res_heywood_ne$case == cases[i] & res_heywood_ne$cors == "mixed" &
    res_heywood_ne$N == N
  ind_h <- res_heywood_ne$case == cases[i] & res_heywood_ne$cors == "strong" &
    res_heywood_ne$N == N
  rows[i] <- paste0(nvec[i], " & ", get_col_h(res_heywood_ne, ind_l),
                    cell_h(res_heywood_ne, ind_l), " & ",
                    get_col_h(res_heywood_ne, ind_m),
                    cell_h(res_heywood_ne, ind_m),
                    " & ", get_col_h(res_heywood_ne, ind_mi),
                    cell_h(res_heywood_ne, ind_mi),
                    " & ", get_col_h(res_heywood_ne, ind_h),
                    cell_h(res_heywood_ne, ind_h), " \\\\")
  
  
}

cat(paste(rows, collapse = "\n"))




