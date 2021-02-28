fl <- list.files("simulation_settings/output/mr/", full.names = TRUE)

dat <- do.call(rbind, lapply(fl, readRDS))
library(EFAtools)
model_control <- expand.grid(
  case = names(population_models$loadings)[1:27],
  cors = names(population_models$phis_3),
  N = c(450), # repeat with 450
  stringsAsFactors = FALSE)
model_control$case_ids <- 1:nrow(model_control)

current_case_ids <- 1:nrow(model_control)

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

dat$dat_id <- rep(rep(1:1000, each = nrow(settings)), nrow(model_control))
library(tibble)
library(tidyr)
library(dplyr)
recovery <- dat %>%
  mutate_if(is.factor, as.character) %>%
  as_tibble() %>%
  left_join(settings, by = c("comm_meth", "criterion_type", "abs_eigen",
                             "conv_crit", "var_type", "p_type", "k")) %>% 
  left_join(model_control, by = c("case", "cors", "N")) %>% 
  filter(case_ids %in% current_case_ids)

saveRDS(recovery, "simulation_settings/output/recovery_450.RDS")
