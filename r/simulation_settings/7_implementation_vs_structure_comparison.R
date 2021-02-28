if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if(!require(EFAtools)) install.packages("EFAtools"); library(EFAtools)

### helper function ============================================================
describe <- function(x, d = 2) {
  mn <- mean(x, na.rm = TRUE)
  mdn <- median(x, na.rm = TRUE)
  low <- min(x, na.rm = TRUE)
  up <- max(x, na.rm = TRUE)
  
  st <- paste0("M = ", EFAtools:::.numformat(mn, d), "; Mdn = ",
               EFAtools:::.numformat(mdn, d), "; [",
               EFAtools:::.numformat(low, d), ", ",
               EFAtools:::.numformat(up, d), "]")
  
  cat(st)
}

### data preparation ===========================================================

model_control <- expand.grid(
  case = names(population_models$loadings)[1:27],
  cors = names(population_models$phis_3),
  N = c(450),
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


# read in data
recovery <- readRDS("output/simulation_settings/recovery_450.RDS")

### fit values between settings vs across cases ================================
dat <- recovery %>% 
  group_by(case_ids, setting_id) %>% 
  summarise(
    p_fc = mean(diff_factor_corres > 0, na.rm = TRUE),
    m_g = mean(g, na.rm = TRUE),
    p_hey = mean(heywood, na.rm = TRUE)
  ) %>% 
  ungroup()

### RMSE
# across implementations
dat %>% 
  group_by(case_ids) %>% 
  arrange(m_g) %>% slice(1, n()) %>%
  mutate(diff_g = m_g - lag(m_g)) %>% 
  drop_na() %>% 
  pull(diff_g) %>% 
  describe()

# across structures
dat %>% 
  group_by(setting_id) %>% 
  arrange(m_g) %>% slice(1, n()) %>%
  mutate(diff_g = m_g - lag(m_g)) %>% 
  drop_na() %>% 
  pull(diff_g) %>% 
  describe()

### Heywood cases
# across implementations
dat %>% 
  group_by(case_ids) %>% 
  arrange(p_hey) %>% slice(1, n()) %>%
  mutate(diff_hey = p_hey - lag(p_hey)) %>% 
  drop_na() %>% 
  pull(diff_hey) %>% 
  describe()

# across structures
dat %>% 
  group_by(setting_id) %>%
  arrange(p_hey) %>% slice(1, n()) %>%
  mutate(diff_hey = p_hey - lag(p_hey)) %>% 
  drop_na() %>% 
  pull(diff_hey) %>% 
  describe()

### factor correspondences
# across implementations
dat %>% 
  group_by(case_ids) %>% 
  arrange(p_fc) %>% slice(1, n()) %>%
  mutate(diff_p = p_fc - lag(p_fc)) %>% 
  drop_na() %>% 
  pull(diff_p) %>% 
  describe()

# across structures
dat %>% 
  group_by(setting_id) %>%
  arrange(p_fc) %>% slice(1, n()) %>%
  mutate(diff_p = p_fc - lag(p_fc)) %>% 
  drop_na() %>% 
  pull(diff_p) %>% 
  describe()

