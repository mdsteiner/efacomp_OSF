##### Replication of original results from psych and SPSS with EFAtools

if (!require(psych)) install.packages("psych"); library(psych)
if (!require(devtools)) install.packages("devtools"); library(devtools)
if (!require(EFAtools)) install.packages("EFAtools"); library(EFAtools)
if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)

# Print decimal instead of exponential notation of numbers
options(scipen = 999)

# prepare data specifications
dat_list <- list(test_models$baseline, test_models$case_1a, test_models$case_6b,
                 test_models$case_11b,
                 list(cormat = DOSPERT$cormat, n_factors = 10, N = DOSPERT$N),
                 list(cormat = IDS2_R, n_factors = 5, N = 1991),
                 list(cormat = WJIV_ages_3_5$cormat, n_factors = 7,
                      N = WJIV_ages_3_5$N),
                 list(cormat = WJIV_ages_20_39$cormat, n_factors = 7,
                      N = WJIV_ages_20_39$N))
names_dat <- c(names(test_models), "DOSPERT", "IDS_2",
               "WJIV_3_5", "WJIV_20_39")
names(dat_list) <- names_dat
init_h2_psych <- c("baseline" = TRUE, "case_1a" = FALSE, "case_6b" = TRUE,
                   "case_11b" = TRUE, "DOSPERT" = TRUE, "IDS_2" = TRUE,
                   "WJIV_3_5" = TRUE, "WJIV_20_39" = TRUE)
init_h2_EFA_psych <- c("baseline" = "smc", "case_1a" = "unity", "case_6b" = "smc",
                       "case_11b" = "smc", "DOSPERT" = "smc", "IDS_2" = "smc",
                       "WJIV_3_5" = "smc", "WJIV_20_39" = "smc")

####  ======================================================
#### Principal Axis Factoring ==============================
####  ======================================================

# Implementations --------

# psych::fa and EFA with type = "psych", type = "SPSS"
psych_paf <- list()
EFA_psych_paf <- list()
EFA_SPSS_paf <- list()

for(i in names_dat){
  # psych::fa without rotation
  psych_paf[[i]] <- psych::fa(dat_list[[i]][["cormat"]],
                              nfactors = dat_list[[i]][["n_factors"]],
                              fm = "pa", rotate = "none", max.iter = 300,
                              SMC = init_h2_psych[i])
  # EFA with type = "psych" without rotation
  EFA_psych_paf[[i]] <- EFA(dat_list[[i]][["cormat"]],
                            n_factors = dat_list[[i]][["n_factors"]],
                            type = "psych", max_iter = 300,
                            init_comm = init_h2_EFA_psych[i])
  # EFA with type = "SPSS" without rotation
  EFA_SPSS_paf[[i]] <- EFA(dat_list[[i]][["cormat"]],
                           n_factors = dat_list[[i]][["n_factors"]],
                           type = "SPSS", max_iter = 300)
}

### Comparison of original results with results from EFA function ----------

paf_psych_rep_load <- list()
paf_SPSS_rep_load <- list()

for(i in names_dat){
  # Compare loadings from psych::fa and EFA with type = "psych"
  paf_psych_rep_load[[i]] <- COMPARE(EFA_psych_paf[[i]]$unrot_loadings,
                                     psych_paf[[i]]$loadings)
  # Compare loadings from SPSS and EFA with type = "SPSS"
  paf_SPSS_rep_load[[i]] <- COMPARE(EFA_SPSS_paf[[i]]$unrot_loadings,
                                    SPSS[[i]]$paf_load)
}


####  =============================================
#### Promax Rotation ==============================
####  =============================================

#### Varimax ==============================

### Implementations ----------

# psych::fa and EFA with type = "psych" and type = "SPSS" with varimax rotation
psych_var <- list()
EFA_psych_var <- list()
EFA_SPSS_var <- list()
Var_SPSS <- list()

for(i in names_dat){
  # psych::fa with varimax rotation
  psych_var[[i]] <- psych::fa(dat_list[[i]][["cormat"]],
                              nfactors = dat_list[[i]][["n_factors"]],
                              fm = "pa", rotate = "varimax", max.iter = 300,
                              SMC = init_h2_psych[i])

  # EFA function with psych configuration and varimax rotation
  EFA_psych_var[[i]] <- EFA(dat_list[[i]][["cormat"]],
                            n_factors = dat_list[[i]][["n_factors"]],
                            type = "psych", rotation = "varimax",
                            max_iter = 300, init_comm = init_h2_EFA_psych[i])

  # EFA function with SPSS configuration and varimax rotation
  EFA_SPSS_var[[i]] <- EFA(dat_list[[i]][["cormat"]],
                           n_factors = dat_list[[i]][["n_factors"]],
                           type = "SPSS", rotation = "varimax", max_iter = 300)
}

### Comparison of original results with results from EFA function -----

var_psych_rep_load <- list()
var_SPSS_rep_load <- list()
var_SPSS_new_load <- list()

for(i in names_dat){
  # Compare loadings from psych::fa and EFA with type = "psych"
  var_psych_rep_load[[i]] <- COMPARE(EFA_psych_var[[i]]$rot_loadings,
                                     psych_var[[i]]$loadings)

  # Compare loadings from SPSS and EFA with type = "SPSS"
  var_SPSS_rep_load[[i]] <- COMPARE(EFA_SPSS_var[[i]]$rot_loadings,
                                    SPSS[[i]]$var_load)
}

#### Promax ==============================

### Implementations ----------

# psych::fa and EFA with type = "psych" and type = "SPSS" with promax rotation
psych_pro <- list()
EFA_psych_pro <- list()
EFA_SPSS_pro <- list()

for(i in names_dat){
  # psych::fa with promax rotation
  psych_pro[[i]] <- psych::fa(dat_list[[i]][["cormat"]],
                              nfactors = dat_list[[i]][["n_factors"]],
                              fm = "pa", rotate = "Promax", max.iter = 300,
                              SMC = init_h2_psych[i])
  # EFA with type = "psych" with promax rotation
  EFA_psych_pro[[i]] <- EFA(dat_list[[i]][["cormat"]],
                            n_factors = dat_list[[i]][["n_factors"]],
                            type = "psych", rotation = "promax",
                            max_iter = 300, init_comm = init_h2_EFA_psych[i])
  # EFA with type = "SPSS" with promax rotation
  EFA_SPSS_pro[[i]] <- EFA(dat_list[[i]][["cormat"]],
                           n_factors = dat_list[[i]][["n_factors"]],
                           type = "SPSS", rotation = "promax", max_iter = 300)
}

### Comparison of original results with results from EFA function ----------
pro_psych_rep_load <- list()
pro_SPSS_rep_load <- list()

for(i in names_dat){
  # Compare pattern matrix from psych::fa and EFA with type = "psych"
  pro_psych_rep_load[[i]] <- COMPARE(EFA_psych_pro[[i]]$rot_loadings,
                                     psych_pro[[i]]$loadings)
  
  # Compare pattern matrix from SPSS and EFA with type = "SPSS"
  pro_SPSS_rep_load[[i]] <- COMPARE(EFA_SPSS_pro[[i]]$rot_loadings,
                                    SPSS[[i]]$pro_load)
}


####  ===========================
#### Save results ===============
####  ===========================

# Replications of original results with functions from EFAtools -----------

Rep_orig <- data.frame("Statistic" = c("paf_load", "", "", "", "", "", "", "",
                                       "var_load", "", "", "", "", "", "", "",
                                       "pro_load", "", "", "", "", "", "", ""),
                       "Data" = rep(c("baseline", "case_1a", "case_6b",
                                      "case_11b", "DOSPERT", "IDS_2",
                                      "WJIV_3_5", "WJIV_20_39"), 3),
                       "psych_results" = rep("", 24),
                       "Md_diff" = c(paf_psych_rep_load$baseline$median_abs_diff,
                                     paf_psych_rep_load$case_1a$median_abs_diff,
                                     paf_psych_rep_load$case_6b$median_abs_diff,
                                     paf_psych_rep_load$case_11b$median_abs_diff,
                                     paf_psych_rep_load$DOSPERT$median_abs_diff,
                                     paf_psych_rep_load$IDS_2$median_abs_diff,
                                     paf_psych_rep_load$WJIV_3_5$median_abs_diff,
                                     paf_psych_rep_load$WJIV_20_39$median_abs_diff,
                                     var_psych_rep_load$baseline$median_abs_diff,
                                     var_psych_rep_load$case_1a$median_abs_diff,
                                     var_psych_rep_load$case_6b$median_abs_diff,
                                     var_psych_rep_load$case_11b$median_abs_diff,
                                     var_psych_rep_load$DOSPERT$median_abs_diff,
                                     var_psych_rep_load$IDS_2$median_abs_diff,
                                     var_psych_rep_load$WJIV_3_5$median_abs_diff,
                                     var_psych_rep_load$WJIV_20_39$median_abs_diff,
                                     pro_psych_rep_load$baseline$median_abs_diff,
                                     pro_psych_rep_load$case_1a$median_abs_diff,
                                     pro_psych_rep_load$case_6b$median_abs_diff,
                                     pro_psych_rep_load$case_11b$median_abs_diff,
                                     pro_psych_rep_load$DOSPERT$median_abs_diff,
                                     pro_psych_rep_load$IDS_2$median_abs_diff,
                                     pro_psych_rep_load$WJIV_3_5$median_abs_diff,
                                     pro_psych_rep_load$WJIV_20_39$median_abs_diff),
                       "Max_diff" = c(paf_psych_rep_load$baseline$max_abs_diff,
                                      paf_psych_rep_load$case_1a$max_abs_diff,
                                      paf_psych_rep_load$case_6b$max_abs_diff,
                                      paf_psych_rep_load$case_11b$max_abs_diff,
                                      paf_psych_rep_load$DOSPERT$max_abs_diff,
                                      paf_psych_rep_load$IDS_2$max_abs_diff,
                                      paf_psych_rep_load$WJIV_3_5$max_abs_diff,
                                      paf_psych_rep_load$WJIV_20_39$max_abs_diff,
                                      var_psych_rep_load$baseline$max_abs_diff,
                                      var_psych_rep_load$case_1a$max_abs_diff,
                                      var_psych_rep_load$case_6b$max_abs_diff,
                                      var_psych_rep_load$case_11b$max_abs_diff,
                                      var_psych_rep_load$DOSPERT$max_abs_diff,
                                      var_psych_rep_load$IDS_2$max_abs_diff,
                                      var_psych_rep_load$WJIV_3_5$max_abs_diff,
                                      var_psych_rep_load$WJIV_20_39$max_abs_diff,
                                      pro_psych_rep_load$baseline$max_abs_diff,
                                      pro_psych_rep_load$case_1a$max_abs_diff,
                                      pro_psych_rep_load$case_6b$max_abs_diff,
                                      pro_psych_rep_load$case_11b$max_abs_diff,
                                      pro_psych_rep_load$DOSPERT$max_abs_diff,
                                      pro_psych_rep_load$IDS_2$max_abs_diff,
                                      pro_psych_rep_load$WJIV_3_5$max_abs_diff,
                                      pro_psych_rep_load$WJIV_20_39$max_abs_diff),
                       "Dec_equ" = c(paf_psych_rep_load$baseline$are_equal,
                                     paf_psych_rep_load$case_1a$are_equal,
                                     paf_psych_rep_load$case_6b$are_equal,
                                     paf_psych_rep_load$case_11b$are_equal,
                                     paf_psych_rep_load$DOSPERT$are_equal,
                                     paf_psych_rep_load$IDS_2$are_equal,
                                     paf_psych_rep_load$WJIV_3_5$are_equal,
                                     paf_psych_rep_load$WJIV_20_39$are_equal,
                                     var_psych_rep_load$baseline$are_equal,
                                     var_psych_rep_load$case_1a$are_equal,
                                     var_psych_rep_load$case_6b$are_equal,
                                     var_psych_rep_load$case_11b$are_equal,
                                     var_psych_rep_load$DOSPERT$are_equal,
                                     var_psych_rep_load$IDS_2$are_equal,
                                     var_psych_rep_load$WJIV_3_5$are_equal,
                                     var_psych_rep_load$WJIV_20_39$are_equal,
                                     pro_psych_rep_load$baseline$are_equal,
                                     pro_psych_rep_load$case_1a$are_equal,
                                     pro_psych_rep_load$case_6b$are_equal,
                                     pro_psych_rep_load$case_11b$are_equal,
                                     pro_psych_rep_load$DOSPERT$are_equal,
                                     pro_psych_rep_load$IDS_2$are_equal,
                                     pro_psych_rep_load$WJIV_3_5$are_equal,
                                     pro_psych_rep_load$WJIV_20_39$are_equal),
                       "Max_dec" = c(paf_psych_rep_load$baseline$max_dec,
                                     paf_psych_rep_load$case_1a$max_dec,
                                     paf_psych_rep_load$case_6b$max_dec,
                                     paf_psych_rep_load$case_11b$max_dec,
                                     paf_psych_rep_load$DOSPERT$max_dec,
                                     paf_psych_rep_load$IDS_2$max_dec,
                                     paf_psych_rep_load$WJIV_3_5$max_dec,
                                     paf_psych_rep_load$WJIV_20_39$max_dec,
                                     var_psych_rep_load$baseline$max_dec,
                                     var_psych_rep_load$case_1a$max_dec,
                                     var_psych_rep_load$case_6b$max_dec,
                                     var_psych_rep_load$case_11b$max_dec,
                                     var_psych_rep_load$DOSPERT$max_dec,
                                     var_psych_rep_load$IDS_2$max_dec,
                                     var_psych_rep_load$WJIV_3_5$max_dec,
                                     var_psych_rep_load$WJIV_20_39$max_dec,
                                     pro_psych_rep_load$baseline$max_dec,
                                     pro_psych_rep_load$case_1a$max_dec,
                                     pro_psych_rep_load$case_6b$max_dec,
                                     pro_psych_rep_load$case_11b$max_dec,
                                     pro_psych_rep_load$DOSPERT$max_dec,
                                     pro_psych_rep_load$IDS_2$max_dec,
                                     pro_psych_rep_load$WJIV_3_5$max_dec,
                                     pro_psych_rep_load$WJIV_20_39$max_dec),
                       "SPSS_results" = rep("", 24),
                       "Md_diff" = c(paf_SPSS_rep_load$baseline$median_abs_diff,
                                     paf_SPSS_rep_load$case_1a$median_abs_diff,
                                     paf_SPSS_rep_load$case_6b$median_abs_diff,
                                     paf_SPSS_rep_load$case_11b$median_abs_diff,
                                     paf_SPSS_rep_load$DOSPERT$median_abs_diff,
                                     paf_SPSS_rep_load$IDS_2$median_abs_diff,
                                     paf_SPSS_rep_load$WJIV_3_5$median_abs_diff,
                                     paf_SPSS_rep_load$WJIV_20_39$median_abs_diff,
                                     var_SPSS_rep_load$baseline$median_abs_diff,
                                     var_SPSS_rep_load$case_1a$median_abs_diff,
                                     var_SPSS_rep_load$case_6b$median_abs_diff,
                                     var_SPSS_rep_load$case_11b$median_abs_diff,
                                     var_SPSS_rep_load$DOSPERT$median_abs_diff,
                                     var_SPSS_rep_load$IDS_2$median_abs_diff,
                                     var_SPSS_rep_load$WJIV_3_5$median_abs_diff,
                                     var_SPSS_rep_load$WJIV_20_39$median_abs_diff,
                                     pro_SPSS_rep_load$baseline$median_abs_diff,
                                     pro_SPSS_rep_load$case_1a$median_abs_diff,
                                     pro_SPSS_rep_load$case_6b$median_abs_diff,
                                     pro_SPSS_rep_load$case_11b$median_abs_diff,
                                     pro_SPSS_rep_load$DOSPERT$median_abs_diff,
                                     pro_SPSS_rep_load$IDS_2$median_abs_diff,
                                     pro_SPSS_rep_load$WJIV_3_5$median_abs_diff,
                                     pro_SPSS_rep_load$WJIV_20_39$median_abs_diff),
                       "Max_diff" = c(paf_SPSS_rep_load$baseline$max_abs_diff,
                                      paf_SPSS_rep_load$case_1a$max_abs_diff,
                                      paf_SPSS_rep_load$case_6b$max_abs_diff,
                                      paf_SPSS_rep_load$case_11b$max_abs_diff,
                                      paf_SPSS_rep_load$DOSPERT$max_abs_diff,
                                      paf_SPSS_rep_load$IDS_2$max_abs_diff,
                                      paf_SPSS_rep_load$WJIV_3_5$max_abs_diff,
                                      paf_SPSS_rep_load$WJIV_20_39$max_abs_diff,
                                      var_SPSS_rep_load$baseline$max_abs_diff,
                                      var_SPSS_rep_load$case_1a$max_abs_diff,
                                      var_SPSS_rep_load$case_6b$max_abs_diff,
                                      var_SPSS_rep_load$case_11b$max_abs_diff,
                                      var_SPSS_rep_load$DOSPERT$max_abs_diff,
                                      var_SPSS_rep_load$IDS_2$max_abs_diff,
                                      var_SPSS_rep_load$WJIV_3_5$max_abs_diff,
                                      var_SPSS_rep_load$WJIV_20_39$max_abs_diff,
                                      pro_SPSS_rep_load$baseline$max_abs_diff,
                                      pro_SPSS_rep_load$case_1a$max_abs_diff,
                                      pro_SPSS_rep_load$case_6b$max_abs_diff,
                                      pro_SPSS_rep_load$case_11b$max_abs_diff,
                                      pro_SPSS_rep_load$DOSPERT$max_abs_diff,
                                      pro_SPSS_rep_load$IDS_2$max_abs_diff,
                                      pro_SPSS_rep_load$WJIV_3_5$max_abs_diff,
                                      pro_SPSS_rep_load$WJIV_20_39$max_abs_diff),
                       "Dec_equ" = c(paf_SPSS_rep_load$baseline$are_equal,
                                     paf_SPSS_rep_load$case_1a$are_equal,
                                     paf_SPSS_rep_load$case_6b$are_equal,
                                     paf_SPSS_rep_load$case_11b$are_equal,
                                     paf_SPSS_rep_load$DOSPERT$are_equal,
                                     paf_SPSS_rep_load$IDS_2$are_equal,
                                     paf_SPSS_rep_load$WJIV_3_5$are_equal,
                                     paf_SPSS_rep_load$WJIV_20_39$are_equal,
                                     var_SPSS_rep_load$baseline$are_equal,
                                     var_SPSS_rep_load$case_1a$are_equal,
                                     var_SPSS_rep_load$case_6b$are_equal,
                                     var_SPSS_rep_load$case_11b$are_equal,
                                     var_SPSS_rep_load$DOSPERT$are_equal,
                                     var_SPSS_rep_load$IDS_2$are_equal,
                                     var_SPSS_rep_load$WJIV_3_5$are_equal,
                                     var_SPSS_rep_load$WJIV_20_39$are_equal,
                                     pro_SPSS_rep_load$baseline$are_equal,
                                     pro_SPSS_rep_load$case_1a$are_equal,
                                     pro_SPSS_rep_load$case_6b$are_equal,
                                     pro_SPSS_rep_load$case_11b$are_equal,
                                     pro_SPSS_rep_load$DOSPERT$are_equal,
                                     pro_SPSS_rep_load$IDS_2$are_equal,
                                     pro_SPSS_rep_load$WJIV_3_5$are_equal,
                                     pro_SPSS_rep_load$WJIV_20_39$are_equal),
                       "Max_dec" = c(paf_SPSS_rep_load$baseline$max_dec,
                                     paf_SPSS_rep_load$case_1a$max_dec,
                                     paf_SPSS_rep_load$case_6b$max_dec,
                                     paf_SPSS_rep_load$case_11b$max_dec,
                                     paf_SPSS_rep_load$DOSPERT$max_dec,
                                     paf_SPSS_rep_load$IDS_2$max_dec,
                                     paf_SPSS_rep_load$WJIV_3_5$max_dec,
                                     paf_SPSS_rep_load$WJIV_20_39$max_dec,
                                     var_SPSS_rep_load$baseline$max_dec,
                                     var_SPSS_rep_load$case_1a$max_dec,
                                     var_SPSS_rep_load$case_6b$max_dec,
                                     var_SPSS_rep_load$case_11b$max_dec,
                                     var_SPSS_rep_load$DOSPERT$max_dec,
                                     var_SPSS_rep_load$IDS_2$max_dec,
                                     var_SPSS_rep_load$WJIV_3_5$max_dec,
                                     var_SPSS_rep_load$WJIV_20_39$max_dec,
                                     pro_SPSS_rep_load$baseline$max_dec,
                                     pro_SPSS_rep_load$case_1a$max_dec,
                                     pro_SPSS_rep_load$case_6b$max_dec,
                                     pro_SPSS_rep_load$case_11b$max_dec,
                                     pro_SPSS_rep_load$DOSPERT$max_dec,
                                     pro_SPSS_rep_load$IDS_2$max_dec,
                                     pro_SPSS_rep_load$WJIV_3_5$max_dec,
                                     pro_SPSS_rep_load$WJIV_20_39$max_dec))

write_csv(Rep_orig, path = "output/rep_orig/Rep_orig.csv")

