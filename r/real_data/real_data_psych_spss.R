# Compare R psych vs. SPSS on real datasets ----

# Load packages where datasets are in
if (!require(psych)) install.packages("psych"); library(psych)
if (!require(EFAtools)) install.packages("EFAtools"); library(EFAtools)
future::plan(future::multisession) # to parallelize parallel analysis
if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if(!require(qgraph)) install.packages("qgraph"); library(qgraph)
if(!require(parallel)) install.packages("parallel"); library(parallel)
if(!require(stringr)) install.packages("stringr"); library(stringr)

# Load data --------

# WISC-V data (copy-right protected and therefore not made available)
WISCV_filenames <- list.files("data/wisc/", pattern = ".RData", 
                              full.names = TRUE)
for (i in WISCV_filenames) load(i)

# Human Cognitive Abilities project data (available at  http://www.iapsych.com/wmfhcaarchive/wmfhcaindex.html)
HCA_filenames <- list.files("data/HCA/", pattern = ".RData", 
                         full.names = TRUE)
for (i in HCA_filenames) load(i)

# Other datasets from R (in EFAtools structure)
other_filenames <- list.files("data/other_datasets/", pattern = ".RData", 
                          full.names = TRUE)
for (i in other_filenames) load(i)

# Get names of datasets 
WISCV_names <- str_remove_all(WISCV_filenames, c("data/wisc//|.RData"))

HCA_names <- str_remove_all(HCA_filenames, c("data/HCA//|.RData"))

other_names <- str_remove_all(other_filenames, 
                              c("data/other_datasets//|.RData"))

# Create vector of dataset names and number of theoretical factors -----

# Names of datasets
names_dat <- c(WISCV_names, "WJIV_ages_3_5", "WJIV_ages_6_8", "WJIV_ages_9_13", 
               "WJIV_ages_14_19", "WJIV_ages_20_39", "WJIV_ages_40_90", 
               HCA_names, other_names, "DOSPERT")

# Theoretical number of factors (number for HCA data from Carroll (1993); 
# NA means unknown) 
nfac_dat <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, # WISC-V
              8, 9, 9, 9, 9, 9, # WJ-IV
              3, 8, 10, 6, 6, 3, 3, 3, 3, 4, 5, NA, 3, 3, NA, 2, 4, 10, 3, 3, 3, 
              5, 8, 3, 4, 7, 9, 5, 8, 5, 4, 1, NA, 4, 4, 5, 8, 3, 2, 9, 4, 7, 2, 
              7, 2, 3, 5, 7, 2, 10, 8, 4, 6, 7, 5, 4, 4, 2, 4, 5, 3, 9, 8, 8, 4, 
              3, 1, 4, 3, 6, 5, NA, 7, 3, 6, 7, 2, 3, 7, 2, 4, 2, 5, 4, NA, 5, 4, 
              5, 4, 4, 4, 4, 3, 7, 3, 3, 7, 4, 5, 7, 5, 4, NA, 6, 8, NA, 3, 7, 4, 
              4, 2, NA, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 2, 2, 4, 4, 2, 5, 5, 7, 
              11, 3, 4, 5, 3, 6, 5, NA, NA, 3, 3, 2, 2, 4, 4, 9, 5, 5, 4, 2, 
              4, 3, 3, 10, 10, 3, 3, 3, 5, 7, 5, 6, 5, 5, 6, 8, 7, 6, 3, 4, 4, 
              NA, 2, 3, 3, 7, 7, 8, 11, 6, 3, 7, 2, 2, 3, 5, 8, 9, 1, 1, 1, 
              7, 4, 11, # HCA (Carroll, 1993)
              3, 6, 6, 5, 30, 4, 3, 1, 1, 5, 5, 5, 3, 1, 1, 2, 4, 5, 4, 7, 1,
              30, 45, 1, 1, 5, 4, 4, 4, 5, NA, 5, 4, 5, 5# other
              ) 
names(nfac_dat) <- names_dat

# Determine (initial) number of factors for PAF ----------

# # Perform a parallel analysis on each dataset (only run this once)
# nfac_pa <- c()
# tryerror <- c()
# 
# for (i in 1:length(names_dat)){
# 
#   temp_dat <- get(names_dat[i])
# 
#   temp_pa <- try(PARALLEL(temp_dat[[1]], N = temp_dat[[2]], eigen_type = "SMC",
#                           n_datasets = 1000, decision_rule = "crawford"),
#                  silent = TRUE)
# 
#   if(class(temp_pa) == "try-error"){
# 
#     nfac_pa[i] <- NA
#     tryerror[i] <- 1
# 
#   } else {
# 
#     nfac_pa[i] <- temp_pa$n_fac_SMC
#     tryerror[i] <- 0
# 
#   }
# 
# }
# 
# saveRDS(nfac_pa, file = "output/real_data/nfac_pa_real_data.RDS")
# saveRDS(tryerror, file = "output/real_data/tryerror_pa_real_data.RDS")

# Load number of factors from parallel analysis
nfac_pa <- readRDS("output/real_data/nfac_pa_real_data.RDS")

# If both the theoretical number of factors and the number of factors from PA are
# unknown, used the number of indicators divided by two to have an upper bound
# for the number of factors to extract
half_ind <- c()

for(i in 1:length(names_dat)){
  half_ind[i] <- round(nrow(get(names_dat[i])$cormat) / 2)
}

sum(is.na(nfac_dat) & is.na(nfac_pa)) # How often are both unknown?

nfac_dat <- ifelse(is.na(nfac_dat) & is.na(nfac_pa), half_ind, 
                   ifelse(is.na(nfac_dat) & !is.na(nfac_pa), 0, nfac_dat))

# Set NAs in nfac_pa to 0 for further computations
nfac_pa <- replace_na(nfac_pa, 0)
  
# Take the larger one of the theoretical and empirically determined number of factors
nfac_max <- ifelse(nfac_dat > nfac_pa, nfac_dat, nfac_pa)
names(nfac_max) <- names_dat

# Check if there are positive numbers only and no NAs
any(nfac_max <= 0)
sum(is.na(nfac_max))

# Prinicpal axis factor analyses with psych and SPSS settings ----------------

# Perform PAF with promax with type = psych and type = SPSS 
# (both with maximum 1e^6 iterations) for each dataset, starting from the determined 
# number of factors and continuing until a solution is found that works with 
# both types (convergence, no Heywood cases and at least two salient loadings 
# per factor)

fin_sol <- list()
res_psych_spss <- list()
fac_con_dat <- list()
diff_load <- list()

for(i in 1:length(names_dat)){

  temp_dat <- get(names_dat[i])
  
  # Set all variables for final dataset to NA
  nfac_admiss_psych <- NA
  nfac_admiss_spss <- NA
  min_salient_psych <- NA
  min_salient_spss <- NA
  nfac_final <- NA

  for (nfact in nfac_max[i]:1){

    # For each of these PAFs, first take smc as communality method, and if this     
    # doesn’t work, unity (for psych) or mac (for SPSS)

    # fit the three models (with increased maximum iteration for psych and SPSS)
    efa_psych <- try(EFA(temp_dat[[1]], n_factors = nfact, type = "psych", 
                         max_iter = 1e6, rotation = "promax", 
                         order_type = "eigen"), silent = TRUE)
    init_psych <- "smc"
    
    # if smc as initial communality estimates fails, rerun with unity
    if (class(efa_psych) == "try-error") {
      
      efa_psych <- try(EFA(temp_dat[[1]], n_factors = nfact, type = "psych", 
                           max_iter = 1e6, rotation = "promax", 
                           init_comm = "unity", order_type = "eigen"), 
                       silent = TRUE)
      init_psych <- "unity"
    }
    
    efa_spss <- try(EFA(temp_dat[[1]], n_factors = nfact, type = "SPSS", 
                        max_iter = 1e6, rotation = "promax", 
                        order_type = "eigen"), silent = TRUE)
    init_spss <- "smc"
    
    # if smc as initial communality estimates fails, rerun with mac
    if (class(efa_spss) == "try-error") {
      
      efa_spss <- try(EFA(temp_dat[[1]], n_factors = nfact, type = "SPSS", 
                          max_iter = 1e6, rotation = "promax", init_comm = "mac",
                          order_type = "eigen"), silent = TRUE)
      init_spss <- "mac"
    }
    
    if(class(efa_psych) != "try-error"){
      
    # Check if there are at least two salient loadings per factor
    salient_psych <- NULL
    for(j in 1:ncol(efa_psych$rot_loadings)) {
      salient_psych[j] <- sum(abs(efa_psych$rot_loadings[, j]) >= .20)
    }
    
    min_salient_psych <- ifelse(any(salient_psych < 2), FALSE, TRUE)
    
    }
    
    if(class(efa_spss) != "try-error"){
      
    salient_spss <- NULL
    for(j in 1:ncol(efa_spss$rot_loadings)){
      salient_spss[j] <- sum(abs(efa_spss$rot_loadings[, j]) >= .20)
    }
    
    min_salient_spss <- ifelse(any(salient_spss < 2), FALSE, TRUE)
    
    }
    
    # Save number of factors for first admissible solution in psych
    if(class(efa_psych) != "try-error" && min_salient_psych == TRUE && 
       all(efa_psych$h2 < .998) &&  all(abs(ifelse(nfact > 1, efa_psych$Structure,
                                                   efa_psych$unrot_loadings)) < .998) && 
       is.na(nfac_admiss_psych)){
        
        nfac_admiss_psych <- ncol(efa_psych$rot_loadings)
        
      }
    
    # Save number of factors for first admissible solution in SPSS
    if(class(efa_spss) != "try-error" && min_salient_spss == TRUE && 
       all(efa_spss$h2 < .998) && all(abs(ifelse(nfact > 1, efa_spss$Structure,
                                                  efa_spss$unrot_loadings)) < .998) && 
       is.na(nfac_admiss_spss)){
        
        nfac_admiss_spss <- ncol(efa_spss$rot_loadings)
        
      }
    
    # Proceed when first solution without try-error in both types is found
    if(class(efa_psych) != "try-error" & class(efa_spss) != "try-error") {
      
      if(all(efa_psych$h2 < .998) & all(efa_spss$h2 < .998) &
         all(abs(ifelse(nfact > 1, efa_psych$Structure, efa_psych$unrot_loadings)) 
             < .998) & all(abs(ifelse(nfact > 1, efa_spss$Structure,
                                      efa_spss$unrot_loadings)) < .998) &
         min_salient_psych == TRUE & min_salient_spss == TRUE) {
          
        nfac_final <- ncol(efa_psych$rot_loadings)
        
        if(is.na(nfac_admiss_psych)){
          nfac_admiss_psych <- nfac_final
        }
  
        if(is.na(nfac_admiss_spss)){
          nfac_admiss_spss <- nfac_final
        }
        
        break
        
      }
        
    }
  
  }

  diff_psych_spss_paf <- NA
  diff_psych_spss_var <- NA
  diff_psych_spss_pro <- NA
  m_diff_paf <- NA
  md_diff_paf <- NA
  min_diff_paf <- NA
  max_diff_paf <- NA
  m_diff_var <- NA
  md_diff_var <- NA
  min_diff_var <- NA
  max_diff_var <- NA
  m_diff_pro <- NA
  md_diff_pro <- NA
  min_diff_pro <- NA
  max_diff_pro <- NA
  diff_corres <- NA
  fac_con_paf <- NA
  fac_con_var <- NA
  fac_con_pro <- NA
  m_fac_con_paf <- NA
  m_fac_con_var <- NA
  m_fac_con_pro <- NA
  min_fac_con_paf <- NA
  min_fac_con_var <- NA
  min_fac_con_pro <- NA
  var_fac_ratio <- NA
  N_fac_ratio <- NA
  m_h2 <- NA
  m_fac_cor <- NA
  cross_psych <- NA
  cross_spss <- NA
  niter_psych <- NA
  niter_spss <- NA

  # After final solutions are found --------
  if(class(efa_psych) != "try-error" & class(efa_spss) != "try-error") {
    
    # Calculate varimax solution
    var_psych <- EFAtools:::.VARIMAX(efa_psych, type = "psych")
    var_spss <- EFAtools:::.VARIMAX(efa_spss, type = "SPSS")
  
    # Calculate differences in loadings and variable-to-factor-correspondences
    diff_psych_spss_paf <- COMPARE(efa_psych$unrot_loadings, efa_spss$unrot_loadings,
                                   thresh = .20)
    diff_psych_spss_var <- COMPARE(var_psych$rot_loadings, var_spss$rot_loadings,
                                   thresh = .20)
    diff_psych_spss_pro <- COMPARE(efa_psych$rot_loadings, efa_spss$rot_loadings,
                                   thresh = .20)

    m_diff_paf <- diff_psych_spss_paf$mean_abs_diff
    md_diff_paf <- diff_psych_spss_paf$median_abs_diff
    min_diff_paf <- diff_psych_spss_paf$min_abs_diff
    max_diff_paf <- diff_psych_spss_paf$max_abs_diff
    
    m_diff_var <- diff_psych_spss_var$mean_abs_diff
    md_diff_var <- diff_psych_spss_var$median_abs_diff
    min_diff_var <- diff_psych_spss_var$min_abs_diff
    max_diff_var <- diff_psych_spss_var$max_abs_diff
    
    m_diff_pro <- diff_psych_spss_pro$mean_abs_diff
    md_diff_pro <- diff_psych_spss_pro$median_abs_diff
    min_diff_pro <- diff_psych_spss_pro$min_abs_diff
    max_diff_pro <- diff_psych_spss_pro$max_abs_diff
    
    diff_corres <- diff_psych_spss_pro$diff_corres
    
    diff_paf_vec <- as.vector(diff_psych_spss_paf$diff)
    diff_var_vec <- as.vector(diff_psych_spss_var$diff)
    diff_pro_vec <- as.vector(diff_psych_spss_pro$diff)
    
    # Calculate factor congruence
    fac_con_paf <- apply(EFAtools:::.factor_congruence(efa_psych$unrot_loadings,
                                                      efa_spss$unrot_loadings),
                         1, max)
    fac_con_var <- apply(EFAtools:::.factor_congruence(var_psych$rot_loadings,
                                                      var_spss$rot_loadings),
                         1, max)
    fac_con_pro <- apply(EFAtools:::.factor_congruence(efa_psych$rot_loadings, 
                                                      efa_spss$rot_loadings),
                         1, max)
    
    m_fac_con_paf <- mean(fac_con_paf)
    m_fac_con_var <- mean(fac_con_var)
    m_fac_con_pro <- mean(fac_con_pro)
    min_fac_con_paf <- min(fac_con_paf)
    min_fac_con_var <- min(fac_con_var)
    min_fac_con_pro <- min(fac_con_pro)
    
    # Calculate variable-to-factor ratio, N-to-fac ratio, mean communality, and
    # mean factor intercorrelation
    var_fac_ratio <- nrow(temp_dat[[1]]) / nfac_final
    N_fac_ratio <- temp_dat[[2]] / nfac_final
    m_h2 <- mean(c(efa_psych$h2, efa_spss$h2))
    m_fac_cor <- mean(c(efa_psych$Phi[lower.tri(efa_psych$Phi)], 
                        efa_spss$Phi[lower.tri(efa_spss$Phi)]))
    
    # Count of cross-loadings for each solution
    temp_cross_psych <- rowSums(abs(efa_psych$rot_loadings) >= .20)
    cross_psych <- sum(ifelse(temp_cross_psych > 0, temp_cross_psych - 1, 
                              temp_cross_psych))
    
    temp_cross_spss <- rowSums(abs(efa_spss$rot_loadings) >= .20)
    cross_spss <- sum(ifelse(temp_cross_spss > 0, temp_cross_spss - 1, 
                              temp_cross_spss))
    
    # Number of iterations for Psych and SPSS
    niter_psych <- efa_psych$iter
    niter_spss <- efa_spss$iter

  }
  
  
  if(class(efa_psych) == "try-error") {
    
    efa_psych <- NA
    
  }
  
  if(class(efa_spss) == "try-error") {
    
    efa_spss <- NA
    
  }
  
  # Create dataframe for loadings differences 
  diff_load[[i]] <- data.frame(datname = rep(names_dat[i], 
                                             (length(diff_paf_vec)*3 + 6)),
                               values = c(diff_paf_vec, m_diff_paf,
                                          max_diff_paf,
                                          diff_var_vec, m_diff_var,
                                          max_diff_var,
                                          diff_pro_vec, m_diff_pro,
                                          max_diff_pro),
                               procedure = c(rep("PAF", (length(diff_paf_vec) 
                                                         + 2)),
                                             rep("Varimax", (length(diff_var_vec) 
                                                             + 2)),
                                             rep("Promax", (length(diff_pro_vec) 
                                                            + 2))),
                               group = rep(c(rep("Overall", length(diff_paf_vec)),
                                             "Mean", "Maximum"), 3))

  # Create dataframe for factor congruences 
  fac_con_dat[[i]] <- data.frame(datname = rep(names_dat[i], 
                                               (length(fac_con_paf)*3 + 6)),
                                 values = c(fac_con_paf, m_fac_con_paf,
                                            min_fac_con_paf,
                                            fac_con_var, m_fac_con_var,
                                            min_fac_con_var,
                                            fac_con_pro, m_fac_con_pro,
                                            min_fac_con_pro),
                                 procedure = c(rep("PAF", (length(fac_con_paf) 
                                                           + 2)),
                                               rep("Varimax", (length(fac_con_var) 
                                                               + 2)),
                                               rep("Promax", (length(fac_con_var) 
                                                              + 2))),
                                 group = rep(c(rep("Overall", length(fac_con_paf)),
                                               "Mean", "Minimum"), 3))
  
  # Create a list for each dataset (names of datasets as names) with
  fin_sol[[names_dat[i]]] <- list(psych = efa_psych, 
                                  spss = efa_spss, 
                                  diff_psych_spss_paf = diff_psych_spss_paf,
                                  diff_psych_spss_var = diff_psych_spss_var,
                                  diff_psych_spss_pro = diff_psych_spss_pro)

  # Create dataframe ------------
  
  # diff_corres: number of variables with different factor correspondence
  res_psych_spss[[i]] <- data.frame(
                              datname = names_dat[i],
                              N = temp_dat[[2]],
                              nfac_theor = nfac_dat[i],
                              nfac_pa = nfac_pa[i],
                              nfac_admiss_psych = nfac_admiss_psych,
                              nfac_admiss_spss = nfac_admiss_spss,
                              nfac_final = nfac_final,
                              nind = nrow(temp_dat[[1]]),
                              var_fac_ratio = var_fac_ratio,
                              N_fac_ratio = N_fac_ratio,
                              m_h2 = m_h2,
                              m_fac_cor = m_fac_cor,
                              comm_meth_psych = init_psych,
                              comm_meth_spss = init_spss,
                              niter_psych = niter_psych,
                              niter_spss = niter_spss,
                              cross_psych = cross_psych,
                              cross_spss = cross_spss,
                              diff_corres = diff_corres,
                              m_diff_paf = m_diff_paf,
                              m_diff_var = m_diff_var,
                              m_diff_pro = m_diff_pro,
                              md_diff_paf = md_diff_paf,
                              md_diff_var = md_diff_var,
                              md_diff_pro = md_diff_pro,
                              min_diff_paf = min_diff_paf,
                              min_diff_var = min_diff_var,
                              min_diff_pro = min_diff_pro,
                              max_diff_paf = max_diff_paf,
                              max_diff_var = max_diff_var,
                              max_diff_pro = max_diff_pro,
                              m_fac_con_paf = m_fac_con_paf,
                              m_fac_con_var = m_fac_con_var,
                              m_fac_con_pro = m_fac_con_pro,
                              min_fac_con_paf = min_fac_con_paf,
                              min_fac_con_var = min_fac_con_var,
                              min_fac_con_pro = min_fac_con_pro)

}

real_data_results <- do.call("rbind", res_psych_spss)
fac_con_dat <- do.call("rbind", fac_con_dat)
diff_load <- do.call("rbind", diff_load)
  
# Save output -----
saveRDS(real_data_results, file = "output/real_data/real_data_results.RDS")
saveRDS(fac_con_dat, file = "output/real_data/fac_con_dat.RDS")
saveRDS(fin_sol, file = "output/real_data/real_data_fin_sol.RDS")
saveRDS(diff_load, file = "output/real_data/diff_load.RDS")