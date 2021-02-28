# Model recovery analysis
args <- commandArgs(TRUE)
eval(parse(text = args))
t1 <- Sys.time()

options(show.error.messages = FALSE, warn = -1)

library(EFAtools)
library(psych)
library(MASS)

# setwd("/Users/msteiner/Documents/Doktorat/efacomp_OSF/r/scicore/simulation_psych_spss/")

# ---------------------------
# Set up control variables
# ---------------------------

model_control <- expand.grid(
  case = names(population_models$loadings)[1:27],
  cors = names(population_models$phis_3),
  N = c(180, 450),
  stringsAsFactors = FALSE)

# maximum number of iterations to be used in PAF
paf_max_iter <- 5e4

# if on server with old R version
# isTRUE <- function(x) is.logical(x) && length(x) == 1L && !is.na(x) && x
# isFALSE <- function(x) is.logical(x) && length(x) == 1L && !is.na(x) && !x


# ---------------------------
# Set up and run cluster
# ----------------------------

# extract the loadings
case <- model_control$case[i]
temp_loadings <- population_models$loadings[[case]]
n_factors <- ncol(temp_loadings)

# extract the intercorrelations
cors <- model_control$cors[i]
temp_phi <- population_models[[paste0("phis_", n_factors)]][[cors]]

# compute variance proportions
vars <- diag(temp_phi %*% t(temp_loadings) %*% temp_loadings)
var_total <- nrow(temp_loadings)
true_var_expl <- sum(vars / var_total)

# extract N
N <- model_control$N[i]

# set up vectors with means and sds
mus <- rep(0, nrow(temp_loadings))
sds <- rep(1, nrow(temp_loadings))

# get correlation matrix from loadings and factor intercorrelations
cormat <- temp_loadings %*% temp_phi %*% t(temp_loadings)
diag(cormat) <- 1

# create covariance matrix for multivariate distribution
covmat <- sds %*% t(sds) * cormat

results_list <- list()

for (sim in 1:1000) {
  
  # simulate data
  sim_dat <- mvrnorm(n = N, mu = mus, Sigma = covmat)
  
  # obtain correlation matrix
  sim_cor <- cor(sim_dat)
  
  # init objects
  nfac_admiss_spss <- NA
  nfac_admiss_psych <- NA
  
  # get first admissible solution for psych and for spss
  for (nfact_i in (n_factors + 4):1){
    
    # For each of these PAFs, first take smc as communality method, and if this doesn’t work, unity (for psych) or mac (for SPSS)
    
    # fit the three models (with increased maximum iteration for psych and SPSS)
    efa_psych <- try(EFA(sim_cor, n_factors = nfact_i, type = "psych",
                         max_iter = paf_max_iter, rotation = "promax"),
                     silent = TRUE)
    
    # if smc as initial communality estimates fails, rerun with unity
    if (class(efa_psych) == "try-error") {
      
      efa_psych <- try(EFA(sim_cor, n_factors = nfact_i, type = "psych",
                           max_iter = paf_max_iter, rotation = "promax",
                           init_comm = "unity"), 
                       silent = TRUE)
    }
    
    efa_spss <- try(EFA(sim_cor, n_factors = nfact_i, type = "SPSS",
                        max_iter = paf_max_iter, rotation = "promax"),
                    silent = TRUE)
    
    # if smc as initial communality estimates fails, rerun with mac
    if (class(efa_spss) == "try-error") {
      
      efa_spss <- try(EFA(sim_cor, n_factors = nfact_i, type = "SPSS",
                          max_iter = paf_max_iter, rotation = "promax",
                          init_comm = "mac"), silent = TRUE)
    }
    
    if (class(efa_psych) != "try-error") {
      
      # Check if there are at least two salient loadings per factor
      salient_psych <- ifelse(nfact_i > 1,
                              colSums(abs(efa_psych$rot_loadings) >= .20),
                              colSums(abs(efa_psych$unrot_loadings) >= .20))
      
      min_salient_psych <- ifelse(any(salient_psych < 2), FALSE, TRUE)
      
    } else {
      min_salient_psych <- FALSE
    }
    
    if (class(efa_spss) != "try-error") {
      salient_spss <- ifelse(nfact_i > 1, colSums(abs(efa_spss$rot_loadings) >= .20),
                             colSums(abs(efa_spss$unrot_loadings) >= .20))
      
      min_salient_spss <- ifelse(any(salient_spss < 2), FALSE, TRUE)
    } else {
      min_salient_spss <- FALSE
    }
    
    # Save number of factors for first admissible solution in psych
    if(class(efa_psych) != "try-error" && min_salient_psych == TRUE && 
       all(efa_psych$h2 < .998) &&
       all(abs(ifelse(nfact_i > 1, efa_psych$Structure, efa_psych$unrot_loadings)) < .998) && 
       is.na(nfac_admiss_psych)){
      
      nfac_admiss_psych <- ncol(efa_psych$unrot_loadings)
      
    }
    
    # Save number of factors for first admissible solution in SPSS
    if(class(efa_spss) != "try-error" && min_salient_spss == TRUE && 
       all(efa_spss$h2 < .998) &&
       all(abs(ifelse(nfact_i > 1, efa_spss$Structure, efa_spss$unrot_loadings)) < .998) && 
       is.na(nfac_admiss_spss)){
      
      nfac_admiss_spss <- ncol(efa_spss$unrot_loadings)
      
    }
    
    # break when a solution was found for both types
    if(!is.na(nfac_admiss_psych) && !is.na(nfac_admiss_spss)) {
      
      break
      
    }
    
  }
  
  # fit the three models (with increased maximum iteration for psych and SPSS)
  efa_psych <- try(EFA(sim_cor, n_factors = n_factors, type = "psych",
                       max_iter = paf_max_iter, rotation = "promax"))
  init_psych <- "smc"
  
  # if smc as initial communality estimates fails, rerun with unity
  if (class(efa_psych) == "try-error") {
    
    efa_psych <- try(EFA(sim_cor, n_factors = n_factors, type = "psych",
                         max_iter = paf_max_iter, rotation = "promax",
                         init_comm = "unity"))
    init_psych <- "unity"
  }
  
  efa_spss <- try(EFA(sim_cor, n_factors = n_factors, type = "SPSS",
                      max_iter = paf_max_iter, rotation = "promax"))
  init_spss <- "smc"
  
  # if smc as initial communality estimates fails, rerun with mac
  if (class(efa_spss) == "try-error") {
    
    efa_spss <- try(EFA(sim_cor, n_factors = n_factors, type = "SPSS",
                        max_iter = paf_max_iter, rotation = "promax",
                        init_comm = "mac"))
    init_spss <- "mac"
  }
  
  heywood_spss <- any(efa_spss$h2 >= .998) || any(abs(ifelse(nfact_i > 1, efa_spss$Structure, efa_spss$unrot_loadings)) >= .998)
  heywood_psych <- any(efa_psych$h2 >= .998) || any(abs(ifelse(nfact_i > 1, efa_psych$Structure, efa_psych$unrot_loadings)) >= .998)
  
  if(class(efa_psych) == "try-error" || heywood_psych){
    delta_var_expl_psych <- NA
    g_psych <- NA
    m_delta_load_psych <- NA
    diff_factor_corres_psych <- NA
    diff_factor_corres_cross_psych <- NA
    if (class(efa_psych) == "try-error") {
      heywood_psych <- NA
      error_psych <- efa_psych[[1]]
      iter_psych <- NA
    } else {
      heywood_psych <- 1
      error_psych <- NA
      iter_psych <- efa_psych$iter
    }
    
  } else {
    temp_psych <- COMPARE(temp_loadings, efa_psych$rot_load, thresh = .2)
    delta_var_expl_psych <- efa_psych$vars_accounted[3, n_factors] - true_var_expl
    g_psych <- temp_psych$g 
    m_delta_load_psych <- temp_psych$mean_abs_diff
    diff_factor_corres_psych <- temp_psych$diff_corres
    diff_factor_corres_cross_psych <- temp_psych$diff_corres_cross
    heywood_psych <- 0
    error_psych <- NA
    iter_psych <- efa_psych$iter
  }
  
  if(class(efa_spss) == "try-error" || heywood_spss){
    delta_var_expl_spss <- NA
    g_spss <- NA
    m_delta_load_spss <- NA
    diff_factor_corres_spss <- NA
    diff_factor_corres_cross_spss <- NA
    if (class(efa_psych) == "try-error") {
      heywood_spss <- NA
      error_spss <- efa_spss[[1]]
      iter_spss <- NA
    } else {
      heywood_spss <- 1
      error_spss <- NA
      iter_spss <- efa_spss$iter
    }
  } else {
    temp_spss <- COMPARE(temp_loadings, efa_spss$rot_load, thresh = .2)
    delta_var_expl_spss <- efa_spss$vars_accounted[3, n_factors] - true_var_expl
    g_spss <- temp_spss$g
    m_delta_load_spss <- temp_spss$mean_abs_diff
    diff_factor_corres_spss <- temp_spss$diff_corres
    diff_factor_corres_cross_spss <- temp_spss$diff_corres_cross
    heywood_spss <- 0
    error_spss <- NA
    iter_spss <- efa_spss$iter
  }
  
  # check difference between psych and SPSS fit
  if(class(efa_spss) != "try-error" && !heywood_spss &&
     class(efa_psych) != "try-error" && !heywood_psych) {
    
    # store differences between spss and psych result
    temp_diff <- COMPARE(efa_spss$rot_load, efa_psych$rot_load, thresh = .2)
    m_delta_load_spss_psych <- temp_diff$mean_abs_diff
    min_delta_load_spss_psych <- temp_diff$min_abs_diff
    max_delta_load_spss_psych <- temp_diff$max_abs_diff
    string_delta_load_spss_psych <- paste(c(temp_diff$diff), collapse = "; ")
    diff_factor_corres_spss_psych <- temp_diff$diff_corres
    diff_factor_corres_cross_spss_psych <- temp_diff$diff_corres_cross
  } else {
    
    m_delta_load_spss_psych <- NA
    min_delta_load_spss_psych <- NA
    max_delta_load_spss_psych <- NA
    string_delta_load_spss_psych <- NA
    diff_factor_corres_spss_psych <- NA
    diff_factor_corres_cross_spss_psych <- NA
    
  }
  
  results_list[[sim]] <- data.frame(
    case = case,
    cors = cors,
    N = N,
    true_var_expl = true_var_expl,
    delta_var_expl_psych = delta_var_expl_psych,
    delta_var_expl_spss = delta_var_expl_spss,
    g_psych = g_psych,
    g_spss = g_spss,
    m_delta_load_psych = m_delta_load_psych,
    m_delta_load_spss = m_delta_load_spss,
    heywood_psych = heywood_psych,
    heywood_spss = heywood_spss,
    error_psych = error_psych,
    error_spss = error_spss,
    iter_psych = iter_psych,
    iter_spss = iter_spss,
    diff_factor_corres_psych = diff_factor_corres_psych,
    diff_factor_corres_spss = diff_factor_corres_spss,
    diff_factor_corres_cross_psych = diff_factor_corres_cross_psych,
    diff_factor_corres_cross_spss = diff_factor_corres_cross_spss,
    init_psych = init_psych,
    init_spss = init_spss,
    m_delta_load_spss_psych = m_delta_load_spss_psych,
    min_delta_load_spss_psych = min_delta_load_spss_psych,
    max_delta_load_spss_psych = max_delta_load_spss_psych,
    string_delta_load_spss_psych = string_delta_load_spss_psych,
    diff_factor_corres_spss_psych = diff_factor_corres_spss_psych,
    diff_factor_corres_cross_spss_psych = diff_factor_corres_cross_spss_psych,
    nfac_admiss_spss = nfac_admiss_spss,
    nfac_admiss_psych = nfac_admiss_psych
  )
  
}

res <- do.call(rbind, results_list)

# save
saveRDS(res, file = paste0("output/mr/model_recovery_psych_spss_results_",
                           i, ".RDS"))

ov_t2 <- Sys.time()
ov_t2 - t1
