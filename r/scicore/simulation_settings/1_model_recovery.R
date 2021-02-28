# Model recovery analysis
args <- commandArgs(TRUE)
eval(parse(text = args))
t1 <- Sys.time()

options(show.error.messages = FALSE, warn = -1)

library(EFAtools)
library(psych)
library(MASS)

# setwd("/Users/msteiner/Documents/Doktorat/efacomp_OSF/r/scicore/simulation_settings/")
sampsize <- 450


# ---------------------------
# Set up control variables
# ---------------------------

model_control <- expand.grid(
  case = names(population_models$loadings)[1:27],
  cors = names(population_models$phis_3),
  N = c(sampsize), # repeat with 450
  stringsAsFactors = FALSE)

settings <- expand.grid(
  comm_meth = c("unity", "mac", "smc"),  # initial communality estimation methods
  criterion_type = c("max_individual", "sums"), # citerion types
  abs_eigen = c(TRUE, FALSE), # absolute eigenvalues yes or no
  conv_crit = c(1e-3, 1e-6), # convergence criterion
  var_type = c("svd", "kaiser"), # varimax type
  p_type = c("unnorm", "norm"), # normalization of target matrix
  k = c(3, 4),
  stringsAsFactors = FALSE)

# maximum number of iterations to be used in PAF
paf_max_iter <- 1e4

# precision to use for varimax
precision = 1e-5

# how to order factors in promax
order_type = "eigen"

# define k according to earlier simulation studies
settings$k[settings$k == 3 & settings$p_type == "norm"] <- 2

# due to old R version on server
# isTRUE <- function(x) is.logical(x) && length(x) == 1L && !is.na(x) && x
# isFALSE <- function(x) is.logical(x) && length(x) == 1L && !is.na(x) && !x

# extract the loadings
case <- model_control$case[i]
temp_loadings <- population_models$loadings[[case]]
n_factors <- ncol(temp_loadings)

# extract the intercorrelations
cors <- model_control$cors[i]
temp_phi <- population_models[[paste0("phis_", n_factors)]][[cors]]

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
  
  set_list <- list()
  # run efa with the different settings
  for (set_i in 1:nrow(settings)) {
    
    comm_i <- settings$comm_meth[set_i]
    crit_type_i <- settings$criterion_type[set_i]
    conv_crit_i <- settings$conv_crit[set_i]
    abs_eigen_i <- settings$abs_eigen[set_i]
    var_type_i <- settings$var_type[set_i]
    p_type_i <- settings$p_type[set_i]
    k_i <- settings$k[set_i]
    
    # fit the models (with increased maximum iteration for psych and SPSS)
    efa_i <- try(EFA(sim_cor, n_factors = n_factors, type = "none", rotation = "promax",
                     init_comm = comm_i, criterion_type = crit_type_i,
                     criterion = conv_crit_i, abs_eigen = abs_eigen_i,
                     max_iter = paf_max_iter, use_cpp = TRUE,
                     signed_loadings = TRUE, P_type = p_type_i, k = k_i,
                     precision = precision, order_type = order_type,
                     varimax_type = var_type_i))
    
    if(class(efa_i) == "try-error" || max(abs(efa_i$h2)) >= .998 ||
       max(abs(efa_i$Structure)) >= .998){
      m_delta_load <- NA
      m_delta_cors <- NA
      diff_factor_corres <- NA
      diff_factor_corres_cross <- NA
      g <- NA
      if (class(efa_i) == "try-error") {
        heywood <- NA
        error <- efa_i[[1]]
        iter <- NA
      } else {
        heywood <- 1
        error <- NA
        iter <- efa_i$iter
      }
      
    } else {
      
      temp_i <- COMPARE(temp_loadings, efa_i$rot_loadings, thresh = .2)
      m_delta_load <- temp_i$mean_abs_diff
      g <- temp_i$g
      diff_factor_corres <- temp_i$diff_corres
      diff_factor_corres_cross <- temp_i$diff_corres_cross
      heywood <- 0
      error <- NA
      iter <- efa_i$iter
      
    }
    
    set_list[[set_i]] <- data.frame(
      case = case,
      cors = cors,
      N = N,
      n_factors = n_factors,
      g = g,
      m_delta_load = m_delta_load,
      heywood = heywood,
      diff_factor_corres = diff_factor_corres,
      diff_factor_corres_cross = diff_factor_corres_cross,
      comm_meth = comm_i,
      criterion_type = crit_type_i,
      abs_eigen = abs_eigen_i,
      conv_crit = conv_crit_i,
      var_type = var_type_i,
      p_type = p_type_i,
      k = k_i,
      iter = iter,
      error = error
    )
    
  }
  
  results_list[[sim]] <- do.call(rbind, set_list)
  
}

res <- do.call(rbind, results_list)

# save
saveRDS(res, file = paste0("output/mr/model_recovery_settings_results_",
                           sampsize, "_", i, ".RDS"))

ov_t2 <- Sys.time()
ov_t2 - t1
