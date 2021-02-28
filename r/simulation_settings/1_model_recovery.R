# Model recovery analysis

if(!require(EFAtools)) install.packages("EFAtools"); library(EFAtools)
if(!require(psych)) install.packages("psych"); library(psych)
if(!require(parallel)) install.packages("parallel"); library(parallel)
if(!require(MASS)) install.packages("MASS"); library(MASS)


setwd("/Users/msteiner/Documents/Doktorat/efacomp_OSF/")

source("r/model_recovery_function.R")

# ---------------------------
# Set up control variables
# ---------------------------

model_control <- expand.grid(
  case = names(population_models$loadings)[1:27],
  cors = names(population_models$phis_3),
  N = c(450), # repeat with 450
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

# ---------------------------
# Set up and run cluster
# ----------------------------

ov_t1 <- Sys.time()

# Set up cluster
cores_n <- detectCores()  # Number of cores to run on
cl <- makeCluster(cores_n)

# Send libraries to cluster
clusterEvalQ(cl, {
  library(EFAtools)
  library(psych)
  library(MASS)
})

# Export objects to cluster
clusterExport(cl, list("model_control", "settings", "paf_max_iter", "precision",
                       "order_type", "mean_abs_diff", "population_models", "g_rmse"))

# Run cluster!
cluster_result_ls <- parallel::parLapply(cl = cl,
                                         fun = model_recovery_settings,
                                         X = 1:nrow(model_control))

stopCluster(cl)  # Stop cluster

# ---------------------------
# Organise cluster results
# ----------------------------

# Convert to dataframe
cluster_result_df <- do.call(what = rbind,
                             args = cluster_result_ls)

# run twice with 500 each (append suffix  _2), otherwise the working memory won't suffice
saveRDS(cluster_result_df, file = "output/simulation_settings/model_recovery_settings_results_450.RDS")

ov_t2 <- Sys.time()
ov_t2 - ov_t1
