# Model recovery analysis

if(!require(EFAtools)) install.packages("EFAtools"); library(EFAtools)
if(!require(psych)) install.packages("psych"); library(psych)
if(!require(parallel)) install.packages("parallel"); library(parallel)
if(!require(MASS)) install.packages("MASS"); library(MASS)

setwd("/Users/msteiner/Documents/Doktorat/AJJ_OSF")

source("r/model_recovery_function.R")

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
clusterExport(cl, list("model_control", "mean_abs_diff", "population_models",
                       "paf_max_iter", "abs_diff", "mean_abs_diff", "g_rmse",
                       "paf_max_iter")) # , "isTRUE", "isFALSE" # if on server with old R version

# Run cluster!
cluster_result_ls <- parallel::parLapply(cl = cl,
                                         fun = model_recovery,
                                         X = 1:nrow(model_control))

stopCluster(cl)  # Stop cluster

# ---------------------------
# Organise cluster results
# ----------------------------

# Convert to dataframe
cluster_result_df <- do.call(what = rbind,
                             args = cluster_result_ls)

saveRDS(cluster_result_df, file = "output/simulation_psych_spss/model_recovery_results.RDS")

ov_t2 <- Sys.time()
ov_t2 - ov_t1