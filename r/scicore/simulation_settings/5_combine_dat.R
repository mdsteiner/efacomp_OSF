library(tidyverse)

nam <- c("results_deviance_settings_180", "results_deviance_NE_settings_180",
         "results_heywood_settings_180", "results_heywood_NE_settings_180",
         "results_factor_corres_settings_180", "results_factor_corres_NE_settings_180",
         "results_deviance_settings_450", "results_deviance_NE_settings_450",
         "results_heywood_settings_450", "results_heywood_NE_settings_450",
         "results_factor_corres_settings_450", "results_factor_corres_NE_settings_450")

paths <- c("r/scicore/simulation_settings/output/ad/", "r/scicore/simulation_settings/output/ad/",
           "r/scicore/simulation_settings/output/ah/", "r/scicore/simulation_settings/output/ah/",
           "r/scicore/simulation_settings/output/af/", "r/scicore/simulation_settings/output/af/",
           "r/scicore/simulation_settings/output/ad/", "r/scicore/simulation_settings/output/ad/",
           "r/scicore/simulation_settings/output/ah/", "r/scicore/simulation_settings/output/ah/",
           "r/scicore/simulation_settings/output/af/", "r/scicore/simulation_settings/output/af/")

for (i in seq_along(nam)) {
  
  dat_temp <- do.call(rbind,
                      lapply(list.files(paths[i], pattern = nam[i], full.names = TRUE),
                             readRDS))
  
  saveRDS(dat_temp, paste0(nam[i], ".RDS"))
  
  
}
