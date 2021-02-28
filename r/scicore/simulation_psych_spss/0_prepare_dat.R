fl <- list.files("simulation_psych_spss/output/mr/", full.names = TRUE)

dat <- do.call(rbind, lapply(fl, readRDS))

saveRDS(dat, "simulation_psych_spss/output/model_recovery_results.RDS")