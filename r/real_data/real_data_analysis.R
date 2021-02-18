# Analysis of comparisons of R psych vs. SPSS with real data sets------

# Load packages
if (!require(psych)) install.packages("psych"); library(psych)
if (!require(devtools)) install.packages("devtools"); library(devtools)
if (!require(EFAtools)) install.packages("EFAtools"); library(EFAtools)
if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if (!require(see)) install.packages("see"); library(see)
if (!require(viridis)) install.packages("viridis"); library(viridis)
if (!require(rstanarm)) install.packages("rstanarm"); library(rstanarm)
options(mc.cores = parallel::detectCores())
if (!require(bayestestR)) install.packages("bayestestR"); library(bayestestR)
if (!require(patchwork)) install.packages("bayestestR"); library(patchwork)

# Load results from comparisons --------
real_dat_res <- readRDS("output/real_data/real_data_results.RDS")
fin_sol <- readRDS("output/real_data/real_data_fin_sol.RDS")
fac_con_dat <- readRDS("output/real_data/fac_con_dat.RDS")
diff_load <- readRDS("output/real_data/diff_load.RDS")

# First check if a final solution was found for all datasets ----
no_sol <- which(is.na(real_dat_res$nfac_final)); no_sol

real_dat_res$datname[no_sol]

# See if this differs for psych and SPSS ---
no_sol_psych <- which(is.na(real_dat_res$nfac_admiss_psych)); no_sol_psych
no_sol_spss <- which(is.na(real_dat_res$nfac_admiss_spss)); no_sol_spss

# Delete cases for which no solution was found ----
fin_sol <- fin_sol[!(names(fin_sol) %in% real_dat_res$datname[no_sol])]
real_dat_res <- real_dat_res[-no_sol, ]

### Analysis of differences in resulting solutions from R psych and SPSS -----

# How often does the first admissible solution differ for R psych vs. SPSS?----
sum(real_dat_res$nfac_admiss_psych != 
      real_dat_res$nfac_admiss_spss) / nrow(real_dat_res)

# Does one systematically result in higher number of factor?----
diff_admiss <- real_dat_res$nfac_admiss_psych - real_dat_res$nfac_admiss_spss

mean(diff_admiss)
table(diff_admiss)

diff_admiss_spss <- sum(diff_admiss < 0) / nrow(real_dat_res)
diff_admiss_psych <- sum(diff_admiss > 0) / nrow(real_dat_res)

diff_admiss_spss / sum(diff_admiss_spss, diff_admiss_psych)
diff_admiss_psych / sum(diff_admiss_spss, diff_admiss_psych)

# Do communality methods sometimes differ between psych and SPSS?
real_dat_res[real_dat_res$comm_meth_psych != real_dat_res$comm_meth_spss, ]
# -> No, no further analyses here.


# Differences in variable-to-factor correspondences -----

# Absolute number and relative to all datasets with admissible solutions
sum(real_dat_res$diff_corres > 0) / nrow(real_dat_res)

table(real_dat_res$diff_corres) / nrow(real_dat_res)

# Include only solutions with minimum 2 factors
nrow(real_dat_res[real_dat_res$nfac_final > 1, ]) # How many solutions?

sum(real_dat_res$diff_corres[real_dat_res$nfac_final > 1] > 0) / nrow(real_dat_res[real_dat_res$nfac_final > 1, ])

table(real_dat_res$diff_corres[real_dat_res$nfac_final > 1]) / nrow(real_dat_res[real_dat_res$nfac_final > 1, ])

# Relative to total number of indicators (= proportion of indicators with 
# different factor correspondences per dataset)
hist(real_dat_res$diff_corres / real_dat_res$nind)


# Crossloadings ------

# Psych
hist(real_dat_res$cross_psych, breaks = 100)
sum(real_dat_res$cross_psych > 0) / nrow(real_dat_res)
hist(real_dat_res$cross_psych / real_dat_res$nind)

# SPSS
hist(real_dat_res$cross_spss, breaks = 100)
sum(real_dat_res$cross_spss > 0) / nrow(real_dat_res)
hist(real_dat_res$cross_spss / real_dat_res$nind)


# Differences in loadings ----- 

# PAF (unrotated loadings)
paf_load_diff <- real_dat_res %>% 
  select(md_diff_paf, m_diff_paf, max_diff_paf) %>% 
  summarize(mean_md = mean(md_diff_paf),
            mean_m = mean(m_diff_paf),
            sd_m = sd(m_diff_paf),
            min_m = min(m_diff_paf),
            max_m = max(m_diff_paf),
            mean_max = mean(max_diff_paf),
            sd_max = sd(max_diff_paf),
            min_max = min(max_diff_paf),
            max_max = max(max_diff_paf))

min(abs(diff_load$values[diff_load$procedure == "PAF"]), na.rm = TRUE)

var_load_diff <- real_dat_res %>% 
  select(md_diff_var, m_diff_var, max_diff_var) %>% 
  summarize(mean_md = mean(md_diff_var),
            mean_m = mean(m_diff_var),
            sd_m = sd(m_diff_var),
            min_m = min(m_diff_var),
            max_m = max(m_diff_var),
            mean_max = mean(max_diff_var),
            sd_max = sd(max_diff_var),
            min_max = min(max_diff_var),
            max_max = max(max_diff_var))

min(abs(diff_load$values[diff_load$procedure == "Varimax"]), na.rm = TRUE)

pro_load_diff <- real_dat_res %>% 
  select(md_diff_pro, m_diff_pro, max_diff_pro) %>% 
  summarize(mean_md = mean(md_diff_pro),
            mean_m = mean(m_diff_pro),
            sd_m = sd(m_diff_pro),
            min_m = min(m_diff_pro),
            max_m = max(m_diff_pro),
            mean_max = mean(max_diff_pro),
            sd_max = sd(max_diff_pro),
            min_max = min(max_diff_pro),
            max_max = max(max_diff_pro))

min(abs(diff_load$values[diff_load$procedure == "Promax"]), na.rm = TRUE)

real_dat_res$datname[which(real_dat_res$max_diff_pro > .3)]


# Plot differences in loadings

# Median differences
hist(real_dat_res$md_diff_pro, breaks = 100)

# Mean and maximum differences
diff_plot <- diff_load %>%
  mutate(values = abs(values),
         statistic = factor(group, levels = c("Overall", "Mean", "Maximum")),
         procedure = factor(procedure, levels = c("PAF", "Varimax", "Promax"))) %>%
  ggplot(aes(x = statistic, y = values, fill = "A")) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5, color = viridis(1)) +
  geom_violin(size = 0.5, color = viridis(1), fill = alpha(viridis(1), 0.3)) +
  scale_fill_manual(values = viridis(1, begin = .2, alpha = .7)) +
  scale_y_continuous(breaks = seq(0, 0.7, 0.05)) +
  facet_wrap(vars(procedure), ncol = 3) +
  stat_summary(fun = "mean", geom = "point", shape = 23, fill = "white") +
  labs(x = "", y = "Differences in Loadings") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major.x = element_blank())

## Factor congruence ----

# Calculate mean and min factor congruence
fac_con_paf <- real_dat_res %>% 
  select(m_fac_con_paf, min_fac_con_paf) %>% 
  summarize(mean_m = mean(m_fac_con_paf),
            min_m = min(m_fac_con_paf),
            mean_min = mean(min_fac_con_paf),
            min_min = min(min_fac_con_paf))

max(abs(fac_con_dat$values[fac_con_dat$procedure == "PAF"]), na.rm = TRUE)

fac_con_var <- real_dat_res %>% 
  select(m_fac_con_var, min_fac_con_var) %>% 
  summarize(mean_m = mean(m_fac_con_var),
            min_m = min(m_fac_con_var),
            mean_min = mean(min_fac_con_var),
            min_min = min(min_fac_con_var))

max(abs(fac_con_dat$values[fac_con_dat$procedure == "Varimax"]), na.rm = TRUE)

fac_con_pro <- real_dat_res %>% 
  select(m_fac_con_pro, min_fac_con_pro) %>% 
  summarize(mean_m = mean(m_fac_con_pro),
            min_m = min(m_fac_con_pro),
            mean_min = mean(min_fac_con_pro),
            min_min = min(min_fac_con_pro))

max(abs(fac_con_dat$values[fac_con_dat$procedure == "Promax"]), na.rm = TRUE)


fac_con_plot <- fac_con_dat %>%
  mutate(statistic = factor(group, levels = c("Overall", "Mean", "Minimum")),
         procedure = factor(procedure, levels = c("PAF", "Varimax", "Promax"))) %>%
  ggplot(aes(x = statistic, y = values, fill = "A")) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5, color = viridis(1)) +
  geom_violin(size = 0.5, color = viridis(1), fill = alpha(viridis(1), 0.3)) +
  scale_fill_manual(values = viridis(1, begin = .2, alpha = .7)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  facet_wrap(vars(procedure), ncol = 3) +
  stat_summary(fun = "mean", geom = "point", shape = 23, fill = "white") +
  labs(x = "", y = "Factor Congruence") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major.x = element_blank())

# Combine plots
pdf("plots/diffload_faccon_real_data.pdf", height = 7, width = 8)
diff_plot / fac_con_plot
dev.off()

# Differences only for those with variable-to-factor ratio < .3

# Mean differences
ggplot(real_dat_res[real_dat_res$var_fac_ratio >= 3, ], aes(x = m_diff_pro)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), binwidth = .001) +
  geom_segment(x = mean(real_dat_res$m_diff_pro), xend = mean(real_dat_res$m_diff_pro),
               y = 0, yend = .2, col = "red", linetype = 2) +
  scale_x_continuous(name = "Mean Differences in Loadings",
                     breaks = seq(0, .2, .01)) +
  scale_y_continuous(name = "Proportion of Data Sets",
                     breaks = seq(0, .2, .01)) +
  theme_bw()

# Maximum differences
ggplot(real_dat_res[real_dat_res$var_fac_ratio >= 3, ], aes(x = max_diff_pro)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), binwidth = .003) +
  geom_segment(x = mean(real_dat_res$max_diff_pro), xend = mean(real_dat_res$max_diff_pro),
               y = 0, yend = .2, col = "red", linetype = 2) +
  scale_x_continuous(name = "Maximum Differences in Loadings",
                     breaks = seq(0, .5, .05)) +
  scale_y_continuous(name = "Proportion of Data Sets",
                     breaks = seq(0, .5, .01)) +
  theme_bw()

# See if differences are related to some properties of the data sets -----

# Preliminary regression analyses

# Variable-to-factor ratio

# Mean differences
summary(lm(m_diff_pro ~ var_fac_ratio, real_dat_res))
plot(lm(m_diff_pro ~ var_fac_ratio, real_dat_res))
cor(real_dat_res$var_fac_ratio, real_dat_res$m_diff_pro,
    method = "spearman")

plot(real_dat_res$var_fac_ratio, real_dat_res$m_diff_pro)

# Maximum differences
summary(lm(max_diff_pro ~ var_fac_ratio, real_dat_res))
plot(lm(max_diff_pro ~ var_fac_ratio, real_dat_res))

cor(real_dat_res$var_fac_ratio, real_dat_res$max_diff_pro,
    method = "spearman")

plot(real_dat_res$var_fac_ratio, real_dat_res$max_diff_pro)

# Sample size-to-factor ratio

# log transform N-to-fac ratio due to skewed distribution
hist_n_fac_ratio <- hist(real_dat_res$N_fac_ratio)
real_dat_res$N_fac_ratio_log <- log10(real_dat_res$N_fac_ratio)
hist_n_fac_ratio_log <- hist(real_dat_res$N_fac_ratio_log)

# Mean differences
summary(lm(m_diff_pro ~ N_fac_ratio_log, real_dat_res))
plot(lm(m_diff_pro ~ N_fac_ratio_log, real_dat_res))

cor(real_dat_res$N_fac_ratio_log, real_dat_res$m_diff_pro,
    method = "spearman", use = "pairwise")

plot(real_dat_res$N_fac_ratio_log, real_dat_res$m_diff_pro)

# Maximum differences
summary(lm(max_diff_pro ~ N_fac_ratio_log, real_dat_res))
plot(lm(max_diff_pro ~ N_fac_ratio_log, real_dat_res))

cor(real_dat_res$N_fac_ratio_log, real_dat_res$max_diff_pro,
    method = "spearman", use = "pairwise")

plot(real_dat_res$N_fac_ratio_log, real_dat_res$max_diff_pro)

# Communalities

# Mean differences
summary(lm(m_diff_pro ~ m_h2, real_dat_res))
plot(lm(m_diff_pro ~ m_h2, real_dat_res))

cor(real_dat_res$m_h2, real_dat_res$m_diff_pro,
    method = "spearman", use = "pairwise")

plot(real_dat_res$m_h2, real_dat_res$m_diff_pro)

# Maximum differences
summary(lm(max_diff_pro ~ m_h2, real_dat_res))
plot(lm(max_diff_pro ~ m_h2, real_dat_res))

cor(real_dat_res$m_h2, real_dat_res$max_diff_pro,
    method = "spearman", use = "pairwise")

plot(real_dat_res$m_h2, real_dat_res$max_diff_pro)

# Factor intercorrelations

# Mean differences
summary(lm(m_diff_pro ~ m_fac_cor, real_dat_res))
plot(lm(m_diff_pro ~ m_fac_cor, real_dat_res))

cor(real_dat_res$m_fac_cor, real_dat_res$m_diff_pro,
    method = "spearman", use = "pairwise")

plot(real_dat_res$m_fac_cor, real_dat_res$m_diff_pro)

# Maximum differences
summary(lm(max_diff_pro ~ m_fac_cor, real_dat_res))
plot(lm(max_diff_pro ~ m_fac_cor, real_dat_res))

cor(real_dat_res$m_fac_cor, real_dat_res$max_diff_pro,
    method = "spearman", use = "pairwise")

plot(real_dat_res$m_fac_cor, real_dat_res$max_diff_pro)


# Bayesian regression analyses

# Variable-to-factor ratio

# Mean differences
mod_var_fac_m <- stan_glm(m_diff_pro + 0.0001 ~ var_fac_ratio,
                data = real_dat_res, cores = 1,
                chains = 4, iter = 15000, family = Gamma(link = "log"),
                prior_intercept = normal(location = 0, scale = 10),
                prior = normal(location = 0, scale = 2.5))

describe_posterior(mod_var_fac_m, ci = .95, rope_ci = 1)

# Maximum differences
mod_var_fac_max <- stan_glm(max_diff_pro + 0.0001 ~ var_fac_ratio,
                          data = real_dat_res, cores = 1,
                          chains = 4, iter = 15000, family = Gamma(link = "log"),
                          prior_intercept = normal(location = 0, scale = 10),
                          prior = normal(location = 0, scale = 2.5))

describe_posterior(mod_var_fac_max, ci = .95, rope_ci = 1)

# Sample size-to-factor ratio

# Mean differences
mod_N_fac_m <- stan_glm(m_diff_pro + 0.0001 ~ N_fac_ratio_log,
                        data = real_dat_res, cores = 1,
                        chains = 4, iter = 15000, family = Gamma(link = "log"),
                        prior_intercept = normal(location = 0, scale = 10),
                        prior = normal(location = 0, scale = 2.5))

describe_posterior(mod_N_fac_m, ci = .95, rope_ci = 1)

# Maximum differences
mod_N_fac_max <- stan_glm(max_diff_pro + 0.0001 ~ N_fac_ratio_log,
                          data = real_dat_res, cores = 1,
                          chains = 4, iter = 15000, family = Gamma(link = "log"),
                          prior_intercept = normal(location = 0, scale = 10),
                          prior = normal(location = 0, scale = 2.5))

describe_posterior(mod_N_fac_max, ci = .95, rope_ci = 1)

# Communalities

# Mean differences
mod_h2_m <- stan_glm(m_diff_pro + 0.0001 ~ m_h2,
                        data = real_dat_res, cores = 1,
                        chains = 4, iter = 15000, family = Gamma(link = "log"),
                        prior_intercept = normal(location = 0, scale = 10),
                        prior = normal(location = 0, scale = 2.5))

describe_posterior(mod_h2_m, ci = .95, rope_ci = 1)

# Maximum differences
mod_h2_max <- stan_glm(max_diff_pro + 0.0001 ~ m_h2,
                          data = real_dat_res, cores = 1,
                          chains = 4, iter = 15000, family = Gamma(link = "log"),
                          prior_intercept = normal(location = 0, scale = 10),
                          prior = normal(location = 0, scale = 2.5))

describe_posterior(mod_h2_max, ci = .95, rope_ci = 1)

# Factor intercorrelations

# Mean differences
mod_fac_cor_m <- stan_glm(m_diff_pro + 0.0001 ~ m_fac_cor,
                     data = real_dat_res, cores = 1,
                     chains = 4, iter = 15000, family = Gamma(link = "log"),
                     prior_intercept = normal(location = 0, scale = 10),
                     prior = normal(location = 0, scale = 2.5))

describe_posterior(mod_fac_cor_m, ci = .95, rope_ci = 1)

# Maximum differences
mod_fac_cor_max <- stan_glm(max_diff_pro + 0.0001 ~ m_fac_cor,
                       data = real_dat_res, cores = 1,
                       chains = 4, iter = 15000, family = Gamma(link = "log"),
                       prior_intercept = normal(location = 0, scale = 10),
                       prior = normal(location = 0, scale = 2.5))

describe_posterior(mod_fac_cor_max, ci = .95, rope_ci = 1)


### Create plots for variable to factor correspondences ------------

## Box plots for factor intercorrelations
phi_list <- list()

for(i in names(fin_sol)){
  if(!is.na(fin_sol[[i]]$psych) && !is.na(fin_sol[[i]]$spss) && 
     ncol(fin_sol[[i]]$psych$unrot_loadings) > 1){
  temp <- fin_sol[[i]]$psych$Phi
  temp <- as.vector(temp[lower.tri(temp)])
  phi_list[[i]] <- data.frame(data = i, Phi = temp, diff_corres = 
                                fin_sol[[i]]$diff_psych_spss_pro$diff_corres)
  }
}

phi_dat <- do.call(rbind, phi_list)

mds_phi <- phi_dat %>%  
  mutate(data = as.character(data)) %>% 
  group_by(data) %>% 
  summarise(md_phi = median(abs(Phi))) %>% 
  ungroup() %>% 
  arrange(md_phi) %>% 
  pull(data)

max_phi_ind <- phi_dat %>%  
  mutate(data = as.character(data)) %>% 
  group_by(data) %>% 
  summarise(md_phi = median(abs(Phi)),
            max_phi = max(abs(Phi))) %>% 
  ungroup() %>% 
  mutate(col_ind = case_when(max_phi < .5 ~ 0,
                              TRUE ~ 1))
  
phi_dat <- phi_dat %>% 
  left_join(max_phi_ind, by = "data")

# How many % of the datasets have at least one correlation > .50?
prop.table(table(phi_dat$col_ind))

# Plot the phis
pdf("plots/phi_plot.pdf", height = 7, width = 11)
phi_dat %>% 
  mutate(data = factor(data, levels = mds_phi)) %>% 
  ggplot(aes(x = data, y = abs(Phi))) +
  geom_hline(yintercept = 0.5, lty = 2, color = "red") +
  geom_boxplot(aes(color = as.factor(col_ind))) +
  scale_color_manual(values = c("darkgray", "black")) +
  labs(x = "Dataset", y = "Absolute Factor Intercorrelations") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
dev.off()


