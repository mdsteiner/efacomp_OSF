## Analyses simulation data settings

if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if(!require(EFAtools)) install.packages("EFAtools"); library(EFAtools)
if(!require(stringr)) install.packages("stringr"); library(stringr)

names_cases <- c("18|3|6", "6|3|6", "9|3|6", "12|3|6", "15|3|6", "18|3|3",
                 "18|3|9", "18|3|369b", "18|3|369w", "18|3|46|1c", "18|3|46|3c",
                 "12|3m|6", "18|3|6n", "6|3|369wb", "9|3|369wb", "12|3|369wb",
                 "15|3|369wb", "12|6|6", "18|6|6", "24|6|6", "30|6|6", "36|6|6",
                 "12|6|369wb", "18|6|369wb", "24|6|369wb", "30|6|369wb", "36|6|369w")

model_control <- expand.grid(
  case = names(population_models$loadings)[1:27],
  cors = names(population_models$phis_3),
  N = c(180),
  stringsAsFactors = FALSE)
model_control$case_ids <- 1:nrow(model_control)


# setting ids:
# psych: 108
# psych with unity: 106
# SPSS: 171
settings <- expand.grid(
  comm_meth = c("unity", "mac", "smc"),  # initial communality estimation methods
  criterion_type = c("max_individual", "sums"), # citerion types
  abs_eigen = c(TRUE, FALSE), # absolute eigenvalues yes or no
  conv_crit = c(1e-3, 1e-6), # convergence criterion
  var_type = c("svd", "kaiser"), # varimax type
  p_type = c("unnorm", "norm"),
  k = c(3, 4),
  stringsAsFactors = FALSE) %>% 
  as_tibble()

# define k according to earlier simulation studies
settings$k[settings$k == 3 & settings$p_type == "norm"] <- 2

settings$setting_id <- 1:nrow(settings)

dev_180 <- readRDS("output/simulation_settings/results_deviance_settings_180.RDS")
dev_180_NE <- readRDS("output/simulation_settings/results_deviance_NE_settings_180.RDS")
dev_450 <- readRDS("output/simulation_settings/results_deviance_settings_450.RDS")
dev_450_NE <- readRDS("output/simulation_settings/results_deviance_NE_settings_450.RDS")

hey_180 <- readRDS("output/simulation_settings/results_heywood_settings_180.RDS")
hey_180_NE <- readRDS("output/simulation_settings/results_heywood_NE_settings_180.RDS")
hey_450 <- readRDS("output/simulation_settings/results_heywood_settings_450.RDS")
hey_450_NE <- readRDS("output/simulation_settings/results_heywood_NE_settings_450.RDS")

fc_180 <- readRDS("output/simulation_settings/results_factor_corres_settings_180.RDS")
fc_180_NE <- readRDS("output/simulation_settings/results_factor_corres_NE_settings_180.RDS")
fc_450 <- readRDS("output/simulation_settings/results_factor_corres_settings_450.RDS")
fc_450_NE <- readRDS("output/simulation_settings/results_factor_corres_NE_settings_450.RDS")

datnames <- c("dev_180", "dev_180_NE", "dev_450", "dev_450_NE", "hey_180",
              "hey_180_NE", "hey_450", "hey_450_NE", "fc_180",
              "fc_180_NE", "fc_450", "fc_450_NE")

datlist <- list(
  dev_180 = dev_180,
  dev_180_NE = dev_180_NE,
  dev_450 = dev_450,
  dev_450_NE = dev_450_NE,
  hey_180 = hey_180,
  hey_180_NE = hey_180_NE,
  hey_450 = hey_450,
  hey_450_NE = hey_450_NE,
  fc_180 = fc_180,
  fc_180_NE = fc_180_NE,
  fc_450 = fc_450,
  fc_450_NE = fc_450_NE
)

res_list <- list()
dev_ids <- c()
hey_ids <- c()
fac_ids <- c()
for (dat_i in 1:length(datlist)){
  
  pp <- matrix(0, ncol = 2, nrow = 192)
  pp[, 1] <- 1:192
  temp <- datlist[[dat_i]] %>% 
    mutate(best_ids = as.character(best_ids)) #%>% 
  #filter(evidence == "Different")
  for (i in 1:192) {
    pp[i, 2] <-  mean(sapply(temp$best_ids[temp$N_datasets > 10], 
                             function(j) i %in% 
                               as.numeric(str_split(j, ";", simplify = TRUE))))
  }
  
  colnames(pp) <- c("id", "percent")
  
  pp <- as_tibble(pp) %>% 
    arrange(desc(percent))
  ids <- pp$id[pp$percent == max(pp$percent)]
  
  res_list[[datnames[dat_i]]] <- list(
    all = pp,
    best = ids
  )
  
  if (dat_i %in% 1:4) {
    dev_ids <- c(dev_ids, ids)
  } else if (dat_i %in% 5:8) {
    hey_ids <- c(hey_ids, ids)
  } else if (dat_i %in% 9:12){
    fac_ids <- c(fac_ids, ids)
  }
  
}


table(dev_ids) %>% sort() # 5, 6, 11, 29, 30, 35, 51, 54, 75, 78,  174, 180
table(hey_ids) %>% sort() %>% {names(.)[. == 4]} %>% paste0(., collapse = ", ")  # 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 16, 22, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 40, 46, 49, 50, 51, 52, 53, 54, 55, 56, 58, 59, 64, 70, 73, 74, 75, 76, 77, 78, 79, 80, 82, 83, 88, 94, 97, 98, 99, 100, 101, 102, 103, 104, 106, 107, 112, 118, 121, 122, 123, 124, 125, 126, 127, 128, 130, 131, 136, 142, 145, 146, 147, 148, 149, 150, 151, 152, 154, 155, 160, 166, 169, 170, 171, 172, 173, 174, 175, 176, 178, 179, 184, 190
table(fac_ids) %>% sort() %>% names() %>% paste0(., collapse = ", ") # 50, 53, 56, 59, 74, 80, 99, 149, 181, 187, 6, 30, 51, 54, 75, 77, 78, 83, 102, 126, 147, 150, 171, 174

settings[settings$setting_id %in% c(5, 6, 11, 29, 30, 35, 57, 60, 81, 84, 174,
                                    180),] # dev_ids
settings[settings$setting_id %in% c(5, 6, 11, 29, 30, 35, 150, 174),] # hey_ids
settings[settings$setting_id %in% c(6, 12, 30, 33, 36, 51, 53, 54, 57, 59, 60, 75, 77, 78, 81, 83, 84, 99, 102, 105, 108, 126, 129, 132, 147, 149, 150, 153, 156, 159, 165, 171, 173, 174, 177, 179, 180, 181, 187),] # fac_ids

# psych: 108
# psych with unity: 106
# SPSS: 171

id_vec <- c(5, 6, 11, 29, 30, 35, 50, 51, 53, 54, 56, 59, 74, 75, 77, 78, 80, 83,
            99, 102, 106, 108, 126, 147, 149, 150, 171, 174, 180, 181, 187)

settings[settings$setting_id %in% id_vec,] %>% pull(comm_meth) %>% 
  toupper() %>% paste0(collapse = ",")
settings[settings$setting_id %in% id_vec,] %>%
  mutate(criterion_type = ifelse(criterion_type == "sums", "sum", "max individual")) %>% 
  pull(criterion_type) %>% paste0(collapse = ",")
settings[settings$setting_id %in% id_vec,] %>% pull(abs_eigen) %>% 
  as.character() %>% paste0(collapse = ",")
settings[settings$setting_id %in% id_vec,] %>% pull(conv_crit) %>% 
  paste0(collapse = ",")
settings[settings$setting_id %in% id_vec,] %>% pull(var_type) %>% 
  paste0(collapse = ",")
settings[settings$setting_id %in% id_vec,] %>% 
  mutate(p_type = case_when(p_type == "norm" ~ "normalized",
                            TRUE ~ "unnormalized")) %>% 
  pull(p_type) %>% paste0(collapse = ",")
settings[settings$setting_id %in% id_vec,] %>% pull(k) %>% 
  paste0(collapse = ",")


# results deviances
res_list[["dev_180"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")
res_list[["dev_180_NE"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")
res_list[["dev_450"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")
res_list[["dev_450_NE"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")

res_list[["dev_180"]][["all"]] %>%  mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))
res_list[["dev_180_NE"]][["all"]] %>% mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))
res_list[["dev_450"]][["all"]] %>%  mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))
res_list[["dev_450_NE"]][["all"]] %>% mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))


res_list[["dev_180"]][["best"]]
res_list[["dev_180_NE"]][["best"]]
res_list[["dev_450"]][["best"]]
res_list[["dev_450_NE"]][["best"]]

# results heywood cases
res_list[["hey_180"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")
res_list[["hey_180_NE"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")
res_list[["hey_450"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")
res_list[["hey_450_NE"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")

res_list[["hey_180"]][["all"]]  %>% mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))
res_list[["hey_180_NE"]][["all"]] %>% mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))
res_list[["hey_450"]][["all"]]  %>% mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))
res_list[["hey_450_NE"]][["all"]]  %>% mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))


res_list[["hey_180"]][["best"]]
res_list[["hey_180_NE"]][["best"]]
res_list[["hey_450"]][["best"]]
res_list[["hey_450_NE"]][["best"]]

# results factor correspondences
res_list[["fc_180"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")
res_list[["fc_180_NE"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")
res_list[["fc_450"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")
res_list[["fc_450_NE"]][["all"]] %>% filter(id %in% id_vec) %>%
  mutate(percent = round(percent, 2)) %>% arrange(id) %>% pull(percent) %>% 
  paste0(collapse = ",")

res_list[["fc_180"]][["all"]]  %>% mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))
res_list[["fc_180_NE"]][["all"]]  %>% mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))
res_list[["fc_450"]][["all"]]  %>% mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))
res_list[["fc_450_NE"]][["all"]]  %>% mutate(percent = round(percent, 2)) %>%
  arrange(desc(percent))


res_list[["fc_180"]][["best"]]
res_list[["fc_180_NE"]][["best"]]
res_list[["fc_450"]][["best"]]
res_list[["fc_450_NE"]][["best"]]


# setting ids:
# psych: 108
# psych with unity: 106
# SPSS: 171

# best: 174
recovery <- readRDS("output/simulation_settings/recovery_180.RDS")
# to decrease size of recovery object
recovery <- recovery %>%
  select(case:g, heywood, diff_factor_corres:error, -iter) %>% 
  mutate(
    dat_id2 = rep(1:(1000 * nrow(model_control)), each = nrow(settings))
  ) %>% 
  mutate_if(is.factor, as.character) %>%
  as_tibble() %>%
  left_join(settings, by = c("comm_meth", "criterion_type", "abs_eigen",
                             "conv_crit", "var_type", "p_type", "k")) %>% 
  left_join(model_control, by = c("case", "cors", "N")) %>% 
  filter(setting_id %in% c(174, 108, 106, 171))

unity_ids <- recovery %>% 
  filter(setting_id == 108 & !is.na(error)) %>% 
  pull(dat_id2)

recovery <- recovery %>% 
  filter(!(setting_id == 106 & !(dat_id2 %in% unity_ids)),
         is.na(error)) %>% 
  mutate(setting_id = case_when(setting_id %in% c(108, 106) ~ "R psych",
                                setting_id == 171 ~"SPSS",
                                setting_id == 174 ~ "Best"))


# 450

recovery_450 <- readRDS("output/simulation_settings/recovery_450.RDS")
# to decrease size of recovery object
recovery_450 <- recovery_450 %>%
  select(case:g, heywood, diff_factor_corres:error, -iter) %>% 
  mutate(
    dat_id2 = rep(1:(1000 * nrow(model_control)), each = nrow(settings))
  ) %>% 
  mutate_if(is.factor, as.character) %>%
  as_tibble() %>%
  left_join(settings, by = c("comm_meth", "criterion_type", "abs_eigen",
                             "conv_crit", "var_type", "p_type", "k")) %>% 
  left_join(model_control, by = c("case", "cors")) %>% 
  filter(setting_id %in% c(174, 108, 106, 171))

unity_ids <- recovery %>% 
  filter(setting_id == 108 & !is.na(error)) %>% 
  pull(dat_id2)

recovery_450 <- recovery_450 %>% 
  filter(!(setting_id == 106 & !(dat_id2 %in% unity_ids)),
         is.na(error)) %>% 
  mutate(setting_id = case_when(setting_id %in% c(108, 106) ~ "R psych",
                                setting_id == 171 ~"SPSS",
                                setting_id == 174 ~ "Best"))


temp <- recovery %>% 
  mutate(case = gsub("case_", "", case)) %>% 
  group_by(setting_id, case_ids, case, cors) %>% 
  summarise(
    m_h = mean(heywood, na.rm = TRUE),
    m_g = mean(g, na.rm = TRUE),
    sd_g = sd(g, na.rm = TRUE),
    max_g = mean(g[g >= quantile(g, prob = .95, na.rm = TRUE)], na.rm = TRUE),
    m_fc = mean(diff_factor_corres, na.rm = TRUE),
    sd_fc = sd(diff_factor_corres, na.rm = TRUE),
    max_fc = mean(diff_factor_corres[diff_factor_corres >= quantile(diff_factor_corres, prob = .95, na.rm = TRUE)], na.rm = TRUE),
    p_fc = mean(diff_factor_corres > 0, na.rm = TRUE)
  )

temp_450 <- recovery_450 %>% 
  mutate(case = gsub("case_", "", case)) %>% 
  group_by(setting_id, case_ids, case, cors) %>% 
  summarise(
    m_h = mean(heywood, na.rm = TRUE),
    m_g = mean(g, na.rm = TRUE),
    sd_g = sd(g, na.rm = TRUE),
    max_g = mean(g[g >= quantile(g, prob = .95, na.rm = TRUE)], na.rm = TRUE),
    m_fc = mean(diff_factor_corres, na.rm = TRUE),
    sd_fc = sd(diff_factor_corres, na.rm = TRUE),
    max_fc = mean(diff_factor_corres[diff_factor_corres >= quantile(diff_factor_corres, prob = .95, na.rm = TRUE)], na.rm = TRUE),
    p_fc = mean(diff_factor_corres > 0, na.rm = TRUE)
  )

### MRMSE plots
pdf("plots/MRMSE_implementations_180.pdf", height = 7, width = 12)
temp %>% 
  ungroup() %>% 
  mutate(case = factor(case, levels = gsub("case_", "", names(population_models$loadings)[1:27]),
                       labels = names_cases),
         cors = factor(cors, levels = c("zero", "moderate", "mixed", "strong"),
                       labels = c("Zero", "Moderate", "Mixed", "Strong"))) %>% 
  ggplot(aes(case, m_g, fill = factor(setting_id))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = m_g - sd_g, ymax = m_g + sd_g), width = .2,
                position=position_dodge(.9)) +
  geom_point(aes(y = max_g, col = factor(setting_id)),
             position = position_dodge(.9), show.legend = FALSE,
             size = .75) +
  labs(y = expression(bold(M[RMSE])), x = "Population Model") +
  theme_light() +
  scale_fill_manual(values = viridis::viridis(3)) +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_grid(cors~.) +
  guides(fill = guide_legend(title = "Implementation")) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
dev.off()

pdf("plots/MRMSE_implementations_450.pdf", height = 7, width = 12)
temp_450 %>% 
  ungroup() %>% 
  mutate(case = factor(case, levels = gsub("case_", "", names(population_models$loadings)[1:27]),
                       labels = names_cases),
         cors = factor(cors, levels = c("zero", "moderate", "mixed", "strong"),
                       labels = c("Zero", "Moderate", "Mixed", "Strong"))) %>% 
  ggplot(aes(case, m_g, fill = factor(setting_id))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = m_g - sd_g, ymax = m_g + sd_g), width = .2,
                position=position_dodge(.9)) +
  geom_point(aes(y = max_g, col = factor(setting_id)),
             position = position_dodge(.9), show.legend = FALSE,
             size = .75) +
  labs(y = expression(bold(M[RMSE])), x = "Population Model") +
  theme_light() +
  scale_fill_manual(values = viridis::viridis(3)) +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_grid(cors~.) +
  guides(fill = guide_legend(title = "Implementation")) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
dev.off()

### p(Heywood)
pdf("plots/heywood_implementations_180.pdf", height = 7, width = 12)

temp %>% 
  ungroup() %>% 
  mutate(case = factor(case, levels = gsub("case_", "", names(population_models$loadings)[1:27]),
                       labels = names_cases),
         cors = factor(cors, levels = c("zero", "moderate", "mixed", "strong"),
                       labels = c("Zero", "Moderate", "Mixed", "Strong"))) %>% 
  ggplot(aes(case, m_h, fill = factor(setting_id))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "p(Heywood Case)", x = "Population Model") +
  theme_light() +
  scale_fill_manual(values = viridis::viridis(3)) +
  facet_grid(cors~.) +
  guides(fill = guide_legend(title = "Implementation")) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
dev.off()

pdf("plots/heywood_implementations_450.pdf", height = 7, width = 12)
temp_450 %>% 
  ungroup() %>% 
  mutate(case = factor(case, levels = gsub("case_", "", names(population_models$loadings)[1:27]),
                       labels = names_cases),
         cors = factor(cors, levels = c("zero", "moderate", "mixed", "strong"),
                       labels = c("Zero", "Moderate", "Mixed", "Strong"))) %>% 
  ggplot(aes(case, m_h, fill = factor(setting_id))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "p(Heywood Case)", x = "Population Model") +
  theme_light() +
  scale_fill_manual(values = viridis::viridis(3)) +
  facet_grid(cors~.) +
  guides(fill = guide_legend(title = "Implementation")) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
dev.off()


### var to fa corres
pdf("plots/VFC_implementations_180.pdf", height = 7, width = 12)
temp %>% 
  ungroup() %>% 
  mutate(case = factor(case, levels = gsub("case_", "", names(population_models$loadings)[1:27]),
                       labels = names_cases),
         cors = factor(cors, levels = c("zero", "moderate", "mixed", "strong"),
                       labels = c("Zero", "Moderate", "Mixed", "Strong"))) %>% 
  ggplot(aes(case, m_fc, fill = factor(setting_id))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = m_fc - sd_fc, ymax = m_fc + sd_fc), width = .2,
                position=position_dodge(.9)) +
  geom_point(aes(y = max_fc, col = factor(setting_id)),
             position = position_dodge(.9), show.legend = FALSE,
             size = .75) +
  labs(y = expression(bold(M[Delta~Incorrect~Indicator-to-Factor~Correspondences])), x = "Population Model") +
  theme_light() +
  scale_fill_manual(values = viridis::viridis(3)) +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_grid(cors~.) +
  guides(fill = guide_legend(title = "Implementation")) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
dev.off()

pdf("plots/VFC_implementations_450.pdf", height = 7, width = 12)
temp_450 %>% 
  ungroup() %>% 
  mutate(case = factor(case, levels = gsub("case_", "", names(population_models$loadings)[1:27]),
                       labels = names_cases),
         cors = factor(cors, levels = c("zero", "moderate", "mixed", "strong"),
                       labels = c("Zero", "Moderate", "Mixed", "Strong"))) %>% 
  ggplot(aes(case, m_fc, fill = factor(setting_id))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = m_fc - sd_fc, ymax = m_fc + sd_fc), width = .2,
                position=position_dodge(.9)) +
  geom_point(aes(y = max_fc, col = factor(setting_id)),
             position = position_dodge(.9), show.legend = FALSE,
             size = .75) +
  labs(y = expression(bold(M[Delta~Incorrect~Indicator-to-Factor~Correspondences])), x = "Population Model") +
  theme_light() +
  scale_fill_manual(values = viridis::viridis(3)) +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_grid(cors~.) +
  guides(fill = guide_legend(title = "Implementation")) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
dev.off()

### p(Var to Fac Correspondences)
pdf("plots/VFC_p_implementations_180.pdf", height = 7, width = 12)

temp %>% 
  ungroup() %>% 
  mutate(case = factor(case, levels = gsub("case_", "", names(population_models$loadings)[1:27]),
                       labels = names_cases),
         cors = factor(cors, levels = c("zero", "moderate", "mixed", "strong"),
                       labels = c("Zero", "Moderate", "Mixed", "Strong"))) %>% 
  ggplot(aes(case, p_fc, fill = factor(setting_id))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "p(Inc. Indicator-to-Factor Correspondences)", x = "Population Model") +
  theme_light() +
  scale_fill_manual(values = viridis::viridis(3)) +
  facet_grid(cors~.) +
  guides(fill = guide_legend(title = "Implementation")) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
dev.off()

pdf("plots/VFC_p_implementations_450.pdf", height = 7, width = 12)
temp_450 %>% 
  ungroup() %>% 
  mutate(case = factor(case, levels = gsub("case_", "", names(population_models$loadings)[1:27]),
                       labels = names_cases),
         cors = factor(cors, levels = c("zero", "moderate", "mixed", "strong"),
                       labels = c("Zero", "Moderate", "Mixed", "Strong"))) %>% 
  ggplot(aes(case, p_fc, fill = factor(setting_id))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "p(Inc. Indicator-to-Factor Correspondences)", x = "Population Model") +
  theme_light() +
  scale_fill_manual(values = viridis::viridis(3)) +
  facet_grid(cors~.) +
  guides(fill = guide_legend(title = "Implementation")) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
dev.off()


### in which structures does it matter and in which does it not?
# tables with evidence

# no NE
table(dev_180$evidence)
table(dev_450$evidence)
table(hey_180$evidence)
table(hey_450$evidence)
table(fc_180$evidence)
table(fc_450$evidence)


# NE
table(dev_180_NE$evidence)
table(dev_450_NE$evidence)
table(hey_180_NE$evidence)
table(hey_450_NE$evidence)
table(fc_180_NE$evidence)
table(fc_450_NE$evidence)

