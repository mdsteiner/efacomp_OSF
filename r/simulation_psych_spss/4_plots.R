if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if(!require(patchwork)) install.packages("patchwork"); library(patchwork)
if(!require(viridis)) install.packages("viridis"); library(viridis)
if(!require(EFAtools)) install.packages("EFAtools"); library(EFAtools)
if(!require(stringr)) install.packages("stringr"); library(stringr)

recovery <- readRDS("output/simulation_psych_spss/model_recovery_results.RDS")

case_colors <- viridis(12, end = .85)
case_colors <- case_colors[c(1, rep(2, 4), 3:7, 7:9, rep(10, 4),
                             rep(11:12, each = 5))]
names(case_colors) <- gsub("case_", " ", names(population_models$loadings[1:27]))

# plot the mean absolute differences in the loadings
p1 <- recovery %>%
  filter(N == 180) %>%
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong"))) %>%
  ggplot(aes(case, m_delta_load_spss_psych, col = case)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .1, alpha = .2) +
  facet_wrap(~ cors, ncol = 1) +
  scale_color_manual(values = case_colors) +
  labs(x = "Population Model", y = expression(M[paste("|", Delta, Pattern~Coefficients, "|")]~SPSS~vs.~R~psych)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black")
  )

p2 <- recovery %>%
  filter(N == 180) %>%
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong")),
         delta_g = g_spss - g_psych) %>%
  ggplot(aes(case, delta_g, col = case)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .1, alpha = .2) +
  facet_wrap(~ cors, ncol = 1) +
  scale_color_manual(values = case_colors) +
  labs(x = "Population Model", y = expression(Delta~RMSE~SPSS~-~R~psych)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black")
  )
pdf("plots/psych_vs_spss_simulations_180.pdf", height = 6, width = 15)
p1 + p2
dev.off()

p1 <- recovery %>%
  filter(N == 450) %>%
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong"))) %>%
  ggplot(aes(case, m_delta_load_spss_psych, col = case)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .1, alpha = .2) +
  facet_wrap(~ cors, ncol = 1) +
  scale_color_manual(values = case_colors) +
  labs(x = "Population Model", y = expression(M[paste("|", Delta, Pattern~Coefficients, "|")]~SPSS~vs.~R~psych)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black")
  )

p2 <- recovery %>%
  filter(N == 450) %>%
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong")),
         delta_g =  g_spss - g_psych) %>%
  ggplot(aes(case, delta_g, col = case)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .1, alpha = .2) +
  facet_wrap(~ cors, ncol = 1) +
  scale_color_manual(values = case_colors) +
  labs(x = "Population Model", y = expression(Delta~RMSE~SPSS~-~R~psych)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black")
  )

pdf("plots/psych_vs_spss_simulations_450.pdf", height = 6, width = 15)
p1 + p2
dev.off()


# plot factor correspondencses ============================================

p1 <- recovery %>%
  filter(N == 180) %>%
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong"))) %>%
  ggplot(aes(case, diff_factor_corres_spss_psych, col = case)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .1, alpha = .2) +
  facet_wrap(~ cors, ncol = 1) +
  scale_color_manual(values = case_colors) +
  labs(x = "Population Model", y = expression(Delta~Ind-to-Fac~Corres~SPSS~vs.~R~psych)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black")
  )

p2 <- recovery %>%
  filter(N == 180) %>%
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong")),
         delta_factor_corres = diff_factor_corres_spss - diff_factor_corres_psych) %>%
  ggplot(aes(case, delta_factor_corres, col = case)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .1, alpha = .2) +
  facet_wrap(~ cors, ncol = 1) +
  scale_color_manual(values = case_colors) +
  labs(x = "Population Model", y = expression(Delta~Acc~Ind-to-Fac~Corres~SPSS~-~R~psych)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black")
  )

pdf("plots/factor_corres_spss_psych_180.pdf", height = 6, width = 15)
p1 + p2
dev.off()

p1 <- recovery %>%
  filter(N == 450) %>%
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong"))) %>%
  ggplot(aes(case, diff_factor_corres_spss_psych, col = case)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .1, alpha = .2) +
  facet_wrap(~ cors, ncol = 1) +
  scale_color_manual(values = case_colors) +
  labs(x = "Population Model", y = expression(Delta~Ind-to-Fac~Corres~SPSS~vs.~R~psych)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black")
  )

p2 <- recovery %>%
  filter(N == 450) %>%
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong")),
         delta_factor_corres = diff_factor_corres_spss - diff_factor_corres_psych) %>%
  ggplot(aes(case, delta_factor_corres, col = case)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .1, alpha = .2) +
  facet_wrap(~ cors, ncol = 1) +
  scale_color_manual(values = case_colors) +
  labs(x = "Population Model", y = expression(Delta~Acc~Ind-to-Fac~Corres~SPSS~-~R~psych)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black")
  )

pdf("plots/factor_corres_spss_psych_450.pdf", height = 6, width = 15)
p1 + p2
dev.off()

# plot heywood cases ============================================

p1 <- recovery %>%
  filter(N == 180) %>%
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong"))) %>%
  group_by(case, cors) %>% 
  summarise(
    m_heywood_psych = mean(heywood_psych),
    m_heywood_spss = mean(heywood_spss)
  ) %>% 
  ungroup() %>% 
  mutate(delta_heywood = m_heywood_spss - m_heywood_psych) %>% 
  ggplot(aes(case, delta_heywood, fill = case)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ cors, ncol = 1) +
  scale_fill_manual(values = case_colors) +
  labs(x = "Population Model", y = expression(Delta~"p(Heywood)"~SPSS~-~R~psych)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black")
  )

p2 <- recovery %>%
  filter(N == 450) %>%
  mutate(case = factor(case, levels = names(population_models$loadings)[1:27],
                       labels = gsub("case_", " ", names(population_models$loadings)[1:27])),
         cors = str_to_title(cors),
         cors = factor(cors, levels = c("Zero", "Moderate", "Mixed", "Strong"))) %>%
  group_by(case, cors) %>% 
  summarise(
    m_heywood_psych = mean(heywood_psych),
    m_heywood_spss = mean(heywood_spss)
  ) %>% 
  ungroup() %>% 
  mutate(delta_heywood = m_heywood_spss - m_heywood_psych) %>% 
  ggplot(aes(case, delta_heywood, fill = case)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ cors, ncol = 1) +
  scale_fill_manual(values = case_colors) +
  labs(x = "Population Model", y = expression(Delta~"p(Heywood)"~SPSS~-~R~psych)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1,
                               colour = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = 0, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "#f2f2f2", colour = "black"),
    axis.line = element_line(colour = "black")
  )
pdf("plots/heywood_spss_psych.pdf", height = 6, width = 15)
p1 + p2
dev.off()


