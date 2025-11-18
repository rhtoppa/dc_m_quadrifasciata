# ============================================
# TITLE: GLMM Analysis of Return Probability and Distance Threshold Detection
# AUTHOR: Toppa RH, et al.
# AFFILIATION: Federal University of São Carlos, Campus Sorocaba
# CONTACT: toppa@ufscar.br
# DESCRIPTION:
# This script analyzes Melipona quadrifasciata return probability using GLMM 
# binomial models with distance as main predictor. Includes threshold detection 
# using hinge functions and comprehensive visualization of distance effects
# on homing success in fragmented landscapes.
# ============================================

# Install & load
base_pkgs <- c("readxl","car","DHARMa","ggplot2","ggeffects","sjPlot","MuMIn",
               "dplyr","glmmTMB","binom","broom.mixed","tidyr")
new <- base_pkgs[!(base_pkgs %in% rownames(installed.packages()))]
if (length(new)) install.packages(new, dependencies = TRUE)

library(readxl); library(car); library(DHARMa); library(ggplot2)
library(ggeffects); library(sjPlot); library(MuMIn)
library(dplyr); library(glmmTMB); library(binom); library(broom.mixed); library(tidyr)

# Import data
data <- read_excel("G:/Meu Drive/UFSCar/projetos/UNIVERSAL_CNPq_2023/publicacoes/new_distance/nova_versao_artigo/submeter_discover/revisao/sheet_all_data_statistic_discovery.xlsx")

# Preparing the dataset ----------------------------------------------
required_cols <- c("event_use","distance_km","colony","release_batch","release_point")
stopifnot(all(required_cols %in% names(data)))

data <- data %>%
  mutate(
    distance_km   = as.numeric(distance_km),          # km (como pedido)
    colony        = factor(colony),
    release_point = factor(release_point),
    release_batch = as.Date(release_batch),           # garante Date
    period        = factor(period, levels = sort(unique(period)))
  )

# Remove NAs only in the model variables
df <- data %>%
  select(all_of(required_cols)) %>%
  drop_na()

# GLMM binomial (logit): distance + period + REs ----------------------
m_glmm <- glmmTMB(
  event_use ~ distance_km +
    (1 | colony) + (1 | release_batch) + (1 | release_point),
  family = binomial(link = "logit"),
  data   = df
)
print(summary(m_glmm))

# DHARMa Diagnosis -----------------------------------------------------
sim <- DHARMa::simulateResiduals(m_glmm, n = 1000)
plot(sim)  # verifique uniformidade / dispersão / outliers (mensagens do qgam são OK)

# ORs per 1 km and 0,5 km + IC95% --------------------------------------
fx <- broom.mixed::tidy(m_glmm, effects = "fixed", conf.int = TRUE)
beta_dist <- fx %>% filter(term == "distance_km")
stopifnot(nrow(beta_dist) == 1)

OR_1km   <- exp(beta_dist$estimate)
ORlo_1   <- exp(beta_dist$conf.low)
ORhi_1   <- exp(beta_dist$conf.high)
OR_0_5   <- exp(beta_dist$estimate * 0.5)
ORlo_0_5 <- exp(beta_dist$conf.low  * 0.5)
ORhi_0_5 <- exp(beta_dist$conf.high * 0.5)

cat("\n# Odds Ratios (distância)\n")
cat(sprintf("OR por 1 km:  %.3f (IC95%% %.3f–%.3f)\n", OR_1km, ORlo_1, ORhi_1))
cat(sprintf("OR por 0,5 km: %.3f (IC95%% %.3f–%.3f)\n\n", OR_0_5, ORlo_0_5, ORhi_0_5))

gg_pred <- ggpredict(
  m_glmm,
  terms = c("distance_km [0:7.5 by=0.25]")
) %>% as.data.frame()

keyd <- c(0, 2, 4, 6, 7.5)
tab_preds <- gg_pred %>%
  filter(x %in% keyd) %>%
  transmute(distance_km = x,
            pred = predicted,
            lo = conf.low,
            hi = conf.high)

cat("# Probabilidades previstas (IC95%) – período ref =", ref_period, "\n")
print(tab_preds, row.names = FALSE)

# Figure: points (binomial CI) + GLMM curve (±CI) ---------------------
# (a) observed proportion + binomial CI by distance
bydist <- df %>%
  group_by(distance_km) %>%
  summarise(n = n(), x = sum(event_use == 1), .groups = "drop")
ci <- binom::binom.confint(bydist$x, bydist$n, methods = "wilson")
bydist <- bind_cols(bydist, ci[, c("lower","upper")])

# (b) projected curve (and CI range) for the reference period
cur <- gg_pred %>% rename(y = predicted, ymin = conf.low, ymax = conf.high)

p <- ggplot(bydist, aes(x = distance_km, y = x/n)) +
  geom_point(aes(size = n), alpha = 0.85) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.08) +
  geom_ribbon(data = cur, aes(x = x, y = y, ymin = ymin, ymax = ymax), alpha = 0.15) +
  geom_line(data = cur, aes(x = x, y = y), linewidth = 1) +
  scale_size_continuous(
    name = "N per distance (sample size)",
    breaks = c(30,60,90,120),
    labels = c("30","60","90","120")) +
  labs(x = "Distance (km)",
       y = "Return probability (binomial 95% CI)",
       title = "Return vs. distance: observed (Wilson CI) and GLMM (±95% CI)",
       subtitle = "Inference from distance-only GLMM; points = descriptive data") +
  theme_minimal(base_size = 12)
print(p)

# Threshold (hinge) – DISTANCE-ONLY, with fair comparison ------------
# - Simple distance-only model (for AICc reference)
m_dist_only <- glmmTMB(
  event_use ~ distance_km +
    (1 | colony) + (1 | release_batch) + (1 | release_point),
  family = binomial(link = "logit"),
  data   = df
)

fit_with_k <- function(k){
  df$hinge_k <- pmax(0, df$distance_km - k)
  m_k <- try(glmmTMB(
    event_use ~ distance_km + hinge_k +
      (1|colony) + (1|release_batch) + (1|release_point),
    family = binomial(link="logit"), data = df
  ), silent = TRUE)
  if (inherits(m_k, "try-error")) return(Inf)
  MuMIn::AICc(m_k)
}

k_seq  <- seq(3, 6, by = 0.1)
aicc_v <- sapply(k_seq, fit_with_k)

best_k <- k_seq[which.min(aicc_v)]
delta  <- min(aicc_v) - MuMIn::AICc(m_dist_only)  # ΔAICc vs distância-only simples

cat(sprintf("\n# Threshold (hinge, distance-only) – melhor k = %.1f km; ΔAICc (hinge - simples) = %.2f\n",
            best_k, delta))
if (is.finite(delta) && delta < -2) {
  cat("→ Evidência a favor de um breakpoint (ΔAICc < -2). Reporte k e IC.\n")
} else {
  cat("→ Sem melhora clara de AICc. Interprete como declínio contínuo.\n")
}



