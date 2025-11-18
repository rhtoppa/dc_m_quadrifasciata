# =============================================================================
# TITLE: GLMM Model Selection for Corridor Landscape Metrics Using Dredge Algorithm
# AUTHOR: Toppa RH, et al.
# AFFILIATION: Federal University of São Carlos, Campus Sorocaba
# CONTACT: toppa@ufscar.br
# DESCRIPTION:
# This script performs automated model selection using the dredge algorithm
# on GLMM binomial models with corridor landscape metrics. Distance is forced
# in all models, with random effects for hierarchical data structure in
# Melipona quadrifasciata homing analysis.
# ==========================================

# Packages ---------------------------------------------------------------
req <- c(
  "readxl","dplyr","tidyr","glmmTMB","MuMIn","DHARMa",
  "broom.mixed","performance","ggplot2","stringr","ggeffects","scales"
)
new <- req[!(req %in% rownames(installed.packages()))]
if (length(new)) install.packages(new, dependencies = TRUE)

library(readxl); library(dplyr); library(tidyr); library(glmmTMB); library(MuMIn)
library(DHARMa); library(broom.mixed); library(performance); library(ggplot2)
library(stringr); library(ggeffects); library(scales)

# Adjust the parameters for data input and output --------------------------------------------------------------
xlsx   <- "G:/Meu Drive/UFSCar/projetos/UNIVERSAL_CNPq_2023/publicacoes/new_distance/nova_versao_artigo/submeter_discover/revisao/sheet_all_data_statistic_discovery.xlsx"
OUTDIR <- "G:/Meu Drive/UFSCar/projetos/UNIVERSAL_CNPq_2023/publicacoes/new_distance/nova_versao_artigo/submeter_discover/revisao/outputs"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Read data and minimum typing -------------------------------------------
dat <- read_excel(xlsx)

df <- dat %>%
  mutate(
    event_use     = as.integer(event_use),
    distance_km   = suppressWarnings(as.numeric(distance_km)),
    colony        = factor(colony),
    release_point = factor(release_point),
    release_batch = as.Date(release_batch)
  )

# Define the CORR candidates (edit here if you want to reduce/enlarge) ------
CAND_CORR <- c("MEAN_NDVI_CORR","MEAN_NDWI_CORR","MEAN_MR_CORR","MAX_NDVI_CORR","MAX_NDWI_CORR","MAX_MR_CORR", "RANGE_NDVI_CORR","RANGE_NDWI_CORR","RANGE_MR_CORR")

CAND_CORR <- CAND_CORR[CAND_CORR %in% names(df)]           # mantém só os que existem
if (length(CAND_CORR) == 0) stop("Nenhuma variável *_CORR encontrada nas colunas.")

# Set up the "full" formula (distance + all candidates) -----------------
rhs_fixed <- paste(c("distance_km", CAND_CORR), collapse = " + ")
rhs_rand  <- "(1|colony) + (1|release_batch) + (1|release_point)"
FORM_FULL <- as.formula(paste0("event_use ~ ", rhs_fixed, " + ", rhs_rand))

cat("\nFórmula FULL (para dredge):\n"); print(FORM_FULL)

# Subset without NA in the used columns
vars_used <- unique(c("event_use","distance_km", CAND_CORR, "colony","release_batch","release_point"))
dsub <- df %>% select(all_of(vars_used)) %>% tidyr::drop_na()
cat("N linhas após remoção de NAs:", nrow(dsub), "\n")

# Adjust FULL and wheel dredge (forced distance) -------------------------
m_full <- glmmTMB(FORM_FULL, family = binomial(link="logit"), data = dsub)

options(na.action = "na.fail")  # exigido pelo dredge
dd <- dredge(m_full, fixed = ~ distance_km, rank = "AICc")
# salva tabela de modelos
write.csv(as.data.frame(dd), file.path(OUTDIR,"AICc_table_CORR_dredge.csv"), row.names = FALSE)

cat("\nTop-5 modelos por AICc:\n")
print(head(dd, 5))

# Select the top-1 model and also mark the models with ΔAICc ≤ 2
best      <- get.models(dd, 1)[[1]]
delta2_ix <- which(dd$delta <= 2)
write.csv(as.data.frame(dd[delta2_ix,]),
          file.path(OUTDIR,"AICc_models_delta2_CORR2.csv"), row.names = FALSE)

# Main outputs of the best model------------------------------------
sink(file.path(OUTDIR,"best_model_summary_CORR.txt"))
cat("BEST MODEL (CORRIDOR)\n\n")
print(summary(best))
cat(sprintf("\nAICc = %.2f\n", MuMIn::AICc(best)))
print(performance::r2(best))   # R2 marginal/condicional
sink()

cat("\nResumo do melhor modelo gravado em: best_model_summary_CORR2.txt\n")
cat(sprintf("AICc (best) = %.2f\n", MuMIn::AICc(best)))
print(performance::r2(best))

# 7) Odds Ratios (distance and indices present in the model) -----------------
fx <- broom.mixed::tidy(best, effects="fixed", conf.int=TRUE)

# distance
if ("distance_km" %in% fx$term) {
  rr <- fx[fx$term=="distance_km", ]
  cat(sprintf("\nDistance OR per 1 km = %.3f (%.3f–%.3f)\n",
              exp(rr$estimate), exp(rr$conf.low), exp(rr$conf.high)))
  cat(sprintf("Distance OR per 0.5 km = %.3f (%.3f–%.3f)\n",
              exp(0.5*rr$estimate), exp(0.5*rr$conf.low), exp(0.5*rr$conf.high)))
}

# LULC indices (only those that entered the best model)
idx_terms <- fx$term[grepl("NDVI|NDWI|NDBI|MR", fx$term, ignore.case = TRUE)]
if (length(idx_terms)) {
  cat("\nIndex effects (OR per +0.1):\n")
  or_table <- lapply(idx_terms, function(t) {
    r <- fx[fx$term==t,]
    c(term=t,
      OR = exp(r$estimate*0.1),
      LCL=exp(r$conf.low*0.1),
      UCL=exp(r$conf.high*0.1))
  }) %>% dplyr::bind_rows()
  print(or_table)
  write.csv(or_table, file.path(OUTDIR,"OR_index_plus0.1_best_CORR.csv"), row.names = FALSE)
}

# 8) Forest plot of the fixed elements -------------------------------------------------
fx_plot <- fx %>%
  mutate(
    OR  = exp(estimate),
    LCL = exp(conf.low),
    UCL = exp(conf.high),
    label = dplyr::case_when(
      term == "distance_km" ~ "Distance (per 1 km)",
      TRUE ~ stringr::str_replace_all(term, "_", " ")
    )
  ) %>% filter(term != "(Intercept)")

p_forest <- ggplot(fx_plot, aes(x = reorder(label, OR), y = OR)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = .15) +
  coord_flip() +
  labs(x = NULL, y = "Odds ratio (95% CI)",
       title = "Fixed effects on return probability (CORRIDOR)") +
  theme_minimal(base_size = 12)

ggsave(file.path(OUTDIR,"forest_OR_best_CORR.png"), p_forest, width=7, height=5, dpi=300)

# Predictions by distance (0–7.5) under P25/P75 scenarios ------------------
# Creates P25/P75 scenarios for EACH index present in the best model, one at a time.
pred_list <- list()
vars_in_best <- intersect(idx_terms, names(dsub))

for (v in vars_in_best) {
  p25 <- quantile(dsub[[v]], 0.25, na.rm=TRUE)
  p75 <- quantile(dsub[[v]], 0.75, na.rm=TRUE)
  
  g_lo <- ggpredict(best, terms = c(
    "distance_km [0:7.5 by=0.25]",
    paste0(v," [", round(p25,3), "]")
  )) %>% as.data.frame() %>% mutate(scenario = paste(v, "low (P25)"))
  
  g_hi <- ggpredict(best, terms = c(
    "distance_km [0:7.5 by=0.25]",
    paste0(v," [", round(p75,3), "]")
  )) %>% as.data.frame() %>% mutate(scenario = paste(v, "high (P75)"))
  
  pred_list[[v]] <- bind_rows(g_lo, g_hi)
}

if (length(pred_list)) {
  pred_df <- bind_rows(pred_list, .id = "index")
  keyd <- c(0,2,4,6,7.5)
  tab_preds <- pred_df %>%
    filter(x %in% keyd) %>%
    transmute(index,
              scenario,
              distance_km = x,
              prob = predicted,
              lo = conf.low,
              hi = conf.high) %>%
    arrange(index, scenario, distance_km)
  
  readr::write_csv(tab_preds, file.path(OUTDIR, "predictions_best_CORR_P25P75.csv"))
  print(tab_preds, n = Inf)
}

# DHARMa Diagnosis ---------------------------------------------------
png(file.path(OUTDIR,"DHARMa_best_CORR.png"), width=900, height=700, res=120)
plot(DHARMa::simulateResiduals(best, n = 1000))
dev.off()

# Collinearity -----------------------------
suppressWarnings({
  colin <- performance::check_collinearity(best)
  capture.output(colin, file = file.path(OUTDIR,"collinearity_best_CORR.txt"))
})

cat("\nArquivos gerados em:\n", OUTDIR, "\n")


