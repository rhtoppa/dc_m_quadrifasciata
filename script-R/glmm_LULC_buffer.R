# =============================================================================
# TITLE: GLMM Analysis of Land Use/Land Cover Buffer Metrics on 
#        Melipona quadrifasciata Return Probability
# AUTHOR: Toppa RH, et al.
# AFFILIATION: Federal University of São Carlos, Campus Sorocaba
# CONTACT: toppa@ufscar.br
# DESCRIPTION:
# This script fits GLMM binomial models with z-score standardized buffer metrics
# to assess their effects on Melipona quadrifasciata return probability. 
# Includes odds ratio calculations, forest plots, and model diagnostics for
# landscape conservation applications.
# ==========================================

# Packages ---------------------------------------------------------------
req <- c("readxl","dplyr","tidyr","glmmTMB","MuMIn","DHARMa",
         "broom.mixed","performance","ggplot2","stringr","scales","readr")
new <- req[!(req %in% rownames(installed.packages()))]
if (length(new)) install.packages(new, dependencies = TRUE)

library(readxl); library(dplyr); library(tidyr); library(glmmTMB); library(MuMIn)
library(DHARMa); library(broom.mixed); library(performance); library(ggplot2)
library(stringr); library(scales); library(readr)

# Adjust the parameters for data input and output --------------------------------------------------------------
xlsx   <- "C:/r_studio/statistic_data.xlsx"
OUTDIR <- "C:/r_studio/buffer_glmm"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Read data and minimum typing -------------------------------------------
dat <- read_excel(xlsx)

df <- dat %>%
  mutate(
    event_use     = as.integer(event_use),
    distance_km   = suppressWarnings(as.numeric(distance_km)),  # Mantemos em km (não padronizado)
    colony        = factor(colony),
    release_point = factor(release_point),
    release_batch = as.Date(release_batch)
  )

# Detects continuous BUFFER variables and standardizes (z-score) ----------
# Important: the names remain the SAME; the values become z-scores.
is_cont_buffer <- function(nm) {
  grepl("_BUFFER$", nm) && !grepl("^period$|^colony$|^release_", nm)
}
BUFF_VARS <- names(df)[sapply(names(df), is_cont_buffer)]

if (length(BUFF_VARS)) {
  df <- df %>%
    mutate(across(all_of(BUFF_VARS), ~ as.numeric(scale(.x)), .names = "{.col}"))
}

# WRITE YOUR FORMULA (use the names from the spreadsheet) ------------------------
FORM <- event_use ~ distance_km + MAX_NDVI_BUFFER +
  (1|colony) + (1|release_batch) + (1|release_point)

cat("\nFórmula usada:\n"); print(FORM)

# Subset without NA in the used columns --------------------------------
vars_used <- unique(c("event_use", all.vars(FORM)))
missing <- setdiff(vars_used, names(df))
if (length(missing)) stop("Variáveis ausentes na planilha: ", paste(missing, collapse=", "))

dsub <- df %>% select(all_of(vars_used)) %>% drop_na()
cat("N linhas após remoção de NAs:", nrow(dsub), "\n")

# GLMM Adjustment --------------------------------------------------------
m <- glmmTMB(FORM, family = binomial(link="logit"), data = dsub)

# Summary, AICc, R2 and DHARMa -------------------------------------------
sink(file.path(OUTDIR, "model_summary_BUFFER.txt"))
print(summary(m))
cat(sprintf("\nAICc = %.2f\n", MuMIn::AICc(m)))
print(performance::r2(m))   # R2 Nakagawa
sink()

print(summary(m))
cat(sprintf("\nAICc = %.2f\n", MuMIn::AICc(m)))
print(performance::r2(m))

png(file.path(OUTDIR,"DHARMa_residuals_BUFFER.png"), width=900, height=700, res=120)
plot(DHARMa::simulateResiduals(m, n = 1000))
dev.off()

suppressWarnings({
  colin <- performance::check_collinearity(m)
  capture.output(colin, file = file.path(OUTDIR,"collinearity_BUFFER.txt"))
})

# Table of ORs (interpretable) ----------------------------------------
# distance_km: OR by +1 km
# variables BUFFER (z-score): OR by +1 standard deviation
fx <- broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE)

or_tbl <- fx %>%
  mutate(
    OR   = exp(estimate),
    LCL  = exp(conf.low),
    UCL  = exp(conf.high),
    scale_note = dplyr::case_when(
      term == "distance_km" ~ "per +1 km",
      term != "(Intercept)" & term %in% BUFF_VARS ~ "per +1 SD (z-score)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(term != "(Intercept)") %>%
  transmute(term, OR, LCL, UCL, scale_note)

write_csv(or_tbl, file.path(OUTDIR, "OR_table_BUFFER.csv"))
print(or_tbl)

# Forest plot on log(OR) scale ----------------------------------------
# Log scale avoids OR "explosions" and facilitates reading at various orders.
fx_plot <- fx %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR  = exp(estimate),
    LCL = exp(conf.low),
    UCL = exp(conf.high),
    label = dplyr::case_when(
      term == "distance_km" ~ "Distance (per 1 km)",
      TRUE ~ stringr::str_replace_all(term, "_", " ")
    )
  )

p_forest <- ggplot(fx_plot, aes(x = OR, y = reorder(label, OR))) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey40") +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = LCL, xmax = UCL, y = reorder(label, OR)),
                orientation = "y", height = .15) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    x = "Odds ratio (95% CI)  [log scale]",
    y = NULL,
    title = "Fixed effects on return probability (BUFFER)",
    subtitle = "Distance: per +1 km | BUFFER indices: per +1 SD (z-score)"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUTDIR,"forest_OR_BUFFER_logscale.png"), p_forest, width=7.5, height=5.5, dpi=300)
print(p_forest)

cat("\nArquivos salvos em:\n", OUTDIR, "\n")


