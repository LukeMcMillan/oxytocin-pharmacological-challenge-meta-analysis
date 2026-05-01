library(readxl)
library(janitor)
library(dplyr)
library(stringr)
library(metafor)
library(clubSandwich)
library(ggplot2)

# =========================================================
# 1) READ DATA
# =========================================================
file_path <- "data/OT_meta_analysis_clean_v5.xlsx"

raw <- read_excel(file_path, sheet = 1) |>
  clean_names()

# Expected columns:
# study, group, n, mean_pre, sd_pre, mean_post, sd_post,
# sample_id, drug_class, percent_female, extraction_status, assay_type

# =========================================================
# 2) CLEAN + AUDITABLE SOURCE DATA
# =========================================================
dat <- raw |>
  mutate(
    study_label       = str_squish(as.character(study)),
    group             = str_to_lower(str_squish(as.character(group))),
    sample_id         = str_squish(as.character(sample_id)),
    drug_class        = str_to_lower(str_squish(as.character(drug_class))),
    percent_female    = suppressWarnings(as.numeric(percent_female)),
    extraction_status = factor(
      str_to_lower(str_squish(as.character(extraction_status))),
      levels = c("not_extracted", "extracted", "not_reported")
    ),
    assay_type = factor(
      str_to_lower(str_squish(as.character(assay_type))),
      levels = c("ria", "elisa", "not_reported")
    )
  ) |>
  filter(!is.na(study_label), study_label != "")

# Explicit hypertonic saline rule:
# hypertonic saline is ACTIVE osmotic challenge, not placebo
dat <- dat |>
  mutate(
    hypertonic_flag = str_detect(group, "\\bhypertonic\\b"),
    placebo_flag = if_else(
      !hypertonic_flag & (
        drug_class %in% c("placebo", "control") |
          str_detect(group, "\\b(placebo|control|vehicle|sham)\\b") |
          str_detect(group, "\\b(isotonic saline|normal saline|0\\.9% saline)\\b")
      ),
      1L, 0L
    ),
    drug_class = if_else(
      hypertonic_flag & drug_class == "placebo",
      "osmotic",
      drug_class
    )
  )

# Active dataset
dat_active <- dat |>
  filter(placebo_flag == 0) |>
  mutate(
    report_id   = str_replace(study_label, "\\s+study\\s*\\d+$", "") |> str_squish(),
    sample_id_u = paste(study_label, sample_id, sep = "__")
  )

# Audit checks
cat("Rows after removing blank studies:", nrow(dat), "\n")
cat("Placebo/control excluded:", sum(dat$placebo_flag == 1), "\n")
cat("Active k:", nrow(dat_active), "\n")
cat("Distinct study labels:", n_distinct(dat_active$study_label), "\n")
cat("Distinct reports:", n_distinct(dat_active$report_id), "\n")

# Missing-data audit
audit_missing <- dat_active |>
  summarise(
    miss_n            = sum(is.na(n)),
    miss_mean_pre     = sum(is.na(mean_pre)),
    miss_sd_pre       = sum(is.na(sd_pre)),
    miss_mean_post    = sum(is.na(mean_post)),
    miss_sd_post      = sum(is.na(sd_post)),
    miss_extraction   = sum(is.na(extraction_status)),
    miss_assay        = sum(is.na(assay_type))
  )

print(audit_missing)

# Drug-class counts
class_counts <- dat_active |>
  group_by(drug_class) |>
  summarise(
    k_effects = n(),
    n_studies = n_distinct(study_label),
    n_samples = n_distinct(sample_id_u),
    .groups = "drop"
  ) |>
  arrange(desc(k_effects))

print(class_counts)

cat("\nExtraction status frequencies:\n")
print(table(dat_active$extraction_status, useNA = "ifany"))

cat("\nAssay type frequencies:\n")
print(table(dat_active$assay_type, useNA = "ifany"))

# Verify extraction and assay labels look correct
cat("\nRaw extraction values in Excel:\n")
print(table(str_to_lower(str_squish(as.character(raw$extraction_status))), useNA = "ifany"))
cat("\nRaw assay_type values in Excel:\n")
print(table(str_to_lower(str_squish(as.character(raw$assay_type))), useNA = "ifany"))

# =========================================================
# 3) EFFECT SIZES
# =========================================================
make_es <- function(dat_in, r_assumed) {
  escalc(
    measure = "SMCC",
    m1i = mean_post, m2i = mean_pre,
    sd1i = sd_post,  sd2i = sd_pre,
    ni  = n,
    ri  = r_assumed,
    data = dat_in
  ) |>
    mutate(
      effect_id = row_number(),
      r_assumed = r_assumed
    )
}

es_r50 <- make_es(dat_active, 0.50)

# =========================================================
# 4) PRIMARY MODEL
# =========================================================
m_primary <- rma.mv(
  yi, vi,
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data = es_r50
)

V_primary    <- vcovCR(m_primary, cluster = es_r50$study_label, type = "CR2")
ct_primary   <- coef_test(m_primary, vcov = V_primary, test = "Satterthwaite")
ci_primary   <- conf_int(m_primary, vcov = V_primary)
pred_primary <- predict(m_primary)

print(m_primary)
print(ct_primary)
print(ci_primary)
print(pred_primary)

# --- I² for primary multilevel model ---
W         <- diag(1 / es_r50$vi)
X         <- model.matrix(m_primary)
P         <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
typical_v <- (m_primary$k - m_primary$p) / sum(diag(P))

sigma2_study  <- m_primary$sigma2[1]
sigma2_sample <- m_primary$sigma2[2]
sigma2_effect <- m_primary$sigma2[3]
total_var     <- sigma2_study + sigma2_sample + sigma2_effect + typical_v

cat("\n=== I² FOR PRIMARY MULTILEVEL MODEL ===\n")
cat("I2 total:        ", round((sigma2_study + sigma2_sample + sigma2_effect) / total_var * 100, 1), "%\n")
cat("I2 study-level:  ", round(sigma2_study  / total_var * 100, 1), "%\n")
cat("I2 sample-level: ", round(sigma2_sample / total_var * 100, 1), "%\n")
cat("I2 effect-level: ", round(sigma2_effect / total_var * 100, 1), "%\n")

# =========================================================
# 5) SENSITIVITY TO r
# =========================================================
run_sensitivity <- function(r_assumed) {
  es <- make_es(dat_active, r_assumed)

  m <- rma.mv(
    yi, vi,
    random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
    method = "REML",
    data = es
  )

  V  <- vcovCR(m, cluster = es$study_label, type = "CR2")
  ct <- coef_test(m, vcov = V, test = "Satterthwaite")
  ci <- conf_int(m, vcov = V)
  pr <- predict(m)

  data.frame(
    r     = r_assumed,
    est   = ct$beta[1],
    se    = ct$SE[1],
    df    = ct$df_Satt[1],
    p     = ct$p_Satt[1],
    ci_lb = ci$CI_L[1],
    ci_ub = ci$CI_U[1],
    pi_lb = pr$pi.lb[1],
    pi_ub = pr$pi.ub[1]
  )
}

sens_table <- bind_rows(
  run_sensitivity(0.30),
  run_sensitivity(0.50),
  run_sensitivity(0.70)
)

print(sens_table)

# =========================================================
# 5b) FISHER (1987) SENSITIVITY ANALYSIS
# =========================================================
es_no_fisher <- es_r50[es_r50$study_label != "Fisher 1987", ]

m_no_fisher <- rma.mv(
  yi, vi,
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data   = es_no_fisher
)

V_nf  <- vcovCR(m_no_fisher, cluster = es_no_fisher$study_label, type = "CR2")
ct_nf <- coef_test(m_no_fisher, vcov = V_nf, test = "Satterthwaite")
ci_nf <- conf_int(m_no_fisher,  vcov = V_nf)
pr_nf <- predict(m_no_fisher)

cat("\n=== FISHER (1987) SENSITIVITY ANALYSIS ===\n")
cat("Pooled SMCC:", round(ct_nf$beta[1], 3), "\n")
cat("95% CI:     ", round(ci_nf$CI_L[1], 3), "to", round(ci_nf$CI_U[1], 3), "\n")
cat("PI:         ", round(pr_nf$pi.lb, 3), "to", round(pr_nf$pi.ub, 3), "\n")

# =========================================================
# 6) DRUG CLASS MODEL
# =========================================================
m_class <- rma.mv(
  yi, vi,
  mods = ~ factor(drug_class) - 1,
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data = es_r50
)

V_class  <- vcovCR(m_class, cluster = es_r50$study_label, type = "CR2")
ct_class <- coef_test(m_class, vcov = V_class, test = "Satterthwaite")
ci_class <- conf_int(m_class, vcov = V_class)

class_results <- data.frame(
  drug_class = gsub("^factor\\(drug_class\\)", "", rownames(ct_class)),
  est   = ct_class$beta,
  se    = ct_class$SE,
  df    = ct_class$df_Satt,
  p     = ct_class$p_Satt,
  ci_lb = ci_class$CI_L,
  ci_ub = ci_class$CI_U
)

print(class_results)

# Osmotic reference model (dynamic: works for any number of drug classes)
n_drug_classes <- n_distinct(es_r50$drug_class)

m_ref <- rma.mv(
  yi, vi,
  mods = ~ relevel(factor(drug_class), ref = "osmotic"),
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data = es_r50
)

V_ref  <- vcovCR(m_ref, cluster = es_r50$study_label, type = "CR2")
ct_ref <- coef_test(m_ref, vcov = V_ref, test = "Satterthwaite")
ci_ref <- conf_int(m_ref, vcov = V_ref)

omnibus_ref <- Wald_test(
  m_ref,
  constraints = constrain_zero(2:n_drug_classes),
  vcov = V_ref,
  test = "HTZ"
)

print(ct_ref)
print(ci_ref)
print(omnibus_ref)

# =========================================================
# 6b) BH-CORRECTED PAIRWISE DRUG CLASS CONTRASTS
# =========================================================
coef_names <- names(coef(m_class))
levels_dc  <- gsub("^factor\\(drug_class\\)", "", coef_names)
n_lev      <- length(levels_dc)
beta_class <- coef(m_class)

cat("\nDrug class levels in model order:\n")
print(levels_dc)

pairs_idx      <- combn(seq_len(n_lev), 2)
contrasts_list <- vector("list", ncol(pairs_idx))

for (j in seq_len(ncol(pairs_idx))) {
  i1 <- pairs_idx[1, j]
  i2 <- pairs_idx[2, j]

  con     <- rep(0, n_lev)
  con[i1] <-  1
  con[i2] <- -1

  delta  <- sum(con * beta_class)
  se_val <- sqrt(as.numeric(t(con) %*% V_class %*% con))

  wt <- Wald_test(
    m_class,
    constraints = matrix(con, nrow = 1),
    vcov        = V_class,
    test        = "HTZ"
  )

  contrasts_list[[j]] <- data.frame(
    contrast = paste(levels_dc[i1], "vs", levels_dc[i2]),
    lev1     = levels_dc[i1],
    lev2     = levels_dc[i2],
    delta    = delta,
    se       = se_val,
    t        = delta / se_val,
    df       = wt$df_denom,
    p_raw    = wt$p_val,
    stringsAsFactors = FALSE
  )
}

pw      <- bind_rows(contrasts_list)
pw$p_BH <- p.adjust(pw$p_raw, method = "BH")
pw$sig  <- ifelse(pw$p_BH < .05, "*", "")
pw      <- pw[order(pw$p_raw), ]

pw_print <- pw
pw_print[, c("delta","se","t","df","p_raw","p_BH")] <-
  round(pw_print[, c("delta","se","t","df","p_raw","p_BH")], 3)

cat("\n=== ALL PAIRWISE DRUG CLASS CONTRASTS (BH-corrected) ===\n")
print(pw_print[, c("contrast","delta","se","t","df","p_raw","p_BH","sig")],
      row.names = FALSE)

# Key a priori: monoaminergic vs osmotic
mono_osm <- pw[
  (pw$lev1 == "monoaminergic" & pw$lev2 == "osmotic") |
  (pw$lev1 == "osmotic"       & pw$lev2 == "monoaminergic"), ]

cat("\n=== KEY A PRIORI: MONOAMINERGIC vs OSMOTIC ===\n")
print(round(mono_osm[, c("delta","se","t","df","p_raw","p_BH")], 3), row.names = FALSE)
if (mono_osm$p_BH < .05) {
  cat(">> SURVIVES BH correction (p_BH =", round(mono_osm$p_BH, 3), ")\n")
} else {
  cat(">> Does NOT survive BH correction (p_BH =", round(mono_osm$p_BH, 3), ")\n")
}

# Osmotic contrasts for Table 2
osm_rows <- pw[pw$lev1 == "osmotic" | pw$lev2 == "osmotic", ]
osm_rows <- osm_rows[order(osm_rows$p_raw), ]
cat("\n=== CONTRASTS AGAINST OSMOTIC (Table 2) ===\n")
print(osm_rows[, c("contrast","delta","se","t","df","p_raw","p_BH")],
      digits = 3, row.names = FALSE)

# =========================================================
# 7) STUDY-COLLAPSED BIAS CHECK
# =========================================================
collapsed <- es_r50 |>
  group_by(study_label) |>
  summarise(
    yi = weighted.mean(yi, 1 / vi),
    vi = 1 / sum(1 / vi),
    .groups = "drop"
  )

m_coll <- rma(yi, vi, data = collapsed, method = "REML")
egger  <- regtest(m_coll, model = "lm")

print(m_coll)
print(egger)

# =========================================================
# 8) SEX MODERATOR
# =========================================================
es_sex <- es_r50 |>
  mutate(
    pf_num = suppressWarnings(as.numeric(percent_female)),
    prop_female = case_when(
      pf_num > 1  ~ pf_num / 100,
      pf_num >= 0 ~ pf_num,
      TRUE        ~ NA_real_
    )
  ) |>
  filter(!is.na(prop_female))

cat("Sex moderator dataset:", nrow(es_sex), "effects from",
    n_distinct(es_sex$study_label), "studies\n")

m_sex <- rma.mv(
  yi, vi,
  mods = ~ prop_female,
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data = es_sex
)

V_sex  <- vcovCR(m_sex, cluster = es_sex$study_label, type = "CR2")
ct_sex <- coef_test(m_sex, vcov = V_sex, test = "Satterthwaite")
ci_sex <- conf_int(m_sex, vcov = V_sex)

print(ct_sex)
print(ci_sex)

# =========================================================
# 9) METHODOLOGICAL MODERATORS: EXTRACTION + ASSAY TYPE
# =========================================================

# --- Confounding check cross-tabs ---
cat("\n=== CROSS-TAB: drug_class x extraction_status ===\n")
print(table(es_r50$drug_class, es_r50$extraction_status, useNA = "ifany"))

cat("\n=== CROSS-TAB: drug_class x assay_type ===\n")
print(table(es_r50$drug_class, es_r50$assay_type, useNA = "ifany"))

cat("\n=== CROSS-TAB: extraction_status x assay_type ===\n")
print(table(es_r50$extraction_status, es_r50$assay_type, useNA = "ifany"))

# Helper: get non-NA cluster vector for a given moderator
cluster_nonNA <- function(data, moderator_col) {
  data$study_label[!is.na(data[[moderator_col]])]
}

# =========================================================
# 9a) EXTRACTION STATUS: MAIN EFFECTS
# =========================================================
m_extr <- rma.mv(
  yi, vi,
  mods   = ~ extraction_status,
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data   = es_r50
)

V_extr  <- vcovCR(m_extr, cluster = cluster_nonNA(es_r50, "extraction_status"), type = "CR2")
ct_extr <- coef_test(m_extr, vcov = V_extr, test = "Satterthwaite")
ci_extr <- conf_int(m_extr,  vcov = V_extr)

cat("\n=== EXTRACTION MODERATOR: coefficient tests ===\n")
print(ct_extr)
print(ci_extr)

n_extr_coefs <- length(coef(m_extr))
omni_extr <- Wald_test(
  m_extr,
  constraints = constrain_zero(2:n_extr_coefs),
  vcov        = V_extr,
  test        = "HTZ"
)
cat("\n=== EXTRACTION OMNIBUS WALD TEST (HTZ) ===\n")
print(omni_extr)

# No-intercept means by extraction
m_extr_means <- rma.mv(
  yi, vi,
  mods   = ~ extraction_status - 1,
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data   = es_r50
)

V_extr_means  <- vcovCR(m_extr_means,
                         cluster = cluster_nonNA(es_r50, "extraction_status"),
                         type = "CR2")
ct_extr_means <- coef_test(m_extr_means, vcov = V_extr_means, test = "Satterthwaite")
ci_extr_means <- conf_int(m_extr_means,  vcov = V_extr_means)

extr_means_tbl <- data.frame(
  extraction = gsub("^extraction_status", "", rownames(ct_extr_means)),
  est   = round(ct_extr_means$beta, 3),
  se    = round(ct_extr_means$SE,   3),
  df    = round(ct_extr_means$df_Satt, 2),
  p     = round(ct_extr_means$p_Satt,  3),
  ci_lb = round(ci_extr_means$CI_L, 3),
  ci_ub = round(ci_extr_means$CI_U, 3)
)

cat("\n=== EXTRACTION STATUS MEANS (no-intercept model) ===\n")
print(extr_means_tbl, row.names = FALSE)

# =========================================================
# 9b) ASSAY TYPE: MAIN EFFECTS
# =========================================================
m_assay <- rma.mv(
  yi, vi,
  mods   = ~ assay_type,
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data   = es_r50
)

V_assay  <- vcovCR(m_assay, cluster = cluster_nonNA(es_r50, "assay_type"), type = "CR2")
ct_assay <- coef_test(m_assay, vcov = V_assay, test = "Satterthwaite")
ci_assay <- conf_int(m_assay,  vcov = V_assay)

cat("\n=== ASSAY TYPE MODERATOR: coefficient tests ===\n")
print(ct_assay)
print(ci_assay)

n_assay_coefs <- length(coef(m_assay))
omni_assay <- Wald_test(
  m_assay,
  constraints = constrain_zero(2:n_assay_coefs),
  vcov        = V_assay,
  test        = "HTZ"
)
cat("\n=== ASSAY TYPE OMNIBUS WALD TEST (HTZ) ===\n")
print(omni_assay)

# No-intercept means by assay type
m_assay_means <- rma.mv(
  yi, vi,
  mods   = ~ assay_type - 1,
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data   = es_r50
)

V_assay_means  <- vcovCR(m_assay_means,
                          cluster = cluster_nonNA(es_r50, "assay_type"),
                          type = "CR2")
ct_assay_means <- coef_test(m_assay_means, vcov = V_assay_means, test = "Satterthwaite")
ci_assay_means <- conf_int(m_assay_means,  vcov = V_assay_means)

assay_means_tbl <- data.frame(
  assay = gsub("^assay_type", "", rownames(ct_assay_means)),
  est   = round(ct_assay_means$beta, 3),
  se    = round(ct_assay_means$SE,   3),
  df    = round(ct_assay_means$df_Satt, 2),
  p     = round(ct_assay_means$p_Satt,  3),
  ci_lb = round(ci_assay_means$CI_L, 3),
  ci_ub = round(ci_assay_means$CI_U, 3)
)

cat("\n=== ASSAY TYPE MEANS (no-intercept model) ===\n")
print(assay_means_tbl, row.names = FALSE)

# =========================================================
# 9c) COMBINED MODEL: drug class + extraction + assay type
# =========================================================

# Subset to rows with both moderators non-missing
es_meth <- es_r50[!is.na(es_r50$extraction_status) & !is.na(es_r50$assay_type), ]
cat("\nCombined moderator dataset: k =", nrow(es_meth), "effects from",
    n_distinct(es_meth$study_label), "studies\n")

m_combined <- rma.mv(
  yi, vi,
  mods   = ~ factor(drug_class) + extraction_status + assay_type,
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data   = es_meth
)

V_comb  <- vcovCR(m_combined, cluster = es_meth$study_label, type = "CR2")
ct_comb <- coef_test(m_combined, vcov = V_comb, test = "Satterthwaite")
ci_comb <- conf_int(m_combined,  vcov = V_comb)

cat("\n=== COMBINED MODEL (drug class + extraction + assay): coefficients ===\n")
print(ct_comb)
print(ci_comb)

# Omnibus for extraction in combined model (dynamic index)
coef_nms   <- names(coef(m_combined))
extr_idx   <- grep("^extraction_status", coef_nms)
assay_idx  <- grep("^assay_type",        coef_nms)

if (length(extr_idx) > 0) {
  omni_extr_comb <- Wald_test(
    m_combined,
    constraints = constrain_zero(extr_idx),
    vcov        = V_comb,
    test        = "HTZ"
  )
  cat("\n=== OMNIBUS: EXTRACTION IN COMBINED MODEL (HTZ) ===\n")
  print(omni_extr_comb)
}

if (length(assay_idx) > 0) {
  omni_assay_comb <- Wald_test(
    m_combined,
    constraints = constrain_zero(assay_idx),
    vcov        = V_comb,
    test        = "HTZ"
  )
  cat("\n=== OMNIBUS: ASSAY TYPE IN COMBINED MODEL (HTZ) ===\n")
  print(omni_assay_comb)
}

# =========================================================
# 9d) VARIANCE COMPONENT COMPARISON
# =========================================================
cat("\n=== VARIANCE COMPONENTS: drug-class-only vs combined ===\n")

# Refit drug-class model on same subset for fair comparison
m_class_sub <- rma.mv(
  yi, vi,
  mods   = ~ factor(drug_class) - 1,
  random = list(~1 | study_label, ~1 | sample_id_u, ~1 | effect_id),
  method = "REML",
  data   = es_meth
)

vc_tbl <- data.frame(
  model  = rep(c("drug class only (matched subset)",
                 "drug class + extraction + assay"), each = 3),
  level  = rep(c("study","sample","effect"), 2),
  sigma2 = round(c(m_class_sub$sigma2, m_combined$sigma2), 5)
)
print(vc_tbl, row.names = FALSE)

pct_red_study  <- (m_class_sub$sigma2[1] - m_combined$sigma2[1]) / m_class_sub$sigma2[1] * 100
pct_red_effect <- (m_class_sub$sigma2[3] - m_combined$sigma2[3]) / m_class_sub$sigma2[3] * 100

cat(sprintf("\nBetween-study variance reduction (methodological moderators): %.1f%%\n",
            pct_red_study))
cat(sprintf("Within-study variance reduction (methodological moderators): %.1f%%\n",
            pct_red_effect))

# =========================================================
# 10) EXPORT AUDIT FILES
# =========================================================
write.csv(dat_active,    "dat_active_clean.csv",         row.names = FALSE)
write.csv(class_counts,  "class_counts.csv",              row.names = FALSE)
write.csv(class_results, "drug_class_results.csv",        row.names = FALSE)
write.csv(sens_table,    "sensitivity_results.csv",       row.names = FALSE)
write.csv(collapsed,     "study_collapsed_es.csv",        row.names = FALSE)
write.csv(pw,            "pairwise_contrasts_BH.csv",     row.names = FALSE)
write.csv(extr_means_tbl,"extraction_means.csv",          row.names = FALSE)
write.csv(assay_means_tbl,"assay_means.csv",              row.names = FALSE)

save.image("OT_metaanalysis_clean_workspace.RData")

