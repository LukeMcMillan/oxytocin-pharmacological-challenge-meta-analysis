# =========================================================
# TABLES AND FIGURES
# All statistics drawn from existing workspace objects only.
# Do NOT recalculate. If an object is missing the script
# will stop and tell you which one.
# =========================================================

# --- Install/load packages ---
for (pkg in c("gt", "webshot2", "ggplot2", "dplyr", "stringr", "metafor",
              "scales", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(gt)
library(webshot2)
library(ggplot2)
library(dplyr)
library(stringr)
library(metafor)
library(scales)

# --- Check required objects ---
required_objects <- c(
  "sens_table", "class_results", "class_counts", "pw", "osm_rows",
  "extr_means_tbl", "assay_means_tbl", "collapsed", "m_coll", "egger",
  "es_r50", "ct_primary", "ci_primary", "pred_primary"
)
missing_objects <- required_objects[!sapply(required_objects, exists)]
if (length(missing_objects) > 0) {
  stop(
    "The following required objects are missing from the workspace:\n  ",
    paste(missing_objects, collapse = ", "),
    "\nLoad OT_metaanalysis_clean_workspace.RData before running this script."
  )
}
cat("All required objects found.\n")

# --- Create output directories ---
dir.create("tables",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# =========================================================
# TABLE 1: SENSITIVITY ANALYSES
# =========================================================
tbl1_data <- sens_table |>
  mutate(
    r_label = case_when(
      r == 0.50 ~ ".50 (primary)",
      r == 0.30 ~ ".30",
      r == 0.70 ~ ".70",
      TRUE      ~ as.character(r)
    ),
    smcc   = round(est, 2),
    ci_str = paste0("[", round(ci_lb, 2), ", ", round(ci_ub, 2), "]"),
    pi_str = paste0("[", round(pi_lb, 2), ", ", round(pi_ub, 2), "]")
  ) |>
  select(r_label, smcc, ci_str, pi_str)

tbl1 <- tbl1_data |>
  gt() |>
  cols_label(
    r_label = html("<i>r</i>"),
    smcc    = "SMCC",
    ci_str  = "95% CR2 CI",
    pi_str  = "95% PI"
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = r_label == ".50 (primary)")
  ) |>
  tab_source_note(
    source_note = md(
      "CR2 = cluster-robust variance estimator with small-sample correction; PI = prediction interval. Primary analysis used *r* = .50. Sensitivity analyses at *r* = .30 and *r* = .70 were preregistered."
    )
  ) |>
  tab_options(table.font.size = 11, table.width = pct(70)) |>
  opt_horizontal_padding(scale = 2)

gtsave(tbl1, "tables/Table1_sensitivity.png", expand = 20)
write.csv(tbl1_data, "tables/Table1_sensitivity.csv", row.names = FALSE)
cat("Table 1 saved.\n")


# =========================================================
# TABLE 2: DRUG CLASS MEAN EFFECTS
# =========================================================
dc_order  <- c("metabolic","monoaminergic","other","osmotic","hormonal","opioid")
dc_labels <- c("Metabolic","Monoaminergic","Other","Osmotic (ref.)","Hormonal","Opioid")

# Merge means + k
tbl2_base <- class_results |>
  left_join(class_counts |> select(drug_class, k_effects), by = "drug_class") |>
  mutate(
    drug_class_lower = drug_class,
    dc_display = dc_labels[match(drug_class, dc_order)],
    smcc_fmt   = round(est, 2),
    ci_str     = paste0("[", round(ci_lb, 2), ", ", round(ci_ub, 2), "]"),
    t_str      = round(est / se, 2),
    df_str     = round(df, 2),
    t_df_str   = paste0(round(est / se, 2), " (", round(df, 2), ")"),
    p_fmt      = round(p, 3),
    sig_stars  = case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "")
  )

# Extract contrasts vs osmotic from osm_rows
contrast_lookup <- osm_rows |>
  mutate(
    non_osm = if_else(lev1 == "osmotic", lev2, lev1),
    delta_display = if_else(lev1 == "osmotic", -delta, delta)
  ) |>
  select(non_osm, delta_display, p_raw, p_BH)

tbl2_data <- tbl2_base |>
  left_join(contrast_lookup, by = c("drug_class_lower" = "non_osm")) |>
  mutate(
    delta_fmt = if_else(!is.na(delta_display), sprintf("%.2f", delta_display), "—"),
    p_raw_c   = if_else(!is.na(p_raw),  sprintf("%.3f", p_raw),  "—"),
    p_BH_c    = if_else(!is.na(p_BH),   sprintf("%.3f", p_BH),   "—"),
    order_col = match(drug_class_lower, dc_order)
  ) |>
  arrange(order_col) |>
  select(dc_display, k_effects, smcc_fmt, ci_str, t_df_str, p_fmt, sig_stars,
         delta_fmt, p_raw_c, p_BH_c)

tbl2 <- tbl2_data |>
  gt() |>
  cols_label(
    dc_display = "Drug Class",
    k_effects  = html("<i>k</i>"),
    smcc_fmt   = "SMCC",
    ci_str     = "95% CI",
    t_df_str   = html("<i>t</i> (df)"),
    p_fmt      = html("<i>p</i>"),
    sig_stars  = "",
    delta_fmt  = html("&Delta; vs Osmotic"),
    p_raw_c    = html("<i>p</i><sub>raw</sub>"),
    p_BH_c     = html("<i>p</i><sub>BH</sub>")
  ) |>
  tab_style(
    style     = cell_fill(color = "#EFEFEF"),
    locations = cells_body(rows = dc_display == "Osmotic (ref.)")
  ) |>
  tab_style(
    style     = cell_text(style = "italic"),
    locations = cells_body(columns = sig_stars)
  ) |>
  tab_source_note(
    source_note = md(
      "SMCC = standardised mean change coefficient. All CIs are 95% CR2 confidence intervals. \u0394 vs Osmotic = pairwise contrast against osmotic reference. *p*-values for \u0394 vs Osmotic are uncorrected; BH-corrected *p*-values reported from the full 15-contrast correction. No contrasts survived BH correction at \u03b1 = .05. * *p* < .05, ** *p* < .01, *** *p* < .001."
    )
  ) |>
  tab_options(table.font.size = 11) |>
  opt_horizontal_padding(scale = 2)

gtsave(tbl2, "tables/Table2_drugclass.png", expand = 20)
write.csv(tbl2_data, "tables/Table2_drugclass.csv", row.names = FALSE)
cat("Table 2 saved.\n")


# =========================================================
# TABLE 3: BH-CORRECTED PAIRWISE CONTRASTS
# =========================================================
tbl3_data <- pw |>
  arrange(p_raw) |>
  mutate(
    delta_fmt = sprintf("%.2f", delta),
    se_fmt    = sprintf("%.2f", se),
    t_fmt     = sprintf("%.2f", t),
    df_fmt    = sprintf("%.2f", df),
    p_raw_fmt = sprintf("%.3f", p_raw),
    p_BH_fmt  = sprintf("%.3f", p_BH),
    sig_raw   = if_else(p_raw < .05, "*", ""),
    sig_BH    = if_else(p_BH  < .05, "*", "")
  ) |>
  select(contrast, delta_fmt, se_fmt, t_fmt, df_fmt,
         p_raw_fmt, sig_raw, p_BH_fmt, sig_BH)

tbl3 <- tbl3_data |>
  gt() |>
  cols_label(
    contrast  = "Contrast",
    delta_fmt = html("&Delta;SMCC"),
    se_fmt    = "SE",
    t_fmt     = html("<i>t</i>"),
    df_fmt    = "df",
    p_raw_fmt = html("<i>p</i><sub>raw</sub>"),
    sig_raw   = html("<i>p</i><sub>raw</sub><br>sig."),
    p_BH_fmt  = html("<i>p</i><sub>BH</sub>"),
    sig_BH    = html("<i>p</i><sub>BH</sub><br>sig.")
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = as.numeric(p_raw_fmt) < .05)
  ) |>
  tab_source_note(
    source_note = md(
      "All contrasts computed from the no-intercept class means model with CR2 cluster-robust standard errors and Satterthwaite-adjusted degrees of freedom. *p*<sub>raw</sub> = uncorrected *p*-value; *p*<sub>BH</sub> = Benjamini\u2013Hochberg FDR-corrected *p*-value across all pairwise contrasts. No contrasts survived BH correction at \u03b1 = .05. Contrasts sorted by raw *p*-value."
    )
  ) |>
  tab_options(table.font.size = 11) |>
  opt_horizontal_padding(scale = 2)

gtsave(tbl3, "tables/Table3_pairwise.png", expand = 20)
write.csv(tbl3_data, "tables/Table3_pairwise.csv", row.names = FALSE)
cat("Table 3 saved.\n")


# =========================================================
# TABLE 3a: EXTRACTION STATUS MEANS
# =========================================================
# k values supplied by user (not_extracted=8, extracted=40, not_reported=42)
k_lookup <- c(not_extracted = 8, extracted = 40, not_reported = 42)

tbl3a_data <- extr_means_tbl |>
  mutate(
    k_val      = k_lookup[trimws(extraction)],
    extr_label = case_when(
      trimws(extraction) == "not_extracted" ~ "Not Extracted",
      trimws(extraction) == "extracted"     ~ "Extracted",
      trimws(extraction) == "not_reported"  ~ "Not Reported",
      TRUE ~ str_to_title(trimws(extraction))
    ),
    smcc_fmt = sprintf("%.2f", est),
    se_fmt   = sprintf("%.2f", se),
    t_df_str = paste0(sprintf("%.2f", est / se), " (", sprintf("%.2f", df), ")"),
    p_fmt    = sprintf("%.3f", p),
    ci_str   = paste0("[", sprintf("%.2f", ci_lb), ", ", sprintf("%.2f", ci_ub), "]")
  ) |>
  select(extr_label, k_val, smcc_fmt, se_fmt, t_df_str, p_fmt, ci_str)

tbl3a <- tbl3a_data |>
  gt() |>
  cols_label(
    extr_label = "Extraction Status",
    k_val      = html("<i>k</i>"),
    smcc_fmt   = "SMCC",
    se_fmt     = "SE",
    t_df_str   = html("<i>t</i> (df)"),
    p_fmt      = html("<i>p</i>"),
    ci_str     = "95% CI"
  ) |>
  tab_source_note(
    source_note = md(
      "Two effect sizes with missing extraction information excluded. SMCC = standardised mean change coefficient. All CIs are 95% CR2 confidence intervals with Satterthwaite-adjusted degrees of freedom. Omnibus Wald test for extraction status: *F*(2, 10.3) = 3.55, *p* = .067. Combined model (drug class + extraction + assay type): extraction *F*(2, 7.9) = 1.37, *p* = .308; assay type *F*(2, 2.65) = 1.74, *p* = .329."
    )
  ) |>
  tab_options(table.font.size = 11, table.width = pct(80)) |>
  opt_horizontal_padding(scale = 2)

gtsave(tbl3a, "tables/Table3a_extraction.png", expand = 20)
write.csv(tbl3a_data, "tables/Table3a_extraction.csv", row.names = FALSE)
cat("Table 3a saved.\n")


# =========================================================
# FIGURE 2: DRUG CLASS HORIZONTAL BAR CHART
# =========================================================
fig2_order <- c("Metabolic","Monoaminergic","Other","Osmotic","Hormonal","Opioid")

fig2_data <- class_results |>
  left_join(class_counts |> select(drug_class, k_effects), by = "drug_class") |>
  mutate(
    dc_label  = case_when(
      drug_class == "metabolic"     ~ "Metabolic",
      drug_class == "monoaminergic" ~ "Monoaminergic",
      drug_class == "other"         ~ "Other",
      drug_class == "osmotic"       ~ "Osmotic",
      drug_class == "hormonal"      ~ "Hormonal",
      drug_class == "opioid"        ~ "Opioid",
      TRUE ~ str_to_title(drug_class)
    ),
    sig_fill  = if_else(p < .05, "sig", "nonsig"),
    sig_stars = case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ ""),
    bar_label = paste0("k = ", k_effects, if_else(sig_stars != "", paste0("  ", sig_stars), "")),
    ci_lb_plot = pmax(ci_lb, -1.5),
    ci_ub_plot = pmin(ci_ub,  2.5),
    is_opioid = drug_class == "opioid",
    dc_label  = factor(dc_label, levels = rev(fig2_order))
  )

# Label x position: just past the clipped CI upper bound
fig2_data <- fig2_data |>
  mutate(label_x = ci_ub_plot + 0.07)

fig2 <- ggplot(fig2_data, aes(x = est, y = dc_label, fill = sig_fill)) +
  geom_col(width = 0.55, colour = NA) +
  geom_errorbarh(
    aes(xmin = ci_lb_plot, xmax = ci_ub_plot),
    height = 0.2, linewidth = 0.7, colour = "grey25"
  ) +
  # Opioid double-headed arrow to indicate clipping
  geom_segment(
    data = filter(fig2_data, is_opioid),
    aes(x = 2.3, xend = 2.5, y = dc_label, yend = dc_label),
    arrow = arrow(ends = "both", length = unit(0.15, "cm"), type = "open"),
    colour = "grey40", inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.6) +
  geom_text(
    aes(x = label_x, label = bar_label),
    hjust = 0, size = 3.2, colour = "grey20"
  ) +
  scale_fill_manual(
    values = c(sig = "#1B4F8A", nonsig = "#7EC8C8"),
    guide  = "none"
  ) +
  scale_x_continuous(
    limits = c(-1.5, 3.2),
    breaks = seq(-1.5, 2.5, 0.5),
    labels = seq(-1.5, 2.5, 0.5),
    expand = c(0, 0)
  ) +
  labs(
    x = "Standardised Mean Change Coefficient (SMCC)",
    y = NULL,
    caption = paste0(
      "Figure 2. Drug class standardised mean change coefficients (SMCC) with 95% CR2 confidence intervals\n",
      "from the no-intercept means model (k = ", nrow(es_r50), "). Dark blue bars indicate classes with individually\n",
      "significant effects (p < .05); pale teal bars indicate non-significant classes. Double-headed arrow\n",
      "for Opioid indicates CI extends beyond axis limits [\u22122.72, 2.64]."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y    = element_blank(),
    panel.grid.minor      = element_blank(),
    panel.grid.major.x    = element_line(colour = "grey90", linewidth = 0.4),
    axis.text             = element_text(size = 11, colour = "grey15"),
    plot.caption          = element_text(hjust = 0, size = 8.5, colour = "grey35",
                                          margin = margin(t = 10)),
    plot.caption.position = "plot",
    plot.margin           = margin(10, 20, 10, 10)
  )

ggsave("figures/Figure2_drugclass.png", fig2,
       width = 18, height = 12, units = "cm", dpi = 300)
ggsave("figures/Figure2_drugclass.pdf", fig2,
       width = 18, height = 12, units = "cm")
cat("Figure 2 saved.\n")


# =========================================================
# FIGURE 3: CONTOUR-ENHANCED FUNNEL PLOT
# =========================================================
# Pull Egger stats from egger object
egger_t   <- round(egger$zval, 2)
egger_df  <- egger$dfs
egger_p   <- egger$pval
egger_p_fmt <- if (egger_p < .001) "< .001" else
  paste0("= ", sub("^0", "", sprintf("%.3f", egger_p)))

# Identify Fisher (1987)
fisher_idx <- which(grepl("Fisher", collapsed$study_label, ignore.case = TRUE) &
                      grepl("1987",   collapsed$study_label))

for (ext in c("png", "pdf")) {
  if (ext == "png") {
    png("figures/Figure3_funnel.png",
        width = 18, height = 16, units = "cm", res = 300, bg = "white")
  } else {
    pdf("figures/Figure3_funnel.pdf",
        width = 18 / 2.54, height = 16 / 2.54)
  }

  funnel(
    m_coll,
    contour     = c(0.90, 0.95, 0.99),
    col.contour = c("grey85", "grey65", "grey40"),
    xlab        = "Standardised Mean Change Coefficient (SMCC)",
    ylab        = "Standard Error",
    back        = "white",
    shade       = NULL,
    legend      = FALSE,
    main        = ""
  )

  # Overlay Fisher (1987) as open triangle
  if (length(fisher_idx) > 0) {
    points(
      collapsed$yi[fisher_idx],
      sqrt(collapsed$vi[fisher_idx]),
      pch = 2, cex = 1.4, lwd = 1.6, col = "black"
    )
    text(
      collapsed$yi[fisher_idx],
      sqrt(collapsed$vi[fisher_idx]),
      labels = "Fisher (1987)",
      pos = 3, cex = 0.72, offset = 0.5
    )
  }

  # Egger annotation
  legend(
    "topleft",
    legend = paste0("Egger's: t(", egger_df, ") = ", egger_t,
                    ", p ", egger_p_fmt),
    bty = "n", cex = 0.82
  )

  # Contour legend
  legend(
    "topright",
    legend = c("90%", "95%", "99%"),
    fill   = c("grey85", "grey65", "grey40"),
    title  = "Significance\ncontour",
    bty    = "n", cex = 0.78
  )

  dev.off()
}
cat("Figure 3 saved.\n")


# =========================================================
# FIGURE 4: FOREST PLOT
# =========================================================
# Sort studies by effect size ascending
ord          <- order(collapsed$yi)
sorted_yi    <- collapsed$yi[ord]
sorted_vi    <- collapsed$vi[ord]
sorted_slab  <- collapsed$study_label[ord]
n_k          <- length(sorted_yi)

# Drug class lookup: one canonical class per study label
drug_class_lookup <- c(
  "Amico 1985"             = "Hormonal",
  "Asla 2025"              = "Other",
  "Atila 2023"             = "Monoaminergic",
  "Atila 2024"             = "Metabolic",
  "Chiodera 1991"          = "Metabolic",
  "Chiodera 1992 (ITT)"    = "Metabolic",
  "Chiodera 1994"          = "Other",
  "Chiodera 1995"          = "Hormonal",
  "Chiodera 1996"          = "Metabolic",
  "Chiodera 1998"          = "Osmotic",
  "Coiro 1991"             = "Metabolic",
  "Coiro 1991 (ANG)"       = "Osmotic",
  "Corbett 2016"           = "Hormonal",
  "Dolder 2017"            = "Osmotic",
  "Dolder 2018"            = "Monoaminergic",
  "Dumont 2009"            = "Monoaminergic",
  "EKSTRÖM 1992"           = "Hormonal",
  "Fisher 1987"            = "Metabolic",
  "Forsling 1999"          = "Other",
  "Forsling 2002"          = "Osmotic",
  "Gabay 2019"             = "Monoaminergic",
  "Galbiati 2025"          = "Hormonal",
  "Johnson 1990"           = "Metabolic",
  "Johnston 2014"          = "Osmotic",
  "Kirkpatrick 2014"       = "Monoaminergic",
  "Kuypers 2014"           = "Monoaminergic",
  "Laczi 1998"             = "Osmotic",
  "Lee 2003"               = "Monoaminergic",
  "Ohlsson 2002"           = "Metabolic",
  "Ottesen 1988"           = "Other",
  "Radant 1992"            = "Monoaminergic",
  "Sailer 2021 study 1"    = "Osmotic",
  "Sailer 2021 study 2"    = "Metabolic",
  "Sailer 2021 study 3"    = "Metabolic",
  "Schmid 2015"            = "Monoaminergic",
  "Srinavasa/Aulinas 2018" = "Osmotic",
  "Vizeli 2022"            = "Monoaminergic",
  "Williams 1985 study 1"  = "Osmotic",
  "Williams 1985 study 2"  = "Osmotic",
  "Williams 1985 study 3"  = "Other"
)
sorted_drug <- drug_class_lookup[sorted_slab]
sorted_drug[is.na(sorted_drug)] <- ""

# Pooled estimate and PI from m_coll
pooled_est  <- as.numeric(coef(m_coll))
pooled_lb   <- m_coll$ci.lb
pooled_ub   <- m_coll$ci.ub
pred_coll   <- predict(m_coll)
pi_lb_val   <- pred_coll$pi.lb
pi_ub_val   <- pred_coll$pi.ub

for (ext in c("png", "pdf")) {
  if (ext == "png") {
    png("figures/Figure4_forest.png",
        width = 30, height = 32, units = "cm", res = 300, bg = "white")
  } else {
    pdf("figures/Figure4_forest.pdf",
        width = 30 / 2.54, height = 32 / 2.54)
  }

  op <- par(mar = c(4, 0, 1, 0))

  forest(
    x          = sorted_yi,
    vi         = sorted_vi,
    slab       = sorted_slab,
    ilab       = sorted_drug,        # drug name column
    ilab.xpos  = -3.2,               # position: between study label and forest
    ilab.pos   = 4,                  # left-aligned text
    xlim       = c(-13, 9),          # extra left room for study + drug columns
    alim       = c(-2, 5),
    at         = seq(-2, 5, 1),
    cex        = 0.68,
    xlab       = "Standardised Mean Change Coefficient (SMCC)",
    header     = c("Study", "SMCC [95% CI]"),
    refline    = 0,
    addfit     = FALSE,
    top        = 3
  )

  # Drug class column header (manually placed to match ilab.xpos)
  text(-3.2, n_k + 3, "Drug Class", font = 2, cex = 0.68, pos = 4)

  # Add pooled RE diamond
  addpoly(
    x     = pooled_est,
    ci.lb = pooled_lb,
    ci.ub = pooled_ub,
    rows  = -1.5,
    mlab  = paste0("Pooled RE (k = ", n_k, ", SMCC = ",
                   round(pooled_est, 2), ")"),
    cex   = 0.72,
    col   = "black"
  )

  # Prediction interval dotted lines
  par(xpd = TRUE)
  segments(
    x0  = pi_lb_val, y0 = -2,
    x1  = pi_lb_val, y1 = n_k + 3,
    lty = 3, lwd = 1.2, col = "grey45"
  )
  segments(
    x0  = pi_ub_val, y0 = -2,
    x1  = pi_ub_val, y1 = n_k + 3,
    lty = 3, lwd = 1.2, col = "grey45"
  )

  # Fisher (1987) upper CI arrowhead if it exceeds axis
  if (length(fisher_idx) > 0) {
    fisher_row_plot <- which(sorted_slab == collapsed$study_label[fisher_idx])
    if (length(fisher_row_plot) > 0 &&
        collapsed$yi[fisher_idx] + 1.96 * sqrt(collapsed$vi[fisher_idx]) > 5) {
      arrows(
        x0  = 4.9,
        y0  = fisher_row_plot,
        x1  = 5.1,
        y1  = fisher_row_plot,
        length = 0.07, col = "grey30"
      )
    }
  }

  par(op)
  dev.off()
}
cat("Figure 4 saved.\n")


# =========================================================
# OPTIONAL FIGURE: DRUG CLASS + EXTRACTION OVERLAY
# =========================================================
ct_raw <- as.data.frame(
  table(es_r50$drug_class, es_r50$extraction_status),
  stringsAsFactors = FALSE
) |>
  rename(drug_class = Var1, extraction = Var2, count = Freq) |>
  filter(!is.na(extraction), extraction != "") |>
  group_by(drug_class) |>
  mutate(prop = count / sum(count)) |>
  ungroup() |>
  mutate(
    dc_label = case_when(
      drug_class == "metabolic"     ~ "Metabolic",
      drug_class == "monoaminergic" ~ "Monoaminergic",
      drug_class == "other"         ~ "Other",
      drug_class == "osmotic"       ~ "Osmotic",
      drug_class == "hormonal"      ~ "Hormonal",
      drug_class == "opioid"        ~ "Opioid",
      TRUE ~ str_to_title(drug_class)
    ),
    dc_label  = factor(dc_label, levels = rev(fig2_order)),
    extraction = factor(
      extraction,
      levels = c("not_extracted", "extracted", "not_reported"),
      labels = c("Not extracted", "Extracted", "Not reported")
    )
  )

p_main_opt <- ggplot(fig2_data, aes(x = est, y = dc_label, fill = sig_fill)) +
  geom_col(width = 0.5, colour = NA) +
  geom_errorbarh(aes(xmin = ci_lb_plot, xmax = ci_ub_plot),
                 height = 0.2, linewidth = 0.6, colour = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  scale_fill_manual(values = c(sig = "#1B4F8A", nonsig = "#7EC8C8"), guide = "none") +
  scale_x_continuous(limits = c(-1.5, 2.5), breaks = seq(-1.5, 2.5, 0.5)) +
  labs(x = "SMCC", y = NULL, title = "Mean effect by drug class") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

p_extract_opt <- ggplot(ct_raw, aes(x = prop, y = dc_label, fill = extraction)) +
  geom_col(width = 0.5, position = "stack") +
  scale_fill_manual(
    values = c("Not extracted" = "#D73027",
               "Extracted"     = "#4DAC26",
               "Not reported"  = "#BDBDBD"),
    name = "Extraction"
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = c(0, 0.5, 1)) +
  labs(x = "Proportion", y = NULL, title = "Extraction status by class") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_blank()
  )

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig_opt <- p_main_opt + p_extract_opt + plot_layout(widths = c(2, 1))
} else {
  fig_opt <- p_main_opt
  message("Install 'patchwork' for the combined optional figure. Saving main panel only.")
}

ggsave("figures/FigureOpt_drugclass_extraction.png", fig_opt,
       width = 24, height = 12, units = "cm", dpi = 300)
cat("Optional figure saved.\n")


# =========================================================
# FINAL SUMMARY
# =========================================================
cat("\n========================================\n")
cat("All exports complete.\n")
cat("Tables in: tables/\n")
cat("  Table1_sensitivity.png / .csv\n")
cat("  Table2_drugclass.png / .csv\n")
cat("  Table3_pairwise.png / .csv\n")
cat("  Table3a_extraction.png / .csv\n")
cat("Figures in: figures/\n")
cat("  Figure2_drugclass.png / .pdf\n")
cat("  Figure3_funnel.png / .pdf\n")
cat("  Figure4_forest.png / .pdf\n")
cat("  FigureOpt_drugclass_extraction.png\n")
cat("========================================\n")

