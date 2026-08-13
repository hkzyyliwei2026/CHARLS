# Final reporting package for SCI manuscript.
# This script reads the completed SCI analysis outputs and produces publication-facing
# main tables, sensitivity summaries, survey-design notes, and an explicit-knot
# restricted cubic spline figure/specification. It does not overwrite Chinese
# submission analyses.

set.seed(42)

`%||%` <- function(a, b) if (length(a) && !is.na(a) && nzchar(a)) a else b

cmd_args <- commandArgs(FALSE)
file_arg <- cmd_args[grep("--file=", cmd_args)][1] %||% ""
script_path <- sub("^--file=", "", file_arg)
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = FALSE)) else getwd()
if (!dir.exists(file.path(script_dir, "tables"))) script_dir <- getwd()

data_dir <- file.path(script_dir, "data")
table_dir <- file.path(script_dir, "tables")
figure_dir <- file.path(script_dir, "figures")
log_dir <- file.path(script_dir, "reproducibility_logs")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

group_order <- c(
  "persistent_nonhypertensive",
  "incident_hypertension",
  "persistent_unaware_or_untreated",
  "new_awareness_or_treatment_uncontrolled",
  "persistent_aware_or_treated_uncontrolled",
  "gained_control",
  "persistent_controlled",
  "control_loss",
  "other_transition"
)

group_labels <- c(
  persistent_nonhypertensive = "Persistent non-hypertensive",
  incident_hypertension = "Incident hypertension",
  persistent_unaware_or_untreated = "Persistent unaware/untreated hypertension",
  new_awareness_or_treatment_uncontrolled = "New awareness/treatment but uncontrolled",
  persistent_aware_or_treated_uncontrolled = "Persistent awareness or treatment but uncontrolled",
  gained_control = "Gained control",
  persistent_controlled = "Persistent controlled",
  control_loss = "Loss of control",
  other_transition = "Other transitions"
)

fmt <- function(x, digits = 1) sprintf(paste0("%.", digits, "f"), as.numeric(x))
fmt_ci <- function(est, lo, hi, digits = 2) sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), est, lo, hi)
read <- function(name) read.csv(file.path(table_dir, name), stringsAsFactors = FALSE, check.names = FALSE)

event_summary <- read("table_transition_group_event_summary.csv")
event_summary <- event_summary[match(group_order, event_summary$transition_group), ]
event_summary$transition_group_label <- unname(group_labels[event_summary$transition_group])

table1_long <- read("table1_baseline_by_transition_group_long.csv")
table1_long$group <- gsub("Persistent aware/treated but uncontrolled", "Persistent awareness or treatment but uncontrolled", table1_long$group, fixed = TRUE)
table1_long$group <- gsub("Control loss", "Loss of control", table1_long$group, fixed = TRUE)

row_map <- c(
  N = "Participants, n",
  cvd_events = "Incident CVD events, n",
  crude_risk = "Crude CVD risk, %",
  age_2011 = "Age, years",
  "sex: female" = "Female sex, n (%)",
  sbp_2011 = "2011 SBP, mm Hg",
  dbp_2011 = "2011 DBP, mm Hg",
  sbp_2015 = "2015 SBP, mm Hg",
  dbp_2015 = "2015 DBP, mm Hg",
  bmi_2011 = "BMI, kg/m2",
  current_smoking = "Current smoking, n (%)",
  alcohol_last_year = "Alcohol drinking in the past year, n (%)",
  diabetes = "Diabetes, n (%)",
  dyslipidemia = "Dyslipidemia, n (%)",
  "hukou: agricultural" = "Agricultural hukou, n (%)",
  "education: no_formal" = "No formal education, n (%)",
  "marital: married_partnered" = "Married or partnered, n (%)"
)

extra_rows <- do.call(rbind, lapply(seq_len(nrow(event_summary)), function(i) {
  data.frame(
    variable = c("cvd_events", "crude_risk"),
    group = rep(event_summary$transition_group_label[i], 2),
    value = c(
      as.character(as.integer(event_summary$cvd_events[i])),
      fmt(event_summary$event_rate_pct[i], 1)
    ),
    stringsAsFactors = FALSE
  )
}))
table1_main_long <- rbind(table1_long, extra_rows)
table1_main_long <- table1_main_long[table1_main_long$variable %in% names(row_map), ]
table1_main_long$row_label <- unname(row_map[table1_main_long$variable])
table1_main_long$row_order <- match(table1_main_long$variable, names(row_map))
table1_main_long$group <- factor(table1_main_long$group, levels = unname(group_labels[group_order]))
table1_main_long <- table1_main_long[order(table1_main_long$row_order, table1_main_long$group), ]
table1_wide <- reshape(
  table1_main_long[, c("row_label", "group", "value")],
  idvar = "row_label",
  timevar = "group",
  direction = "wide"
)
names(table1_wide) <- sub("^value\\.", "", names(table1_wide))
table1_wide <- table1_wide[match(unname(row_map), table1_wide$row_label), ]
names(table1_wide)[1] <- "Characteristic"
write.csv(table1_wide, file.path(table_dir, "main_table1_baseline_by_9_transition_groups.csv"), row.names = FALSE)

logistic <- read("table2_logistic_or_all_models.csv")
logistic$transition_group_label <- unname(group_labels[logistic$transition_group])
model_wide <- logistic[logistic$outcome == "incident_cvd_2015_2018", c("transition_group", "transition_group_label", "model", "estimate_ci")]
model_wide <- reshape(model_wide, idvar = c("transition_group", "transition_group_label"), timevar = "model", direction = "wide")
names(model_wide) <- sub("^estimate_ci\\.", "", names(model_wide))
model_wide <- model_wide[match(group_order[-1], model_wide$transition_group), ]
event_lookup <- event_summary[, c("transition_group", "n", "cvd_events")]
table2_main <- merge(model_wide, event_lookup, by = "transition_group", all.x = TRUE, sort = FALSE)
table2_main <- table2_main[match(group_order[-1], table2_main$transition_group), ]
table2_main$`Events/participants` <- sprintf("%d/%d", as.integer(table2_main$cvd_events), as.integer(table2_main$n))
table2_main <- table2_main[, c("transition_group_label", "Events/participants", "model1_crude", "model2_demographic", "model3_full")]
names(table2_main) <- c("Transition group", "Events/participants", "Crude OR (95% CI)", "Model 2 OR (95% CI)", "Fully adjusted OR (95% CI)")
write.csv(table2_main, file.path(table_dir, "main_table2_logistic_three_models.csv"), row.names = FALSE)

rr <- read("table3_modified_poisson_rr.csv")
rr$transition_group_label <- unname(group_labels[rr$transition_group])
ar <- read("table3_adjusted_absolute_risk_and_rd.csv")
ar$transition_group_label <- unname(group_labels[ar$transition_group])
table3_main <- merge(
  rr[, c("transition_group", "estimate_ci")],
  ar[, c("transition_group", "adjusted_risk_pct", "adjusted_risk_ci_lower_pct", "adjusted_risk_ci_upper_pct", "risk_difference_per_1000", "rd_ci_lower_per_1000", "rd_ci_upper_per_1000")],
  by = "transition_group",
  all.y = TRUE,
  sort = FALSE
)
table3_main$transition_group_label <- unname(group_labels[table3_main$transition_group])
table3_main <- table3_main[match(group_order, table3_main$transition_group), ]
table3_main$`Adjusted RR (95% CI)` <- ifelse(is.na(table3_main$estimate_ci), "1.00 (reference)", table3_main$estimate_ci)
table3_main$`Adjusted risk per 1000 (95% CI)` <- sprintf(
  "%.1f (%.1f-%.1f)",
  table3_main$adjusted_risk_pct * 10,
  table3_main$adjusted_risk_ci_lower_pct * 10,
  table3_main$adjusted_risk_ci_upper_pct * 10
)
table3_main$`Risk difference per 1000 (95% CI)` <- sprintf(
  "%.1f (%.1f-%.1f)",
  table3_main$risk_difference_per_1000,
  table3_main$rd_ci_lower_per_1000,
  table3_main$rd_ci_upper_per_1000
)
table3_main <- table3_main[, c("transition_group_label", "Adjusted RR (95% CI)", "Adjusted risk per 1000 (95% CI)", "Risk difference per 1000 (95% CI)")]
names(table3_main)[1] <- "Transition group"
write.csv(table3_main, file.path(table_dir, "main_table3_rr_absolute_risk.csv"), row.names = FALSE)

event_flags <- event_summary[, c("transition_group_label", "n", "cvd_events", "heart_events", "stroke_events")]
event_flags$cvd_event_interpretation <- ifelse(
  event_flags$cvd_events < 10, "Do not model independently",
  ifelse(event_flags$cvd_events < 20, "Exploratory interpretation", "Retained in fully adjusted model")
)
write.csv(event_flags, file.path(table_dir, "supp_table_exposure_event_count_flags.csv"), row.names = FALSE)

core <- "persistent_aware_or_treated_uncontrolled"
pull_effect <- function(file, analysis, measure = "OR", model = NULL) {
  x <- read(file)
  if (!is.null(model) && "model" %in% names(x)) x <- x[x$model == model, ]
  row <- x[x$transition_group == core, ][1, ]
  data.frame(
    analysis = analysis,
    measure = measure,
    effect_95ci = row$estimate_ci,
    p = row$p_text,
    n_model = row$n_model,
    events_model = row$events_model,
    stringsAsFactors = FALSE
  )
}
sens <- rbind(
  pull_effect("table2_logistic_or_all_models.csv", "Complete-case fully adjusted logistic", "OR", model = "model3_full"),
  pull_effect("table3_modified_poisson_rr.csv", "Modified Poisson with robust standard errors", "RR"),
  {
    x <- read("supp_table_multiple_imputation_model3_or.csv")
    row <- x[x$transition_group == core, ][1, ]
    data.frame(analysis = "Multiple imputation by chained equations", measure = "OR", effect_95ci = row$estimate_ci, p = row$p_text, n_model = NA, events_model = NA)
  },
  pull_effect("supp_table_sensitivity_130_80_or.csv", "Alternative hypertension threshold of 130/80 mm Hg", "OR"),
  pull_effect("supp_table_sensitivity_strict_2016_2018_or.csv", "Strict incident CVD definition restricted to 2016-2018", "OR"),
  pull_effect("supp_table_cluster_robust_psu_or.csv", "PSU cluster-robust standard errors", "OR"),
  pull_effect("supp_table_survey_design_weighted_or.csv", "Survey-design model with baseline biomarker weight", "OR"),
  pull_effect("supp_table_ipcw_or.csv", "Inverse-probability-of-censoring weighting", "OR")
)
write.csv(sens, file.path(table_dir, "supp_table_core_exposure_sensitivity_summary.csv"), row.names = FALSE)

survey_details <- data.frame(
  item = c(
    "Sampling weight",
    "Rationale for weight choice",
    "Weight scaling",
    "Weight truncation",
    "PSU",
    "Strata",
    "Finite population correction",
    "Survey design specification",
    "IPCW combination",
    "Interpretation"
  ),
  value = c(
    "2011 baseline biomarker weight: bio_weight2, stored as baseline_biomarker_weight",
    "The exposure definition used measured blood pressure from the biomarker examination; therefore the biomarker/examination weight was more relevant than household or interview-only weights.",
    "Weights were standardized by dividing by their analytic-sample mean before modeling.",
    "Survey sensitivity analyses used weights truncated at the 1st and 99th percentiles.",
    "communityID from the CHARLS 2011 PSU/weight files, stored as cluster_psu.",
    "Proxy strata constructed as province_eng combined with urban_nbs.",
    "No finite population correction was specified.",
    "survey::svydesign(id = ~cluster_psu, strata = ~strata_proxy, weights = ~w_baseline_biomarker_trunc, nest = TRUE), with lonely PSUs adjusted.",
    "IPCW was analyzed as a separate sensitivity analysis and was not multiplied by the baseline sampling weight.",
    "Because no validated three-wave longitudinal analysis weight was available, this analysis should be described as survey-weighted sensitivity analysis rather than as a fully nationally representative estimate."
  ),
  stringsAsFactors = FALSE
)
write.csv(survey_details, file.path(table_dir, "supp_table_survey_weight_design_details.csv"), row.names = FALSE)

# Explicit-knot restricted cubic spline for cumulative mean SBP.
df <- read.csv(file.path(data_dir, "sci_landmark_analysis_dataset.csv"), stringsAsFactors = FALSE)
df$transition_group <- factor(df$transition_group, levels = group_order)
df$sex <- factor(df$sex)
df$hukou <- factor(df$hukou)
df$education <- factor(df$education)
df$marital <- factor(df$marital)
vars <- c("incident_cvd_2015_2018", "mean_sbp_2011_2015", "age_per10", "sex", "hukou", "education", "marital", "bmi_2011", "current_smoking", "alcohol_last_year", "diabetes", "dyslipidemia")
dat <- df[complete.cases(df[, vars]), vars]
kn <- as.numeric(quantile(dat$mean_sbp_2011_2015, probs = c(0.05, 0.35, 0.65, 0.95), na.rm = TRUE))
names(kn) <- c("p5", "p35", "p65", "p95")
rhs_cov <- "age_per10 + sex + hukou + education + marital + bmi_2011 + current_smoking + alcohol_last_year + diabetes + dyslipidemia"
base_fit <- glm(as.formula(paste("incident_cvd_2015_2018 ~", rhs_cov)), data = dat, family = binomial())
linear_fit <- glm(as.formula(paste("incident_cvd_2015_2018 ~ mean_sbp_2011_2015 +", rhs_cov)), data = dat, family = binomial())
spline_fit <- glm(
  as.formula(paste("incident_cvd_2015_2018 ~ splines::ns(mean_sbp_2011_2015, knots = c(", kn[2], ",", kn[3], "), Boundary.knots = c(", kn[1], ",", kn[4], ")) +", rhs_cov)),
  data = dat,
  family = binomial()
)
overall_p <- anova(base_fit, spline_fit, test = "Chisq")$`Pr(>Chi)`[2]
nonlinear_p <- anova(linear_fit, spline_fit, test = "Chisq")$`Pr(>Chi)`[2]

template <- dat[1, , drop = FALSE]
template$age_per10 <- mean(dat$age_per10, na.rm = TRUE)
template$bmi_2011 <- mean(dat$bmi_2011, na.rm = TRUE)
for (v in c("current_smoking", "alcohol_last_year", "diabetes", "dyslipidemia")) template[[v]] <- round(mean(dat[[v]], na.rm = TRUE))
for (v in c("sex", "hukou", "education", "marital")) template[[v]] <- names(sort(table(dat[[v]]), decreasing = TRUE))[1]

grid <- data.frame(mean_sbp_2011_2015 = seq(kn[1], kn[4], length.out = 160))
pred_data <- template[rep(1, nrow(grid)), ]
pred_data$mean_sbp_2011_2015 <- grid$mean_sbp_2011_2015
ref_data <- template
ref_data$mean_sbp_2011_2015 <- 120
pred <- predict(spline_fit, newdata = pred_data, type = "link", se.fit = TRUE)
ref <- predict(spline_fit, newdata = ref_data, type = "link", se.fit = TRUE)
plot_data <- data.frame(
  mean_sbp_2011_2015 = grid$mean_sbp_2011_2015,
  or = exp(pred$fit - as.numeric(ref$fit)),
  lower = exp(pred$fit - as.numeric(ref$fit) - 1.96 * pred$se.fit),
  upper = exp(pred$fit - as.numeric(ref$fit) + 1.96 * pred$se.fit)
)
write.csv(plot_data, file.path(table_dir, "figure3_restricted_cubic_spline_plot_data.csv"), row.names = FALSE)

spline_spec <- data.frame(
  item = c("Exposure", "Formula", "Spline type", "Knots", "Reference value", "Extreme values", "Adjustment set", "Overall association P", "Nonlinearity P"),
  value = c(
    "Cumulative mean SBP from 2011 and 2015",
    "mean_sbp_2011_2015 = (SBP_2011 + SBP_2015) / 2",
    "Natural cubic spline/restricted cubic spline with four knots",
    sprintf("5th, 35th, 65th, and 95th percentiles: %.1f, %.1f, %.1f, and %.1f mm Hg", kn[1], kn[2], kn[3], kn[4]),
    "120 mm Hg",
    "No participants were excluded only for the spline model; the curve is plotted over the 5th to 95th percentile range.",
    "Same as the fully adjusted model: age, sex, hukou, education, marital status, BMI, current smoking, alcohol drinking in the past year, diabetes, and dyslipidemia.",
    ifelse(overall_p < 0.001, "<0.001", sprintf("%.3f", overall_p)),
    ifelse(nonlinear_p < 0.001, "<0.001", sprintf("%.3f", nonlinear_p))
  ),
  stringsAsFactors = FALSE
)
write.csv(spline_spec, file.path(table_dir, "supp_table_spline_specification.csv"), row.names = FALSE)
write.csv(data.frame(overall_p = overall_p, nonlinear_p = nonlinear_p, reference_sbp = 120, knot_p5 = kn[1], knot_p35 = kn[2], knot_p65 = kn[3], knot_p95 = kn[4]), file.path(table_dir, "supp_table_mean_sbp_spline_tests_explicit_knots.csv"), row.names = FALSE)

png(file.path(figure_dir, "figure3_cumulative_mean_sbp_restricted_cubic_spline.png"), width = 1800, height = 1300, res = 220)
par(mar = c(4.6, 4.8, 1.2, 1.0))
plot(plot_data$mean_sbp_2011_2015, plot_data$or, type = "n", ylim = range(c(plot_data$lower, plot_data$upper), finite = TRUE), xlab = "Cumulative mean SBP, mm Hg", ylab = "Odds ratio for incident CVD")
polygon(c(plot_data$mean_sbp_2011_2015, rev(plot_data$mean_sbp_2011_2015)), c(plot_data$lower, rev(plot_data$upper)), col = rgb(0.1, 0.35, 0.65, 0.18), border = NA)
lines(plot_data$mean_sbp_2011_2015, plot_data$or, lwd = 2.4, col = "#1F4E79")
abline(h = 1, lty = 2, col = "gray45")
abline(v = 120, lty = 3, col = "gray45")
dev.off()

pdf(file.path(figure_dir, "figure3_cumulative_mean_sbp_restricted_cubic_spline.pdf"), width = 7.0, height = 5.0)
par(mar = c(4.6, 4.8, 1.2, 1.0))
plot(plot_data$mean_sbp_2011_2015, plot_data$or, type = "n", ylim = range(c(plot_data$lower, plot_data$upper), finite = TRUE), xlab = "Cumulative mean SBP, mm Hg", ylab = "Odds ratio for incident CVD")
polygon(c(plot_data$mean_sbp_2011_2015, rev(plot_data$mean_sbp_2011_2015)), c(plot_data$lower, rev(plot_data$upper)), col = rgb(0.1, 0.35, 0.65, 0.18), border = NA)
lines(plot_data$mean_sbp_2011_2015, plot_data$or, lwd = 2.4, col = "#1F4E79")
abline(h = 1, lty = 2, col = "gray45")
abline(v = 120, lty = 3, col = "gray45")
dev.off()

cat("Final SCI reporting package generated.\n")
cat("Main tables: main_table1_baseline_by_9_transition_groups.csv; main_table2_logistic_three_models.csv; main_table3_rr_absolute_risk.csv\n")
cat("Core sensitivity summary: supp_table_core_exposure_sensitivity_summary.csv\n")
cat("Explicit spline: figure3_cumulative_mean_sbp_restricted_cubic_spline.[png/pdf]\n")
