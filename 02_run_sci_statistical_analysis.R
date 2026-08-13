# Analysis: SCI statistical analysis for CHARLS hypertension transition landmark cohort
# Date: 2026-08-11
# Random seed: 42
# R: see runtime
# Key packages: base stats, sandwich, survey, splines

set.seed(42)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_args[grep("--file=", cmd_args)][1]
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)
script_dir <- ifelse(is.na(script_path) || !nzchar(script_path) || script_path == ".", getwd(), dirname(script_path))

data_dir <- file.path(script_dir, "data")
table_dir <- file.path(script_dir, "tables")
figure_dir <- file.path(script_dir, "figures")
log_dir <- file.path(script_dir, "reproducibility_logs")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

need_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)
has_sandwich <- need_pkg("sandwich")
has_survey <- need_pkg("survey")
has_mice <- need_pkg("mice")

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
  persistent_aware_or_treated_uncontrolled = "Persistent aware/treated but uncontrolled",
  gained_control = "Gained control",
  persistent_controlled = "Persistent controlled",
  control_loss = "Control loss",
  other_transition = "Other transitions"
)

format_p <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

or_ci <- function(beta, se) {
  c(exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se))
}

prep_factors <- function(df, exposure = "transition_group") {
  df[[exposure]] <- factor(df[[exposure]], levels = group_order)
  df$sex <- factor(df$sex, levels = c("male", "female"))
  df$hukou <- factor(df$hukou, levels = c("agricultural", "non_agricultural", "other"))
  df$education <- factor(df$education, levels = c("no_formal", "primary_or_below", "middle_school", "high_school_plus"))
  df$marital <- factor(df$marital, levels = c("married_partnered", "not_married_partnered"))
  if ("age_group" %in% names(df)) {
    df$age_group <- factor(df$age_group, levels = c("45-59", ">=60"))
  }
  df
}

model_terms <- list(
  model1_crude = "transition_group",
  model2_demographic = "transition_group + age_per10 + sex + hukou + education + marital",
  model3_full = "transition_group + age_per10 + sex + hukou + education + marital + bmi_2011 + current_smoking + alcohol_last_year + diabetes + dyslipidemia"
)

auc_value <- function(y, pred) {
  ok <- is.finite(y) & is.finite(pred)
  y <- y[ok]
  pred <- pred[ok]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(pred, ties.method = "average")
  (sum(ranks[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

extract_exposure <- function(fit, model_data, outcome, model_name, exposure = "transition_group", vc = NULL, scale_name = "OR") {
  sm <- summary(fit)$coefficients
  if (is.null(vc)) {
    beta <- sm[, "Estimate"]
    se <- sm[, "Std. Error"]
    p <- if ("Pr(>|z|)" %in% colnames(sm)) sm[, "Pr(>|z|)"] else sm[, ncol(sm)]
  } else {
    beta <- coef(fit)
    se <- sqrt(diag(vc))
    p <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  }
  terms <- names(beta)
  keep <- startsWith(terms, exposure)
  kept <- terms[keep]
  if (!length(kept)) return(data.frame())
  groups <- sub(paste0("^", exposure), "", kept)
  est <- exp(beta[kept])
  lo <- exp(beta[kept] - 1.96 * se[kept])
  hi <- exp(beta[kept] + 1.96 * se[kept])
  data.frame(
    outcome = outcome,
    model = model_name,
    scale = scale_name,
    reference = group_labels[["persistent_nonhypertensive"]],
    transition_group = groups,
    transition_group_label = unname(group_labels[groups]),
    n_model = nrow(model_data),
    events_model = sum(model_data[[outcome]] == 1),
    estimate = est,
    ci_lower = lo,
    ci_upper = hi,
    p_value = p[kept],
    estimate_ci = sprintf("%.2f (%.2f-%.2f)", est, lo, hi),
    p_text = vapply(p[kept], format_p, character(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

fit_glm_complete <- function(df, outcome, rhs, family, weights = NULL) {
  form <- as.formula(paste(outcome, "~", rhs))
  vars <- all.vars(form)
  if (!is.null(weights)) vars <- unique(c(vars, weights))
  dat <- df[complete.cases(df[, vars]), vars, drop = FALSE]
  if (is.null(weights)) {
    fit <- glm(form, data = dat, family = family)
  } else {
    fit <- glm(form, data = dat, family = family, weights = dat[[weights]])
  }
  list(fit = fit, data = dat)
}

make_md <- function(df, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(paste0("| ", paste(names(df), collapse = " | "), " |"), con)
  writeLines(paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"), con)
  invisible(apply(df, 1, function(row) writeLines(paste0("| ", paste(row, collapse = " | "), " |"), con)))
}

summ_cont <- function(x) sprintf("%.1f +/- %.1f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
summ_bin <- function(x) sprintf("%d (%.1f)", sum(x == 1, na.rm = TRUE), mean(x == 1, na.rm = TRUE) * 100)

df <- read.csv(
  file.path(data_dir, "sci_landmark_analysis_dataset.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(ID = "character", ID_w1 = "character", cluster_psu = "character")
)
landmark_all <- read.csv(
  file.path(data_dir, "sci_landmark_all_before_2018_followup.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(ID = "character", ID_w1 = "character", cluster_psu = "character")
)
df <- prep_factors(df)
landmark_all <- prep_factors(landmark_all)

# Table 1 by transition group.
cont_vars <- c("age_2011", "bmi_2011", "sbp_2011", "dbp_2011", "sbp_2015", "dbp_2015")
bin_vars <- c("current_smoking", "alcohol_last_year", "diabetes", "dyslipidemia")
cat_vars <- c("sex", "hukou", "education", "marital")
table1_rows <- list()
for (g in group_order) {
  sub <- df[df$transition_group == g, ]
  table1_rows[[length(table1_rows) + 1]] <- data.frame(variable = "N", group = group_labels[g], value = as.character(nrow(sub)))
  for (v in cont_vars) {
    table1_rows[[length(table1_rows) + 1]] <- data.frame(variable = v, group = group_labels[g], value = summ_cont(sub[[v]]))
  }
  for (v in bin_vars) {
    table1_rows[[length(table1_rows) + 1]] <- data.frame(variable = v, group = group_labels[g], value = summ_bin(sub[[v]]))
  }
  for (v in cat_vars) {
    tab <- table(sub[[v]], useNA = "no")
    for (lev in names(tab)) {
      table1_rows[[length(table1_rows) + 1]] <- data.frame(
        variable = paste(v, lev, sep = ": "),
        group = group_labels[g],
        value = sprintf("%d (%.1f)", tab[[lev]], tab[[lev]] / nrow(sub) * 100)
      )
    }
  }
}
table1_long <- do.call(rbind, table1_rows)
write.csv(table1_long, file.path(table_dir, "table1_baseline_by_transition_group_long.csv"), row.names = FALSE)

# Sparse event checks.
sparse <- aggregate(
  cbind(n = rep(1, nrow(df)), cvd = df$incident_cvd_2015_2018, heart = df$incident_heart_2015_2018, stroke = df$incident_stroke_2015_2018) ~ transition_group,
  data = df,
  FUN = sum
)
sparse$transition_group_label <- group_labels[as.character(sparse$transition_group)]
sparse$cvd_sparse_flag <- sparse$cvd < 15
sparse$heart_sparse_flag <- sparse$heart < 10
sparse$stroke_sparse_flag <- sparse$stroke < 10
write.csv(sparse, file.path(table_dir, "supp_table_sparse_event_check.csv"), row.names = FALSE)

# Main logistic models.
logistic_rows <- list()
diagnostic_rows <- list()
for (outcome in c("incident_cvd_2015_2018", "incident_heart_2015_2018", "incident_stroke_2015_2018")) {
  for (mn in names(model_terms)) {
    obj <- fit_glm_complete(df, outcome, model_terms[[mn]], binomial())
    logistic_rows[[length(logistic_rows) + 1]] <- extract_exposure(obj$fit, obj$data, outcome, mn)
    diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
      outcome = outcome,
      model = mn,
      n = nrow(obj$data),
      events = sum(obj$data[[outcome]] == 1),
      predictors_df = length(coef(obj$fit)) - 1,
      epv = sum(obj$data[[outcome]] == 1) / (length(coef(obj$fit)) - 1),
      auc = auc_value(obj$data[[outcome]], fitted(obj$fit))
    )
  }
}
logistic_results <- do.call(rbind, logistic_rows)
logistic_results$p_bh <- ave(logistic_results$p_value, paste(logistic_results$outcome, logistic_results$model), FUN = function(p) p.adjust(p, "BH"))
logistic_results$p_bh_text <- vapply(logistic_results$p_bh, format_p, character(1))
diagnostics <- do.call(rbind, diagnostic_rows)
write.csv(logistic_results, file.path(table_dir, "table2_logistic_or_all_models.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(table_dir, "supp_table_model_diagnostics.csv"), row.names = FALSE)

primary_m3 <- logistic_results[logistic_results$outcome == "incident_cvd_2015_2018" & logistic_results$model == "model3_full", ]
table2 <- primary_m3[, c("transition_group_label", "n_model", "events_model", "estimate_ci", "p_text", "p_bh_text")]
names(table2) <- c("Transition group", "Model N", "Events", "Model 3 OR (95% CI)", "P", "BH-adjusted P")
write.csv(table2, file.path(table_dir, "table2_primary_cvd_model3_or.csv"), row.names = FALSE)
make_md(table2, file.path(table_dir, "table2_primary_cvd_model3_or.md"))

# Modified Poisson RR.
poisson_obj <- fit_glm_complete(df, "incident_cvd_2015_2018", model_terms[["model3_full"]], poisson(link = "log"))
poisson_vc <- if (has_sandwich) sandwich::vcovHC(poisson_obj$fit, type = "HC0") else vcov(poisson_obj$fit)
rr <- extract_exposure(poisson_obj$fit, poisson_obj$data, "incident_cvd_2015_2018", "modified_poisson_model3", vc = poisson_vc, scale_name = "RR")
write.csv(rr, file.path(table_dir, "table3_modified_poisson_rr.csv"), row.names = FALSE)

# Adjusted absolute risks and risk differences from Model 3 logistic.
logit_m3 <- fit_glm_complete(df, "incident_cvd_2015_2018", model_terms[["model3_full"]], binomial())
fit <- logit_m3$fit
model_data <- logit_m3$data
risk_one <- function(g, beta = coef(fit)) {
  nd <- model_data
  nd$transition_group <- factor(g, levels = group_order)
  mm <- model.matrix(delete.response(terms(fit)), data = nd)
  mm <- mm[, names(beta), drop = FALSE]
  mean(plogis(as.numeric(mm %*% beta)))
}
risk_point <- sapply(group_order, risk_one)
if (requireNamespace("MASS", quietly = TRUE)) {
  beta_draws <- MASS::mvrnorm(1000, mu = coef(fit), Sigma = vcov(fit))
  risk_draws <- apply(beta_draws, 1, function(b) sapply(group_order, risk_one, beta = b))
  risk_lo <- apply(risk_draws, 1, quantile, 0.025, na.rm = TRUE)
  risk_hi <- apply(risk_draws, 1, quantile, 0.975, na.rm = TRUE)
  rd_draws <- sweep(risk_draws, 2, risk_draws["persistent_nonhypertensive", ], "-")
  rd_lo <- apply(rd_draws, 1, quantile, 0.025, na.rm = TRUE)
  rd_hi <- apply(rd_draws, 1, quantile, 0.975, na.rm = TRUE)
} else {
  risk_lo <- rep(NA_real_, length(group_order))
  risk_hi <- rep(NA_real_, length(group_order))
  rd_lo <- rep(NA_real_, length(group_order))
  rd_hi <- rep(NA_real_, length(group_order))
}
rd <- risk_point - risk_point[["persistent_nonhypertensive"]]
absrisk <- data.frame(
  transition_group = group_order,
  transition_group_label = unname(group_labels[group_order]),
  adjusted_risk_pct = risk_point * 100,
  adjusted_risk_ci_lower_pct = risk_lo * 100,
  adjusted_risk_ci_upper_pct = risk_hi * 100,
  risk_difference_per_1000 = rd * 1000,
  rd_ci_lower_per_1000 = rd_lo * 1000,
  rd_ci_upper_per_1000 = rd_hi * 1000,
  row.names = NULL
)
write.csv(absrisk, file.path(table_dir, "table3_adjusted_absolute_risk_and_rd.csv"), row.names = FALSE)

# Secondary outcomes table.
secondary <- logistic_results[logistic_results$outcome %in% c("incident_heart_2015_2018", "incident_stroke_2015_2018") & logistic_results$model == "model3_full", ]
write.csv(secondary, file.path(table_dir, "table4_secondary_outcomes_model3_or.csv"), row.names = FALSE)

# Sensitivity: 130/80 mmHg.
df130 <- df
df130$transition_group <- factor(df130$transition_group_130_80, levels = group_order)
sens130 <- fit_glm_complete(df130, "incident_cvd_2015_2018", model_terms[["model3_full"]], binomial())
sens130_out <- extract_exposure(sens130$fit, sens130$data, "incident_cvd_2015_2018", "sensitivity_130_80_model3")
write.csv(sens130_out, file.path(table_dir, "supp_table_sensitivity_130_80_or.csv"), row.names = FALSE)

# Sensitivity: strict known 2016-2018 events only.
strict_df <- df[!is.na(df$strict_cvd_2016_2018), ]
strict_obj <- fit_glm_complete(strict_df, "strict_cvd_2016_2018", model_terms[["model3_full"]], binomial())
strict_out <- extract_exposure(strict_obj$fit, strict_obj$data, "strict_cvd_2016_2018", "sensitivity_strict_2016_2018_model3")
write.csv(strict_out, file.path(table_dir, "supp_table_sensitivity_strict_2016_2018_or.csv"), row.names = FALSE)

# Sensitivity: additional adjustment for baseline SBP and DBP.
bp_rhs <- paste(model_terms[["model3_full"]], "+ sbp_2011 + dbp_2011")
bp_obj <- fit_glm_complete(df, "incident_cvd_2015_2018", bp_rhs, binomial())
bp_out <- extract_exposure(bp_obj$fit, bp_obj$data, "incident_cvd_2015_2018", "sensitivity_baseline_bp_adjusted")
write.csv(bp_out, file.path(table_dir, "supp_table_sensitivity_baseline_bp_adjusted_or.csv"), row.names = FALSE)

# Cluster-robust and survey-weighted sensitivity analyses.
cluster_obj <- fit_glm_complete(df, "incident_cvd_2015_2018", model_terms[["model3_full"]], binomial())
cluster_vars <- unique(c(all.vars(as.formula(paste("incident_cvd_2015_2018 ~", model_terms[["model3_full"]]))), "cluster_psu"))
cluster_dat <- df[complete.cases(df[, cluster_vars]), cluster_vars, drop = FALSE]
cluster_fit <- glm(as.formula(paste("incident_cvd_2015_2018 ~", model_terms[["model3_full"]])), data = cluster_dat, family = binomial())
cluster_vc <- if (has_sandwich) sandwich::vcovCL(cluster_fit, cluster = cluster_dat$cluster_psu, type = "HC0") else vcov(cluster_fit)
cluster_out <- extract_exposure(cluster_fit, cluster_dat, "incident_cvd_2015_2018", "cluster_robust_psu_model3", vc = cluster_vc)
write.csv(cluster_out, file.path(table_dir, "supp_table_cluster_robust_psu_or.csv"), row.names = FALSE)

weighted_df <- df
weighted_df$w_baseline_biomarker_scaled <- weighted_df$baseline_biomarker_weight / mean(weighted_df$baseline_biomarker_weight, na.rm = TRUE)
q <- quantile(weighted_df$w_baseline_biomarker_scaled, c(0.01, 0.99), na.rm = TRUE)
weighted_df$w_baseline_biomarker_trunc <- pmin(pmax(weighted_df$w_baseline_biomarker_scaled, q[1]), q[2])
weighted_obj <- fit_glm_complete(weighted_df, "incident_cvd_2015_2018", model_terms[["model3_full"]], quasibinomial(), weights = "w_baseline_biomarker_trunc")
weighted_out <- extract_exposure(weighted_obj$fit, weighted_obj$data, "incident_cvd_2015_2018", "baseline_biomarker_weighted_model3")
write.csv(weighted_out, file.path(table_dir, "supp_table_baseline_weighted_or.csv"), row.names = FALSE)

if (has_survey) {
  options(survey.lonely.psu = "adjust")
  svy_vars <- unique(c(all.vars(as.formula(paste("incident_cvd_2015_2018 ~", model_terms[["model3_full"]]))), "cluster_psu", "strata_proxy", "w_baseline_biomarker_trunc"))
  svy_dat <- weighted_df[complete.cases(weighted_df[, svy_vars]), svy_vars, drop = FALSE]
  des <- survey::svydesign(id = ~cluster_psu, strata = ~strata_proxy, weights = ~w_baseline_biomarker_trunc, data = svy_dat, nest = TRUE)
  svy_fit <- survey::svyglm(as.formula(paste("incident_cvd_2015_2018 ~", model_terms[["model3_full"]])), design = des, family = quasibinomial())
  svy_out <- extract_exposure(svy_fit, svy_dat, "incident_cvd_2015_2018", "survey_design_baseline_weight_model3")
  write.csv(svy_out, file.path(table_dir, "supp_table_survey_design_weighted_or.csv"), row.names = FALSE)
}

# Loss to follow-up comparison and IPCW sensitivity.
landmark_all$observed_binary <- as.numeric(as.character(landmark_all$observed_2018_outcome) %in% c("TRUE", "True", "true", "1"))
loss_vars <- c("age_2011", "sex", "education", "hukou", "bmi_2011", "sbp_2011", "dbp_2011", "diabetes", "dyslipidemia", "transition_group")
loss_binary_vars <- c("diabetes", "dyslipidemia")
loss_rows <- list()
for (v in loss_vars) {
  if (v %in% loss_binary_vars) {
    loss_rows[[length(loss_rows) + 1]] <- data.frame(
      variable = v,
      observed_2018 = summ_bin(landmark_all[landmark_all$observed_binary == 1, v]),
      lost_2018 = summ_bin(landmark_all[landmark_all$observed_binary == 0, v])
    )
  } else if (is.numeric(landmark_all[[v]])) {
    loss_rows[[length(loss_rows) + 1]] <- data.frame(variable = v, observed_2018 = summ_cont(landmark_all[landmark_all$observed_binary == 1, v]), lost_2018 = summ_cont(landmark_all[landmark_all$observed_binary == 0, v]))
  } else {
    tab <- table(landmark_all[[v]], landmark_all$observed_binary)
    for (lev in rownames(tab)) {
      loss_rows[[length(loss_rows) + 1]] <- data.frame(
        variable = paste(v, lev, sep = ": "),
        observed_2018 = if ("1" %in% colnames(tab)) sprintf("%d", tab[lev, "1"]) else "0",
        lost_2018 = if ("0" %in% colnames(tab)) sprintf("%d", tab[lev, "0"]) else "0"
      )
    }
  }
}
loss_table <- do.call(rbind, loss_rows)
write.csv(loss_table, file.path(table_dir, "supp_table_loss_to_followup_comparison.csv"), row.names = FALSE)

ipcw_rhs <- "transition_group + age_per10 + sex + hukou + education + marital + bmi_2011 + current_smoking + alcohol_last_year + diabetes + dyslipidemia + sbp_2011 + dbp_2011 + sbp_2015 + dbp_2015"
obs_obj <- fit_glm_complete(landmark_all, "observed_binary", ipcw_rhs, binomial())
obs_data <- obs_obj$data
obs_data$ID <- landmark_all[rownames(obs_data), "ID"]
obs_data$incident_cvd_2015_2018 <- landmark_all[rownames(obs_data), "incident_cvd_2015_2018"]
obs_data$pred_observed <- pmin(pmax(fitted(obs_obj$fit), 0.02), 0.98)
stab <- mean(obs_data$observed_binary == 1)
obs_data$ipcw <- ifelse(obs_data$observed_binary == 1, stab / obs_data$pred_observed, NA_real_)
iq <- quantile(obs_data$ipcw, c(0.01, 0.99), na.rm = TRUE)
obs_data$ipcw_trunc <- pmin(pmax(obs_data$ipcw, iq[1]), iq[2])
ipcw_weights <- obs_data[obs_data$observed_binary == 1, c("ID", "ipcw_trunc")]
ipcw_analysis <- merge(df, ipcw_weights, by = "ID")
ipcw_analysis <- prep_factors(ipcw_analysis)
ipcw_obj <- fit_glm_complete(ipcw_analysis, "incident_cvd_2015_2018", model_terms[["model3_full"]], quasibinomial(), weights = "ipcw_trunc")
ipcw_fit <- ipcw_obj$fit
ipcw_analysis <- ipcw_obj$data
ipcw_out <- extract_exposure(ipcw_fit, ipcw_analysis, "incident_cvd_2015_2018", "ipcw_model3")
write.csv(ipcw_out, file.path(table_dir, "supp_table_ipcw_or.csv"), row.names = FALSE)
write.csv(data.frame(n_landmark = nrow(landmark_all), n_observed = sum(landmark_all$observed_binary == 1), n_lost = sum(landmark_all$observed_binary == 0), ipcw_p1 = iq[1], ipcw_p99 = iq[2]), file.path(table_dir, "supp_table_ipcw_weight_summary.csv"), row.names = FALSE)

# Multiple imputation sensitivity. Prefer mice if installed; fallback is custom
# stochastic chained imputation because covariate missingness is low.
impute_vars <- c("age_per10", "sex", "hukou", "education", "marital", "bmi_2011", "current_smoking", "alcohol_last_year", "diabetes", "dyslipidemia", "transition_group", "incident_cvd_2015_2018", "sbp_2011", "dbp_2011", "sbp_2015", "dbp_2015", "incident_heart_2015_2018", "incident_stroke_2015_2018")
mi_base <- df[, impute_vars]
mi_base <- prep_factors(mi_base)
mi_results <- list()
if (has_mice) {
  imp <- mice::mice(mi_base, m = 20, maxit = 20, seed = 42, printFlag = FALSE)
  for (i in 1:20) {
    comp <- mice::complete(imp, i)
    comp <- prep_factors(comp)
    obj <- fit_glm_complete(comp, "incident_cvd_2015_2018", model_terms[["model3_full"]], binomial())
    b <- coef(obj$fit)
    vc <- vcov(obj$fit)
    terms <- names(b)[startsWith(names(b), "transition_group")]
    mi_results[[i]] <- data.frame(imputation = i, term = terms, beta = b[terms], var = diag(vc)[terms], row.names = NULL)
  }
  mi_method <- "mice_package_m20_maxit20"
} else {
  for (i in 1:20) {
    comp <- mi_base
    for (v in c("sex", "hukou", "education", "marital")) {
      miss <- is.na(comp[[v]])
      if (any(miss)) comp[[v]][miss] <- sample(comp[[v]][!is.na(comp[[v]])], sum(miss), replace = TRUE)
    }
    for (v in c("current_smoking", "alcohol_last_year", "diabetes", "dyslipidemia")) {
      miss <- is.na(comp[[v]])
      if (any(miss)) {
        p <- mean(comp[[v]] == 1, na.rm = TRUE)
        comp[[v]][miss] <- rbinom(sum(miss), 1, p)
      }
    }
    miss_bmi <- is.na(comp$bmi_2011)
    if (any(miss_bmi)) {
      fit_bmi <- lm(bmi_2011 ~ age_per10 + sex + hukou + education + marital + transition_group + sbp_2011 + dbp_2011 + sbp_2015 + dbp_2015, data = comp)
      pred <- predict(fit_bmi, newdata = comp[miss_bmi, ])
      resid_pool <- residuals(fit_bmi)
      comp$bmi_2011[miss_bmi] <- pred + sample(resid_pool, sum(miss_bmi), replace = TRUE)
    }
    comp <- prep_factors(comp)
    obj <- fit_glm_complete(comp, "incident_cvd_2015_2018", model_terms[["model3_full"]], binomial())
    b <- coef(obj$fit)
    vc <- vcov(obj$fit)
    terms <- names(b)[startsWith(names(b), "transition_group")]
    mi_results[[i]] <- data.frame(imputation = i, term = terms, beta = b[terms], var = diag(vc)[terms], row.names = NULL)
  }
  mi_method <- "fallback_stochastic_chained_imputation_m20_due_to_missing_mice_package"
}
mi_long <- do.call(rbind, mi_results)
mi_pool <- do.call(rbind, lapply(split(mi_long, mi_long$term), function(x) {
  m <- nrow(x)
  qbar <- mean(x$beta)
  ubar <- mean(x$var)
  bvar <- var(x$beta)
  tvar <- ubar + (1 + 1 / m) * bvar
  se <- sqrt(tvar)
  grp <- sub("^transition_group", "", unique(x$term))
  data.frame(
    method = mi_method,
    transition_group = grp,
    transition_group_label = group_labels[grp],
    or = exp(qbar),
    ci_lower = exp(qbar - 1.96 * se),
    ci_upper = exp(qbar + 1.96 * se),
    p_value = 2 * pnorm(abs(qbar / se), lower.tail = FALSE),
    estimate_ci = sprintf("%.2f (%.2f-%.2f)", exp(qbar), exp(qbar - 1.96 * se), exp(qbar + 1.96 * se)),
    p_text = format_p(2 * pnorm(abs(qbar / se), lower.tail = FALSE)),
    row.names = NULL
  )
}))
write.csv(mi_pool, file.path(table_dir, "supp_table_multiple_imputation_model3_or.csv"), row.names = FALSE)

# Subgroups by sex and age.
subgroup_rows <- list()
for (sg in c("sex", "age_group")) {
  base_rhs <- model_terms[["model3_full"]]
  no_sg_rhs <- gsub(paste0(" \\+ ", sg), "", base_rhs, fixed = TRUE)
  no_sg_rhs <- gsub(paste0(sg, " \\+ "), "", no_sg_rhs)
  int_full <- as.formula(paste("incident_cvd_2015_2018 ~", no_sg_rhs, "+ transition_group:", sg))
  int_reduced <- as.formula(paste("incident_cvd_2015_2018 ~", base_rhs))
  int_vars <- unique(c(all.vars(int_full), all.vars(int_reduced)))
  int_dat <- df[complete.cases(df[, int_vars]), int_vars, drop = FALSE]
  p_int <- tryCatch(anova(glm(int_reduced, data = int_dat, family = binomial()), glm(int_full, data = int_dat, family = binomial()), test = "LRT")$`Pr(>Chi)`[2], error = function(e) NA_real_)
  for (lev in levels(df[[sg]])) {
    sub <- df[df[[sg]] == lev, ]
    rhs <- base_rhs
    rhs <- gsub(paste0(" + ", sg), "", rhs, fixed = TRUE)
    rhs <- gsub(paste0(sg, " + "), "", rhs, fixed = TRUE)
    obj <- fit_glm_complete(sub, "incident_cvd_2015_2018", rhs, binomial())
    out <- extract_exposure(obj$fit, obj$data, "incident_cvd_2015_2018", paste0("subgroup_", sg, "_", lev))
    out$subgroup_variable <- sg
    out$subgroup_level <- lev
    out$interaction_p <- p_int
    subgroup_rows[[length(subgroup_rows) + 1]] <- out
  }
}
subgroups <- do.call(rbind, subgroup_rows)
write.csv(subgroups, file.path(table_dir, "supp_table_subgroup_sex_age_or.csv"), row.names = FALSE)

# Dose-response support: cumulative average SBP using natural cubic spline.
rcs_rhs <- "splines::ns(mean_sbp_2011_2015, df = 3) + age_per10 + sex + hukou + education + marital + bmi_2011 + current_smoking + alcohol_last_year + diabetes + dyslipidemia"
rcs_obj <- fit_glm_complete(df, "incident_cvd_2015_2018", rcs_rhs, binomial())
linear_obj <- fit_glm_complete(df, "incident_cvd_2015_2018", gsub("splines::ns\\(mean_sbp_2011_2015, df = 3\\)", "mean_sbp_2011_2015", rcs_rhs), binomial())
null_obj <- fit_glm_complete(df, "incident_cvd_2015_2018", "age_per10 + sex + hukou + education + marital + bmi_2011 + current_smoking + alcohol_last_year + diabetes + dyslipidemia", binomial())
nonlinear_p <- anova(linear_obj$fit, rcs_obj$fit, test = "LRT")$`Pr(>Chi)`[2]
overall_p <- anova(null_obj$fit, rcs_obj$fit, test = "LRT")$`Pr(>Chi)`[2]
sbp_grid <- seq(quantile(rcs_obj$data$mean_sbp_2011_2015, 0.05), quantile(rcs_obj$data$mean_sbp_2011_2015, 0.95), length.out = 100)
ref <- 120
nd <- rcs_obj$data[rep(1, length(sbp_grid)), ]
for (v in names(nd)) {
  if (is.numeric(nd[[v]])) nd[[v]] <- mean(rcs_obj$data[[v]], na.rm = TRUE)
}
for (v in c("sex", "hukou", "education", "marital")) {
  nd[[v]] <- factor(names(sort(table(rcs_obj$data[[v]]), decreasing = TRUE))[1], levels = levels(rcs_obj$data[[v]]))
}
nd$mean_sbp_2011_2015 <- sbp_grid
pred <- predict(rcs_obj$fit, newdata = nd, type = "link", se.fit = TRUE)
nd_ref <- nd[1, ]
nd_ref$mean_sbp_2011_2015 <- ref
ref_lp <- predict(rcs_obj$fit, newdata = nd_ref, type = "link")
rcs_plot <- data.frame(mean_sbp = sbp_grid, or = exp(pred$fit - ref_lp), lo = exp(pred$fit - 1.96 * pred$se.fit - ref_lp), hi = exp(pred$fit + 1.96 * pred$se.fit - ref_lp))
write.csv(data.frame(overall_p = overall_p, nonlinear_p = nonlinear_p, reference_sbp = ref), file.path(table_dir, "supp_table_mean_sbp_spline_tests.csv"), row.names = FALSE)
write.csv(rcs_plot, file.path(table_dir, "figure3_mean_sbp_spline_plot_data.csv"), row.names = FALSE)

# Figures.
forest <- primary_m3
forest$transition_group_label <- factor(forest$transition_group_label, levels = rev(forest$transition_group_label))
draw_forest <- function(path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") pdf(path, width = 8, height = 5)
  if (device == "png") png(path, width = 2400, height = 1500, res = 300)
  op <- par(mar = c(4, 12, 2, 2))
  y <- seq_len(nrow(forest))
  plot(forest$estimate, y, xlim = c(0.4, max(forest$ci_upper, na.rm = TRUE) * 1.15), yaxt = "n", xlab = "OR (95% CI)", ylab = "", pch = 19, log = "x")
  segments(forest$ci_lower, y, forest$ci_upper, y)
  abline(v = 1, lty = 2, col = "gray40")
  axis(2, at = y, labels = forest$transition_group_label, las = 1, cex.axis = 0.75)
  title("Changes in hypertension management status and incident CVD")
  par(op)
  dev.off()
}
draw_forest(file.path(figure_dir, "figure2_primary_model3_forest.pdf"), "pdf")
draw_forest(file.path(figure_dir, "figure2_primary_model3_forest.png"), "png")

draw_spline <- function(path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") pdf(path, width = 6, height = 4)
  if (device == "png") png(path, width = 1800, height = 1200, res = 300)
  plot(rcs_plot$mean_sbp, rcs_plot$or, type = "l", lwd = 2, xlab = "Cumulative average SBP, mmHg", ylab = "OR for incident CVD")
  lines(rcs_plot$mean_sbp, rcs_plot$lo, lty = 2, col = "gray40")
  lines(rcs_plot$mean_sbp, rcs_plot$hi, lty = 2, col = "gray40")
  abline(h = 1, lty = 3)
  dev.off()
}
draw_spline(file.path(figure_dir, "figure3_mean_sbp_spline.pdf"), "pdf")
draw_spline(file.path(figure_dir, "figure3_mean_sbp_spline.png"), "png")

# Manifest and manuscript-ready summary.
summary_lines <- c(
  "# SCI Analysis Outputs",
  paste0("Generated: ", Sys.Date()),
  "Study type: prospective landmark cohort analysis",
  "",
  "## Main data",
  "- `data/sci_landmark_analysis_dataset.csv` -- 2015 landmark analysis dataset with 2018 outcomes",
  "- `data/sci_landmark_all_before_2018_followup.csv` -- 2015 landmark eligible sample before 2018 follow-up restriction",
  "",
  "## Main tables",
  "- `tables/supp_table_5x5_transition_matrix_details.csv` -- full 5 x 5 transition matrix with events, event rates, and mean SBP/DBP",
  "- `tables/table1_baseline_by_transition_group_long.csv` -- baseline characteristics by SCI transition group",
  "- `tables/table2_primary_cvd_model3_or.csv` -- primary fully adjusted ORs",
  "- `tables/table3_modified_poisson_rr.csv` -- modified Poisson RR sensitivity",
  "- `tables/table3_adjusted_absolute_risk_and_rd.csv` -- marginal adjusted risk and risk difference",
  "- `tables/table4_secondary_outcomes_model3_or.csv` -- secondary heart disease and stroke outcomes",
  "",
  "## Sensitivity and supplementary analyses",
  "- `tables/supp_table_missingness.csv` -- variable missingness",
  "- `tables/supp_table_event_year_quality.csv` -- CVD diagnosis year quality",
  "- `tables/supp_table_sensitivity_130_80_or.csv` -- 130/80 mmHg threshold sensitivity",
  "- `tables/supp_table_sensitivity_strict_2016_2018_or.csv` -- strict known 2016-2018 event sensitivity",
  "- `tables/supp_table_multiple_imputation_model3_or.csv` -- multiple imputation sensitivity",
  "- `tables/supp_table_cluster_robust_psu_or.csv` -- PSU cluster-robust sensitivity",
  "- `tables/supp_table_baseline_weighted_or.csv` -- baseline biomarker weighted sensitivity",
  "- `tables/supp_table_ipcw_or.csv` -- inverse probability of censoring weighted sensitivity",
  "- `tables/supp_table_subgroup_sex_age_or.csv` -- sex and age subgroup analyses",
  "- `tables/supp_table_mean_sbp_spline_tests.csv` -- cumulative average SBP spline tests",
  "",
  "## Figures",
  "- `figures/figure2_primary_model3_forest.pdf` / `.png` -- primary forest plot",
  "- `figures/figure3_mean_sbp_spline.pdf` / `.png` -- cumulative average SBP spline plot"
)
writeLines(summary_lines, file.path(script_dir, "_analysis_outputs.md"))

res <- merge(primary_m3[, c("transition_group", "transition_group_label", "estimate_ci", "p_text")], absrisk[, c("transition_group", "adjusted_risk_pct", "risk_difference_per_1000")], by = "transition_group")
res <- res[match(group_order[-1], res$transition_group), ]
methods_text <- c(
  "Statistical methods snippet:",
  "We conducted a prospective landmark cohort analysis using 2015 as the landmark time point. Hypertension management status transitions from 2011 to 2015 were modeled as the exposure, and incident CVD from 2015 to 2018 was the primary outcome.",
  "Logistic regression was used as the primary model because follow-up was assessed over a fixed interval. Model 1 was unadjusted, Model 2 adjusted for age, sex, hukou, education and marital status, and Model 3 further adjusted for BMI, current smoking, alcohol drinking in the past year, diabetes and dyslipidemia.",
  "Sensitivity analyses included modified Poisson regression with robust standard errors, marginal standardized absolute risks, multiple imputation for missing covariates, 130/80 mmHg hypertension threshold, strict 2016-2018 event definition, baseline biomarker weighting, PSU cluster-robust standard errors, IPCW for loss to follow-up, and sex/age subgroup analyses."
)
results_text <- c(
  "Results snippet:",
  paste0("The 2015 landmark analysis included ", nrow(df), " participants without CVD at the landmark, among whom ", sum(df$incident_cvd_2015_2018 == 1), " developed incident CVD by 2018."),
  paste0("In the fully adjusted model, the largest association was observed for ", primary_m3$transition_group_label[which.max(primary_m3$estimate)], " (OR=", primary_m3$estimate_ci[which.max(primary_m3$estimate)], "; P=", primary_m3$p_text[which.max(primary_m3$estimate)], ")."),
  paste0("Adjusted absolute risk differences ranged from ", sprintf("%.1f", min(absrisk$risk_difference_per_1000[-1], na.rm = TRUE)), " to ", sprintf("%.1f", max(absrisk$risk_difference_per_1000[-1], na.rm = TRUE)), " additional CVD events per 1,000 persons over approximately 3 years compared with persistent non-hypertension.")
)
writeLines(c(methods_text, "", results_text), file.path(script_dir, "sci_methods_and_results_snippets.txt"))

cat("SCI statistical analysis complete\n")
cat("Primary Model 3 ORs:\n")
print(table2, row.names = FALSE)
cat("\nModified Poisson RR:\n")
print(rr[, c("transition_group_label", "estimate_ci", "p_text")], row.names = FALSE)
cat("\nAdjusted absolute risks and risk differences per 1000:\n")
print(absrisk[, c("transition_group_label", "adjusted_risk_pct", "risk_difference_per_1000")], row.names = FALSE)
cat("\nMultiple imputation method:", mi_method, "\n")
