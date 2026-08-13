# BMC-oriented RR-primary reanalysis for the CHARLS hypertension transition study.
# Outputs are written under SCI submission tables and do not modify Chinese analyses.

set.seed(42)

cmd_args <- commandArgs(FALSE)
file_arg <- cmd_args[grep("--file=", cmd_args)][1]
script_dir <- if (length(file_arg) && !is.na(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)) else getwd()
if (!dir.exists(file.path(script_dir, "tables"))) script_dir <- getwd()

data_dir <- file.path(script_dir, "data")
table_dir <- file.path(script_dir, "tables")
figure_dir <- file.path(script_dir, "figures")
raw_root <- Sys.getenv("CHARLS_RAW_ROOT", unset = file.path(script_dir, "raw", "_extracted"))
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

need_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)
if (!need_pkg("sandwich")) stop("The sandwich package is required for cluster-robust modified Poisson models.")

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
  persistent_nonhypertensive = "Persistent non-hypertension",
  incident_hypertension = "Incident hypertension",
  persistent_unaware_or_untreated = "Persistent unawareness/no treatment",
  new_awareness_or_treatment_uncontrolled = "New awareness/treatment without control",
  # Label must match the manuscript text and Tables 1-4; do not reintroduce
  # "but uncontrolled", which is the five-state name, not the transition group.
  persistent_aware_or_treated_uncontrolled = "Persistent awareness or treatment without control",
  gained_control = "Gained control",
  persistent_controlled = "Persistent control",
  control_loss = "Loss of control",
  other_transition = "Other transitions"
)

format_p <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

yes <- function(x, code = 1) !is.na(x) & x == code

mean_valid <- function(dat, vars, lower, upper) {
  vals <- as.data.frame(lapply(dat[, vars], function(x) suppressWarnings(as.numeric(x))))
  vals[vals < lower | vals > upper] <- NA
  apply(vals, 1, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
}

convert_wave1_to_wave2_id <- function(id) {
  id <- as.character(id)
  paste0(substr(id, 1, 9), "0", substr(id, 10, nchar(id)))
}

pad_charls_id <- function(x) {
  x <- trimws(as.character(x))
  missing <- is.na(x) | x == "NA"
  x[!missing] <- ifelse(nchar(x[!missing]) < 12, paste0(strrep("0", 12 - nchar(x[!missing])), x[!missing]), x[!missing])
  x[missing] <- NA_character_
  x
}

classify_htn_state <- function(sbp, dbp, diagnosis, awareness, chinese_med, western_med, western_only = FALSE) {
  bp_valid <- !is.na(sbp) & !is.na(dbp)
  high <- bp_valid & (sbp >= 140 | dbp >= 90)
  treated <- if (western_only) yes(western_med, 2) else (yes(chinese_med, 1) | yes(western_med, 2))
  aware <- yes(diagnosis, 1) | yes(awareness, 1) | treated
  # Branch order must match add_hypertension_state() in
  # 01_build_sci_landmark_dataset.py: "normal" requires the participant to be
  # unaware AND untreated. Testing (!high & !treated) first would move aware but
  # untreated participants with normal blood pressure into the non-hypertensive
  # reference group.
  ifelse(
    !bp_valid, NA_character_,
    ifelse(!aware & !high, "normal",
      ifelse(!aware & high, "unaware_htn",
        ifelse(!treated, "aware_untreated",
          ifelse(high, "treated_uncontrolled", "treated_controlled")
        )
      )
    )
  )
}

reconstructed_2015_status_r <- function(dat, disease_index) {
  previous_yes <- yes(dat[[paste0("zda007_", disease_index, "_")]], 1)
  if (disease_index == 1) previous_yes <- previous_yes | yes(dat$zda008_1_, 1)
  confirmation <- dat[[paste0("da007_w2_1_", disease_index, "_")]]
  new_diagnosis <- dat[[paste0("da007_w2_2_", disease_index, "_")]]
  raw_status <- dat[[paste0("da007_", disease_index, "_")]]
  out <- rep(NA_integer_, nrow(dat))
  has_confirmation <- !is.na(confirmation)
  out[has_confirmation & previous_yes] <- as.integer(yes(confirmation[has_confirmation & previous_yes], 1) | (yes(confirmation[has_confirmation & previous_yes], 2) & yes(new_diagnosis[has_confirmation & previous_yes], 1)))
  out[has_confirmation & !previous_yes] <- as.integer(yes(confirmation[has_confirmation & !previous_yes], 2) | (yes(confirmation[has_confirmation & !previous_yes], 1) & yes(new_diagnosis[has_confirmation & !previous_yes], 1)))
  no_confirmation <- is.na(out) & !is.na(raw_status)
  out[no_confirmation] <- as.integer(yes(raw_status[no_confirmation], 1))
  out
}

map_transition_group <- function(state_2011, state_2015) {
  key <- paste(state_2011, state_2015, sep = "__")
  out <- rep(NA_character_, length(key))
  out[key == "normal__normal"] <- "persistent_nonhypertensive"
  out[state_2011 == "normal" & state_2015 != "normal" & !is.na(state_2015)] <- "incident_hypertension"
  out[key %in% c("unaware_htn__unaware_htn", "unaware_htn__aware_untreated", "aware_untreated__unaware_htn", "aware_untreated__aware_untreated")] <- "persistent_unaware_or_untreated"
  out[key == "unaware_htn__treated_uncontrolled"] <- "new_awareness_or_treatment_uncontrolled"
  out[key %in% c("aware_untreated__treated_uncontrolled", "treated_uncontrolled__aware_untreated", "treated_uncontrolled__treated_uncontrolled")] <- "persistent_aware_or_treated_uncontrolled"
  out[key %in% c("unaware_htn__treated_controlled", "aware_untreated__treated_controlled", "treated_uncontrolled__treated_controlled")] <- "gained_control"
  out[key == "treated_controlled__treated_controlled"] <- "persistent_controlled"
  out[key %in% c("treated_controlled__unaware_htn", "treated_controlled__aware_untreated", "treated_controlled__treated_uncontrolled")] <- "control_loss"
  out[is.na(out) & !is.na(state_2011) & !is.na(state_2015)] <- "other_transition"
  out
}

prep_factors <- function(df, exposure = "transition_group") {
  df[[exposure]] <- factor(df[[exposure]], levels = group_order)
  df$sex <- factor(df$sex, levels = c("male", "female"))
  df$hukou <- factor(df$hukou, levels = c("agricultural", "non_agricultural", "other"))
  df$education <- factor(df$education, levels = c("no_formal", "primary_or_below", "middle_school", "high_school_plus"))
  df$marital <- factor(df$marital, levels = c("married_partnered", "not_married_partnered"))
  if ("age_group" %in% names(df)) df$age_group <- factor(df$age_group, levels = c("45-59", ">=60"))
  df
}

fit_complete <- function(df, outcome, rhs, family = poisson(link = "log"), weights = NULL, extra_vars = NULL) {
  form <- as.formula(paste(outcome, "~", rhs))
  vars <- all.vars(form)
  if (!is.null(weights)) vars <- unique(c(vars, weights))
  if (!is.null(extra_vars)) vars <- unique(c(vars, extra_vars))
  dat <- df[complete.cases(df[, vars]), vars, drop = FALSE]
  if (is.null(weights)) {
    fit <- glm(form, data = dat, family = family)
  } else {
    fit <- glm(form, data = dat, family = family, weights = dat[[weights]])
  }
  list(fit = fit, data = dat)
}

vcov_cluster <- function(fit, dat) {
  sandwich::vcovCL(fit, cluster = dat$cluster_psu, type = "HC0")
}

extract_transition <- function(fit, dat, outcome, model_name, vc, scale = "RR") {
  beta <- coef(fit)
  se <- sqrt(diag(vc))
  p <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  terms <- names(beta)
  keep <- startsWith(terms, "transition_group")
  kept <- terms[keep]
  groups <- sub("^transition_group", "", kept)
  est <- exp(beta[kept])
  lo <- exp(beta[kept] - 1.96 * se[kept])
  hi <- exp(beta[kept] + 1.96 * se[kept])
  data.frame(
    outcome = outcome,
    model = model_name,
    scale = scale,
    reference = group_labels[["persistent_nonhypertensive"]],
    transition_group = groups,
    transition_group_label = unname(group_labels[groups]),
    n_model = nrow(dat),
    events_model = sum(dat[[outcome]] == 1),
    estimate = est,
    ci_lower = lo,
    ci_upper = hi,
    p_value = p[kept],
    estimate_ci = sprintf("%.2f (%.2f-%.2f)", est, lo, hi),
    p_text = vapply(p[kept], format_p, character(1)),
    row.names = NULL
  )
}

global_wald <- function(fit, vc) {
  beta <- coef(fit)
  keep <- startsWith(names(beta), "transition_group")
  b <- beta[keep]
  v <- vc[keep, keep, drop = FALSE]
  stat <- as.numeric(t(b) %*% solve(v, b))
  df <- length(b)
  p <- pchisq(stat, df = df, lower.tail = FALSE)
  data.frame(test = "Global Wald test for transition_group", statistic = stat, df = df, p_value = p, p_text = format_p(p))
}

contrast_rr <- function(fit, vc, group_a, group_b, label) {
  beta <- coef(fit)
  coef_name <- function(g) paste0("transition_group", g)
  L <- rep(0, length(beta))
  names(L) <- names(beta)
  if (group_a != "persistent_nonhypertensive") L[coef_name(group_a)] <- L[coef_name(group_a)] + 1
  if (group_b != "persistent_nonhypertensive") L[coef_name(group_b)] <- L[coef_name(group_b)] - 1
  log_rr <- sum(L * beta)
  se <- sqrt(as.numeric(t(L) %*% vc %*% L))
  p <- 2 * pnorm(abs(log_rr / se), lower.tail = FALSE)
  data.frame(
    contrast = label,
    numerator = unname(group_labels[group_a]),
    denominator = unname(group_labels[group_b]),
    rr = exp(log_rr),
    ci_lower = exp(log_rr - 1.96 * se),
    ci_upper = exp(log_rr + 1.96 * se),
    p_value = p,
    rr_95ci = sprintf("%.2f (%.2f-%.2f)", exp(log_rr), exp(log_rr - 1.96 * se), exp(log_rr + 1.96 * se)),
    p_text = format_p(p),
    row.names = NULL
  )
}

standardized_risk_draws <- function(df, rhs, outcome = "incident_cvd_2015_2018", ndraw = 1000) {
  obj <- fit_complete(df, outcome, rhs, binomial())
  fit <- obj$fit
  dat <- obj$data
  risk_one <- function(g, beta = coef(fit)) {
    nd <- dat
    nd$transition_group <- factor(g, levels = group_order)
    mm <- model.matrix(delete.response(terms(fit)), data = nd)
    mm <- mm[, names(beta), drop = FALSE]
    mean(plogis(as.numeric(mm %*% beta)))
  }
  risk_point <- sapply(group_order, risk_one)
  beta_draws <- MASS::mvrnorm(ndraw, mu = coef(fit), Sigma = vcov(fit))
  risk_draws <- apply(beta_draws, 1, function(b) sapply(group_order, risk_one, beta = b))
  list(point = risk_point, draws = risk_draws)
}

gcomp_point <- function(fit, dat) {
  beta <- coef(fit)
  risk_one <- function(g) {
    nd <- dat
    nd$transition_group <- factor(g, levels = group_order)
    mm <- model.matrix(delete.response(terms(fit)), data = nd)
    mm <- mm[, names(beta), drop = FALSE]
    mean(plogis(as.numeric(mm %*% beta)))
  }
  risks <- sapply(group_order, risk_one)
  ref_risk <- risks[["persistent_nonhypertensive"]]
  data.frame(
    transition_group = group_order,
    adjusted_risk = as.numeric(risks),
    marginal_rr = as.numeric(risks / ref_risk),
    risk_difference = as.numeric(risks - ref_risk),
    row.names = NULL
  )
}

cluster_bootstrap_gcomp <- function(df, rhs, outcome = "incident_cvd_2015_2018", n_boot = 1000) {
  obj <- fit_complete(df, outcome, rhs, binomial(), extra_vars = "cluster_psu")
  dat <- obj$data
  point <- gcomp_point(obj$fit, dat)
  clusters <- unique(dat$cluster_psu)
  cluster_rows <- split(seq_len(nrow(dat)), dat$cluster_psu)
  draws <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
    boot_idx <- unlist(cluster_rows[sampled_clusters], use.names = FALSE)
    boot_dat <- dat[boot_idx, , drop = FALSE]
    boot_fit <- tryCatch(glm(as.formula(paste(outcome, "~", rhs)), data = boot_dat, family = binomial(), model = FALSE, y = FALSE), error = function(e) NULL)
    if (is.null(boot_fit) || any(!is.finite(coef(boot_fit)))) next
    boot_point <- tryCatch(gcomp_point(boot_fit, boot_dat), error = function(e) NULL)
    if (is.null(boot_point)) next
    boot_point$bootstrap <- b
    draws[[b]] <- boot_point
  }
  draws <- do.call(rbind, draws)
  if (is.null(draws) || nrow(draws) == 0) stop("No valid bootstrap replicates for g-computation.")
  summarize_metric <- function(metric, prefix) {
    pieces <- lapply(split(draws[[metric]], draws$transition_group), function(x) {
      qs <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
      c(lower = unname(qs[1]), upper = unname(qs[2]))
    })
    out <- data.frame(transition_group = names(pieces), do.call(rbind, pieces), row.names = NULL)
    names(out)[2:3] <- paste0(prefix, c("_ci_lower", "_ci_upper"))
    out
  }
  risk_ci <- summarize_metric("adjusted_risk", "risk")
  rr_ci <- summarize_metric("marginal_rr", "rr")
  rd_ci <- summarize_metric("risk_difference", "rd")
  out <- Reduce(function(x, y) merge(x, y, by = "transition_group", all = TRUE), list(point, risk_ci, rr_ci, rd_ci))
  out$valid_bootstrap_replicates <- length(unique(draws$bootstrap))
  out <- out[, c(
    "transition_group", "adjusted_risk", "risk_ci_lower", "risk_ci_upper",
    "marginal_rr", "rr_ci_lower", "rr_ci_upper",
    "risk_difference", "rd_ci_lower", "rd_ci_upper", "valid_bootstrap_replicates"
  )]
  attr(out, "draws") <- draws
  out
}

gcomp_pair_contrast <- function(gcomp, group_a, group_b, label) {
  a <- gcomp[gcomp$transition_group == group_a, ]
  b <- gcomp[gcomp$transition_group == group_b, ]
  draws <- attr(gcomp, "draws")
  draw_a <- draws[draws$transition_group == group_a, c("bootstrap", "adjusted_risk")]
  draw_b <- draws[draws$transition_group == group_b, c("bootstrap", "adjusted_risk")]
  pair_draws <- merge(draw_a, draw_b, by = "bootstrap", suffixes = c("_a", "_b"))
  rr_draw <- pair_draws$adjusted_risk_a / pair_draws$adjusted_risk_b
  rd_draw <- (pair_draws$adjusted_risk_a - pair_draws$adjusted_risk_b) * 1000
  data.frame(
    contrast = label,
    numerator = unname(group_labels[group_a]),
    denominator = unname(group_labels[group_b]),
    marginal_rr = a$adjusted_risk / b$adjusted_risk,
    rr_ci_lower = unname(quantile(rr_draw, 0.025, na.rm = TRUE)),
    rr_ci_upper = unname(quantile(rr_draw, 0.975, na.rm = TRUE)),
    risk_difference_per_1000 = (a$adjusted_risk - b$adjusted_risk) * 1000,
    rd_ci_lower_per_1000 = unname(quantile(rd_draw, 0.025, na.rm = TRUE)),
    rd_ci_upper_per_1000 = unname(quantile(rd_draw, 0.975, na.rm = TRUE)),
    valid_bootstrap_replicates = length(unique(pair_draws$bootstrap)),
    row.names = NULL
  )
}

model_terms <- list(
  model1_crude = "transition_group",
  model2_demographic = "transition_group + age_per10 + sex + hukou + education + marital",
  model3_full = "transition_group + age_per10 + sex + hukou + education + marital + bmi_2011 + current_smoking + alcohol_last_year + diabetes + dyslipidemia"
)

df <- read.csv(file.path(data_dir, "sci_landmark_analysis_dataset.csv"), stringsAsFactors = FALSE, colClasses = c(ID = "character", ID_w1 = "character", cluster_psu = "character"))
landmark_all <- read.csv(file.path(data_dir, "sci_landmark_all_before_2018_followup.csv"), stringsAsFactors = FALSE, colClasses = c(ID = "character", ID_w1 = "character", cluster_psu = "character"))
df <- prep_factors(df)
landmark_all <- prep_factors(landmark_all)

# Rebuild selected 2015 landmark covariates from the raw 2015 health file because
# some chronic-disease fields in the first derived SCI dataset were incomplete.
raw_2015_health <- file.path(raw_root, "2015", "Health_Status_and_Functioning.dta")
if (requireNamespace("haven", quietly = TRUE) && file.exists(raw_2015_health)) {
  h15 <- as.data.frame(haven::read_dta(
    raw_2015_health,
    col_select = c("ID", "zda007_2_", "zda007_3_", "da007_2_", "da007_3_", "da007_w2_2_2_", "da007_w2_2_3_", "zda059", "da059")
  ))
  h15$ID <- as.character(h15$ID)
  is_yes <- function(x) !is.na(x) & x == 1
  is_no <- function(x) !is.na(x) & x == 2
  h15$bmc_dyslipidemia_2015 <- as.integer(is_yes(h15$zda007_2_) | is_yes(h15$da007_2_) | is_yes(h15$da007_w2_2_2_))
  h15$bmc_diabetes_2015 <- as.integer(is_yes(h15$zda007_3_) | is_yes(h15$da007_3_) | is_yes(h15$da007_w2_2_3_))
  h15$bmc_current_smoking_2015 <- ifelse(
    is_yes(h15$zda059) | is_yes(h15$da059), 1,
    ifelse(is_no(h15$zda059) | is_no(h15$da059), 0, NA)
  )
  h15 <- h15[, c("ID", "bmc_dyslipidemia_2015", "bmc_diabetes_2015", "bmc_current_smoking_2015")]
  df <- merge(df, h15, by = "ID", all.x = TRUE)
} else {
  df$bmc_dyslipidemia_2015 <- if ("dyslipidemia_2015" %in% names(df)) df$dyslipidemia_2015 else NA
  df$bmc_diabetes_2015 <- if ("diabetes_2015" %in% names(df)) df$diabetes_2015 else NA
  df$bmc_current_smoking_2015 <- if ("current_smoking_2015" %in% names(df)) df$current_smoking_2015 else NA
}

# Main RR models: modified Poisson with PSU/community cluster-robust SEs.
rr_rows <- list()
fits <- list()
for (mn in names(model_terms)) {
  obj <- fit_complete(df, "incident_cvd_2015_2018", model_terms[[mn]], poisson(link = "log"), extra_vars = "cluster_psu")
  vc <- vcov_cluster(obj$fit, obj$data)
  fits[[mn]] <- list(fit = obj$fit, data = obj$data, vc = vc)
  rr_rows[[length(rr_rows) + 1]] <- extract_transition(obj$fit, obj$data, "incident_cvd_2015_2018", mn, vc, "RR")
}
rr_all <- do.call(rbind, rr_rows)
write.csv(rr_all, file.path(table_dir, "bmc_table2_poisson_rr_all_models_cluster.csv"), row.names = FALSE)

event_summary <- read.csv(file.path(table_dir, "table_transition_group_event_summary.csv"), stringsAsFactors = FALSE)
event_summary <- event_summary[match(group_order, event_summary$transition_group), ]

model_wide <- rr_all[, c("transition_group", "transition_group_label", "model", "estimate_ci")]
model_wide <- reshape(model_wide, idvar = c("transition_group", "transition_group_label"), timevar = "model", direction = "wide")
names(model_wide) <- sub("^estimate_ci\\.", "", names(model_wide))
model_wide <- model_wide[match(group_order[-1], model_wide$transition_group), ]
table2 <- merge(model_wide, event_summary[, c("transition_group", "n", "cvd_events")], by = "transition_group", all.x = TRUE, sort = FALSE)
table2 <- table2[match(group_order[-1], table2$transition_group), ]
table2$`Events/participants` <- sprintf("%d/%d", as.integer(table2$cvd_events), as.integer(table2$n))
table2 <- table2[, c("transition_group_label", "Events/participants", "model1_crude", "model2_demographic", "model3_full")]
names(table2) <- c("Transition group", "Events/participants", "Crude RR (95% CI)", "Model 2 RR (95% CI)", "Fully adjusted RR (95% CI)")
write.csv(table2, file.path(table_dir, "bmc_main_table2_poisson_rr_three_models.csv"), row.names = FALSE)

gcomp <- cluster_bootstrap_gcomp(df, model_terms$model3_full, n_boot = 1000)
write.csv(gcomp, file.path(table_dir, "bmc_gcomp_cluster_bootstrap_absolute_risk.csv"), row.names = FALSE)
table3 <- gcomp[match(group_order, gcomp$transition_group), ]
table3$transition_group_label <- unname(group_labels[table3$transition_group])
table3$`Marginal RR (95% CI)` <- ifelse(
  table3$transition_group == "persistent_nonhypertensive",
  "1.00 (reference)",
  sprintf("%.2f (%.2f-%.2f)", table3$marginal_rr, table3$rr_ci_lower, table3$rr_ci_upper)
)
table3$`Adjusted risk per 1000 (95% CI)` <- sprintf("%.1f (%.1f-%.1f)", table3$adjusted_risk * 1000, table3$risk_ci_lower * 1000, table3$risk_ci_upper * 1000)
table3$`Risk difference per 1000 (95% CI)` <- ifelse(
  table3$transition_group == "persistent_nonhypertensive",
  "0.0 (reference)",
  sprintf("%.1f (%.1f-%.1f)", table3$risk_difference * 1000, table3$rd_ci_lower * 1000, table3$rd_ci_upper * 1000)
)
table3 <- table3[, c("transition_group_label", "Marginal RR (95% CI)", "Adjusted risk per 1000 (95% CI)", "Risk difference per 1000 (95% CI)")]
names(table3)[1] <- "Transition group"
write.csv(table3, file.path(table_dir, "bmc_main_table3_rr_absolute_risk.csv"), row.names = FALSE)

global <- global_wald(fits$model3_full$fit, fits$model3_full$vc)
write.csv(global, file.path(table_dir, "bmc_global_wald_transition_group.csv"), row.names = FALSE)

# Prespecified clinical contrasts among hypertension states.
poisson_contrasts <- rbind(
  contrast_rr(fits$model3_full$fit, fits$model3_full$vc, "persistent_aware_or_treated_uncontrolled", "persistent_controlled", "Persistent awareness/treatment without control vs persistent controlled"),
  contrast_rr(fits$model3_full$fit, fits$model3_full$vc, "persistent_aware_or_treated_uncontrolled", "gained_control", "Persistent awareness/treatment without control vs gained control"),
  contrast_rr(fits$model3_full$fit, fits$model3_full$vc, "control_loss", "persistent_controlled", "Loss of control vs persistent controlled")
)
poisson_contrasts$p_holm <- p.adjust(poisson_contrasts$p_value, method = "holm")
poisson_contrasts$p_holm_text <- vapply(poisson_contrasts$p_holm, format_p, character(1))
poisson_contrasts$modified_poisson_rr_95ci <- poisson_contrasts$rr_95ci

gcomp_contrasts <- rbind(
  gcomp_pair_contrast(gcomp, "persistent_aware_or_treated_uncontrolled", "persistent_controlled", "Persistent awareness/treatment without control vs persistent controlled"),
  gcomp_pair_contrast(gcomp, "persistent_aware_or_treated_uncontrolled", "gained_control", "Persistent awareness/treatment without control vs gained control"),
  gcomp_pair_contrast(gcomp, "control_loss", "persistent_controlled", "Loss of control vs persistent controlled")
)
gcomp_contrasts$marginal_rr_95ci <- sprintf("%.2f (%.2f-%.2f)", gcomp_contrasts$marginal_rr, gcomp_contrasts$rr_ci_lower, gcomp_contrasts$rr_ci_upper)
gcomp_contrasts$rd_95ci_per_1000 <- sprintf("%.1f (%.1f-%.1f)", gcomp_contrasts$risk_difference_per_1000, gcomp_contrasts$rd_ci_lower_per_1000, gcomp_contrasts$rd_ci_upper_per_1000)
contrasts <- merge(
  gcomp_contrasts,
  poisson_contrasts[, c("contrast", "modified_poisson_rr_95ci", "p_value", "p_text", "p_holm", "p_holm_text")],
  by = "contrast",
  all.x = TRUE,
  sort = FALSE
)
write.csv(contrasts, file.path(table_dir, "bmc_supp_table_prespecified_clinical_contrasts.csv"), row.names = FALSE)

# 2015 landmark-covariate sensitivity model.
df$age_2015_per10 <- (df$age_2011 + 4) / 10
landmark_rhs <- "transition_group + age_2015_per10 + sex + hukou + education + marital + bmi_2015 + bmc_current_smoking_2015 + alcohol_last_year_2015 + bmc_diabetes_2015 + bmc_dyslipidemia_2015"
landmark_vars <- unique(c(all.vars(as.formula(paste("incident_cvd_2015_2018 ~", landmark_rhs))), "cluster_psu"))
landmark_cc <- df[complete.cases(df[, landmark_vars]), landmark_vars, drop = FALSE]
if (nrow(landmark_cc) > 0) {
  landmark_obj <- fit_complete(df, "incident_cvd_2015_2018", landmark_rhs, poisson(link = "log"), extra_vars = "cluster_psu")
  landmark_vc <- vcov_cluster(landmark_obj$fit, landmark_obj$data)
  landmark_out <- extract_transition(landmark_obj$fit, landmark_obj$data, "incident_cvd_2015_2018", "landmark_2015_covariate_poisson_cluster", landmark_vc, "RR")
} else {
  landmark_out <- data.frame(
    note = "2015 landmark-covariate model could not be estimated because 2015 covariates were unavailable after complete-case filtering."
  )
}
write.csv(landmark_out, file.path(table_dir, "bmc_supp_table_2015_landmark_covariate_rr.csv"), row.names = FALSE)

# Western-medicine-only treatment definition sensitivity analysis.
raw_2011_health <- file.path(raw_root, "2011", "health_status_and_functioning.dta")
raw_2011_bio <- file.path(raw_root, "2011", "biomarkers.dta")
raw_2015_bio <- file.path(raw_root, "2015", "Biomarker.dta")
if (requireNamespace("haven", quietly = TRUE) && file.exists(raw_2011_health) && file.exists(raw_2011_bio) && file.exists(raw_2015_health) && file.exists(raw_2015_bio)) {
  h11_west <- as.data.frame(haven::read_dta(raw_2011_health, col_select = c("ID", "da007_1_", "da008_1_", "da011s1", "da011s2")))
  b11_west <- as.data.frame(haven::read_dta(raw_2011_bio, col_select = c("ID", "qa003", "qa004", "qa007", "qa008", "qa011", "qa012")))
  h15_west <- as.data.frame(haven::read_dta(
    raw_2015_health,
    col_select = c("ID", "zda007_1_", "zda008_1_", "da007_1_", "da007_w2_1_1_", "da007_w2_2_1_", "da011s1", "da011s2")
  ))
  b15_west <- as.data.frame(haven::read_dta(raw_2015_bio, col_select = c("ID", "qa003", "qa004", "qa007", "qa008", "qa011", "qa012")))
  h11_west$ID <- as.character(h11_west$ID)
  b11_west$ID <- as.character(b11_west$ID)
  h15_west$ID <- as.character(h15_west$ID)
  b15_west$ID <- as.character(b15_west$ID)
  w11_west <- merge(h11_west, b11_west, by = "ID")
  w11_west$ID <- pad_charls_id(convert_wave1_to_wave2_id(w11_west$ID))
  # Plausibility bounds must match valid_sbp()/valid_dbp() in
  # 01_build_sci_landmark_dataset.py (SBP 50-260, DBP 30-160).
  w11_west$sbp_2011_western <- mean_valid(w11_west, c("qa003", "qa007", "qa011"), 50, 260)
  w11_west$dbp_2011_western <- mean_valid(w11_west, c("qa004", "qa008", "qa012"), 30, 160)
  diagnosis_2011_western <- w11_west$da007_1_
  awareness_2011_western <- w11_west$da008_1_
  chinese_med_2011 <- w11_west$da011s1
  western_med_2011 <- w11_west$da011s2
  w11_west$htn_state_2011_western_only <- classify_htn_state(
    w11_west$sbp_2011_western,
    w11_west$dbp_2011_western,
    diagnosis_2011_western,
    awareness_2011_western,
    chinese_med_2011,
    western_med_2011,
    TRUE
  )
  w15_west <- merge(h15_west, b15_west, by = "ID")
  w15_west$ID <- pad_charls_id(w15_west$ID)
  w15_west$sbp_2015_western <- mean_valid(w15_west, c("qa003", "qa007", "qa011"), 50, 260)
  w15_west$dbp_2015_western <- mean_valid(w15_west, c("qa004", "qa008", "qa012"), 30, 160)
  w15_west$htn_diag_2015_western <- reconstructed_2015_status_r(w15_west, 1)
  awareness_2015_western <- w15_west$zda008_1_
  chinese_med_2015 <- w15_west$da011s1
  western_med_2015 <- w15_west$da011s2
  w15_west$htn_state_2015_western_only <- classify_htn_state(
    w15_west$sbp_2015_western,
    w15_west$dbp_2015_western,
    w15_west$htn_diag_2015_western,
    awareness_2015_western,
    chinese_med_2015,
    western_med_2015,
    TRUE
  )
  # Self-check: with western_only = FALSE this classifier must reproduce the
  # primary state assignment produced by 01_build_sci_landmark_dataset.py.
  # Any mismatch means the sensitivity analysis is changing more than the
  # treatment definition, so it must fail loudly rather than be reported.
  w11_west$htn_state_2011_maindef <- classify_htn_state(
    w11_west$sbp_2011_western,
    w11_west$dbp_2011_western,
    diagnosis_2011_western,
    awareness_2011_western,
    chinese_med_2011,
    western_med_2011,
    FALSE
  )
  w15_west$htn_state_2015_maindef <- classify_htn_state(
    w15_west$sbp_2015_western,
    w15_west$dbp_2015_western,
    w15_west$htn_diag_2015_western,
    awareness_2015_western,
    chinese_med_2015,
    western_med_2015,
    FALSE
  )
  disagree <- function(a, b) (is.na(a) != is.na(b)) | (!is.na(a) & !is.na(b) & a != b)
  check_state <- function(py_state, r_ids, r_state) {
    merged <- merge(
      data.frame(ID_join = pad_charls_id(df$ID), py = py_state, stringsAsFactors = FALSE),
      data.frame(ID_join = pad_charls_id(r_ids), r = r_state, stringsAsFactors = FALSE),
      by = "ID_join"
    )
    sum(disagree(merged$py, merged$r))
  }
  mismatch_2011 <- check_state(df$htn_state_2011, w11_west$ID, w11_west$htn_state_2011_maindef)
  mismatch_2015 <- check_state(df$htn_state_2015, w15_west$ID, w15_west$htn_state_2015_maindef)
  cat(sprintf(
    "State-classifier self-check (western_only = FALSE): 2011 mismatches = %d, 2015 mismatches = %d\n",
    mismatch_2011, mismatch_2015
  ))
  if (mismatch_2011 > 0 || mismatch_2015 > 0) {
    stop("classify_htn_state() does not reproduce the primary state assignment; western-medicine-only analysis aborted.")
  }

  western_states <- merge(w11_west[, c("ID", "htn_state_2011_western_only")], w15_west[, c("ID", "htn_state_2015_western_only")], by = "ID")
  western_states$transition_group_western_only <- map_transition_group(western_states$htn_state_2011_western_only, western_states$htn_state_2015_western_only)
  df_west_base <- df
  df_west_base$ID_join <- pad_charls_id(df_west_base$ID)
  western_states$ID_join <- pad_charls_id(western_states$ID)
  df_west <- merge(df_west_base, western_states[, c("ID_join", "transition_group_western_only")], by = "ID_join", all.x = TRUE)
  df_west$transition_group_original <- df_west$transition_group
  df_west$transition_group <- factor(df_west$transition_group_western_only, levels = group_order)
  western_obj <- fit_complete(df_west, "incident_cvd_2015_2018", model_terms$model3_full, poisson(link = "log"), extra_vars = "cluster_psu")
  western_vc <- vcov_cluster(western_obj$fit, western_obj$data)
  western_out <- extract_transition(western_obj$fit, western_obj$data, "incident_cvd_2015_2018", "western_medicine_only_treatment_definition", western_vc, "RR")
  western_count_dat <- df_west[!is.na(df_west$transition_group), c("transition_group", "incident_cvd_2015_2018")]
  western_count_dat$n <- 1
  western_count_dat$events <- western_count_dat$incident_cvd_2015_2018
  western_counts <- aggregate(cbind(n, events) ~ transition_group, data = western_count_dat, FUN = sum)
  write.csv(western_out, file.path(table_dir, "bmc_supp_table_western_medicine_only_rr.csv"), row.names = FALSE)
  write.csv(western_counts, file.path(table_dir, "bmc_supp_table_western_medicine_only_counts.csv"), row.names = FALSE)
}

# Logistic OR as supplementary model.
logit_obj <- fit_complete(df, "incident_cvd_2015_2018", model_terms$model3_full, binomial(), extra_vars = "cluster_psu")
logit_vc <- vcov_cluster(logit_obj$fit, logit_obj$data)
logit_or <- extract_transition(logit_obj$fit, logit_obj$data, "incident_cvd_2015_2018", "supplementary_logistic_cluster_model3", logit_vc, "OR")
write.csv(logit_or, file.path(table_dir, "bmc_supp_table_logistic_or_cluster_model3.csv"), row.names = FALSE)

# RR subgroup analyses and interaction tests.
subgroup_rows <- list()
base_terms_vec <- c("transition_group", "age_per10", "sex", "hukou", "education", "marital", "bmi_2011", "current_smoking", "alcohol_last_year", "diabetes", "dyslipidemia")
for (sg in c("sex", "age_group")) {
  cat("Fitting interaction subgroup model:", sg, "\n")
  base_rhs <- paste(base_terms_vec, collapse = " + ")
  no_sg_terms <- if (sg %in% base_terms_vec) setdiff(base_terms_vec, sg) else base_terms_vec
  no_sg_rhs <- paste(no_sg_terms, collapse = " + ")
  int_full <- as.formula(paste("incident_cvd_2015_2018 ~", paste(unique(c(no_sg_terms, sg, paste0("transition_group:", sg))), collapse = " + ")))
  int_reduced <- as.formula(paste("incident_cvd_2015_2018 ~", base_rhs))
  vars <- unique(c(all.vars(int_full), all.vars(int_reduced), "cluster_psu"))
  dat_int <- df[complete.cases(df[, vars]), vars]
  full_fit <- glm(int_full, data = dat_int, family = poisson(link = "log"))
  red_fit <- glm(int_reduced, data = dat_int, family = poisson(link = "log"))
  p_int <- anova(red_fit, full_fit, test = "Chisq")$`Pr(>Chi)`[2]
  for (lev in levels(df[[sg]])) {
    cat("Fitting stratified subgroup model:", sg, lev, "\n")
    sub <- df[df[[sg]] == lev, ]
    rhs <- no_sg_rhs
    obj <- fit_complete(sub, "incident_cvd_2015_2018", rhs, poisson(link = "log"), extra_vars = "cluster_psu")
    vc <- vcov_cluster(obj$fit, obj$data)
    out <- extract_transition(obj$fit, obj$data, "incident_cvd_2015_2018", paste0("subgroup_", sg, "_", lev), vc, "RR")
    out$subgroup_variable <- sg
    out$subgroup_level <- lev
    out$interaction_p <- p_int
    out$interaction_p_text <- format_p(p_int)
    subgroup_rows[[length(subgroup_rows) + 1]] <- out
  }
}
subgroups <- do.call(rbind, subgroup_rows)
write.csv(subgroups, file.path(table_dir, "bmc_supp_table_subgroup_sex_age_rr_interaction.csv"), row.names = FALSE)

# IPCW diagnostics.
landmark_all$observed_binary <- landmark_all$observed_2018_outcome %in% c(TRUE, "TRUE", "True", "true", 1, "1")
ipcw_rhs <- "transition_group + age_per10 + sex + hukou + education + marital + bmi_2011 + current_smoking + alcohol_last_year + diabetes + dyslipidemia + sbp_2011 + dbp_2011 + sbp_2015 + dbp_2015"
obs_obj <- fit_complete(landmark_all, "observed_binary", ipcw_rhs, binomial())
obs_data <- obs_obj$data
obs_data$pred_observed <- predict(obs_obj$fit, type = "response")
stab <- mean(obs_data$observed_binary == 1)
obs_data$ipcw_raw <- ifelse(obs_data$observed_binary == 1, stab / obs_data$pred_observed, NA_real_)
iq <- quantile(obs_data$ipcw_raw, c(0.01, 0.99), na.rm = TRUE)
obs_data$ipcw_trunc <- pmin(pmax(obs_data$ipcw_raw, iq[1]), iq[2])
obs_obs <- obs_data[obs_data$observed_binary == 1, ]
ipcw_diag <- data.frame(
  n_landmark_model = nrow(obs_data),
  n_observed_model = nrow(obs_obs),
  raw_mean = mean(obs_obs$ipcw_raw),
  raw_median = median(obs_obs$ipcw_raw),
  raw_min = min(obs_obs$ipcw_raw),
  raw_max = max(obs_obs$ipcw_raw),
  trunc_p1 = iq[1],
  trunc_p99 = iq[2],
  trunc_mean = mean(obs_obs$ipcw_trunc),
  trunc_median = median(obs_obs$ipcw_trunc),
  trunc_min = min(obs_obs$ipcw_trunc),
  trunc_max = max(obs_obs$ipcw_trunc),
  effective_sample_size = sum(obs_obs$ipcw_trunc)^2 / sum(obs_obs$ipcw_trunc^2)
)
write.csv(ipcw_diag, file.path(table_dir, "bmc_supp_table_ipcw_diagnostics.csv"), row.names = FALSE)

# Death-related loss to follow-up diagnostics and conservative crude scenarios.
sample_info_2018 <- file.path(raw_root, "2018", "Sample_Infor.dta")
exit_module_2020 <- file.path(raw_root, "2020", "Exit_Module.dta")
if (requireNamespace("haven", quietly = TRUE) && file.exists(sample_info_2018)) {
  si <- as.data.frame(haven::read_dta(sample_info_2018, col_select = c("ID", "died", "crosssection")))
  si$ID_join <- pad_charls_id(si$ID)
  death_base <- landmark_all
  death_base$ID_join <- pad_charls_id(death_base$ID)
  death_df <- merge(death_base, si[, c("ID_join", "died", "crosssection")], by = "ID_join", all.x = TRUE)
  death_df$observed_binary <- death_df$observed_2018_outcome %in% c(TRUE, "TRUE", "True", "true", 1, "1")
  death_df$lost <- !death_df$observed_binary
  death_df$loss_status <- ifelse(
    !death_df$lost, "Observed 2018 CVD outcome",
    ifelse(!is.na(death_df$died) & death_df$died == 1, "Died before 2018 outcome assessment",
      ifelse(is.na(death_df$died) & is.na(death_df$crosssection), "No 2018 sample-information status", "In 2018 roster but no CVD outcome")
    )
  )
  loss_summary <- as.data.frame(table(death_df$loss_status), stringsAsFactors = FALSE)
  names(loss_summary) <- c("loss_status", "n")
  write.csv(loss_summary, file.path(table_dir, "bmc_supp_table_death_loss_status.csv"), row.names = FALSE)

  death_by_group <- aggregate(
    cbind(
      landmark_n = rep(1, nrow(death_df)),
      observed = as.integer(death_df$observed_binary),
      lost = as.integer(death_df$lost),
      died = as.integer(death_df$lost & !is.na(death_df$died) & death_df$died == 1),
      no_2018_status = as.integer(death_df$lost & is.na(death_df$died) & is.na(death_df$crosssection)),
      roster_no_outcome = as.integer(death_df$lost & !is.na(death_df$died) & death_df$died != 1)
    ),
    by = list(transition_group = death_df$transition_group, transition_group_label = death_df$transition_group_label),
    FUN = sum
  )
  death_by_group$death_rate_pct <- death_by_group$died / death_by_group$landmark_n * 100
  write.csv(death_by_group, file.path(table_dir, "bmc_supp_table_death_by_transition_group.csv"), row.names = FALSE)

  core <- "persistent_aware_or_treated_uncontrolled"
  ref <- "persistent_nonhypertensive"
  scenario_rr <- function(dat, scenario, include_rows, event_var) {
    sub <- dat[include_rows, ]
    one <- sub[sub$transition_group == core, ]
    zero <- sub[sub$transition_group == ref, ]
    risk_core <- mean(one[[event_var]])
    risk_ref <- mean(zero[[event_var]])
    data.frame(
      scenario = scenario,
      reference_group = ref,
      comparison_group = core,
      reference_n = nrow(zero),
      reference_events = sum(zero[[event_var]]),
      comparison_n = nrow(one),
      comparison_events = sum(one[[event_var]]),
      crude_rr = risk_core / risk_ref
    )
  }
  death_df$event_observed <- ifelse(death_df$observed_binary, death_df$incident_cvd_2015_2018, NA)
  death_df$event_death_non_event <- ifelse(death_df$observed_binary, death_df$incident_cvd_2015_2018, 0)
  death_df$event_death_event <- ifelse(death_df$observed_binary, death_df$incident_cvd_2015_2018, as.integer(!is.na(death_df$died) & death_df$died == 1))
  death_df$event_all_lost_non_event <- ifelse(death_df$observed_binary, death_df$incident_cvd_2015_2018, 0)
  death_df$event_all_lost_event <- ifelse(death_df$observed_binary, death_df$incident_cvd_2015_2018, 1)
  scenarios <- rbind(
    scenario_rr(death_df, "Main analysis: participants without 2018 outcome excluded", death_df$observed_binary, "event_observed"),
    scenario_rr(death_df, "Deaths counted as non-events; other missing outcomes excluded", death_df$observed_binary | (!is.na(death_df$died) & death_df$died == 1), "event_death_non_event"),
    scenario_rr(death_df, "Deaths counted as events; other missing outcomes excluded", death_df$observed_binary | (!is.na(death_df$died) & death_df$died == 1), "event_death_event"),
    scenario_rr(death_df, "All missing 2018 outcomes counted as non-events", rep(TRUE, nrow(death_df)), "event_all_lost_non_event"),
    scenario_rr(death_df, "All missing 2018 outcomes counted as events", rep(TRUE, nrow(death_df)), "event_all_lost_event")
  )
  scenarios$crude_rr <- round(scenarios$crude_rr, 2)
  write.csv(scenarios, file.path(table_dir, "bmc_supp_table_death_extreme_scenarios.csv"), row.names = FALSE)
}
if (requireNamespace("haven", quietly = TRUE) && file.exists(exit_module_2020)) {
  ex <- as.data.frame(haven::read_dta(exit_module_2020, col_select = c("ID", "exb001_1", "exb001_2", "exb001_3")))
  ex$ID_join <- pad_charls_id(ex$ID)
  ex$death_year <- as.numeric(ex$exb001_1)
  exit_base <- landmark_all[, c("ID", "transition_group", "transition_group_label")]
  exit_base$ID_join <- pad_charls_id(exit_base$ID)
  ex_join <- merge(exit_base, ex[, c("ID_join", "death_year")], by = "ID_join", all.x = TRUE)
  death_year_summary <- aggregate(
    rep(1, nrow(ex_join[!is.na(ex_join$death_year), ])),
    by = list(
      transition_group = ex_join$transition_group[!is.na(ex_join$death_year)],
      transition_group_label = ex_join$transition_group_label[!is.na(ex_join$death_year)],
      death_year = ex_join$death_year[!is.na(ex_join$death_year)]
    ),
    FUN = sum
  )
  names(death_year_summary)[names(death_year_summary) == "x"] <- "n"
  write.csv(death_year_summary, file.path(table_dir, "bmc_supp_table_exit_module_death_years.csv"), row.names = FALSE)
}

mice_details <- data.frame(
  item = c("Number of imputations", "Iterations", "Imputed variables", "Variables not imputed", "Variables included in imputation model", "Combination rule"),
  value = c(
    "20",
    "20",
    "Missing baseline covariates: sex, hukou, BMI, current smoking, alcohol drinking in the past year, diabetes, and dyslipidemia.",
    "The primary exposure transition group and incident CVD outcome were not imputed.",
    "The imputation model included exposure group, incident CVD, age, sex, hukou, education, marital status, BMI, smoking, alcohol drinking, diabetes, dyslipidemia, 2011/2015 SBP and DBP, and heart disease and stroke outcomes.",
    "Model estimates were pooled using Rubin's rules."
  )
)
write.csv(mice_details, file.path(table_dir, "bmc_supp_table_mice_details.csv"), row.names = FALSE)

cat("BMC RR-primary reanalysis complete.\n")
cat("Main RR table: bmc_main_table2_poisson_rr_three_models.csv\n")
cat("Global Wald test: bmc_global_wald_transition_group.csv\n")
cat("Prespecified contrasts: bmc_supp_table_prespecified_clinical_contrasts.csv\n")
