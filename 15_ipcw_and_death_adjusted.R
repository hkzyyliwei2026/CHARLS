# Separate death from true non-response in the IPCW model, and run covariate-adjusted
# death/missing-outcome scenario analyses.
#
# The original IPCW model treated every participant without an observed 2018 outcome
# as censored, mixing 314 deaths with 391 true non-respondents. Death is not loss to
# follow-up: weighting deaths back in estimates a counterfactual in which those
# participants had not died. Here the censoring model is restricted to participants
# known to be alive, and deaths are handled separately by scenario analyses that use
# the same covariate set as the primary model.

set.seed(42)

cmd_args <- commandArgs(FALSE)
file_arg <- cmd_args[grep("--file=", cmd_args)][1]
script_dir <- if (length(file_arg) && !is.na(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)) else getwd()
data_dir <- file.path(script_dir, "data")
table_dir <- file.path(script_dir, "tables")
raw_root <- Sys.getenv("CHARLS_RAW_ROOT", unset = file.path(script_dir, "raw", "_extracted"))
if (!requireNamespace("sandwich", quietly = TRUE)) stop("sandwich is required.")
if (!requireNamespace("haven", quietly = TRUE)) stop("haven is required.")

group_order <- c(
  "persistent_nonhypertensive", "incident_hypertension", "persistent_unaware_or_untreated",
  "new_awareness_or_treatment_uncontrolled", "persistent_aware_or_treated_uncontrolled",
  "gained_control", "persistent_controlled", "control_loss", "other_transition"
)
group_labels <- c(
  persistent_nonhypertensive = "Persistent non-hypertension",
  incident_hypertension = "Incident hypertension",
  persistent_unaware_or_untreated = "Persistent unawareness/no treatment",
  new_awareness_or_treatment_uncontrolled = "New awareness/treatment without control",
  persistent_aware_or_treated_uncontrolled = "Persistent awareness or treatment without control",
  gained_control = "Gained control", persistent_controlled = "Persistent control",
  control_loss = "Loss of control", other_transition = "Other transitions"
)

pad_charls_id <- function(x) {
  x <- trimws(as.character(x))
  miss <- is.na(x) | x == "NA"
  x[!miss] <- ifelse(nchar(x[!miss]) < 12, paste0(strrep("0", 12 - nchar(x[!miss])), x[!miss]), x[!miss])
  x[miss] <- NA_character_
  x
}

prep <- function(d) {
  d$transition_group <- factor(d$transition_group, levels = group_order)
  d$sex <- factor(d$sex, levels = c("male", "female"))
  d$hukou <- factor(d$hukou, levels = c("agricultural", "non_agricultural", "other"))
  d$education <- factor(d$education, levels = c("no_formal", "primary_or_below", "middle_school", "high_school_plus"))
  d$marital <- factor(d$marital, levels = c("married_partnered", "not_married_partnered"))
  d
}

COVS <- "age_per10 + sex + hukou + education + marital + bmi_2011 + current_smoking + alcohol_last_year + diabetes + dyslipidemia"
RHS <- paste("transition_group +", COVS)

fit_extract <- function(dat, outcome, model_name, weights = NULL) {
  form <- as.formula(paste(outcome, "~", RHS))
  vars <- unique(c(all.vars(form), "cluster_psu", weights))
  d <- dat[complete.cases(dat[, vars]), vars, drop = FALSE]
  fit <- if (is.null(weights)) {
    glm(form, data = d, family = poisson(link = "log"))
  } else {
    suppressWarnings(glm(form, data = d, family = poisson(link = "log"), weights = d[[weights]]))
  }
  vc <- sandwich::vcovCL(fit, cluster = d$cluster_psu, type = "HC0")
  b <- coef(fit); se <- sqrt(diag(vc))
  keep <- startsWith(names(b), "transition_group")
  g <- sub("^transition_group", "", names(b)[keep])
  data.frame(
    model = model_name, transition_group = g, transition_group_label = unname(group_labels[g]),
    n_model = nrow(d), events_model = sum(d[[outcome]] == 1),
    estimate = exp(b[keep]), ci_lower = exp(b[keep] - 1.96 * se[keep]), ci_upper = exp(b[keep] + 1.96 * se[keep]),
    p_value = 2 * pnorm(abs(b[keep] / se[keep]), lower.tail = FALSE),
    estimate_ci = sprintf("%.2f (%.2f-%.2f)", exp(b[keep]), exp(b[keep] - 1.96 * se[keep]), exp(b[keep] + 1.96 * se[keep])),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

landmark <- read.csv(file.path(data_dir, "sci_landmark_all_before_2018_followup.csv"),
                     stringsAsFactors = FALSE,
                     colClasses = c(ID = "character", ID_w1 = "character", cluster_psu = "character"))
landmark$ID_join <- pad_charls_id(landmark$ID)
landmark <- prep(landmark)

si_path <- file.path(raw_root, "2018", "Sample_Infor.dta")
si <- as.data.frame(haven::read_dta(si_path, col_select = c("ID", "died")))
si$ID_join <- pad_charls_id(si$ID)
landmark <- merge(landmark, si[, c("ID_join", "died")], by = "ID_join", all.x = TRUE)

landmark$observed <- as.integer(landmark$observed_2018_outcome %in% c(TRUE, "True", "TRUE", 1, "1"))
landmark$dead <- as.integer(!is.na(landmark$died) & landmark$died == 1)
cat(sprintf("Landmark eligible %d = observed %d + died %d + true non-response %d\n",
            nrow(landmark), sum(landmark$observed), sum(landmark$dead),
            sum(landmark$observed == 0 & landmark$dead == 0)))
stopifnot(sum(landmark$observed) == 6975, sum(landmark$dead) == 314)

## ---- 1. IPCW restricted to participants known to be alive -------------------
alive <- landmark[landmark$dead == 0, ]
cens_rhs <- paste("transition_group +", COVS, "+ sbp_2011 + dbp_2011 + sbp_2015 + dbp_2015")
cens_vars <- unique(c(all.vars(as.formula(paste("observed ~", cens_rhs))), "cluster_psu"))
alive_cc <- alive[complete.cases(alive[, cens_vars]), ]
cat(sprintf("Censoring model among survivors: %d of %d with complete predictors; %d observed\n",
            nrow(alive_cc), nrow(alive), sum(alive_cc$observed)))
cens_fit <- glm(as.formula(paste("observed ~", cens_rhs)), data = alive_cc, family = binomial())
alive_cc$p_obs <- predict(cens_fit, type = "response")
alive_cc$ipcw <- 1 / alive_cc$p_obs
obs <- alive_cc[alive_cc$observed == 1, ]
lo <- quantile(obs$ipcw, 0.01); hi <- quantile(obs$ipcw, 0.99)
obs$ipcw_trunc <- pmin(pmax(obs$ipcw, lo), hi)
ess <- sum(obs$ipcw_trunc)^2 / sum(obs$ipcw_trunc^2)
cat(sprintf("Survivor-only IPCW: mean %.3f, median %.3f, range %.3f-%.3f, ESS %.0f (n=%d)\n",
            mean(obs$ipcw_trunc), median(obs$ipcw_trunc), min(obs$ipcw_trunc), max(obs$ipcw_trunc), ess, nrow(obs)))
ipcw_out <- fit_extract(obs, "incident_cvd_2015_2018", "ipcw_survivors_only", weights = "ipcw_trunc")

diag <- data.frame(
  quantity = c("Landmark eligible", "Died before the 2018 interview", "Alive, eligible for the censoring model",
               "Alive with complete censoring-model predictors", "Observed outcome among those",
               "Weight mean", "Weight median", "Weight minimum", "Weight maximum", "Effective sample size"),
  value = c(nrow(landmark), sum(landmark$dead), nrow(alive), nrow(alive_cc), sum(alive_cc$observed),
            round(mean(obs$ipcw_trunc), 3), round(median(obs$ipcw_trunc), 3),
            round(min(obs$ipcw_trunc), 3), round(max(obs$ipcw_trunc), 3), round(ess)),
  stringsAsFactors = FALSE
)

## ---- 2. Covariate-adjusted death / missing-outcome scenarios ----------------
scen <- list()
d_a <- landmark[landmark$observed == 1 | landmark$dead == 1, ]
d_a$y <- ifelse(d_a$observed == 1, d_a$incident_cvd_2015_2018, 0)
scen[["deaths_as_non_events"]] <- d_a

d_b <- d_a
d_b$y <- ifelse(d_b$dead == 1, 1, d_b$incident_cvd_2015_2018)
scen[["cvd_or_death_composite"]] <- d_b

d_c <- landmark
d_c$y <- ifelse(d_c$observed == 1, d_c$incident_cvd_2015_2018, 0)
scen[["all_missing_as_non_events"]] <- d_c

d_d <- landmark
d_d$y <- ifelse(d_d$observed == 1, d_d$incident_cvd_2015_2018, 1)
scen[["all_missing_as_events"]] <- d_d

adj <- do.call(rbind, lapply(names(scen), function(nm) fit_extract(scen[[nm]], "y", nm)))
primary <- fit_extract(landmark[landmark$observed == 1, ], "incident_cvd_2015_2018", "primary_replication")

out <- rbind(primary, adj, ipcw_out)
write.csv(out, file.path(table_dir, "bmc_supp_table_death_scenarios_adjusted_rr.csv"), row.names = FALSE)
write.csv(diag, file.path(table_dir, "bmc_supp_table_ipcw_survivors_diagnostics.csv"), row.names = FALSE)

cat("\nFully adjusted RR for persistent awareness or treatment without control:\n")
for (m in unique(out$model)) {
  r <- out[out$model == m & out$transition_group == "persistent_aware_or_treated_uncontrolled", ]
  cat(sprintf("  %-28s n=%5d events=%4d  RR %s\n", m, r$n_model, r$events_model, r$estimate_ci))
}
cat("\nWritten: bmc_supp_table_death_scenarios_adjusted_rr.csv, bmc_supp_table_ipcw_survivors_diagnostics.csv\n")
