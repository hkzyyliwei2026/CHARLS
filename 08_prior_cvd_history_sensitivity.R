# Prior-CVD-history sensitivity analysis.
#
# The primary cohort excludes participants reporting heart disease or stroke in
# 2011 or at the 2015 landmark. Self-reported chronic conditions are unstable
# across CHARLS waves, so a participant who reported CVD in 2013 but not in 2015
# and again in 2018 would be counted as an incident case. This script flags such
# participants using the 2013 wave and refits the fully adjusted primary model
# after excluding them.

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
  gained_control = "Gained control",
  persistent_controlled = "Persistent control",
  control_loss = "Loss of control",
  other_transition = "Other transitions"
)

yes <- function(x, code = 1) !is.na(x) & x == code

pad_charls_id <- function(x) {
  x <- trimws(as.character(x))
  missing <- is.na(x) | x == "NA"
  x[!missing] <- ifelse(nchar(x[!missing]) < 12, paste0(strrep("0", 12 - nchar(x[!missing])), x[!missing]), x[!missing])
  x[missing] <- NA_character_
  x
}

# Same reconstruction rule used for the 2015 wave in 05_bmc_rr_primary_reanalysis.R.
reconstructed_wave2_status <- function(dat, disease_index) {
  previous_yes <- yes(dat[[paste0("zda007_", disease_index, "_")]], 1)
  confirmation <- dat[[paste0("da007_w2_1_", disease_index, "_")]]
  new_diagnosis <- dat[[paste0("da007_w2_2_", disease_index, "_")]]
  raw_status <- dat[[paste0("da007_", disease_index, "_")]]
  out <- rep(NA_integer_, nrow(dat))
  has_conf <- !is.na(confirmation)
  i <- has_conf & previous_yes
  out[i] <- as.integer(yes(confirmation[i], 1) | (yes(confirmation[i], 2) & yes(new_diagnosis[i], 1)))
  j <- has_conf & !previous_yes
  out[j] <- as.integer(yes(confirmation[j], 2) | (yes(confirmation[j], 1) & yes(new_diagnosis[j], 1)))
  k <- is.na(out) & !is.na(raw_status)
  out[k] <- as.integer(yes(raw_status[k], 1))
  out
}

prep_factors <- function(df) {
  df$transition_group <- factor(df$transition_group, levels = group_order)
  df$sex <- factor(df$sex, levels = c("male", "female"))
  df$hukou <- factor(df$hukou, levels = c("agricultural", "non_agricultural", "other"))
  df$education <- factor(df$education, levels = c("no_formal", "primary_or_below", "middle_school", "high_school_plus"))
  df$marital <- factor(df$marital, levels = c("married_partnered", "not_married_partnered"))
  df
}

RHS <- paste(
  "transition_group + age_per10 + sex + hukou + education + marital + bmi_2011",
  "+ current_smoking + alcohol_last_year + diabetes + dyslipidemia"
)

fit_and_extract <- function(dat, model_name) {
  form <- as.formula(paste("incident_cvd_2015_2018 ~", RHS))
  vars <- unique(c(all.vars(form), "cluster_psu"))
  d <- dat[complete.cases(dat[, vars]), vars, drop = FALSE]
  fit <- glm(form, data = d, family = poisson(link = "log"))
  vc <- sandwich::vcovCL(fit, cluster = d$cluster_psu, type = "HC0")
  beta <- coef(fit); se <- sqrt(diag(vc))
  keep <- startsWith(names(beta), "transition_group")
  groups <- sub("^transition_group", "", names(beta)[keep])
  data.frame(
    model = model_name,
    transition_group = groups,
    transition_group_label = unname(group_labels[groups]),
    n_model = nrow(d),
    events_model = sum(d$incident_cvd_2015_2018 == 1),
    estimate = exp(beta[keep]),
    ci_lower = exp(beta[keep] - 1.96 * se[keep]),
    ci_upper = exp(beta[keep] + 1.96 * se[keep]),
    p_value = 2 * pnorm(abs(beta[keep] / se[keep]), lower.tail = FALSE),
    estimate_ci = sprintf("%.2f (%.2f-%.2f)", exp(beta[keep]),
                          exp(beta[keep] - 1.96 * se[keep]), exp(beta[keep] + 1.96 * se[keep])),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

df <- read.csv(file.path(data_dir, "sci_landmark_analysis_dataset.csv"),
               stringsAsFactors = FALSE,
               colClasses = c(ID = "character", ID_w1 = "character", cluster_psu = "character"))
df$ID_join <- pad_charls_id(df$ID)
df <- prep_factors(df)

h13_path <- file.path(raw_root, "2013", "Health_Status_and_Functioning.dta")
if (!file.exists(h13_path)) stop("2013 health status file not found: ", h13_path)
h13 <- as.data.frame(haven::read_dta(h13_path, col_select = c(
  "ID",
  "zda007_7_", "da007_w2_1_7_", "da007_w2_2_7_", "da007_7_",
  "zda007_8_", "da007_w2_1_8_", "da007_w2_2_8_", "da007_8_"
)))
h13$heart_2013 <- reconstructed_wave2_status(h13, 7)
h13$stroke_2013 <- reconstructed_wave2_status(h13, 8)
h13$ID_join <- pad_charls_id(h13$ID)
h13$cvd_2013_reported <- as.integer(yes(h13$heart_2013, 1) | yes(h13$stroke_2013, 1))

merged <- merge(df, h13[, c("ID_join", "heart_2013", "stroke_2013", "cvd_2013_reported")],
                by = "ID_join", all.x = TRUE)
merged$linked_2013 <- as.integer(!is.na(merged$cvd_2013_reported))
merged$prior_cvd_2013 <- as.integer(yes(merged$cvd_2013_reported, 1))

cat(sprintf("Primary cohort: %d\n", nrow(merged)))
cat(sprintf("Linked to the 2013 wave: %d (%.1f%%)\n",
            sum(merged$linked_2013), 100 * mean(merged$linked_2013)))
cat(sprintf("Reported heart disease or stroke in 2013 despite being CVD-free at the 2015 landmark: %d\n",
            sum(merged$prior_cvd_2013)))
cat(sprintf("  of whom counted as incident CVD in 2018: %d of %d total events (%.1f%%)\n",
            sum(merged$prior_cvd_2013 == 1 & merged$incident_cvd_2015_2018 == 1, na.rm = TRUE),
            sum(merged$incident_cvd_2015_2018 == 1, na.rm = TRUE),
            100 * sum(merged$prior_cvd_2013 == 1 & merged$incident_cvd_2015_2018 == 1, na.rm = TRUE) /
              sum(merged$incident_cvd_2015_2018 == 1, na.rm = TRUE)))

by_group <- do.call(rbind, lapply(group_order, function(g) {
  d <- merged[merged$transition_group == g, ]
  data.frame(
    transition_group = g,
    transition_group_label = unname(group_labels[g]),
    n = nrow(d),
    linked_2013 = sum(d$linked_2013),
    prior_cvd_2013 = sum(d$prior_cvd_2013),
    prior_cvd_2013_pct = round(100 * sum(d$prior_cvd_2013) / max(sum(d$linked_2013), 1), 1),
    events = sum(d$incident_cvd_2015_2018 == 1, na.rm = TRUE),
    events_with_prior_cvd_2013 = sum(d$prior_cvd_2013 == 1 & d$incident_cvd_2015_2018 == 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
print(by_group, row.names = FALSE)

primary <- fit_and_extract(merged, "primary_replication")
restricted <- fit_and_extract(merged[merged$prior_cvd_2013 == 0, ], "excluding_prior_cvd_2013")
linked_only <- fit_and_extract(merged[merged$linked_2013 == 1, ], "restricted_to_2013_linked")
# Isolates the effect of the exclusion itself: same denominator as
# restricted_to_2013_linked, minus the participants flagged in 2013.
linked_clean <- fit_and_extract(
  merged[merged$linked_2013 == 1 & merged$prior_cvd_2013 == 0, ],
  "linked_2013_and_no_prior_cvd"
)

out <- rbind(primary, restricted, linked_only, linked_clean)
write.csv(out, file.path(table_dir, "bmc_supp_table_prior_cvd_2013_sensitivity_rr.csv"), row.names = FALSE)
write.csv(by_group, file.path(table_dir, "bmc_supp_table_prior_cvd_2013_by_group.csv"), row.names = FALSE)

cat("\nFully adjusted RR for persistent awareness or treatment without control:\n")
for (m in unique(out$model)) {
  row <- out[out$model == m & out$transition_group == "persistent_aware_or_treated_uncontrolled", ]
  cat(sprintf("  %-28s n=%d events=%d  RR %s\n", m, row$n_model, row$events_model, row$estimate_ci))
}
cat("\nWritten: bmc_supp_table_prior_cvd_2013_sensitivity_rr.csv, bmc_supp_table_prior_cvd_2013_by_group.csv\n")
