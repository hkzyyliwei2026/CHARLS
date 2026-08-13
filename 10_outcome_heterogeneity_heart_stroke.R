# Formal test of effect heterogeneity between incident heart disease and incident stroke.
#
# The two secondary outcomes can occur in the same participant, so the coefficients
# from two separate models cannot simply be compared. Participant records are stacked
# (one row per outcome type) and a modified Poisson model with a full outcome-type
# interaction is fitted, with community-PSU cluster-robust standard errors. Because
# both rows of a participant fall in the same cluster, this also accounts for the
# within-participant correlation. Heterogeneity is assessed with a Wald test on the
# transition group by outcome type interaction terms.

set.seed(42)

cmd_args <- commandArgs(FALSE)
file_arg <- cmd_args[grep("--file=", cmd_args)][1]
script_dir <- if (length(file_arg) && !is.na(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)) else getwd()
data_dir <- file.path(script_dir, "data")
table_dir <- file.path(script_dir, "tables")
if (!requireNamespace("sandwich", quietly = TRUE)) stop("sandwich is required.")

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

df <- read.csv(file.path(data_dir, "sci_landmark_analysis_dataset.csv"),
               stringsAsFactors = FALSE,
               colClasses = c(ID = "character", ID_w1 = "character", cluster_psu = "character"))
df$transition_group <- factor(df$transition_group, levels = group_order)
df$sex <- factor(df$sex, levels = c("male", "female"))
df$hukou <- factor(df$hukou, levels = c("agricultural", "non_agricultural", "other"))
df$education <- factor(df$education, levels = c("no_formal", "primary_or_below", "middle_school", "high_school_plus"))
df$marital <- factor(df$marital, levels = c("married_partnered", "not_married_partnered"))

COVS <- c("age_per10", "sex", "hukou", "education", "marital", "bmi_2011",
          "current_smoking", "alcohol_last_year", "diabetes", "dyslipidemia")
keep <- c("ID", "cluster_psu", "transition_group", COVS,
          "incident_heart_2015_2018", "incident_stroke_2015_2018")
d <- df[complete.cases(df[, keep]), keep]

heart <- d; heart$outcome_type <- "heart"; heart$y <- heart$incident_heart_2015_2018
stroke <- d; stroke$outcome_type <- "stroke"; stroke$y <- stroke$incident_stroke_2015_2018
stacked <- rbind(heart, stroke)
stacked$outcome_type <- factor(stacked$outcome_type, levels = c("heart", "stroke"))

cat(sprintf("Participants: %d | stacked rows: %d | heart events: %d | stroke events: %d\n",
            nrow(d), nrow(stacked), sum(heart$y), sum(stroke$y)))

rhs <- paste("outcome_type * (transition_group +", paste(COVS, collapse = " + "), ")")
fit <- glm(as.formula(paste("y ~", rhs)), data = stacked, family = poisson(link = "log"))
vc <- sandwich::vcovCL(fit, cluster = stacked$cluster_psu, type = "HC0")
beta <- coef(fit)

inter <- grep("^outcome_typestroke:transition_group", names(beta), value = TRUE)
idx <- match(inter, names(beta))
ok <- !is.na(beta[inter])
inter <- inter[ok]; idx <- idx[ok]
stat <- as.numeric(t(beta[inter]) %*% solve(vc[idx, idx, drop = FALSE]) %*% beta[inter])
dfree <- length(inter)
p_global <- pchisq(stat, df = dfree, lower.tail = FALSE)
cat(sprintf("\nGlobal Wald test for transition group by outcome type interaction: chi2 = %.2f, df = %d, P = %.3f\n",
            stat, dfree, p_global))

se <- sqrt(diag(vc))
rows <- lapply(group_order[-1], function(g) {
  b_heart <- paste0("transition_group", g)
  b_int <- paste0("outcome_typestroke:transition_group", g)
  if (!(b_heart %in% names(beta)) || !(b_int %in% names(beta))) return(NULL)
  i1 <- match(b_heart, names(beta)); i2 <- match(b_int, names(beta))
  v_stroke <- vc[i1, i1] + vc[i2, i2] + 2 * vc[i1, i2]
  b_stroke <- beta[i1] + beta[i2]
  data.frame(
    transition_group = g,
    transition_group_label = unname(group_labels[g]),
    heart_rr = sprintf("%.2f (%.2f-%.2f)", exp(beta[i1]),
                       exp(beta[i1] - 1.96 * se[i1]), exp(beta[i1] + 1.96 * se[i1])),
    stroke_rr = sprintf("%.2f (%.2f-%.2f)", exp(b_stroke),
                        exp(b_stroke - 1.96 * sqrt(v_stroke)), exp(b_stroke + 1.96 * sqrt(v_stroke))),
    ratio_of_rr = sprintf("%.2f (%.2f-%.2f)", exp(beta[i2]),
                          exp(beta[i2] - 1.96 * se[i2]), exp(beta[i2] + 1.96 * se[i2])),
    p_interaction = 2 * pnorm(abs(beta[i2] / se[i2]), lower.tail = FALSE),
    row.names = NULL, stringsAsFactors = FALSE
  )
})
out <- do.call(rbind, rows)
out$p_interaction_text <- ifelse(out$p_interaction < 0.001, "<0.001", sprintf("%.3f", out$p_interaction))
print(out[, c("transition_group_label", "heart_rr", "stroke_rr", "ratio_of_rr", "p_interaction_text")],
      row.names = FALSE)

global <- data.frame(
  test = "Transition group x outcome type interaction (stacked modified Poisson, PSU cluster-robust)",
  chi_square = stat, df = dfree, p_value = p_global,
  participants = nrow(d), heart_events = sum(heart$y), stroke_events = sum(stroke$y),
  stringsAsFactors = FALSE
)
write.csv(out, file.path(table_dir, "bmc_supp_table_outcome_heterogeneity_by_group.csv"), row.names = FALSE)
write.csv(global, file.path(table_dir, "bmc_supp_table_outcome_heterogeneity_global.csv"), row.names = FALSE)
cat("\nWritten: bmc_supp_table_outcome_heterogeneity_{by_group,global}.csv\n")
