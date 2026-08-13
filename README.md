# CHARLS hypertension management transitions and incident CVD

This repository contains the reproducible analysis code and non-identifiable summary outputs for the manuscript:

**Hypertension management transitions and incident cardiovascular disease among middle-aged and older adults in China: a CHARLS landmark cohort study**

## Data availability

The China Health and Retirement Longitudinal Study (CHARLS) individual-level data are not redistributed here. CHARLS data are available to registered users through the CHARLS project website after application and approval:

<https://charls.pku.edu.cn/>

Do not commit raw CHARLS files, extracted `.dta` files, or participant-level derived datasets to this repository.

## Repository contents

- `01_build_sci_landmark_dataset.py`: builds the 2015 landmark cohort, hypertension management transition variables, sample-flow table, 5 x 5 transition matrix, missingness table, and event-year checks.
- `02_run_sci_statistical_analysis.R`: runs the original SCI statistical analysis, including logistic models, modified Poisson models, absolute risks, MICE, alternative thresholds, IPCW, subgroup analyses, and spline analysis.
- `05_bmc_rr_primary_reanalysis.R`: final BMC-oriented RR-primary analysis, including modified Poisson models, marginal risk estimates, clinically motivated contrasts, survey/IPCW/death scenarios, Western-medicine-only sensitivity analysis, and software-version outputs.
- `08_prior_cvd_history_sensitivity.R`, `10_outcome_heterogeneity_heart_stroke.R`, `15_ipcw_and_death_adjusted.R`: additional sensitivity and secondary analyses.
- `tables/`: non-identifiable aggregate tables and model outputs.
- `figures/`: final analysis figures.
- `reproducibility_logs/`: logs from the main analysis runs.

The files in `reproducibility_logs/` are unedited records of the original analysis runs. The participant-flow labels and starting counts were subsequently expanded in the manuscript figure to show the 17,705 to 13,965 linkage steps explicitly.

## Expected local data layout

After obtaining CHARLS data, create a local data directory outside version control and point the scripts to it with environment variables.

Recommended layout:

```text
raw/
  _extracted/
    2011/
      demographic_background.dta
      health_status_and_functioning.dta
      biomarkers.dta
    2013/
      Health_Status_and_Functioning.dta
    2015/
      Health_Status_and_Functioning.dta
      Biomarker.dta
    2018/
      Health_Status_and_Functioning.dta
      Sample_Infor.dta
    2020/
      Exit_Module.dta
  2011/
    PSU.zip
```

The 2011 `weight.dta` file should be extracted from the CHARLS 2011 `weight.rar` archive and placed in `_tmp_extract/weight.dta`, or in another directory supplied through `CHARLS_TMP_ROOT`.

## Running the analysis

Set paths for your local CHARLS files:

```bash
export CHARLS_RAW_ROOT="/path/to/raw/_extracted"
export CHARLS_DESIGN_ROOT="/path/to/raw"
export CHARLS_TMP_ROOT="/path/to/local/tmp_extract"
```

Then run the scripts in order:

```bash
python3 01_build_sci_landmark_dataset.py
Rscript 02_run_sci_statistical_analysis.R
Rscript 03_finalize_sci_reporting_package.R
Rscript 05_bmc_rr_primary_reanalysis.R
Rscript 08_prior_cvd_history_sensitivity.R
Rscript 10_outcome_heterogeneity_heart_stroke.R
Rscript 15_ipcw_and_death_adjusted.R
python3 13_redraw_figures.py
```

The core statistical reproducibility chain is `01`, `02`, `03`, `05`, `08`, `10`, and `15`. Script `13` regenerates the final figures from aggregate outputs.

## Software

The final software versions used for the submitted analysis are recorded in:

- `tables/bmc_supp_table_session_info.csv`
- `reproducibility_logs/`

## Privacy and licensing note

This repository is intended for analysis transparency only. It does not grant access to CHARLS data and does not change the CHARLS data-use terms. The analysis code is released under the MIT License; third-party CHARLS data remain subject to CHARLS data-use terms.
