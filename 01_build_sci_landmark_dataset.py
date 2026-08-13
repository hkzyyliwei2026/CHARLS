"""
Analysis: CHARLS hypertension management transition landmark cohort for SCI submission
Date: 2026-08-11
Random seed: 42
Python: see runtime
Key packages: numpy, pandas, pyreadstat

This script builds the SCI-analysis landmark dataset from locally obtained
CHARLS files. Raw CHARLS files are read-only; derived individual-level data are
written to the local data/ directory, which is intentionally excluded from Git.
"""

from __future__ import annotations

import math
import os
import shutil
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadstat


np.random.seed(42)

SCRIPT_DIR = Path(__file__).resolve().parent
CHARLS_ROOT = Path(os.environ.get("CHARLS_ROOT", SCRIPT_DIR / "raw")).expanduser()
RAW_ROOT = Path(os.environ.get("CHARLS_RAW_ROOT", CHARLS_ROOT / "_extracted")).expanduser()
DESIGN_ROOT = Path(os.environ.get("CHARLS_DESIGN_ROOT", CHARLS_ROOT)).expanduser()
TMP_DIR = Path(os.environ.get("CHARLS_TMP_ROOT", SCRIPT_DIR / "_tmp_extract")).expanduser()
DATA_DIR = SCRIPT_DIR / "data"
TABLE_DIR = SCRIPT_DIR / "tables"

STATE_ORDER = [
    "normal",
    "unaware_htn",
    "aware_untreated",
    "treated_uncontrolled",
    "treated_controlled",
]

STATE_LABELS = {
    "normal": "Non-hypertensive",
    "unaware_htn": "Unaware hypertension",
    "aware_untreated": "Aware but untreated",
    "treated_uncontrolled": "Treated but uncontrolled",
    "treated_controlled": "Treated and controlled",
}

GROUP_ORDER = [
    "persistent_nonhypertensive",
    "incident_hypertension",
    "persistent_unaware_or_untreated",
    "new_awareness_or_treatment_uncontrolled",
    "persistent_aware_or_treated_uncontrolled",
    "gained_control",
    "persistent_controlled",
    "control_loss",
    "other_transition",
]

GROUP_LABELS = {
    "persistent_nonhypertensive": "Persistent non-hypertension",
    "incident_hypertension": "Incident hypertension",
    "persistent_unaware_or_untreated": "Persistent unawareness/no treatment",
    "new_awareness_or_treatment_uncontrolled": "New awareness/treatment without control",
    "persistent_aware_or_treated_uncontrolled": "Persistent aware/treated but uncontrolled",
    "gained_control": "Gained control",
    "persistent_controlled": "Persistent control",
    "control_loss": "Control loss",
    "other_transition": "Other transitions",
}


def read_dta(relative_path: str, usecols: list[str] | None = None) -> pd.DataFrame:
    df, _ = pyreadstat.read_dta(
        str(RAW_ROOT / relative_path),
        usecols=usecols,
        apply_value_formats=False,
    )
    if "ID" in df.columns:
        df["ID"] = df["ID"].astype(str)
    return df


def finite_number(value: object) -> bool:
    return isinstance(value, (int, float, np.integer, np.floating)) and math.isfinite(float(value))


def yes(value: object, code: int = 1) -> bool:
    return finite_number(value) and int(value) == code


def no(value: object, code: int = 2) -> bool:
    return finite_number(value) and int(value) == code


def valid_sbp(value: object) -> bool:
    return finite_number(value) and 50 <= float(value) <= 260


def valid_dbp(value: object) -> bool:
    return finite_number(value) and 30 <= float(value) <= 160


def mean_valid(row: pd.Series, columns: list[str], validator) -> float:
    values = [float(row[col]) for col in columns if validator(row[col])]
    return float(np.mean(values)) if values else np.nan


def convert_wave1_to_wave2_id(id_w1: str) -> str:
    return id_w1[:9] + "0" + id_w1[9:]


def add_bp_means(df: pd.DataFrame, suffix: str, threshold_sbp: int = 140, threshold_dbp: int = 90) -> pd.DataFrame:
    out = df.copy()
    out[f"sbp_{suffix}"] = out.apply(lambda row: mean_valid(row, ["qa003", "qa007", "qa011"], valid_sbp), axis=1)
    out[f"dbp_{suffix}"] = out.apply(lambda row: mean_valid(row, ["qa004", "qa008", "qa012"], valid_dbp), axis=1)
    out[f"bp_valid_{suffix}"] = out[f"sbp_{suffix}"].notna() & out[f"dbp_{suffix}"].notna()
    out[f"bp_high_{suffix}"] = (out[f"sbp_{suffix}"] >= threshold_sbp) | (out[f"dbp_{suffix}"] >= threshold_dbp)
    return out


def reconstructed_2015_status(row: pd.Series, disease_index: int) -> int | float:
    previous_yes = yes(row.get(f"zda007_{disease_index}_"))
    if disease_index == 1 and yes(row.get("zda008_1_")):
        previous_yes = True

    confirmation = row.get(f"da007_w2_1_{disease_index}_")
    new_diagnosis = row.get(f"da007_w2_2_{disease_index}_")
    raw_status = row.get(f"da007_{disease_index}_")

    if finite_number(confirmation):
        if previous_yes:
            return int(yes(confirmation) or (no(confirmation) and yes(new_diagnosis)))
        return int(no(confirmation) or (yes(confirmation) and yes(new_diagnosis)))
    if finite_number(raw_status):
        return int(yes(raw_status))
    return np.nan


def add_hypertension_state(df: pd.DataFrame, suffix: str, diagnosis_col: str, awareness_col: str | None) -> pd.DataFrame:
    out = df.copy()
    treated = out["da011s1"].apply(lambda value: yes(value, 1)) | out["da011s2"].apply(lambda value: yes(value, 2))
    aware = out[diagnosis_col].apply(lambda value: yes(value, 1))
    if awareness_col and awareness_col in out.columns:
        aware = aware | out[awareness_col].apply(lambda value: yes(value, 1))
    aware = aware | treated
    out[f"htn_treated_{suffix}"] = treated.astype(int)
    out[f"htn_aware_{suffix}"] = aware.astype(int)

    def classify(row: pd.Series) -> str | float:
        if not bool(row[f"bp_valid_{suffix}"]):
            return np.nan
        high = bool(row[f"bp_high_{suffix}"])
        is_aware = bool(row[f"htn_aware_{suffix}"])
        is_treated = bool(row[f"htn_treated_{suffix}"])
        if not is_aware and not is_treated and not high:
            return "normal"
        if not is_aware and not is_treated and high:
            return "unaware_htn"
        if is_aware and not is_treated:
            return "aware_untreated"
        if is_treated and high:
            return "treated_uncontrolled"
        return "treated_controlled"

    out[f"htn_state_{suffix}"] = out.apply(classify, axis=1)
    return out


def yes_code(value: object, code: int = 1) -> float:
    if not finite_number(value):
        return np.nan
    return float(int(value) == code)


def education_group(value: object) -> str | float:
    if not finite_number(value):
        return np.nan
    code = int(value)
    if code == 1:
        return "no_formal"
    if code in {2, 3, 4}:
        return "primary_or_below"
    if code == 5:
        return "middle_school"
    if code >= 6:
        return "high_school_plus"
    return np.nan


def marital_group(value: object) -> str | float:
    if not finite_number(value):
        return np.nan
    return "married_partnered" if int(value) in {1, 2} else "not_married_partnered"


def hukou_group(value: object) -> str | float:
    if not finite_number(value):
        return np.nan
    code = int(value)
    if code == 1:
        return "agricultural"
    if code == 2:
        return "non_agricultural"
    if code in {3, 4}:
        return "other"
    return np.nan


def sex_group(value: object) -> str | float:
    if not finite_number(value):
        return np.nan
    code = int(value)
    if code == 1:
        return "male"
    if code == 2:
        return "female"
    return np.nan


def drink_group(value: object) -> float:
    if not finite_number(value):
        return np.nan
    return float(int(value) in {1, 2})


def valid_height_cm(value: object) -> bool:
    return finite_number(value) and 100 <= float(value) <= 220


def valid_weight_kg(value: object) -> bool:
    return finite_number(value) and 25 <= float(value) <= 250


def transition_group(row: pd.Series, start_col: str = "htn_state_2011", end_col: str = "htn_state_2015") -> str:
    start = row[start_col]
    end = row[end_col]
    uncontrolled = {"unaware_htn", "aware_untreated", "treated_uncontrolled"}
    untreated_or_unaware = {"unaware_htn", "aware_untreated"}
    aware_or_treated_uncontrolled = {"aware_untreated", "treated_uncontrolled"}

    if start == "normal" and end == "normal":
        return "persistent_nonhypertensive"
    if start == "normal" and end != "normal":
        return "incident_hypertension"
    if start in untreated_or_unaware and end in untreated_or_unaware:
        return "persistent_unaware_or_untreated"
    if start == "unaware_htn" and end in aware_or_treated_uncontrolled:
        return "new_awareness_or_treatment_uncontrolled"
    if start in aware_or_treated_uncontrolled and end in aware_or_treated_uncontrolled:
        return "persistent_aware_or_treated_uncontrolled"
    if start in uncontrolled and end == "treated_controlled":
        return "gained_control"
    if start == "treated_controlled" and end == "treated_controlled":
        return "persistent_controlled"
    if start == "treated_controlled" and end in uncontrolled:
        return "control_loss"
    return "other_transition"


def classify_state_from_existing(row: pd.Series, suffix: str, sbp_threshold: int, dbp_threshold: int) -> str | float:
    if not bool(row[f"bp_valid_{suffix}"]):
        return np.nan
    high = bool((row[f"sbp_{suffix}"] >= sbp_threshold) or (row[f"dbp_{suffix}"] >= dbp_threshold))
    aware = bool(row[f"htn_aware_{suffix}"] == 1)
    treated = bool(row[f"htn_treated_{suffix}"] == 1)
    if not aware and not treated and not high:
        return "normal"
    if not aware and not treated and high:
        return "unaware_htn"
    if aware and not treated:
        return "aware_untreated"
    if treated and high:
        return "treated_uncontrolled"
    return "treated_controlled"


def extract_design_files() -> None:
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    psu_path = TMP_DIR / "PSU.dta"
    if not psu_path.exists():
        psu_zip = DESIGN_ROOT / "2011" / "PSU.zip"
        if not psu_zip.exists():
            raise FileNotFoundError(
                f"Missing {psu_zip}. Set CHARLS_DESIGN_ROOT to the directory containing the original 2011/PSU.zip file."
            )
        with zipfile.ZipFile(psu_zip) as zf:
            zf.extract("PSU.dta", TMP_DIR)
    weight_path = TMP_DIR / "weight.dta"
    if not weight_path.exists():
        weight_source = DESIGN_ROOT / "2011" / "weight.dta"
        if weight_source.exists():
            shutil.copy2(weight_source, weight_path)
        else:
            raise FileNotFoundError(
                f"Missing {weight_path}. Extract CHARLS 2011 weight.rar and place weight.dta in {TMP_DIR}, "
                "or set CHARLS_TMP_ROOT to a directory containing weight.dta."
            )


def add_design_variables(df: pd.DataFrame) -> pd.DataFrame:
    extract_design_files()
    psu, _ = pyreadstat.read_dta(str(TMP_DIR / "PSU.dta"), apply_value_formats=False)
    psu["communityID"] = psu["communityID"].astype(str)
    weight11, _ = pyreadstat.read_dta(str(TMP_DIR / "weight.dta"), apply_value_formats=False)
    weight11["ID_w1"] = weight11["ID"].astype(str)
    weight11 = weight11[["ID_w1", "communityID", "householdID", "bio_weight2", "ind_weight_ad2"]]
    weight11["communityID"] = weight11["communityID"].astype(str)
    weight11["householdID"] = weight11["householdID"].astype(str)
    out = df.merge(weight11, on="ID_w1", how="left").merge(psu[["communityID", "province_eng", "city_eng", "urban_nbs", "areatype"]], on="communityID", how="left")
    out["cluster_psu"] = out["communityID"]
    out["strata_proxy"] = out["province_eng"].astype(str) + "_" + out["urban_nbs"].astype(str)
    out["baseline_biomarker_weight"] = out["bio_weight2"]
    return out


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    demo11 = read_dta("2011/demographic_background.dta", ["ID", "ba002_1", "rgender", "bd001", "be001", "bc001"])
    health11 = read_dta(
        "2011/health_status_and_functioning.dta",
        ["ID", "da007_1_", "da008_1_", "da007_2_", "da007_3_", "da007_7_", "da007_8_", "da011s1", "da011s2", "da059", "da067"],
    )
    bio11 = read_dta("2011/biomarkers.dta", ["ID", "qa003", "qa004", "qa007", "qa008", "qa011", "qa012", "qi002", "ql002"])

    w1_demo_health = demo11.merge(health11, on="ID", how="inner")
    w1 = w1_demo_health.merge(bio11, on="ID", how="inner")
    w1["ID_w1"] = w1["ID"]
    w1["ID"] = w1["ID_w1"].apply(convert_wave1_to_wave2_id)
    w1["age_2011"] = 2011 - pd.to_numeric(w1["ba002_1"], errors="coerce")
    w1["baseline_cvd_2011"] = w1["da007_7_"].apply(lambda value: yes(value, 1)) | w1["da007_8_"].apply(lambda value: yes(value, 1))
    w1 = add_bp_means(w1, "2011")
    w1 = add_hypertension_state(w1, "2011", "da007_1_", "da008_1_")
    w1["sex"] = w1["rgender"].apply(sex_group)
    w1["education"] = w1["bd001"].apply(education_group)
    w1["marital"] = w1["be001"].apply(marital_group)
    w1["hukou"] = w1["bc001"].apply(hukou_group)
    w1["current_smoking"] = w1["da059"].apply(lambda value: yes_code(value, 1))
    w1["alcohol_last_year"] = w1["da067"].apply(drink_group)
    w1["diabetes"] = w1["da007_3_"].apply(lambda value: yes_code(value, 1))
    w1["dyslipidemia"] = w1["da007_2_"].apply(lambda value: yes_code(value, 1))
    w1["height_m"] = w1["qi002"].where(w1["qi002"].apply(valid_height_cm)) / 100
    w1["weight_kg"] = w1["ql002"].where(w1["ql002"].apply(valid_weight_kg))
    w1["bmi_2011"] = w1["weight_kg"] / w1["height_m"] ** 2

    w1_keep = w1[
        [
            "ID",
            "ID_w1",
            "age_2011",
            "sex",
            "education",
            "marital",
            "hukou",
            "current_smoking",
            "alcohol_last_year",
            "diabetes",
            "dyslipidemia",
            "bmi_2011",
            "baseline_cvd_2011",
            "sbp_2011",
            "dbp_2011",
            "bp_valid_2011",
            "bp_high_2011",
            "htn_aware_2011",
            "htn_treated_2011",
            "htn_state_2011",
        ]
    ]

    health15 = read_dta(
        "2015/Health_Status_and_Functioning.dta",
        [
            "ID",
            "zda007_1_",
            "zda007_7_",
            "zda007_8_",
            "zda008_1_",
            "da007_1_",
            "da007_w2_1_1_",
            "da007_w2_2_1_",
            "da007_w2_1_7_",
            "da007_w2_2_7_",
            "da007_w2_1_8_",
            "da007_w2_2_8_",
            "da011s1",
            "da011s2",
            "da059",
            "da067",
            "da007_2_",
            "da007_3_",
        ],
    )
    bio15 = read_dta("2015/Biomarker.dta", ["ID", "qa003", "qa004", "qa007", "qa008", "qa011", "qa012", "qi002", "ql002"])
    health15["htn_diag_2015"] = health15.apply(lambda row: reconstructed_2015_status(row, 1), axis=1)
    health15["heart_2015"] = health15.apply(lambda row: reconstructed_2015_status(row, 7), axis=1)
    health15["stroke_2015"] = health15.apply(lambda row: reconstructed_2015_status(row, 8), axis=1)
    w3 = health15.merge(bio15, on="ID", how="inner")
    w3["cvd_2015"] = w3["heart_2015"].apply(lambda value: yes(value, 1)) | w3["stroke_2015"].apply(lambda value: yes(value, 1))
    w3 = add_bp_means(w3, "2015")
    w3 = add_hypertension_state(w3, "2015", "htn_diag_2015", "zda008_1_")
    w3["current_smoking_2015"] = w3["da059"].apply(lambda value: yes_code(value, 1))
    w3["alcohol_last_year_2015"] = w3["da067"].apply(drink_group)
    w3["dyslipidemia_2015"] = w3["da007_2_"].apply(lambda value: yes_code(value, 1))
    w3["height_m_2015"] = w3["qi002"].where(w3["qi002"].apply(valid_height_cm)) / 100
    w3["weight_kg_2015"] = w3["ql002"].where(w3["ql002"].apply(valid_weight_kg))
    w3["bmi_2015"] = w3["weight_kg_2015"] / w3["height_m_2015"] ** 2

    w3_keep = w3[
        [
            "ID",
            "heart_2015",
            "stroke_2015",
            "cvd_2015",
            "sbp_2015",
            "dbp_2015",
            "bp_valid_2015",
            "bp_high_2015",
            "htn_aware_2015",
            "htn_treated_2015",
            "htn_state_2015",
            "current_smoking_2015",
            "alcohol_last_year_2015",
            "dyslipidemia_2015",
            "bmi_2015",
        ]
    ]

    follow18 = read_dta(
        "2018/Health_Status_and_Functioning.dta",
        ["ID", "da007_7_", "da007_8_", "da009_1_7_", "da009_1_8_", "da007_w2_5", "da019_w2_1"],
    )
    follow = follow18.copy()
    follow["heart_year_2018_report"] = pd.to_numeric(follow["da009_1_7_"], errors="coerce")
    follow["stroke_year_2018_report"] = pd.to_numeric(follow["da009_1_8_"], errors="coerce")
    follow["incident_heart_2015_2018"] = (
        follow["da007_w2_5"].apply(lambda value: yes(value, 1))
        | (follow["da007_7_"].apply(lambda value: yes(value, 1)) & (follow["heart_year_2018_report"] > 2015))
    ).astype(int)
    follow["incident_stroke_2015_2018"] = (
        follow["da019_w2_1"].apply(lambda value: yes(value, 1))
        | (follow["da007_8_"].apply(lambda value: yes(value, 1)) & (follow["stroke_year_2018_report"] > 2015))
    ).astype(int)
    follow["incident_cvd_2015_2018"] = ((follow["incident_heart_2015_2018"] == 1) | (follow["incident_stroke_2015_2018"] == 1)).astype(int)

    heart_known_strict = follow["heart_year_2018_report"].between(2016, 2018)
    stroke_known_strict = follow["stroke_year_2018_report"].between(2016, 2018)
    heart_ambiguous = (follow["incident_heart_2015_2018"] == 1) & ~heart_known_strict
    stroke_ambiguous = (follow["incident_stroke_2015_2018"] == 1) & ~stroke_known_strict
    follow["strict_cvd_2016_2018"] = np.where(
        heart_ambiguous | stroke_ambiguous,
        np.nan,
        (heart_known_strict | stroke_known_strict).astype(int),
    )
    event_years = follow[["heart_year_2018_report", "stroke_year_2018_report"]]
    valid_event_years = event_years.where(event_years.apply(lambda col: col.between(2010, 2025)))
    follow["cvd_event_year_min"] = valid_event_years.min(axis=1)

    follow_keep = follow[
        [
            "ID",
            "incident_heart_2015_2018",
            "incident_stroke_2015_2018",
            "incident_cvd_2015_2018",
            "strict_cvd_2016_2018",
            "heart_year_2018_report",
            "stroke_year_2018_report",
            "cvd_event_year_min",
        ]
    ]

    linked = w1_keep.merge(w3_keep, on="ID", how="inner").merge(follow_keep, on="ID", how="left")
    linked = add_design_variables(linked)
    linked["age_eligible_2011"] = linked["age_2011"] >= 45
    linked["complete_transition"] = linked["htn_state_2011"].notna() & linked["htn_state_2015"].notna()
    linked["observed_2018_outcome"] = linked["incident_cvd_2015_2018"].notna()
    linked["cvd_free_landmark_2015"] = linked["cvd_2015"] == 0
    linked["landmark_eligible_before_followup"] = (
        linked["age_eligible_2011"] & ~linked["baseline_cvd_2011"] & linked["complete_transition"] & linked["cvd_free_landmark_2015"]
    )
    linked["analysis_eligible"] = linked["landmark_eligible_before_followup"] & linked["observed_2018_outcome"]

    landmark_all = linked[linked["landmark_eligible_before_followup"]].copy()
    analysis = linked[linked["analysis_eligible"]].copy()
    for df in (landmark_all, analysis):
        df["transition_group"] = df.apply(transition_group, axis=1)
        df["transition_group_label"] = df["transition_group"].map(GROUP_LABELS)
        df["age_per10"] = df["age_2011"] / 10
        df["age_group"] = np.where(df["age_2011"] < 60, "45-59", ">=60")
        df["mean_sbp_2011_2015"] = (df["sbp_2011"] + df["sbp_2015"]) / 2
        df["mean_dbp_2011_2015"] = (df["dbp_2011"] + df["dbp_2015"]) / 2
        df["htn_state_2011_130_80"] = df.apply(lambda row: classify_state_from_existing(row, "2011", 130, 80), axis=1)
        df["htn_state_2015_130_80"] = df.apply(lambda row: classify_state_from_existing(row, "2015", 130, 80), axis=1)
        df["transition_group_130_80"] = df.apply(lambda row: transition_group(row, "htn_state_2011_130_80", "htn_state_2015_130_80"), axis=1)

    flow = pd.DataFrame(
        [
            {"step": "2011 demographic module", "n": len(demo11)},
            {"step": "Linked with 2011 health status module", "n": len(w1_demo_health)},
            {"step": "Linked with 2011 physical measurement module (biomarkers), including measured blood pressure", "n": len(w1_keep)},
            {"step": "Linked with 2015 health status and physical measurement modules", "n": len(linked)},
            {"step": "Age >=45 years in 2011", "n": int(linked["age_eligible_2011"].sum())},
            {"step": "No heart disease or stroke in 2011", "n": int((linked["age_eligible_2011"] & ~linked["baseline_cvd_2011"]).sum())},
            {"step": "Complete 2011 and 2015 hypertension management states", "n": int((linked["age_eligible_2011"] & ~linked["baseline_cvd_2011"] & linked["complete_transition"]).sum())},
            {"step": "CVD-free at the 2015 landmark", "n": len(landmark_all)},
            {"step": "Observed 2018 CVD outcome", "n": len(analysis)},
        ]
    )

    matrix = (
        analysis.groupby(["htn_state_2011", "htn_state_2015"], dropna=False)
        .agg(
            n=("ID", "size"),
            cvd_events=("incident_cvd_2015_2018", "sum"),
            heart_events=("incident_heart_2015_2018", "sum"),
            stroke_events=("incident_stroke_2015_2018", "sum"),
            event_rate_pct=("incident_cvd_2015_2018", lambda x: x.mean() * 100),
            mean_sbp_2011=("sbp_2011", "mean"),
            mean_dbp_2011=("dbp_2011", "mean"),
            mean_sbp_2015=("sbp_2015", "mean"),
            mean_dbp_2015=("dbp_2015", "mean"),
        )
        .reset_index()
    )
    matrix["state_2011_label"] = matrix["htn_state_2011"].map(STATE_LABELS)
    matrix["state_2015_label"] = matrix["htn_state_2015"].map(STATE_LABELS)

    counts_matrix = pd.crosstab(analysis["htn_state_2011"], analysis["htn_state_2015"]).reindex(index=STATE_ORDER, columns=STATE_ORDER, fill_value=0)
    group_summary = (
        analysis.groupby("transition_group")
        .agg(
            n=("ID", "size"),
            cvd_events=("incident_cvd_2015_2018", "sum"),
            heart_events=("incident_heart_2015_2018", "sum"),
            stroke_events=("incident_stroke_2015_2018", "sum"),
            event_rate_pct=("incident_cvd_2015_2018", lambda x: x.mean() * 100),
            mean_sbp_2011=("sbp_2011", "mean"),
            mean_dbp_2011=("dbp_2011", "mean"),
            mean_sbp_2015=("sbp_2015", "mean"),
            mean_dbp_2015=("dbp_2015", "mean"),
        )
        .reindex(GROUP_ORDER)
        .reset_index()
    )
    group_summary["transition_group_label"] = group_summary["transition_group"].map(GROUP_LABELS)
    group_summary["sparse_event_flag"] = group_summary["cvd_events"] < 15

    event_year_quality = pd.DataFrame(
        [
            {"metric": "incident_cvd_events", "n": int(analysis["incident_cvd_2015_2018"].sum())},
            {"metric": "events_with_any_reported_diagnosis_year", "n": int((analysis.loc[analysis["incident_cvd_2015_2018"] == 1, "cvd_event_year_min"].notna()).sum())},
            {"metric": "events_with_year_2015", "n": int((analysis.loc[analysis["incident_cvd_2015_2018"] == 1, "cvd_event_year_min"] == 2015).sum())},
            {"metric": "events_with_year_2016", "n": int((analysis.loc[analysis["incident_cvd_2015_2018"] == 1, "cvd_event_year_min"] == 2016).sum())},
            {"metric": "events_with_year_2017", "n": int((analysis.loc[analysis["incident_cvd_2015_2018"] == 1, "cvd_event_year_min"] == 2017).sum())},
            {"metric": "events_with_year_2018", "n": int((analysis.loc[analysis["incident_cvd_2015_2018"] == 1, "cvd_event_year_min"] == 2018).sum())},
            {"metric": "events_missing_valid_year_or_since-last-interview-only", "n": int(((analysis["incident_cvd_2015_2018"] == 1) & analysis["strict_cvd_2016_2018"].isna()).sum())},
            {"metric": "events_with_year_before_2015", "n": int((analysis.loc[analysis["incident_cvd_2015_2018"] == 1, "cvd_event_year_min"] < 2015).sum())},
            {"metric": "events_with_year_after_2018", "n": int((analysis.loc[analysis["incident_cvd_2015_2018"] == 1, "cvd_event_year_min"] > 2018).sum())},
        ]
    )

    missing_vars = [
        "age_2011",
        "sex",
        "hukou",
        "education",
        "marital",
        "bmi_2011",
        "current_smoking",
        "alcohol_last_year",
        "diabetes",
        "dyslipidemia",
        "sbp_2011",
        "dbp_2011",
        "sbp_2015",
        "dbp_2015",
        "baseline_biomarker_weight",
        "cluster_psu",
    ]
    missing = pd.DataFrame(
        {
            "variable": missing_vars,
            "n_total": len(analysis),
            "n_missing": [int(analysis[v].isna().sum()) for v in missing_vars],
        }
    )
    missing["missing_pct"] = missing["n_missing"] / missing["n_total"] * 100

    DATA_DIR.joinpath("sci_landmark_all_before_2018_followup.csv").write_text("")
    landmark_all.to_csv(DATA_DIR / "sci_landmark_all_before_2018_followup.csv", index=False)
    analysis.to_csv(DATA_DIR / "sci_landmark_analysis_dataset.csv", index=False)
    flow.to_csv(TABLE_DIR / "sci_sample_flow.csv", index=False)
    matrix.to_csv(TABLE_DIR / "supp_table_5x5_transition_matrix_details.csv", index=False)
    counts_matrix.to_csv(TABLE_DIR / "supp_table_5x5_transition_matrix_counts.csv")
    group_summary.to_csv(TABLE_DIR / "table_transition_group_event_summary.csv", index=False)
    event_year_quality.to_csv(TABLE_DIR / "supp_table_event_year_quality.csv", index=False)
    missing.to_csv(TABLE_DIR / "supp_table_missingness.csv", index=False)

    print("SCI landmark dataset built")
    print(flow.to_string(index=False))
    print("\nMain transition groups")
    print(group_summary[["transition_group_label", "n", "cvd_events", "event_rate_pct", "sparse_event_flag"]].to_string(index=False))
    print("\nEvent year quality")
    print(event_year_quality.to_string(index=False))
    print("\nMissingness")
    print(missing.to_string(index=False))


if __name__ == "__main__":
    main()
