#!/usr/bin/env python3
"""Redraw Figure 1 (STROBE-style flow with exclusions) and Figure 2 (forest plot).

Figure 1 previously showed only the retained count at each step, which does not
meet STROBE item 13. Figure 2 was generated before the transition-group labels
were harmonised and carried four stale labels, no reference row, and no numeric
estimates.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle

ROOT = Path(__file__).resolve().parent
FIGDIR = ROOT / "figures"
MSDIR = ROOT / "manuscript"
TABLES = ROOT / "tables"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 9,
    "axes.linewidth": 0.8,
})

NAVY = "#1f3f66"

# --- Figure 1 ----------------------------------------------------------------
# (retained box text, n retained, exclusion text, n excluded)
FLOW = [
    ("2011 CHARLS demographic module", 17705, None, None),
    ("Linked with the 2011 health status module", 17594,
     "No 2011 health status module record", 111),
    ("Linked with the 2011 physical examination module\n"
     "containing measured blood pressure", 13965,
     "No 2011 physical examination,\nso blood pressure was not measured", 3629),
    ("Linked with the 2015 health status and\nphysical measurement modules", 10051,
     "Not interviewed or measured in 2015", 3914),
    ("Aged ≥45 years in 2011", 9850, "Aged <45 years in 2011", 201),
    ("Free of heart disease and stroke in 2011", 8590,
     "Heart disease or stroke reported in 2011", 1260),
    ("Complete 2011 and 2015 hypertension\nmanagement states", 8438,
     "Incomplete 2011 or 2015 blood pressure,\nawareness, or treatment information", 152),
    ("CVD-free at the 2015 landmark", 7680,
     "Heart disease or stroke reported by the\n2015 landmark", 758),
    ("Primary landmark cohort with an observed\n2018 CVD outcome", 6975,
     "No observed 2018 outcome\n(314 died, 386 absent from the 2018\nroster, 5 without outcome data)", 705),
]


def draw_figure1(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 10.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, len(FLOW) * 2.2)
    ax.axis("off")
    box_w, x_left = 4.6, 0.3
    y = len(FLOW) * 2.2 - 1.0
    for i, (label, n, excl, n_excl) in enumerate(FLOW):
        lines = label.count("\n") + 1
        h = 0.5 + 0.28 * lines
        ax.add_patch(Rectangle((x_left, y - h / 2), box_w, h,
                               fill=False, edgecolor=NAVY, linewidth=1.0))
        ax.text(x_left + box_w / 2, y, f"{label}\n(n={n:,})",
                ha="center", va="center", fontsize=8.5, linespacing=1.35)
        if i < len(FLOW) - 1:
            nxt = FLOW[i + 1]
            nxt_h = 0.5 + 0.28 * (nxt[0].count("\n") + 1)
            y_next = y - 2.2
            ax.add_patch(FancyArrowPatch((x_left + box_w / 2, y - h / 2),
                                         (x_left + box_w / 2, y_next + nxt_h / 2),
                                         arrowstyle="-|>", mutation_scale=11,
                                         color=NAVY, linewidth=0.9))
            e_lines = nxt[2].count("\n") + 1
            e_h = 0.45 + 0.26 * e_lines
            y_e = y - 1.1
            ax.add_patch(FancyArrowPatch((x_left + box_w / 2, y_e),
                                         (x_left + box_w + 0.55, y_e),
                                         arrowstyle="-|>", mutation_scale=10,
                                         color="0.35", linewidth=0.8))
            ax.add_patch(Rectangle((x_left + box_w + 0.55, y_e - e_h / 2), 4.3, e_h,
                                   fill=False, edgecolor="0.45", linewidth=0.8,
                                   linestyle=(0, (3, 2))))
            ax.text(x_left + box_w + 0.55 + 2.15, y_e,
                    f"Excluded (n={nxt[3]:,}):\n{nxt[2]}",
                    ha="center", va="center", fontsize=7.6, color="0.2", linespacing=1.3)
            y = y_next
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


# --- Figure 2 ----------------------------------------------------------------
ADJUSTED_COUNTS = {
    "Persistent non-hypertension": "170/3,477",
    "Incident hypertension": "71/915",
    "Persistent unawareness/no treatment": "45/631",
    "New awareness/treatment without control": "18/165",
    "Persistent awareness or treatment without control": "62/428",
    "Gained control": "28/283",
    "Persistent control": "15/164",
    "Loss of control": "17/190",
    "Other transitions": "29/454",
}


def parse_ci(text: str) -> tuple[float, float, float]:
    body = text.replace("–", "-").replace("(", " ").replace(")", " ")
    est, rng = body.split()[0], body.split()[1]
    lo, hi = rng.split("-")
    return float(est), float(lo), float(hi)


def draw_figure2(path: Path) -> None:
    src = pd.read_csv(TABLES / "bmc_main_table2_poisson_rr_three_models.csv")
    label_fix = lambda s: s.replace("but uncontrolled", "without control")
    rows = [("Persistent non-hypertension", 1.0, None, None, "1.00 (reference)")]
    for _, r in src.iterrows():
        text = str(r["Fully adjusted RR (95% CI)"])
        est, lo, hi = parse_ci(text)
        rows.append((label_fix(str(r.iloc[0])), est, lo, hi, text.replace("-", "–")))

    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    ys = list(range(len(rows)))[::-1]
    for y, (_, est, lo, hi, _) in zip(ys, rows):
        if lo is not None:
            ax.plot([lo, hi], [y, y], color=NAVY, linewidth=1.4, solid_capstyle="butt")
            ax.plot([lo, lo], [y - .12, y + .12], color=NAVY, linewidth=1.4)
            ax.plot([hi, hi], [y - .12, y + .12], color=NAVY, linewidth=1.4)
            ax.plot(est, y, "o", color=NAVY, markersize=6)
        else:
            ax.plot(est, y, "D", color=NAVY, markersize=7,
                    markerfacecolor="white", markeredgewidth=1.4)
    ax.axvline(1.0, color="0.4", linestyle="--", linewidth=0.9)
    ax.set_xscale("log")
    ax.set_xlim(0.7, 4.2)
    ax.set_xticks([0.8, 1.0, 1.5, 2.0, 3.0, 4.0])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.6)
    ax.set_xlabel("Fully adjusted risk ratio (95% CI), log scale", fontsize=9.5)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="0.9", linewidth=0.6)
    ax.set_axisbelow(True)

    trans = ax.get_yaxis_transform()
    # get_yaxis_transform() blends axes-fraction x with data-coordinate y, so the
    # header row must be placed just above the topmost data row.
    header_y = len(rows) - 0.55
    ax.text(1.03, header_y, "RR (95% CI)", transform=trans, fontsize=8.6,
            fontweight="bold", ha="left", va="center")
    ax.text(1.31, header_y, "Events/n", transform=trans, fontsize=8.6,
            fontweight="bold", ha="left", va="center")
    for y, (label, _, _, _, text) in zip(ys, rows):
        ax.text(1.03, y, text, transform=trans, fontsize=8.4, va="center")
        ax.text(1.31, y, ADJUSTED_COUNTS[label], transform=trans, fontsize=8.4, va="center")
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    f1 = FIGDIR / "figure1_participant_selection_flowchart.png"
    f2 = FIGDIR / "figure2_adjusted_rr_forest.png"
    for f in (f1, f2):
        if f.exists():
            shutil.copy2(f, f.with_name(f.stem + "_pre_editorial.png"))
    draw_figure1(f1)
    draw_figure2(f2)
    shutil.copy2(f1, MSDIR / "Figure_1.png")
    shutil.copy2(f2, MSDIR / "Figure_2.png")
    print(f"Redrawn: {f1.name}, {f2.name}")
    print(f"Copied to {MSDIR.name}/Figure_1.png and Figure_2.png")


if __name__ == "__main__":
    main()
