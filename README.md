# Carbetocin sensitivity analyses (E-value, multiple imputation, pattern-mixture)

Reproducible analysis code for the sensitivity analyses reported in:

> Mao S-P, et al. *Carbetocin and peri-operative hemoglobin change at cesarean
> delivery in a setting of self-paid use: an 8-year retrospective cohort study.*
> BMC Pregnancy and Childbirth (under review).

This repository contains **only the analysis code**. It reproduces the three
sensitivity analyses added during peer review for the primary cesarean
peri-operative hemoglobin outcome, using the manuscript's fully adjusted model:

```
delta_Hb ~ carbetocin + age + BMI + parity + birth_weight + high_risk + C(calendar_year)
```

## Analyses

| Analysis | Purpose | Expected result |
|---|---|---|
| E-value (VanderWeele & Ding) | Robustness to unmeasured confounding | ~3.2 (point), ~2.8 (CI limit) |
| Multiple imputation (MICE, m=20) | Robustness to missing-at-random data | +0.75 g/dL (95% CI 0.62–0.88) |
| Delta / pattern-mixture (tipping point) | Robustness to not-missing-at-random data | ~2.0 g/dL |

The complete-case fully adjusted difference is +0.80 g/dL (95% CI 0.66–0.94),
reproducing the +0.79 g/dL reported in the manuscript.

## Reproducibility

A fixed random seed (`np.random.seed(2026)`) is set, so the multiple-imputation
and tipping-point results are fully reproducible.

## Usage

```bash
pip install -r requirements.txt
python carbetocin_R2_sensitivity_analyses.py
```

## Data availability

The clinical dataset is **not included** in this repository. It contains
identifiable clinical information extracted from the hospital electronic medical
record and is subject to institutional privacy and ethics restrictions. The data
are available from the corresponding author on reasonable request and with the
approval of the institutional review board. The script expects an Excel file
(`result_v5_with_BT_final.xlsx`) with sheets `cs` and `cs-high risk`; column
names used by the script are documented inline.

## License

MIT License (see `LICENSE`).

## Contact

Shih-Peng Mao, M.D., M.S. — Department of Obstetrics and Gynecology,
Taipei Medical University-Shuang Ho Hospital, New Taipei City, Taiwan.
ORCID: 0000-0003-4914-3179
