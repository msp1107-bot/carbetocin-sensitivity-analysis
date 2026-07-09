"""
carbetocin_R2_sensitivity_analyses.py
-------------------------------------
Reproduces the three sensitivity analyses added in Revision R2 for the cesarean
peri-operative hemoglobin outcome, using the manuscript's FULLY ADJUSTED model:

    delta_Hb ~ carbetocin + age + BMI + parity + birth_weight + high_risk + C(calendar_year)

Outputs (expected, on result_v5_with_BT_final.xlsx):
    - Fully adjusted difference  : +0.80 g/dL (95% CI ~0.66-0.94)  [reproduces reported +0.79]
    - R8.1 E-value               : ~3.2 (point), ~2.8 (CI limit closest to null)
    - R8.2 MICE pooled difference : +0.75 g/dL (95% CI 0.62-0.88)
    - R8.2 MNAR tipping point     : ~2.0 g/dL

Requirements:  pip install pandas numpy scipy statsmodels openpyxl
Data:          result_v5_with_BT_final.xlsx  (sheets 'cs' and 'cs-high risk')
Author:        Revision R2
"""

import re, warnings
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.imputation.mice import MICEData, MICE
warnings.filterwarnings("ignore")
np.random.seed(2026)   # fixed seed: MICE/tipping fully reproducible

DATA = "result_v5_with_BT_final.xlsx"
COVARS = "age + BMI + parity + bw + highrisk + C(year)"   # fully adjusted model
FORMULA = f"delta ~ carb + {COVARS}"

# ---------------------------------------------------------------- load & build
def gest_weeks(x):
    """(not used in the primary model; kept for completeness)"""
    if pd.isna(x):
        return np.nan
    m = re.search(r"(\d+)\D+(\d+)", str(x))
    return int(m.group(1)) + int(m.group(2)) / 7 if m else np.nan

xl = pd.ExcelFile(DATA)
cs = pd.concat(
    [xl.parse("cs").assign(highrisk=0), xl.parse("cs-high risk").assign(highrisk=1)],
    ignore_index=True,
)
cs["carb"]  = cs["Duratocin"].map(lambda x: 1 if str(x).strip() == "+"
                                  else (0 if str(x).strip() == "-" else np.nan))
cs["delta"] = pd.to_numeric(cs["delta_Hgb"], errors="coerce")           # post - pre Hb
cs = cs.rename(columns={"年齡": "age", "產次": "parity", "胎兒出生體重GM": "bw"})
for c in ["age", "BMI", "parity", "bw"]:
    cs[c] = pd.to_numeric(cs[c], errors="coerce")
cs["year"] = pd.to_datetime(cs["嬰兒娩出"], errors="coerce").dt.year     # calendar year

d = cs[cs["carb"].notna()].copy()
print(f"Cesarean cohort with known arm: n={len(d)} "
      f"(carbetocin {int((d.carb==1).sum())}, control {int((d.carb==0).sum())})")

# ============================================================ FULLY ADJUSTED FIT
cc = d.dropna(subset=["delta", "age", "BMI", "parity", "bw", "highrisk", "year"]).copy()
m = smf.ols(FORMULA, data=cc).fit()
est, (lo, hi) = m.params["carb"], m.conf_int().loc["carb"]
sd_raw, sd_res = cc["delta"].std(), np.sqrt(m.scale)
print(f"\n[FULLY ADJUSTED] complete-case n={len(cc)}: "
      f"diff = {est:+.3f} g/dL  95% CI {lo:+.3f} to {hi:+.3f}")

# ============================================================ R8.1  E-VALUE
def evalue(mean_diff, sd):
    d_std = abs(mean_diff) / sd                 # standardized mean difference (Cohen's d)
    rr = np.exp(0.91 * d_std)                   # VanderWeele-Ding continuous-outcome approx.
    rr = max(rr, 1 / rr)
    return rr, rr + np.sqrt(rr * (rr - 1))

print("\n=== R8.1  E-VALUE ===")
for lab, sd in [("SD of outcome (raw)", sd_raw), ("residual SD", sd_res)]:
    rr_pt, e_pt = evalue(est, sd)
    _,     e_ci = evalue(min(abs(lo), abs(hi)), sd)   # CI limit closest to null
    print(f"  [{lab} = {sd:.3f}]  d={abs(est)/sd:.3f}  RR={rr_pt:.2f}  "
          f"E-value(point)={e_pt:.2f}  E-value(CI limit)={e_ci:.2f}")

# ============================================================ R8.2  MISSINGNESS
print("\n=== R8.2  MISSING-DATA MECHANISM ===")
d["miss"] = d["delta"].isna().astype(int)
a, b = d.loc[d.carb == 1, "miss"], d.loc[d.carb == 0, "miss"]
orr, pf = stats.fisher_exact([[(a == 0).sum(), (a == 1).sum()],
                              [(b == 0).sum(), (b == 1).sum()]])
print(f"  availability: carbetocin {(1-a.mean())*100:.1f}% vs control {(1-b.mean())*100:.1f}%  "
      f"(Fisher OR={orr:.2f}, p={pf:.2e})")

dl = d.dropna(subset=["age", "BMI", "parity", "bw", "highrisk", "year"]).copy()
ml = smf.logit("miss ~ carb + age + BMI + parity + bw + highrisk + C(year)",
               data=dl).fit(disp=0)
ll0 = smf.logit("miss ~ 1", data=dl).fit(disp=0).llf
lr = 2 * (ml.llf - ll0)
p_global = stats.chi2.sf(lr, len(ml.params) - 1)
print(f"  global LR test that missingness depends on observed data: "
      f"chi2={lr:.1f}, p={p_global:.2e}  => reject MCAR")
print("  predictors of missingness (p<0.05):")
for v in ml.params.index:
    if v != "Intercept" and ml.pvalues[v] < 0.05:
        print(f"    {v:14s} OR={np.exp(ml.params[v]):.3f}  p={ml.pvalues[v]:.4f}")

# ------------------------------------------------ multiple imputation (MICE, m=20)
print("\n=== R8.2  MULTIPLE IMPUTATION (MICE, m=20, MAR) ===")
mi_cols = ["delta", "carb", "age", "BMI", "parity", "bw", "highrisk", "year"]
mi = MICEData(d[mi_cols].copy())
fit = MICE(FORMULA, sm.OLS, mi).fit(n_burnin=10, n_imputations=20)
row = fit.summary().tables[1].loc["carb"]
print(f"  pooled adjusted difference = {float(row['Coef.']):+.3f} g/dL  "
      f"95% CI {float(row['[0.025']):+.3f} to {float(row['0.975]']):+.3f}")

# ------------------------------------------------ MNAR delta / tipping point
print("\n=== R8.2  MNAR DELTA / TIPPING POINT ===")
miss_mask = d["delta"].isna().values
prev = None
for shift in np.arange(0.0, 4.01, 0.25):
    ests = []
    for _ in range(8):                       # a few imputations per shift
        mi2 = MICEData(d[mi_cols].copy())
        for _ in range(8):
            mi2.update_all()
        dd = mi2.data.copy()
        dd.loc[miss_mask & (dd["carb"] == 1), "delta"] -= shift   # push carb-arm missing down
        ests.append(smf.ols(FORMULA, data=dd).fit().params["carb"])
    mean_est = np.mean(ests)
    if prev is not None and prev > 0 and mean_est <= 0:
        print(f"  tipping point ~ {shift:.2f} g/dL  (adjusted difference {mean_est:+.3f})")
        break
    prev = mean_est

print("\nDone. Compare the printed values with the manuscript "
      "(E-value ~3.2/2.8; MICE +0.75 [0.62-0.88]; tipping point ~2.0 g/dL).")
