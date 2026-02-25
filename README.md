# 📊 Automated Statistical Analysis Pipeline for R

> **Drop in any dataset. Get a complete, publication-ready statistical analysis — automatically.**

A single self-contained R script that handles the full statistical workflow for multi-group comparison studies: from raw data to descriptive stats, normality testing, parametric and non-parametric group tests, post-hoc comparisons, effect sizes, power analysis, correlations, regression, and fully formatted exports — with zero manual steps between them.

---

## ✨ Features

| Category | What's included |
|---|---|
| **Descriptive Statistics** | Mean, SD, SE, median, IQR, Q1/Q3, min/max, range, CV%, skewness, kurtosis, 95% CI — per group |
| **Normality Testing** | Shapiro-Wilk · Kolmogorov-Smirnov (Lilliefors) · Anderson-Darling · Jarque-Bera + QQ plots |
| **Group Comparisons** | Auto-selects the right test: Student/Welch t-test, Mann-Whitney U, One-Way ANOVA, Welch ANOVA, Kruskal-Wallis |
| **Variance Homogeneity** | Levene's test + Bartlett's test — informs parametric test choice |
| **Post-hoc Tests** | Tukey HSD · Games-Howell · Dunn · Pairwise Wilcoxon — all run; best one flagged automatically |
| **Effect Sizes** | Cohen's d · Hedges' g · Glass's Δ · Rank-biserial r · η² · ω² · ε² — with interpretation labels |
| **Power Analysis** | Observed power + minimum n for 80% and 90% power per outcome |
| **Correlation** | Pearson + Spearman with p-values; corrplots + GGally scatter matrix |
| **Regression** | Simple and full OLS per outcome — coefficients, 95% CI, model fit, diagnostic plots |
| **Missing Data** | Per-variable count/% + visual pattern chart (naniar) |
| **Outlier Detection** | Z-score flagging at configurable threshold |
| **ANCOVA** | *(optional)* Covariate-adjusted group means via emmeans |
| **Repeated Measures** | *(optional)* LME4 mixed models with ICC and marginal/conditional R² |
| **Survival Analysis** | *(optional)* Kaplan-Meier curves · log-rank test · Cox proportional hazards |
| **Exports** | Multi-sheet Excel workbook · CSV tables · 300 dpi PNG plots · Word report · RDS cache |

---

## 🚀 Quick Start

```r
# 1. Clone or download statistical_pipeline.R

# 2. Run with built-in example (no data needed)
source("statistical_pipeline.R")

# That's it. Outputs land in pipeline_output/
```

To use your own data, edit the `CONFIG` block at the top of the script before sourcing:

```r
CONFIG$data_source  <- "csv"
CONFIG$data_path    <- "your_data.csv"
CONFIG$group_var    <- "treatment_group"
CONFIG$outcome_vars <- c("score", "weight", "blood_pressure")

source("statistical_pipeline.R")
```

---

## ⚙️ Configuration

All settings live in a single `CONFIG` list. No digging through code.

```r
CONFIG <- list(

  # --- DATA SOURCE ---
  data_source    = "builtin_example",  # "csv" | "excel" | "rds" | "builtin_example"
  data_path      = "",

  # --- VARIABLE ROLES ---
  group_var      = "group",            # categorical grouping column
  outcome_vars   = NULL,               # NULL = auto-detect all numeric columns
  covariate_vars = NULL,               # for ANCOVA (optional)
  id_var         = NULL,               # subject ID for repeated measures
  time_var       = NULL,               # time column for survival / longitudinal
  event_var      = NULL,               # 0/1 event indicator for survival

  # --- TOGGLE ANALYSES ---
  run_descriptives     = TRUE,
  run_normality        = TRUE,
  run_group_tests      = TRUE,
  run_posthoc          = TRUE,
  run_effect_sizes     = TRUE,
  run_power            = TRUE,
  run_correlation      = TRUE,
  run_regression       = TRUE,
  run_ancova           = FALSE,
  run_survival         = FALSE,
  run_repeated         = FALSE,
  run_missing_analysis = TRUE,

  # --- THRESHOLDS ---
  alpha             = 0.05,
  normality_alpha   = 0.05,
  p_adjust_method   = "BH",           # "bonferroni" | "holm" | "BH" | "fdr" | "none"
  outlier_threshold = 3,

  # --- OUTPUT ---
  output_dir    = "pipeline_output",
  export_excel  = TRUE,
  export_csv    = TRUE,
  export_plots  = TRUE,
  export_report = TRUE,               # Word .docx summary
  verbose       = TRUE
)
```

---

## 🧠 Automatic Test Selection

The pipeline reads your data and picks the right test — no guessing required.

```
Is n per group ≥ min_n_per_group?
  └─ YES →  Run normality tests (Shapiro-Wilk primary decision)
              └─ Normal?
                    ├─ YES → Check variance homogeneity (Levene's)
                    │           ├─ Equal var   → Student t / One-Way ANOVA   → Tukey HSD
                    │           └─ Unequal var → Welch t  / Welch ANOVA      → Games-Howell
                    └─ NO  → Mann-Whitney U / Kruskal-Wallis                 → Dunn test
```

Both parametric and non-parametric results are **always stored** — only the *primary* recommendation changes. You can inspect any test result from `RESULTS`.

---

## 📁 Output Structure

```
pipeline_output/
├── RESULTS.rds                     # Complete R object — reload anytime
├── statistical_report.docx         # Auto-generated Word summary report
│
├── plots/
│   ├── dist_{variable}.png         # Violin + box + beeswarm per outcome
│   ├── qq_{variable}.png           # Q-Q normality plots
│   ├── posthoc_{variable}.png      # Post-hoc significance brackets
│   ├── ridge_distributions.png     # All outcomes, scaled, ridge plot
│   ├── mean_se_barplot.png         # Group means ± SE
│   ├── pvalue_heatmap.png          # Raw vs adjusted p-values heatmap
│   ├── correlation_pearson.png     # Pearson corrplot
│   ├── correlation_spearman.png    # Spearman corrplot
│   ├── scatter_matrix.png          # Pairwise scatter with group coloring
│   ├── regression_diagnostics_{variable}.png
│   ├── missing_values.png
│   └── survival_km.png             # (if survival enabled)
│
└── tables/
    ├── full_results.xlsx           # All tables in one multi-sheet workbook
    ├── descriptives.csv
    ├── normality.csv
    ├── group_tests.csv
    ├── effect_sizes.csv
    ├── power.csv
    ├── outliers.csv
    └── missing.csv
```

---

## 📦 Accessing Results in R

```r
# Reload at any time
RESULTS <- readRDS("pipeline_output/RESULTS.rds")

# Key fields
RESULTS$descriptives    # Per-group descriptive stats
RESULTS$normality       # All four normality tests
RESULTS$group_tests     # Test statistics, raw + adjusted p-values
RESULTS$posthoc         # Post-hoc pairwise comparisons (list by variable)
RESULTS$effect_sizes    # All effect size estimates
RESULTS$power           # Observed power + required n
RESULTS$correlation     # Pearson and Spearman rcorr objects
RESULTS$regression      # OLS models, tidy coefficients, glance fit stats
RESULTS$outliers        # Flagged observations
RESULTS$missing         # Missing value summary
```

---

## 📋 Requirements

- **R ≥ 4.2**
- Internet connection on first run (packages auto-install from CRAN)
- No other setup needed

**Packages installed automatically:**

`tidyverse` · `data.table` · `rstatix` · `coin` · `DescTools` · `PMCMRplus` · `effectsize` · `pwr` · `emmeans` · `lme4` · `lmerTest` · `nortest` · `moments` · `corrplot` · `Hmisc` · `psych` · `car` · `MASS` · `broom` · `MuMIn` · `survival` · `survminer` · `ggplot2` · `ggpubr` · `ggbeeswarm` · `ggridges` · `patchwork` · `viridis` · `RColorBrewer` · `openxlsx` · `flextable` · `officer` · `janitor` · `skimr` · `naniar` · `mice`

---

## 🔬 Effect Size Benchmarks

| Measure | Small | Medium | Large |
|---|---|---|---|
| Cohen's d / Hedges' g | 0.2 | 0.5 | 0.8 |
| Eta-squared η² | 0.01 | 0.06 | 0.14 |
| Omega-squared ω² | 0.01 | 0.06 | 0.14 |
| Cohen's f (ANOVA) | 0.10 | 0.25 | 0.40 |
| Rank-biserial r | 0.10 | 0.30 | 0.50 |

---

## 🩺 Design Principles

- **Zero-config default** — sources and runs on a 90-row, 3-group synthetic example immediately, no data needed
- **Resilient** — every section is wrapped in `tryCatch`; one failure never stops the pipeline
- **Transparent** — every test selection decision is logged and stored in `RESULTS`
- **Reproducible** — `set.seed` is used throughout; full state saved as `.rds`
- **Modular** — toggle any analysis on or off with a single `TRUE`/`FALSE` in `CONFIG`
- **Non-destructive** — never overwrites your raw data; only writes to `output_dir`

---

## 🔧 Extending the Pipeline

New analysis sections slot in cleanly following this pattern:

```r
# Add to CONFIG
CONFIG$run_my_analysis <- FALSE

# Add a new section in the script
if (CONFIG$run_my_analysis) {
  log_msg("=== MY ANALYSIS ===")
  RESULTS$my_analysis <- safe_run("my_analysis", {
    # ... your code here ...
  })
}
```

Common extensions: logistic regression (`glm(family=binomial)`), Bayesian tests (`BayesFactor`), machine learning (`tidymodels`), Poisson/negative binomial for count data.

---

## 📝 Methods Paragraph (for papers)

> Data were analysed using R (v4.x) with an automated statistical pipeline (v1.0.0). Normality was assessed using Shapiro-Wilk tests. Group differences were evaluated using [test from `RESULTS$group_tests$test_used`] with p-values adjusted using the Benjamini-Hochberg procedure. Effect sizes were estimated as [Cohen's d / η²] and interpreted per Cohen (1988). Statistical power was computed using the `pwr` package. Statistical significance was defined as adjusted p < 0.05.

---

## 📄 Files

| File | Description |
|---|---|
| `statistical_pipeline.R` | The pipeline — source this file |
| `pipeline_documentation.docx` | Full documentation: all CONFIG options, decision logic, output reference, interpretation guide, troubleshooting |

---

## 📜 License

MIT — free to use, modify, and distribute. Attribution appreciated.