################################################################################
##                                                                            ##
##          AUTOMATED STATISTICAL ANALYSIS PIPELINE FOR R                    ##
##          Version 1.0.0 | Multi-Group Comprehensive Analysis               ##
##                                                                            ##
################################################################################
## DESCRIPTION:
##   Full-spectrum statistical pipeline. Drop in any data frame and get
##   descriptive stats, normality tests, group comparisons (parametric &
##   non-parametric), post-hoc, effect sizes, correlations, regression,
##   power analysis, and publication-ready exports.
##
## USAGE:
##   1. Set CONFIG block below.
##   2. source("statistical_pipeline.R")
##   3. All outputs land in OUTPUT_DIR.
################################################################################

# ============================================================
# 0. DEPENDENCIES
# ============================================================

required_packages <- c(
  # Core
  "tidyverse", "data.table",
  # Statistics
  "rstatix", "coin", "DescTools", "PMCMRplus",
  "effectsize", "pwr", "emmeans", "lme4", "lmerTest",
  # Normality / diagnostics
  "nortest", "moments",
  # Correlation
  "corrplot", "Hmisc", "psych",
  # Regression / modelling
  "car", "MASS", "broom", "MuMIn",
  # Survival (optional)
  "survival", "survminer",
  # Visualization
  "ggplot2", "ggpubr", "ggbeeswarm", "ggridges",
  "patchwork", "viridis", "RColorBrewer",
  # Export
  "openxlsx", "flextable", "officer", "knitr", "kableExtra",
  # Utilities
  "janitor", "skimr", "naniar", "mice"
)

install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) {
    message("Installing missing packages: ", paste(new, collapse = ", "))
    install.packages(new, repos = "https://cran.r-project.org", quiet = TRUE)
  }
}
install_if_missing(required_packages)

suppressPackageStartupMessages({
  lapply(required_packages, function(p) {
    tryCatch(library(p, character.only = TRUE),
             error = function(e) message("Could not load: ", p))
  })
})

# ============================================================
# 1. CONFIGURATION  ← Edit here
# ============================================================

CONFIG <- list(

  ## --- DATA SOURCE ---
  # Options: "builtin_example" | "csv" | "excel" | "rds" | "custom"
  data_source   = "builtin_example",
  data_path     = "",                  # path when data_source != "builtin_example"

  ## --- VARIABLE ROLES ---
  group_var     = "group",             # categorical grouping column
  outcome_vars  = NULL,                # NULL = auto-detect all numeric columns
  covariate_vars = NULL,               # columns to include as ANCOVA covariates
  id_var        = NULL,                # subject ID for repeated-measures
  time_var      = NULL,                # time column for longitudinal / survival
  event_var     = NULL,                # event indicator for survival (0/1)

  ## --- ANALYSIS TOGGLES ---
  run_descriptives    = TRUE,
  run_normality       = TRUE,
  run_group_tests     = TRUE,
  run_posthoc         = TRUE,
  run_effect_sizes    = TRUE,
  run_power           = TRUE,
  run_correlation     = TRUE,
  run_regression      = TRUE,
  run_ancova          = FALSE,         # set TRUE if covariate_vars supplied
  run_nonparametric   = TRUE,
  run_survival        = FALSE,         # requires time_var + event_var
  run_repeated        = FALSE,         # requires id_var + time_var
  run_missing_analysis = TRUE,

  ## --- THRESHOLDS ---
  alpha              = 0.05,
  normality_alpha    = 0.05,
  p_adjust_method    = "BH",           # "bonferroni" | "holm" | "BH" | "fdr" | "none"
  min_n_per_group    = 3,
  outlier_threshold  = 3,              # SD-based outlier flag (|z| > threshold)

  ## --- VISUALIZATION ---
  color_palette      = "Set2",         # RColorBrewer palette
  plot_width         = 10,
  plot_height        = 7,
  plot_dpi           = 300,
  theme_choice       = "bw",           # "bw" | "minimal" | "classic" | "dark"

  ## --- OUTPUT ---
  output_dir         = "pipeline_output",
  export_excel       = TRUE,
  export_csv         = TRUE,
  export_plots       = TRUE,
  export_report      = TRUE,           # Word .docx summary
  verbose            = TRUE
)

# ============================================================
# 2. SETUP HELPERS
# ============================================================

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(CONFIG$output_dir, "tables"), showWarnings = FALSE)

log_msg <- function(...) {
  if (CONFIG$verbose) message(format(Sys.time(), "[%H:%M:%S]"), " ", ...)
}

safe_run <- function(label, expr) {
  log_msg("Running: ", label)
  tryCatch(expr, error = function(e) {
    warning("FAILED [", label, "]: ", conditionMessage(e))
    NULL
  })
}

save_plot <- function(p, name, w = CONFIG$plot_width, h = CONFIG$plot_height) {
  if (!CONFIG$export_plots || is.null(p)) return(invisible(NULL))
  path <- file.path(CONFIG$output_dir, "plots", paste0(name, ".png"))
  ggsave(path, plot = p, width = w, height = h, dpi = CONFIG$plot_dpi,
         bg = "white")
  log_msg("  Saved plot: ", basename(path))
}

pval_fmt <- function(p) {
  ifelse(p < 0.001, "< 0.001", formatC(p, digits = 3, format = "f"))
}

gg_theme <- switch(CONFIG$theme_choice,
  bw      = theme_bw,
  minimal = theme_minimal,
  classic = theme_classic,
  dark    = theme_dark,
  theme_bw
)

# Results collector
RESULTS <- list()

# ============================================================
# 3. LOAD DATA
# ============================================================

log_msg("=== LOADING DATA ===")

load_data <- function(cfg) {
  switch(cfg$data_source,
    builtin_example = {
      set.seed(42)
      n <- 90
      data.frame(
        id     = 1:n,
        group  = rep(c("Control", "Treatment_A", "Treatment_B"), each = 30),
        age    = c(rnorm(30, 45, 8), rnorm(30, 47, 9), rnorm(30, 44, 7)),
        score  = c(rnorm(30, 50, 10), rnorm(30, 62, 11), rnorm(30, 58, 9)),
        weight = c(rnorm(30, 70, 12), rnorm(30, 68, 11), rnorm(30, 72, 13)),
        bp_sys = c(rnorm(30, 120, 15), rnorm(30, 115, 14), rnorm(30, 118, 16)),
        bp_dia = c(rnorm(30, 80, 10), rnorm(30, 76, 9), rnorm(30, 78, 11)),
        time   = rexp(n, 0.05),
        event  = rbinom(n, 1, 0.6),
        sex    = sample(c("M", "F"), n, replace = TRUE),
        stringsAsFactors = FALSE
      )
    },
    csv   = read.csv(cfg$data_path, stringsAsFactors = FALSE),
    excel = openxlsx::read.xlsx(cfg$data_path),
    rds   = readRDS(cfg$data_path),
    stop("Unknown data_source: ", cfg$data_source)
  )
}

raw_data <- load_data(CONFIG)
raw_data <- janitor::clean_names(raw_data)

# Normalise column names in config
CONFIG$group_var <- janitor::make_clean_names(CONFIG$group_var)

# Auto-detect outcome variables
if (is.null(CONFIG$outcome_vars)) {
  exclude <- c(CONFIG$group_var, CONFIG$id_var, CONFIG$time_var,
               CONFIG$event_var, CONFIG$covariate_vars)
  CONFIG$outcome_vars <- names(raw_data)[
    sapply(raw_data, is.numeric) & !names(raw_data) %in% exclude
  ]
}
CONFIG$outcome_vars <- janitor::make_clean_names(CONFIG$outcome_vars)

log_msg("Rows: ", nrow(raw_data), " | Cols: ", ncol(raw_data))
log_msg("Outcome variables: ", paste(CONFIG$outcome_vars, collapse = ", "))
log_msg("Groups: ", paste(unique(raw_data[[CONFIG$group_var]]), collapse = ", "))

df <- raw_data
df[[CONFIG$group_var]] <- as.factor(df[[CONFIG$group_var]])
RESULTS$data <- df

# ============================================================
# 4. MISSING DATA ANALYSIS
# ============================================================

if (CONFIG$run_missing_analysis) {
  log_msg("=== MISSING DATA ===")

  missing_summary <- naniar::miss_var_summary(df)
  RESULTS$missing <- missing_summary

  p_miss <- safe_run("missing plot", {
    naniar::gg_miss_var(df) + gg_theme() +
      labs(title = "Missing Values by Variable")
  })
  save_plot(p_miss, "missing_values")

  # Simple imputation report (do NOT auto-impute; just flag)
  if (any(missing_summary$n_miss > 0)) {
    log_msg("  Missing values detected. Review RESULTS$missing.")
    log_msg("  For imputation, set: df <- complete(mice(df, m=5, method='pmm'))")
  }
}

# ============================================================
# 5. DESCRIPTIVE STATISTICS
# ============================================================

if (CONFIG$run_descriptives) {
  log_msg("=== DESCRIPTIVE STATISTICS ===")

  ## 5a. Overall skim
  RESULTS$skim <- skimr::skim(df)

  ## 5b. Per-group descriptives for each outcome
  desc_all <- lapply(CONFIG$outcome_vars, function(v) {
    df %>%
      group_by(across(all_of(CONFIG$group_var))) %>%
      summarise(
        n         = n(),
        n_missing  = sum(is.na(.data[[v]])),
        mean      = mean(.data[[v]], na.rm = TRUE),
        sd        = sd(.data[[v]], na.rm = TRUE),
        se        = sd / sqrt(n),
        median    = median(.data[[v]], na.rm = TRUE),
        iqr       = IQR(.data[[v]], na.rm = TRUE),
        q25       = quantile(.data[[v]], 0.25, na.rm = TRUE),
        q75       = quantile(.data[[v]], 0.75, na.rm = TRUE),
        min       = min(.data[[v]], na.rm = TRUE),
        max       = max(.data[[v]], na.rm = TRUE),
        range     = max - min,
        cv_pct    = (sd / mean) * 100,
        skewness  = moments::skewness(.data[[v]], na.rm = TRUE),
        kurtosis  = moments::kurtosis(.data[[v]], na.rm = TRUE),
        ci_lower  = mean - qt(0.975, df = n - 1) * se,
        ci_upper  = mean + qt(0.975, df = n - 1) * se,
        .groups   = "drop"
      ) %>%
      mutate(variable = v, .before = 1)
  }) %>% bind_rows()

  RESULTS$descriptives <- desc_all

  ## 5c. Outlier detection
  outliers_all <- lapply(CONFIG$outcome_vars, function(v) {
    df %>%
      mutate(
        z_score  = as.numeric(scale(.data[[v]])),
        is_outlier = abs(z_score) > CONFIG$outlier_threshold
      ) %>%
      filter(is_outlier) %>%
      mutate(variable = v, value = .data[[v]]) %>%
      select(any_of(c(CONFIG$id_var, CONFIG$group_var, "variable", "value", "z_score")))
  }) %>% bind_rows()

  RESULTS$outliers <- outliers_all
  if (nrow(outliers_all) > 0)
    log_msg("  Outliers detected (|z| > ", CONFIG$outlier_threshold, "): ",
            nrow(outliers_all), " observations")

  ## 5d. Visualization: summary table plot
  p_desc <- ggpubr::ggsummary_table <- function() NULL  # placeholder

  # Box + violin + bee-swarm for each outcome
  for (v in CONFIG$outcome_vars) {
    p <- df %>%
      ggplot(aes(x = .data[[CONFIG$group_var]], y = .data[[v]],
                 fill = .data[[CONFIG$group_var]], color = .data[[CONFIG$group_var]])) +
      geom_violin(alpha = 0.3, trim = FALSE) +
      geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA,
                   color = "black", fill = "white") +
      geom_beeswarm(size = 1.5, alpha = 0.6, cex = 2) +
      scale_fill_brewer(palette = CONFIG$color_palette) +
      scale_color_brewer(palette = CONFIG$color_palette) +
      gg_theme() +
      theme(legend.position = "none") +
      labs(title = paste("Distribution of", v, "by Group"),
           x = CONFIG$group_var, y = v)
    save_plot(p, paste0("dist_", v))
  }

  # Ridge plot (all outcomes combined, scaled)
  ridge_data <- df %>%
    select(all_of(c(CONFIG$group_var, CONFIG$outcome_vars))) %>%
    pivot_longer(-all_of(CONFIG$group_var), names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    mutate(value_scaled = as.numeric(scale(value))) %>%
    ungroup()

  p_ridge <- ridge_data %>%
    ggplot(aes(x = value_scaled, y = .data[[CONFIG$group_var]],
               fill = .data[[CONFIG$group_var]])) +
    geom_density_ridges(alpha = 0.6, scale = 0.9) +
    facet_wrap(~variable, scales = "free_x") +
    scale_fill_brewer(palette = CONFIG$color_palette) +
    gg_theme() + theme(legend.position = "none") +
    labs(title = "Scaled Distributions by Group (Ridge Plot)",
         x = "Scaled Value", y = CONFIG$group_var)
  save_plot(p_ridge, "ridge_distributions", w = 14, h = 8)
}

# ============================================================
# 6. NORMALITY TESTS
# ============================================================

if (CONFIG$run_normality) {
  log_msg("=== NORMALITY TESTS ===")

  normality_results <- lapply(CONFIG$outcome_vars, function(v) {
    df %>%
      group_by(across(all_of(CONFIG$group_var))) %>%
      group_modify(~{
        x <- .x[[v]][!is.na(.x[[v]])]
        n <- length(x)
        tibble(
          variable      = v,
          n             = n,
          # Shapiro-Wilk (n < 5000)
          SW_stat   = if (n >= 3 && n <= 5000) shapiro.test(x)$statistic else NA_real_,
          SW_p      = if (n >= 3 && n <= 5000) shapiro.test(x)$p.value   else NA_real_,
          # Kolmogorov-Smirnov (Lilliefors)
          KS_stat   = if (n >= 4) nortest::lillie.test(x)$statistic else NA_real_,
          KS_p      = if (n >= 4) nortest::lillie.test(x)$p.value   else NA_real_,
          # Anderson-Darling
          AD_stat   = if (n >= 7) nortest::ad.test(x)$statistic else NA_real_,
          AD_p      = if (n >= 7) nortest::ad.test(x)$p.value   else NA_real_,
          # Jarque-Bera
          JB_stat   = if (n >= 8) moments::jarque.test(x)$statistic else NA_real_,
          JB_p      = if (n >= 8) moments::jarque.test(x)$p.value   else NA_real_,
          skewness  = moments::skewness(x),
          kurtosis  = moments::kurtosis(x),
          is_normal = (if (n >= 3 && n <= 5000) shapiro.test(x)$p.value else 1) >
                      CONFIG$normality_alpha
        )
      }) %>% ungroup()
  }) %>% bind_rows()

  RESULTS$normality <- normality_results

  # QQ plots
  for (v in CONFIG$outcome_vars) {
    p_qq <- df %>%
      ggplot(aes(sample = .data[[v]])) +
      stat_qq() + stat_qq_line(color = "red") +
      facet_wrap(as.formula(paste("~", CONFIG$group_var))) +
      gg_theme() +
      labs(title = paste("Q-Q Plot:", v), x = "Theoretical", y = "Sample")
    save_plot(p_qq, paste0("qq_", v))
  }
}

# ============================================================
# 7. GROUP COMPARISON TESTS
# ============================================================

if (CONFIG$run_group_tests) {
  log_msg("=== GROUP COMPARISON TESTS ===")

  n_groups <- nlevels(df[[CONFIG$group_var]])

  group_tests <- lapply(CONFIG$outcome_vars, function(v) {
    x   <- df[[v]]
    grp <- df[[CONFIG$group_var]]
    ok  <- !is.na(x) & !is.na(grp)
    x   <- x[ok]; grp <- droplevels(grp[ok])

    # Check group sizes
    g_sizes <- table(grp)
    if (any(g_sizes < CONFIG$min_n_per_group)) {
      return(tibble(variable = v, note = "Insufficient n in at least one group"))
    }

    # Normality decision (majority rule)
    norm_ok <- if (!is.null(RESULTS$normality)) {
      sub <- RESULTS$normality %>% filter(variable == v)
      mean(sub$is_normal, na.rm = TRUE) >= 0.5
    } else TRUE

    # Variance homogeneity
    levene_p <- tryCatch({
      car::leveneTest(x ~ grp)$`Pr(>F)`[1]
    }, error = function(e) NA_real_)
    bartlett_p <- tryCatch(bartlett.test(x ~ grp)$p.value, error = function(e) NA_real_)
    var_equal <- !is.na(levene_p) && levene_p > CONFIG$alpha

    row <- tibble(
      variable          = v,
      n_total           = length(x),
      n_groups          = n_groups,
      normality_assumed = norm_ok,
      levene_p          = levene_p,
      bartlett_p        = bartlett_p,
      var_homogeneous   = var_equal
    )

    if (n_groups == 2) {
      # ---------- TWO-GROUP ----------
      g1 <- x[grp == levels(grp)[1]]
      g2 <- x[grp == levels(grp)[2]]

      # Parametric
      tt <- tryCatch(t.test(g1, g2, var.equal = var_equal), error = function(e) NULL)
      wt <- tryCatch(t.test(g1, g2, var.equal = FALSE), error = function(e) NULL)  # Welch

      # Non-parametric
      mw <- tryCatch(wilcox.test(g1, g2, exact = FALSE, correct = TRUE),
                     error = function(e) NULL)

      row <- row %>% mutate(
        test_used        = if (norm_ok) "Student/Welch t-test" else "Mann-Whitney U",
        t_stat           = if (!is.null(tt)) tt$statistic else NA_real_,
        df_t             = if (!is.null(tt)) tt$parameter else NA_real_,
        t_pval           = if (!is.null(tt)) tt$p.value else NA_real_,
        welch_t          = if (!is.null(wt)) wt$statistic else NA_real_,
        welch_pval       = if (!is.null(wt)) wt$p.value else NA_real_,
        mw_W             = if (!is.null(mw)) mw$statistic else NA_real_,
        mw_pval          = if (!is.null(mw)) mw$p.value else NA_real_,
        primary_pval     = if (norm_ok) wt$p.value else mw$p.value,
        significant      = !is.na(primary_pval) && primary_pval < CONFIG$alpha,
        mean_diff        = mean(g1, na.rm=TRUE) - mean(g2, na.rm=TRUE),
        ci_diff_lower    = if (!is.null(wt)) wt$conf.int[1] else NA_real_,
        ci_diff_upper    = if (!is.null(wt)) wt$conf.int[2] else NA_real_
      )

    } else {
      # ---------- MULTI-GROUP ----------
      fmla <- as.formula(paste(v, "~", CONFIG$group_var))

      # One-Way ANOVA
      aov_res  <- tryCatch(aov(fmla, data = df[ok,]), error = function(e) NULL)
      aov_sum  <- if (!is.null(aov_res)) summary(aov_res)[[1]] else NULL

      # Welch ANOVA (robust to heteroscedasticity)
      welch_res <- tryCatch(oneway.test(fmla, data = df[ok,], var.equal = FALSE),
                            error = function(e) NULL)

      # Brown-Forsythe
      bf_res <- tryCatch(
        PMCMRplus::welchManyOneTTest(fmla, data = df[ok,]),
        error = function(e) NULL
      )

      # Kruskal-Wallis
      kw_res <- tryCatch(kruskal.test(fmla, data = df[ok,]), error = function(e) NULL)

      # Eta-squared (ANOVA)
      eta2 <- if (!is.null(aov_res)) {
        tryCatch(effectsize::eta_squared(aov_res)$Eta2[1], error = function(e) NA_real_)
      } else NA_real_

      row <- row %>% mutate(
        test_used    = if (norm_ok) "One-Way ANOVA" else "Kruskal-Wallis",
        F_stat       = if (!is.null(aov_sum)) aov_sum$`F value`[1] else NA_real_,
        df_between   = if (!is.null(aov_sum)) aov_sum$Df[1] else NA_real_,
        df_within    = if (!is.null(aov_sum)) aov_sum$Df[2] else NA_real_,
        anova_pval   = if (!is.null(aov_sum)) aov_sum$`Pr(>F)`[1] else NA_real_,
        welch_F      = if (!is.null(welch_res)) welch_res$statistic else NA_real_,
        welch_pval   = if (!is.null(welch_res)) welch_res$p.value else NA_real_,
        kw_H         = if (!is.null(kw_res)) kw_res$statistic else NA_real_,
        kw_pval      = if (!is.null(kw_res)) kw_res$p.value else NA_real_,
        eta_squared  = eta2,
        primary_pval = if (norm_ok && var_equal) {
          if (!is.null(aov_sum)) aov_sum$`Pr(>F)`[1] else NA_real_
        } else if (norm_ok && !var_equal) {
          if (!is.null(welch_res)) welch_res$p.value else NA_real_
        } else {
          if (!is.null(kw_res)) kw_res$p.value else NA_real_
        },
        significant  = !is.na(primary_pval) && primary_pval < CONFIG$alpha
      )
    }
    row
  }) %>% bind_rows()

  # Adjust p-values across outcomes
  group_tests$primary_pval_adj <- p.adjust(group_tests$primary_pval,
                                            method = CONFIG$p_adjust_method)
  group_tests$significant_adj  <- group_tests$primary_pval_adj < CONFIG$alpha

  RESULTS$group_tests <- group_tests
}

# ============================================================
# 8. POST-HOC TESTS
# ============================================================

if (CONFIG$run_posthoc && !is.null(RESULTS$group_tests)) {
  log_msg("=== POST-HOC TESTS ===")

  n_groups <- nlevels(df[[CONFIG$group_var]])

  if (n_groups >= 3) {
    posthoc_all <- lapply(CONFIG$outcome_vars, function(v) {
      fmla <- as.formula(paste(v, "~", CONFIG$group_var))

      row_info <- RESULTS$group_tests %>% filter(variable == v)
      norm_ok   <- if (nrow(row_info) > 0) row_info$normality_assumed[1] else TRUE
      var_eq    <- if (nrow(row_info) > 0) row_info$var_homogeneous[1] else TRUE

      out <- list()

      # Tukey HSD (parametric, equal variance)
      out$tukey <- tryCatch({
        aov_fit <- aov(fmla, data = df)
        TukeyHSD(aov_fit)[[CONFIG$group_var]] %>%
          as.data.frame() %>%
          rownames_to_column("comparison") %>%
          mutate(variable = v, method = "Tukey HSD")
      }, error = function(e) NULL)

      # Games-Howell (parametric, unequal variance)
      out$games_howell <- tryCatch({
        rstatix::games_howell_test(df, as.formula(paste(v, "~", CONFIG$group_var))) %>%
          mutate(variable = v, method = "Games-Howell")
      }, error = function(e) NULL)

      # Dunn test (non-parametric)
      out$dunn <- tryCatch({
        rstatix::dunn_test(df, as.formula(paste(v, "~", CONFIG$group_var)),
                           p.adjust.method = CONFIG$p_adjust_method) %>%
          mutate(variable = v, method = "Dunn (KW post-hoc)")
      }, error = function(e) NULL)

      # Wilcoxon pairwise (non-parametric with FDR)
      out$pairwise_wilcox <- tryCatch({
        rstatix::pairwise_wilcox_test(df,
          as.formula(paste(v, "~", CONFIG$group_var)),
          p.adjust.method = CONFIG$p_adjust_method) %>%
          mutate(variable = v, method = "Pairwise Wilcoxon")
      }, error = function(e) NULL)

      # Primary recommended test
      out$recommended <- if (norm_ok && var_eq) out$tukey else
                         if (norm_ok && !var_eq) out$games_howell else
                         out$dunn

      out
    })
    names(posthoc_all) <- CONFIG$outcome_vars
    RESULTS$posthoc <- posthoc_all

    # Significance bracket plots
    for (v in CONFIG$outcome_vars) {
      stat_test <- tryCatch({
        rstatix::tukey_hsd(df, as.formula(paste(v, "~", CONFIG$group_var)))
      }, error = function(e) NULL)

      if (!is.null(stat_test)) {
        p_ph <- ggpubr::ggboxplot(df, x = CONFIG$group_var, y = v,
                                   fill = CONFIG$group_var,
                                   palette = CONFIG$color_palette) +
          ggpubr::stat_pvalue_manual(stat_test, label = "p.adj.signif",
                                     step.increase = 0.1) +
          gg_theme() + theme(legend.position = "none") +
          labs(title = paste("Post-hoc comparisons:", v))
        save_plot(p_ph, paste0("posthoc_", v))
      }
    }
  }
}

# ============================================================
# 9. EFFECT SIZES
# ============================================================

if (CONFIG$run_effect_sizes) {
  log_msg("=== EFFECT SIZES ===")

  n_groups <- nlevels(df[[CONFIG$group_var]])

  effect_sizes <- lapply(CONFIG$outcome_vars, function(v) {
    row <- tibble(variable = v)

    if (n_groups == 2) {
      g1 <- df[[v]][df[[CONFIG$group_var]] == levels(df[[CONFIG$group_var]])[1]]
      g2 <- df[[v]][df[[CONFIG$group_var]] == levels(df[[CONFIG$group_var]])[2]]
      g1 <- g1[!is.na(g1)]; g2 <- g2[!is.na(g2)]

      # Cohen's d
      cd <- tryCatch(effectsize::cohens_d(g1, g2), error = function(e) NULL)
      # Hedges' g
      hg <- tryCatch(effectsize::hedges_g(g1, g2), error = function(e) NULL)
      # Glass's delta
      gd <- tryCatch(effectsize::glass_delta(g1, g2), error = function(e) NULL)
      # Rank-biserial r (for Mann-Whitney)
      rb <- tryCatch(effectsize::rank_biserial(g1, g2), error = function(e) NULL)

      row <- row %>% mutate(
        cohens_d       = if (!is.null(cd)) cd$Cohens_d else NA_real_,
        hedges_g       = if (!is.null(hg)) hg$Hedges_g else NA_real_,
        glass_delta    = if (!is.null(gd)) gd$Glass_delta else NA_real_,
        rank_biserial  = if (!is.null(rb)) rb$r_rank_biserial else NA_real_,
        interpretation = case_when(
          !is.null(cd) && abs(cd$Cohens_d) < 0.2 ~ "negligible",
          !is.null(cd) && abs(cd$Cohens_d) < 0.5 ~ "small",
          !is.null(cd) && abs(cd$Cohens_d) < 0.8 ~ "medium",
          !is.null(cd)                            ~ "large",
          TRUE                                    ~ NA_character_
        )
      )
    } else {
      fmla    <- as.formula(paste(v, "~", CONFIG$group_var))
      aov_fit <- tryCatch(aov(fmla, data = df), error = function(e) NULL)

      eta2   <- tryCatch(effectsize::eta_squared(aov_fit)$Eta2[1], error = function(e) NA_real_)
      omega2 <- tryCatch(effectsize::omega_squared(aov_fit)$Omega2[1], error = function(e) NA_real_)
      eps2   <- tryCatch(effectsize::epsilon_squared(aov_fit)$Epsilon2[1], error = function(e) NA_real_)

      # Epsilon-squared for Kruskal-Wallis
      kw_eps <- tryCatch({
        kw <- kruskal.test(fmla, data = df)
        n  <- sum(!is.na(df[[v]]))
        k  <- nlevels(df[[CONFIG$group_var]])
        (kw$statistic - k + 1) / (n - k)
      }, error = function(e) NA_real_)

      row <- row %>% mutate(
        eta_squared    = eta2,
        omega_squared  = omega2,
        epsilon_squared = eps2,
        kw_epsilon_sq  = kw_eps,
        interpretation = case_when(
          !is.na(eta2) & eta2 < 0.01 ~ "negligible",
          !is.na(eta2) & eta2 < 0.06 ~ "small",
          !is.na(eta2) & eta2 < 0.14 ~ "medium",
          !is.na(eta2)               ~ "large",
          TRUE                       ~ NA_character_
        )
      )
    }
    row
  }) %>% bind_rows()

  RESULTS$effect_sizes <- effect_sizes
}

# ============================================================
# 10. POWER ANALYSIS
# ============================================================

if (CONFIG$run_power) {
  log_msg("=== POWER ANALYSIS ===")

  n_groups <- nlevels(df[[CONFIG$group_var]])

  power_results <- lapply(CONFIG$outcome_vars, function(v) {
    n_per_group <- min(table(df[[CONFIG$group_var]]))

    es_row <- RESULTS$effect_sizes %>% filter(variable == v)

    if (n_groups == 2) {
      d <- if (!is.null(es_row) && nrow(es_row) > 0) es_row$cohens_d[1] else 0.5
      d <- if (is.na(d)) 0.5 else abs(d)

      pw <- tryCatch(pwr::pwr.t.test(n = n_per_group, d = d,
                                      sig.level = CONFIG$alpha, type = "two.sample"),
                     error = function(e) NULL)
      n_needed <- tryCatch(pwr::pwr.t.test(d = d, sig.level = CONFIG$alpha,
                                            power = 0.80, type = "two.sample")$n,
                           error = function(e) NA_real_)

      tibble(variable = v, test = "t-test", effect_size = d,
             observed_n_per_group = n_per_group,
             observed_power       = if (!is.null(pw)) pw$power else NA_real_,
             n_for_80pct_power    = ceiling(n_needed),
             n_for_90pct_power    = ceiling(tryCatch(
               pwr::pwr.t.test(d=d, sig.level=CONFIG$alpha, power=0.90,
                               type="two.sample")$n, error=function(e) NA_real_)))

    } else {
      f_val <- if (!is.null(es_row) && nrow(es_row) > 0 && !is.na(es_row$eta_squared[1])) {
        effectsize::eta2_to_f(es_row$eta_squared[1])
      } else 0.25

      pw <- tryCatch(pwr::pwr.anova.test(k = n_groups, n = n_per_group,
                                          f = f_val, sig.level = CONFIG$alpha),
                     error = function(e) NULL)
      n_needed <- tryCatch(pwr::pwr.anova.test(k = n_groups, f = f_val,
                                                sig.level = CONFIG$alpha, power = 0.80)$n,
                           error = function(e) NA_real_)

      tibble(variable = v, test = "ANOVA", effect_size = f_val,
             observed_n_per_group = n_per_group,
             observed_power       = if (!is.null(pw)) pw$power else NA_real_,
             n_for_80pct_power    = ceiling(n_needed),
             n_for_90pct_power    = ceiling(tryCatch(
               pwr::pwr.anova.test(k=n_groups, f=f_val, sig.level=CONFIG$alpha,
                                   power=0.90)$n, error=function(e) NA_real_)))
    }
  }) %>% bind_rows()

  RESULTS$power <- power_results
}

# ============================================================
# 11. CORRELATION ANALYSIS
# ============================================================

if (CONFIG$run_correlation) {
  log_msg("=== CORRELATION ANALYSIS ===")

  num_df <- df %>% select(all_of(CONFIG$outcome_vars)) %>% drop_na()

  # Pearson
  cor_pearson <- tryCatch(Hmisc::rcorr(as.matrix(num_df), type = "pearson"),
                          error = function(e) NULL)
  # Spearman
  cor_spearman <- tryCatch(Hmisc::rcorr(as.matrix(num_df), type = "spearman"),
                           error = function(e) NULL)

  RESULTS$correlation <- list(pearson = cor_pearson, spearman = cor_spearman)

  # Correlation plots
  if (!is.null(cor_pearson)) {
    png(file.path(CONFIG$output_dir, "plots", "correlation_pearson.png"),
        width = CONFIG$plot_width, height = CONFIG$plot_height,
        units = "in", res = CONFIG$plot_dpi, bg = "white")
    corrplot::corrplot.mixed(cor_pearson$r,
      upper      = "ellipse",
      lower      = "number",
      tl.col     = "black",
      tl.cex     = 0.8,
      p.mat      = cor_pearson$P,
      sig.level  = CONFIG$alpha,
      insig      = "label_sig",
      pch.cex    = 1.2,
      title      = "Pearson Correlation Matrix",
      mar        = c(0, 0, 2, 0)
    )
    dev.off()
  }

  if (!is.null(cor_spearman)) {
    png(file.path(CONFIG$output_dir, "plots", "correlation_spearman.png"),
        width = CONFIG$plot_width, height = CONFIG$plot_height,
        units = "in", res = CONFIG$plot_dpi, bg = "white")
    corrplot::corrplot.mixed(cor_spearman$r,
      upper      = "ellipse",
      lower      = "number",
      tl.col     = "black",
      tl.cex     = 0.8,
      p.mat      = cor_spearman$P,
      sig.level  = CONFIG$alpha,
      insig      = "label_sig",
      pch.cex    = 1.2,
      title      = "Spearman Correlation Matrix",
      mar        = c(0, 0, 2, 0)
    )
    dev.off()
  }

  # Scatter matrix
  if (length(CONFIG$outcome_vars) >= 2 && length(CONFIG$outcome_vars) <= 8) {
    p_pairs <- GGally::ggpairs(
      df, columns = CONFIG$outcome_vars,
      mapping = aes(color = .data[[CONFIG$group_var]], alpha = 0.6),
      upper = list(continuous = "cor"),
      lower = list(continuous = "smooth"),
      diag  = list(continuous = "densityDiag")
    ) + scale_color_brewer(palette = CONFIG$color_palette) +
      gg_theme() +
      labs(title = "Scatter Matrix with Group Coloring")
    save_plot(p_pairs, "scatter_matrix", w = 14, h = 12)
  }
}

# ============================================================
# 12. REGRESSION ANALYSIS
# ============================================================

if (CONFIG$run_regression) {
  log_msg("=== REGRESSION ANALYSIS ===")

  regression_results <- list()

  for (v in CONFIG$outcome_vars) {
    # Simple: outcome ~ group
    fmla_simple <- as.formula(paste(v, "~", CONFIG$group_var))

    # Full: add all other outcomes as predictors
    other_vars <- setdiff(CONFIG$outcome_vars, v)
    fmla_full  <- if (length(other_vars) > 0) {
      as.formula(paste(v, "~", CONFIG$group_var, "+",
                       paste(other_vars, collapse = " + ")))
    } else fmla_simple

    lm_simple <- tryCatch(lm(fmla_simple, data = df), error = function(e) NULL)
    lm_full   <- tryCatch(lm(fmla_full,   data = df), error = function(e) NULL)

    regression_results[[v]] <- list(
      simple_model   = lm_simple,
      simple_tidy    = if (!is.null(lm_simple)) broom::tidy(lm_simple, conf.int = TRUE) else NULL,
      simple_glance  = if (!is.null(lm_simple)) broom::glance(lm_simple) else NULL,
      full_model     = lm_full,
      full_tidy      = if (!is.null(lm_full)) broom::tidy(lm_full, conf.int = TRUE) else NULL,
      full_glance    = if (!is.null(lm_full)) broom::glance(lm_full) else NULL
    )

    # Diagnostic plots
    if (!is.null(lm_full)) {
      png(file.path(CONFIG$output_dir, "plots", paste0("regression_diagnostics_", v, ".png")),
          width = 12, height = 10, units = "in", res = CONFIG$plot_dpi, bg = "white")
      par(mfrow = c(2, 2))
      plot(lm_full, main = paste("Regression Diagnostics:", v))
      dev.off()
    }
  }

  RESULTS$regression <- regression_results
}

# ============================================================
# 13. ANCOVA
# ============================================================

if (CONFIG$run_ancova && !is.null(CONFIG$covariate_vars)) {
  log_msg("=== ANCOVA ===")

  ancova_results <- lapply(CONFIG$outcome_vars, function(v) {
    fmla <- as.formula(paste(v, "~", CONFIG$group_var, "+",
                             paste(CONFIG$covariate_vars, collapse = " + ")))
    fit  <- tryCatch(aov(fmla, data = df), error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    list(
      model     = fit,
      summary   = summary(fit),
      tidy      = broom::tidy(fit),
      emmeans   = tryCatch(emmeans::emmeans(fit, CONFIG$group_var) %>%
                             pairs() %>% as.data.frame(), error = function(e) NULL),
      eta2      = tryCatch(effectsize::eta_squared(fit), error = function(e) NULL)
    )
  })
  names(ancova_results) <- CONFIG$outcome_vars
  RESULTS$ancova <- ancova_results
}

# ============================================================
# 14. REPEATED MEASURES
# ============================================================

if (CONFIG$run_repeated && !is.null(CONFIG$id_var) && !is.null(CONFIG$time_var)) {
  log_msg("=== REPEATED MEASURES / MIXED MODELS ===")

  rm_results <- lapply(CONFIG$outcome_vars, function(v) {
    fmla_lmer <- as.formula(
      paste(v, "~ (1|", CONFIG$id_var, ") +", CONFIG$time_var, "*", CONFIG$group_var)
    )
    fit <- tryCatch(lmerTest::lmer(fmla_lmer, data = df, REML = TRUE),
                   error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    list(
      model   = fit,
      anova   = tryCatch(anova(fit), error = function(e) NULL),
      tidy    = broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE),
      icc     = tryCatch(performance::icc(fit), error = function(e) NULL),
      r2      = tryCatch(MuMIn::r.squaredGLMM(fit), error = function(e) NULL)
    )
  })
  names(rm_results) <- CONFIG$outcome_vars
  RESULTS$repeated_measures <- rm_results
}

# ============================================================
# 15. SURVIVAL ANALYSIS
# ============================================================

if (CONFIG$run_survival && !is.null(CONFIG$time_var) && !is.null(CONFIG$event_var)) {
  log_msg("=== SURVIVAL ANALYSIS ===")

  surv_obj <- tryCatch(
    survival::Surv(df[[CONFIG$time_var]], df[[CONFIG$event_var]]),
    error = function(e) NULL
  )

  if (!is.null(surv_obj)) {
    fmla_km <- as.formula(paste("surv_obj ~", CONFIG$group_var))
    km_fit  <- tryCatch(survival::survfit(fmla_km, data = df), error = function(e) NULL)
    lr_test <- tryCatch(survival::survdiff(fmla_km, data = df), error = function(e) NULL)
    cox_fit <- tryCatch(survival::coxph(fmla_km, data = df), error = function(e) NULL)

    RESULTS$survival <- list(
      km_fit   = km_fit,
      logrank  = lr_test,
      cox_fit  = cox_fit,
      cox_tidy = if (!is.null(cox_fit)) broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) else NULL
    )

    if (!is.null(km_fit)) {
      p_km <- survminer::ggsurvplot(
        km_fit, data = df,
        pval         = TRUE,
        conf.int     = TRUE,
        risk.table   = TRUE,
        palette      = CONFIG$color_palette,
        ggtheme      = gg_theme(),
        title        = "Kaplan-Meier Survival Curves"
      )
      ggsave(file.path(CONFIG$output_dir, "plots", "survival_km.png"),
             plot = p_km$plot, width = CONFIG$plot_width, height = CONFIG$plot_height,
             dpi = CONFIG$plot_dpi, bg = "white")
    }
  }
}

# ============================================================
# 16. COMPREHENSIVE VISUALIZATION PANEL
# ============================================================

log_msg("=== GENERATING SUMMARY VISUALIZATIONS ===")

# Mean ± SE bar chart for all outcomes
if (!is.null(RESULTS$descriptives)) {
  p_bar <- RESULTS$descriptives %>%
    ggplot(aes(x = .data[[CONFIG$group_var]], y = mean,
               fill = .data[[CONFIG$group_var]], color = .data[[CONFIG$group_var]])) +
    geom_col(alpha = 0.7, position = position_dodge()) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.2, color = "black",
                  position = position_dodge(0.9)) +
    facet_wrap(~variable, scales = "free_y") +
    scale_fill_brewer(palette = CONFIG$color_palette) +
    scale_color_brewer(palette = CONFIG$color_palette) +
    gg_theme() + theme(legend.position = "none",
                       axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = "Group Means ± SE", x = NULL, y = "Mean ± SE")
  save_plot(p_bar, "mean_se_barplot", w = 14, h = 10)
}

# P-value heatmap
if (!is.null(RESULTS$group_tests) && nrow(RESULTS$group_tests) > 0) {
  pval_heat_data <- RESULTS$group_tests %>%
    select(variable, primary_pval, primary_pval_adj) %>%
    pivot_longer(-variable, names_to = "type", values_to = "pval") %>%
    mutate(neg_log10_p = -log10(pval),
           sig_label   = case_when(pval < 0.001 ~ "***",
                                   pval < 0.01  ~ "**",
                                   pval < 0.05  ~ "*",
                                   TRUE         ~ "ns"))

  p_pheat <- pval_heat_data %>%
    ggplot(aes(x = type, y = variable, fill = neg_log10_p)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sig_label), size = 5) +
    scale_fill_viridis_c(name = "-log10(p)", option = "magma", direction = -1) +
    gg_theme() +
    labs(title = "P-value Heatmap (raw vs adjusted)",
         x = "P-value Type", y = "Outcome Variable")
  save_plot(p_pheat, "pvalue_heatmap")
}

# ============================================================
# 17. EXPORT ALL RESULTS
# ============================================================

log_msg("=== EXPORTING RESULTS ===")

## 17a. Excel workbook
if (CONFIG$export_excel) {
  wb <- openxlsx::createWorkbook()

  add_sheet <- function(wb, sheet_name, data, ...) {
    if (!is.null(data) && (is.data.frame(data) || is.matrix(data))) {
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeDataTable(wb, sheet_name, as.data.frame(data), ...)
    }
  }

  add_sheet(wb, "Descriptives",   RESULTS$descriptives)
  add_sheet(wb, "Normality",      RESULTS$normality)
  add_sheet(wb, "GroupTests",     RESULTS$group_tests)
  add_sheet(wb, "EffectSizes",    RESULTS$effect_sizes)
  add_sheet(wb, "PowerAnalysis",  RESULTS$power)
  add_sheet(wb, "Missing",        RESULTS$missing)
  add_sheet(wb, "Outliers",       RESULTS$outliers)

  if (!is.null(RESULTS$posthoc)) {
    lapply(names(RESULTS$posthoc), function(v) {
      ph <- RESULTS$posthoc[[v]]$recommended
      if (!is.null(ph) && is.data.frame(ph)) {
        nm <- paste0("PostHoc_", substr(v, 1, 20))
        openxlsx::addWorksheet(wb, nm)
        openxlsx::writeDataTable(wb, nm, ph)
      }
    })
  }

  if (!is.null(RESULTS$regression)) {
    reg_tidy <- lapply(names(RESULTS$regression), function(v) {
      r <- RESULTS$regression[[v]]$full_tidy
      if (!is.null(r)) mutate(r, variable = v)
    }) %>% bind_rows()
    add_sheet(wb, "Regression", reg_tidy)

    reg_glance <- lapply(names(RESULTS$regression), function(v) {
      g <- RESULTS$regression[[v]]$full_glance
      if (!is.null(g)) mutate(g, variable = v)
    }) %>% bind_rows()
    add_sheet(wb, "RegressionFit", reg_glance)
  }

  if (!is.null(RESULTS$survival$cox_tidy)) {
    add_sheet(wb, "Survival_Cox", RESULTS$survival$cox_tidy)
  }

  xlsx_path <- file.path(CONFIG$output_dir, "tables", "full_results.xlsx")
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  log_msg("  Excel: ", xlsx_path)
}

## 17b. CSV exports
if (CONFIG$export_csv) {
  csv_map <- list(
    descriptives = RESULTS$descriptives,
    normality    = RESULTS$normality,
    group_tests  = RESULTS$group_tests,
    effect_sizes = RESULTS$effect_sizes,
    power        = RESULTS$power,
    outliers     = RESULTS$outliers,
    missing      = RESULTS$missing
  )
  lapply(names(csv_map), function(nm) {
    d <- csv_map[[nm]]
    if (!is.null(d) && is.data.frame(d)) {
      write.csv(d,
        file.path(CONFIG$output_dir, "tables", paste0(nm, ".csv")),
        row.names = FALSE
      )
    }
  })
  log_msg("  CSVs saved to: ", file.path(CONFIG$output_dir, "tables"))
}

## 17c. RDS (complete object)
saveRDS(RESULTS, file.path(CONFIG$output_dir, "RESULTS.rds"))
log_msg("  RDS: ", file.path(CONFIG$output_dir, "RESULTS.rds"))

# ============================================================
# 18. WORD REPORT
# ============================================================

if (CONFIG$export_report) {
  log_msg("=== GENERATING WORD REPORT ===")

  doc <- officer::read_docx()

  add_h1 <- function(doc, txt) officer::body_add_par(doc, txt, style = "heading 1")
  add_h2 <- function(doc, txt) officer::body_add_par(doc, txt, style = "heading 2")
  add_p  <- function(doc, txt) officer::body_add_par(doc, txt, style = "Normal")
  add_ft <- function(doc, df_in, caption = "") {
    if (is.null(df_in) || !is.data.frame(df_in) || nrow(df_in) == 0) return(doc)
    ft <- flextable::flextable(head(as.data.frame(df_in), 50)) %>%
      flextable::autofit() %>%
      flextable::theme_vanilla()
    if (nchar(caption) > 0)
      ft <- flextable::set_caption(ft, caption)
    flextable::body_add_flextable(doc, ft)
  }

  # Cover
  doc <- add_h1(doc, "Automated Statistical Analysis Report")
  doc <- add_p(doc, paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  doc <- add_p(doc, paste("Data rows:", nrow(df), "| Variables:", ncol(df)))
  doc <- add_p(doc, paste("Grouping variable:", CONFIG$group_var))
  doc <- add_p(doc, paste("Outcome variables:", paste(CONFIG$outcome_vars, collapse = ", ")))
  doc <- add_p(doc, paste("Alpha:", CONFIG$alpha, "| P-adjust:", CONFIG$p_adjust_method))
  doc <- officer::body_add_par(doc, "", style = "Normal")

  # Descriptives
  doc <- add_h1(doc, "1. Descriptive Statistics")
  doc <- add_ft(doc, RESULTS$descriptives, "Descriptive statistics by group")

  # Normality
  if (!is.null(RESULTS$normality)) {
    doc <- add_h1(doc, "2. Normality Tests")
    doc <- add_ft(doc, RESULTS$normality, "Normality test results (SW = Shapiro-Wilk, KS = Kolmogorov-Smirnov, AD = Anderson-Darling, JB = Jarque-Bera)")
  }

  # Group tests
  if (!is.null(RESULTS$group_tests)) {
    doc <- add_h1(doc, "3. Group Comparison Tests")
    doc <- add_ft(doc, RESULTS$group_tests %>%
                    select(variable, test_used, primary_pval, primary_pval_adj,
                           significant, significant_adj, any_of(c("F_stat", "t_stat",
                             "kw_H", "eta_squared", "mean_diff"))),
                  "Group test results")
  }

  # Effect sizes
  if (!is.null(RESULTS$effect_sizes)) {
    doc <- add_h1(doc, "4. Effect Sizes")
    doc <- add_ft(doc, RESULTS$effect_sizes, "Effect size estimates")
  }

  # Power
  if (!is.null(RESULTS$power)) {
    doc <- add_h1(doc, "5. Power Analysis")
    doc <- add_ft(doc, RESULTS$power, "Statistical power and required sample sizes")
  }

  # Correlation
  doc <- add_h1(doc, "6. Correlation Analysis")
  if (!is.null(RESULTS$correlation$pearson)) {
    doc <- add_p(doc, "Pearson correlation matrix (r values):")
    doc <- add_ft(doc, as.data.frame(round(RESULTS$correlation$pearson$r, 3)),
                  "Pearson r")
    doc <- add_p(doc, "Pearson p-values:")
    doc <- add_ft(doc, as.data.frame(round(RESULTS$correlation$pearson$P, 4)),
                  "Pearson p-values")
  }

  # Regression
  if (!is.null(RESULTS$regression)) {
    doc <- add_h1(doc, "7. Regression Analysis")
    for (v in names(RESULTS$regression)) {
      doc <- add_h2(doc, paste("Outcome:", v))
      doc <- add_ft(doc, RESULTS$regression[[v]]$full_tidy,
                    paste("Regression coefficients:", v))
      doc <- add_ft(doc, RESULTS$regression[[v]]$full_glance,
                    paste("Model fit:", v))
    }
  }

  # Missing
  if (!is.null(RESULTS$missing)) {
    doc <- add_h1(doc, "8. Missing Data")
    doc <- add_ft(doc, RESULTS$missing, "Missing value summary")
  }

  # Outliers
  if (!is.null(RESULTS$outliers) && nrow(RESULTS$outliers) > 0) {
    doc <- add_h1(doc, "9. Outliers")
    doc <- add_ft(doc, RESULTS$outliers, paste0("Outliers (|z| > ", CONFIG$outlier_threshold, ")"))
  }

  report_path <- file.path(CONFIG$output_dir, "statistical_report.docx")
  print(doc, target = report_path)
  log_msg("  Report: ", report_path)
}

# ============================================================
# 19. CONSOLE SUMMARY
# ============================================================

log_msg("")
log_msg("============================================================")
log_msg("  PIPELINE COMPLETE")
log_msg("============================================================")

if (!is.null(RESULTS$group_tests)) {
  sig <- RESULTS$group_tests %>% filter(significant_adj == TRUE)
  log_msg("  Significant outcomes (adjusted p < ", CONFIG$alpha, "): ",
          if (nrow(sig) > 0) paste(sig$variable, collapse = ", ") else "none")
}

log_msg("  Outputs in: ", normalizePath(CONFIG$output_dir))
log_msg("  Plots:      ", file.path(CONFIG$output_dir, "plots"))
log_msg("  Tables:     ", file.path(CONFIG$output_dir, "tables"))
log_msg("  Report:     ", file.path(CONFIG$output_dir, "statistical_report.docx"))
log_msg("  RDS cache:  ", file.path(CONFIG$output_dir, "RESULTS.rds"))
log_msg("")
log_msg("  Access results: RESULTS$descriptives, RESULTS$group_tests, etc.")
log_msg("============================================================")

# Return invisibly for interactive use
invisible(RESULTS)