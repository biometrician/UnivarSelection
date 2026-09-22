################################################################################
# UnivarSelection_V1.R
# Consequences of Applying Univariable Preselection Compared to
# Backward Elimination - Educational Shiny app
#
# September 2026
# Daniela Dunkler, Theresa Ullmann
#
# Note: deployment to shinyapps.io only works if the file is outside saved
#       the R project
################################################################################

library(shiny)
library(bslib)
library(MASS)
library(abe)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

################################################################################
# Helper functions
################################################################################

calc_sigma2 <- function(beta, rho, r2) {
  D_P <- mvrnorm(n = 100000, mu = c(0, 0),
                 Sigma = matrix(c(1, rho, rho, 1), ncol = 2))
  Xb  <- as.matrix(D_P) %*% beta
  var(Xb) * (1 - r2) / r2
}

# Full simulation: returns all estimates, SEs, p-values, selection flags, and BE estimates
run_simulation_full <- function(nsim, n, beta, rho, sigma2, alpha, seed,
                                progress_fn = NULL) {
  set.seed(seed)
  Sig <- matrix(c(1, rho, rho, 1), ncol = 2)

  b1_uv  <- b2_uv  <- b1_mv  <- b2_mv  <- numeric(nsim)
  se1_uv <- se2_uv <- se1_mv <- se2_mv <- numeric(nsim)
  p1_uv  <- p2_uv  <- p1_mv  <- p2_mv  <- numeric(nsim)
  # BE: 0 when variable dropped (unconditional mean = selection bias included)
  b1_be  <- b2_be  <- numeric(nsim)
  se1_be <- se2_be <- numeric(nsim)
  p1_be  <- p2_be  <- rep(NA_real_, nsim)   # NA when dropped
  sel_uv_x1 <- sel_uv_x2 <- sel_be_x1 <- sel_be_x2 <- logical(nsim)
  # Residual df of the refitted BE model (depends on how many predictors BE
  # retained: n - 3 if both x1 & x2 retained, n - 2 if only one retained,
  # n - 1 if neither retained). Needed because the stored BE SE comes from
  # whichever of these refits actually produced it.
  df_be <- numeric(nsim)

  for (i in seq_len(nsim)) {
    if (!is.null(progress_fn) && i %% 50 == 0)
      progress_fn(i / nsim, detail = paste0("Iteration ", i, " / ", nsim))

    D  <- mvrnorm(n = n, mu = c(0, 0), Sigma = Sig)
    y  <- as.numeric(as.matrix(D) %*% beta + rnorm(n, 0, sqrt(sigma2)))
    df <- data.frame(x1 = D[, 1], x2 = D[, 2], y = y)

    m1 <- lm(y ~ x1,      data = df)
    m2 <- lm(y ~ x2,      data = df)
    mf <- lm(y ~ x1 + x2, data = df, x = TRUE, y = TRUE)
    mb <- tryCatch(
      abe(mf, data = df, alpha = alpha, tau = Inf, type.test = "F", verbose = FALSE),
      error = function(e) mf
    )

    cm1 <- coef(summary(m1))
    cm2 <- coef(summary(m2))
    cmf <- coef(summary(mf))
    cmb <- coef(summary(mb))

    b1_uv[i]  <- cm1["x1", "Estimate"];  se1_uv[i] <- cm1["x1", "Std. Error"]
    p1_uv[i]  <- cm1["x1", "Pr(>|t|)"]
    b2_uv[i]  <- cm2["x2", "Estimate"];  se2_uv[i] <- cm2["x2", "Std. Error"]
    p2_uv[i]  <- cm2["x2", "Pr(>|t|)"]
    b1_mv[i]  <- cmf["x1", "Estimate"];  se1_mv[i] <- cmf["x1", "Std. Error"]
    p1_mv[i]  <- cmf["x1", "Pr(>|t|)"]
    b2_mv[i]  <- cmf["x2", "Estimate"];  se2_mv[i] <- cmf["x2", "Std. Error"]
    p2_mv[i]  <- cmf["x2", "Pr(>|t|)"]

    b1_be[i]  <- if ("x1" %in% rownames(cmb)) cmb["x1", "Estimate"]   else 0
    se1_be[i] <- if ("x1" %in% rownames(cmb)) cmb["x1", "Std. Error"] else 0
    b2_be[i]  <- if ("x2" %in% rownames(cmb)) cmb["x2", "Estimate"]   else 0
    se2_be[i] <- if ("x2" %in% rownames(cmb)) cmb["x2", "Std. Error"] else 0
    p1_be[i]  <- if ("x1" %in% rownames(cmb)) cmb["x1", "Pr(>|t|)"] else NA_real_
    p2_be[i]  <- if ("x2" %in% rownames(cmb)) cmb["x2", "Pr(>|t|)"] else NA_real_

    sel_uv_x1[i] <- p1_uv[i] < alpha
    sel_uv_x2[i] <- p2_uv[i] < alpha
    bc <- names(coef(mb))
    sel_be_x1[i] <- "x1" %in% bc
    sel_be_x2[i] <- "x2" %in% bc
    # df_be = n - (number of coefficients in the refitted BE model, incl. intercept)
    df_be[i] <- n - length(bc)
  }

  data.frame(
    b1_uv, b2_uv, b1_mv, b2_mv, b1_be, b2_be,
    se1_uv, se2_uv, se1_mv, se2_mv, se1_be, se2_be,
    p1_uv, p2_uv, p1_mv, p2_mv, p1_be, p2_be,
    sel_uv_x1, sel_uv_x2, sel_be_x1, sel_be_x2, df_be
  )
}

# Sensitivity simulation: selection rates AND mean coefficient estimates (UV, MV, BE)
run_sim_sens <- function(nsim_lite, n, beta, rho, sigma2, alpha, seed) {
  set.seed(seed)
  Sig <- matrix(c(1, rho, rho, 1), ncol = 2)
  uv_ok <- be_ok <- logical(nsim_lite)
  b1_uv <- b2_uv <- b1_mv <- b2_mv <- numeric(nsim_lite)
  b1_be <- b2_be <- rep(NA_real_, nsim_lite)

  p1_vec <- numeric(nsim_lite); p2_vec <- numeric(nsim_lite)
  sel_uv_x1 <- logical(nsim_lite); sel_uv_x2 <- logical(nsim_lite)
  sel_be_x1 <- logical(nsim_lite); sel_be_x2 <- logical(nsim_lite)

  for (i in seq_len(nsim_lite)) {
    D  <- mvrnorm(n = n, mu = c(0, 0), Sigma = Sig)
    y  <- as.numeric(as.matrix(D) %*% beta + rnorm(n, 0, sqrt(sigma2)))
    df <- data.frame(x1 = D[, 1], x2 = D[, 2], y = y)

    m1 <- lm(y ~ x1,      data = df)
    m2 <- lm(y ~ x2,      data = df)
    mf <- lm(y ~ x1 + x2, data = df, x = TRUE, y = TRUE)

    p1 <- coef(summary(m1))["x1", "Pr(>|t|)"]
    p2 <- coef(summary(m2))["x2", "Pr(>|t|)"]
    p1_vec[i] <- p1; p2_vec[i] <- p2
    b1_uv[i] <- coef(summary(m1))["x1", "Estimate"]
    b2_uv[i] <- coef(summary(m2))["x2", "Estimate"]
    b1_mv[i] <- coef(summary(mf))["x1", "Estimate"]
    b2_mv[i] <- coef(summary(mf))["x2", "Estimate"]

    uv_ok[i]     <- (p1 < alpha) & (p2 < alpha)
    sel_uv_x1[i] <- p1 < alpha
    sel_uv_x2[i] <- p2 < alpha

    mb <- tryCatch(
      abe(mf, data = df, alpha = alpha, tau = Inf, type.test = "F", verbose = FALSE),
      error = function(e) mf
    )
    bc       <- names(coef(mb))
    be_ok[i]     <- ("x1" %in% bc) & ("x2" %in% bc)
    sel_be_x1[i] <- "x1" %in% bc
    sel_be_x2[i] <- "x2" %in% bc
    if ("x1" %in% bc) b1_be[i] <- coef(mb)["x1"]
    if ("x2" %in% bc) b2_be[i] <- coef(mb)["x2"]
  }
  c(
    prop_uv         = mean(uv_ok),                prop_be         = mean(be_ok),
    prop_uv_x1_only = mean( sel_uv_x1 & !sel_uv_x2),
    prop_uv_x2_only = mean(!sel_uv_x1 &  sel_uv_x2),
    prop_be_x1_only = mean( sel_be_x1 & !sel_be_x2),
    prop_be_x2_only = mean(!sel_be_x1 &  sel_be_x2),
    mean_b1_uv = mean(b1_uv),              mean_b2_uv = mean(b2_uv),
    mean_b1_mv = mean(b1_mv),              mean_b2_mv = mean(b2_mv),
    mean_b1_be = mean(b1_be, na.rm = TRUE), mean_b2_be = mean(b2_be, na.rm = TRUE)
  )
}

# Selection frequency table ====================================================
summarise_selection <- function(res) {
  cats <- c("Both selected", "Only X1", "Only X2", "Neither")
  mk   <- function(x1, x2) c(
    "Both selected" = mean( x1 &  x2),
    "Only X1"       = mean( x1 & !x2),
    "Only X2"       = mean(!x1 &  x2),
    "Neither"       = mean(!x1 & !x2)
  )
  uv <- mk(res$sel_uv_x1, res$sel_uv_x2)
  be <- mk(res$sel_be_x1, res$sel_be_x2)
  data.frame(
    Method     = factor(
      rep(c("Backward Elimination", "Univariable Preselection"), each = 4),
      levels = c("Backward Elimination", "Univariable Preselection")
    ),
    Category   = factor(rep(cats, 2), levels = cats),
    Proportion = c(be, uv)
  )
}

# Bias / Variance / Coverage table (Univariable, BE, Multivariable) ============
bias_var_table <- function(res, beta, n) {
  qt_uv <- qt(0.975, df = n - 2)
  qt_mv <- qt(0.975, df = n - 3)
  # BE's residual df varies by iteration (n-3 both retained, n-2 one retained,
  # n-1 neither retained); use the per-iteration df actually used to fit each
  # BE refit.
  qt_be <- qt(0.975, df = res$df_be)

  cvg <- function(b, se, qt_df, true) {
    covered <- (b - qt_df * se <= true) & (b + qt_df * se >= true)
    covered[is.na(covered)] <- FALSE
    mean(covered)
  }

  data.frame(
    Model             = c("Univariable Preselection", "Univariable Preselection",
                          "Backward Elimination",     "Backward Elimination",
                          "Multivariable Full Model", "Multivariable Full Model"),
    Variable          = c("X1", "X2", "X1", "X2", "X1", "X2"),
    `True beta`       = c(beta[1], beta[2], beta[1], beta[2], beta[1], beta[2]),
    `Mean beta-hat`   = round(c(mean(res$b1_uv), mean(res$b2_uv),
                                mean(res$b1_be), mean(res$b2_be),
                                mean(res$b1_mv), mean(res$b2_mv)), 4),
    Bias              = round(c(mean(res$b1_uv) - beta[1],
                                mean(res$b2_uv) - beta[2],
                                mean(res$b1_be) - beta[1],
                                mean(res$b2_be) - beta[2],
                                mean(res$b1_mv) - beta[1],
                                mean(res$b2_mv) - beta[2]), 4),
    `Empirical SE`    = round(c(sd(res$b1_uv), sd(res$b2_uv),
                                sd(res$b1_be), sd(res$b2_be),
                                sd(res$b1_mv), sd(res$b2_mv)), 4),
    `Coverage 95% CI` = percent(c(
      cvg(res$b1_uv, res$se1_uv, qt_uv, beta[1]),
      cvg(res$b2_uv, res$se2_uv, qt_uv, beta[2]),
      cvg(res$b1_be, res$se1_be, qt_be, beta[1]),
      cvg(res$b2_be, res$se2_be, qt_be, beta[2]),
      cvg(res$b1_mv, res$se1_mv, qt_mv, beta[1]),
      cvg(res$b2_mv, res$se2_mv, qt_mv, beta[2])
    ), accuracy = 0.1),
    check.names = FALSE
  )
}

# OVB detection ================================================================
is_suppressor <- function(b1, b2, rho) {
  (b1 != 0) & (b2 != 0) & (rho != 0)
}
is_classic_suppressor <- function(b1, b2, rho) {
  (b1 * b2 * rho) < 0
}
# Do β1 and β2 have opposite signs? (determines direction of bias)
signs_differ <- function(b1, b2) {
  (b1 > 0 & b2 < 0) | (b1 < 0 & b2 > 0)
}
# Null-beta scenario: at least one coefficient is exactly zero
is_null_beta <- function(b1, b2) (b1 == 0) | (b2 == 0)

# Unified colour palette (UV=red, BE=green, MV=blue) ===========================
COL_UV  <- "#e74c3c"
COL_BE  <- "#27ae60"
COL_MV  <- "#2980b9"

COL_METHOD <- c("Univariable"          = COL_UV,
                "Multivariable"        = COL_MV,
                "Backward Elimination" = COL_BE)

COL_COEF   <- c("Univariable"          = COL_UV,
                "Multivariable"        = COL_MV,
                "Backward Elimination" = COL_BE)

# Ensures each facet's free y-scale spans at least `min_range` (used for the
# mean beta-hat sensitivity plots so small effect sizes remain readable).
min_yrange_df <- function(df, group_var, y_var, min_range = 0.1) {
  df %>%
    group_by(.data[[group_var]]) %>%
    summarise(ymin = min(.data[[y_var]], na.rm = TRUE),
              ymax = max(.data[[y_var]], na.rm = TRUE),
              .groups = "drop") %>%
    mutate(
      span = ymax - ymin,
      pad  = pmax(0, (min_range - span) / 2),
      ymin = ymin - pad,
      ymax = ymax + pad
    )
}

# Row background colours for HTML tables
ROW_BG <- c("Univariable Preselection" = "#fce8e6",
            "Backward Elimination"    = "#e8f5e9",
            "Multivariable Full Model" = "#e3f2fd")

# Shared ggplot2 theme =========================================================
app_theme <- function() {
  theme_bw(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(colour = "#555", size = 11),
      strip.background = element_rect(fill = "#2c3e50", colour = NA),
      strip.text       = element_text(colour = "white", face = "bold"),
      legend.position  = "bottom"
    )
}

# Info box helper
info_box <- function(border_col, ...) {
  div(
    style = paste0("background:#f8f9fa; border-left:4px solid ", border_col,
                   "; padding:12px 18px; border-radius:4px;"),
    h5("Interpretation", style = "margin-top:0;"),
    tags$ul(...)
  )
}

# HTML table helper with row colouring by Model ================================
coloured_table <- function(df) {
  model_col <- c(
    "Univariable Preselection" = "#fce8e6",
    "Backward Elimination"     = "#e8f5e9",
    "Multivariable Full Model"  = "#e3f2fd"
  )
  hdr <- tags$thead(
    tags$tr(lapply(names(df), function(h)
      tags$th(style = "background:#2c3e50; color:white; padding:6px 10px;", h)))
  )
  rows <- lapply(seq_len(nrow(df)), function(i) {
    mdl <- as.character(df[i, "Model"])
    bg  <- if (!is.na(model_col[mdl])) model_col[mdl] else "#ffffff"
    tags$tr(
      style = paste0("background:", bg, ";"),
      lapply(df[i, ], function(v)
        tags$td(style = "padding:5px 10px; border:1px solid #dee2e6;", v))
    )
  })
  tags$table(
    class = "table table-bordered table-sm",
    style = "font-size:14px; width:auto;",
    hdr,
    tags$tbody(rows)
  )
}

# Format a p-value for display
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) "< 0.001" else sprintf("= %.3f", p)
}

################################################################################
# UI
################################################################################
ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly", base_font = "Source Sans Pro"),

  titlePanel(div(
    h2("Consequences of Applying Univariable Preselection Compared to Backward Elimination",
       style = "margin-bottom:2px;"),
    p("Simulation illustrating why univariable pre-selection of covariates is in some situations problematic.",
      style = "color:#555; font-size:15px; margin-top:4px;")
  ), windowTitle = "univar_selection"),

  sidebarLayout(
    # Sidebar ==================================================================
    sidebarPanel(
      width = 3,

      uiOutput("suppression_alert"),

      h5("True model parameters", class = "text-primary"),
      helpText(HTML("Linear model:")),
      br(),
      helpText(HTML("Y = &beta;<sub>1</sub>X<sub>1</sub> + &beta;<sub>2</sub>X<sub>2</sub> + &epsilon;")),
      br(), br(),

      sliderInput("beta1", HTML("&beta;<sub>1</sub> (X1)"),
                  min = -1, max = 1, value = 0.5, step = 0.05),
      sliderInput("beta2", HTML("&beta;<sub>2</sub> (X2)"),
                  min = -1, max = 1, value = 0.3, step = 0.05),

      sliderInput("rho", HTML("Correlation of X1 and X2 (r)"),
                  min = -0.9, max = 0.9, value = -0.35, step = 0.05),
      div(style = "margin-top:-10px; margin-bottom:4px; font-size:11px; color:#888;
                   display:flex; justify-content:space-between; padding:0 6px;",
          span("negative"), span("|  0  |"), span("positive")),

      sliderInput("r2", HTML("R<sup>2</sup> (population)"),
                  min = 0.1, max = 0.99, value = 0.9, step = 0.01),

      hr(),
      h5("Study & simulation settings", class = "text-primary"),

      numericInput("n", "Sample size (n)",
                   value = 40, min = 10, max = 500, step = 5),
      sliderInput("alpha", HTML("&alpha; (selection threshold)"),
                  min = 0.01, max = 0.20, value = 0.05, step = 0.01),
      numericInput("nsim", "Number of simulations for main simulation",
                   value = 500, min = 100, max = 5000, step = 100),
      numericInput("nsim_sens",
                   "Number of simulations per grid point for the sensitivity analysis (tab 6)",
                   value = 200, min = 50, max = 1000, step = 50),
      numericInput("seed", "Random seed for reproducibility",
                   value = 123, min = 1, max = 99999, step = 1),

      hr(),

      actionButton("run", "RUN MAIN SIMULATION to populate tabs 1-5",
                   class = "btn-primary",
                   style = "width:100%; font-weight:600; margin-bottom:6px;"),
      actionButton("run_sensitivity", "RUN SENSITIVITY ANALYSIS to populate tab 6",
                   class = "btn-warning",
                   style = "width:100%; font-weight:600; margin-bottom:6px;"),
      actionButton("new_dataset", "NEW SINGLE DATASET to populate tab 7",
                   class = "btn-info",
                   style = "width:100%; font-weight:600;"),

      br(), br(),
      div(
        style = "background:#fff3cd; border-left:4px solid #f39c12;
                  padding:8px 12px; border-radius:4px; font-size:12px;",
        HTML("<b>Running times (approx.):</b><br>
              Main simulation (500 runs, n = 40): ~10 s<br>
              Sensitivity analysis: ~2 min
              (19 r values + 20 sample sizes, each with 200 runs)<br>
              <i>Increasing the number of simulations or the sample size increases runtime linearly.</i>")
      ),

      hr(),
      h5(HTML("Derived noise variance &sigma;<sup>2</sup>"), class = "text-muted"),
      verbatimTextOutput("sigma2_display")
    ),

    # Main panel with tabs =====================================================
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        selected = "Introduction",

        # Tab 0: Introduction ==================================================
        tabPanel(
          "Introduction",
          br(),
          div(style = "max-width:860px;",
              h4("What does this app do?", style = "color:#2c3e50;"),
              p("This app uses computer simulation to illustrate a common mistake in
               statistical analysis: testing each predictor variable individually
               (univariable selection) before building a multivariable model.
               We compare this approach with a principled alternative called
               backward elimination."),
              p("The core message is simple:"),
              div(
                style = "background:#eaf4fb; border-left:4px solid #2980b9;
                        padding:12px 18px; border-radius:4px; margin-bottom:16px;",
                HTML("<b>When two predictor variables are correlated with each other,
                    univariable screening can give misleading results.
                    A variable that is truly important may appear unimportant
                    when tested alone. Backward elimination avoids this problem
                    by always considering all variables together.</b>")
              ),
              div(
                style = "background:#fef9e7; border-left:4px solid #f39c12;
                        padding:12px 18px; border-radius:4px; margin-bottom:16px;",
                HTML("<b>Why does this matter in medicine and epidemiology?</b><br>
                    In virtually all real-world medical and epidemiological studies,
                    predictor variables are correlated with each other.
                    For example, age and BMI, smoking and alcohol use, or
                    blood pressure and kidney function rarely vary independently.
                    Correlated predictors are the norm, not the exception.
                    This means that the problems demonstrated in this app are not
                    rare edge cases. They occur routinely whenever
                    univariable screening is used as a variable selection step.")
              ),

              hr(),

              h4("Mathematical notation explained", style = "color:#2c3e50;"),
              p("The simulation is based on a simple linear regression model with
               two predictors. Here is what each symbol means:"),

              tags$table(
                class = "table table-bordered table-sm",
                style = "max-width:720px; font-size:14px;",
                tags$thead(tags$tr(
                  tags$th("Symbol"), tags$th("Name"), tags$th("Meaning in plain language")
                )),
                tags$tbody(
                  tags$tr(tags$td(HTML("Y")),
                          tags$td("Outcome variable"),
                          tags$td("The dependent variable you want to explain or predict (e.g. blood pressure).")),
                  tags$tr(tags$td(HTML("X<sub>1</sub>, X<sub>2</sub>")),
                          tags$td("Predictor variables"),
                          tags$td("Two independent variables (e.g. age and BMI). Both are truly associated with Y in the simulation.")),
                  tags$tr(tags$td(HTML("&beta;<sub>1</sub>, &beta;<sub>2</sub> (beta1, beta2)")),
                          tags$td("True regression coefficients"),
                          tags$td(HTML("The true effect sizes of X<sub>1</sub> and X<sub>2</sub> on Y. These are the values we are trying to recover by fitting a model. Set by the sliders."))),
                  tags$tr(tags$td(HTML("&beta;&#770; (beta-hat)")),
                          tags$td("Estimated regression coefficient"),
                          tags$td(HTML("The coefficient estimated from data by fitting a regression model. It will differ from the true &beta; due to sampling variability."))),
                  tags$tr(tags$td(HTML("&epsilon; (epsilon)")),
                          tags$td("Random error"),
                          tags$td(HTML("The part of Y that cannot be explained by X<sub>1</sub> and X<sub>2</sub>. It represents measurement error and unmeasured factors."))),
                  tags$tr(tags$td(HTML("r")),
                          tags$td(HTML("Correlation between X<sub>1</sub> and X<sub>2</sub>")),
                          tags$td("It quantifies how strongly the two predictors are linearly related. In the simulation it ranges from -0.9 (strong negative) to +0.9 (strong positive). Zero means independent.")),
                  tags$tr(tags$td(HTML("R<sup>2</sup>")),
                          tags$td("Coefficient of determination"),
                          tags$td(HTML("The proportion of the variance in Y explained by X<sub>1</sub> and X<sub>2</sub> together. For example, R<sup>2</sup> = 0.9 means that 90% of the variability in Y is explained by X<sub>1</sub> and X<sub>2</sub> together.
                           A higher R<sup>2</sup> means greater signal-to-noise ratio, which improves power: both methods are more likely to detect and retain the truly relevant variables.
                           A lower R<sup>2</sup> increases noise, reducing power for both approaches and making selection errors more frequent."))),
                  tags$tr(tags$td(HTML("&sigma;<sup>2</sup> (sigma squared)")),
                          tags$td("Error variance"),
                          tags$td(HTML("The variance of the random error &epsilon;. Calculated automatically from &beta;, r and R<sup>2</sup> to achieve the target R<sup>2</sup>."))),
                  tags$tr(tags$td(HTML("&alpha; (alpha)")),
                          tags$td("Significance threshold"),
                          tags$td("The p-value cutoff for selecting a variable. In practice, commonly set to 0.05.")),
                  tags$tr(tags$td("n"),
                          tags$td("Sample size"),
                          tags$td("Number of study participants in each simulated dataset.")),
                  tags$tr(tags$td("nsim"),
                          tags$td("Number of simulations"),
                          tags$td("It quantifies how many independent datasets are generated and analysed. More simulations give more stable results but take longer to compute.")),
                  tags$tr(tags$td(HTML("Bias")),
                          tags$td("Systematic error"),
                          tags$td(HTML("The mean difference between the estimated and true coefficient (mean(&beta;&#770;) &minus; &beta;). A biased variable selection method gives systematically wrong answers even with large samples."))),
                  tags$tr(tags$td("SE"),
                          tags$td("Standard error"),
                          tags$td(HTML("The empirical standard deviation of estimated coefficients (&beta;&#770;) across all simulated datasets. It measures precision."))),
                  tags$tr(tags$td("Coverage"),
                          tags$td("95% CI coverage"),
                          tags$td(HTML("The percentage of simulations in which the 95% confidence interval contains the true coefficient (&beta;). It should be approximately 95% for a correctly specified model.")))
                )
              ),

              hr(),

              h4("The two strategies compared", style = "color:#2c3e50;"),
              fluidRow(
                column(6, div(
                  style = "background:#fdf2f2; border:1px solid #e74c3c;
                            border-radius:6px; padding:14px; height:100%;",
                  h5(style = "color:#e74c3c; margin-top:0;", "Univariable Preselection"),
                  tags$ol(
                    tags$li(HTML("Derive a set of all predictors that are candidates for the multivariable model.")),
                    tags$li(HTML("Fit a separate univariable linear regression model for each predictor: lm(Y ~ X<sub>1</sub>) and lm(Y ~ X<sub>2</sub>).")),
                    tags$li(HTML("Keep only those predictors with p &lt; &alpha;.")),
                    tags$li("Derive the final model using only the selected predictors.")
                  ),
                  p(style = "color:#c0392b; font-size:13px;",
                    HTML("<b>Problem:</b> When X<sub>1</sub> and X<sub>2</sub> are
                          correlated, the univariable test for one variable does not
                          account for the other. This can make a truly important
                          variable look unimportant or vice versa."))
                )),
                column(6, div(
                  style = "background:#eaf4fb; border:1px solid #2980b9;
                            border-radius:6px; padding:14px; height:100%;",
                  h5(style = "color:#2980b9; margin-top:0;", "Backward Elimination"),
                  tags$ol(
                    tags$li(HTML("Derive a set of all predictors that are candidates for the multivariable model.")),
                    tags$li(HTML("Start with a reasonable full model containing all predictors: lm(Y ~ X<sub>1</sub> + X<sub>2</sub>).")),
                    tags$li("At each step, remove the predictor with the largest p-value, provided it exceeds the threshold &alpha;."),
                    tags$li(HTML("Stop when all remaining predictors are significant for the threshold &alpha;."))
                  ),
                  p(style = "color:#1a5276; font-size:13px;",
                    HTML("<b>Advantage:</b> All variables are always evaluated
                          in the context of the others. The correlation between
                          X<sub>1</sub> and X<sub>2</sub> is accounted for at every step."))
                ))
              ),

              hr(),

              h4("The masking phenomenon", style = "color:#2c3e50;"),
              p(HTML("A particularly striking failure of univariable screening occurs when
                   the signs of the two effect coefficients (&beta;<sub>1</sub>,
                   &beta;<sub>2</sub>) and the correlation (r) are arranged so
                   that their product is negative
                   [sign(&beta;<sub>1</sub> &times; &beta;<sub>2</sub> &times; r) &lt; 0].
                   In this <b>masking scenario</b>, one variable \"masks\"
                   the effect of the other in a univariable test, causing its
                   p-value to be much larger than it should be.")),
              p(HTML("A related scenario is <b>amplification</b> (ommitted variable bias), which occurs when the
                   product &beta;<sub>1</sub> &times; &beta;<sub>2</sub> &times; r &gt; 0.
                   Here the univariable estimates are biased <i>away</i> from zero, making
                   variables appear more important than they truly are.
                   Both scenarios arise whenever all three quantities (&beta;<sub>1</sub>,
                   &beta;<sub>2</sub>, r) are nonzero.")),
              p(HTML("The app detects these situations and gives a warning in the sidebar.")),
              p(HTML("<i>Example using the default parameters: &beta;<sub>1</sub> = 0.5, &beta;<sub>2</sub> = 0.3,
                    r = &minus;0.35.</i> The product of those three numbers is negative, so
                    X<sub>2</sub> is suppressed by the negative correlation.
                    In univariable screening X<sub>2</sub> is frequently missed,
                    while backward elimination keeps X<sub>2</sub> reliably.")),

              hr(),

              h4("How to use this app", style = "color:#2c3e50;"),
              tags$ol(
                tags$li(HTML("<b>Set parameters</b> in the left sidebar. Start with the defaults, which show a clear masking scenario.")),
                tags$li(HTML("<b>Click \"RUN MAIN SIMULATION\"</b> to generate results for tabs 1 to 5. A progress bar will appear in the corresponding tabs.")),
                tags$li(HTML("<b>Click \"RUN SENSITIVITY ANALYSIS\"</b> to populate tab 6. This may take 2&ndash;3 minutes. A progress bar will appear in tab 6.")),
                tags$li(HTML("<b>Click \"NEW SINGLE DATASET\"</b> to populate tab 7, which shows one concrete example.")),
                tags$li("Explore the tabs in order; each one adds a layer of understanding.")
              ),

              hr(),

              h4("Guide to the tabs", style = "color:#2c3e50;"),
              tags$table(
                class = "table table-striped table-bordered table-sm",
                style = "font-size:14px;",
                tags$thead(tags$tr(
                  tags$th("Tab"), tags$th("What you see"), tags$th("Key question answered")
                )),
                tags$tbody(
                  tags$tr(
                    tags$td("1. Selection Frequency"),
                    tags$td(HTML("Bar chart showing how often each method selects both, one, or neither variable (X<sub>1</sub>, X<sub>2</sub>).")),
                    tags$td("Does the variable selection method find the variables that truly belong in the model?")),
                  tags$tr(
                    tags$td("2. Coefficient Estimates"),
                    tags$td(HTML("Density plots of estimated coefficients (&beta;&#770;) for all three methods.")),
                    tags$td("Is the estimated effect size correct on average, i.e., unbiased?")),
                  tags$tr(
                    tags$td("3. p-Value Distributions"),
                    tags$td(HTML("Histograms of p-values for X<sub>1</sub> and X<sub>2</sub> for all three methods.")),
                    tags$td("Are p-values well-calibrated, or distorted by correlation?")),
                  tags$tr(
                    tags$td("4. Bias & Variance"),
                    tags$td(HTML("Table with bias, empirical SE, and CI coverage of the estimated coefficients (&beta;&#770;) for all three approaches.")),
                    tags$td("How large is the systematic error, and do confidence intervals work?")),
                  tags$tr(
                    tags$td("5. Scatterplot of Estimates"),
                    tags$td(HTML("Univariable versus multivariable estimated coefficients (&beta;&#770;) per simulated dataset.")),
                    tags$td("How different are the two estimates, and in which direction?")),
                  tags$tr(
                    tags$td("6. Sensitivity Analysis"),
                    tags$td(HTML("Line plots across correlation (r) and sample size (n) visualizing correct selection rate and mean estimated coefficients (&beta;&#770;) for all three approaches.")),
                    tags$td("How does performance change as correlation or sample size varies?")),
                  tags$tr(
                    tags$td("7. Single Dataset"),
                    tags$td("One concrete simulated dataset with all model outputs."),
                    tags$td("What does a single analysis look like? Walk through backward elimination step by step."))
                )
              ),

              hr(),

              h4("Expected running times", style = "color:#2c3e50;"),
              div(
                style = "background:#f8f9fa; border-left:4px solid #7f8c8d;
                        padding:12px 18px; border-radius:4px;",
                tags$ul(
                  tags$li(HTML("<b>Main simulation</b> with 500 simulations and n = 40: approximately 10 seconds.")),
                  tags$li(HTML("<b>Main simulation</b> with 2000 simulations and n = 200: approximately 2 minutes.")),
                  tags$li(HTML("<b>Sensitivity analysis</b> with 200 simulations per grid point:
                              approximately 2&ndash;3 minutes.
                              The r-sweep evaluates 19 values of r (from &minus;0.9 to +0.9 in steps of 0.1);
                              the n-sweep evaluates 20 sample sizes (from 5 to 500).
                              Both sweeps use all remaining sidebar settings.")),
                  tags$li(HTML("<i>Runtime increases roughly linearly with the number of simulations and the sample size. Use fewer simulations for quick exploration, more for stable final results.</i>"))
                )
              ),
              hr(),
              h4("Impressum", style = "color:#2c3e50;"),
              div(
                style = "background:#f8f9fa; border:1px solid #dee2e6;
                        padding:16px 20px; border-radius:6px; font-size:14px;",
                h5("Authors", style = "margin-top:0;"),
                tags$p(HTML(
                  "Daniela Dunkler &amp; Theresa Ullmann<br>
                 Medical University of Vienna, Center for Medical Data Science,
                 Institute of Clinical Biometrics<br>
                 Spitalgasse 23, 1090 Vienna, Austria<br>
                 E-mail: daniela . dunkler @ meduniwien . ac . at"
                )),
                h5("Licence"),
                tags$p(HTML(
                  "Released under the
                 <a href='https://www.gnu.org/licenses/gpl-3.0.html' target='_blank'>
                 GNU General Public License v3.0 (GPL-3)</a>."
                )),
                h5("Citation"),
                tags$p(HTML("Ullmann T., Heinze G., Kappenberg F., Henrion M., Sauerbrei W., Collins G.,
                Leonhardt C.-S., Nold M., and Dunkler D. for TG2 of the STRATOS initiative
                (2026). <i>The Problem with Univariable Selection in Regression
                Modelling&mdash;and What To Do Instead.</i> In preparation.")),
                h5("Source code"),
                tags$p(HTML('This app is implemented in the R Shiny web application framework and is deployed through
               RStudio\'s webservice <a href="https://www.shinyapps.io/">shinyapps.io</a>.
               The code is available in the accompanying
               <a href="https://github.com/biometrician/UnivarSelection/">Github</a> repository.')),
                h5("Disclaimer"),
                tags$p("For educational purposes only. The authors accept no liability for use of results in clinical or policy decisions."),
                h5("Version"),
                tags$p(HTML("Version 1.0 &mdash; [RELEASE September 2026]"), class = "text-muted")
              )

          ) # end max-width div
        ), # end Introduction tab

        # Tab 1: Selection Frequency ===========================================
        tabPanel(
          "1. Selection Frequency",
          br(),
          uiOutput("both_zero_warn_sel"),
          plotOutput("sel_plot", height = "400px"),
          br(),
          div(
            style = "background:#f8f9fa; border-left:4px solid #2c3e50;
                      padding:12px 18px; border-radius:4px;",
            h5("Interpretation", style = "margin-top:0;"),
            tags$ul(
              uiOutput("sel_interp_bullet1"),
              uiOutput("sel_interp_bullet2"),
              uiOutput("sel_interp_bullet3"),
              uiOutput("sel_interp_bullet4"),
              uiOutput("sel_r2_bullet")
            )
          ),
          br(),
          uiOutput("num_results_header"),
          tableOutput("result_table")
        ),

        # Tab 2: Coefficient Estimates =========================================
        tabPanel(
          "2. Coefficient Estimates",
          br(),
          uiOutput("both_zero_warn_coef"),
          plotOutput("coef_density_plot", height = "500px"),
          br(),
          div(
            style = "background:#f8f9fa; border-left:4px solid #e74c3c;
                      padding:12px 18px; border-radius:4px;",
            h5("Interpretation", style = "margin-top:0;"),
            tags$ul(
              tags$li(HTML("The dashed vertical line represents the true coefficient (&beta;).")),
              uiOutput("coef_uv_bullet"),
              uiOutput("coef_mv_bullet"),
              uiOutput("coef_be_bullet"),
              uiOutput("coef_r2_bullet"),
              uiOutput("coef_zero_r_bullet"),
              uiOutput("coef_scenario_bullet")
            )
          )
        ),

        # Tab 3: p-Value Distributions =========================================
        tabPanel(
          "3. p-Value Distributions",
          br(),
          uiOutput("both_zero_warn_pval"),
          plotOutput("pval_plot", height = "620px"),
          br(),
          div(
            style = "background:#f8f9fa; border-left:4px solid #27ae60;
                      padding:12px 18px; border-radius:4px;",
            h5("Interpretation", style = "margin-top:0;"),
            tags$ul(
              tags$li(HTML("Under H<sub>1</sub> (effect truly present), p-values should pile up near 0, i.e. indicating high power.")),
              tags$li(HTML("Distorted (i.e. right-shifted) univariable p-values mean that the correlation between X<sub>1</sub> and X<sub>2</sub> is masking the true effect.")),
              tags$li(HTML("The vertical dashed line represents the significance threshold (&alpha;). The proportion of p-values to the left of this line equals the empirical power.")),
              uiOutput("pval_be_bullet"),
              uiOutput("pval_r2_bullet"),
              uiOutput("pval_zero_r_bullet"),
              uiOutput("pval_scenario_bullet")
            )
          )
        ),

        # Tab 4: Bias & Variance Table =========================================
        tabPanel(
          "4. Bias & Variance",
          br(),
          uiOutput("both_zero_warn_bv"),
          h5("Bias, Empirical Standard Error (SE), and 95% Confidence Interval (CI) Coverage",
             class = "text-primary"),
          div(
            style = "background:#eaf4fb; border-left:4px solid #2980b9;
                      padding:10px 16px; border-radius:4px; margin-bottom:14px; font-size:16px;",
            HTML(
              "<b>Bias</b> is a systematic deviation of the estimated coefficient from
               the true coefficient (mean(&beta;&#770;) &minus; &beta;). A biased estimator
               gives the wrong answer on average, even with a very large sample. An unbiased estimator has bias close to 0.<br>
               <b>Empirical SE</b> is the standard deviation of the estimated coefficients (&beta;&#770;) computed
               across all simulated datasets. It quantifies the variability (imprecision) of
               the estimates from one dataset to the next.
               A smaller SE means the estimates are more consistent, but consistency around a wrong value (when bias is present) does not help &mdash; it only makes the wrong answer more reproducible.
               In the multivariable model, the SE is generally larger than in the univariable model because adjusting for a correlated predictor increases uncertainty; this is the bias&ndash;variance trade-off.
               Backward elimination inherits the SE of the full model when both variables are retained, but shows a reduced (and misleadingly small) SE when a variable is dropped, because the 0&ndash;imputation reduces variability artificially.<br>
               <b>Coverage 95% CI</b> is the percentage of simulations in which the
               95% confidence interval for the coefficient (&beta;) actually contains the true coefficient (&beta;).
               For a correctly specified model this should be approximately 95%.
               Bias causes under-coverage, i.e., the estimated coefficients miss the truth more than 5% of the time."
            )
          ),
          br(),
          p(HTML("Note: For Backward Elimination, when a variable is dropped, its estimate is set to 0 and its SE to 0.
                  The mean, bias, SE and coverage are computed across <i>all</i> simulations (unconditional), so they reflect both estimation uncertainty and selection uncertainty."),
            style = "font-size:13px; color:#555;"),
          br(),
          uiOutput("bv_table_ui"),
          br(),
          div(
            style = "background:#f8f9fa; border-left:4px solid #8e44ad;
                      padding:12px 18px; border-radius:4px;",
            h5("Interpretation", style = "margin-top:0;"),
            tags$ul(
              tags$li(HTML("<b>True beta</b>: The true coefficient (&beta;) value set by the slider. The estimator should recover this value.")),
              tags$li(HTML("<b>Mean beta-hat (&beta;&#770;)</b>: The mean estimated coefficient across all simulations.")),
              tags$li(HTML("A lower R<sup>2</sup> (more noise) leads to a larger empirical SE for all methods, reflecting reduced precision. The bias is not directly affected by R<sup>2</sup>, but low R<sup>2</sup> reduces power, increasing the chance that backward elimination drops a variable and thereby increases its unconditional bias.")),
              uiOutput("bv_scenario_bullet")
            )
          )
        ),

        # Tab 5: Scatterplot of Estimates ======================================
        tabPanel(
          "5. Scatterplot of Estimates",
          br(),
          uiOutput("both_zero_warn_scatter"),
          plotOutput("scatter_plot", height = "460px"),
          br(),
          div(
            style = "background:#f8f9fa; border-left:4px solid #e67e22;
                      padding:12px 18px; border-radius:4px;",
            h5("Interpretation", style = "margin-top:0;"),
            tags$ul(
              tags$li("Each point represents one simulated dataset."),
              tags$li(HTML("The red dashed diagonal represents perfect agreement between univariable and multivariable estimated coefficients (&beta;&#770;).")),
              tags$li(HTML("The blue dotted lines visualize the true coefficients (&beta;<sub>1</sub> and &beta;<sub>2</sub>).")),
              uiOutput("scatter_scenario_bullet")
            )
          )
        ),

        # Tab 6: Sensitivity Analysis ==========================================
        tabPanel(
          "6. Sensitivity Analysis",
          br(),
          uiOutput("both_zero_warn_sens"),
          fluidRow(
            column(6, uiOutput("hdr_rho_sel"), plotOutput("sens_rho_plot", height = "340px")),
            column(6, uiOutput("hdr_n_sel"),   plotOutput("sens_n_plot",   height = "340px"))
          ),
          br(),
          uiOutput("sens_interp_row1"),
          br(),
          hr(),
          fluidRow(
            column(6, uiOutput("hdr_rho_coef"), plotOutput("sens_rho_coef_plot", height = "360px")),
            column(6, uiOutput("hdr_n_coef"),   plotOutput("sens_n_coef_plot",   height = "360px"))
          ),
          br(),
          div(
            style = "background:#f8f9fa; border-left:4px solid #16a085;
                      padding:12px 18px; border-radius:4px;",
            h5("Interpretation", style = "margin-top:0;"),
            tags$ul(
              tags$li(HTML("<b>Left plot:</b> Mean estimated coefficient (&beta;&#770;) from univariable (<span style='color:#e74c3c;'>red</span>), backward elimination (<span style='color:#27ae60;'>green</span>), and multivariable (<span style='color:#2980b9;'>blue</span>, reference) models at each value of r. The horizontal dashed lines give the true coefficients (&beta;<sub>1</sub> and &beta;<sub>2</sub>). Deviation from the dashed line represents bias. The dotted vertical line marks the currently selected r.")),
              tags$li(HTML("<b>Right plot:</b> Same as the left plot, but across sample sizes n. The multivariable &beta;&#770; (blue) converges to the truth as n grows. Backward elimination (green) converges toward the multivariable estimate as n increases (because at large n BE almost always retains both variables). The dotted vertical line marks the currently selected n.")),
              tags$li(HTML("When X<sub>1</sub> and X<sub>2</sub> are uncorrelated (r = 0), all three approaches coincide (no omitted-variable bias). As |r| increases, the univariable estimates diverge from the truth while the multivariable estimates remain unbiased.")),
              uiOutput("sens_coef_scenario_bullet")
            )
          )
        ),

        # Tab 7: Single Dataset Inspector ======================================
        tabPanel(
          "7. Single Dataset",
          br(),
          p("One simulated dataset generated with the current parameters.
             Click 'NEW SINGLE DATASET' in the sidebar to draw a fresh example.",
            style = "color:#555; font-size:13px;"),
          br(),
          fluidRow(
            column(4, plotOutput("single_x1y",  height = "240px")),
            column(4, plotOutput("single_x2y",  height = "240px")),
            column(4, plotOutput("single_x1x2", height = "240px"))
          ),
          fluidRow(
            column(4, plotOutput("single_dens_x1", height = "180px")),
            column(4, plotOutput("single_dens_x2", height = "180px")),
            column(4)
          ),
          br(),
          uiOutput("single_plot_interp"),

          hr(),

          # Row 3: Combined unadjusted & adjusted prediction plots (one legend each)
          h5("Unadjusted & Adjusted Prediction Plots", class = "text-primary", style = "margin-top:10px;"),
          p(HTML("Each panel overlays the unadjusted (univariable, ignoring the other predictor) and the
                  adjusted (multivariable, other predictor held fixed at its mean) regression line
                  for the same predictor, with a single shared legend."),
            style = "color:#555; font-size:13px;"),
          fluidRow(
            column(6, plotOutput("single_pred_x1", height = "300px")),
            column(6, plotOutput("single_pred_x2", height = "300px"))
          ),
          uiOutput("single_pred_interp"),

          # ---------------------------------------------------------------------
          # The two blocks below (Partial Residual Plots and Added-Variable
          # Plots) are intentionally disabled
          # The matching server calculations/plots are commented out as well
          # (see "Tab 7" server section) so that nothing is visible, but the
          # code remains available in case it should be re-enabled later.
          # ---------------------------------------------------------------------
          # # Row 5: Partial residual (component + residual) plots
          # h5("Partial Residual Plots (Component + Residual)", class = "text-primary"),
          # p(HTML("Shows X<sub>j</sub> vs. (fitted component of X<sub>j</sub> + model residual).
          #         The slope of the line equals the multivariable coefficient &beta;&#770;<sub>j</sub>.
          #         <b>Use for diagnostics:</b> The scatter around the line indicates residual variability after adjustment."),
          #   style = "color:#555; font-size:13px;"),
          # fluidRow(
          #   column(6, plotOutput("single_partres_x1", height = "240px")),
          #   column(6, plotOutput("single_partres_x2", height = "240px"))
          # ),
          # br(),
          # # Row 6: Added-variable (partial regression) plots
          # h5("Added-Variable Plots (Partial Regression)", class = "text-primary"),
          # p(HTML("For each predictor X<sub>j</sub>: plot (residuals of X<sub>j</sub> on all other predictors)
          #         vs. (residuals of Y on all other predictors). The slope of the line equals
          #         the multivariable coefficient &beta;&#770;<sub>j</sub> exactly.
          #         <b>Use for exact coefficient interpretation:</b> the plot shows the unique contribution
          #         of X<sub>j</sub> after removing the shared variation with the other predictor.
          #         Influential points visible here have a strong effect on the estimated coefficient."),
          #   style = "color:#555; font-size:13px;"),
          # fluidRow(
          #   column(6, plotOutput("single_avplot_x1", height = "240px")),
          #   column(6, plotOutput("single_avplot_x2", height = "240px"))
          # ),

          hr(),

          p(tags$b("Output from R"), style = "font-size:15px;"),
          fluidRow(
            column(4,
                   tags$b("Univariable model: lm(y ~ x1)", style = "color:#e74c3c;"),
                   verbatimTextOutput("single_m1")),
            column(4,
                   tags$b("Univariable model: lm(y ~ x2)", style = "color:#e67e22;"),
                   verbatimTextOutput("single_m2")),
            column(4,
                   tags$b("Full multivariable model: lm(y ~ x1 + x2)", style = "color:#2980b9;"),
                   verbatimTextOutput("single_mf"))
          ),
          uiOutput("single_model_interp"),
          br(),
          p(tags$b("Output from R"), style = "font-size:15px;"),
          tags$b("Backward Elimination: step-by-step output", style = "color:#27ae60;"),
          verbatimTextOutput("single_be"),
          uiOutput("single_be_interp"),

          hr(),

          h5("Summary: Coefficient Estimates from All Models", class = "text-primary"),
          uiOutput("single_results_table_ui")
        )
      ) # end tabsetPanel
    ) # end mainPanel
  ) # end sidebarLayout
) # end fluidPage

################################################################################
# Server
################################################################################
server <- function(input, output, session) {

  # Derived sigma^2 =============================================================
  sigma2_val <- reactive({
    req(input$beta1, input$beta2, input$rho, input$r2)
    calc_sigma2(c(input$beta1, input$beta2), input$rho, input$r2)
  })

  output$sigma2_display <- renderText({
    sprintf("sigma^2 = %.4f\n(noise variance to\nachieve target R^2)", sigma2_val())
  })

  # OVB checks (reactives used throughout) =====================================
  ovb_active <- reactive({ is_suppressor(input$beta1, input$beta2, input$rho) })
  classic_suppressor <- reactive({ is_classic_suppressor(input$beta1, input$beta2, input$rho) })
  betas_opposite_sign <- reactive({ signs_differ(input$beta1, input$beta2) })
  # Null-beta scenario reactives
  null_beta  <- reactive({ is_null_beta(input$beta1, input$beta2) })
  b1_is_zero <- reactive({ input$beta1 == 0 })
  b2_is_zero <- reactive({ input$beta2 == 0 })
  both_zero  <- reactive({ input$beta1 == 0 & input$beta2 == 0 })
  # r = 0 with both betas nonzero: no OVB but MV has smaller residual variance
  zero_r_both_active <- reactive({
    abs(input$rho) < 1e-9 & input$beta1 != 0 & input$beta2 != 0
  })

  # Suppression annotation (sidebar) amber bg for classic suppressor ==========
  output$suppression_alert <- renderUI({
    # Null-beta scenario overrides OVB alert
    if (null_beta()) {
      if (both_zero()) {
        lbl <- "Null model: both coefficients are zero"
        msg <- HTML("Both &beta;<sub>1</sub> = 0 and &beta;<sub>2</sub> = 0.
                     Neither X<sub>1</sub> nor X<sub>2</sub> has any true association with Y.
                     The correct variable selection decision is to drop both predictors.
                     At significance level &alpha;, both methods will include a variable by chance
                     approximately 100&alpha;% of the time (type I error rate).")
        bg  <- "#e8f5e9"; col <- "#27ae60"
      } else if (b1_is_zero()) {
        if (abs(input$rho) > 0 & input$beta2 != 0) {
          lbl <- "Null coefficient scenario: β₁ = 0, r ≠ 0"
          msg <- HTML("&beta;<sub>1</sub> = 0: X<sub>1</sub> has no true effect on Y.
                       The correct decision is to drop X<sub>1</sub> and retain X<sub>2</sub>.
                       <b>Because r &ne; 0</b>, X<sub>1</sub> is correlated with X<sub>2</sub>,
                       which does have a true effect. Therefore, univariable screening will
                       show a spurious association for X<sub>1</sub>
                       (false positive), leading to over-selection of X<sub>1</sub>.
                       Backward elimination, starting from the full model, correctly
                       distinguishes the null from the active variable.")
        } else {
          lbl <- "Null coefficient scenario: β₁ = 0, r = 0"
          msg <- HTML("&beta;<sub>1</sub> = 0: X<sub>1</sub> has no true effect on Y.
                       Because r = 0, there is no correlation-induced bias.
                       Both methods will select X<sub>1</sub> at approximately the nominal rate &alpha;
                       (type I error), and retain X<sub>2</sub> at the appropriate power.")
        }
        bg  <- "#fff8e1"; col <- "#f9a825"
      } else {
        # b2_is_zero
        if (abs(input$rho) > 0 & input$beta1 != 0) {
          lbl <- "Null coefficient scenario: β₂ = 0, r ≠ 0"
          msg <- HTML("&beta;<sub>2</sub> = 0: X<sub>2</sub> has no true effect on Y.
                       The correct decision is to drop X<sub>2</sub> and retain X<sub>1</sub>.
                       <b>Because r ≠ 0</b>, X<sub>2</sub> is correlated with X<sub>1</sub>,
                       which does have a true effect. Therefore, univariable screening will
                       show a spurious association for X<sub>2</sub>
                       (false positive), leading to over-selection of X<sub>2</sub>.
                       Backward elimination correctly identifies the active variable.")
        } else {
          lbl <- "Null coefficient scenario: β₂ = 0, r = 0"
          msg <- HTML("&beta;<sub>2</sub> = 0: X<sub>2</sub> has no true effect on Y.
                       Because r = 0, there is no correlation-induced bias.
                       Both methods will select X<sub>2</sub> at approximately the nominal rate &alpha;
                       (type I error), and retain X<sub>1</sub> at the appropriate power.")
        }
        bg  <- "#fff8e1"; col <- "#f9a825"
      }
      return(div(
        style = paste0("background:", bg, "; border-left:4px solid ", col,
                       "; padding:10px 14px; border-radius:4px; margin-bottom:14px;"),
        tags$b(lbl),
        tags$p(style = "margin:5px 0 0; font-size:12px;", msg)
      ))
    }

    if (!ovb_active()) return(NULL)

    if (classic_suppressor()) {
      bg  <- "#fff3cd"   # amber, same hue as the left border (#f39c12)
      col <- "#f39c12"
      lbl <- "Masking scenario"
      msg <- HTML("One variable masks the effect of the other because the two
          effect coefficients (&beta;<sub>1</sub>, &beta;<sub>2</sub>) and
          the correlation (r) are arranged so that their product is
          <i>negative</i> [sign(&beta;<sub>1</sub> &times; &beta;<sub>2</sub>
          &times; r) &lt; 0]. At least one univariable p-value is inflated relative to the multivariable model.
          Univariable preselection is likely to drop it even if
          it truly belongs in the model. See Introduction for details.")
    } else {
      bg  <- "#f8f0fc"
      col <- "#8e44ad"
      lbl <- "Amplification scenario"
      msg <- HTML("The signs of the two
          effect coefficients (&beta;<sub>1</sub>, &beta;<sub>2</sub>) and
          the correlation (r) are arranged so that their
          product is <i>positive</i>.
          Univariable estimates are biased <i>away</i> from zero, making
          variables appear more important than they truly are.
          Both variables will tend to be over-selected. Estimation is
          biased even though selection itself may appear satisfactory.
          See Introduction for details.")
    }
    div(
      style = paste0("background:", bg, "; border-left:4px solid ", col,
                     "; padding:10px 14px; border-radius:4px; margin-bottom:14px;"),
      tags$b(lbl),
      tags$p(style = "margin:5px 0 0; font-size:12px;", msg)
    )
  })

  # Main simulation =============================================================
  sim_res <- eventReactive(input$run, {
    req(sigma2_val())
    withProgress(message = "Running main simulation: ", value = 0, {
      run_simulation_full(
        nsim     = input$nsim,
        n        = input$n,
        beta     = c(input$beta1, input$beta2),
        rho      = input$rho,
        sigma2   = sigma2_val(),
        alpha    = input$alpha,
        seed     = input$seed,
        progress_fn = function(p, detail) setProgress(p, detail = detail)
      )
    })
  })

  sel_df <- reactive({ req(sim_res()); summarise_selection(sim_res()) })

  # Tab 1: Selection frequency ==================================================
  output$sel_plot <- renderPlot({
    req(sel_df())
    cat_col <- c("Both selected" = "#27ae60", "Only X1" = "#e67e22",
                 "Only X2"       = "#e74c3c", "Neither" = "#7f8c8d")
    ggplot(sel_df(), aes(x = Category, y = Proportion, fill = Category)) +
      geom_col(width = 0.6, colour = "white") +
      geom_text(aes(label = percent(Proportion, accuracy = 0.1)),
                vjust = -0.4, size = 3.8, fontface = "bold") +
      facet_wrap(~ Method, ncol = 2) +
      scale_y_continuous(labels = percent_format(accuracy = 1),
                         limits = c(0, 1.1), expand = c(0, 0)) +
      scale_fill_manual(values = cat_col, guide = "none") +
      labs(
        title    = "Variable Selection Frequency",
        subtitle = sprintf(
          "nsim = %d  |  n = %d  |  beta1 = %.2f  |  beta2 = %.2f  |  r = %.2f  |  R2 = %.2f  |  alpha = %.2f",
          input$nsim, input$n, input$beta1, input$beta2,
          input$rho, input$r2, input$alpha),
        x = NULL, y = "Proportion of simulations"
      ) +
      app_theme() +
      theme(legend.position = "none", panel.grid.major.x = element_blank())
  })

  output$sel_interp_bullet4 <- renderUI({
    if (null_beta()) {
      if (both_zero()) {
        tags$li(HTML(
          "<b>Null model scenario</b> (&beta;<sub>1</sub> = &beta;<sub>2</sub> = 0):
           Neither variable has a true effect on Y.
           The <b>correct</b> decision is <b>'Neither selected'</b>.
           The proportion labelled 'Both selected' represents the type&nbsp;I error rate
           (selecting variables that have no effect)."
        ))
      } else if (b1_is_zero()) {
        if (abs(input$rho) > 1e-9 & input$beta2 != 0) {
          tags$li(HTML(
            "<b>Null coefficient scenario</b> (&beta;<sub>1</sub> = 0, r &ne; 0):
             X<sub>1</sub> has no true effect, so the <b>correct</b> decision is
             <b>'Only X<sub>2</sub> selected'</b>.
             Because r &ne; 0, univariable screening picks up a spurious association
             for X<sub>1</sub> (false positive), so <b>'Only X<sub>1</sub>'</b> and
             <b>'Both selected'</b> represent selection errors.
             Backward elimination by evaluating both variables jointly
             is better at identifying that X<sub>1</sub> has no independent effect."
          ))
        } else {
          tags$li(HTML(
            "<b>Null coefficient scenario</b> (&beta;<sub>1</sub> = 0, r = 0):
             X<sub>1</sub> has no true effect and is uncorrelated with X<sub>2</sub>.
             The <b>correct</b> decision is <b>'Only X<sub>2</sub> selected'</b>.
             Both methods should select X<sub>1</sub> at approximately the nominal type&nbsp;I
             error rate (&alpha; &asymp; 5%)."
          ))
        }
      } else {
        if (abs(input$rho) > 1e-9 & input$beta1 != 0) {
          tags$li(HTML(
            "<b>Null coefficient scenario</b> (&beta;<sub>2</sub> = 0, r &ne; 0):
             X<sub>2</sub> has no true effect, so the <b>correct</b> decision is
             <b>'Only X<sub>1</sub> selected'</b>.
             Because r &ne; 0, univariable screening picks up a spurious association
             for X<sub>2</sub> (false positive).
             Backward elimination is better at identifying that X<sub>2</sub> has
             no independent effect."
          ))
        } else {
          tags$li(HTML(
            "<b>Null coefficient scenario</b> (&beta;<sub>2</sub> = 0, r = 0):
             X<sub>2</sub> has no true effect and is uncorrelated with X<sub>1</sub>.
             The <b>correct</b> decision is <b>'Only X<sub>1</sub> selected'</b>.
             Both methods should select X<sub>2</sub> at approximately the nominal type&nbsp;I
             error rate (&alpha; &asymp; 5%)."
          ))
        }
      }
    } else if (!ovb_active()) {
      tags$li(HTML(
        "In this scenario r = 0 or at least one coefficient is 0, so there is
         no omitted-variable bias. Univariable and multivariable estimates coincide
         and both variable selection methods perform equivalently."
      ))
    } else if (classic_suppressor()) {
      tags$li(HTML(
        "This is a <b>masking scenario</b>
         [sign(&beta;<sub>1</sub> &times; &beta;<sub>2</sub> &times; r) &lt; 0]:
         one variable masks the effect of the other in a univariable test,
         inflating its p-value. Univariable preselection is particularly likely
         to drop a variable that genuinely belongs in the model."
      ))
    } else {
      tags$li(HTML(
        "This is an <b>amplification scenario</b>, i.e., the signs of the two coefficients (&beta;<sub>1</sub>, &beta;<sub>2</sub>) and
        the correlation (r) are arranged so that their product is <i>positive</i>
         [sign(&beta;<sub>1</sub> &times; &beta;<sub>2</sub> &times; r) &gt; 0].
         Univariable estimates are biased <i>away</i> from zero. Both variables
         tend to appear even more significant than they truly are. Variable selection
         may look satisfactory, but the estimated coefficients are inflated. The
         multivariable model correctly adjusts for the shared variance."
      ))
    }
  })

  output$num_results_header <- renderUI({
    req(sim_res())
    h5(paste0("Numerical results (", input$nsim, " simulations)"), class = "text-primary")
  })

  output$result_table <- renderTable({
    req(sel_df())
    sel_df() |>
      mutate(Proportion = percent(Proportion, accuracy = 0.1)) |>
      pivot_wider(names_from = Method, values_from = Proportion) |>
      rename(`Selection outcome` = Category)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # Tab 2: Coefficient density plots (UV + MV + BE) ============================
  output$coef_density_plot <- renderPlot({
    req(sim_res())
    validate(
      need(!both_zero(),
           "Coefficient density plot is not meaningful when both β1 = 0 and β2 = 0 (Y has no systematic variation). Please set at least one coefficient to a non-zero value.")
    )
    res <- sim_res()
    df  <- bind_rows(
      data.frame(Variable = "X1", Method = "Univariable",          bhat = res$b1_uv),
      data.frame(Variable = "X1", Method = "Multivariable",        bhat = res$b1_mv),
      data.frame(Variable = "X1", Method = "Backward Elimination", bhat = res$b1_be),
      data.frame(Variable = "X2", Method = "Univariable",          bhat = res$b2_uv),
      data.frame(Variable = "X2", Method = "Multivariable",        bhat = res$b2_mv),
      data.frame(Variable = "X2", Method = "Backward Elimination", bhat = res$b2_be)
    ) |>
      mutate(
        Variable = recode(Variable,
                          X1 = paste0("X1   (true \u03b2\u2081 = ", input$beta1, ")"),
                          X2 = paste0("X2   (true \u03b2\u2082 = ", input$beta2, ")"))
      )
    refs <- data.frame(
      Variable  = c(paste0("X1   (true \u03b2\u2081 = ", input$beta1, ")"),
                    paste0("X2   (true \u03b2\u2082 = ", input$beta2, ")")),
      true_beta = c(input$beta1, input$beta2)
    )
    ggplot(df, aes(x = bhat, fill = Method, colour = Method)) +
      geom_density(alpha = 0.25, linewidth = 0.9) +
      geom_vline(data = refs, aes(xintercept = true_beta),
                 linetype = "dashed", colour = "black", linewidth = 1.1,
                 inherit.aes = FALSE) +
      facet_wrap(~ Variable, scales = "free", ncol = 2) +
      scale_fill_manual(values = COL_METHOD) +
      scale_colour_manual(values = COL_METHOD) +
      labs(
        title    = "Distribution of Coefficient Estimates",
        subtitle = "Black dashed line = true \u03b2.  Shift from dashed line = omitted-variable bias.",
        x        = "Estimated \u03b2\u0302",
        y        = "Density",
        fill     = "Model", colour = "Model"
      ) +
      app_theme() +
      theme(strip.text = element_text(colour = "white", face = "bold", size = 15))
  })

  # Scenario bullet for Tab 2
  output$coef_scenario_bullet <- renderUI({
    if (null_beta()) {
      if (both_zero()) {
        tags$li(HTML(
          "<b>Null model</b> (&beta;<sub>1</sub> = &beta;<sub>2</sub> = 0): All coefficient estimates are centred around 0 (the true value). The density curves for all three methods are nearly identical and centred on 0."
        ))
      } else if (b1_is_zero()) {
        tags$li(HTML(
          "<b>Null coefficient (&beta;<sub>1</sub> = 0):</b> The multivariable estimate for X<sub>1</sub> is unbiased (centred on 0). If r &ne; 0, the univariable estimate for X<sub>1</sub> is biased away from 0 (E[&beta;&#770;<sub>1,uv</sub>] = 0 + &beta;<sub>2</sub> &times; r = &beta;<sub>2</sub>r), so the univariable density for X<sub>1</sub> is shifted. For X<sub>2</sub>, the univariable estimate is unbiased since &beta;<sub>1</sub> = 0 (E[&beta;&#770;<sub>2,uv</sub>] = &beta;<sub>2</sub> + 0 &times; r = &beta;<sub>2</sub>). Backward elimination tends to drop X<sub>1</sub> from the full model, so its green density for X<sub>1</sub> is a mixture at 0 (dropped) and the true multivariable distribution (retained)."
        ))
      } else {
        tags$li(HTML(
          "<b>Null coefficient (&beta;<sub>2</sub> = 0):</b> The multivariable estimate for X<sub>2</sub> is unbiased (centred on 0). If r &ne; 0, the univariable estimate for X<sub>2</sub> is biased (E[&beta;&#770;<sub>2,uv</sub>] = 0 + &beta;<sub>1</sub> &times; r = &beta;<sub>1</sub>r), so the univariable density for X<sub>2</sub> is shifted. For X<sub>1</sub>, the univariable estimate is unbiased since &beta;<sub>2</sub> = 0. Backward elimination tends to drop X<sub>2</sub> from the full model."
        ))
      }
    } else if (!ovb_active()) {
      tags$li(HTML(
        "No omitted-variable bias is present (r = 0 or a coefficient is 0). All three density curves overlap and are centred on the true &beta;. All three methods are equivalent."
      ))
    } else if (classic_suppressor()) {
      if (betas_opposite_sign()) {
        tags$li(HTML(
          "<b>Suppressor scenario with opposite-sign coefficients:</b> Univariable estimates for X<sub>1</sub> and X<sub>2</sub> are biased in opposite directions relative to their true values (one is over-estimated, the other is under-estimated). The multivariable model remains unbiased for both."
        ))
      } else {
        tags$li(HTML(
          "<b>Suppressor scenario with same-sign coefficients:</b> Both univariable estimates are biased toward zero (i.e. attenuated). The effect of each variable is masked by the negative correlation with the other. The multivariable model correctly recovers both true coefficients."
        ))
      }
    } else {
      if (betas_opposite_sign()) {
        tags$li(HTML(
          "<b>Amplification scenario with opposite-sign coefficients:</b> Univariable estimates for X<sub>1</sub> and X<sub>2</sub> are biased in opposite directions (one is inflated positively, the other negatively). The multivariable model remains unbiased."
        ))
      } else {
        tags$li(HTML(
          "<b>Amplification scenario with same-sign coefficients:</b> Both univariable estimates are biased away from zero (i.e. inflated). Each variable appears more strongly associated with Y than it truly is. The multivariable model correctly recovers both true coefficients."
        ))
      }
    }
  })

  # R² bullet for Tab 2
  output$coef_r2_bullet <- renderUI({
    tags$li(HTML(paste0(
      "The current R<sup>2</sup> = ", input$r2,
      ". A lower R<sup>2</sup> produces wider density curves (larger SE) for all methods, reflecting greater estimation uncertainty. The position of the univariable peak relative to the true &beta; (bias) is determined by the correlation between X<sub>1</sub> and X<sub>2</sub> (r) and the coefficients, not by R<sup>2</sup>."
    )))
  })

  # Tab 3: p-value histograms (UV + MV + BE) ====================================
  output$pval_plot <- renderPlot({
    req(sim_res())
    validate(
      need(!both_zero(),
           "p-value distribution is not meaningful when both β1 = 0 and β2 = 0. All p-values follow a uniform distribution. Please set at least one coefficient to a non-zero value.")
    )
    res <- sim_res()
    # For BE, only include simulations where the variable was retained (p not NA)
    df  <- bind_rows(
      data.frame(Variable = "X1", Method = "Univariable",          pval = res$p1_uv),
      data.frame(Variable = "X1", Method = "Multivariable",        pval = res$p1_mv),
      data.frame(Variable = "X1", Method = "Backward Elimination", pval = res$p1_be),
      data.frame(Variable = "X2", Method = "Univariable",          pval = res$p2_uv),
      data.frame(Variable = "X2", Method = "Multivariable",        pval = res$p2_mv),
      data.frame(Variable = "X2", Method = "Backward Elimination", pval = res$p2_be)
    ) |> filter(!is.na(pval)) |>
      mutate(Method = factor(Method,
                             levels = c("Univariable", "Multivariable", "Backward Elimination")))

    ggplot(df, aes(x = pval, fill = Method)) +
      geom_histogram(bins = 25, alpha = 0.75, position = "identity",
                     colour = "white", linewidth = 0.2) +
      geom_vline(xintercept = input$alpha, linetype = "dashed",
                 colour = "black", linewidth = 0.9) +
      facet_grid(Method ~ Variable, scales = "free_y") +
      scale_fill_manual(values = COL_METHOD, guide = "none") +
      labs(
        title    = "p-Value Distributions",
        subtitle = paste0("alpha = ", input$alpha,
                          " (dashed).  Pile-up near 0 = high power.  ",
                          "Right-skew = inflated p-values due to omitted-variable bias."),
        x = "p-value", y = "Count"
      ) +
      app_theme() +
      theme(legend.position = "none",
            strip.text = element_text(colour = "white", face = "bold", size = 15))
  })

  output$pval_r2_bullet <- renderUI({
    tags$li(HTML(paste0(
      "R<sup>2</sup> = ", input$r2,
      ". A lower R<sup>2</sup> shifts all p-value histograms to the right, which corresponds to less power. Fewer simulations reach significance, and the pile-up near 0 becomes less pronounced for all three methods."
    )))
  })

  output$pval_scenario_bullet <- renderUI({
    if (null_beta()) {
      if (both_zero()) {
        tags$li(HTML(
          "<b>Null model</b>: Neither variable has any true effect. Under H<sub>0</sub>, p-values follow a uniform distribution on [0, 1] for both methods. The proportion of p-values below &alpha; equals the type&nbsp;I error rate (should be &asymp; &alpha; = 5%)."
        ))
      } else if (b1_is_zero()) {
        if (abs(input$rho) > 1e-9 & input$beta2 != 0) {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>1</sub> = 0, r &ne; 0):</b> For X<sub>1</sub> (the null variable), the univariable p-values pile up near 0 (false positives due to correlation with X<sub>2</sub>), while the multivariable p-values are roughly uniform (correctly non-significant). For X<sub>2</sub> (the true predictor), both univariable and multivariable p-values pile up near 0 (high power, correctly significant)."
          ))
        } else {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>1</sub> = 0, r = 0):</b> For X<sub>1</sub>, p-values are approximately uniformly distributed for both methods (correct type&nbsp;I error control). For X<sub>2</sub>, p-values pile up near 0 (high power)."
          ))
        }
      } else {
        if (abs(input$rho) > 1e-9 & input$beta1 != 0) {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>2</sub> = 0, r &ne; 0):</b> For X<sub>2</sub> (the null variable), the univariable p-values pile up near 0 (false positives due to correlation with X<sub>1</sub>), while the multivariable p-values are roughly uniform. For X<sub>1</sub> (the true predictor), both approaches show high power."
          ))
        } else {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>2</sub> = 0, r = 0):</b> For X<sub>2</sub>, p-values are approximately uniformly distributed (correct type&nbsp;I error control). For X<sub>1</sub>, p-values pile up near 0 (high power)."
          ))
        }
      }
    } else if (!ovb_active()) {
      tags$li(HTML(
        "No omitted-variable bias: the univariable and multivariable p-value distributions are identical. No distortion from correlation."
      ))
    } else if (classic_suppressor()) {
      tags$li(HTML(
        "<b>Suppressor scenario:</b> The univariable p-value distribution for the suppressed variable is right-shifted (inflated) compared to the multivariable distribution. Univariable p-values can exceed &alpha; even though the variable is truly important. The backward elimination distribution mirrors the multivariable distribution when the variable is retained, but is empty for simulations where it was dropped."
      ))
    } else {
      tags$li(HTML(
        "<b>Amplification scenario:</b> The univariable p-value distributions are shifted toward zero (deflated) compared to the multivariable distributions, because variables appear more strongly associated than they truly are. Both univariable and backward elimination methods have high power (selection rates), but the coefficient estimates are biased."
      ))
    }
  })

  # Tab 4: Bias-variance-coverage table (coloured HTML) ========================
  output$bv_table_ui <- renderUI({
    if (both_zero()) return(NULL)
    req(sim_res())
    df <- bias_var_table(sim_res(), c(input$beta1, input$beta2), input$n)
    coloured_table(df)
  })

  # Scenario-specific bullet for Tab 4
  output$bv_scenario_bullet <- renderUI({
    if (null_beta()) {
      if (both_zero()) {
        tags$li(HTML(
          "<b>Null model</b>: Both true coefficients are 0. Mean beta-hat is close to 0 for all methods (unbiased). Coverage is approximately 95% for all methods. This represents ideal behaviour under the null hypothesis."
        ))
      } else if (b1_is_zero()) {
        if (abs(input$rho) > 1e-9 & input$beta2 != 0) {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>1</sub> = 0, r &ne; 0):</b> For X<sub>1</sub>, the multivariable mean beta-hat is close to 0 (unbiased), but the univariable mean beta-hat is biased toward &beta;<sub>2</sub> &times; r (spurious association). Coverage for the univariable CI of &beta;<sub>1</sub> is below 95% because the CI is centred on the wrong value. For X<sub>2</sub>, both univariable and multivariable estimates are unbiased."
          ))
        } else {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>1</sub> = 0, r = 0):</b> For X<sub>1</sub>, all methods give unbiased estimates near 0. Coverage is approximately 95% for all. For X<sub>2</sub>, all three methods are unbiased."
          ))
        }
      } else {
        if (abs(input$rho) > 1e-9 & input$beta1 != 0) {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>2</sub> = 0, r &ne; 0):</b> For X<sub>2</sub>, the multivariable mean beta-hat is close to 0 (unbiased), but the univariable mean beta-hat is biased toward &beta;<sub>1</sub> &times; r. Coverage for the univariable CI of &beta;<sub>2</sub> is below 95%. For X<sub>1</sub>, both methods are unbiased."
          ))
        } else {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>2</sub> = 0, r = 0):</b> For X<sub>2</sub>, all methods give unbiased estimates near 0 with approximately 95% coverage. For X<sub>1</sub>, all methods are unbiased."
          ))
        }
      }
    } else if (!ovb_active()) {
      tags$li(HTML(
        "No omitted-variable bias: bias is close to 0 for all three methods. Coverage of the 95% CI is approximately 95% for all. The Empirical SE is slightly larger for the multivariable model due to adjusting for a correlated predictor."
      ))
    } else if (classic_suppressor()) {
      tags$li(HTML(
        "<b>Suppressor scenario:</b> Bias in the univariable model leads to under-coverage. The 95% CI in the univariable model misses the true coefficient (&beta;) more often than 5% of the time, because the CI is centred on a biased estimate. The multivariable model maintains coverage close to 95%."
      ))
    } else {
      tags$li(HTML(
        "<b>Amplification scenario:</b> The univariable model over-estimates the coefficients (bias away from zero). The 95% CI is shifted away from the true &beta; and also shows under-coverage, just in the opposite direction compared to the masking scenario. The multivariable model maintains approximately 95% coverage."
      ))
    }
  })

  # Tab 5: Scatterplot of Estimates =============================================
  output$scatter_plot <- renderPlot({
    req(sim_res())
    validate(
      need(!both_zero(),
           "Scatterplot is not meaningful when both β1 = 0 and β2 = 0 (all estimates are at zero). Please set at least one coefficient to a non-zero value.")
    )
    res <- sim_res()
    df  <- bind_rows(
      data.frame(Variable = "X1", b_uv = res$b1_uv, b_mv = res$b1_mv, true_b = input$beta1),
      data.frame(Variable = "X2", b_uv = res$b2_uv, b_mv = res$b2_mv, true_b = input$beta2)
    )
    ax_min <- min(df$b_uv, df$b_mv)
    ax_max <- max(df$b_uv, df$b_mv)
    pad    <- (ax_max - ax_min) * 0.04
    ggplot(df, aes(x = b_uv, y = b_mv)) +
      geom_point(alpha = 0.25, size = 1.4, colour = "#444") +
      geom_abline(slope = 1, intercept = 0, colour = COL_UV,
                  linewidth = 1.0, linetype = "dashed") +
      geom_vline(aes(xintercept = true_b), linetype = "dotted",
                 colour = COL_MV, linewidth = 0.9) +
      geom_hline(aes(yintercept = true_b), linetype = "dotted",
                 colour = COL_MV, linewidth = 0.9) +
      facet_wrap(~ Variable) +
      coord_equal(xlim = c(ax_min - pad, ax_max + pad),
                  ylim = c(ax_min - pad, ax_max + pad)) +
      labs(
        title    = "Univariable versus Multivariable \u03b2\u0302 per Simulation",
        subtitle = "Red dashed = perfect agreement.  Blue dotted = true \u03b2.  Deviation = omitted-variable bias.",
        x = "\u03b2\u0302 from univariable model",
        y = "\u03b2\u0302 from multivariable model"
      ) +
      app_theme() +
      theme(legend.position = "none",
            strip.text = element_text(colour = "white", face = "bold", size = 15))
  })

  # Scenario bullet for Tab 5
  output$scatter_scenario_bullet <- renderUI({
    if (null_beta()) {
      if (both_zero()) {
        tags$li(HTML(
          "<b>Null model:</b> Both univariable and multivariable estimates are centred on 0. Points cluster around the origin along the diagonal. No systematic deviation; both estimates are unbiased."
        ))
      } else if (b1_is_zero()) {
        if (abs(input$rho) > 1e-9 & input$beta2 != 0) {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>1</sub> = 0, r &ne; 0):</b> For X<sub>1</sub>, the univariable estimates are systematically shifted away from 0 (false association), while multivariable estimates cluster near 0 (true value). The cloud of points for X<sub>1</sub> is displaced horizontally from the diagonal &mdash; this is the spurious univariable association caused by r &ne; 0. For X<sub>2</sub>, points lie along the diagonal (both methods unbiased)."
          ))
        } else {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>1</sub> = 0, r = 0):</b> For X<sub>1</sub>, points cluster near the origin (both methods estimate ~0, correctly). For X<sub>2</sub>, points lie along the diagonal."
          ))
        }
      } else {
        if (abs(input$rho) > 1e-9 & input$beta1 != 0) {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>2</sub> = 0, r &ne; 0):</b> For X<sub>2</sub>, the univariable estimates are systematically shifted from 0 (spurious association), while multivariable estimates cluster near 0. The horizontal displacement for X<sub>2</sub> reveals the false positive induced by r &ne; 0. For X<sub>1</sub>, points lie along the diagonal."
          ))
        } else {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>2</sub> = 0, r = 0):</b> For X<sub>2</sub>, points cluster near the origin (both methods correctly estimate ~0). For X<sub>1</sub>, points lie along the diagonal."
          ))
        }
      }
    } else if (!ovb_active()) {
      tags$li(HTML(
        "No omitted-variable bias: points cluster tightly along the red diagonal, indicating that univariable and multivariable estimates agree. There is no systematic deviation."
      ))
    } else if (classic_suppressor()) {
      if (betas_opposite_sign()) {
        tags$li(HTML(
          "<b>Suppressor scenario (opposite-sign coefficients):</b> Points are systematically displaced from the diagonal in opposite directions for X<sub>1</sub> and X<sub>2</sub>. One variable's univariable estimate is too small, the other's is too large relative to the multivariable estimate. The larger |r|, the greater the displacement."
        ))
      } else {
        tags$li(HTML(
          "<b>Suppressor scenario (same-sign coefficients):</b> Points are systematically displaced from the diagonal in the same direction for both X<sub>1</sub> and X<sub>2</sub>: univariable estimates are smaller (closer to zero) than the corresponding multivariable estimates, because the negative correlation attenuates both effects. The larger |r|, the greater the attenuation."
        ))
      }
    } else {
      if (betas_opposite_sign()) {
        tags$li(HTML(
          "<b>Amplification scenario (opposite-sign coefficients):</b> Points are systematically displaced from the diagonal in opposite directions for X<sub>1</sub> and X<sub>2</sub>. Univariable estimates diverge from the multivariable estimates, each in the direction away from zero for its respective variable."
        ))
      } else {
        tags$li(HTML(
          "<b>Amplification scenario (same-sign coefficients):</b> Points are systematically displaced from the diagonal in the same direction: univariable estimates are larger (further from zero) than the multivariable estimates for both variables. The positive correlation amplifies both univariable effects. The larger |r|, the greater the inflation."
        ))
      }
    }
  })

  # Tab 6: Sensitivity analysis =================================================
  sens_res <- eventReactive(input$run_sensitivity, {
    req(sigma2_val())
    rho_grid <- seq(-0.9, 0.9, by = 0.1)
    n_grid   <- round(seq(10, 500, length.out = 20))
    n_steps  <- length(rho_grid) + length(n_grid)

    withProgress(message = "Running sensitivity analysis: ", value = 0, {
      rho_df <- bind_rows(lapply(seq_along(rho_grid), function(j) {
        rho_j <- rho_grid[j]
        s2_j  <- calc_sigma2(c(input$beta1, input$beta2), rho_j, input$r2)
        r     <- run_sim_sens(
          nsim_lite = input$nsim_sens, n = input$n,
          beta = c(input$beta1, input$beta2),
          rho = rho_j, sigma2 = s2_j, alpha = input$alpha, seed = input$seed
        )
        incProgress(1 / n_steps, detail = sprintf("r sweep: r = %.2f", rho_j))
        data.frame(rho = rho_j,
                   Univariable = unname(r["prop_uv"]),
                   `Backward Elimination` = unname(r["prop_be"]),
                   prop_uv_x1_only = unname(r["prop_uv_x1_only"]),
                   prop_uv_x2_only = unname(r["prop_uv_x2_only"]),
                   prop_be_x1_only = unname(r["prop_be_x1_only"]),
                   prop_be_x2_only = unname(r["prop_be_x2_only"]),
                   b1_uv = unname(r["mean_b1_uv"]), b2_uv = unname(r["mean_b2_uv"]),
                   b1_mv = unname(r["mean_b1_mv"]), b2_mv = unname(r["mean_b2_mv"]),
                   b1_be = unname(r["mean_b1_be"]), b2_be = unname(r["mean_b2_be"]),
                   check.names = FALSE)
      }))

      n_df <- bind_rows(lapply(seq_along(n_grid), function(j) {
        n_j <- n_grid[j]
        r   <- run_sim_sens(
          nsim_lite = input$nsim_sens, n = n_j,
          beta = c(input$beta1, input$beta2),
          rho = input$rho, sigma2 = sigma2_val(), alpha = input$alpha, seed = input$seed
        )
        incProgress(1 / n_steps, detail = sprintf("n sweep: n = %d", n_j))
        data.frame(n = n_j,
                   Univariable = unname(r["prop_uv"]),
                   `Backward Elimination` = unname(r["prop_be"]),
                   prop_uv_x1_only = unname(r["prop_uv_x1_only"]),
                   prop_uv_x2_only = unname(r["prop_uv_x2_only"]),
                   prop_be_x1_only = unname(r["prop_be_x1_only"]),
                   prop_be_x2_only = unname(r["prop_be_x2_only"]),
                   b1_uv = unname(r["mean_b1_uv"]), b2_uv = unname(r["mean_b2_uv"]),
                   b1_mv = unname(r["mean_b1_mv"]), b2_mv = unname(r["mean_b2_mv"]),
                   b1_be = unname(r["mean_b1_be"]), b2_be = unname(r["mean_b2_be"]),
                   check.names = FALSE)
      }))
      list(rho = rho_df, n = n_df)
    })
  })

  # Reactive plot headers for Tab 6
  output$hdr_rho_sel <- renderUI({
    lbl <- if (null_beta()) {
      if (both_zero()) "Type I error rate (both selected) vs. r"
      else if (b1_is_zero()) "Power (X2) and type I error (X1) vs. r"
      else "Power (X1) and type I error (X2) vs. r"
    } else "Correct selection rate vs. r"
    h5(HTML(paste0(lbl, " &mdash; for sample size n = ", input$n)), class = "text-primary")
  })
  output$hdr_n_sel <- renderUI({
    lbl <- if (null_beta()) {
      if (both_zero()) "Type I error rate (both selected) vs. n"
      else if (b1_is_zero()) "Power (X2) and type I error (X1) vs. n"
      else "Power (X1) and type I error (X2) vs. n"
    } else "Correct selection rate vs. n"
    h5(HTML(paste0(lbl, " &mdash; for correlation r = ", input$rho)), class = "text-primary")
  })
  output$hdr_rho_coef <- renderUI({ h5(HTML(paste0("Mean &beta;&#770; vs. r &mdash; for sample size n = ", input$n)), class = "text-primary") })
  output$hdr_n_coef   <- renderUI({ h5(HTML(paste0("Mean &beta;&#770; vs. n &mdash; for correlation r = ", input$rho)), class = "text-primary") })

  output$sens_interp_row1 <- renderUI({
    li_last <- if (null_beta()) {
      if (both_zero()) {
        tags$li(HTML(
          "<b>Null model:</b> Both variables have no true effect. The correct selection rate for 'both selected' is 0 &mdash; selecting either variable is an error. The displayed proportion represents the joint type&nbsp;I error rate. Both methods should show low 'both selected' rates, ideally near &alpha;<sup>2</sup>."
        ))
      } else if (b1_is_zero()) {
        tags$li(HTML(
          "<b>Null coefficient (&beta;<sub>1</sub> = 0):</b> Note that the plot shows P(both X<sub>1</sub> &amp; X<sub>2</sub> selected), but here the <b>correct</b> selection is <b>only X<sub>2</sub></b>. Selecting both is an error (false positive for X<sub>1</sub>). If r &ne; 0, univariable screening will spuriously include X<sub>1</sub> more often than &alpha;; backward elimination &mdash; jointly evaluating both &mdash; provides better type&nbsp;I error control for the null variable."
        ))
      } else {
        tags$li(HTML(
          "<b>Null coefficient (&beta;<sub>2</sub> = 0):</b> Note that the plot shows P(both X<sub>1</sub> &amp; X<sub>2</sub> selected), but here the <b>correct</b> selection is <b>only X<sub>1</sub></b>. Selecting both is an error (false positive for X<sub>2</sub>). If r &ne; 0, univariable screening will spuriously include X<sub>2</sub>; backward elimination provides better type&nbsp;I error control."
        ))
      }
    } else if (classic_suppressor()) {
      tags$li(HTML(
        "In this <b>masking scenario</b>, the performance gap persists
         even at large sample sizes (n), because additional data cannot fix a systematic bias.
         It only makes the biased estimate more precise."
      ))
    } else if (ovb_active()) {
      tags$li(HTML(
        "In this <b>amplification scenario</b>, both variables tend to appear significant
         even in univariable testing. Both methods show high correct selection rates.
         However, the estimated coefficients are inflated (see the row two plots below),
         and the multivariable model is still needed for unbiased estimation."
      ))
    } else {
      tags$li(HTML(
        "No omitted-variable bias is present (r = 0 or at least one coefficient is 0).
         Both methods perform equivalently. Performance improves with sample size (n)
         for both methods."
      ))
    }
    div(
      style = "background:#f8f9fa; border-left:4px solid #16a085;
                padding:12px 18px; border-radius:4px;",
      h5("Interpretation", style = "margin-top:0;"),
      tags$ul(
        uiOutput("sens_yaxis_bullet"),
        tags$li(HTML(paste0(
          "<b>Left plot:</b> Correlation (r) swept from &minus;0.9 to +0.9 by steps of 0.1.
           All other parameters remain fixed (n = ", input$n, ").
           The dotted vertical line gives the currently selected r."
        ))),
        tags$li(HTML(
          "Selection based on univariate preselection is affected by the correlation between X<sub>1</sub> and X<sub>2</sub> (r). Selection based on backward elimination is more robust to the correlation between X<sub>1</sub> and X<sub>2</sub> (r)."
        )),
        tags$li(HTML(paste0(
          "<b>Right plot:</b> Sample size (n) increased from 10 to 500 (20 evenly spaced values).
           All other parameters remain fixed (r = ", input$rho, ").
           The upper x-axis shows n per variable (npv = n / 2).
           The dotted vertical line visualizes the currently selected sample size (n)."
        ))),
        li_last
      )
    )
  })

  # Scenario-specific bullet for Tab 6 row 2 interpretation ====================
  output$sens_coef_scenario_bullet <- renderUI({
    if (null_beta()) {
      if (both_zero()) {
        tags$li(HTML(
          "<b>Null model:</b> All three curves lie near 0 (the true value) for both variables, regardless of r or n. No bias; all methods are correct on average."
        ))
      } else if (b1_is_zero()) {
        if (abs(input$rho) > 1e-9 & input$beta2 != 0) {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>1</sub> = 0, r &ne; 0):</b> For X<sub>1</sub>, the univariable curve (red) is biased away from 0 and tracks &beta;<sub>2</sub> &times; r as r changes; this is the spurious association formula E[&beta;&#770;<sub>1,uv</sub>] = &beta;<sub>2</sub>r. The multivariable curve (blue) stays near 0 (unbiased). Backward elimination (green) converges to the multivariable curve at large n, because it drops X<sub>1</sub> more reliably as sample size grows. For X<sub>2</sub>, all three methods are unbiased."
          ))
        } else {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>1</sub> = 0, r = 0):</b> All curves for X<sub>1</sub> lie near 0 (no bias). For X<sub>2</sub>, all methods are unbiased."
          ))
        }
      } else {
        if (abs(input$rho) > 1e-9 & input$beta1 != 0) {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>2</sub> = 0, r &ne; 0):</b> For X<sub>2</sub>, the univariable curve (red) is biased and tracks &beta;<sub>1</sub> &times; r (the spurious association formula E[&beta;&#770;<sub>2,uv</sub>] = &beta;<sub>1</sub>r). The multivariable curve (blue) stays near 0. Backward elimination (green) converges to the multivariable curve at large n. For X<sub>1</sub>, all methods are unbiased."
          ))
        } else {
          tags$li(HTML(
            "<b>Null coefficient (&beta;<sub>2</sub> = 0, r = 0):</b> All curves for X<sub>2</sub> lie near 0 (no bias). For X<sub>1</sub>, all methods are unbiased."
          ))
        }
      }
    } else if (!ovb_active()) {
      tags$li(HTML(
        "No omitted-variable bias: all three curves lie on top of each other and on the true &beta; dashed line at all values of r and n. There is no direction of bias to describe."
      ))
    } else if (classic_suppressor()) {
      if (betas_opposite_sign()) {
        tags$li(HTML(
          "<b>Suppressor scenario (opposite-sign coefficients):</b> The direction of bias in the univariable estimates differs for X<sub>1</sub> and X<sub>2</sub>: one is over-estimated and the other is under-estimated relative to the truth. This is because the bias formula E[&beta;&#770;<sub>1,uv</sub>] = &beta;<sub>1</sub> + &beta;<sub>2</sub>r gives a bias with a sign that depends on whether &beta;<sub>2</sub>r is positive or negative, and vice versa for X<sub>2</sub>. With opposite-sign coefficients and r of the same sign as one but opposite to the other, the biases go in opposite directions. Backward elimination (green) tracks between the univariable and multivariable curves depending on how often it retains each variable."
        ))
      } else {
        tags$li(HTML(
          "<b>Suppressor scenario (same-sign coefficients):</b> Both univariable estimates are biased in the same direction &mdash; toward zero. The omitted-variable bias formula E[&beta;&#770;<sub>1,uv</sub>] = &beta;<sub>1</sub> + &beta;<sub>2</sub>r gives a negative bias for both variables when &beta;<sub>1</sub>, &beta;<sub>2</sub>, and r lead to same-direction attenuation. Backward elimination (green) converges to the multivariable (blue) curve as n increases, because at large n BE nearly always retains both variables. At small n, BE is intermediate between univariable and multivariable estimates."
        ))
      }
    } else {
      if (betas_opposite_sign()) {
        tags$li(HTML(
          "<b>Amplification scenario (opposite-sign coefficients):</b> The univariable estimates for X<sub>1</sub> and X<sub>2</sub> are inflated in opposite directions (one more positive, one more negative than the true value). Backward elimination (green) typically retains both variables and tracks close to the multivariable (blue) curve, especially at larger n."
        ))
      } else {
        tags$li(HTML(
          "<b>Amplification scenario (same-sign coefficients):</b> Both univariable estimates are inflated away from zero in the same direction. The larger |r|, the greater the inflation. Backward elimination (green) converges toward the multivariable (blue) estimate as n grows, because it retains both variables more reliably at larger sample sizes."
        ))
      }
    }
  })

  # Helpers for sensitivity plots
  sens_x_theme <- function() {
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  }
  empty_sec_axis <- sec_axis(~ ., name = " ", labels = NULL, breaks = NULL)
  npv_axis_theme <- function() {
    theme(
      axis.title.x.top = element_text(margin = margin(b = 14, t = 6)),
      axis.text.x.top  = element_text(margin = margin(b = 8))
    )
  }

  sel_y_label <- reactive({
    if (b1_is_zero() & !b2_is_zero())
      "P(X2 selected & X1 not selected)"
    else if (b2_is_zero() & !b1_is_zero())
      "P(X1 selected & X2 not selected)"
    else if (both_zero())
      "P(both X1 & X2 selected) [error rate]"
    else
      "P(both X1 & X2 selected)"
  })

  # Prepare sensitivity df for selection plots, extracting correct metric
  sens_sel_df_rho <- reactive({
    req(sens_res())
    raw <- sens_res()$rho
    if (b1_is_zero() & !b2_is_zero()) {
      data.frame(rho = raw$rho,
                 Univariable = raw[["prop_uv_x2_only"]],
                 `Backward Elimination` = raw[["prop_be_x2_only"]],
                 check.names = FALSE)
    } else if (b2_is_zero() & !b1_is_zero()) {
      data.frame(rho = raw$rho,
                 Univariable = raw[["prop_uv_x1_only"]],
                 `Backward Elimination` = raw[["prop_be_x1_only"]],
                 check.names = FALSE)
    } else {
      raw[, c("rho", "Univariable", "Backward Elimination")]
    }
  })
  sens_sel_df_n <- reactive({
    req(sens_res())
    raw <- sens_res()$n
    if (b1_is_zero() & !b2_is_zero()) {
      data.frame(n = raw$n,
                 Univariable = raw[["prop_uv_x2_only"]],
                 `Backward Elimination` = raw[["prop_be_x2_only"]],
                 check.names = FALSE)
    } else if (b2_is_zero() & !b1_is_zero()) {
      data.frame(n = raw$n,
                 Univariable = raw[["prop_uv_x1_only"]],
                 `Backward Elimination` = raw[["prop_be_x1_only"]],
                 check.names = FALSE)
    } else {
      raw[, c("n", "Univariable", "Backward Elimination")]
    }
  })

  # Sensitivity: selection rate vs rho (BE=green) --------------------------------
  output$sens_rho_plot <- renderPlot({
    req(sens_sel_df_rho())
    df <- sens_sel_df_rho() |>
      pivot_longer(c("Univariable", "Backward Elimination"),
                   names_to = "Method", values_to = "Rate")
    ggplot(df, aes(x = rho, y = Rate, colour = Method, group = Method)) +
      geom_line(linewidth = 1.3) + geom_point(size = 2.8) +
      geom_vline(xintercept = input$rho, linetype = "dotted",
                 colour = "#555", linewidth = 0.8) +
      annotate("text", x = input$rho + 0.04, y = 0.04,
               label = paste0("r = ", input$rho),
               size = 3.5, colour = "#555", hjust = 0) +
      scale_x_continuous(
        breaks       = seq(-0.9, 0.9, by = 0.3),
        minor_breaks = seq(-0.9, 0.9, by = 0.1),
        labels       = label_number(accuracy = 0.1),
        sec.axis     = empty_sec_axis
      ) +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
      scale_colour_manual(values = c("Univariable"          = COL_UV,
                                     "Backward Elimination" = COL_BE)) +
      labs(x = "Correlation of X1 and X2 (r)",
           y = sel_y_label(), colour = NULL) +
      app_theme() + sens_x_theme()
  })

  # Sensitivity: selection rate vs n (BE=green) ----------------------------------
  output$sens_n_plot <- renderPlot({
    req(sens_sel_df_n())
    df <- sens_sel_df_n() |>
      pivot_longer(c("Univariable", "Backward Elimination"),
                   names_to = "Method", values_to = "Rate")
    n_vals <- sort(unique(df$n))
    ggplot(df, aes(x = n, y = Rate, colour = Method, group = Method)) +
      geom_line(linewidth = 1.3) + geom_point(size = 2.8) +
      geom_vline(xintercept = input$n, linetype = "dotted",
                 colour = "#555", linewidth = 0.8) +
      annotate("text", x = input$n + 5, y = 0.04,
               label = paste0("n = ", input$n),
               size = 3.5, colour = "#555", hjust = 0) +
      scale_x_continuous(
        breaks   = pretty(n_vals, n = 8),
        labels   = label_number(accuracy = 1),
        sec.axis = sec_axis(~ . / 2, name = "n per variable (npv)",
                            breaks = pretty(n_vals / 2, n = 8),
                            labels = label_number(accuracy = 1))
      ) +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
      scale_colour_manual(values = c("Univariable"          = COL_UV,
                                     "Backward Elimination" = COL_BE)) +
      labs(x = "Sample size n",
           y = sel_y_label(), colour = NULL) +
      app_theme() + sens_x_theme() + npv_axis_theme()
  })

  # Sensitivity: mean beta-hat vs rho (UV=red, BE=green, MV=blue) ---------------
  output$sens_rho_coef_plot <- renderPlot({
    req(sens_res())
    df <- sens_res()$rho |>
      pivot_longer(cols = c(b1_uv, b2_uv, b1_mv, b2_mv, b1_be, b2_be),
                   names_to = "Series", values_to = "Mean_bhat") |>
      mutate(
        Variable = ifelse(grepl("b1", Series), "X1", "X2"),
        Method   = case_when(
          grepl("_uv", Series) ~ "Univariable",
          grepl("_mv", Series) ~ "Multivariable",
          grepl("_be", Series) ~ "Backward Elimination"
        )
      )
    refs     <- data.frame(Variable = c("X1", "X2"), true_beta = c(input$beta1, input$beta2))
    vline_df <- data.frame(xval = input$rho, lbl = paste0("Selected r = ", input$rho))
    yrange_df <- min_yrange_df(
      bind_rows(df %>% select(Variable, Mean_bhat),
                refs %>% rename(Mean_bhat = true_beta)),
      "Variable", "Mean_bhat", min_range = 0.01)
    ggplot(df, aes(x = rho, y = Mean_bhat, colour = Method, group = Method)) +
      geom_blank(data = yrange_df, aes(x = input$rho, y = ymin),
                 inherit.aes = FALSE) +
      geom_blank(data = yrange_df, aes(x = input$rho, y = ymax),
                 inherit.aes = FALSE) +
      geom_line(linewidth = 1.2) + geom_point(size = 2.4) +
      geom_hline(data = refs, aes(yintercept = true_beta, linetype = "True \u03b2"),
                 colour = "black", linewidth = 0.85) +
      geom_vline(data = vline_df, aes(xintercept = xval, linetype = lbl),
                 colour = "#555", linewidth = 0.8) +
      facet_wrap(~ Variable, ncol = 2, scales = "free_y") +
      scale_x_continuous(
        breaks       = seq(-0.9, 0.9, by = 0.3),
        minor_breaks = seq(-0.9, 0.9, by = 0.1),
        labels       = label_number(accuracy = 0.1),
        sec.axis     = empty_sec_axis
      ) +
      scale_colour_manual(name = "Model", values = COL_COEF) +
      scale_linetype_manual(
        name   = NULL,
        values = setNames(c("dashed", "dotted"), c("True \u03b2", unique(vline_df$lbl))),
        guide  = guide_legend(
          override.aes = list(colour = c("black", "#555"), linewidth = c(0.85, 0.80)))
      ) +
      labs(x = "Correlation of X1 and X2 (r)", y = "Mean \u03b2\u0302") +
      app_theme() +
      theme(strip.text = element_text(colour = "white", face = "bold", size = 13)) +
      sens_x_theme()
  })

  # Sensitivity: mean beta-hat vs n (UV=red, BE=green, MV=blue) -----------------
  output$sens_n_coef_plot <- renderPlot({
    req(sens_res())
    df <- sens_res()$n |>
      pivot_longer(cols = c(b1_uv, b2_uv, b1_mv, b2_mv, b1_be, b2_be),
                   names_to = "Series", values_to = "Mean_bhat") |>
      mutate(
        Variable = ifelse(grepl("b1", Series), "X1", "X2"),
        Method   = case_when(
          grepl("_uv", Series) ~ "Univariable",
          grepl("_mv", Series) ~ "Multivariable",
          grepl("_be", Series) ~ "Backward Elimination"
        )
      )
    refs    <- data.frame(Variable = c("X1", "X2"), true_beta = c(input$beta1, input$beta2))
    n_vals  <- sort(unique(df$n))
    vline_df <- data.frame(xval = input$n, lbl = paste0("Selected n = ", input$n))
    yrange_df <- min_yrange_df(
      bind_rows(df %>% select(Variable, Mean_bhat),
                refs %>% rename(Mean_bhat = true_beta)),
      "Variable", "Mean_bhat", min_range = 0.01)
    ggplot(df, aes(x = n, y = Mean_bhat, colour = Method, group = Method)) +
      geom_blank(data = yrange_df, aes(x = input$n, y = ymin),
                 inherit.aes = FALSE) +
      geom_blank(data = yrange_df, aes(x = input$n, y = ymax),
                 inherit.aes = FALSE) +
      geom_line(linewidth = 1.2) + geom_point(size = 2.4) +
      geom_hline(data = refs, aes(yintercept = true_beta, linetype = "True \u03b2"),
                 colour = "black", linewidth = 0.85) +
      geom_vline(data = vline_df, aes(xintercept = xval, linetype = lbl),
                 colour = "#555", linewidth = 0.8) +
      facet_wrap(~ Variable, ncol = 2, scales = "free_y") +
      scale_x_continuous(
        breaks   = pretty(n_vals, n = 8),
        labels   = label_number(accuracy = 1),
        sec.axis = sec_axis(~ . / 2, name = "n per variable (npv)",
                            breaks = pretty(n_vals / 2, n = 8),
                            labels = label_number(accuracy = 1))
      ) +
      scale_colour_manual(name = "Model", values = COL_COEF) +
      scale_linetype_manual(
        name   = NULL,
        values = setNames(c("dashed", "dotted"), c("True \u03b2", unique(vline_df$lbl))),
        guide  = guide_legend(
          override.aes = list(colour = c("black", "#555"), linewidth = c(0.85, 0.80)))
      ) +
      labs(x = "Sample size n", y = "Mean \u03b2\u0302") +
      app_theme() +
      theme(strip.text = element_text(colour = "white", face = "bold", size = 13)) +
      sens_x_theme() + npv_axis_theme()
  })

  # Tab 7: Single dataset inspector =============================================
  dataset_counter <- reactiveVal(0)
  observeEvent(input$new_dataset, { dataset_counter(dataset_counter() + 1) })
  observeEvent(input$run,         { dataset_counter(dataset_counter() + 1) })

  one_data <- eventReactive(dataset_counter(), {
    req(sigma2_val(), dataset_counter() > 0)
    set.seed(input$seed + 9999 + dataset_counter())
    Sig <- matrix(c(1, input$rho, input$rho, 1), ncol = 2)
    D   <- mvrnorm(n = input$n, mu = c(0, 0), Sigma = Sig)
    y   <- as.numeric(as.matrix(D) %*% c(input$beta1, input$beta2) +
                        rnorm(input$n, 0, sqrt(sigma2_val())))
    df  <- data.frame(x1 = D[, 1], x2 = D[, 2], y = y)

    m1 <- lm(y ~ x1,      data = df)
    m2 <- lm(y ~ x2,      data = df)
    mf <- lm(y ~ x1 + x2, data = df, x = TRUE, y = TRUE)
    mb <- tryCatch(
      abe(mf, data = df, alpha = input$alpha, tau = Inf, type.test = "F", verbose = FALSE),
      error = function(e) mf
    )
    be_out <- capture.output(
      tryCatch(
        abe(mf, data = df, alpha = input$alpha, tau = Inf, type.test = "F", verbose = TRUE),
        error = function(e) cat("BE error:", conditionMessage(e))
      )
    )
    list(df = df, m1 = m1, m2 = m2, mf = mf, mb = mb, be_out = be_out)
  })

  scatter_lm <- function(df, xcol, ycol, title, col) {
    ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
      geom_point(colour = col, alpha = 0.7, size = 2.2) +
      geom_smooth(method = "lm", se = TRUE, colour = col,
                  fill = col, alpha = 0.15, linewidth = 0.9) +
      labs(title = title, x = xcol, y = ycol) +
      app_theme() + theme(legend.position = "none")
  }

  output$single_x1y  <- renderPlot({ req(one_data()); scatter_lm(one_data()$df, "x1", "y", "X1 vs Y  (univariable)", COL_UV) })
  output$single_x2y  <- renderPlot({ req(one_data()); scatter_lm(one_data()$df, "x2", "y", "X2 vs Y  (univariable)", COL_UV) })
  output$single_x1x2 <- renderPlot({
    req(one_data())
    r_hat <- round(cor(one_data()$df$x1, one_data()$df$x2), 3)
    scatter_lm(one_data()$df, "x1", "x2", paste0("X1 vs X2  (r-hat = ", r_hat, ")"), "#27ae60")
  })
  output$single_dens_x1 <- renderPlot({
    req(one_data())
    ggplot(one_data()$df, aes(x = x1)) +
      geom_histogram(aes(y = after_stat(density)), bins = 15,
                     fill = COL_UV, alpha = 0.4, colour = "white") +
      geom_density(colour = COL_UV, linewidth = 1.0) +
      labs(title = "Distribution of X1", x = "X1", y = "Density") +
      app_theme() + theme(legend.position = "none")
  })
  output$single_dens_x2 <- renderPlot({
    req(one_data())
    ggplot(one_data()$df, aes(x = x2)) +
      geom_histogram(aes(y = after_stat(density)), bins = 15,
                     fill = COL_UV, alpha = 0.4, colour = "white") +
      geom_density(colour = COL_UV, linewidth = 1.0) +
      labs(title = "Distribution of X2", x = "X2", y = "Density") +
      app_theme() + theme(legend.position = "none")
  })

  output$single_plot_interp <- renderUI({
    req(one_data())
    r_hat <- round(cor(one_data()$df$x1, one_data()$df$x2), 3)
    div(
      style = "background:#f8f9fa; border-left:4px solid #2c3e50;
                padding:12px 18px; border-radius:4px;",
      h5("What the plots show", style = "margin-top:0;"),
      p(HTML(paste0(
        "The scatter plots above separately show the relationship of X<sub>1</sub> and X<sub>2</sub>
         with the outcome Y, and the correlation between X<sub>1</sub> and X<sub>2</sub>.
         In this dataset (n = ", input$n, "), the observed correlation between
         X<sub>1</sub> and X<sub>2</sub> is <b>", r_hat, "</b> (true r = ",
        input$rho, "). As expected from the simulation design, the density plots show that X<sub>1</sub> and X<sub>2</sub> follow an approximately normal distribution. This may not be obvious if the sample size is small."
      ))),
      uiOutput("single_clinical_para")
    )
  })

  output$single_m1 <- renderPrint({ req(one_data()); summary(one_data()$m1) })
  output$single_m2 <- renderPrint({ req(one_data()); summary(one_data()$m2) })
  output$single_mf <- renderPrint({ req(one_data()); summary(one_data()$mf) })
  output$single_be <- renderText({ req(one_data()); paste(one_data()$be_out, collapse = "\n") })

  output$single_model_interp <- renderUI({
    req(one_data())
    d   <- one_data()
    cm1 <- coef(summary(d$m1));  cm2 <- coef(summary(d$m2));  cmf <- coef(summary(d$mf))
    b1u <- round(cm1["x1", "Estimate"], 2);  p1u <- fmt_p(cm1["x1", "Pr(>|t|)"])
    b2u <- round(cm2["x2", "Estimate"], 2);  p2u <- fmt_p(cm2["x2", "Pr(>|t|)"])
    b1m <- round(cmf["x1", "Estimate"], 2);  p1m <- fmt_p(cmf["x1", "Pr(>|t|)"])
    b2m <- round(cmf["x2", "Estimate"], 2);  p2m <- fmt_p(cmf["x2", "Pr(>|t|)"])
    div(
      style = "background:#f8f9fa; border-left:4px solid #2980b9;
                padding:12px 18px; border-radius:4px; margin-top:10px;",
      h5("How to interpret the R output", style = "margin-top:0;"),
      p(HTML(paste0(
        "<b>Univariable model for X<sub>1</sub>:</b> Ignoring X<sub>2</sub>,
         a one-unit increase in X<sub>1</sub> is associated with a change of
         <b>", b1u, "</b> units in Y (p ", p1u, ").
         In the clinical example: a 1-SD increase in age is associated with a
         ", b1u, " SD change in systolic blood pressure, <i>without accounting for BMI</i>."))),
      p(HTML(paste0(
        "<b>Univariable model for X<sub>2</sub>:</b> Ignoring X<sub>1</sub>,
         a one-unit increase in X<sub>2</sub> is associated with a change of
         <b>", b2u, "</b> units in Y (p ", p2u, ").
         In the clinical example: a 1-SD increase in BMI is associated with a
         ", b2u, " SD change in systolic blood pressure, <i>without accounting for age</i>."))),
      p(HTML(paste0(
        "<b>Full multivariable model:</b> Adjusting for each other simultaneously,
         X<sub>1</sub> has a coefficient of <b>", b1m, "</b> (p ", p1m, ")
         and X<sub>2</sub> has a coefficient of <b>", b2m, "</b> (p ", p2m, ").
         These are the estimates <i>adjusted for the other predictor</i>."))),
      uiOutput("single_compare_para")
    )
  })

  output$single_be_interp <- renderUI({
    req(one_data())
    be_coefs <- names(coef(one_data()$mb))
    x1_kept  <- "x1" %in% be_coefs;  x2_kept <- "x2" %in% be_coefs
    outcome  <- case_when(
      x1_kept & x2_kept  ~ "Both X1 and X2 were retained in the final model.",
      x1_kept & !x2_kept ~ "Only X1 was retained; X2 was removed.",
      !x1_kept & x2_kept ~ "Only X2 was retained; X1 was removed.",
      TRUE               ~ "Neither X1 nor X2 was retained (only the intercept remains)."
    )
    # "correct" depends on whether each beta is zero
    should_have_x1 <- input$beta1 != 0
    should_have_x2 <- input$beta2 != 0
    correct <- (x1_kept == should_have_x1) & (x2_kept == should_have_x2)

    correct_outcome_str <- if (should_have_x1 & should_have_x2) {
      "both X<sub>1</sub> and X<sub>2</sub> retained"
    } else if (should_have_x1 & !should_have_x2) {
      "only X<sub>1</sub> retained and X<sub>2</sub> dropped"
    } else if (!should_have_x1 & should_have_x2) {
      "only X<sub>2</sub> retained and X<sub>1</sub> dropped"
    } else {
      "neither variable retained"
    }

    div(
      style = paste0("background:#f8f9fa; border-left:4px solid ",
                     if (correct) COL_BE else COL_UV,
                     "; padding:12px 18px; border-radius:4px; margin-top:10px;"),
      h5("How to interpret the R output", style = "margin-top:0;"),
      p(HTML(
        "The step-by-step output shows which predictor (if any) was removed at each step,
         based on its p-value in the current model. The algorithm starts with the full model. Elimination stops when all
         remaining predictors have p &lt; &alpha;."
      )),
      p(HTML(paste0(
        "<b>Result for this dataset:</b> ", outcome, " ",
        if (correct)
          paste0("This is the <b>correct decision</b>: the correct outcome for these parameters is ",
                 correct_outcome_str, ".")
        else
          paste0("This is an <b>incorrect decision</b>: for these parameters the correct outcome is ",
                 correct_outcome_str,
                 " (&beta;<sub>1</sub> = ", input$beta1, ", &beta;<sub>2</sub> = ", input$beta2, ").",
                 " This illustrates the kind of error that can occur in a single dataset,
                 even with a principled method.")
      )))
    )
  })

  # Summary table: UV → BE → MV, with True Estimate + Bias as last two italic columns
  output$single_results_table_ui <- renderUI({
    req(one_data())
    d   <- one_data()
    cm1 <- coef(summary(d$m1));  cm2 <- coef(summary(d$m2))
    cmf <- coef(summary(d$mf));  cmb <- coef(summary(d$mb))

    get_row <- function(cm, var) {
      if (var %in% rownames(cm))
        unname(cm[var, c("Estimate", "Std. Error", "Pr(>|t|)")])
      else
        c(NA_real_, NA_real_, NA_real_)
    }

    r1b <- get_row(cmb, "x1");  r2b <- get_row(cmb, "x2")

    rows_data <- list(
      list(model = "Univariable preselection", var = "X1",
           est = cm1["x1","Estimate"], se = cm1["x1","Std. Error"],
           pv  = cm1["x1","Pr(>|t|)"], true = input$beta1),
      list(model = "Univariable preselection", var = "X2",
           est = cm2["x2","Estimate"], se = cm2["x2","Std. Error"],
           pv  = cm2["x2","Pr(>|t|)"], true = input$beta2),
      list(model = "Backward elimination", var = "X1",
           est = r1b[1], se = r1b[2], pv = r1b[3], true = input$beta1),
      list(model = "Backward elimination", var = "X2",
           est = r2b[1], se = r2b[2], pv = r2b[3], true = input$beta2),
      list(model = "Multivariable full model", var = "X1",
           est = cmf["x1","Estimate"], se = cmf["x1","Std. Error"],
           pv  = cmf["x1","Pr(>|t|)"], true = input$beta1),
      list(model = "Multivariable full model", var = "X2",
           est = cmf["x2","Estimate"], se = cmf["x2","Std. Error"],
           pv  = cmf["x2","Pr(>|t|)"], true = input$beta2)
    )

    row_bg <- c("Univariable preselection" = "#fce8e6",
                "Backward elimination"    = "#e8f5e9",
                "Multivariable full model"  = "#e3f2fd")

    fmt_pv <- function(p) {
      if (is.na(p)) "<i>removed</i>"
      else if (p < 0.001) "&lt; 0.001"
      else sprintf("%.3f", p)
    }

    col_style <- "padding:5px 10px; border:1px solid #dee2e6;"
    italic_style <- paste0(col_style, " font-style:italic;")

    hdr <- tags$thead(tags$tr(
      lapply(c("Model","Variable","Estimate","SE","p-value",
               "<i>True Estimate</i>","<i>Bias</i>"), function(h)
                 tags$th(style = "background:#2c3e50; color:white; padding:6px 10px;",
                         HTML(h)))
    ))

    trs <- lapply(rows_data, function(r) {
      bg  <- row_bg[r$model]
      est <- if (is.na(r$est)) "\u2014" else sprintf("%.4f", r$est)
      se  <- if (is.na(r$se))  "\u2014" else sprintf("%.4f", r$se)
      bia <- if (is.na(r$est)) "\u2014" else sprintf("%.4f", r$est - r$true)
      tags$tr(
        style = paste0("background:", bg, ";"),
        tags$td(style = col_style,     r$model),
        tags$td(style = col_style,     r$var),
        tags$td(style = col_style,     est),
        tags$td(style = col_style,     se),
        tags$td(style = col_style,     HTML(fmt_pv(r$pv))),
        tags$td(style = italic_style,  sprintf("%.4f", r$true)),
        tags$td(style = italic_style,  bia)
      )
    })

    tags$table(
      class = "table table-bordered table-sm",
      style = "font-size:14px; width:auto;",
      hdr,
      tags$tbody(trs)
    )
  })


  # Tab 2: r=0 bullet ==========================================================
  output$coef_zero_r_bullet <- renderUI({
    if (!zero_r_both_active()) return(NULL)
    tags$li(HTML(
      "<b>Note (r = 0):</b> When the two predictors are uncorrelated, there is no
       omitted-variable bias: all three density curves are centred on the true &beta;
       (no shift). Hence, the densities of univariable and multivariable estimates are
       centred at the same point. However, the <b>multivariable densities are narrower</b>
       (smaller SE) because fitting both variables jointly reduces the residual variance
       and increases estimation precision."
    ))
  })

  # Tab 1: scenario-specific first three bullets ===============================
  output$sel_interp_bullet1 <- renderUI({
    if (both_zero()) {
      tags$li(HTML(
        "<b>Both selected</b> represents a <b>type I error</b> for both variables (selecting variables that have no effect). A lower proportion is better here."
      ))
    } else if (null_beta()) {
      if (b1_is_zero()) {
        tags$li(HTML(
          "<b>'Only X<sub>2</sub> selected'</b> is the correct decision. A higher proportion for this outcome is better."
        ))
      } else {
        tags$li(HTML(
          "<b>'Only X<sub>1</sub> selected'</b> is the correct decision. A higher proportion for this outcome is better."
        ))
      }
    } else {
      tags$li(HTML(
        "<b>Both selected</b> corresponds to the correct decision. A higher proportion in the simulation is better."
      ))
    }
  })
  output$sel_interp_bullet2 <- renderUI({
    if (both_zero()) {
      tags$li(HTML(
        "<b>Only X<sub>1</sub> / Only X<sub>2</sub></b> each represent a single false positive (selecting one variable that has no effect)."
      ))
    } else if (null_beta()) {
      if (b1_is_zero()) {
        tags$li(HTML(
          "<b>'Both selected'</b> is a false positive for X<sub>1</sub> (type I error). <b>'Only X<sub>1</sub>'</b> is also an error: X<sub>1</sub> is selected but X<sub>2</sub> is missed."
        ))
      } else {
        tags$li(HTML(
          "<b>'Both selected'</b> is a false positive for X<sub>2</sub> (type I error). <b>'Only X<sub>2</sub>'</b> is also an error: X<sub>2</sub> is selected but X<sub>1</sub> is missed."
        ))
      }
    } else {
      tags$li(HTML(
        "<b>Only X<sub>1</sub> / Only X<sub>2</sub></b> means that one true predictor was missed (false negative)."
      ))
    }
  })
  output$sel_interp_bullet3 <- renderUI({
    if (both_zero()) {
      tags$li(HTML(
        "<b>Neither selected</b> is the <b>correct decision</b>. A higher proportion for this outcome is better."
      ))
    } else if (null_beta()) {
      if (b1_is_zero()) {
        tags$li(HTML(
          "<b>'Neither selected'</b> is an error: the true predictor X<sub>2</sub> was missed."
        ))
      } else {
        tags$li(HTML(
          "<b>'Neither selected'</b> is an error: the true predictor X<sub>1</sub> was missed."
        ))
      }
    } else {
      tags$li(HTML(
        "<b>Neither</b> corresponds to both predictors dropped by the variable selection method. This is the worst outcome."
      ))
    }
  })

  # Tab 3: BE excluded count bullet ============================================
  output$pval_be_bullet <- renderUI({
    req(sim_res())
    res <- sim_res()
    n_excl_x1  <- sum(is.na(res$p1_be))
    n_excl_x2  <- sum(is.na(res$p2_be))
    pct_x1     <- round(100 * n_excl_x1 / input$nsim, 1)
    pct_x2     <- round(100 * n_excl_x2 / input$nsim, 1)
    tags$li(HTML(paste0(
      "<span style='color:#27ae60;'><b>Backward Elimination (green)</b></span>: ",
      "p-values are shown only for simulations where the variable was retained in the final model. ",
      "They are always below &alpha; by construction (that is the stopping criterion of BE).<br>In this simulation:<br>",
      "For X<sub>1</sub>: <b>", n_excl_x1, " simulations (", pct_x1, "%)</b> were excluded ",
      "because X<sub>1</sub> was dropped by BE. ",
      "<br>For X<sub>2</sub>: <b>", n_excl_x2, " simulations (", pct_x2, "%)</b> were excluded ",
      "because X<sub>2</sub> was dropped by BE."
    )))
  })

  # Tab 3: r=0 bullet (only shown when r=0 and both betas active) ==============
  output$pval_zero_r_bullet <- renderUI({
    if (!zero_r_both_active()) return(NULL)
    tags$li(HTML(
      "<b>Note (r = 0):</b> When the two predictors are uncorrelated, there is no
       omitted-variable bias and the univariable estimates are unbiased.
       However, the <b>residual variance in the multivariable model is smaller</b> than
       in the univariable models, because both X<sub>1</sub> and X<sub>2</sub> explain
       independent portions of the variance in Y.
       As a result, the standard errors (and therefore p-values) for both variables are
       <i>smaller</i> in the multivariable model than in the univariable models.
       This effect is particularly pronounced for the 'weaker' variable (the one with
       the smaller |&beta;|): its univariable p-value may exceed &alpha; (not significant),
       while its multivariable p-value falls below &alpha; (significant),
       because the 'stronger' variable reduces the residual variance and thereby
       increases the precision of all coefficient estimates.
       This is an argument <i>in favour</i> of fitting the multivariable model even when r = 0."
    ))
  })

  # Tab 6 row 1: scenario-aware Y-axis description =============================
  output$sens_yaxis_bullet <- renderUI({
    lbl <- sel_y_label()
    if (both_zero()) {
      tags$li(HTML(paste0(
        "<b>Y axis of both plots:</b> ", lbl, ".
         In this null model scenario, selecting both variables is an <b>error</b> (joint type&nbsp;I error).
         A <i>lower</i> proportion is better here."
      )))
    } else if (b1_is_zero() & !b2_is_zero()) {
      tags$li(HTML(paste0(
        "<b>Y axis of both plots:</b> ", lbl, ".
         X<sub>1</sub> has no true effect (&beta;<sub>1</sub> = 0) and X<sub>2</sub> does.
         The correct decision is to select X<sub>2</sub> and <i>not</i> select X<sub>1</sub>.
         A <i>higher</i> proportion represents better performance."
      )))
    } else if (b2_is_zero() & !b1_is_zero()) {
      tags$li(HTML(paste0(
        "<b>Y axis of both plots:</b> ", lbl, ".
         X<sub>2</sub> has no true effect (&beta;<sub>2</sub> = 0) and X<sub>1</sub> does.
         The correct decision is to select X<sub>1</sub> and <i>not</i> select X<sub>2</sub>.
         A <i>higher</i> proportion represents better performance."
      )))
    } else {
      tags$li(HTML(paste0(
        "<b>Y axis of both plots:</b> ", lbl, " is
         the proportion of simulations in which the variable selection method correctly
         keeps both X<sub>1</sub> and X<sub>2</sub>. A <i>higher</i> proportion is better."
      )))
    }
  })

  # Tab 7: dynamic clinical illustration paragraph =============================
  output$single_clinical_para <- renderUI({
    r_val <- input$rho
    r_desc <- if (abs(r_val) < 1e-9) {
      "No correlation (r = 0) between age and BMI means the two variables are independent in this cohort."
    } else if (r_val < 0) {
      paste0("A negative correlation between age and BMI (r = ", r_val,
             ") reflects a hypothetical cohort where older participants tend to have slightly lower BMI (e.g. a geriatric cohort).")
    } else {
      paste0("A positive correlation between age and BMI (r = ", r_val,
             ") reflects a hypothetical cohort where older participants tend to have slightly higher BMI.")
    }
    misleading_sent <- if (ovb_active()) {
      " In this setting, testing age and BMI separately (univariable) can be misleading because each variable partially absorbs or masks the effect of the other."
    } else if (zero_r_both_active()) {
      " Because r = 0 here, each univariable test gives an unbiased estimate of its own coefficient. However, fitting both predictors jointly (multivariable) reduces the residual variance and thereby increases power for each variable."
    } else {
      ""
    }
    p(HTML(paste0(
      "<b>Clinical illustration:</b>
       Imagine X<sub>1</sub> represents <i>standardised age</i>
       (1 unit &asymp; 1 SD above the population mean age, e.g. ~10 years above the mean)
       and X<sub>2</sub> represents <i>standardised BMI</i>
       (1 unit &asymp; 1 SD above mean BMI, e.g. ~4 kg/m<sup>2</sup> above the mean).
       The outcome Y could be a standardised measure of <i>systolic blood pressure</i>.
       ", r_desc, misleading_sent
    )))
  })

  # Tab 7: scenario-specific compare sentence ==================================
  output$single_compare_para <- renderUI({
    req(one_data())
    true_str <- paste0("True &beta;<sub>1</sub> = ", input$beta1,
                       " and true &beta;<sub>2</sub> = ", input$beta2, ". ")
    compare_str <- if (ovb_active()) {
      "<i>Compare the univariable and multivariable estimates: if they differ noticeably, this reflects the omitted-variable bias caused by the correlation between X<sub>1</sub> and X<sub>2</sub>.</i>"
    } else if (zero_r_both_active()) {
      "<i>Because r = 0, the univariable and multivariable coefficient estimates should be similar on average (no OVB). However, the standard errors and p-values will typically be smaller in the multivariable model, because adding the second predictor reduces the residual variance.</i>"
    } else if (null_beta()) {
      if (b1_is_zero()) {
        "<i>X<sub>1</sub> has no true effect (&beta;<sub>1</sub> = 0). The multivariable estimate for X<sub>1</sub> should be close to 0; any deviation is sampling variability. If r &ne; 0, the univariable estimate for X<sub>1</sub> will be biased toward &beta;<sub>2</sub> &times; r.</i>"
      } else if (b2_is_zero()) {
        "<i>X<sub>2</sub> has no true effect (&beta;<sub>2</sub> = 0). The multivariable estimate for X<sub>2</sub> should be close to 0; any deviation is sampling variability. If r &ne; 0, the univariable estimate for X<sub>2</sub> will be biased toward &beta;<sub>1</sub> &times; r.</i>"
      } else {
        "<i>Neither variable has a true effect. All estimates should be close to 0; deviations are sampling variability.</i>"
      }
    } else {
      "<i>r = 0 and both coefficients are nonzero. The univariable and multivariable estimates should agree on average; any difference is sampling variability.</i>"
    }
    p(HTML(paste0(true_str, compare_str)))
  })


  # both_zero banner helper ====================================================
  both_zero_banner <- function() {
    div(
      style = "background:#fce4ec; border-left:4px solid #c0392b;
                padding:12px 18px; border-radius:4px; margin-bottom:16px;",
      tags$b("⚠ Both β₁ = 0 and β₂ = 0"),
      tags$p(
        style = "margin:6px 0 0; font-size:13px;",
        "Coefficient density plot is not meaningful when both β₁ = 0 and β₂ = 0 (Y has no systematic variation). Please set at least one coefficient to a non-zero value."
      )
    )
  }
  output$both_zero_warn_sel     <- renderUI({ if (both_zero()) both_zero_banner() })
  output$both_zero_warn_coef    <- renderUI({ if (both_zero()) both_zero_banner() })
  output$both_zero_warn_pval    <- renderUI({ if (both_zero()) both_zero_banner() })
  output$both_zero_warn_bv      <- renderUI({ if (both_zero()) both_zero_banner() })
  output$both_zero_warn_scatter <- renderUI({ if (both_zero()) both_zero_banner() })
  output$both_zero_warn_sens    <- renderUI({ if (both_zero()) both_zero_banner() })



  # Tab 1: scenario-specific R2 bullet =========================================
  output$sel_r2_bullet <- renderUI({
    if (both_zero()) {
      tags$li(HTML(paste0(
        "R<sup>2</sup> = ", input$r2, ". When both coefficients are 0, R<sup>2</sup>
         controls the noise level, but it does not affect the type&nbsp;I error rate
         (which should remain at &asymp; &alpha; regardless of R<sup>2</sup>)."
      )))
    } else if (null_beta()) {
      if (b1_is_zero()) {
        tags$li(HTML(paste0(
          "R<sup>2</sup> = ", input$r2, ". A higher R<sup>2</sup> increases the power to detect X<sub>2</sub> (the true predictor) for both methods. It does not reduce the type&nbsp;I error rate for X<sub>1</sub> (the null variable)."
        )))
      } else {
        tags$li(HTML(paste0(
          "R<sup>2</sup> = ", input$r2, ". A higher R<sup>2</sup> increases the power to detect X<sub>1</sub> (the true predictor) for both methods. It does not reduce the type&nbsp;I error rate for X<sub>2</sub> (the null variable)."
        )))
      }
    } else {
      tags$li(HTML(paste0(
        "R<sup>2</sup> = ", input$r2, ". A higher R<sup>2</sup> increases statistical power, so both methods are more likely to select both X<sub>1</sub> and X<sub>2</sub> correctly. A lower R<sup>2</sup> makes correct selection harder for both methods, but the advantage of backward elimination over univariable preselection is typically maintained."
      )))
    }
  })

  # Tab 2: adaptive coefficient bullets ========================================
  output$coef_uv_bullet <- renderUI({
    if (null_beta()) {
      if (b1_is_zero() & !b2_is_zero()) {
        tags$li(HTML("<span style='color:#e74c3c;'><b>Univariable (red)</b></span>: For X<sub>1</sub> (&beta;<sub>1</sub> = 0), a non-zero r causes a spurious shift in the univariable estimate (E[&beta;&#770;<sub>1,uv</sub>] = &beta;<sub>2</sub>r). For X<sub>2</sub>, the univariable estimate is unbiased."))
      } else if (!b1_is_zero() & b2_is_zero()) {
        tags$li(HTML("<span style='color:#e74c3c;'><b>Univariable (red)</b></span>: For X<sub>2</sub> (&beta;<sub>2</sub> = 0), a non-zero r causes a spurious shift in the univariable estimate (E[&beta;&#770;<sub>2,uv</sub>] = &beta;<sub>1</sub>r). For X<sub>1</sub>, the univariable estimate is unbiased."))
      } else {
        tags$li(HTML("<span style='color:#e74c3c;'><b>Univariable (red)</b></span>: Both coefficients are 0; all density curves are centred near 0 (no systematic variation)."))
      }
    } else {
      tags$li(HTML("<span style='color:#e74c3c;'><b>Univariable (red)</b></span>: A shifted density indicates omitted-variable bias. The estimated coefficient is systematically wrong."))
    }
  })
  output$coef_mv_bullet <- renderUI({
    if (null_beta()) {
      tags$li(HTML("<span style='color:#2980b9;'><b>Multivariable (blue)</b></span>: Estimates centred on the true value for both variables (0 for the null variable, true &beta; for the active variable). The multivariable model is unbiased regardless of r."))
    } else {
      tags$li(HTML("<span style='color:#2980b9;'><b>Multivariable (blue)</b></span>: Estimates centred on the true value. Hence, the estimates are unbiased."))
    }
  })
  output$coef_be_bullet <- renderUI({
    if (null_beta()) {
      if (b1_is_zero() & !b2_is_zero()) {
        tags$li(HTML("<span style='color:#27ae60;'><b>Backward Elimination (green)</b></span>: BE tends to drop X<sub>1</sub> (the null variable) from the model. When it does, the estimate is set to 0. The green density for X<sub>1</sub> is a mixture at 0 (dropped) and the multivariable estimate (retained). BE correctly concentrates on X<sub>2</sub> as n grows."))
      } else if (!b1_is_zero() & b2_is_zero()) {
        tags$li(HTML("<span style='color:#27ae60;'><b>Backward Elimination (green)</b></span>: BE tends to drop X<sub>2</sub> (the null variable) from the model. When it does, the estimate is set to 0. The green density for X<sub>2</sub> is a mixture at 0 (dropped) and the multivariable estimate (retained). BE correctly concentrates on X<sub>1</sub> as n grows."))
      } else {
        tags$li(HTML("<span style='color:#27ae60;'><b>Backward Elimination (green)</b></span>: BE should drop both variables (null model). Its estimates are a mixture at 0 (dropped, correct) and the multivariable estimates (when it fails to drop, type I error)."))
      }
    } else {
      tags$li(HTML("<span style='color:#27ae60;'><b>Backward Elimination (green)</b></span>: When it retains both variables, its estimates match the multivariable model. When it drops a variable (set to 0), the unconditional mean is biased toward zero."))
    }
  })

  # Tab 7: Combined unadjusted & adjusted prediction plots (one legend each) ===
  # Replaces the previously separate "unadjusted" and "adjusted" plots: both
  # the unadjusted (univariable) and adjusted (multivariable, other predictor
  # fixed at its mean) regression lines are now drawn in the same panel, with
  # one shared legend distinguishing the two lines.
  single_pred_plot <- function(var_x, var_other) {
    req(one_data())
    d   <- one_data()
    df  <- d$df
    mf  <- d$mf
    cf  <- coef(mf)
    x   <- df[[var_x]]
    x_other_mean <- mean(df[[var_other]])
    x_seq   <- seq(min(x), max(x), length.out = 200)
    if (var_x == "x1") {
      y_pred    <- cf["(Intercept)"] + cf["x1"] * x_seq + cf["x2"] * x_other_mean
      b_mv      <- cf["x1"]
      b_uv      <- coef(d$m1)["x1"]
      title_str <- paste0("X1 vs Y (unadjusted & adjusted, X2 = ", round(x_other_mean, 2), ")")
    } else {
      y_pred    <- cf["(Intercept)"] + cf["x1"] * x_other_mean + cf["x2"] * x_seq
      b_mv      <- cf["x2"]
      b_uv      <- coef(d$m2)["x2"]
      title_str <- paste0("X2 vs Y (unadjusted & adjusted, X1 = ", round(x_other_mean, 2), ")")
    }
    pred_df <- data.frame(x = x_seq, y = y_pred)
    lbl_uv  <- "Unadjusted (univariable)"
    lbl_mv  <- "Adjusted (multivariable, other = mean)"
    ggplot(df, aes(x = .data[[var_x]], y = y)) +
      geom_point(colour = "grey40", alpha = 0.6, size = 2) +
      geom_smooth(aes(colour = lbl_uv, fill = lbl_uv),
                  method = "lm", se = TRUE, linewidth = 1.1, alpha = 0.15,
                  formula = y ~ x) +
      geom_line(data = pred_df, aes(x = x, y = y, colour = lbl_mv),
                linewidth = 1.1, inherit.aes = FALSE) +
      scale_colour_manual(name = NULL,
                          values = setNames(c(COL_UV, COL_MV), c(lbl_uv, lbl_mv))) +
      scale_fill_manual(name = NULL, values = setNames(COL_UV, lbl_uv), guide = "none") +
      labs(title = title_str, x = var_x, y = "y",
           subtitle = paste0("Unadjusted slope = ", round(b_uv, 3),
                             "   |   Adjusted slope = ", round(b_mv, 3))) +
      app_theme()
  }
  output$single_pred_x1 <- renderPlot({ single_pred_plot("x1", "x2") })
  output$single_pred_x2 <- renderPlot({ single_pred_plot("x2", "x1") })

  # Tab 7: Explanation box for the combined prediction plots ===================
  output$single_pred_interp <- renderUI({
    req(one_data())
    d  <- one_data()
    cf <- coef(d$mf)
    div(
      style = "background:#f8f9fa; border-left:4px solid #16a085;
                padding:12px 18px; border-radius:4px; margin-top:10px;",
      h5("How to read the prediction plots", style = "margin-top:0;"),
      p(HTML(paste0(
        "Each panel overlays two regression lines for the same predictor: the ",
        "<b><span style='color:", COL_UV, ";'>unadjusted (univariable)</span></b> line, fitted while ",
        "ignoring the other predictor, and the ",
        "<b><span style='color:", COL_MV, ";'>adjusted (multivariable)</span></b> line, which holds the ",
        "other predictor fixed at its mean. If the correlation between X<sub>1</sub> and X<sub>2</sub> (r) ",
        "is non-zero and at least one true coefficient is non-zero, the two lines will diverge in slope ",
        "&mdash; this divergence <i>is</i> the omitted-variable bias. If the lines nearly coincide, ",
        "adjustment makes little difference for this predictor."
      ))),
      p(HTML(paste0(
        "<b>In this dataset:</b> unadjusted slope for X<sub>1</sub> = ", round(coef(d$m1)["x1"], 3),
        ", adjusted slope for X<sub>1</sub> = ", round(cf["x1"], 3),
        "; unadjusted slope for X<sub>2</sub> = ", round(coef(d$m2)["x2"], 3),
        ", adjusted slope for X<sub>2</sub> = ", round(cf["x2"], 3), "."
      )))
    )
  })

  # -----------------------------------------------------------------------------
  # The following two plot blocks (Partial Residual Plots and Added-Variable
  # Plots) are intentionally disabled (commented out) per request. They are
  # kept here, fully functional, in case they should be re-enabled later.
  # The matching plotOutput() calls in the UI are commented out as well.
  # -----------------------------------------------------------------------------
  # # Tab 7: Partial residual (component + residual) plots =======================
  # single_partres_plot <- function(var_x) {
  #   req(one_data())
  #   d   <- one_data()
  #   df  <- d$df
  #   mf  <- d$mf
  #   cf  <- coef(mf)
  #   x   <- df[[var_x]]
  #   b   <- if (var_x == "x1") cf["x1"] else cf["x2"]
  #   partial_resid <- b * x + residuals(mf)
  #   prd <- data.frame(x = x, pr = partial_resid)
  #   title_str <- if (var_x == "x1") "Partial Residual: X1" else "Partial Residual: X2"
  #   ggplot(prd, aes(x = x, y = pr)) +
  #     geom_point(colour = COL_MV, alpha = 0.6, size = 2) +
  #     geom_smooth(method = "lm", se = TRUE, colour = COL_MV,
  #                 fill = COL_MV, alpha = 0.15, linewidth = 1.0,
  #                 formula = y ~ x) +
  #     labs(title = title_str, x = var_x,
  #          y = paste0("β&#x302;*", var_x, " + residual"),
  #          subtitle = paste0("Slope = ", round(b, 3),
  #                            " = multivariable coefficient")) +
  #     app_theme() + theme(legend.position = "none")
  # }
  # output$single_partres_x1 <- renderPlot({ single_partres_plot("x1") })
  # output$single_partres_x2 <- renderPlot({ single_partres_plot("x2") })
  #
  # # Tab 7: Added-variable (partial regression) plots ===========================
  # single_avplot <- function(var_x, var_other) {
  #   req(one_data())
  #   d   <- one_data()
  #   df  <- d$df
  #   mf  <- d$mf
  #   cf  <- coef(mf)
  #   # Residuals of Y on other predictor
  #   e_y  <- residuals(lm(as.formula(paste("y ~", var_other)), data = df))
  #   # Residuals of X_j on other predictor
  #   e_xj <- residuals(lm(as.formula(paste(var_x, "~", var_other)), data = df))
  #   b    <- if (var_x == "x1") cf["x1"] else cf["x2"]
  #   avd  <- data.frame(e_xj = e_xj, e_y = e_y)
  #   title_str <- if (var_x == "x1") "Added-Variable: X1 | X2" else "Added-Variable: X2 | X1"
  #   ggplot(avd, aes(x = e_xj, y = e_y)) +
  #     geom_point(colour = COL_BE, alpha = 0.6, size = 2) +
  #     geom_smooth(method = "lm", se = TRUE, colour = COL_BE,
  #                 fill = COL_BE, alpha = 0.15, linewidth = 1.0,
  #                 formula = y ~ x) +
  #     labs(title = title_str,
  #          x = paste0("e(", var_x, " | ", var_other, ")"),
  #          y = paste0("e(Y | ", var_other, ")"),
  #          subtitle = paste0("Slope = ", round(b, 3),
  #                            " = multivariable coefficient of ", var_x)) +
  #     app_theme() + theme(legend.position = "none")
  # }
  # output$single_avplot_x1 <- renderPlot({ single_avplot("x1", "x2") })
  # output$single_avplot_x2 <- renderPlot({ single_avplot("x2", "x1") })

}

shinyApp(ui = ui, server = server)
