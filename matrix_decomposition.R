# Matrix spectral decomposition for phenotypic data
# References:
# Nyholt DR (2004). Am J Hum Genet 74(4):765-769.
# Li and Ji (2005). Heredity 95:221-227.
# Original implementation: http://gump.qimr.edu.au/general/daleN/matSpDlite/

#' Estimate the effective number of independent phenotypic variables
#'
#' @param input_file CSV file with participant IDs in the first column.
#' @param correlation_file Destination for the correlation matrix RDS file.
#' @param report_file Destination for the text report.
#' @param alpha Experiment-wide significance level.
#' @return Invisibly returns the correlation matrix, eigenvalues, effective
#'   variable counts, and significance thresholds.
run_matrix_decomposition <- function(
    input_file,
    correlation_file,
    report_file = "matSpDlite.out",
    alpha = 0.05
) {
  # Read phenotypes and exclude the participant ID column.
  phenotypes <- data.table::fread(input_file, header = TRUE, sep = ",")
  phenotype_data <- phenotypes[, -1, with = FALSE]

  correlation_matrix <- cor(
    phenotype_data,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  saveRDS(as.data.frame(correlation_matrix), file = correlation_file)

  # Decompose the absolute correlation matrix.
  eigenvalues <- eigen(
    t(abs(correlation_matrix)),
    symmetric = TRUE,
    only.values = TRUE
  )$values
  adjusted_eigenvalues <- pmax(eigenvalues, 0)
  n_variables <- length(eigenvalues)

  original_variance <- var(eigenvalues)
  adjusted_variance <- var(adjusted_eigenvalues)

  if (all(eigenvalues == 1)) {
    original_variance <- 0
    adjusted_variance <- 0
  }

  # Nyholt estimates, before and after setting negative eigenvalues to zero.
  effective_original <- n_variables * (
    1 - (n_variables - 1) * original_variance / n_variables^2
  )
  effective_nyholt <- n_variables * (
    1 - (n_variables - 1) * adjusted_variance / n_variables^2
  )

  # Li and Ji estimate (Equation 5).
  effective_li_ji <- sum(
    (adjusted_eigenvalues >= 1) +
      (adjusted_eigenvalues - floor(adjusted_eigenvalues))
  )

  thresholds <- c(
    bonferroni_nyholt = alpha / effective_nyholt,
    sidak_nyholt = 1 - (1 - alpha)^(1 / effective_nyholt),
    sidak_li_ji = 1 - (1 - alpha)^(1 / effective_li_ji)
  )

  # Write the report.
  eigenvalue_table <- data.frame(
    Factor = seq_len(n_variables),
    Eigenvalue = round(adjusted_eigenvalues, 4)
  )

  report <- c(
    "Matrix spectral decomposition: Phenotypic data",
    sprintf("Number of variables: %d", n_variables),
    "",
    "Eigenvalues (negative values set to zero):",
    capture.output(print(eigenvalue_table, row.names = FALSE)),
    "",
    "Nyholt (2004; negative eigenvalues set to zero)",
    sprintf("Variance of eigenvalues: %.4f", adjusted_variance),
    sprintf("Effective number of independent variables: %.4f", effective_nyholt),
    sprintf(
      "Bonferroni significance threshold (alpha / Veff; alpha = %g): %.8g",
      alpha, thresholds[["bonferroni_nyholt"]]
    ),
    "",
    "Li and Ji (2005; Equation 5)",
    sprintf("Effective number of independent variables: %.4f", effective_li_ji),
    sprintf(
      "Experiment-wide significance threshold (alpha = %g): %.8g",
      alpha, thresholds[["sidak_li_ji"]]
    ),
    "",
    "Original matSpDlite guidance: use Li and Ji unless Veff < VeffLi."
  )
  writeLines(report, con = report_file)

  invisible(list(
    correlation_matrix = correlation_matrix,
    eigenvalues = eigenvalues,
    adjusted_eigenvalues = adjusted_eigenvalues,
    n_variables = n_variables,
    effective_variables = c(
      nyholt_original = effective_original,
      nyholt_adjusted = effective_nyholt,
      li_ji = effective_li_ji
    ),
    thresholds = thresholds
  ))
}

# Example usage:
# results <- run_matrix_decomposition(
#   input_file = "/Users/amir/Desktop/test/All_regional_HCP_ready.csv",
#   correlation_file = "/Users/amir/Desktop/test/correlationbetweenROIS.rds",
#   report_file = "/Users/amir/Desktop/test/matSpDlite.out"
# )
