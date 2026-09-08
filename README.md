# Matrix Spectral Decomposition for Phenotypic Data

Estimate the effective number of independent phenotypic variables from their correlation structure and calculate corresponding multiple-testing significance thresholds in R.

The `run_matrix_decomposition()` function reads participant-level measurements from a CSV file, excludes the participant ID column, and calculates estimates based on [Nyholt (2004)](https://doi.org/10.1086/383251) and [Li and Ji (2005)](https://doi.org/10.1038/sj.hdy.6800717). It saves a correlation matrix and a text report, and returns the results for further analysis.

This README describes the function implementation, including its use of absolute correlations and its handling of negative eigenvalues.

## Requirements

- R
- The `data.table` package

Install the dependency from an R session:

```r
install.packages("data.table")
```

All other functions used by the script are included with R. RStudio is optional.

## Input format

Provide a comma-separated CSV file with a header row and one participant per row:

```csv
participant_id,ROI_1,ROI_2,ROI_3
sub_001,2.45,3.12,1.87
sub_002,2.61,3.05,1.92
sub_003,2.38,3.21,NA
sub_004,2.72,3.18,2.04
```

- The first column contains participant IDs and is excluded by position, regardless of its name or data type.
- Every remaining column is treated as a phenotype and must contain numeric measurements or missing values. Include at least one phenotype.
- Each phenotype must have nonzero variance, with enough overlapping observations for every pairwise correlation to be defined.
- Missing measurements are handled using `use = "pairwise.complete.obs"`: each correlation uses participants with observations for both variables.

Prepare any transformations, covariate adjustment, or variable filtering before running the function. Additional ID or metadata columns are not automatically excluded.

## Quick start

Save the function definition in `matrix_decomposition.R`. The examples below assume this script is in your working directory and your input file is `data/phenotypes.csv`.

```r
source("matrix_decomposition.R")

dir.create("results", recursive = TRUE, showWarnings = FALSE)

results <- run_matrix_decomposition(
  input_file = "data/phenotypes.csv",
  correlation_file = "results/phenotype_correlations.rds",
  report_file = "results/matSpDlite.out",
  alpha = 0.05
)

# Inspect the effective variable counts and significance thresholds.
results$effective_variables
results$thresholds
```

Sourcing the function definition alone does not run the analysis. Call `run_matrix_decomposition()` to generate the outputs.

## Function arguments

```r
run_matrix_decomposition(
  input_file,
  correlation_file,
  report_file = "matSpDlite.out",
  alpha = 0.05
)
```

| Argument | Default | Description |
| --- | --- | --- |
| `input_file` | Required | Path to the input CSV, with participant IDs in its first column. |
| `correlation_file` | Required | Destination for the correlation matrix saved as an RDS file. |
| `report_file` | `"matSpDlite.out"` | Destination for the plain-text analysis report. |
| `alpha` | `0.05` | Target experiment-wide significance level; supply a single number strictly between 0 and 1. |

Paths can be absolute or relative to the R working directory. Output directories must already exist; the quick-start example creates its output directory explicitly. Existing files at the specified output paths are overwritten. The function does not explicitly validate `alpha`.

## Outputs

### Saved files

| Output | Contents |
| --- | --- |
| Correlation RDS | The signed Pearson correlation matrix, stored as an R data frame with phenotype names. Participant IDs are excluded. |
| Text report | Number of phenotypes, adjusted eigenvalues, adjusted eigenvalue variance, Nyholt and Li–Ji effective counts, a Bonferroni-style Nyholt threshold, and a Šidák-style Li–Ji threshold. |

The report displays eigenvalues, variance, and effective counts to four decimal places. The returned results retain numerical precision.

Read the saved correlation data frame with:

```r
correlation_df <- readRDS("results/phenotype_correlations.rds")
correlation_matrix <- as.matrix(correlation_df)
```

### Returned R object

The function invisibly returns a list. Assign the call to an object, as in the quick-start example, to retain it.

| Element | Contents |
| --- | --- |
| `correlation_matrix` | Signed Pearson correlation matrix. |
| `eigenvalues` | Eigenvalues of the absolute correlation matrix, before clipping. |
| `adjusted_eigenvalues` | Eigenvalues after negative values are set to zero. |
| `n_variables` | Number of phenotype columns, excluding the first column. |
| `effective_variables` | Named vector containing `nyholt_original`, `nyholt_adjusted`, and `li_ji`. |
| `thresholds` | Named vector containing `bonferroni_nyholt`, `sidak_nyholt`, and `sidak_li_ji`. |

For example:

```r
results$n_variables
results$effective_variables[["li_ji"]]
results$thresholds[["sidak_li_ji"]]

# Inspect the spectrum before clipping negative eigenvalues.
range(results$eigenvalues)
sum(results$eigenvalues < 0)
```

The original Nyholt estimate and the Šidák-style Nyholt threshold are available in the returned object but are not printed in the text report.

## Statistical calculations

### Correlation matrix and eigenvalues

Let $R$ be the signed Pearson correlation matrix and $M$ the number of phenotype columns. The function decomposes the elementwise absolute matrix:

$$
A = \lvert R \rvert.
$$

For eigenvalues $\lambda_1, \ldots, \lambda_M$ of $A$, the adjusted values are:

$$
\lambda_i^{+} = \max(\lambda_i, 0).
$$

The saved correlation matrix retains its original signs. Taking absolute values and clipping negative eigenvalues affect the decomposition calculations only.

### Nyholt estimate

Following the variance-based approach described by [Nyholt (2004)](https://doi.org/10.1086/383251), the implementation calculates:

$$
M_{\mathrm{eff,Nyholt}}
= M\left(1 - \frac{(M-1)s_{\lambda}^{2}}{M^{2}}\right),
$$

where $s_{\lambda}^{2}$ is the sample variance calculated by R's `var()` function.

The function evaluates this expression twice:

- `nyholt_original` uses the eigenvalues of $A$ before clipping.
- `nyholt_adjusted` uses the adjusted eigenvalues $\lambda_i^{+}$.

The Nyholt thresholds use `nyholt_adjusted`. If all original eigenvalues equal 1, the function sets both eigenvalue variances to zero, yielding an effective count of $M$.

### Li–Ji estimate

The implementation applies the Equation 5 expression attributed to [Li and Ji (2005)](https://doi.org/10.1038/sj.hdy.6800717) to the adjusted eigenvalues:

$$
M_{\mathrm{eff,LiJi}} = \sum_{i=1}^{M} \left[\mathbf{1}(\lambda_i^{+} \geq 1) + \lambda_i^{+} - \left\lfloor \lambda_i^{+} \right\rfloor\right].
$$

Here, $\mathbf{1}(\lambda_i^{+} \geq 1)$ equals 1 when the adjusted eigenvalue is at least 1 and 0 otherwise. Effective counts can be non-integer values.

### Significance thresholds

Let $N$ be `nyholt_adjusted` and $L$ be `li_ji`. The returned thresholds are:

| Returned name | Formula | Type |
| --- | --- | --- |
| `bonferroni_nyholt` | $\alpha / N$ | Bonferroni-style threshold using the adjusted Nyholt count. |
| `sidak_nyholt` | $1 - (1-\alpha)^{1/N}$ | Šidák-style threshold using the adjusted Nyholt count. |
| `sidak_li_ji` | $1 - (1-\alpha)^{1/L}$ | Šidák-style threshold using the Li–Ji count. |

These formulas substitute an estimated effective count for the total number of tests. They provide approximate multiplicity adjustments; an estimated count alone does not establish family-wise error control for every study design.

Apply a selected cutoff to p-values calculated separately for the corresponding family of phenotype tests. The function returns all three thresholds without automatically selecting one. Its report retains the supplied matSpDlite guidance to use Li–Ji unless the adjusted Nyholt effective count is smaller.

## Implementation notes

- **Missing and constant data:** Nonnumeric phenotype columns, constant columns, or insufficient paired observations can prevent correlation calculation or leave undefined entries that cause eigendecomposition to fail. The function does not automatically clean these inputs.
- **Negative eigenvalues:** Pairwise missing-data handling can produce a correlation matrix that is not positive semidefinite, as described in the [R correlation documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html). The function clips negative eigenvalues of the absolute matrix to zero, without reconstructing a valid correlation matrix or renormalizing the spectrum. Review the original eigenvalues when interpreting the results; the report does not emit a negative-eigenvalue warning.
- **Redundant variables:** Duplicate and collinear phenotype columns are retained. `n_variables` is the input phenotype count, rather than a count after redundancy removal.
- **Numerical boundaries:** The Li–Ji calculation uses direct floating-point comparisons and `floor()` without a tolerance near integer eigenvalues.
- **Scope of correction:** The estimates describe dependence among the supplied phenotypes. Account separately for additional testing dimensions, such as multiple predictors or analysis models, when defining the family of tests.

## References and attribution

If using these methods in research, cite the original methodological papers:

1. Nyholt, D. R. (2004). [A simple correction for multiple testing for single-nucleotide polymorphisms in linkage disequilibrium with each other](https://doi.org/10.1086/383251). *The American Journal of Human Genetics*, **74**(4), 765–769.
2. Li, J., & Ji, L. (2005). [Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix](https://doi.org/10.1038/sj.hdy.6800717). *Heredity*, **95**, 221–227.

The supplied source credits Dale Nyholt's [matSpDlite implementation](http://gump.qimr.edu.au/general/daleN/matSpDlite/).
