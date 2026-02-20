# ============================================================
# R/qtl: GENERIC QTL ANALYSIS PIPELINE
# Works for any backcross (bc) or intercross (f2) dataset
# ============================================================

library(qtl)


# ============================================================
# USER SETTINGS — edit this section for your dataset
# ============================================================

# Path to your data file (put "." if file is in working directory)
DATA_DIR  <- "."
DATA_FILE <- "qtl_data.csv"       # your input CSV file name

# Cross type: "bc" for backcross, "f2" for intercross
CROSS_TYPE <- "bc"

# Genotype codes in your file:
#   backcross (bc) : c("A", "H")
#   intercross(f2) : c("A", "H", "B")
GENO_CODES <- c("A", "H")

# Phenotype column number to analyse (1 = first phenotype)
PHENO_COL <- 1

# Significance level for QTL detection
ALPHA <- 0.05

# Number of permutations for significance threshold
# (1000 recommended for final analysis, 100 for quick testing)
N_PERM_SINGLE <- 1000
N_PERM_TWO    <- 100

# LOD drop for confidence intervals (1.5 is standard)
LOD_DROP <- 1.5

# Bayes interval probability (0.95 is standard)
BAYES_PROB <- 0.95

# Step size (cM) for genotype probability calculations
STEP_SIZE <- 2

# Genotyping error probability
ERROR_PROB <- 0.01


# ============================================================
# STEP 1: READ DATA
# ============================================================

cat("Reading cross data...\n")

cross <- read.cross(
  format    = "csv",
  dir       = DATA_DIR,
  file      = DATA_FILE,
  genotypes = GENO_CODES,
  crosstype = CROSS_TYPE
)

# Check your data loaded correctly:
#   - correct number of individuals?
#   - correct number of markers?
#   - correct number of phenotypes?
#   - percent genotyped should be high (>90%)
summary(cross)
plot(cross)


# ============================================================
# STEP 2: CALCULATE GENOTYPE PROBABILITIES
# ============================================================
# Uses Hidden Markov Model to estimate genotype probabilities
# at evenly spaced positions (every STEP_SIZE cM)
# Must be done before genome scanning

cross <- calc.genoprob(cross,
                       step       = STEP_SIZE,
                       error.prob = ERROR_PROB)


# ============================================================
# STEP 3: ESTIMATE GENETIC MAP
# ============================================================
# Re-estimates marker positions from your genotype data
# Compare original map (cross) vs re-estimated map (newmap)
# Large differences may indicate data quality issues

newmap <- est.map(cross, error.prob = ERROR_PROB, verbose = FALSE)
plot.map(cross, newmap,
         main = "Original Map vs Re-estimated Map")


# ============================================================
# STEP 4: CHECK FOR GENOTYPING ERRORS
# ============================================================
# Calculates error LOD scores for each genotype call
# High scores (>3) suggest likely genotyping errors
# These are not automatically removed — inspect and decide

cross <- calc.errorlod(cross, error.prob = ERROR_PROB)

cat("\nPotential genotyping errors (errorlod > 3):\n")
top.errorlod(cross, cutoff = 3)


# ============================================================
# STEP 5: SINGLE-QTL GENOME SCAN
# ============================================================
# Scans each position across genome for QTL signal
# Haley-Knott (hk) is fast; EM is slower but more accurate
# pheno.col selects which phenotype to analyse

out.hk <- scanone(cross, method = "hk", pheno.col = PHENO_COL)
out.em <- scanone(cross, method = "em", pheno.col = PHENO_COL)

# Plot both methods — they should largely agree
# Large disagreements may indicate model fit issues
plot(out.hk, out.em,
     col  = c("blue", "red"),
     ylab = "LOD score",
     main = "Single-QTL Genome Scan")
legend("topright",
       legend = c("Haley-Knott", "EM"),
       col    = c("blue", "red"),
       lwd    = 2)


# ============================================================
# STEP 6: PERMUTATION TEST — significance threshold
# ============================================================
# Randomly shuffles phenotype labels to build null distribution
# Determines LOD threshold for genome-wide significance
# More permutations = more accurate threshold

set.seed(123)  # for reproducibility
perm.hk <- scanone(cross,
                   method    = "hk",
                   pheno.col = PHENO_COL,
                   n.perm    = N_PERM_SINGLE,
                   verbose   = FALSE)

# Print LOD thresholds at 5% and 10% significance
cat("\nLOD significance thresholds:\n")
summary(perm.hk, alpha = c(0.05, 0.10))

# Plot scan with threshold lines
plot(out.hk,
     ylab = "LOD score",
     main = "Genome Scan with Significance Thresholds")
add.threshold(out.hk, perms = perm.hk,
              alpha = 0.05, col = "red",  lty = 2)
add.threshold(out.hk, perms = perm.hk,
              alpha = 0.10, col = "blue", lty = 2)
legend("topright",
       legend = c("LOD", "5% threshold", "10% threshold"),
       col    = c("black", "red", "blue"),
       lty    = c(1, 2, 2))


# ============================================================
# STEP 7: SUMMARIZE SIGNIFICANT QTLs
# ============================================================
# Lists all QTL peaks above the significance threshold
# These chromosomes and positions guide all remaining steps

cat("\nSignificant QTLs at alpha =", ALPHA, ":\n")
sig_qtls <- summary(out.hk,
                    perms    = perm.hk,
                    alpha    = ALPHA,
                    pvalues  = TRUE)
print(sig_qtls)

# Extract significant chromosomes automatically
# (used in Steps 8, 10, 11 below)
sig_chrs <- unique(sig_qtls$chr)
sig_pos  <- sig_qtls$pos
cat("\nSignificant chromosomes:", sig_chrs, "\n")
cat("Positions (cM):", sig_pos, "\n")


# ============================================================
# STEP 8: CONFIDENCE INTERVALS
# ============================================================
# Calculated automatically for all significant QTLs
# lodint  : region where LOD drops by LOD_DROP from peak
# bayesint: Bayesian credible interval

cat("\nConfidence intervals for significant QTLs:\n")

for (ch in sig_chrs) {

  cat("\n--- Chromosome", ch, "---\n")

  cat("LOD drop interval (drop =", LOD_DROP, "):\n")
  print(lodint(out.hk, chr = ch, drop = LOD_DROP))

  cat("Bayes interval (prob =", BAYES_PROB, "):\n")
  print(bayesint(out.hk, chr = ch, prob = BAYES_PROB))

}


# ============================================================
# STEP 9: TWO-QTL GENOME SCAN
# ============================================================
# Scans all pairs of positions for two-QTL models
# Detects epistasis (interaction between two QTLs)
# Takes longer than single-QTL scan

cat("\nRunning two-QTL genome scan...\n")
out2 <- scantwo(cross,
                method    = "hk",
                pheno.col = PHENO_COL,
                verbose   = FALSE)
plot(out2, main = "Two-QTL Genome Scan")

# Permutation test for two-QTL significance
set.seed(456)
perm2 <- scantwo(cross,
                 method    = "hk",
                 pheno.col = PHENO_COL,
                 n.perm    = N_PERM_TWO,
                 verbose   = FALSE)

cat("\nTwo-QTL permutation thresholds:\n")
summary(perm2, alpha = ALPHA)

cat("\nSignificant QTL pairs:\n")
summary(out2, perms = perm2, alpha = ALPHA)


# ============================================================
# STEP 10: FIT MULTIPLE-QTL MODEL
# ============================================================
# Only runs if 2+ significant QTLs were found in Step 7
# IMPORTANT: if you know better peak positions from Step 7,
# update sig_pos accordingly before running this step

if (length(sig_chrs) >= 2) {

  cat("\nFitting multiple-QTL model...\n")

  qtl_fit <- makeqtl(cross,
                     chr  = sig_chrs,
                     pos  = sig_pos,
                     what = "prob")

  # Build formula automatically: y ~ Q1 + Q2 + ... + Qn
  qtl_terms    <- paste0("Q", seq_along(sig_chrs))
  add_formula  <- as.formula(paste("y ~", paste(qtl_terms, collapse = " + ")))

  # Build interaction formula: y ~ Q1 + Q2 + Q1:Q2 etc.
  int_terms   <- combn(qtl_terms, 2, FUN = function(x) paste(x, collapse = ":"))
  int_formula <- as.formula(paste("y ~",
                                  paste(c(qtl_terms, int_terms), collapse = " + ")))

  # Additive model (no interactions)
  cat("\nAdditive model:\n")
  fit1 <- fitqtl(cross,
                 qtl      = qtl_fit,
                 formula  = add_formula,
                 method   = "hk",
                 get.ests = TRUE)
  summary(fit1)

  # Interaction model
  cat("\nInteraction model:\n")
  fit2 <- fitqtl(cross,
                 qtl      = qtl_fit,
                 formula  = int_formula,
                 method   = "hk",
                 get.ests = TRUE)
  summary(fit2)

  # Refine QTL peak positions using the additive model
  cat("\nRefined QTL positions:\n")
  rqtl <- refineqtl(cross,
                    qtl     = qtl_fit,
                    formula = add_formula,
                    method  = "hk",
                    verbose = FALSE)
  summary(rqtl)

} else if (length(sig_chrs) == 1) {

  cat("\nOnly one significant QTL found — skipping multiple-QTL model.\n")

} else {

  cat("\nNo significant QTLs found — skipping multiple-QTL model.\n")

}


# ============================================================
# STEP 11: EFFECT PLOTS
# ============================================================
# Shows mean phenotype value for each genotype class at each QTL
# Also shows epistasis plot for pairs of QTLs

# Run sim.geno first to avoid warnings in effectplot
cross <- sim.geno(cross,
                  step       = STEP_SIZE,
                  n.draws    = 64,
                  error.prob = ERROR_PROB)

if (length(sig_chrs) >= 1) {

  # Individual QTL effect plots
  par(mfrow = c(1, length(sig_chrs)))

  for (i in seq_along(sig_chrs)) {
    mname <- paste0(sig_chrs[i], "@", sig_pos[i])
    effectplot(cross,
               mname1 = mname,
               main   = paste("QTL Effect: Chr", sig_chrs[i],
                              "@", sig_pos[i], "cM"))
  }

  # Pairwise epistasis plots for all QTL pairs
  if (length(sig_chrs) >= 2) {

    par(mfrow = c(1, 1))
    qtl_pairs <- combn(seq_along(sig_chrs), 2)

    for (k in seq_len(ncol(qtl_pairs))) {
      i <- qtl_pairs[1, k]
      j <- qtl_pairs[2, k]
      m1 <- paste0(sig_chrs[i], "@", sig_pos[i])
      m2 <- paste0(sig_chrs[j], "@", sig_pos[j])
      effectplot(cross,
                 mname1 = m1,
                 mname2 = m2,
                 main   = paste("Epistasis: Chr", sig_chrs[i],
                                "x Chr", sig_chrs[j]))
    }
  }
}

cat("\nAnalysis complete.\n")
