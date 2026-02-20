# ============================================
# R/qtl: QTL ANALYSIS FROM EXISTING FILES
# ============================================

library(qtl)

# ── Set your working directory if needed ───────────────────────────────────
# setwd("C:/your/path/here")   # uncomment and edit if required
# getwd()                       # confirm current working directory


# ============================================
# STEP 1: READ DATA
# ============================================

cross <- read.cross(
  format    = "csv",
  dir       = ".",
  file      = "qtl_data.csv",
  genotypes = c("A", "H"),
  crosstype = "bc"
)

summary(cross)
plot(cross)


# ============================================
# STEP 2: CALCULATE GENOTYPE PROBABILITIES
# ============================================

cross <- calc.genoprob(cross, step = 2, error.prob = 0.01)


# ============================================
# STEP 3: ESTIMATE GENETIC MAP
# ============================================

newmap <- est.map(cross, error.prob = 0.01, verbose = FALSE)
plot.map(cross, newmap)


# ============================================
# STEP 4: CHECK FOR GENOTYPING ERRORS
# ============================================

cross <- calc.errorlod(cross, error.prob = 0.01)
top.errorlod(cross, cutoff = 3)


# ============================================
# STEP 5: SINGLE-QTL GENOME SCAN
# ============================================

out.hk <- scanone(cross, method = "hk")
out.em <- scanone(cross, method = "em")

plot(out.hk, out.em,
     col  = c("blue", "red"),
     ylab = "LOD score",
     main = "Single-QTL Genome Scan")
legend("topright",
       legend = c("Haley-Knott", "EM"),
       col    = c("blue", "red"),
       lwd    = 2)


# ============================================
# STEP 6: PERMUTATION TEST
# ============================================

set.seed(456)
perm.hk <- scanone(cross, method = "hk", n.perm = 1000, verbose = FALSE)

summary(perm.hk, alpha = c(0.05, 0.10))

plot(out.hk,
     ylab = "LOD score",
     main = "Genome Scan with Significance Thresholds")
add.threshold(out.hk, perms = perm.hk, alpha = 0.05, col = "red",  lty = 2)
add.threshold(out.hk, perms = perm.hk, alpha = 0.10, col = "blue", lty = 2)
legend("topright",
       legend = c("LOD", "5% threshold", "10% threshold"),
       col    = c("black", "red", "blue"),
       lty    = c(1, 2, 2))


# ============================================
# STEP 7: SUMMARIZE SIGNIFICANT QTLs
# ============================================

summary(out.hk, perms = perm.hk, alpha = 0.05, pvalues = TRUE)


# ============================================
# STEP 8: CONFIDENCE INTERVALS
# ============================================

lodint(out.hk,   chr = 1, drop = 1.5)
lodint(out.hk,   chr = 5, drop = 1.5)

bayesint(out.hk, chr = 1, prob = 0.95)
bayesint(out.hk, chr = 5, prob = 0.95)


# ============================================
# STEP 9: TWO-QTL GENOME SCAN
# ============================================

out2 <- scantwo(cross, method = "hk", verbose = FALSE)
plot(out2, main = "Two-QTL Genome Scan")

set.seed(789)
perm2 <- scantwo(cross, method = "hk", n.perm = 100, verbose = FALSE)
summary(perm2, alpha = 0.05)
summary(out2,  perms = perm2, alpha = 0.05)


# ============================================
# STEP 10: FIT MULTIPLE-QTL MODEL
# ============================================

qtl_fit <- makeqtl(cross, chr = c(1, 5), pos = c(50, 30), what = "prob")

# Additive model
fit1 <- fitqtl(cross,
               qtl      = qtl_fit,
               formula  = y ~ Q1 + Q2,
               method   = "hk",
               get.ests = TRUE)
summary(fit1)

# Interaction model
fit2 <- fitqtl(cross,
               qtl      = qtl_fit,
               formula  = y ~ Q1 + Q2 + Q1:Q2,
               method   = "hk",
               get.ests = TRUE)
summary(fit2)

# Refine QTL positions
rqtl <- refineqtl(cross,
                  qtl     = qtl_fit,
                  formula = y ~ Q1 + Q2,
                  method  = "hk",
                  verbose = FALSE)
summary(rqtl)


# ============================================
# STEP 11: EFFECT PLOTS
# ============================================

cross <- sim.geno(cross, step = 2, n.draws = 64, error.prob = 0.01)

par(mfrow = c(1, 2))
effectplot(cross, mname1 = "1@50", main = "QTL Effect: Chr 1 @ 50 cM")
effectplot(cross, mname1 = "5@30", main = "QTL Effect: Chr 5 @ 30 cM")

par(mfrow = c(1, 1))
effectplot(cross,
           mname1 = "1@50",
           mname2 = "5@30",
           main   = "Epistasis: Chr1 x Chr5")
