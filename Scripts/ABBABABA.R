# Read ABBABABA2 output (has header with CHR BLOCKstart BLOCKend Numer Denom ...)
dat <- read.table("Dgr_ABBABABA2_P1Mog_P2Sky_P3Mex_out.abbababa2",
                  header = TRUE)

# Check column names to confirm Numer / Denom
colnames(dat)[1:10]

num <- dat$Numer
den <- dat$Denom

# Full-data D
D_full <- sum(num) / sum(den)

# Jackknife over blocks
B <- length(num)
D_jk <- numeric(B)

for (i in 1:B) {
  D_jk[i] <- (sum(num) - num[i]) / (sum(den) - den[i])
}

D_bar <- mean(D_jk)
var_D <- (B - 1) / B * sum((D_jk - D_bar)^2)
SE_D <- sqrt(var_D)
Z <- D_full / SE_D

cat("D =", D_full, "SE =", SE_D, "Z =", Z, "\n")
