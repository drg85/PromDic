# Parameters (change as desired)
N <- 1000
p <- 0.15 # proportion of 1's
imax <- 1000    # number of simulations
h <- 3.841      # chi-square(1) threshold

nOnes  <- round(p * N)
nZeros <- N - nOnes

# Track whether each run has a FP
hasFP   <- logical(imax)   # any G > 0
hasFPEx <- logical(imax)   # any G > 0 & all counts >= 5

for (sim in seq_len(imax)) {
  # Generate sequential random vector (biased walk)
  x <- numeric(N + 1)  # cumulative zeros
  y <- numeric(N + 1)  # cumulative ones
  remainingZeros <- nZeros
  remainingOnes  <- nOnes
  
  for (k in seq_len(N)) {
    px <- remainingZeros / (remainingZeros + remainingOnes)
    if (runif(1) < px) {
      # place a zero
      x[k + 1] <- x[k] + 1
      y[k + 1] <- y[k]
      remainingZeros <- remainingZeros - 1
    } else {
      # place a one
      x[k + 1] <- x[k]
      y[k + 1] <- y[k] + 1
      remainingOnes <- remainingOnes - 1
    }
  }
  
  # Evaluate G-statistic at all split points
  for (k in seq_len(N - 1)) {
    a <- y[k + 1]        # ones on the left
    b <- x[k + 1]        # zeros on the left
    c <- nOnes  - a      # ones on the right
    d <- nZeros - b      # zeros on the right
    
    g <- (N * (a * d - b * c)^2) / ((a + b) * (c + d) * (a + c) * (b + d)) - h
    minv <- min(c(a, b, c, d))
    
    if (g > 0) {
      hasFP[sim] <- TRUE
      if (minv >= 5) {
        hasFPEx[sim] <- TRUE
      }
    }
  }
}

# Results
totalFP   <- sum(hasFP)
totalFPEx <- sum(hasFPEx)

cat(sprintf("Out of %d runs:\n", imax))
cat(sprintf(" (a) Runs with any G>0: %d (%.2f%%)\n", 
            totalFP, 100 * totalFP / imax))
cat(sprintf(" (b) Runs with any G>0 & all counts>=5: %d (%.2f%%)\n", 
            totalFPEx, 100 * totalFPEx / imax))
