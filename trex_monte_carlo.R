# Input values to generate a distribution of expected number of T. rex that ever lived

# ----------------------
# --- Body Mass (g) ----
# ----------------------
BMmean <- 6790400
BMsdev <- 665753.2823

AsymBMmean <- 7090
AsymBMsdev <- 1012.755

# ------------------------------
# ---- Pop density (ind/km^2)---
# ---- value for intercept ----- 
# ---- for Damuth's rule eq ----
# ------------------------------

logA_mid <- 2.991
logA_sd <- 0.6060698

# -------------------------------
# ---- Sexual maturity (yrs) ----
# -------------------------------
SMmean <- 15.5 
SMsdev <- 0.7653

# --------------------------------------
# --- Life history curve (gompertz) ----
# --------------------------------------
#lx <- c(0.998, 0.995, 0.992, 0.987, 0.982, 0.975, 0.967, 0.957,
#        0.944, 0.929, 0.910, 0.887, 0.859, 0.826, 0.786, 0.739,
#        0.684, 0.621, 0.550, 0.473, 0.393, 0.311, 0.232, 0.161,
#        0.102, 0.058, 0.029, 0.012)

  # Weibull lx values
  # lx <-   c(0.993, 0.987, 0.980, 0.974, 0.967, 0.961, 0.954, 0.947,
  #           0.939, 0.930, 0.919, 0.904, 0.884, 0.858, 0.823, 0.776,
  #           0.717, 0.643, 0.556, 0.458, 0.355, 0.255, 0.166, 0.097,
  #           0.049, 0.021, 0.007, 0.002)

# --------------------------------
# ---- Temporal range (M.yrs) ----
# --------------------------------

TRmin <- 1.2
TRmax <- 3.6

# --------------------------------
# --- Geographic range (M.km^2) --
# --------------------------------

GRmean <- 2.3
GRsdev <- 0.449


s#ource("functions.R")
source("functions_2.R")

n <- 1000000

#list_trex <- trex_ever_lived(n, BMmean, BMsdev, logA_mid, logA_sd, SMmean, SMsdev, lx, TRmin, TRmax, GRmean, GRsdev)
timestamp()
list_trex2 <- trex_ever_lived2(n, AsymBMmean, AsymBMsdev, logA_mid, logA_sd, SMmean, SMsdev, TRmin, TRmax, GRmean, GRsdev)
timestamp()

trex_ranges <- t(sapply(list_trex2, quantile, c(0.025, 0.5,0.975), na.rm = T))
lapply(list_trex2, function(x){sum(is.na(x))})

save(file = "Monte_carlo_trex.RData", list_trex2, trex_ranges)

timestamp()
sensitivity <- trex_ever_lived2(damuth_slope = -2/3,n, AsymBMmean, AsymBMsdev, logA_mid, logA_sd, SMmean, SMsdev, TRmin, TRmax, GRmean, GRsdev)
timestamp()
sensitivity_trex_ranges <- t(sapply(sensitivity, quantile, c(0.025, 0.5,0.975), na.rm = T))
lapply(sensitivity, function(x){sum(is.na(x))})
save(file = "Monte_carlo_trex_damuth_sensibility.RData", sensitivity_trex_ranges, sensitivity)
