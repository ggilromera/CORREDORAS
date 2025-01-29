# Set working directory dynamically to the script's location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Ensure the directory path is correctly set where the .csv file is stored

## Install and load required packages
required_packages <- c("clam", "rbacon", "dplyr", "tidyverse")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, library, character.only = TRUE)

# Define parameters for rbacon
params <- list(
  thick = 50,         # Thickness of sediment slices (in cm) used in the age-depth model
  d.by = 1,          # Depth interval (in cm) at which ages are estimated
  core.name = "bsm24_v1",  # Name of the core being analyzed
  acc.mean = 10,      # Mean accumulation rate (cm per year)
  acc.shape = 1.5,    # Shape parameter of the accumulation rate distribution
  mem.strength = 0.8, # Strength of memory in accumulation rate changes
  mem.mean = 0.5      # Mean memory (how much accumulation rate depends on previous depth)
)

# Run Bacon model
x <- Bacon(
  params$core.name,   # Core name should be the first argument
  coredir = getwd(),
  acc.mean = params$acc.mean,
  acc.shape = params$acc.shape,
  mem.mean = params$mem.mean,
  mem.strength = params$mem.strength,
  thick = params$thick,
  d.by = params$d.by,
  ask = FALSE,
  suggest = FALSE,
  plot.pdf = TRUE,
  depths.file = FALSE,
  normal = FALSE,
  rotate.axes = TRUE
)


## Generate accumulation rate plots
accrate.age.ghost()   # Estimates yr/cm in an age scale
accrate.depth.ghost() # Estimates yr/cm in a depth scale

## Calibrating original dates
# To calibrate radiocarbon dates one by one in Rbacon, use the calibrate() function.
# Example: 
# calibration_result <- calibrate(5000, 30, calCurves = "IntCal20")
# This will calibrate a radiocarbon date of 5000 years BP with an uncertainty of 30 years using the IntCal20 calibration curve.
# You can then plot the calibration:
# plot(calibration_result)

# Define functions for weighted mean and standard deviation
wt.mean <- function(x, wt) {
  valid <- which(is.finite(x * wt))
  wt <- wt[valid]
  x <- x[valid]
  sum(wt * x) / sum(wt)
}

wt.sd <- function(x, wt) {
  valid <- which(is.finite(x + wt))
  wt <- wt[valid]
  x <- x[valid]
  xbar <- wt.mean(x, wt)
  sqrt(sum(wt * (x - xbar)^2) * (sum(wt) / (sum(wt)^2 - sum(wt^2))))
}

# Extract calibration data
calibration <- info$calib$probs

# Retrieve ages and probabilities from the first calibration matrix
ages <- calibration[[1]][, 1]
probs <- calibration[[1]][, 2]

# Plot calibration curve
plot(ages, probs, type = "l")

# Compute weighted mean and standard deviation for each date
calibration.mean <- sapply(calibration, function(x) wt.mean(x[, 1], x[, 2]))
calibration.sd <- sapply(calibration, function(x) wt.sd(x[, 1], x[, 2]))

# Read dates data frame and add calibration values
df <- read.csv(
  file = file.path(getwd(), "bsm24_v1", "bsm24_v1.csv"),
  header = TRUE
) %>% 
  mutate(
    calibration.mean = calibration.mean,
    calibration.sd = calibration.sd
  )

# Save the updated dataframe
output_dir <- file.path(getwd(), "bsm24_v1")
setwd(output_dir)
write.csv(df, file = "bsm24_v1_cal.csv", row.names = FALSE)

