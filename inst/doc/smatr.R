## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE
)

## ----setup, eval = FALSE------------------------------------------------------
# # Install via CRAN or github
# install.packages("smatr")
# devtools::install_github("dfalster/smatr3")

## ----message=FALSE------------------------------------------------------------
# Load library
library(smatr)

## -----------------------------------------------------------------------------
data(leaflife)
data(leafmeas)

#For more details on each data set
?leaflife
?leafmeas

## -----------------------------------------------------------------------------
# Extract only low-nutrient, low-rainfall data from leaflife dataset:
leaf.low <- subset(leaflife, soilp == 'low' & rain == 'low')

# Fit a single MA for log(leaf longevity) vs log(leaf mass per area):
ma(longev ~ lma, log='xy', data=leaflife)

## -----------------------------------------------------------------------------
# Test if the MA slope is not significantly different from 1 for longevity and leaf mass per area (LMA):
ma_obj <- ma(longev ~ lma, log='xy', slope.test=1, data=leaflife) 

summary(ma_obj)
plot(ma_obj)

## -----------------------------------------------------------------------------
# First we subset the data that is relevant
leaf.low.soilp <- subset(leaflife, soilp == 'low')

# Then using * we can fit SMAs separately at each of high and low rainfall sites and test for whether sites with different rainfall share a common slope
rain_sma_obj <- sma(longev~lma*rain, log="xy", data=leaf.low.soilp)

# Lets plot longevity vs LMA separately for each group first:
plot(rain_sma_obj)

## ----message=FALSE------------------------------------------------------------
library(ggplot2)

# ggplot2 version of the grouped SMA fit
ggplot(rain_sma_obj)

# The result is a normal ggplot object, so extend it however you like:
ggplot(rain_sma_obj) +
  theme_minimal() +
  labs(title = "Leaf longevity vs LMA", colour = "Rainfall")

## -----------------------------------------------------------------------------
# Fit SMAs for each group and tests if several SMA lines share a common slope. 
rain_sma_obj <- sma(longev~lma*rain, log="xy", data=leaf.low.soilp)

summary(rain_sma_obj)

## -----------------------------------------------------------------------------
# Fit SMAs separately at each of high and low rainfall sites,
# and test if there is a common slope and whether it is equal to 1:
rain_slope_obj <- sma(longev~lma*rain, log="xy", slope.test=1, data=leaf.low.soilp)

summary(rain_slope_obj)

## -----------------------------------------------------------------------------
# Fit SMAs separately at each of high and low rainfall sites, and test whether sites differ in the elevation of their SMA
rain_elev_obj <- sma(longev~lma*rain, log="xy", type = "elevation", data=leaf.low.soilp)

summary(rain_elev_obj)

## -----------------------------------------------------------------------------
# Fit SMAs with common slope across each of high and low rainfall sites, and test for no shift along common SMA:
rain_shift_obj <- sma(longev~lma+rain, log="xy", type="shift", data=leaf.low.soilp)
summary(rain_shift_obj)

