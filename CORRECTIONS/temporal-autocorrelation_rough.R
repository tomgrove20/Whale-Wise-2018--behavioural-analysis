##########################################
### EXPLORING TEMPORAL AUTOCORRELATION ###
##########################################
### 21/05/2023
### Tom Grove
### tomgrove20@yahoo.co.uk

### PACKAGES
packages <- c("tidyverse", "mgcv", "itsadug")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

#---------------- DATA --------------------

# dive time foc model
mod <- readRDS("intermediate-products/gam/dive-time/foc/gam_dive-time_foc_final.rds")

acf(residuals(mod),main="raw residual ACF")

# dive time ais model
mod <- readRDS("intermediate-products/gam/dive-time/ais/gam_dive-time_ais_final.rds")

acf(resid(mod), lag.max = 20, main = "pACF")
pacf(resid(mod), lag.max = 20, main = "pACF")

# ibi ais model
mod <- readRDS("intermediate-products/gam/ibi-mean/ais/gam_ibi-mean_ais_final.rds")

acf(resid(mod), lag.max = 20, main = "pACF")
pacf(resid(mod), lag.max = 20, main = "pACF")

# swim speed ais
mod <- readRDS("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_final.rds")

acf(resid(mod)[tot$folnum.unique==565], lag.max = 15, main = "pACF")
pacf(resid(mod)[tot$folnum.unique==565], lag.max = 20, main = "pACF")

# di ais
mod <- readRDS("intermediate-products/gam/di/ais/gam_di_ais_final.rds")

acf(resid(mod), lag.max = 20, main = "pACF")
pacf(resid(mod), lag.max = 20, main = "pACF")

# surface feeding
mod <- readRDS("intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_final.rds")
acf(resid(mod), lag.max = 20, main = "pACF")
pacf(resid(mod), lag.max = 20, main = "pACF")
