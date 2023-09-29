### Demo of secr
### Tom Grove
### 14/02/2020

library(secr)

# Link: https://cran.r-project.org/web/packages/secr/secr.pdf#page=214&zoom=100,133,553
# Another example: https://www.otago.ac.nz/density/pdfs/secr-tutorial.pdf

## navigate to folder with raw data files
olddir <- setwd (system.file("extdata", package="secr"))

# Look at trap data
traps <- read.table("trap.txt")
traps # 100 traps, each with associated coordinates
# col 1 = trapID. col 2 = X. col 3 = Y

# Look at capture data
capt <- read.table("capt.txt")
capt
# col 1 = session. col 2 = ID. col 3 = occasion. col 4 = detector. Unsure of 5/6
# Cols you need: session, ID, occasion, detector

## construct capthist object from raw data
captdata <- read.capthist ("capt.txt", "trap.txt", fmt = "XY", detector = "single")

## Let's look at the data
summary(captdata)
# n: number of distinct individuals detected on each occasion t
# u: number of individuals detected for the first time on each occasion t
# f: number of individuals detected on exactly t occasions
# M(t+1): cumulative number of detected individuals on each occasion t

plot (captdata, tracks = TRUE)

# Let's look at movement in a graph
m <- unlist(moves(captdata))
hist(m, breaks = seq(-61/2, 500,61), xlab = "Movement m", main = "")

# For a quick and biased estimate of sigma (spatial scale)
initialsigma <- RPSV(captdata, CC = TRUE)
cat("Quick and biased estimate of sigma =", initialsigma, "m\n")


## generate demonstration fits
secrdemo.0 <- secr.fit (captdata) # Null model 
secrdemo.CL <- secr.fit (captdata, CL = TRUE)
secrdemo.b <- secr.fit (captdata, model = list(g0 ~ b)) # fits a learned trap response. b = behaviour

## restore previous setting
setwd(olddir)

## display the null model fit, using the print method for secr
secrdemo.0
# Last table tells us that estimated density is 5.47 animals per hectare
# 95% CI: 4.35 - 6.90 animals per hectare
# For this simple model, there is one beta estimate for each parameter

# g0 and sigma jointly determine the detection function that you can easily plot with 95% CIs
# Plotting the null model
plot(secrdemo.0, limits = TRUE)


## compare fit of models
AIC(secrdemo.0, secrdemo.b)

## display estimates for the two models (single session)
collate(secrdemo.0, secrdemo.b)[1,,,]
