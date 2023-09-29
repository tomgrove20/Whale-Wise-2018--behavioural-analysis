### EXPLORING TEMPORAL AUTOCORRELATION
### DIVE TIME AIS EXAMPLE

### 21/05/2023
### Tom Grove
### tomgrove20@yahoo.co.uk

### PACKAGES
packages <- c("tidyverse", "mgcv", "itsadug")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

#---------------- FUNCTIONS + THEME --------------------

source("code/functions.R")
source("code/themes.R")

#---------------- DATA --------------------

# dive time foc model
mod <- readRDS("intermediate-products/gam-v2/dive-time/ais/gam_dive-time_ais_final_v2.rds")

# total data set
tot <- read.csv("intermediate-products/response-var-dfs/dive-time_ais.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding.follow", "surface.active.follow"), as.factor) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(rib.num.30.1500<=4) # only 1/1243 rows had 5 ribs in 30.1500

# no transformations needed

#---------------- GLOBAL ACF PLOTS --------------------

# ACF
lag = 15
auto <- acf(resid(mod), lag.max = lag, plot = FALSE)
dat <- data.frame(val = auto$acf, lag = auto$lag)
alpha <- 0.95; conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(auto$n.used) # extracting significance interval

dat %>% ggplot(aes(x=lag, y = val)) +
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  geom_hline(yintercept = 0, col = "black") +
  labs(y="ACF", x="Lag", title= "All data") +
  geom_segment(aes(xend=lag, yend=0)) +geom_point()

# PACF
lag = 15
auto <- pacf(resid(mod), plot = FALSE, lag.max = lag)
dat <- data.frame(val = auto$acf, lag = auto$lag)
alpha <- 0.95; conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(auto$n.used) # extracting significance interval

dat %>% ggplot(aes(x=lag, y = val)) +
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  geom_hline(yintercept = 0, col = "black") +
  labs(y="PACF", x="Lag", title= "All data") +
  geom_segment(aes(xend=lag, yend=0)) +geom_point()


#---------------- FOLLOW-SPECIFIC ACF PLOTS --------------------

# now we need to calculate the four longest follows
fols <- as.character((tot %>% group_by(folnum.unique) %>% summarise(n = n()) %>% 
                        arrange(-n) %>%  slice(1:8))$folnum.unique) 

# ACF
# extract data for each fol
dat.tot <- data.frame()
for (i in fols) {
  lag = 15
  auto <- acf(resid(mod)[tot$folnum.unique==i], lag.max = lag, plot = FALSE)
  alpha <- 0.95; conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(auto$n.used) # extracting significance interval
  dat <- data.frame(val = auto$acf, lag = auto$lag, fol = i, lwr = conf.lims[1], upr = conf.lims[2])
  dat.tot <- bind_rows(dat.tot, dat)
}

dat.tot %>% ggplot(aes(x=lag, y = val)) +
  geom_hline(aes(yintercept=lwr), lty=2, col='blue') +
  geom_hline(aes(yintercept=upr), lty=2, col='blue') +
  geom_hline(aes(yintercept=0), col = "black") +
  labs(y="ACF", x="Lag") +
  geom_segment(aes(xend=lag, yend=0)) +geom_point() +
  facet_wrap(~fol, ncol = 2)

# PACF
# extract data for each fol
dat.tot <- data.frame()
for (i in fols) {
  lag = 15
  auto <- pacf(resid(mod)[tot$folnum.unique==i], lag.max = lag, plot = FALSE)
  alpha <- 0.95; conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(auto$n.used) # extracting significance interval
  dat <- data.frame(val = auto$acf, lag = auto$lag, fol = i, lwr = conf.lims[1], upr = conf.lims[2])
  dat.tot <- bind_rows(dat.tot, dat)
}

dat.tot %>% ggplot(aes(x=lag, y = val)) +
  geom_hline(aes(yintercept=lwr), lty=2, col='blue') +
  geom_hline(aes(yintercept=upr), lty=2, col='blue') +
  geom_hline(aes(yintercept=0), col = "black") +
  labs(y="PACF", x="Lag") +
  geom_segment(aes(xend=lag, yend=0)) +geom_point() +
  facet_wrap(~fol, ncol = 2)


