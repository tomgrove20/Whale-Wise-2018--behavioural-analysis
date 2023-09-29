##################################################
### PLOTS FOR HUSAVIK WHALE-WATCHING COMPANIES ###
##################################################
### 04/10/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

# 1. transformations
# 2. Collinearity
# 3. time frame selection
# 4. Final GAM (GCV + backward selection)
# 5. Results + plotting


### PACKAGES
packages <- c("tidyverse", "ppcor","RColorBrewer", "scales","survMisc", "Metrics", "corrplot", "lme4", "MASS", "mgcv", "tidymv", "mgcViz", "gridExtra", "gratia", "ggcorrplot", "cowplot", "ggpubr")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)


#---------------- FUNCTIONS + THEME --------------------

source("code/functions.R")
source("code/themes.R")
isnt.na = function(x){!is.na(x)}

#############
### FOCAL ###
#############

#---------------- BREATHS PER SURFACING: FOCAL --------------------

# focal vessel
tot <- read.csv("intermediate-products/response-var-dfs/surface-interval_focal.csv") %>%
  filter(!is.na(ibi.mean) & !is.na(dist.mean)) %>%
  mutate_at(vars(contains("dt")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate.last", "year", "folnum.unique"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1)  %>% # important dummy variable to factor out random effect when predicting response
  # dplyr::select(-contains(".10")) %>%
  # dplyr::select(-contains(c("0_end", "0_start"))) %>%
  rename_at(vars(contains("_mean")), funs(str_replace(.,"_mean",""))) %>%
  filter(vess.accel.min.300>-10) %>% # temporary filtering until boat stuff properly processes
  filter(vess.accel.max.300<0.5) %>% # filtering out until boat stuff properly processes
  # creating new transformed columns
  mutate(across(contains(c("dive.time","dist")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # sqrt mean dist
         across(contains("speed.var"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt dist max
         across(contains("vess.di."), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))

gam <- readRDS("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_final.rds")

# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.num = c("day.of.year","log.dive.time.before",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "vess.accel.max.60", "vess.accel.60") # numeric 
var.rand = c("folnum.unique") # random 


### VESSEL SPEED (60 sec)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "vess.speed.60", var.fac, var.num, var.rand) %>% 
    mutate(vess.speed.60 = vess.speed.60*2.23694) %>%
    ggplot(aes(x = vess.speed.60, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel speed (mph, previous 1 minute)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin)
ggsave("intermediate-products/presentation-plots/companies-202210/speed/companies_breath-num_speed.png", height = 4.5, width = 6.5, bg = "transparent")

### VESSEL SPEED SD (60 sec)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "sqrt.vess.speed.var.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.vess.speed.var.60^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Variation in boat speed (previous 1 minute)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt"))
ggsave("intermediate-products/presentation-plots/companies-202210/speed-var/companies_breath-num_speed-var.png", height = 4.5, width = 6.5, bg = "transparent")

### VESSEL DI (60 sec)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "arcsin.vess.di.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sin(arcsin.vess.di.60), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Boat path straightness (previous 1 minute)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "asn", breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99)))
ggsave("intermediate-products/presentation-plots/companies-202210/di/companies_breath-num_di.png", height = 4.5, width = 6.5, bg = "transparent")


#---------------- IBI: FOCAL --------------------

# focal vessel
tot <- read.csv("intermediate-products/response-var-dfs/surface-interval_focal.csv") %>%
  filter(!is.na(ibi.mean) & !is.na(dist.mean)) %>%
  mutate_at(vars(contains("dt")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate.last", "year", "folnum.unique"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1)  %>% # important dummy variable to factor out random effect when predicting response
  # dplyr::select(-contains(".10")) %>%
  # dplyr::select(-contains(c("0_end", "0_start"))) %>%
  rename_at(vars(contains("_mean")), funs(str_replace(.,"_mean",""))) %>%
  filter(vess.accel.min.300>-10) %>% # temporary filtering until boat stuff properly processes
  filter(vess.accel.max.300<0.5) %>% # filtering out until boat stuff properly processes
  # creating new transformed columns
  mutate(across(contains(c("dive.time","dist")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # sqrt mean dist
         across(contains("speed.var"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt dist max
         across(contains("vess.di."), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))

gam <- readRDS("intermediate-products/gam/ibi-mean/foc/gam_ibi-mean_foc_final.rds")

# defining variable combinations
var.num = c("log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.300")
var.rand = c("folnum.unique")
var.fac = NA

### DISTANCE
(pr <- pred.var(gam, tot, resp = "ibi.mean", var = "log.dist.mean", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.dist.mean), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Distance to whale (meters)", y = "Time between breaths (seconds)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(50,100,200,500,1000, 2000)) +
    scale_y_continuous(breaks = seq(0,40,2)))
ggsave("intermediate-products/presentation-plots/companies-202210/dist/companies_ibi-mean_dist.mean.png", height = 4.5, width = 6.5, bg = "transparent")

### SPEED
(pr <- pred.var(gam, tot, resp = "ibi.mean", var = "vess.speed.60", var.fac, var.num, var.rand) %>% 
    mutate(vess.speed.60 = vess.speed.60*2.23694) %>%
    ggplot(aes(x = vess.speed.60, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Boat speed (mph, previous 1 minute)", y = "Time between breaths (seconds)") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,40,2)))
ggsave("intermediate-products/presentation-plots/companies-202210/speed/companies_ibi-mean_speed.png", height = 4.5, width = 6.5, bg = "transparent")


### SPEED VARIATION 
(pr <- pred.var(gam, tot, resp = "ibi.mean", var = "sqrt.vess.speed.var.300", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.vess.speed.var.300^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Variation in boat speed (previous 5 minutes)", y = "Time between breaths (seconds)") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,40,2)) + scale_x_continuous(trans = "sqrt"))
ggsave("intermediate-products/presentation-plots/companies-202210/speed-var/companies_ibi-mean_speed-var.png", height = 4.5, width = 6.5, bg = "transparent")


#---------------- DIVE TIME: FOCAL --------------------

# focal vessel
tot <- read.csv("intermediate-products/response-var-dfs/dive-time_focal.csv") %>%
  dplyr::select(-c("lat.whale", "lon.whale")) %>% drop_na() %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding.follow", "surface.active.follow"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(vess.accel.max.300<0.5) %>% # filtering out until boat stuff properly processes
  mutate(across(contains("speed.var"), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # log speed var
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"), # arcsin DI
         log.dist = log(dist))

gam <- readRDS("intermediate-products/gam/dive-time/foc/gam_dive-time_foc_final.rds")

# final model
var.fac = c("year") # factor
var.num = c("day.of.year", "log.dist", "vess.speed.300",  "log.vess.speed.var.300", "vess.accel.60") # numeric 
var.rand = c("folnum.unique") # random 

### DISTANCE
(pr <- pred.var(gam, tot, resp = "dive.time", var = "log.dist", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.dist), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Distance (m)", y = "Dive time (seconds)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(20,50,100,200,400,1000)))
ggsave("intermediate-products/presentation-plots/companies-202210/dist/companies_dive-time_dist.png", height = 4.5, width = 6.5, bg = "transparent")

### SPEED
(pr <- pred.var(gam, tot, resp = "dive.time", var = "vess.speed.300", var.fac, var.num, var.rand) %>% 
    mutate(vess.speed.300 = vess.speed.300*2.23694) %>%
    ggplot(aes(x = vess.speed.300, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel speed (mph, previous 5 minutes)", y = "Dive time (sec)") + partial.theme + partial.grid.margin)
ggsave("intermediate-products/presentation-plots/companies-202210/speed/companies_dive-time_speed.png", height = 4.5, width = 6.5, bg = "transparent")

### SPEED VARIATION (5 mins)
(pr <- pred.var(gam, tot, resp = "dive.time", var = "log.vess.speed.var.300", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.vess.speed.var.300), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Variation in boat speed (previous 5 minutes)", y = "Dive time (seconds)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(0.1, 0.2, 0.5, 1, 2)))
ggsave("intermediate-products/presentation-plots/companies-202210/speed-var/companies_dive-time_speed.var.png", height = 4.5, width = 6.5, bg = "transparent")


#---------------- DI: FOCAL --------------------

tot <- read.csv("intermediate-products/response-var-dfs/di_focal.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding", "surface.active"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(difftime.next<2000) %>% # until properly processed
  filter(vess.accel.max.300 < 2) %>% # filtering until properly processed
  filter(vess.accel.60 >-0.2) %>% # filtering until properly processed
  filter(vess.speed.var.60<3) %>% # filtering until properly processed
  mutate(arcsin.DI = asin(DI)) %>%
  mutate(across(contains(c("difftime","vess.speed.var")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(ends_with(c("speed.60", "speed.300", "dist")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))

gam <- readRDS("intermediate-products/gam/di/foc/gam_di_foc_final.rds")

var.fac = c("group", "seastate")
var.num = c("day.of.year", "sqrt.dist", "log.difftime.prev", "log.difftime.next", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.300", "vess.accel.max.60")
var.rand = c("folnum.unique")


### DISTANCE
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "sqrt.dist", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.dist^2, y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Distance to whale (meters)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(20,50,100,200,400,800)) +
    scale_y_continuous(trans = asin.trans, breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)))
ggsave("intermediate-products/presentation-plots/companies-202210/dist/companies_di_dist.png", height = 4.5, width = 6.5, bg = "transparent")


### SPEED
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "sqrt.vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = (sqrt.vess.speed.60^2)*2.23694, y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Boat speed (mph, previous 1 minute)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.9, 0.95, 0.98, 0.99, 1)) +
    scale_x_continuous(trans = "sqrt", breaks = c(0.2, 1,2,4,8,12)))
ggsave("intermediate-products/presentation-plots/companies-202210/speed/companies_di_speed.png", height = 4.5, width = 6.5, bg = "transparent")

### ACCELERATION
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "vess.accel.max.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = vess.accel.max.60, y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = bquote("Maximum boat accleration (previous 1 minute)"), y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.7, 0.8, 0.9, 0.95, 0.99, 1)))
ggsave("intermediate-products/presentation-plots/companies-202210/accel/companies_di_accel.png", height = 4.5, width = 6.5, bg = "transparent")




#---------------- SPEED --------------------

# foc
tot <- read.csv("intermediate-products/response-var-dfs/speed_focal.csv") %>%
  mutate_at(vars(contains("datetime")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.active", "surface.feeding"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  filter_at(vars(contains("vess")), isnt.na) %>%
  mutate(dummy = 1) %>%
  # important dummy variable to factor out random effect when predicting response
  filter(vess.accel.max.300 < 5) %>% # filtering until properly processed
  filter(vess.accel.60 >-0.2) %>% # filtering until properly processed
  filter(vess.speed.var.60<3) %>% # filtering until properly processed
  mutate(across(contains(c("diff.time","vess.speed.var")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(ends_with(c("speed.60", "speed.300", "dist")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))

gam <- readRDS("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_final.rds")

## PARTIALS + PREDICTED
var.fac = c("seastate", "year", "surface.feeding", "surface.active")
var.num = c("day.of.year", "sqrt.dist", "log.diff.time", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.60", "vess.accel.60")
var.rand = c("folnum.unique")


### DISTANCE
(pr <- pred.var(gam, tot, resp = "speed", var = "sqrt.dist", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.dist^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Distance to whale (meters)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(20,50,100,200,400,800)) +
    scale_y_continuous(breaks = seq(1,10, by = 0.5)))
ggsave("intermediate-products/presentation-plots/companies-202210/dist/companies_speed_dist.png", height = 4.5, width = 6.5, bg = "transparent")


### VESSEL SPEED
(pr <- pred.var(gam, tot, resp = "speed", var = "sqrt.vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = (sqrt.vess.speed.60^2)*2.23694, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel speed (mph, previous 1 minute)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(0.1, 1, 2, 4, 8, 12)) +
    scale_y_continuous(breaks = seq(1,10, by = 0.5)))
ggsave("intermediate-products/presentation-plots/companies-202210/speed/companies_speed_speed.png", height = 4.5, width = 6.5, bg = "transparent")


### VESSEL SPEED VARIATION
(pr <- pred.var(gam, tot, resp = "speed", var = "log.vess.speed.var.300", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.vess.speed.var.300), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Variation in boat speed (previous 5 minutes)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2)) +
    scale_y_continuous(breaks = seq(1,10, by = 0.5)))
ggsave("intermediate-products/presentation-plots/companies-202210/speed-var/companies_speed-var_speed.png", height = 4.5, width = 6.5, bg = "transparent")

### VESSEL DI
(pr <- pred.var(gam, tot, resp = "speed", var = "arcsin.vess.di.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sin(arcsin.vess.di.60), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Boat path straightness (previous 1 minute)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = asin.trans, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99,1)))
ggsave("intermediate-products/presentation-plots/companies-202210/di/companies_di_speed.png", height = 4.5, width = 6.5, bg = "transparent")



#---------------- SURFACE ACTIVITY: FOCAL --------------------

tot <- read.csv("intermediate-products/response-var-dfs/surface-active+feeding_focal.csv") %>%
  dplyr::select(-c("lat.whale", "lon.whale")) %>% drop_na() %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(dist<1000) %>% # more likely to see surfacings that are breaches at longer distances
  filter(vess.accel.min.300>-10) %>% # temporary filtering until boat stuff properly processes
  filter(vess.accel.max.300<0.5) %>%
  mutate(across(contains("speed.var"), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # log speed var
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"), # arcsin DI
         log.dist = log(dist))

gam <- readRDS("intermediate-products/gam/surface-active/foc/gam_surface-active_foc_final.rds")

var.fac = c("group", "seastate", "year") # factor
var.num = c("day.of.year", "log.dist", "vess.speed.60",  "log.vess.speed.var.60", "vess.accel.max.300")


### DISTANCE
(pr <- pred.var(gam, tot, resp = "surface.active", var = "log.dist", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.dist), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Distance to whale (meters)", y = "Rate of surface activity") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(10,20,50,100,200,500,1000)) +
    coord_cartesian(ylim = c(0,0.02)))
ggsave("intermediate-products/presentation-plots/companies-202210/dist/companies_surface-active_dist.png", height = 4.5, width = 6.5, bg = "transparent")


#---------------- SURFACE FEEDING: FOCAL --------------------

# focal vessel
tot <- read.csv("intermediate-products/response-var-dfs/surface-active+feeding_focal.csv") %>%
  dplyr::select(-c("lat.whale", "lon.whale")) %>% drop_na() %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(dist<1000) %>% # more likely to see surfacings that are breaches at longer distances
  filter(vess.accel.min.300>-10) %>% # temporary filtering until boat stuff properly processes
  filter(vess.accel.max.300<0.5) %>% # filtering out until boat stuff properly processes
  mutate(across(contains("speed.var"), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # log speed var
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"), # arcsin DI
         log.dist = log(dist))

gam <- readRDS("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_final.rds")

var.fac = c("seastate") 
var.num = c("day.of.year", "vess.speed.60",  "arcsin.vess.di.300")
var.rand = c("folnum.unique")

### VESSEL DI
(pr <- pred.var(gam, tot, resp = "surface.feeding", var = "arcsin.vess.di.300", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sin(arcsin.vess.di.300), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Boat path straightness (previous 1 minute)", y = "Rate of surface feeding") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = asin.trans, breaks = c(0, 0.4,0.6, 0.8, 0.9, 0.99)) + coord_cartesian(ylim = c(0,0.2)))
ggsave("intermediate-products/presentation-plots/companies-202210/di/companies_surface-active_di.png", height = 4.5, width = 6.5, bg = "transparent")


###########
### AIS ###
###########


#---------------- SURFACE ACTIVITY: AIS --------------------

# AIS
tot <- read.csv("intermediate-products/response-var-dfs/surface-active+feeding_ais.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique"), as.factor) %>%
  mutate(dummy = 1)  %>% # important dummy variable to factor out random effect when predicting response
  mutate(across(contains("meandist"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt mean dist
         across(contains("dist.max"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}")) # sqrt dist max

gam <- readRDS("intermediate-products/gam/surface-active/gam_surface-active_ais_final.rds")

var.fac = c("group")
var.num = c("oak.num.10.1500", "rib.num.30.1500", "sqrt.oak.meandist.10.1500", "sqrt.rib.meandist.30.1500")

### OAK NUMBER
(pr <- pred.var(gam, tot, resp = "surface.active", var = "oak.num.10.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = oak.num.10.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Number of oak boats (within 1500 m and 10 minutes)", y = "Rate of surface activity") + partial.theme + partial.grid.margin +
    coord_cartesian(ylim = c(0,0.02)))
ggsave("intermediate-products/presentation-plots/companies-202210/oak-num/companies_surface-active_oak-num.png", height = 4.5, width = 6.5, bg = "transparent")

### RIB NUMBER
(pr <- pred.var(gam, tot, resp = "surface.active", var = "rib.num.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = rib.num.30.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Number of RIBs (within 1500 m and 30 minutes)", y = "Rate of surface activity") + partial.theme + partial.grid.margin +
    coord_cartesian(ylim = c(0,0.3)))
ggsave("intermediate-products/presentation-plots/companies-202210/rib-num/companies_surface-active_rib-num.png", height = 4.5, width = 6.5, bg = "transparent")


#---------------- IBI: AIS --------------------

# AIS
tot <- read.csv("intermediate-products/response-var-dfs/surface-interval_ais.csv") %>%
  filter(!is.na(resp.rate) & !is.na(lon.end)) %>%
  mutate_at(vars(contains("dt")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate.last", "year", "folnum.unique", "surface.feeding.follow", "surface.active.follow"), as.factor) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
# filter(rib.num.30.1500<=4) # only 1/1243 rows had 5 ribs in 30.1500
  mutate(across(contains("meandist"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt mean dist
         across(contains("dist.max"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}")) # sqrt dist max

gam <- readRDS("intermediate-products/gam/ibi-mean/ais/gam_ibi-mean_ais_final.rds")

var.fac = c("surface.feeding.follow")
var.num = c("day.of.year", "sqrt.rib.meandist.10.1500", "oak.num.30.1500", "rib.num.10.1500")
var.rand = c("folnum.unique")

### OAK NUMBER
(pr <- pred.var(gam, tot, resp = "ibi.mean", var = "oak.num.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = oak.num.30.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Number of oak boats (within 1500 m and 30 minutes)", y = "Time between breaths (seconds)") + partial.theme + partial.grid.margin +
    scale_x_continuous(breaks = c(0:10)))
ggsave("intermediate-products/presentation-plots/companies-202210/oak-num/companies_ibi_oak-num.png", height = 4.5, width = 6.5, bg = "transparent")

### RIB NUMBER
(pr <- pred.var(gam, tot, resp = "ibi.mean", var = "rib.num.10.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = rib.num.10.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Number of RIBs (within 1500 m and 10 minutes)", y = "Time between breaths (seconds)") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,20,2)) +
    scale_x_continuous(breaks = c(0:10)))
ggsave("intermediate-products/presentation-plots/companies-202210/rib-num/companies_ibi_rib-num.png", height = 4.5, width = 6.5, bg = "transparent")


#---------------- DIVE TIME: AIS --------------------

tot <- read.csv("intermediate-products/response-var-dfs/dive-time_ais.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding.follow", "surface.active.follow"), as.factor) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(rib.num.30.1500<=4) # only 1/1243 rows had 5 ribs in 30.1500

gam <- readRDS("intermediate-products/gam/dive-time/ais/gam_dive-time_ais_final.rds")

var.fac = c("year", "surface.feeding.follow")
var.num = c("day.of.year", "rib.meandist.10.1500", "oak.num.30.1500", "rib.num.30.1500")
var.rand = c("folnum.unique")

### OAK NUMBER
(pr <- pred.var(gam, tot, resp = "dive.time", var = "oak.num.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = oak.num.30.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Number of oak boats (within 1500 m and 30 minutes)", y = "Dive time (seconds)") + partial.theme + partial.grid.margin +
    scale_x_continuous(breaks = c(0:10)))
ggsave("intermediate-products/presentation-plots/companies-202210/oak-num/companies_dive-time_oak-num.png", height = 4.5, width = 6.5, bg = "transparent")


### RIB NUMBER
(pr <- pred.var(gam, tot, resp = "dive.time", var = "rib.num.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = rib.num.30.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Number of RIBs (within 1500 m and 30 minutes)", y = "Dive time (seconds)") + partial.theme + partial.grid.margin +
    scale_x_continuous(breaks = c(0:10)))
ggsave("intermediate-products/presentation-plots/companies-202210/rib-num/companies_dive-time_rib-num.png", height = 4.5, width = 6.5, bg = "transparent")


##################
### COMPLIANCE ###
##################

# total!
tot <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv") %>%  
  filter(!is.na(focal.vess.speed.10)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>%
  dplyr::select(datetime, dist,focal.vess.speed.10) %>%
  mutate(focal.vess.speed.10 = focal.vess.speed.10*2.23694) %>% # m/s to mph
  mutate(compliant.300 = ifelse(focal.vess.speed.10>6,">6 mph", "<6 mph"),
         compliant.50 = ifelse(focal.vess.speed.10>3,">3 mph", "<3 mph"))
tot.join <- tot %>% dplyr::select(-dist, -focal.vess.speed.10)


#---------------- APPROACHING ZONE --------------------
# dive stuff
dive <- read.csv("intermediate-products/response-var-dfs/dive-time_ais.csv") %>% 
  filter(!is.na(focal.vess.speed.10)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>%  left_join(tot.join)

# speed
speed <- read.csv("intermediate-products/response-var-dfs/speed_ais.csv") %>%
  filter(!is.na(focal.vess.speed.10)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>%  left_join(tot.join)

# how about dive times?
df <- dive %>% filter(dist<= 300 & dist>50 & !is.na(compliant.300)) # 198 dives only, 38 = not compliant
(p.dive <- ggplot(data = df, aes(x = compliant.300, y = exp(log(dive.time)))) + geom_boxplot(width = 0.5)+
    scale_y_continuous(trans = "log", breaks = c(10,20,50,100,200,400,800)) +
    labs(x = "", y = "Dive time (seconds)") + plot.theme + theme(plot.margin =  unit(c(10,40,10,10), "pt")))

# speed
df <- speed %>% filter(dist<= 300 & dist>50 & !is.na(compliant.300)) # 931 rows, 293 not compliant
(p.speed <- ggplot(data = df, aes(x = compliant.300, y = speed)) + geom_boxplot(width = 0.5)+
    labs(x = "", y = "Swim speed (km/h)") + plot.theme + theme(plot.margin =  unit(c(10,10,10,40), "pt")))

# bring it all together
plot_grid(p.dive,  p.speed,align = "hv")
ggsave("intermediate-products/presentation-plots/companies-202210/compliance/companies_approaching-zeons_no-yes.png", dpi = 600, height = 4, width = 9.5, bg = "white")


#######################
### ENCOUNTER STAGE ###
#######################

#---------------- BEFORE/DURING/AFTER --------------------

# total
tot <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv") %>%
  mutate_at(.vars = c("surface.active", "surface.feeding"), funs(recode(., `1` = "yes", `0` = "no"))) %>%
  mutate_at(.vars = c("surface.active", "surface.feeding"), funs(as.factor)) %>%
  group_by(folnum.unique) %>% mutate(dive.time = lead(dive.time)) %>% ungroup()
# decided to keep dive time attached to second surfacing (allows us to look at after!)

surfacing <- read.csv("intermediate-products/response-var-dfs/surface-interval_ais.csv") %>%
  rename(datetime = dt.end) %>% dplyr::select(datetime, ibi.mean, breath.num)

speed <- read.csv("intermediate-products/response-var-dfs/speed_ais.csv") %>%
  dplyr::select(datetime, speed)

di <- read.csv("intermediate-products/response-var-dfs/di_ais.csv") %>%
  dplyr::select(datetime, DI)

tot <- tot %>%
  left_join(surfacing, by = "datetime") %>% left_join(speed, by = "datetime") %>% left_join(di, by = "datetime")

# before/during/after
tot.stage <- tot %>% ungroup() %>%
  filter(!is.na(encounter.stage)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>% # going for 10 minutes, 1500 m
  mutate(encounter.stage = fct_relevel(encounter.stage, "before", "during", "after"))

# approach/encounter/depart
tot.enc <- tot %>% ungroup() %>%
  filter(!is.na(vess.encounter.stage)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>% # going for 10 minutes, 1500 m
  mutate(vess.encounter.stage = recode(vess.encounter.stage, `-1` = "approach", `0` = "encounter", `1` = "departure"),
         vess.encounter.stage = fct_relevel(vess.encounter.stage, "approach", "encounter", "departure"))

#---------------- BEFORE/DURING/AFTER --------------------

# speed
df = tot.stage %>% filter(!is.na(encounter.stage) & !is.na(speed))# 1266
(p.speed <- ggplot(data = df, aes(x = encounter.stage, y = speed)) +
    geom_boxplot(size = 0.8, width = 0.4, fill = "transparent") + labs(x = "", y = "Swim speed (km/h)") +
    scale_y_continuous(trans = "sqrt", breaks = c(0.5, 2, 5, 10, 15, 20)) + plot.theme)
ggsave("intermediate-products/presentation-plots/companies-202210/encounter/companies_enc-stage_speed.png", dpi = 600, height = 3.5, width = 5.5, bg = "white")

#---------------- APPROACH/ENCOUNTER --------------------

# speed
df = tot.enc %>% filter(!is.na(speed) & !is.na(vess.encounter.stage))# 488
df <- df %>% filter(vess.encounter.stage != "departure")
(p.speed <- ggplot(data = df, aes(x = vess.encounter.stage, y = speed)) +
    geom_boxplot(size = 0.8, width = 0.4, fill = "transparent") + labs(x = "", y = "Swim speed (km/h)") +
    scale_y_continuous(trans = "sqrt", breaks = c(0.5, 2, 5, 10, 15, 20)) + plot.theme)
ggsave("intermediate-products/presentation-plots/companies-202210/encounter/companies_vess-enc-stage_speed.png", dpi = 600, height = 3.5, width = 4.5, bg = "white")


########################
### ENCOUNTER MINUTE ###
########################

#---------------- DATA --------------------

tot <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv") %>%
  mutate(trip.id = paste0(as.Date(datetime),"_",vessel)) %>%
  group_by(trip.id, id.focal.catalogue) %>%
  # we only want the first follow for a trip for each animal
  filter(folnum.unique == first(folnum.unique)) %>%
  filter(!is.na(encounter.minute)) %>%
  group_by(folnum.unique) %>%
  filter(first(boat.ids.30.1500) == vessel | first(boat.ids.30.1500) == "") %>%
  mutate(sqrt.encounter.minute = sqrt(encounter.minute),
         folnum.unique = as.factor(folnum.unique)) %>%
  mutate(dummy = 1)

dive <- read.csv("intermediate-products/response-var-dfs/surface-interval_ais.csv") %>%
  rename(datetime = dt.end) %>%
  dplyr::select(datetime, ibi.mean, breath.num)

speed <- read.csv("intermediate-products/response-var-dfs/speed_ais.csv") %>%
  dplyr::select(datetime, speed, move.id.tot)

di <- read.csv("intermediate-products/response-var-dfs/di_ais.csv") %>%
  dplyr::select(datetime, DI)

tot <- tot %>% 
  left_join(dive) %>%
  left_join(speed) %>% 
  left_join(di)


#---------------- SPEED --------------------

df <- tot
mod.speed  <- gam(speed ~ s(sqrt.encounter.minute, bs = "cs", k = 10) +
                    s(folnum.unique, bs = "re", by = dummy), data = df, family = gaussian(link = "sqrt"))
(pr <- plot.folmin(mod.speed, df, resp = "speed", var = "sqrt.encounter.minute", var.rand = "folnum.unique") %>%
    ggplot(aes(x = sqrt.encounter.minute^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Encounter minute", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(0.1,1,2,5,10,20,30)))
ggsave("intermediate-products/presentation-plots/companies-202210/enc-minute/companies_speed_enc-minute.png", height = 4.5, width = 6.5, bg = "transparent")


#---------------- DIVE TIME --------------------

df <- tot %>% filter(encounter.minute>=10)
mod.divetime <- gam(dive.time ~ s(sqrt.encounter.minute, bs = "cs", k = 10) +
                              s(folnum.unique, bs = "re", by = dummy), data = df, family = gaussian(link = "log"))
(pr <- plot.folmin(mod.divetime, df, resp = "dive.time", var = "sqrt.encounter.minute", var.rand = "folnum.unique") %>%
    ggplot(aes(x = sqrt.encounter.minute^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Encounter minute", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(0.1,1,2,5,10,20,30)))
ggsave("intermediate-products/presentation-plots/companies-202210/enc-minute/companies_dive-time_enc-minute.png", height = 4.5, width = 6.5, bg = "transparent")


#####################
### EXAMPLE PLOTS ###
#####################

# creating eg data frame
df.eg <- data.frame(boat = 1:10, whale = 1:10) %>%
  mutate(lcl.small = whale-0.4, ucl.small = whale+0.4,
         lcl.big = whale-8, ucl.big = whale+8)

# small CI
ggplot(data = df.eg, aes(x = boat)) +
  geom_line(aes(y = whale), size = 1) +
  geom_ribbon(aes(ymin = lcl.small, ymax = ucl.small), 
              linetype = "dashed", color = "black", fill = "transparent") +
  partial.theme + partial.grid.margin +
  labs(x = "Boat behaviour", y = "Whale behaviour")
ggsave("intermediate-products/presentation-plots/companies-202210/CI-eg_small.png", height = 4.5, width = 6.5, bg = "transparent")


# big CI
ggplot(data = df.eg, aes(x = boat)) +
  geom_line(aes(y = whale), size = 1) +
  geom_ribbon(aes(ymin = lcl.big, ymax = ucl.big), 
              linetype = "dashed", color = "black", fill = "transparent") +
  partial.theme + partial.grid.margin +
  labs(x = "Boat behaviour", y = "Whale behaviour")
ggsave("intermediate-products/presentation-plots/companies-202210/CI-eg_big.png", height = 4.5, width = 6.5, bg = "transparent")
