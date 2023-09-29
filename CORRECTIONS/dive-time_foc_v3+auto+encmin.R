### V3: GAM: DIVE TIME, FOCAL

### 25/05/2023
### Tom Grove
### tomgrove20@yahoo.co.uk

# 1. transformations
# 2. Collinearity
# 4. Final GAM (GCV + backward selection)
# 5. Results + plotting


### PACKAGES
packages <- c("tidyverse", "ppcor","RColorBrewer", "scales","survMisc", "Metrics", "corrplot", "lme4", "MASS", "mgcv", "tidymv", "mgcViz", "gridExtra", "gratia", "ggcorrplot", "cowplot", "ggpubr", "gamm4", "forecast")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)


#---------------- FUNCTIONS + THEME --------------------

source("code/functions.R")
source("code/themes.R")


#---------------- DATA --------------------

# conversion rate from m/s to mph
conv <- 2.23694 # for m/s to mph
conv.sq <- 8052.9692102775 # for m/s2 to mph2

# focal vessel
tot <- read.csv("intermediate-products/response-var-dfs/dive-time_focal.csv") %>%
  dplyr::select(-c("lat.whale", "lon.whale")) %>% drop_na() %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding.follow", "surface.active.follow"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(vess.accel.max.300<0.5 & log(vess.speed.var.60)>-3) %>% # filtering out until boat stuff properly processes
  # we also need to convert m/s to mph!
  mutate(across(contains("speed"), ~.*conv), # speed and sd speed
         across(contains("accel"), ~.*conv.sq)) # acceleration

# quick look at the data set
colSums(is.na(tot)) 

# now we need to add encounter minute to this (rough!!)
follows.encmin <- read.csv("intermediate-products/follows/follows_start+end_consecutive_rough.csv") %>%
  mutate(folnum.unique = as.factor(folnum.unique))

tot <- tot %>% left_join(follows.encmin)  %>%
  group_by(enc.num) %>%
  mutate(encounter.minute = as.numeric(difftime(datetime, start, unit = "mins"))) %>% ungroup()


#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = dive.time)) + geom_histogram()
ggplot(data = tot, aes(x = log(dive.time))) + geom_histogram()
# based on later results (wonky QQ). Log-transforming to reduce skew!

# define explanatory variables, then view distribution
vars <- c("dist", "seastate", "group", "surface.feeding.follow", "surface.active.follow", "day.of.year", "vess.di.60", "vess.accel.60", "vess.speed.60", "vess.speed.var.60", "vess.accel.max.300")
for (i in vars) {
  if(class(tot[,i]) == "numeric"){ # if numeric, provide a histogram
    print(ggplot(data = tot, aes(x = get(i))) + geom_histogram() + xlab(i))
  } else { # if not numeric, provide a bar plot
    print(ggplot(data = tot, aes(x = get(i))) + geom_bar() + xlab(i))
  }
  readline(prompt = "next plot")
}

# accel  max not ideal but difficult to transform since crossing 0. DI is left-skewed and speed var is right skewed

# speed.var: log transformation
sort(tot$vess.speed.var.60); ggplot(data = tot, aes(x = log(vess.speed.var.60))) + geom_histogram()
# di: arcsin transformation
ggplot(data = tot, aes(x = asin(vess.di.60))) + geom_histogram()
ggplot(data = tot, aes(x = vess.di.60)) + geom_histogram()
# dist: log transformation
ggplot(data = tot, aes(x = log(dist))) + geom_histogram()
ggplot(data = tot, aes(x = (vess.accel.max.60)^(1/3))) + geom_histogram()

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains("speed.var"), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # log speed var
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"), # arcsin DI
         log.dist = log(dist),
         across(contains(c("encounter.minute")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains("accel.max"), .fns = list(cubrt = ~sign(.)*(abs(.)^(1/3))), .names = "{fn}.{col}"))


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group", "surface.feeding.follow", "surface.active.follow", "year", "day.of.year", "log.dist", "arcsin.vess.di.60",  "vess.speed.60", "log.vess.speed.var.60", "sqrt.encounter.minute")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_collinear.png", scale = 1.2)
# actually no troubling correlations at all! Let's keep them in!


#---------------- FULL MODEL + AUTOCORRELATION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.num = c("day.of.year", "log.dist", "vess.speed.60",  "log.vess.speed.var.60", "arcsin.vess.di.60", "sqrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("dive.time ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "log"), select = TRUE))
summary(gam) # dev = 91.6%, gcv = 3484

# View ACF plots
acf.pacf.plot(acf.structured(gam, 4, "gam"), 4)
ggsave("intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_acf.png", width = 11, height = 3.6)

# trying with a bam
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
# for bam we need to log-transform
tot <- tot %>% mutate(log.dive.time = log(dive.time))
f <- as.formula(paste("log.dive.time ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam1 <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "identity"),
                            rho = -0.4, start = ar.start))

acf.pacf.plot(acf.structured(gam1, 4, "bam"),4)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.num = c("day.of.year", "log.dist", "vess.speed.60",  "log.vess.speed.var.60", "arcsin.vess.di.60", "sqrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random

# formula
f <- as.formula(paste("log.dive.time ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# gam!
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = -0.4, start = ar.start))
summary(gam); AIC(gam) # dev = 78.3%, -REML = 246.9, AIC = 390.7


# final model (manual selection)
var.fac = c("seastate", "year", "surface.feeding.follow")
var.num = c("day.of.year", "log.dist", "vess.speed.60",  "sqrt.encounter.minute")
f <- as.formula(paste("log.dive.time ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = -0.4, start = ar.start))
summary(gam); AIC(gam) # dev = 78.6%, -REML = 247, AIC = 382.5


# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_final_v3.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_final_v3.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_results_v3.csv", row.names = FALSE)

var.fac = c("seastate", "year", "surface.feeding.follow")
var.num = c("day.of.year", "log.dist", "vess.speed.60",  "sqrt.encounter.minute")
var.rand = c("folnum.unique") # random 

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
print(plot(b, allTerms = T), pages = 1)

# day.of.year
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "log.dive.time", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    coord_cartesian(ylim = c(0,2000)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_julian_v3.png", height = 4.5, width = 13)


# log.dist
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(distance (m))", y = "s(log(distance))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "log.dive.time", var = "log.dist", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.dist), y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Distance (m)", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(20,50,100,200,400,1000)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_dist_v3.png", height = 4.5, width = 13)


# vess.speed.60
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Vessel speed (mph)", y = "s(vessel speed)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "log.dive.time", var = "vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = vess.speed.60, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Speed (mph)", y = "Dive time (sec)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_vess.speed.60_v3.png", height = 4.5, width = 13)

# sqrt.encounter.minute
(sm <- plot(sm(b, 4)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(encounter minute)", y = "s(sqrt(encounter minute))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "dive.time", var = "sqrt.encounter.minute", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.encounter.minute^2, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Encounter minute", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(0,1,5,10,20,30)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_encounter.minute_v3.png", height = 4.5, width = 13)

# seastate
v <- "seastate"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "log.dive.time", var = "seastate", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = seastate, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    coord_cartesian(ylim = c(30,370)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_seastate_v3.png", height = 4.5, width = 13)

# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 7)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "log.dive.time", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    coord_cartesian(ylim = c(30,370)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_year_v3.png", height = 4.5, width = 13)

# surface.feeding.follow
v <- "surface.feeding.follow"
plotdata <- plot(b, allTerms = T, select = 8)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface feeding in follow", y = "Partial effect of surface feeding"))
(pr <- pred.var.fac(gam, tot, resp = "log.dive.time", var = "surface.feeding.follow", var.fac, var.num, var.rand) %>%
    ggplot(aes(x = surface.feeding.follow, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Surface feeding in follow", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")) +
    coord_cartesian(ylim = c(30,370)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/foc/gam_dive-time_foc_surface.feeding.follow_v3.png", height = 4.5, width = 13)
