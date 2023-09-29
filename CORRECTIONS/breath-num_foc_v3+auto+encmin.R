### V2: GAM: SURFACING BREATHS, FOCAL 

### 11/03/2023
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


#---------------- DATA --------------------

# conversion rate from m/s to mph
conv <- 2.23694 # for m/s to mph
conv.sq <- 8052.9692102775 # for m/s2 to mph2

# focal vessel
tot <- read.csv("intermediate-products/response-var-dfs/surface-interval_focal.csv") %>%
  filter(!is.na(dist.mean)) %>%
  mutate_at(vars(contains("dt")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate.last", "year", "folnum.unique", "surface.active.follow", "surface.feeding.follow"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1)  %>% # important dummy variable to factor out random effect when predicting response
  dplyr::select(-contains(".10")) %>%
  dplyr::select(-contains(c("0_end", "0_start"))) %>%
  rename_at(vars(contains("_mean")), funs(str_replace(.,"_mean",""))) %>%
  filter(vess.accel.max.60^(1/3)>-2) %>%
  # we also need to convert m/s to mph!
  mutate(across(contains("speed"), ~.*conv), # speed and sd speed
         across(contains("accel"), ~.*conv.sq)) # acceleration
# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 142 surfacing intervals, 56 follows

# now we need to add encounter minute to this (rough!!)
follows.encmin <- read.csv("intermediate-products/follows/follows_start+end_consecutive_rough.csv") %>%
  mutate(folnum.unique = as.factor(folnum.unique))

tot <- tot %>% left_join(follows.encmin) %>%
  group_by(enc.num) %>%
  mutate(encounter.minute = as.numeric(difftime(dt.start, start, unit = "mins"))) %>% ungroup()

#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = breath.num)) + geom_histogram()
ggplot(data = tot, aes(x = log(breath.num))) + geom_histogram() # so Poisson with log link

# define explanatory variables, then view distribution
vars <- c("dist.mean", "seastate.last", "dive.time.before", "group", "day.of.year", "vess.di.60", "vess.accel.60", "vess.speed.60", "vess.speed.var.300", "vess.accel.max.60")
for (i in vars) {
  if(class(tot[,i]) == "numeric"){ # if numeric, provide a histogram
    print(ggplot(data = tot, aes(x = get(i))) + geom_histogram() + xlab(i))
  } else { # if not numeric, provide a bar plot
    print(ggplot(data = tot, aes(x = get(i))) + geom_bar() + xlab(i))
  }
  readline(prompt = "next plot")
}

# dist.mean (log)
ggplot(data = tot, aes(x = log(dist.mean))) + geom_histogram()
# dive time before (log)
ggplot(data = tot, aes(x = log(dive.time.before))) + geom_histogram()
# DI (arcsin)
ggplot(data = tot, aes(x = asin(vess.di.60))) + geom_histogram()
# speed var (sqrt)
ggplot(data = tot, aes(x = sqrt(vess.speed.var.60))) + geom_histogram()
# acceleration
ggplot(data = tot, aes(x = (vess.accel.max.60)^(1/3))) + geom_histogram()
# cannot transform accel.max (either side of 0!)

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("dive.time","dist")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(contains("speed.var"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains("vess.di."), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"),
         across(contains(c("accel.max", "encounter.minute")), .fns = list(cubrt = ~sign(.)*(abs(.)^(1/3))), .names = "{fn}.{col}"))


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate.last", "group", "surface.feeding.follow", "surface.active.follow", "year", "day.of.year", "log.dive.time.before", "log.dist.mean", "arcsin.vess.di.60", "vess.speed.60", "sqrt.vess.speed.var.60", "cubrt.encounter.minute")

# correlation matrix
x <- tot[,vars] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_collinear.png", scale = 1.2)

# actually no troubling correlations at all! Let's keep them in!


#---------------- FULL MODEL + AUTOCORRELATION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.num = c("day.of.year","log.dive.time.before",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "cubrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = quasipoisson(link = "log")))
summary(gam) # dev = 91.1%, gcv = 0.347

# View ACF plots
acf.pacf.plot(acf.structured(gam, 5, "gam"), 5)
ggsave("intermediate-products/gam-v3/breaths-surfacing/foc/gam_breaths-surfacing_foc_acf.png", width = 11, height = 3.6)

# trying with a bam
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
# for bam we need to log-transform
tot <- tot %>% mutate(log.breath.num = log(breath.num))
f <- as.formula(paste("log.breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam1 <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                            family = gaussian(link = "identity"), rho = -0.3, start = ar.start))

acf.pacf.plot(acf.structured(gam1, 5, "bam"),5)

#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.num = c("day.of.year","log.dive.time.before",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "cubrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("log.breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                            family = gaussian(link = "identity"), rho = -0.3, start = ar.start))
summary(gam); AIC(gam) # dev = 80.2%, -REML = 95.4, AIC = 147.2


# remove log.dist.mean, vess.speed.60
var.fac = c("group", "year", "surface.active.follow")
var.num = c("log.dive.time.before", "sqrt.vess.speed.var.60", "cubrt.encounter.minute")
var.rand = c("folnum.unique")
f <- as.formula(paste("log.breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = -0.3, start = ar.start))
summary(gam); AIC(gam) # dev = 80.3%, -REML = 93.2, AIC = 143.6

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam
saveRDS(gam.final, "intermediate-products/gam-v3/breaths-surfacing/foc/gam_breaths-surfacing_foc_final_v3.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam-v3/breaths-surfacing/foc/gam_breaths-surfacing_foc_final_v3.rds")

# defining variable combinations
var.fac = c("group", "year", "surface.active.follow")
var.num = c("log.dive.time.before", "sqrt.vess.speed.var.60", "cubrt.encounter.minute")
var.rand = c("folnum.unique")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/breaths-surfacing/foc/gam_breaths-surfacing_foc_results_v3.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
print(plot(b, allTerms = T), pages = 1)

# log.dive.time.before
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(previous time dive (sec))", y = "s(log(previous dive time))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "log.dive.time.before", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.dive.time.before), y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Previous dive time (sec)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(50,100,200,400,800)) + 
    scale_y_continuous(breaks = seq(0,20,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/foc/gam_breaths-surfacing_foc_dive.time.before_v3.png", height = 4.5, width = 13)


# sqrt.vess.speed.var.60
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(mean vessel speed SD (mph))", y = "s(sqrt(mean vessel speed SD))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "sqrt.vess.speed.var.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.vess.speed.var.60^2, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Mean vessel speed SD (mph)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt"))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/foc/gam_breaths-surfacing_foc_vess.speed.var.60_v3.png", height = 4.5, width = 13)


# cubrt.encounter.minute
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Encounter minute", y = "s(encounter minute)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "cubrt.encounter.minute", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = cubrt.encounter.minute^3, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Encounter minute", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = cubrt.trans, breaks = c(0,1,5,10,20,30)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/foc/gam_breaths-surfacing_foc_encounter.minute_v3.png", height = 4.5, width = 13)


# group
v <- "group"
plotdata <- plot(b, allTerms = T, select = 5)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "breath.num", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/foc/gam_breaths-surfacing_foc_group_v3.png", height = 4.5, width = 13)


# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "breath.num", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/foc/gam_breaths-surfacing_foc_year_v3.png", height = 4.5, width = 13)


# surface.active.follow
v <- "surface.active.follow"
plotdata <- plot(b, allTerms = T, select = 7)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface activity in follow", y = "Partial effect of surface activity"))
(pr <- pred.var.fac(gam, tot, resp = "ibi.mean", var = "surface.active.follow", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.active.follow, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Surface activity in follow", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/foc/gam_breaths-surfacing_foc_surface.active.follow_v3.png", height = 4.5, width = 13)



