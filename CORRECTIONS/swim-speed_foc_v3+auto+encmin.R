
### V3: GAM: SWIM SPEED, FOCAL 
### 24/05/2023
### Tom Grove
### tomgrove20@yahoo.co.uk

# 1. transformations
# 2. Collinearity
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


#---------------- DATA --------------------

# conversion rate from m/s to mph
conv <- 2.23694 # for m/s to mph
conv.sq <- 8052.9692102775 # for m/s2 to mph2

# foc
tot <- read.csv("intermediate-products/response-var-dfs/speed_focal.csv") %>%
  mutate_at(vars(contains("datetime")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.active", "surface.feeding"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  filter_at(vars(contains("vess")), isnt.na) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  # we also need to convert m/s to mph!
  mutate(across(contains("vess.speed"), ~.*conv), # speed and sd speed
         across(contains("accel"), ~.*conv.sq)) %>% # acceleration
  filter(vess.accel.max.60<20000, vess.speed.var.60 < 9)

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 737 surfacing pairs, 126 follows

# now we need to add encounter minute to this (rough!!)
follows.encmin <- read.csv("intermediate-products/follows/follows_start+end_consecutive_rough.csv") %>%
  mutate(folnum.unique = as.factor(folnum.unique))

tot <- tot %>% left_join(follows.encmin) %>%
  group_by(enc.num) %>%
  mutate(encounter.minute = as.numeric(difftime(datetime, start, unit = "mins"))) %>% ungroup()


#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = speed)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(speed))) + geom_histogram() # so Gaussian with log link?

# define explanatory variables, then view distribution
vars <- c("dist", "seastate", "surface.feeding", "surface.active", "group", "diff.time", "day.of.year", "vess.di.60", "vess.accel.60", "vess.speed.60", "vess.speed.var.60", "vess.accel.max.300")
for (i in vars) {
  if(class(tot[,i]) == "numeric"){ # if numeric, provide a histogram
    print(ggplot(data = tot, aes(x = get(i))) + geom_histogram() + xlab(i))
  } else { # if not numeric, provide a bar plot
    print(ggplot(data = tot, aes(x = get(i))) + geom_bar() + xlab(i))
  }
  readline(prompt = "next plot")
}

# dist
ggplot(data = tot, aes(x = sqrt(dist))) + geom_histogram()
# diff.time 
ggplot(data = tot, aes(x = log(diff.time))) + geom_histogram()
# di
ggplot(data = tot, aes(x = vess.di.60)) + geom_histogram()
ggplot(data = tot, aes(x = asin(vess.di.60))) + geom_histogram()
# vess.accel
ggplot(data = tot, aes(x = vess.accel.60)) + geom_histogram() # can't really do anything
# vess.speed
ggplot(data = tot, aes(x = sqrt(vess.speed.60))) + geom_histogram()
# vess.speed.var
ggplot(data = tot, aes(x = sqrt(vess.speed.var.60))) + geom_histogram()
# vess.accel.max
ggplot(data = tot, aes(x = (vess.accel.max.60)^(1/3))) + geom_histogram()

# speed.var: log transformation
sort(tot$vess.speed.var.60); ggplot(data = tot, aes(x = log(vess.speed.var.60))) + geom_histogram()
# di: arcsin transformation
ggplot(data = tot, aes(x = asin(vess.di.60))) + geom_histogram()
# dist: log transformation
ggplot(data = tot, aes(x = log(dist))) + geom_histogram()

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("diff.time")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(contains(c("vess.speed", "encounter.minute")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains(c("accel", "dist")), .fns = list(cubrt = ~sign(.)*(abs(.)^(1/3))), .names = "{fn}.{col}"),
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group","surface.feeding", "surface.active", "diff.time", "year", "day.of.year", "cubrt.dist", "arcsin.vess.di.60",  "sqrt.vess.speed.60", "sqrt.vess.speed.var.60","sqrt.encounter.minute")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))


#---------------- FULL MODEL + AUTOCORRELATION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") # factor
var.num = c("day.of.year", "cubrt.dist", "log.diff.time", "sqrt.vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "sqrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random

# formula
f <- as.formula(paste("speed ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "log")))
summary(gam) # dev = 42.3%, GCV = 4.33

# View ACF plots
acf.pacf.plot(acf.structured(gam, 10, "gam"), 10)
ggsave("intermediate-products/gam-v3/swim-speed/foc/gam_swim-speed_foc_acf.png", width = 11, height = 3.6)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") # factor
var.num = c("day.of.year", "cubrt.dist", "log.diff.time", "sqrt.vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "sqrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("speed ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "log")))
summary(gam); AIC(gam) # dev = 39.7%, -REML = 1623.6, AIC = 3184.3

# final model (manual selection)
var.fac = c("surface.feeding") 
var.num = c("log.diff.time", "sqrt.vess.speed.60", "arcsin.vess.di.60", "sqrt.encounter.minute")
var.rand = c("folnum.unique")
f <- as.formula(paste("speed ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "log")))
summary(gam); AIC(gam) # dev = 39.5%, -REML = 1618.4, AIC = 3181.7

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam-v3/swim-speed/foc/gam_swim-speed_foc_final_v3.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam-v3/swim-speed/foc/gam_swim-speed_foc_final_v3.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/swim-speed/foc/gam_swim-speed_foc_results_v3.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
var.fac = c("surface.feeding", "surface.active")
var.num = c("day.of.year", "cubrt.dist", "log.diff.time", "sqrt.vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "sqrt.encounter.minute") # numeric 
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)

# log.diff.time
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(inter-breath interval (sec))", y = "s(log(inter-breath interval))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "log.diff.time", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.diff.time), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Inter-breath interval (sec)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(5,10,20,50,100,200,400,800)) +
    scale_y_continuous(breaks = seq(1,10, by = 1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/swim-speed/foc/gam_swim-speed_foc_ibi_v3.png", height = 4.5, width = 13)


# vess.speed
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(vessel speed (mph))", y = "s(sqrt(vessel speed))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "sqrt.vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.vess.speed.60^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel speed (mph)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(0.5, 1, 2, 5, 10)) +
    scale_y_continuous(breaks = seq(1,10, by = 0.5)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/swim-speed/foc/gam_swim-speed_foc_vess.speed.60_v3.png", height = 4.5, width = 13)

# DI
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "arcsin(vessel DI)", y = "s(arcsin(vessel DI))") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "arcsin.vess.di.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sin(arcsin.vess.di.60), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel DI", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = asin.trans, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99,1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/swim-speed/foc/gam_swim-speed_foc_vess.di.60_v3.png", height = 4.5, width = 13)

# encounter.minute
(sm <- plot(sm(b, 4)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = bquote("sqrt(encounter minute)"), y = "s(sqrt(encounter minute))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "sqrt.encounter.minute", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.encounter.minute^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Encounter minute", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin + scale_y_continuous(breaks = seq(0,16, by = 2)) +
    scale_x_continuous(trans = "sqrt", breaks = c(0, 1, 2, 5, 10, 20, 50)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/swim-speed/foc/gam_swim-speed_foc_enc.min_v3.png", height = 4.5, width = 13)

# surface.feeding
v <- "surface.feeding"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +  
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface feeding", y = "Partial effect of surface feeding"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "surface.feeding", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.feeding, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Surface feeding", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")) + scale_y_continuous(breaks = seq(0,20,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/swim-speed/foc/gam_swim-speed_foc_surface.feeding_v3.png", height = 4.5, width = 13)
