######################
### V2: GAM: DI, FOCAL ###
######################
### 12/01/2023
### Tom Grove
### tomgrove20@yahoo.co.uk

# 1. transformations
# 2. Collinearity
# 4. Final GAM (GCV + backward selection)
# 5. Results + plotting


### PACKAGES
packages <- c("tidyverse", "ppcor","RColorBrewer", "scales","survMisc", "Metrics", "corrplot", "lme4", "MASS", "mgcv", "tidymv", "mgcViz", "gridExtra", "gratia", "ggcorrplot", "cowplot", "ggpubr", "ungeviz")
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

# foc
tot <- read.csv("intermediate-products/response-var-dfs/di_focal.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding", "surface.active"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  # we also need to convert m/s to mph!
  mutate(across(contains("vess.speed"), ~.*conv), # speed and sd speed
         across(contains("accel"), ~.*conv.sq)) %>% # acceleration
  filter(vess.accel.max.60<20000, vess.speed.var.60 < 9)

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 598 surfacing pairs, 100 follows

# now we need to add encounter minute to this (rough!!)
follows.encmin <- read.csv("intermediate-products/follows/follows_start+end_consecutive_rough.csv") %>%
  mutate(folnum.unique = as.factor(folnum.unique))

tot <- tot %>% left_join(follows.encmin) %>%
  group_by(enc.num) %>%
  mutate(encounter.minute = as.numeric(difftime(datetime, start, unit = "mins"))) %>% ungroup()


#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = DI)) + geom_histogram()
ggplot(data = tot, aes(x = asin(DI))) + geom_histogram() # going to arcsin transform this simply to reduce skew


# define explanatory variables, then view distribution
vars <- c("dist", "seastate", "surface.feeding", "surface.active", "group", "difftime.prev", "difftime.next", "day.of.year", "vess.di.60", "vess.accel.60", "vess.speed.60", "vess.speed.var.60", "vess.accel.max.300")
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
ggplot(data = tot, aes(x = log(difftime.prev))) + geom_histogram()
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
ggplot(data = tot, aes(x = cube_root(vess.accel.max.60))) + geom_histogram()


# creating new transformed columns
tot <- tot %>%
  mutate(arcsin.DI = asin(DI),
         across(contains(c("difftime")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(contains(c("vess.speed", "dist")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains(c("accel", "encounter.minute")), .fns = list(cubrt = ~sign(.)*(abs(.)^(1/3))), .names = "{fn}.{col}"),
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group","surface.feeding", "surface.active", "log.difftime.prev", "log.difftime.next", "year", "day.of.year", "sqrt.dist", "arcsin.vess.di.60",  "sqrt.vess.speed.60", "sqrt.vess.speed.var.60", "cubrt.encounter.minute")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_collinear.png", scale = 1.2)

# as expected, correlated accel and accel.max but only 0.6. Keeping in for now!


#---------------- FULL MODEL + AUTOCORRELATION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") # factor
var.num = c("day.of.year", "sqrt.dist", "log.difftime.prev", "log.difftime.next", "sqrt.vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "cubrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "identity")))
summary(gam)

# View ACF plots
acf.pacf.plot(acf.structured(gam, 10, "gam"), 10)
ggsave("intermediate-products/gam-v3/di/foc/gam_di_foc_acf.png", width = 11, height = 3.6)

# trying with a bam
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
# for bam we need to log-transform
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam1 <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                            family = gaussian(link = "identity"), rho = 0.1, start = ar.start))

acf.pacf.plot(acf.structured(gam1, 10, "bam"),10)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") # factor
var.num = c("day.of.year", "sqrt.dist", "log.difftime.prev", "log.difftime.next", "sqrt.vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "cubrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = 0.1, start = ar.start))
summary(gam); AIC(gam) # dev = 14.2%, -REML = 113.5, AIC = 222.0 

# final model (manual selection)
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") 
var.num = c("log.difftime.prev", "sqrt.vess.speed.60",  "arcsin.vess.di.60", "cubrt.encounter.minute")
var.rand = c("folnum.unique")
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = 0.1, start = ar.start))
summary(gam); AIC(gam) # dev = 14.3%, -REML = 113.5, AIC = 222

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam-v3/di/foc/gam_di_foc_final_v3.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam-v3/di/foc/gam_di_foc_final_v3.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/di/foc/gam_di_foc_results_v3.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") 
var.num = c("log.difftime.prev", "sqrt.vess.speed.60",  "arcsin.vess.di.60", "cubrt.encounter.minute")
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)

# difftime.prev
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(previous IBI (sec))", y = "s(log(previous IBI))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "log.difftime.prev", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.difftime.prev), y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Previous IBI (sec)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.8, 0.9, 0.95, 0.98, 0.99, 1)) +
    scale_x_continuous(trans = "log", breaks = c(5, 10, 20, 50, 100, 200, 400, 800)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/foc/gam_di_foc_difftime.prev_v3.png", height = 4.5, width = 13)

# vess.speed.60
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(vessel speed (mph))", y = "s(sqrt(vessel speed))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "sqrt.vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.vess.speed.60^2, y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel speed (mph)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.9, 0.95, 0.98, 0.99, 1)) +
    scale_x_continuous(trans = "sqrt", breaks = c(0.1, 0.5, 1,2,4,6,8,10)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/foc/gam_di_foc_vess.speed.60_v3.png", height = 4.5, width = 13)

# vess.di.60
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "arcsin(vessel DI)", y = "s(arcsin(vessel DI))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "arcsin.vess.di.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sin(arcsin.vess.di.60), y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel DI", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.7, 0.8, 0.85, 0.9, 0.95, 0.99, 1)) +
    scale_x_continuous(trans = asin.trans, breaks = c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/foc/gam_di_foc_vess.di.60_v3.png", height = 4.5, width = 13)


# cubrt.encounter.minute
(sm <- plot(sm(b, 4)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "cubrt(encounter minute)", y = "s(cubrt(encounter minute))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "cubrt.encounter.minute", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = cubrt.encounter.minute^(3), y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Encounter minute", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 1)) +
    scale_x_continuous(trans = cubrt.trans, breaks = c(0,1,5,10,20,30,50)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/foc/gam_di_foc_encounter.minute_v3.png", height = 4.5, width = 13)


# group
v <- "group"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.5, 0.7,0.8, 0.9, 0.95, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/foc/gam_di_foc_group_v3.png", height = 4.5, width = 13)

# sea state
v <- "seastate"
plotdata <- plot(b, allTerms = T, select = 7)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "seastate", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = seastate, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.5, 0.7,0.8, 0.9, 0.95, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/foc/gam_di_foc_seastate_v3.png", height = 4.5, width = 13)

# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 8)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Directness index") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/foc/gam_di_foc_year_v3.png", height = 4.5, width = 13)

# surface.feeding
v <- "surface.feeding"
plotdata <- plot(b, allTerms = T, select = 9)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +  
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface feeding", y = "Partial effect of surface feeding"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "surface.feeding", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.feeding, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Surface feeding", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")) +
    scale_y_continuous(trans = asin.trans, breaks = c(0.5, 0.7,0.8, 0.9, 0.95, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/foc/gam_di_foc_surface.feeding_v3.png", height = 4.5, width = 13)

# surface.feeding
v <- "surface.active"
plotdata <- plot(b, allTerms = T, select = 10)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +  
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface activity", y = "Partial effect of surface activity"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "surface.active", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.active, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Surface active", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")) +
    scale_y_continuous(trans = asin.trans, breaks = c(0.5, 0.7,0.8, 0.9, 0.95, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/foc/gam_di_foc_surface.active_v3.png", height = 4.5, width = 13)
