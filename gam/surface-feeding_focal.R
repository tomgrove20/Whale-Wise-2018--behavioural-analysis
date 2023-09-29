####################################
### GAM: SURFACE ACTIVITY, FOCAL ###
####################################
### 16/02/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

# 1. transformations
# 2. Collinearity
# 3. time frame selection
# 4. Final GAM (GCV + backward selection)
# 5. Results + plotting

# note: PQL is used in gamm or gamm4, which performs poorly with binary data. So we're sticking with gam for surface activity and surface feeding



### PACKAGES
packages <- c("tidyverse", "ppcor","RColorBrewer", "scales","survMisc", "Metrics", "corrplot", "lme4", "MASS", "mgcv", "tidymv", "mgcViz", "gridExtra", "gratia", "ggcorrplot", "cowplot", "ggpubr")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)


#---------------- FUNCTIONS + THEME --------------------

source("code/functions.R")
source("code/themes.R")


#---------------- DATA --------------------

# focal vessel
tot <- read.csv("intermediate-products/response-var-dfs/surface-active+feeding_focal.csv") %>%
  dplyr::select(-c("lat.whale", "lon.whale")) %>% drop_na() %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(dist<1000) %>% # more likely to see surfacings that are breaches at longer distances
  filter(vess.accel.min.300>-10) %>% # temporary filtering until boat stuff properly processes
  filter(vess.accel.max.300<0.5) # filtering out until boat stuff properly processes

# quick look at the data set
colSums(is.na(tot)) 


#---------------- TRANSFORMATIONS --------------------

# define explanatory variables, then view distribution
vars <- c("dist", "seastate", "group", "day.of.year", "vess.di.60", "vess.accel.60", "vess.speed.60", "vess.speed.var.60", "vess.accel.max.300")
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
# dist: log transformation
ggplot(data = tot, aes(x = log(dist))) + geom_histogram()

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains("speed.var"), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # log speed var
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"), # arcsin DI
         log.dist = log(dist))


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group", "year", "day.of.year", "log.dist", "arcsin.vess.di.60",  "vess.speed.60", "log.vess.speed.var.60","vess.accel.60", "vess.accel.max.60")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_collinear.png", scale = 1.2)

# as expected, correlated accel and accel.max
# going for accel max for now


#---------------- TIME/DISTANCE FRAME --------------------

# first defining explanatory variables (not running folnum.unique as random in this)
var.fac = c("group", "seastate", "year") # factor
var.rand = c("folnum.unique") # random 
var.num = c("day.of.year", "log.dist") # numeric 
# specifying time/distance frame
ts = c(60, 300) # time/distance frame
var.t = c("vess.speed", "log.vess.speed.var", "arcsin.vess.di", "vess.accel.max") # varying vars!

# Combining together into a lagged data frame
t.df <- data.frame(var = rep(var.t, each = length(ts)), t = rep(ts,length(var.t)), t.orig = 60) %>%
  mutate(var.orig = paste0(var,".",t.orig), var.target = paste0(var,".",t))

(var.dy <- c(var.num, unique(t.df$var.orig)))

# now running the loop!
for (i in 1:nrow(t.df)) {
  
  # variable list, different time/distance frame each time
  v <- var.dy %>% recode(!!t.df[i,"var.orig"] := t.df[i,"var.target"])
  
  f1 <- as.formula(paste("surface.feeding ~", 
                         paste("s(",v,", bs = 'cs', k=10)", collapse= "+"),"+", # smooths
                         paste(var.fac, collapse = "+")))
  
  # run binomial GAM
  gam <- do.call("gam", list(as.formula(f1), data=as.name("tot"), method = "GCV.Cp", family = binomial))
  
  # and add AIC to df!
  t.df[i,"AIC"] = AIC(gam)
  
  print (t.df[i,"var.target"]) # progress checker
}

t.df # di: 300, speed.var: 60, accel.max: 300, speed: 60
ggplot(data = t.df, aes(x = as.factor(t), y = AIC, group = var, color = var)) +
  geom_point(size = 2, alpha = 0.5) + geom_line(alpha = 0.5) + 
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Time frame (seconds)", color = "Variable") + plot.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_time_aic.png", dpi = 600, height = 5, width = 8)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year") # factor
var.num = c("day.of.year", "log.dist", "vess.speed.60",  "log.vess.speed.var.60", "arcsin.vess.di.300", "vess.accel.max.300") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 65.4%, ubre = -0.8075. Remove year


# Remove year
var.fac = c("group", "seastate")
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 67.9%, ubre = -0.8121. good and good. removing

# Remove group
var.fac = c("seastate")
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 69%, ubre = -0.8122. good and good. removing

# Remove sea state
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 63.4%, ubre = -0.8032. bad and bad. not removing
var.fac = c("seastate")

# Remove vess.accel.max.300
var.num = c("day.of.year", "log.dist", "vess.speed.60",  "log.vess.speed.var.60", "arcsin.vess.di.300")
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 69%, ubre = -0.8122. no change. removing

# Remove log.dist
var.num = c("day.of.year", "vess.speed.60",  "log.vess.speed.var.60", "arcsin.vess.di.300")
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 69%, ubre = -0.8122. no change. removing

# Remove log.vess.speed.var.60
var.num = c("day.of.year", "vess.speed.60",  "arcsin.vess.di.300")
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 69%, ubre = -0.8122. no change. removing

# final model
var.fac = c("seastate") # factor
var.num = c("day.of.year", "vess.speed.60",  "arcsin.vess.di.300") # numeric 
var.rand = c("folnum.unique") # random 

f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 69%, ubre = -0.8122

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_final.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_final.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_results.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
print(plot(b, allTerms = T), pages = 1)
var.fac = c("seastate") 
var.num = c("day.of.year", "vess.speed.60",  "arcsin.vess.di.300")
var.rand = c("folnum.unique")

# day.of.year
(sm <- plot( sm(b, 1) ) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +
    partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "surface.feeding", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Rate of surface feeding") + partial.theme + partial.grid.margin + 
    coord_cartesian(ylim = c(0,1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_julian.png", height = 4.5, width = 13)


# speed
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Vessel speed (m/s, 60 sec)", y = "s(vessel speed)") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "surface.feeding", var = "vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = vess.speed.60, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel speed (m/s, 60 sec)", y = "Rate of surface feeding") + partial.theme + partial.grid.margin +
    coord_cartesian(ylim = c(0,0.7)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_vess.speed.60.png", height = 4.5, width = 13)


# DI
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "arcsin(vessel DI (300 sec))", y = "s(arcsin(vessel DI))") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "surface.feeding", var = "arcsin.vess.di.300", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sin(arcsin.vess.di.300), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel DI (300 sec)", y = "Rate of surface feeding") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = asin.trans, breaks = c(0, 0.4,0.6, 0.8, 0.9, 0.99)) + coord_cartesian(ylim = c(0,0.4)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_vess.di.300.png", height = 4.5, width = 13)


# seastate
v <- "seastate"
plotdata <- plot(b, allTerms = T, select = 5)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "surface.feeding", var = "seastate", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = seastate, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Rate of surface feeding") + partial.theme + partial.grid.margin + 
    coord_cartesian(ylim = c(0,0.04)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_seastate.png", height = 4.5, width = 13)
