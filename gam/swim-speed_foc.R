##############################
### GAM: SWIM SPEED, FOCAL ###
##############################
### 18/02/2022
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


#---------------- DATA --------------------

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
  filter(vess.speed.var.60<3) # filtering until properly processed

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 737 surfacing pairs, 126 follows


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
ggplot(data = tot, aes(x = log(vess.speed.var.60))) + geom_histogram()
# vess.accel.max
ggplot(data = tot, aes(x = vess.accel.max.60)) + geom_histogram()

# accel  max not ideal but difficult to transform since crossing 0. DI is left-skewed and speed var is right skewed

# speed.var: log transformation
sort(tot$vess.speed.var.60); ggplot(data = tot, aes(x = log(vess.speed.var.60))) + geom_histogram()
# di: arcsin transformation
ggplot(data = tot, aes(x = asin(vess.di.60))) + geom_histogram()
# dist: log transformation
ggplot(data = tot, aes(x = log(dist))) + geom_histogram()

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("diff.time","vess.speed.var")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(ends_with(c("speed.60", "speed.300", "dist")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group","surface.feeding", "surface.active", "diff.time", "year", "day.of.year", "sqrt.dist", "arcsin.vess.di.60",  "sqrt.vess.speed.60", "log.vess.speed.var.60","vess.accel.60", "vess.accel.max.60")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_collinear.png", scale = 1.2)


#---------------- TIME/DISTANCE FRAME --------------------

# first defining explanatory variables (not running folnum.unique as random in this)
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") # factor
var.rand = c("folnum.unique") # random 
var.num = c("day.of.year", "sqrt.dist", "diff.time", "vess.accel.60") # numeric 
# specifying time/distance frame
ts = c(60, 300) # time/distance frame
var.t = c("sqrt.vess.speed", "log.vess.speed.var", "arcsin.vess.di", "vess.accel.max") # varying vars!

# Combining together into a lagged data frame
t.df <- data.frame(var = rep(var.t, each = length(ts)), t = rep(ts,length(var.t)), t.orig = 60) %>%
  mutate(var.orig = paste0(var,".",t.orig), var.target = paste0(var,".",t))

(var.dy <- c(var.num, unique(t.df$var.orig)))

# now running the loop!
for (i in 1:nrow(t.df)) {
  
  # variable list, different time/distance frame each time
  v <- var.dy %>% recode(!!t.df[i,"var.orig"] := t.df[i,"var.target"])
  
  f1 <- as.formula(paste("speed ~", 
                         paste("s(",v,", bs = 'cs', k=10)", collapse= "+"),"+", # smooths
                         paste(var.fac, collapse = "+")))
  
  # run binomial GAM
  gam <- do.call("gam", list(as.formula(f1), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
  
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
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_time_aic.png", dpi = 600, height = 5, width = 8)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") # factor
var.num = c("day.of.year", "sqrt.dist", "log.diff.time", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.60", "vess.accel.max.300", "vess.accel.60") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("speed ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 45.4%, GCV = 4.137

# remove year
var.fac = c("group", "seastate", "surface.feeding", "surface.active") # factor
f <- as.formula(paste("speed ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 44.6%, GCV = 4.138. bad and bad, don't remove.
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active")



# remove accel.max
var.num = c("day.of.year", "sqrt.dist", "log.diff.time", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.60", "vess.accel.60")
f <- as.formula(paste("speed ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 45.4%, GCV = 4.137. No change, so remove!


# remove group
var.fac = c("seastate", "year", "surface.feeding", "surface.active")
f <- as.formula(paste("speed ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 45.3%, GCV = 4.129. tiny bit bad and quite good, so removing

# final model
var.fac = c("seastate", "year", "surface.feeding", "surface.active")
var.num = c("day.of.year", "sqrt.dist", "log.diff.time", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.60", "vess.accel.60")
var.rand = c("folnum.unique")

f <- as.formula(paste("speed ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 43.4%, ubre = 11.06.
# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_final.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_final.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_results.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
var.fac = c("seastate", "year", "surface.feeding", "surface.active")
var.num = c("day.of.year", "sqrt.dist", "log.diff.time", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.60", "vess.accel.60")
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)

# day.of.year
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(1,10, by = 0.5)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_julian.png", height = 4.5, width = 13)


# dist
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(distance (m))", y = "s(sqrt(distance))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "sqrt.dist", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.dist^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Distance (m)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(20,50,100,200,400,800)) +
    scale_y_continuous(breaks = seq(1,10, by = 0.5)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_dist.png", height = 4.5, width = 13)


# log.diff.time
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(inter-breath interval (sec))", y = "s(inter-breath interval)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "log.diff.time", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.diff.time), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Inter-breath interval (sec)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(5,10,20,50,100,200,400,800)) +
    scale_y_continuous(breaks = seq(1,10, by = 1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_ibi.png", height = 4.5, width = 13)


# vess.speed
(sm <- plot(sm(b, 4)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(vessel speed (m/s, 60 sec))", y = "s(sqrt(vessel speed))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "sqrt.vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.vess.speed.60^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel speed (m/s, 60 sec)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(0.2, 0.5, 1, 2, 3,4, 5)) +
    scale_y_continuous(breaks = seq(1,10, by = 0.5)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_vess.speed.60.png", height = 4.5, width = 13)


# log.vess.speed.var.300
(sm <- plot(sm(b, 5)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(vessel speed SD (m/s, 300 sec))", y = "s(log(vessel speed SD))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "log.vess.speed.var.300", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.vess.speed.var.300), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel speed SD (m/s, 300 sec)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2)) +
    scale_y_continuous(breaks = seq(1,10, by = 0.5)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_vess.speed.var.300.png", height = 4.5, width = 13)


# DI
(sm <- plot(sm(b, 6)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "arcsin(vessel DI (60 sec))", y = "s(arcsin(vessel DI))") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "arcsin.vess.di.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sin(arcsin.vess.di.60), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel DI (60 sec)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = asin.trans, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99,1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_vess.di.60.png", height = 4.5, width = 13)


# vess.accel.60
(sm <- plot(sm(b, 7)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = bquote("Vessel accleration ("*m/s^2*", 60 sec)"), y = "s(vessel acceleration)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "vess.accel.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = vess.accel.60, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = bquote("Vessel accleration ("*m/s^2*", 60 sec)"), y = "Swim speed (km/h)") + partial.theme + partial.grid.margin + scale_y_continuous(breaks = seq(0,16, by = 2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_vess.accel.60.png", height = 4.5, width = 13)


# seastate
v <- "seastate"
plotdata <- plot(b, allTerms = T, select = 9)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "seastate", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = seastate, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_seastate.png", height = 4.5, width = 13)


# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 10)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +  
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,16, 2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_year.png", height = 4.5, width = 13)

# surface.feeding
v <- "surface.feeding"
plotdata <- plot(b, allTerms = T, select = 11)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
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
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_surface.feeding.png", height = 4.5, width = 13)


# surface.active
v <- "surface.active"
plotdata <- plot(b, allTerms = T, select = 12)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +  
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface activity", y = "Partial effect of surface activity"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "surface.active", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.active, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Surface activity", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")) + scale_y_continuous(breaks = seq(0,20,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_surface.active.png", height = 4.5, width = 13)
