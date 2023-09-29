#####################################
### GAM: SURFACING BREATHS, FOCAL ###
#####################################
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


#---------------- DATA --------------------

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
  filter(vess.accel.max.300 < 5) %>% # filtering until properly processed
  filter(vess.speed.var.300 < 3) # filtering until properly processed
# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 142 surfacing intervals, 56 follows


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
ggplot(data = tot, aes(x = sqrt(vess.speed.var.300))) + geom_histogram()
# cannot transform accel.max (either side of 0!)

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("dive.time","dist")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(contains("speed.var"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains("vess.di."), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate.last", "group", "surface.feeding.follow", "surface.active.follow", "year", "day.of.year", "log.dive.time.before", "log.dist.mean", "arcsin.vess.di.60", "vess.speed.60", "sqrt.vess.speed.var.60","vess.accel.60", "vess.accel.max.60")

# correlation matrix
x <- tot[,vars] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_collinear.png", scale = 1.2)

# actually no troubling correlations at all! Let's keep them in!


#---------------- TIME/DISTANCE FRAME --------------------

# first defining explanatory variables (not running folnum.unique as random in this)
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.rand = c("folnum.unique") # random 
var.num = c("day.of.year", "log.dist.mean", "vess.accel.60", "log.dive.time.before") # numeric 
# specifying time/distance frame
ts = c(60, 300) # time/distance frame
var.t = c("sqrt.vess.speed.var", "arcsin.vess.di", "vess.accel.max", "vess.speed") # varying vars!

# Combining together into a lagged data frame
t.df <- data.frame(var = rep(var.t, each = length(ts)), t = rep(ts,length(var.t)), t.orig = 60) %>%
  mutate(var.orig = paste0(var,".",t.orig), var.target = paste0(var,".",t))

(var.dy <- c(var.num, unique(t.df$var.orig)))

# now running the loop!
for (i in 1:nrow(t.df)) {
  
  # variable list, different time/distance frame each time
  v <- var.dy %>% recode(!!t.df[i,"var.orig"] := t.df[i,"var.target"])
  
  f1 <- as.formula(paste("ibi.mean ~", 
                         paste("s(",v,", bs = 'cs', k=5)", collapse= "+"),"+", # smooths
                         paste(var.fac, collapse = "+")))
  
  # run binomial GAM
  gam <- do.call("gam", list(as.formula(f1), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
  
  # and dev and gcv to df!
  t.df[i,"Deviance explained"] = summary(gam)$dev.expl
  t.df[i,"GCV score"] = gam$gcv.ubre
  
  print (t.df[i,"var.target"]) # progress checker
}

t.df # di: 60, speed.var: 60, accel.max: 60, speed: 60
t.df %>% gather(key = "type", value = "score", 'Deviance explained', 'GCV score') %>%
  ggplot(aes(x = as.factor(t), y = score, group = var, color = var)) +
  geom_point(size = 2, alpha = 0.5) + geom_line() + facet_wrap(~type, scales = "free_y", ncol = 1) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Time frame (seconds)", color = "Variable") + plot.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_time_dev+gcv.png", dpi = 600, height = 8, width = 8)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.num = c("day.of.year","log.dive.time.before",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "vess.accel.max.60", "vess.accel.60") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
summary(gam) # dev = 92.6%, gcv = 0.3474


# Remove seastate.last
var.fac = c("group", "year", "surface.feeding.follow", "surface.active.follow")
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
summary(gam) # dev = 92.5%, gcv = 0.3448. bit bad and good. removing


# Remove surface.feeding.follow
var.fac = c("group", "year", "surface.active.follow")
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
summary(gam) # dev = 93.1%, gcv = 0.34071. good and very good. removing


# Remove year
var.fac = c("group", "surface.active.follow")
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
summary(gam) # dev = 92.2%, gcv = 0.3195. quite bad and very good. not removing
var.fac = c("group", "year", "surface.active.follow")


# Remove day.of.year
var.num = c("log.dive.time.before",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "vess.accel.max.60", "vess.accel.60")
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
summary(gam) # dev = 92.4%, gcv = 0.34106. bad and bad. not removing
var.num = c("log.dive.time.before", "day.of.year",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "vess.accel.max.60", "vess.accel.60")


# final model
var.fac = c("group", "year",  "surface.active.follow")
var.num = c("log.dive.time.before", "day.of.year",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "vess.accel.max.60", "vess.accel.60")
var.rand = c("folnum.unique")

f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
summary(gam) # dev = 93.1%, gcv = 0.34071

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_final.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_final.rds")

# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.num = c("day.of.year","log.dive.time.before",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "vess.accel.max.60", "vess.accel.60") # numeric 
var.rand = c("folnum.unique") # random 

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_results.csv", row.names = FALSE)

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
    ggplot(aes(x = exp(log.dive.time.before), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Previous dive time (sec)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(50,100,200,400,800)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_dive.time.before.png", height = 4.5, width = 13)


# day.of.year
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_julian.png", height = 4.5, width = 13)

# log.dist.mean
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(mean distance(m))", y = "s(log(mean distance))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "log.dist.mean", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.dist.mean), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Mean distance (m)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(20, 50,100,200,500,1000)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_dist.mean.png", height = 4.5, width = 13)


# vess.speed.60
(sm <- plot(sm(b, 4)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Mean vessel speed (m/s, 60 sec)", y = "s(mean vessel speed)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = vess.speed.60, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Mean vessel speed (m/s, 60 sec)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_vess.speed.60.png", height = 4.5, width = 13)


# sqrt.vess.speed.var.60
(sm <- plot(sm(b, 5)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(mean vessel speed SD (m/s, 60 sec))", y = "s(sqrt(mean vessel speed SD))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "sqrt.vess.speed.var.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.vess.speed.var.60^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Mean vessel speed SD (m/s, 60 sec)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt"))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_vess.speed.var.60.png", height = 4.5, width = 13)


# arcsin.vess.di.60
(sm <- plot(sm(b, 6)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "arcsin(mean vessel DI, 60 sec)", y = "s(arcsin(mean vessel DI))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "arcsin.vess.di.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sin(arcsin.vess.di.60), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Mean vessel DI (60 sec)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "asn", breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_vess.di.60.png", height = 4.5, width = 13)


# vess.accel.max.60
(sm <- plot(sm(b, 7)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = bquote("Max vessel accleration ("*m/s^2*", 60 sec)"), y = "s(max vessel acceleration)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "vess.accel.max.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = vess.accel.max.60, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = bquote("Max vessel accleration ("*m/s^2*", 60 sec)"), y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_vess.accel.max.60.png", height = 4.5, width = 13)


# vess.accel.60
(sm <- plot(sm(b, 8)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = bquote("Mean vessel accleration ("*m/s^2*", 60 sec)"), y = "s(mean vessel acceleration)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "vess.accel.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = vess.accel.60, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = bquote("Mean vessel accleration ("*m/s^2*", 60 sec)"), y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_vess.accel.60.png", height = 4.5, width = 13)


# group
plotdata <- plot(b, allTerms = T, select = 10)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- ggplot(data = plotdata, aes(x = x, y = y)) + geom_point(size = 3) + 
    geom_errorbar(aes(ymax = y+ci, ymin = y-ci), linetype = "dashed", width = 0.4) + 
    geom_vpline(data = tot, aes(x = group, y = jitter.pos), position = position_jitter(h = 0), height =jitter.height, size = 0.001) +
    partial.theme + partial.grid.margin +  scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "breath.num", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    coord_cartesian(ylim = c(0,5)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_group.png", height = 4.5, width = 13)


# year
plotdata <- plot(b, allTerms = T, select = 11)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- ggplot(data = plotdata, aes(x = x, y = y)) + geom_point(size = 3) + 
    geom_errorbar(aes(ymax = y+ci, ymin = y-ci), linetype = "dashed", width = 0.4) + 
    geom_vpline(data = tot, aes(x = year, y = jitter.pos), position = position_jitter(h = 0), height =jitter.height, size = 0.001) +
    partial.theme + partial.grid.margin +  scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "breath.num", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    coord_cartesian(ylim = c(0,13)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_year.png", height = 4.5, width = 13)


# surface.active.follow
plotdata <- plot(b, allTerms = T, select = 12)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- ggplot(data = plotdata, aes(x = x, y = y)) + geom_point(size = 3) + 
    geom_errorbar(aes(ymax = y+ci, ymin = y-ci), linetype = "dashed", width = 0.4) + 
    geom_vpline(data = tot, aes(x = surface.active.follow, y = jitter.pos), position = position_jitter(h = 0), height =jitter.height, size = 0.001) +
    partial.theme + partial.grid.margin +  scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface activity in follow", y = "Partial effect of surface activity"))
(pr <- pred.var.fac(gam, tot, resp = "ibi.mean", var = "surface.active.follow", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.active.follow, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Surface activity in follow", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/foc/gam_breaths-surfacing_foc_surface.active.follow.png", height = 4.5, width = 13)



