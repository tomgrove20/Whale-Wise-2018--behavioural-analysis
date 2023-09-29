#################################
### GAM: SURFACING IBI, FOCAL ###
#################################
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
  filter(!is.na(ibi.mean) & !is.na(dist.mean)) %>%
  mutate_at(vars(contains("dt")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate.last", "year", "folnum.unique"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1)  %>% # important dummy variable to factor out random effect when predicting response
  # dplyr::select(-contains(".10")) %>%
  # dplyr::select(-contains(c("0_end", "0_start"))) %>%
  rename_at(vars(contains("_mean")), funs(str_replace(.,"_mean",""))) %>%
  filter(vess.accel.min.300>-10) %>% # temporary filtering until boat stuff properly processes
  filter(vess.accel.max.300<0.5) # filtering out until boat stuff properly processes

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 66 surfacing intervals, 40 follows


#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = ibi.mean)) + geom_histogram()
ggplot(data = tot, aes(x = log(ibi.mean))) + geom_histogram() # perhaps quasi-poisson with log link?

# define explanatory variables, then view distribution
vars <- c("dist.end", "seastate.last", "dive.time.before", "group", "day.of.year", "vess.di.60_mean", "vess.accel.60_mean", "vess.speed.60_mean", "vess.speed.var.60_mean", "vess.accel.max.60_mean")
for (i in vars) {
  if(class(tot[,i]) == "numeric"){ # if numeric, provide a histogram
    print(ggplot(data = tot, aes(x = get(i))) + geom_histogram() + xlab(i))
  } else { # if not numeric, provide a bar plot
    print(ggplot(data = tot, aes(x = get(i))) + geom_bar() + xlab(i))
  }
  readline(prompt = "next plot")
}


# dist.end
ggplot(data = tot, aes(x = log(dist.end))) + geom_histogram()
# dive.time.before
ggplot(data = tot, aes(x = log(dive.time.before))) + geom_histogram()
# di
ggplot(data = tot, aes(x = asin(vess.di.60_mean))) + geom_histogram()
# speed.var?
ggplot(data = tot, aes(x = sqrt(vess.speed.var.60_mean))) + geom_histogram()

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("dive.time","dist")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # sqrt mean dist
         across(contains("speed.var"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt dist max
         across(contains("vess.di."), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))





#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate.last", "group", "surface.feeding.follow", "surface.active.follow", "year", "day.of.year", "log.dive.time.before", "log.dist.end", "arcsin.vess.di.60", "vess.speed.60", "sqrt.vess.speed.var.60","vess.accel.60", "vess.accel.max.60")

# correlation matrix
x <- tot[,vars] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/ibi-mean/foc/gam_ibi-mean_foc_collinear.png", scale = 1.2)

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
                         paste("s(",v,", bs = 'cs', k=3)", collapse= "+"),"+", # smooths
                         paste(var.fac, collapse = "+")))
  
  # run binomial GAM
  gam <- do.call("gam", list(as.formula(f1), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
  
  # and add AIC to df!
  t.df[i,"AIC"] = AIC(gam)
  
  print (t.df[i,"var.target"]) # progress checker
}

t.df # di: 60, speed.var: 300, accel.max: 300, speed: 300
ggplot(data = t.df, aes(x = as.factor(t), y = AIC, group = var, color = var)) +
  geom_point(size = 2, alpha = 0.5) + geom_line(alpha = 0.5) + 
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Time frame (seconds)", color = "Variable") + plot.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("intermediate-products/gam/ibi-mean/foc/gam_ibi-mean_foc_time_dev+gcv.png", dpi = 600, height = 8, width = 8)

#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.num = c("day.of.year","log.dive.time.before",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.300", "arcsin.vess.di.300", "vess.accel.max.60", "vess.accel.60") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("ibi.mean ~", 
                      paste("s(",var.num,", bs = 'cs', k=3)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 93.8%, gcv = 4.507


# remove surface.active.follow, surface.feeding.follow, seastate, group, year
f <- as.formula(paste("ibi.mean ~", 
                      paste("s(",var.num,", bs = 'cs', k=3)", collapse= "+"),"+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 94.5%, gcv = 3.79


# remove speed.var, vess.di, accel.max, accel, day.of.year, log.dive.time.before
var.num = c("log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.300")
f <- as.formula(paste("ibi.mean ~", 
                      paste("s(",var.num,", bs = 'cs', k=3)", collapse= "+"),"+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 94.2%, gcv = 3.79. no change, remove

# final model
var.num = c("log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.300")
var.rand = c("folnum.unique")

f <- as.formula(paste("ibi.mean ~", 
                      paste("s(",var.num,", bs = 'cs', k=3)", collapse= "+"),"+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 94.2%, gcv = 3.791

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam/ibi-mean/foc/gam_ibi-mean_foc_final.rds")



#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam/ibi-mean/foc/gam_ibi-mean_foc_final.rds")

# defining variable combinations
var.num = c("log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.300")
var.rand = c("folnum.unique")
var.fac = NA

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam/ibi-mean/foc/gam_ibi-mean_foc_results.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
print(plot(b, allTerms = T), pages = 1)

# log.dist.mean
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(mean distance (m))", y = "s(log(mean distance))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "ibi.mean", var = "log.dist.mean", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.dist.mean), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Mean distance (m)", y = "IBI mean (sec)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(50,100,200,400,800)) +
    scale_y_continuous(breaks = seq(0,40,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/ibi-mean/foc/gam_ibi-mean_foc_dist.mean.png", height = 4.5, width = 13)

# vess.speed.60
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Mean vessel speed (m/s, 60 sec)", y = "s(mean vessel speed)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "ibi.mean", var = "vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = vess.speed.60, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Mean vessel speed (m/s, 60 sec)", y = "IBI mean (sec)") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,40,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/ibi-mean/foc/gam_ibi-mean_foc_vess.speed.60.png", height = 4.5, width = 13)

# sqrt.vess.speed.var.300
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(mean vessel speed SD (m/s, 300 sec))", y = "s(sqrt(mean vessel speed SD))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "ibi.mean", var = "sqrt.vess.speed.var.300", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.vess.speed.var.300^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Mean vessel speed SD (m/s, 300 sec)", y = "IBI mean (sec)") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,40,2)) + scale_x_continuous(trans = "sqrt"))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/ibi-mean/foc/gam_ibi-mean_foc_vess.speed.var.300.png", height = 4.5, width = 13)
