######################
### GAM: DI, FOCAL ###
######################
### 20/02/2022
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

# foc
tot <- read.csv("intermediate-products/response-var-dfs/di_focal.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding", "surface.active"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(difftime.next<2000) %>% # until properly processed
  filter(vess.accel.max.300 < 2) %>% # filtering until properly processed
  filter(vess.accel.60 >-0.2) %>% # filtering until properly processed
  filter(vess.speed.var.60<3) # filtering until properly processed

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 598 surfacing pairs, 100 follows


#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = DI)) + geom_histogram()
ggplot(data = tot, aes(x = asin(DI))) + geom_histogram() # going to arcsin transform this simply to reduce skew

tot <- tot %>%
  mutate(arcsin.DI = asin(DI))


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
ggplot(data = tot, aes(x = log(vess.speed.var.60))) + geom_histogram()
# vess.accel.max
ggplot(data = tot, aes(x = vess.accel.max.60)) + geom_histogram()


# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("difftime","vess.speed.var")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(ends_with(c("speed.60", "speed.300", "dist")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group","surface.feeding", "surface.active", "log.difftime.prev", "log.difftime.next", "year", "day.of.year", "sqrt.dist", "arcsin.vess.di.60",  "sqrt.vess.speed.60", "log.vess.speed.var.60","vess.accel.60", "vess.accel.max.60")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_collinear.png", scale = 1.2)

# as expected, correlated accel and accel.max but only 0.6. Keeping in for now!


#---------------- TIME/DISTANCE FRAME --------------------

# first defining explanatory variables (not running folnum.unique as random in this)
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") # factor
var.rand = c("folnum.unique") # random 
var.num = c("day.of.year", "sqrt.dist", "log.difftime.prev", "log.difftime.next", "vess.accel.60") # numeric 
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
  
  f1 <- as.formula(paste("arcsin.DI ~", 
                         paste("s(",v,", bs = 'cs', k=10)", collapse= "+"),"+", # smooths
                         paste(var.fac, collapse = "+")))
  
  # run binomial GAM
  gam <- do.call("gam", list(as.formula(f1), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "identity")))
  
  # and add AIC to df!
  t.df[i,"AIC"] = AIC(gam)
  
  print (t.df[i,"var.target"]) # progress checker
}

t.df # di: 60, speed.var: 300, accel.max: 60, speed: 60
ggplot(data = t.df, aes(x = as.factor(t), y = AIC, group = var, color = var)) +
  geom_point(size = 2, alpha = 0.5) + geom_line(alpha = 0.5) + 
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Time frame (seconds)", color = "Variable") + plot.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_time_aic.png", dpi = 600, height = 5, width = 8)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") # factor
var.num = c("day.of.year", "sqrt.dist", "log.difftime.prev", "log.difftime.next", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.300", "vess.accel.max.60", "vess.accel.60") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "identity")))
summary(gam) # dev = 27.5%, GCV = 0.083397


# try removing year
var.fac = c("group", "seastate", "surface.feeding", "surface.active") # factor
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link= "identity")))
summary(gam) # dev = 27.8%, GCV = 0.083256. Remove


# remove vess.accel.60
var.num = c("day.of.year", "sqrt.dist", "log.difftime.prev", "log.difftime.next", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.300", "vess.accel.max.60")
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link= "identity")))
summary(gam) # dev = 27.2%, gcv = 0.083065. Bad and good, but it really doesn't do anything. Remove

# remove surface.active
var.fac = c("group", "seastate", "surface.feeding") # factor
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link= "identity")))
summary(gam) # dev = 27.3%, GCV = 0.082813. Good and good. Remove

# remove surface feeding
var.fac = c("group", "seastate") # factor
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link= "identity")))
summary(gam) # dev = 27.6%, GCV = 0.082942. Good and no real change. Remove

# final model
var.fac = c("group", "seastate")
var.num = c("day.of.year", "sqrt.dist", "log.difftime.prev", "log.difftime.next", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.300", "vess.accel.max.60")
var.rand = c("folnum.unique")

f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link= "identity")))
summary(gam) # dev = 27.6%, GCV = 0.082942


# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam/di/foc/gam_di_foc_final.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam/di/foc/gam_di_foc_final.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam/di/foc/gam_di_foc_results.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
var.fac = c("group", "seastate")
var.num = c("day.of.year", "sqrt.dist", "log.difftime.prev", "log.difftime.next", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.300", "vess.accel.max.60")
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)

# day of year
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_julian-day.png", height = 4.5, width = 13)


# dist
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(distance (m))", y = "s(sqrt(distance))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "sqrt.dist", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.dist^2, y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Distance (m)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(20,50,100,200,400,800)) +
    scale_y_continuous(trans = asin.trans, breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_dist.png", height = 4.5, width = 13)

# difftime.prev
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(previous IBI (sec))", y = "s(log(previous IBI))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "log.difftime.prev", var.fac, var.num, var.rand) %>% 
    mutate(upr = ifelse(exp(log.difftime.prev)>200, asin(1), upr)) %>%
    ggplot(aes(x = exp(log.difftime.prev), y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Previous IBI (sec)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.8, 0.9, 0.95, 0.99, 1)) +
    scale_x_continuous(trans = "log", breaks = c(5, 10, 20, 50, 100, 200, 400, 800)) +
    coord_cartesian(ylim = c(0.79, 0.9994)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_difftime.prev.png", height = 4.5, width = 13)

# difftime.next
(sm <- plot(sm(b, 4)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(following IBI (sec))", y = "s(log(following IBI))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "log.difftime.next", var.fac, var.num, var.rand) %>% 
    mutate(upr = ifelse(exp(log.difftime.next)>320, asin(1), upr)) %>%
    ggplot(aes(x = exp(log.difftime.next), y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Following IBI (sec)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.8, 0.85, 0.9, 0.95, 0.99, 1)) +
    scale_x_continuous(trans = "log", breaks = c(5, 10, 20, 50, 100, 200, 400, 800)) +
    coord_cartesian(ylim = c(0.8, 0.998)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_difftime.next.png", height = 4.5, width = 13)

# vess.speed.60
(sm <- plot(sm(b, 5)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(vessel speed (m/s, 60 sec))", y = "s(sqrt(vessel speed))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "sqrt.vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.vess.speed.60^2, y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel speed (m/s, 60 sec)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.9, 0.95, 0.98, 0.99, 1)) +
    scale_x_continuous(trans = "sqrt", breaks = c(0.1, 0.5, 1,2,3,4,6,8)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_vess.speed.60.png", height = 4.5, width = 13)

# vess.speed.var.300
(sm <- plot(sm(b, 6)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(vessel speed SD (m/s, 300 sec))", y = "s(log(vessel speed SD))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "log.vess.speed.var.300", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.vess.speed.var.300), y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel speed SD (m/s, 300 sec)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.8, 0.85, 0.9, 0.95, 0.99, 1)) +
    scale_x_continuous(trans = "log", breaks = c(0.1, 0.2, 0.5, 1, 2, 5)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_vess.speed.var.300.png", height = 4.5, width = 13)


# vess.di.300
(sm <- plot(sm(b, 7)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "arcsin(vessel DI (300 sec))", y = "s(arcsin(vessel DI))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "arcsin.vess.di.300", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sin(arcsin.vess.di.300), y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel DI (300 sec)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.7, 0.8, 0.85, 0.9, 0.95, 0.99, 1)) +
    scale_x_continuous(trans = asin.trans, breaks = c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_vess.di.300.png", height = 4.5, width = 13)


# vess.accel.max.60
(sm <- plot(sm(b, 8)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = bquote("Max vessel accleration ("*m/s^2*", 60 sec)"), y = "s(max vessel acceleration)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "vess.accel.max.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = vess.accel.max.60, y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = bquote("Max vessel accleration ("*m/s^2*", 60 sec)"), y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.7, 0.8, 0.9, 0.95, 0.99, 1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_vess.accel.max.60.png", height = 4.5, width = 13)


# group
v <- "group"
plotdata <- plot(b, allTerms = T, select = 10)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.5, 0.7,0.8, 0.9, 0.95, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_group.png", height = 4.5, width = 13)

# sea state
v <- "seastate"
plotdata <- plot(b, allTerms = T, select = 11)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "seastate", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = seastate, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.5, 0.7,0.8, 0.9, 0.95, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/di/foc/gam_di_foc_seastate.png", height = 4.5, width = 13)
