############################
### GAM: SWIM SPEED, AIS ###
############################
### 18/02/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

# 1. transformations
# 2. Collinearity
# 3. time frame selection
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

# AIS
tot <- read.csv("intermediate-products/response-var-dfs/speed_ais.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding", "surface.active"), as.factor) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(speed<20) %>% # speeds above 20 km/h are unrealistic
  filter(diff.time<2000) # filtering until data processed properly

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 4155 surfacing pairs, 530 follows


#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = speed)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(speed))) + geom_histogram() # so Gaussian with log link?

# define explanatory variables, then view distribution
vars <- c("seastate", "group","surface.feeding", "day.of.year", "diff.time", "oak.meandist.30.1500", "rib.meandist.30.1500", "oak.dist.max.30.1500", "rib.dist.max.30.1500", "oak.num.30.1500", "rib.num.30.1500")
for (i in vars) {
  if(class(tot[,i]) == "numeric"){ # if numeric, provide a histogram
    print(ggplot(data = tot, aes(x = get(i))) + geom_histogram() + xlab(i))
  } else { # if not numeric, provide a bar plot
    print(ggplot(data = tot, aes(x = get(i))) + geom_bar() + xlab(i))
  }
  readline(prompt = "next plot")
}

# thinking of transforming oak and rib distances, annoying with all the 0s
ggplot(data = tot, aes(x = oak.meandist.30.1500)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(oak.meandist.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = rib.meandist.30.1500)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(rib.meandist.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(oak.dist.max.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(rib.dist.max.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = log(diff.time))) + geom_histogram()

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("meandist", "dist.max")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt
         across(contains(c("diff.time")), .fns = list(log = ~log(.)), .names = "{fn}.{col}")) # sqrt dist max


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group", "surface.feeding", "surface.active", "diff.time", "day.of.year", "oak.num.30.1500", "rib.num.30.1500","sqrt.oak.meandist.30.1500", "sqrt.rib.meandist.30.1500", "sqrt.oak.dist.max.30.1500", "sqrt.rib.dist.max.30.1500",  "year")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 2.5) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_collinear.png")

# max.dist causing most problems
# rib.num ~ rib.meandist = 0.72 but not too bad


#---------------- TIME/DISTANCE FRAME --------------------

# first defining explanatory variables (not running folnum.unique as random in this)
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") # factor
var.rand = c("folnum.unique") # random 
var.num = c("day.of.year", "log.diff.time") # numeric 
# specifying time/distance frame
ts = c("30.1500", "10.500", "10.1000", "10.1500") # time/distance frame
var.num.t = c("sqrt.oak.meandist", "sqrt.rib.meandist") # varying continuous numeric vars
var.int.t = c("oak.num", "rib.num") # varying integer vars

# Combining together into a lagged data frame
t.df <- data.frame(
  var = c(rep(var.num.t, each = length(ts)),rep(var.int.t, each = length(ts))),
  type = c(rep("num",length(ts)*length(var.num.t)), rep("int", length(ts)*length(var.int.t))),
  t = rep(ts, length(var.num.t)+length(var.int.t)), t.orig = "10.1000") %>%
  mutate(var.orig = paste0(var,".",t.orig), var.target = paste0(var,".",t))

# specifying original vars for formula
(var.dy.num <- c(var.num, unique(filter(t.df, type == "num")$var.orig)))
(var.dy.int <- c(unique(filter(t.df, type == "int")$var.orig)))

# now running the loop!
for (i in 1:nrow(t.df)) {
  
  # variable list, different time/distance frame each time
  v.num <- var.dy.num %>% recode(!!t.df[i,"var.orig"] := t.df[i,"var.target"])
  v.int <- var.dy.int %>% recode(!!t.df[i,"var.orig"] := t.df[i,"var.target"])
  
  f1 <- as.formula(paste("speed ~", 
                         paste("s(",v.int,", bs = 'cs', k=3)", collapse= "+"),"+", # integer smooths
                         paste("s(",v.num,", bs = 'cs', k=10)", collapse = "+"), "+", # numeric smooths
                         paste(var.fac, collapse = "+")))
  
  # run binomial GAM
  gam <- do.call("gam", list(as.formula(f1), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
  
  # and add AIC to df!
  t.df[i,"AIC"] = AIC(gam)
  
  print (t.df[i,"var.target"]) # progress checker
}

t.df 
ggplot(data = t.df, aes(x = t, y = AIC, group = var, color = var)) +
  geom_point(size = 2, alpha = 0.5) + geom_line() + 
  scale_color_brewer(palette = "Dark2") +
  scale_x_discrete(labels=c("10 min, 500 m", "10 min, 1000 m", "10 min, 1500 m", "30 min, 1500 m")) +
  labs(x = "", color = "Variable") + plot.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_time-dist_aic.png", dpi = 600, height = 5, width = 8)

# going for 10.1500 for oak and 30.1500 for rib


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active")
var.int = c("oak.num.10.1500", "rib.num.30.1500")
var.num = c("day.of.year", "log.diff.time", "sqrt.oak.meandist.10.1500", "sqrt.rib.meandist.30.1500")
var.rand = c("folnum.unique")


# formula
f <- as.formula(paste("speed ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 30.7%, gcv = 4.7562

# remove year
var.fac = c("group", "seastate", "surface.feeding", "surface.active")
f <- as.formula(paste("speed ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 30.7%, gcv = 4.7546. Nothing and bit good, so removing

# remove seastate
var.fac = c("group", "surface.feeding", "surface.active")
f <- as.formula(paste("speed ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 30.7%, gcv = 4.7614 # no change and bit bad, so keeping

# remove oak.num.10.1500
var.fac = c("group", "seastate", "surface.feeding", "surface.active")
var.int = c("rib.num.30.1500")
f <- as.formula(paste("speed ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 30.7%, gcv = 4.7546 # no change, so removing


# remove rib.num.30.1500
f <- as.formula(paste("speed ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 30.7%, gcv = 4.756 # no change and a bit bad, so not removing

# final model
var.fac = c("group", "seastate", "surface.feeding", "surface.active")
var.int = c("rib.num.30.1500")
var.num = c("day.of.year", "log.diff.time", "sqrt.oak.meandist.10.1500", "sqrt.rib.meandist.30.1500")
var.rand = c("folnum.unique")
f <- as.formula(paste("speed ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 30.7%, gcv = 4.755

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_final.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_final.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_results.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
var.fac = c("group", "seastate", "surface.feeding", "surface.active")
var.num = c("rib.num.30.1500", "day.of.year", "log.diff.time", "sqrt.oak.meandist.10.1500", "sqrt.rib.meandist.30.1500")
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)

# rib.num.30.1500
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "RIB boats (30 min, 1500 m)", y = "s(RIB boats)") + partial.theme + partial.grid.margin + scale_x_continuous(breaks = c(0:10)))
(pr <- pred.var(gam, tot, resp = "speed", var = "rib.num.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = rib.num.30.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "RIB boats (30 min, 1500 m)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_rib.num.30.1500.png", height = 4.5, width = 13)

# day.of.year
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_julian.png", height = 4.5, width = 13)


# log.diff.time
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(inter-breath interval (sec))", y = "s(inter-breath interval)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "log.diff.time", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.diff.time), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Inter-breath interval (sec)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(5,10,20,50,100,200,400,800)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_ibi.png", height = 4.5, width = 13)


# oak.meandist.10.1500
(sm <- plot(sm(b, 4)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(oak boat mean distance (m in 10 min, 1500 m))", y = "s(sqrt(oak boat mean distance))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "sqrt.oak.meandist.10.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.oak.meandist.10.1500^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Oak boat mean distance (m in 10 min, 1500 m)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(20,200,500,1000,2000, 3000)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_oak.meandist.10.1500.png", height = 4.5, width = 13)


# oak.meandist.10.1500
(sm <- plot(sm(b, 5)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(RIB mean distance (m in 30 min, 1500 m))", y = "s(sqrt(RIB mean distance))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "sqrt.rib.meandist.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.rib.meandist.30.1500^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "RIB mean distance (m in 30 min, 1500 m)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(20,100,500, 1000, 2000, 4000, 8000)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_rib.meandist.30.1500.png", height = 4.5, width = 13)


# group
v <- "group"
plotdata <- plot(b, allTerms = T, select = 7)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +  
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_group.png", height = 4.5, width = 13)


# seastate
v <- "seastate"
plotdata <- plot(b, allTerms = T, select = 8)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    partial.theme + partial.grid.margin +  labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "seastate", var.fac, var.num, var.rand) %>% drop_na(seastate) %>%
    ggplot(aes(x = seastate, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_seastate.png", height = 4.5, width = 13)

# surface.feeding
v <- "surface.feeding"
plotdata <- plot(b, allTerms = T, select = 9)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    partial.theme + partial.grid.margin +  
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface feeding", y = "Partial effect of surface feeding"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "surface.feeding", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.feeding, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Surface feeding", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_surface.feeding.png", height = 4.5, width = 13)

# surface.active
v <- "surface.active"
plotdata <- plot(b, allTerms = T, select = 10)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface activity", y = "Partial effect of surface activity"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "surface.active", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.active, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Surface activity", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_surface.active.png", height = 4.5, width = 13)

