###########################
### GAM: DIVE TIME, AIS ###
###########################
### 16/02/2022
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

# AIS
tot <- read.csv("intermediate-products/response-var-dfs/dive-time_ais.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding.follow", "surface.active.follow"), as.factor) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(rib.num.30.1500<=4) # only 1/1243 rows had 5 ribs in 30.1500

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 1243 dives from 387 follows


#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = dive.time)) + geom_histogram()
ggplot(data = tot, aes(x = log(dive.time))) + geom_histogram() # need log link!

# define explanatory variables, then view distribution
vars <- c("seastate", "group","surface.feeding.follow",  "day.of.year", "oak.meandist.30.1500", "rib.meandist.30.1500", "oak.dist.max.30.1500", "rib.dist.max.30.1500", "oak.num.30.1500", "rib.num.30.1500")
for (i in vars) {
  if(class(tot[,i]) == "numeric"){ # if numeric, provide a histogram
    print(ggplot(data = tot, aes(x = get(i))) + geom_histogram() + xlab(i))
  } else { # if not numeric, provide a bar plot
    print(ggplot(data = tot, aes(x = get(i))) + geom_bar() + xlab(i))
  }
  readline(prompt = "next plot")
}

# not transforming rib and oak distances. Log transformed won't work for rib (presence of 0s), want oak and rib to be comparable!
ggplot(data = tot, aes(x = oak.meandist.30.1500)) + geom_histogram()
ggplot(data = tot, aes(x = rib.meandist.30.1500)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(oak.meandist.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(rib.meandist.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(oak.dist.max.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = rib.dist.max.30.1500)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(rib.dist.max.30.1500))) + geom_histogram()

# not transforming any columns actually!


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group", "surface.feeding.follow", "surface.active.follow", "year", "day.of.year", "oak.meandist.30.1500", "rib.meandist.30.1500", "oak.dist.max.30.1500", "rib.dist.max.30.1500", "oak.num.30.1500", "rib.num.30.1500")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/dive-time/ais/gam_dive-time_collinear.png")

# all problems are solved if we remove dist.max. We'd rather go with meandist anyway


#---------------- TIME/DISTANCE FRAME --------------------

# first defining explanatory variables (not running folnum.unique as random in this)
var.fac = c("group", "seastate", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.rand = c("folnum.unique") # random 
var.num = c("day.of.year") # numeric 
# specifying time/distance frame
ts = c("30.1500", "10.500", "10.1000", "10.1500") # time/distance frame
var.num.t = c("oak.meandist", "rib.meandist") # varying continuous numeric vars
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
  
  f1 <- as.formula(paste("dive.time ~", 
                         paste("s(",v.int,", bs = 'cs', k=3)", collapse= "+"),"+", # integer smooths
                         paste("s(",v.num,", bs = 'cs', k=10)", collapse = "+"), "+", # numeric smooths
                         paste(var.fac, collapse = "+")))
  
  # run binomial GAM
  gam <- do.call("gam", list(as.formula(f1), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
  
  # and add AIC to df!
  t.df[i,"AIC"] = AIC(gam)
  
  print (t.df[i,"var.target"]) # progress checker
}

t.df # num = 30.1500, dist = 10.1500

ggplot(data = t.df, aes(x = t, y = AIC, group = var, color = var)) +
  geom_point(size = 2, alpha = 0.5) + geom_line() + 
  scale_color_brewer(palette = "Dark2") +
  scale_x_discrete(labels=c("10 min, 500 m", "10 min, 1000 m", "10 min, 1500 m", "30 min, 1500 m")) +
  labs(x = "", color = "Variable") + plot.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("intermediate-products/gam/dive-time/ais/gam_dive-time_ais_time-dist_aic.png", dpi = 600, height = 5, width = 8)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding.follow", "surface.active.follow")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "oak.meandist.10.1500", "rib.meandist.10.1500")
var.rand = c("folnum.unique")

# formula
f <- as.formula(paste("dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 83%, gcv = 6755.7

# Remove group
var.fac = c("year", "seastate", "surface.feeding.follow", "surface.active.follow")
f <- as.formula(paste("dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 83%, gcv = 6747. no change and good, so overall remove

# remove seastate
var.fac = c("year", "surface.feeding.follow", "surface.active.follow")
f <- as.formula(paste("dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 83%, gcv = 6748. no change and barely change, so overall remove

# Remove oak.meandist.10.1500
var.num = c("day.of.year", "rib.meandist.10.1500")
f <- as.formula(paste("dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 82.9%, gcv = 6744. tiny bit bad and quite good, so removing

# Remove surface.active.follow
var.fac = c("year", "surface.feeding.follow")
f <- as.formula(paste("dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 82.9%, gcv = 6743. no change and good. removing

# Remove rib.num.30.1500
var.int = c("oak.num.30.1500")
f <- as.formula(paste("dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 82.8%, gcv = 6767. bad and bad, removing
var.int = c("oak.num.30.1500", "rib.num.30.1500")


# final model
var.fac = c("year", "surface.feeding.follow")
var.num = c("day.of.year", "rib.meandist.10.1500")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.rand = c("folnum.unique")

f <- as.formula(paste("dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 5)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = gaussian(link = "log")))
summary(gam) # dev = 82.8%, gcv = 6767

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam/dive-time/ais/gam_dive-time_ais_final.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam/dive-time/ais/gam_dive-time_ais_final.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam/dive-time/ais/gam_dive-time_ais_results.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
var.fac = c("year", "surface.feeding.follow")
var.num = c("day.of.year", "rib.meandist.10.1500", "oak.num.30.1500", "rib.num.30.1500")
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)

# oak.num.30.1500
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Oak boats (30 min, 1500 m)", y = "s(oak boats)") + partial.theme + partial.grid.margin + scale_x_continuous(breaks = c(0:10)))
(pr <- pred.var(gam, tot, resp = "dive.time", var = "oak.num.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = oak.num.30.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Oak boats (30 min, 1500 m)", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    scale_x_continuous(breaks = c(0:10)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/dive-time/ais/gam_dive-time_ais_oak.num.30.1500.png", height = 4.5, width = 13)


# rib.num.30.1500
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "RIB boats (30 min, 1500 m)", y = "s(RIB boats)") + partial.theme + partial.grid.margin + scale_x_continuous(breaks = c(0:10)))
(pr <- pred.var(gam, tot, resp = "dive.time", var = "rib.num.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = rib.num.30.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "RIB boats (30 min, 1500 m)", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    scale_x_continuous(breaks = c(0:10)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/dive-time/ais/gam_dive-time_ais_rib.num.30.1500.png", height = 4.5, width = 13)


# day.of.year
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "dive.time", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Dive time (sec)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/dive-time/ais/gam_dive-time_ais_julian.png", height = 4.5, width = 13)

# rib.meandist.10.1500
(sm <- plot(sm(b, 4)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "RIB mean distance (m in 10 min, 1500 m)", y = "s(RIB mean distance)") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "dive.time", var = "rib.meandist.10.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = rib.meandist.10.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "RIB mean distance (m in 10 min, 1500 m)", y = "Dive time (sec)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/dive-time/ais/gam_dive-time_ais_rib.meandist.10.1500.png", height=4.5, width=13)

# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "dive.time", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Dive time (sec)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/dive-time/ais/gam_dive-time_ais_year.png", height = 4.5, width = 13)

# surface.feeding.follow
v <- "surface.feeding.follow"
plotdata <- plot(b, allTerms = T, select = 7)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface feeding in follow", y = "Partial effect of surface feeding"))
(pr <- pred.var.fac(gam, tot, resp = "dive.time", var = "surface.feeding.follow", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.feeding.follow, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Surface feeding in follow", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/dive-time/ais/gam_dive-time_ais_surface.feeding.follow.png", height = 4.5, width = 13)

