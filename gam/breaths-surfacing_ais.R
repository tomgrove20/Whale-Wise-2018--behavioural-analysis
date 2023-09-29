###################################
### GAM: SURFACING BREATHS, AIS ###
###################################
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

# AIS
tot <- read.csv("intermediate-products/response-var-dfs/surface-interval_ais.csv") %>%
  filter(!is.na(lon.end)) %>%
  mutate_at(vars(contains("dt")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate.last", "year", "folnum.unique", "surface.feeding.follow", "surface.active.follow"), as.factor) %>%
  mutate(dummy = 1) # important dummy variable to factor out random effect when predicting response
# filter(rib.num.30.1500<=4) # only 1/1243 rows had 5 ribs in 30.1500

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 690 intervals, 242 follows


#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = breath.num)) + geom_histogram()
ggplot(data = tot, aes(x = log(breath.num))) + geom_histogram() # so Poisson with log link?

# define explanatory variables, then view distribution
vars <- c("seastate.last", "group","surface.feeding.follow", "dive.time.before", "day.of.year", "oak.meandist.30.1500", "rib.meandist.30.1500", "oak.dist.max.30.1500", "rib.dist.max.30.1500", "oak.num.30.1500", "rib.num.30.1500")
for (i in vars) {
  if(class(tot[,i]) == "numeric"){ # if numeric, provide a histogram
    print(ggplot(data = tot, aes(x = get(i))) + geom_histogram() + xlab(i))
  } else { # if not numeric, provide a bar plot
    print(ggplot(data = tot, aes(x = get(i))) + geom_bar() + xlab(i))
  }
  readline(prompt = "next plot")
}

# dive.time.before
ggplot(data = tot, aes(x = log(dive.time.before))) + geom_histogram()
# rib.meandist
ggplot(data = tot, aes(x = oak.meandist.30.1500)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(rib.meandist.30.1500))) + geom_histogram()

# sqrt transform oak and rib meandist and dist.max
# creating new transformed columns
tot <- tot %>%
  mutate(log.dive.time.before = log(dive.time.before),
         across(contains("meandist"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt mean dist
         across(contains("dist.max"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}")) # sqrt dist max


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate.last", "group", "surface.feeding.follow", "surface.active.follow",  "log.dive.time.before", "year", "day.of.year", "oak.meandist.30.1500", "sqrt.rib.meandist.10.1500", "oak.dist.max.30.1500", "sqrt.rib.dist.max.10.1500", "oak.num.30.1500", "rib.num.30.1500")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/ibi-mean/ais/gam_ibi-mean_collinear.png")

# all problems are solved if we remove dist.max


#---------------- TIME/DISTANCE FRAME --------------------

# first defining explanatory variables (not running folnum.unique as random in this)
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.rand = c("folnum.unique") # random 
var.num = c("day.of.year", "log.dive.time.before") # numeric 
# specifying time/distance frame
ts = c("30.1500", "10.500", "10.1000", "10.1500") # time/distance frame
var.num.t = c("oak.meandist", "sqrt.rib.meandist") # varying continuous numeric vars
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
  
  f1 <- as.formula(paste("breath.num ~", 
                         paste("s(",v.int,", bs = 'cs', k=3)", collapse= "+"),"+", # integer smooths
                         paste("s(",v.num,", bs = 'cs', k=10)", collapse = "+"), "+", # numeric smooths
                         paste(var.fac, collapse = "+")))
  
  # run binomial GAM
  gam <- do.call("gam", list(as.formula(f1), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
  
  # and dev and gcv to df!
  t.df[i,"Deviance explained"] = summary(gam)$dev.expl
  t.df[i,"GCV score"] = gam$gcv.ubre
  
  print (t.df[i,"var.target"]) # progress checker
}

t.df # oak = 30.1500, rib = 10.1000

t.df %>% gather(key = "type", value = "score", 'Deviance explained', 'GCV score') %>%
  ggplot(aes(x = t, y = score, group = var, color = var)) +
  geom_point(size = 2, alpha = 0.5) + geom_line() + facet_wrap(~type, scales = "free_y", ncol = 1) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_discrete(labels=c("10 min, 500 m", "10 min, 1000 m", "10 min, 1500 m", "30 min, 1500 m")) +
  labs(x = "", y = "", color = "Variable") + plot.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("intermediate-products/gam/breaths-surfacing/ais/gam_breaths-surfacing_ais_time-dist_dev+gcv.png", dpi = 600, height = 8, width = 8)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow")
var.int = c("oak.num.30.1500", "rib.num.10.1500")
var.num = c("day.of.year", "oak.meandist.30.1500", "sqrt.rib.meandist.10.1500", "log.dive.time.before")
var.rand = c("folnum.unique")

# formula
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
summary(gam) # dev = 85.3%, gcv = 0.4795


# Remove oak.num, rib.num and rib.meandist
var.num = c("day.of.year", "oak.meandist.30.1500", "log.dive.time.before")
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
summary(gam) # dev = 85.3%, gcv = 0.4795. no change, removing


# Remove surface.feeding.follow
var.fac = c("group", "seastate.last", "year", "surface.active.follow")
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
summary(gam) # dev = 85.3%, gcv = 0.4783. no change, removing. All significant


# final model
var.fac = c("group", "seastate.last", "year", "surface.active.follow")
var.num = c("day.of.year", "oak.meandist.30.1500", "log.dive.time.before")
var.rand = c("folnum.unique")
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = quasipoisson(link = "log")))
summary(gam) # dev = 85.3%, gcv = 0.4785

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam/breaths-surfacing/ais/gam_breaths-surfacing_ais_final.rds")


#---------------- MODEL RESULTS --------------------
#gam <- readRDS("intermediate-products/gam/breaths-surfacing/ais/gam_breaths-surfacing_ais_final.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam/breaths-surfacing/ais/gam_breaths-surfacing_ais_results.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED

# vars
var.fac = c("group", "seastate.last", "year", "surface.active.follow")
var.num = c("day.of.year", "oak.meandist.30.1500", "log.dive.time.before")
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)

# day.of.year
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/ais/gam_breaths-surfacing_ais_julian.png", height = 4.5, width = 13)

# oak.meandist.30.1500
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Oak boat mean distance (m in 30 min, 1500 m)", y = "s(oak boat mean distance)") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "oak.meandist.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = oak.meandist.30.1500, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Oak boat mean distance (m in 30 min, 1500 m)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/ais/gam_breaths-surfacing_ais_oak.meandist.30.1500.png", height=4.5, width=13)


# dive time before
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(dive time before (sec))", y = "s(log(dive time before))") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "log.dive.time.before", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.dive.time.before), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Dive time before (sec)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin + scale_x_continuous(trans = "log", breaks = c(20, 50, 100, 200, 400, 800)) +
    scale_y_continuous(breaks = seq(2,16, by = 2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/ais/gam_breaths-surfacing_ais_dive.time.before.png", height=4.5, width=13)


# group
v <- "group"
plotdata <- plot(b, allTerms = T, select = 5)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "breath.num", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,20,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/ais/gam_breaths-surfacing_ais_group.png", height = 4.5, width = 13)


# seastate.last
v <- "seastate.last"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "breath.num", var = "seastate.last", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = seastate.last, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,20,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/ais/gam_breaths-surfacing_ais_seastate.png", height = 4.5, width = 13)


# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 7)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "breath.num", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,20,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/breaths-surfacing/ais/gam_breaths-surfacing_ais_year.png", height = 4.5, width = 13)



  