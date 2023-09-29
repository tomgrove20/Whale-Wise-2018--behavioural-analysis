###################################
### V2: GAM: SURFACING BREATHS, AIS ###
###################################
### 11/01/2023
### Tom Grove
### tomgrove20@yahoo.co.uk

# 1. transformations
# 2. Collinearity
# 3. time frame selection
# 4. Final GAM (GCV + backward selection)
# 5. Results + plotting


### PACKAGES
# installing edited ungeviz
#devtools::install("code/ungeviz_tjg")
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


#---------------- FULL MODEL + AUTOCORRELATION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "oak.meandist.30.1500", "rib.meandist.30.1500", "log.dive.time.before")
var.rand = c("folnum.unique")

# formula
f <- as.formula(paste("breath.num ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = quasipoisson(link = "log")))
summary(gam) # dev = 85.4%, gcv = 0.471

# View ACF plots
acf.pacf.plot(acf.structured(gam, 5, "gam"), 5)
ggsave("intermediate-products/gam-v3/breaths-surfacing/ais/gam_breath-num_ais_acf.png", width = 11, height = 3.6)

# trying with a bam
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
# for bam we need to log-transform
tot <- tot %>% mutate(log.breath.num = log(breath.num))
f <- as.formula(paste("log.breath.num ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam1 <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "identity"),
                            rho = 0.2, start = ar.start))

acf.pacf.plot(acf.structured(gam1, 5, "gam"), 5)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "oak.meandist.30.1500", "rib.meandist.30.1500", "log.dive.time.before")
var.rand = c("folnum.unique")

# formula
f <- as.formula(paste("log.breath.num ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                            family = gaussian(link = "identity"), rho = 0.2, start = ar.start))
summary(gam); AIC(gam) # dev = 72%, -REML = 480.42, AIC = 926.3

# final model after manual selection
var.fac = c("group", "seastate.last", "year", "surface.active.follow")
var.num = c("day.of.year", "log.dive.time.before")
var.rand = c("folnum.unique")
f <- as.formula(paste("log.breath.num ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = 0.2, start = ar.start))
summary(gam); AIC(gam) # dev = 72.1%, -REML = 479.2, AIC = 922.3

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam-v3/breaths-surfacing/ais/gam_breaths-surfacing_ais_final_v3.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam-v3/breaths-surfacing/ais/gam_breaths-surfacing_ais_final_v3.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/breaths-surfacing/ais/gam_breaths-surfacing_ais_results_v3.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED

# vars
var.fac = c("group", "seastate.last", "year", "surface.active.follow")
var.num = c("day.of.year", "log.dive.time.before")
var.rand = c("folnum.unique")
print(plot(b, allTerms = T), pages = 1)

# day.of.year
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/ais/gam_breaths-surfacing_ais_julian_v3.png", height = 4.5, width = 13)


# dive time before
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(dive time before (sec))", y = "s(log(dive time before))") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "breath.num", var = "log.dive.time.before", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.dive.time.before), y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Dive time before (sec)", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin + scale_x_continuous(trans = "log", breaks = c(20, 50, 100, 200, 400, 800)) +
    scale_y_continuous(breaks = seq(2,16, by = 2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/ais/gam_breaths-surfacing_ais_dive.time.before_v3.png", height=4.5, width=13)


# group
v <- "group"
plotdata <- plot(b, allTerms = T, select = 4)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "log.breath.num", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,20,2)) +
    coord_cartesian(ylim = c(1,7.7)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/ais/gam_breaths-surfacing_ais_group_v3.png", comb, height = 4.5, width = 13)


# seastate.last
v <- "seastate.last"
plotdata <- plot(b, allTerms = T, select = 5)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "breath.num", var = "seastate.last", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = seastate.last, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,20,2))+
    coord_cartesian(ylim = c(1,7.7)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/ais/gam_breaths-surfacing_ais_seastate_v3.png", height = 4.5, width = 13)


# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "breath.num", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,20,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/ais/gam_breaths-surfacing_ais_year_v3.png", height = 4.5, width = 13)

# surface activity
v <- "surface.active.follow"
plotdata <- plot(b, allTerms = T, select = 7)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface activity in follow", y = "Partial effect of surface activity"))
(pr <- pred.var.fac(gam, tot, resp = "log.breath.num", var = "surface.active.follow", var.fac, var.num, var.rand) %>%
    ggplot(aes(x = surface.active.follow, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Surface activity in follow", y = "Breaths per surfacing interval") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")) +
    coord_cartesian(ylim = c(1,10)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/breaths-surfacing/ais/gam_breaths-surfacing_ais_surface.active.follow_v3.png", height = 4.5, width = 13)
