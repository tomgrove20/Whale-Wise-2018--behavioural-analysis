###############################
### V2: GAM: SURFACING IBI, AIS ###
###############################
### 11/01/2023
### Tom Grove
### tomgrove20@yahoo.co.uk

# 1. transformations
# 2. Collinearity
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
  filter(!is.na(resp.rate) & !is.na(lon.end)) %>%
  mutate_at(vars(contains("dt")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate.last", "year", "folnum.unique", "surface.feeding.follow", "surface.active.follow"), as.factor) %>%
  mutate(dummy = 1) # important dummy variable to factor out random effect when predicting response
  # filter(rib.num.30.1500<=4) # only 1/1243 rows had 5 ribs in 30.1500

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 359 intervals, 189 follows


#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = ibi.mean)) + geom_histogram()
ggplot(data = tot, aes(x = log(ibi.mean))) + geom_histogram()

# define explanatory variables, then view distribution
vars <- c("seastate.last", "group","surface.feeding.follow",  "day.of.year", "oak.meandist.30.1500", "rib.meandist.30.1500", "oak.dist.max.30.1500", "rib.dist.max.30.1500", "oak.num.30.1500", "rib.num.30.1500")
for (i in vars) {
  if(class(tot[,i]) == "numeric"){ # if numeric, provide a histogram
    print(ggplot(data = tot, aes(x = get(i))) + geom_histogram() + xlab(i))
  } else { # if not numeric, provide a bar plot
    print(ggplot(data = tot, aes(x = get(i))) + geom_bar() + xlab(i))
  }
  readline(prompt = "next plot")
}

# sqrt transform oak and rib meandist and dist.max
# creating new transformed columns
tot <- tot %>%
  mutate(log.dive.time.before = log(dive.time.before),
         across(contains("meandist"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt mean dist
         across(contains("dist.max"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}")) # sqrt dist max


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate.last", "group", "surface.feeding.follow", "surface.active.follow", "log.dive.time.before", "year", "day.of.year", "oak.meandist.30.1500", "sqrt.rib.meandist.30.1500", "sqrt.oak.dist.max.30.1500", "sqrt.rib.dist.max.10.1500", "oak.num.30.1500", "rib.num.30.1500")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/ibi-mean/ais/gam_ibi-mean_collinear.png")

# most problems are solved if we remove dist.max


#---------------- FULL MODEL + AUTOCORRELATION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "oak.meandist.30.1500", "sqrt.rib.meandist.30.1500", "log.dive.time.before")
var.rand = c("folnum.unique")

# formula
f <- as.formula(paste("ibi.mean ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "log")))
summary(gam) # dev = 71.9%, gcv = 8.44. Remove sea state

# View ACF plots
acf.pacf.plot(acf.structured(gam, 5, "gam"), 5)
ggsave("intermediate-products/gam-v3/ibi-mean/ais/gam_ibi-mean_ais_acf.png", width = 11, height = 3.6)

# trying with a bam
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
# for bam we need to log-transform
tot <- tot %>% mutate(log.ibi.mean = log(ibi.mean))
f <- as.formula(paste("log.ibi.mean ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam1 <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "identity"),
                            rho = 0.4, start = ar.start))

acf.pacf.plot(acf.structured(gam1, 5, "gam"), 5)
# rho = 0.2 is perfect!!



#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "oak.meandist.30.1500", "sqrt.rib.meandist.30.1500", "log.dive.time.before")
var.rand = c("folnum.unique")

# formula
f <- as.formula(paste("log.ibi.mean ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = 0.4, start = ar.start))
summary(gam); AIC(gam) # dev = 56.8%, -REML = -2.00, AIC = -69.00

# final model (after manual selection)
var.fac = c("surface.feeding.follow")
var.int = c("oak.num.30.1500")
var.num = c("day.of.year")
var.rand = c("folnum.unique")
f <- as.formula(paste("log.ibi.mean ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = 0.4, start = ar.start))
summary(gam); AIC(gam) # dev = 57.1%, -REML = -11.3, AIC = -71.7

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam-v3/ibi-mean/ais/gam_ibi-mean_ais_final_v3.rds")


#---------------- MODEL RESULTS --------------------

# gam <- readRDS("intermediate-products/gam-v3/ibi-mean/ais/gam_ibi-mean_ais_final_v3.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/ibi-mean/ais/gam_ibi-mean_ais_results_v3.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
var.fac = c("surface.feeding.follow")
var.num = c("oak.num.30.1500", "day.of.year")
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)

# oak.num.30.1500
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Oak boats", y = "s(oak boats)") + partial.theme + partial.grid.margin + scale_x_continuous(breaks = c(0:10)))
(pr <- pred.var(gam, tot, resp = "log.ibi.mean", var = "oak.num.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = oak.num.30.1500, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Oak boats", y = "Mean IBI (sec)") + partial.theme + partial.grid.margin +
    scale_x_continuous(breaks = c(0:10)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/ibi-mean/ais/gam_ibi-mean_ais_oak.num.30.1500_v3.png", height = 4.5, width = 13)

# day.of.year
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "log.ibi.mean", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Mean IBI (sec)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/ibi-mean/ais/gam_ibi-mean_ais_julian_v3.png", height = 4.5, width = 13)


# surface.feeding.follow
v <- "surface.feeding.follow"
plotdata <- plot(b, allTerms = T, select = 4)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface feeding in follow", y = "Partial effect of surface feeding"))
(pr <- pred.var.fac(gam, tot, resp = "log.ibi.mean", var = "surface.feeding.follow", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.feeding.follow, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Surface feeding in follow", y = "Mean IBI (sec)") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/ibi-mean/ais/gam_ibi-mean_ais_surface.feeding.follow_v3.png", height = 4.5, width = 13)

