#################################
### V2: GAM: SURFACING IBI, FOCAL ###
#################################
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

# conversion rate from m/s to mph
conv <- 2.23694 # for m/s to mph
conv.sq <- 8052.9692102775 # for m/s2 to mph2

# focal vessel
tot <- read.csv("intermediate-products/response-var-dfs/surface-interval_focal.csv") %>%
  filter(!is.na(ibi.mean) & !is.na(dist.mean)) %>%
  mutate_at(vars(contains("dt")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate.last", "year", "folnum.unique", "surface.feeding.follow", "surface.active.follow"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1)  %>% # important dummy variable to factor out random effect when predicting response
  # dplyr::select(-contains(".10")) %>%
  # dplyr::select(-contains(c("0_end", "0_start"))) %>%
  rename_at(vars(contains("_mean")), funs(str_replace(.,"_mean",""))) %>%
  # we also need to convert m/s to mph!
  mutate(across(contains("speed"), ~.*conv), # speed and sd speed
         across(contains("accel"), ~.*conv.sq)) # acceleration

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 66 surfacing intervals, 40 follows

# now we need to add encounter minute to this (rough!!)
follows.encmin <- read.csv("intermediate-products/follows/follows_start+end_consecutive_rough.csv") %>%
  mutate(folnum.unique = as.factor(folnum.unique)) 

tot <- tot %>% left_join(follows.encmin) %>%
  group_by(enc.num) %>%
  mutate(encounter.minute = as.numeric(difftime(dt.start, start, unit = "mins"))) %>% ungroup()


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
ggplot(data = tot, aes(x = asin(vess.di.60))) + geom_histogram()
# speed.var?
ggplot(data = tot, aes(x = sqrt(vess.speed.var.60))) + geom_histogram()
# acceleration
ggplot(data = tot, aes(x = sqrt(vess.accel.max.60))) + geom_histogram()

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("dive.time","dist")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # sqrt mean dist
         across(contains(c("speed.var", "accel.max", "encounter.minute")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt dist max
         across(contains("vess.di."), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))



#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate.last", "group", "surface.feeding.follow", "surface.active.follow", "year", "day.of.year", "log.dive.time.before", "log.dist.end", "arcsin.vess.di.60", "vess.speed.60", "sqrt.vess.speed.var.60", "sqrt.encounter.minute")

# correlation matrix
x <- tot[,vars] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/ibi-mean/foc/gam_ibi-mean_foc_collinear.png", scale = 1.2)

# actually no troubling correlations at all! Let's keep them in!


#---------------- FULL MODE + AUTOCORRELATION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.num = c("day.of.year","log.dive.time.before",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "sqrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("ibi.mean ~", 
                      paste("s(",var.num,", bs = 'cs', k=3)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "log")))
summary(gam) # dev = 93.3%, gcv = 4.771

acf.pacf.plot(acf.structured(gam, 5, "gam"), 5)
# nonsensical to look at ACF plots here because we have 64 points and 39 follows, only 15 have >1 IBI value!

# trying with a bam
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
# for bam we need to log-transform
tot <- tot %>% mutate(log.ibi.mean = log(ibi.mean))
f <- as.formula(paste("log.ibi.mean ~", 
                      paste("s(",var.num,", bs = 'cs', k=3)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam1 <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                            family = gaussian(link = "identity"), rho = -0.99, start = ar.start))

acf.pacf.plot(acf.structured(gam1, 3, "bam"),3)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate.last", "year", "surface.feeding.follow", "surface.active.follow") # factor
var.num = c("day.of.year","log.dive.time.before",  "log.dist.mean", "vess.speed.60",  "sqrt.vess.speed.var.60", "arcsin.vess.di.60", "sqrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("ibi.mean ~", 
                      paste("s(",var.num,", bs = 'cs', k=3)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "log")))
summary(gam); AIC(gam) # dev = 91.3%, -REML = 161.0, AIC = 253

# final model (after manual selection)
var.fac = c("year", "surface.active.follow")
var.num = c("log.dist.mean", "vess.speed.60") 
var.rand = c("folnum.unique")
f <- as.formula(paste("ibi.mean ~", 
                      paste("s(",var.num,", bs = 'cs', k=3)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "log")))
summary(gam); AIC(gam) # dev = 91.6%, -REML = 156.3, AIC = 247.8

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam-v3/ibi-mean/foc/gam_ibi-mean_foc_final_v3.rds")



#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam-v3/ibi-mean/foc/gam_ibi-mean_foc_final_v3.rds")

# defining variable combinations
var.fac = c("year", "surface.active.follow")
var.num = c("log.dist.mean", "vess.speed.60") 
var.rand = c("folnum.unique")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/ibi-mean/foc/gam_ibi-mean_foc_results_v3.csv", row.names = FALSE)

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
    labs(x = "Mean distance (m)", y = "Mean IBI (sec)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(50,100,200,500,1000)) +
    scale_y_continuous(breaks = seq(0,40,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/ibi-mean/foc/gam_ibi-mean_foc_dist.mean_v3.png", height = 4.5, width = 13)

# vess.speed.60
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Mean vessel speed (mph)", y = "s(mean vessel speed)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "ibi.mean", var = "vess.speed.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = vess.speed.60, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Mean vessel speed (mph)", y = "Mean IBI (sec)") + partial.theme + partial.grid.margin +
    scale_y_continuous(breaks = seq(0,40,2)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/ibi-mean/foc/gam_ibi-mean_foc_vess.speed.60_v3.png", height = 4.5, width = 13)

# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 4)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "ibi.mean", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Mean IBI (sec)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/ibi-mean/foc/gam_ibi-mean_foc_year_v3.png", height = 4.5, width = 13)

# surface.active.follow
v <- "surface.active.follow"
plotdata <- plot(b, allTerms = T, select = 5)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface activity in follow", y = "Partial effect of surface activity"))
(pr <- pred.var.fac(gam, tot, resp = "ibi.mean", var = "surface.active.follow", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.active.follow, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Surface activity in follow", y = "Mean IBI (sec)") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/ibi-mean/foc/gam_ibi-mean_foc_surface.active.follow_v3.png", height = 4.5, width = 13)
