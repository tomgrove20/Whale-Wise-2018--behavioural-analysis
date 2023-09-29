############################
### V2: GAM: SWIM SPEED, AIS ###
############################
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
vars <- c("seastate", "group", "surface.feeding", "surface.active", "diff.time", "day.of.year", "oak.num.30.1500", "rib.num.30.1500","sqrt.oak.meandist.30.1500", "sqrt.rib.meandist.30.1500", "year")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 2.5) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_collinear.png")

# max.dist causing most problems
# rib.num ~ rib.meandist = 0.72 but not too bad


#---------------- FULL MODEL+ AUTOCORRELATION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "log.diff.time", "sqrt.oak.meandist.30.1500", "rib.meandist.30.1500")
var.rand = c("folnum.unique")


# formula
f <- as.formula(paste("speed ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "log")))
summary(gam) # dev = 30.4%, gcv = 4.767

# View ACF plots
acf.pacf.plot(acf.structured(gam, 10, "gam"), 10)
ggsave("intermediate-products/gam-v3/swim-speed/ais/gam_swim-speed_ais_acf.png", width = 11, height = 3.6)


# trying with a bam
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
# for bam we need to log-transform
tot <- tot %>% mutate(log.speed = log(speed))
f <- as.formula(paste("log.speed ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam1 <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "identity"),
                            rho = -0.05, start = ar.start))
# HOW TO GET THE CORRECT RESIDUALS FROM A BAM? bam$rsd or residuals(bam)????
acf.pacf.plot(acf.structured(gam1, 10, "bam"), 10)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "log.diff.time", "sqrt.oak.meandist.30.1500", "rib.meandist.30.1500")
var.rand = c("folnum.unique")


# formula
f <- as.formula(paste("log.speed ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = -0.05, start = ar.start))
summary(gam) # dev = 26.5%, REML = 2986.1

# final model (after manual selection)
var.num = c("day.of.year", "log.diff.time")
var.fac = c("group", "year", "surface.feeding")
f <- as.formula(paste("log.speed ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = -0.05, start = ar.start))
summary(gam); AIC(gam) # dev = 26.5%, REML = 2984.1, AIC = 5769.4

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam-v3/swim-speed/ais/gam_swim-speed_ais_final_v3.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam-v2/swim-speed/ais/gam_swim-speed_ais_final_v2.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/swim-speed/ais/gam_swim-speed_ais_results_v3.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
var.num = c("day.of.year", "log.diff.time")
var.fac = c("group", "year", "surface.feeding")
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)


# day.of.year
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "log.speed", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/swim-speed/ais/gam_swim-speed_ais_julian_v3.png", height = 4.5, width = 13)


# log.diff.time
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(inter-breath interval (sec))", y = "s(inter-breath interval)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "speed", var = "log.diff.time", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = exp(log.diff.time), y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Inter-breath interval (sec)", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "log", breaks = c(5,10,20,50,100,200,400,800)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/swim-speed/ais/gam_swim-speed_ais_ibi_v3.png", height = 4.5, width = 13)

# group
v <- "group"
plotdata <- plot(b, allTerms = T, select = 4)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +  
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/swim-speed/ais/gam_swim-speed_ais_group_v3.png", height = 4.5, width = 13)


# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 5)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "year", var.fac, var.num, var.rand) %>% drop_na(year) %>%
    ggplot(aes(x = year, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/swim-speed/ais/gam_swim-speed_ais_year_v3.png", height = 4.5, width = 13)

# surface.feeding
v <- "surface.feeding"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) + 
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface feeding", y = "Partial effect of surface feeding"))
(pr <- pred.var.fac(gam, tot, resp = "speed", var = "surface.feeding", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.feeding, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Surface feeding", y = "Swim speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/swim-speed/ais/gam_swim-speed_ais_surface.feeding_v3.png", height = 4.5, width = 13)
