
### V3: GAM: DIVE TIME, AIS 

### 29/05/2023
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


#---------------- FULL MODEL + AUTOCORRELATION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding.follow", "surface.active.follow")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "oak.meandist.30.1500", "rib.meandist.30.1500")
var.rand = c("folnum.unique")

# formula
f <- as.formula(paste("dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "log")))
summary(gam) # dev = 82.2%, gcv = 6991.6

# 313/390 follows have 5 or fewer observations (clear dives), so keeping acf plot to 5
# View ACF plots
acf.pacf.plot(acf.structured(gam, 5, "gam"), 5)
ggsave("intermediate-products/gam-v3/dive-time/ais/gam_dive-time_ais_acf.png", width = 11, height = 3.6)

# trying with a bam
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
# for bam we need to log-transform
tot <- tot %>% mutate(log.dive.time = log(dive.time))
f <- as.formula(paste("log.dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam1 <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "identity"),
                            rho = -0.1, start = ar.start))

acf.pacf.plot(acf.structured(gam1, 4, "bam"), 4)



#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding.follow", "surface.active.follow")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "oak.meandist.30.1500", "rib.meandist.30.1500")
var.rand = c("folnum.unique")

# formula
f <- as.formula(paste("log.dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "identity"),
                            rho = -0.1, start = ar.start))
summary(gam); AIC(gam) # dev = 76.8%, -REML = 1188, AIC = 1984.7

# final model (after manual removal of terms)
var.fac = c("year", "surface.feeding.follow")
var.num = c("day.of.year")
var.int = c("oak.num.30.1500")
var.rand = c("folnum.unique")
f <- as.formula(paste("log.dive.time ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "identity"),
                           rho = -0.1, start = ar.start))
summary(gam); AIC(gam) # dev = 76.8%, =REML = 1184.7, AIC = 1979

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam-v3/dive-time/ais/gam_dive-time_ais_final_v3.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam-v3/dive-time/ais/gam_dive-time_ais_final_v3.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/dive-time/ais/gam_dive-time_ais_results_v3.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
var.fac = c("year", "surface.feeding.follow")
var.num = c("day.of.year", "oak.num.30.1500")
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)

# oak.num.30.1500
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Oak boats", y = "s(oak boats)") + partial.theme + partial.grid.margin + scale_x_continuous(breaks = c(0:10)))
(pr <- pred.var(gam, tot, resp = "dive.time", var = "oak.num.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = oak.num.30.1500, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Oak boats", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    scale_x_continuous(breaks = c(0:10)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/ais/gam_dive-time_ais_oak.num.30.1500_v3.png", height = 4.5, width = 13)

# day.of.year
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "dive.time", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = exp(predict))) + geom_line() +
    geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Dive time (sec)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/ais/gam_dive-time_ais_julian_v3.png", height = 4.5, width = 13)


# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 4)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "dive.time", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Dive time (sec)") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/ais/gam_dive-time_ais_year_v3.png", height = 4.5, width = 13)

# surface.feeding.follow
v <- "surface.feeding.follow"
plotdata <- plot(b, allTerms = T, select = 5)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface feeding in follow", y = "Partial effect of surface feeding"))
(pr <- pred.var.fac(gam, tot, resp = "dive.time", var = "surface.feeding.follow", var.fac, var.num, var.rand) %>%
    ggplot(aes(x = surface.feeding.follow, y = exp(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Surface feeding in follow", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/dive-time/ais/gam_dive-time_ais_surface.feeding.follow_v3.png", height = 4.5, width = 13)

