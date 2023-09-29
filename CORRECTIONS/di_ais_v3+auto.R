####################
### V2: GAM: DI, AIS ###
####################
### 11/01/2023
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
tot <- read.csv("intermediate-products/response-var-dfs/di_ais.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding", "surface.active"), as.factor) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(difftime.next<2000) # until properly processed

# quick look at the data set
colSums(is.na(tot)) 
nrow(tot); length(unique(tot$folnum.unique)) # 3555 surfacing pairs, 452 follows


#---------------- TRANSFORMATIONS --------------------

# response variable
ggplot(data = tot, aes(x = DI)) + geom_histogram()
ggplot(data = tot, aes(x = asin(DI))) + geom_histogram() # going to arcsin transform this simply to reduce skew

# define explanatory variables, then view distribution
vars <- c("seastate", "group","surface.feeding", "surface.active", "day.of.year", "difftime.next", "difftime.prev", "oak.meandist.30.1500", "rib.meandist.30.1500", "oak.dist.max.30.1500", "rib.dist.max.30.1500", "oak.num.30.1500", "rib.num.30.1500")
for (i in vars) {
  if(class(tot[,i]) == "numeric"){ # if numeric, provide a histogram
    print(ggplot(data = tot, aes(x = get(i))) + geom_histogram() + xlab(i))
  } else { # if not numeric, provide a bar plot
    print(ggplot(data = tot, aes(x = get(i))) + geom_bar() + xlab(i))
  }
  readline(prompt = "next plot")
}

# difftime.next
# difftime.prev
# oak.meandist
ggplot(data = tot, aes(x = oak.meandist.30.1500)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(oak.meandist.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = rib.meandist.30.1500)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(rib.meandist.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = oak.dist.max.30.1500)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(oak.dist.max.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(rib.dist.max.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = log(difftime.prev))) + geom_histogram()
ggplot(data = tot, aes(x = log(difftime.next))) + geom_histogram()

# creating new transformed columns
tot <- tot %>%
  mutate(arcsin.DI = asin(DI),
         across(contains(c("meandist", "dist.max")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt
         across(contains(c("difftime")), .fns = list(log = ~log(.)), .names = "{fn}.{col}")) # sqrt dist max


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group", "surface.feeding", "surface.active", "log.difftime.prev", "log.difftime.next", "day.of.year", "oak.num.30.1500", "rib.num.30.1500","sqrt.oak.meandist.30.1500", "rib.meandist.30.1500", "sqrt.oak.dist.max.30.1500", "sqrt.rib.dist.max.30.1500",  "year")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 2.5) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/di/ais/gam_di_ais_collinear.png")

# max.dist causing most problems
# rib.num ~ rib.meandist = 0.72 but not too bad


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "log.difftime.prev", "log.difftime.next", "sqrt.oak.meandist.30.1500", "rib.meandist.30.1500")
var.rand = c("folnum.unique")


# formula
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "identity")))
summary(gam) # dev = 22.4%, gcv = 0.073

acf.pacf.plot(acf.structured(gam, 10, "gam"), 10)
ggsave("intermediate-products/gam-v3/di/ais/gam_di_ais_acf.png", width = 11, height = 3.6)

# trying bam
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
gam1 <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", family = gaussian(link = "identity"),
                            rho = 0.18, start = ar.start))
acf.pacf.plot(acf.structured(gam1, 10, "bam"), 10)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active")
var.int = c("oak.num.30.1500", "rib.num.30.1500")
var.num = c("day.of.year", "log.difftime.prev", "log.difftime.next", "sqrt.oak.meandist.30.1500", "rib.meandist.30.1500")
var.rand = c("folnum.unique")


# formula
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                            family = gaussian(link = "identity"), rho = 0.18, start = ar.start))
summary(gam); AIC(gam) # dev = 16.8%, -REML = 391.4, AIC = 696.1

# final model (after manual selection)
var.fac = c("group", "seastate", "surface.feeding", "surface.active")
var.int = c("rib.num.30.1500")
var.num = c("log.difftime.prev")
var.rand = c("folnum.unique")
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "REML", 
                           family = gaussian(link = "identity"), rho = 0.18, start = ar.start))
summary(gam); AIC(gam) # dev = 16.8%, -REML = 386.8, AIC = 695.5

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam-v3/di/ais/gam_di_ais_final_v3.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam-v3/di/ais/gam_di_ais_final_v3.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/di/ais/gam_di_ais_results_v3.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
var.fac = c("group", "seastate", "surface.feeding", "surface.active")
var.num = c("rib.num.30.1500", "log.difftime.prev")
var.rand = c("folnum.unique")

print(plot(b, allTerms = T), pages = 1)

# rib.num.30.1500
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "RIB boats", y = "s(RIB boats)") + partial.theme + partial.grid.margin + scale_x_continuous(breaks = c(0:10)))
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "rib.num.30.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = rib.num.30.1500, y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "RIB boats", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.9, 0.95, 0.96, 0.97, 0.98, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/ais/gam_di_ais_rib.num.30.1500_v3.png", height = 4.5, width = 13)


# difftime.prev
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "log(previous IBI (sec))", y = "s(log(previous IBI))") + partial.theme +  partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "arcsin.DI", var = "log.difftime.prev", var.fac, var.num, var.rand) %>%
    mutate(upr = ifelse(exp(log.difftime.prev)>550, asin(1), upr)) %>%
    ggplot(aes(x = exp(log.difftime.prev), y = sin(predict))) + geom_line() +
    geom_ribbon(aes(ymin = sin(lwr), ymax = sin(upr)), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Previous IBI (sec)", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.93, 0.95, 0.97, 0.98, 0.99,0.99)) +
    scale_x_continuous(trans = "log", breaks = c(5, 10, 20, 50, 100, 200, 400, 800)) +
    coord_cartesian(ylim = c(0.92, 0.998)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/ais/gam_di_ais_difftime.prev_v3.png", height = 4.5, width = 13)

# group
v <- "group"
plotdata <- plot(b, allTerms = T, select = 4)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +  
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "group", var.fac, var.num, var.rand) %>% drop_na(group) %>%
    ggplot(aes(x = group, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.9, 0.95, 0.98, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/ais/gam_di_ais_group_v3.png", height = 4.5, width = 13)

# seastate
v <- "seastate"
plotdata <- plot(b, allTerms = T, select = 5)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +  
    labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "seastate", var.fac, var.num, var.rand) %>% drop_na(seastate) %>%
    ggplot(aes(x = seastate, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_y_continuous(trans = asin.trans, breaks = c(0.9, 0.95, 0.98, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/ais/gam_di_ais_seastate_v3.png", height = 4.5, width = 13)

# surface.feeding
v <- "surface.feeding"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) + 
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface feeding", y = "Partial effect of surface feeding"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "surface.feeding", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.feeding, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Surface feeding", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")) +
    scale_y_continuous(trans = asin.trans, breaks = c(0.9, 0.95, 0.98, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/ais/gam_di_ais_surface.feeding_v3.png", height = 4.5, width = 13)

# surface active
v <- "surface.active"
plotdata <- plot(b, allTerms = T, select = 7)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
(sm <- eval(var.fac.plot) +
    scale_x_discrete(labels = c("no", "yes")) +
    labs(x = "Surface activity", y = "Partial effect of surface activity"))
(pr <- pred.var.fac(gam, tot, resp = "arcsin.DI", var = "surface.active", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = surface.active, y = sin(predict))) + geom_boxplot(width = 0.4) +
    labs(x = "Surface activity", y = "Directness index") + partial.theme + partial.grid.margin +
    scale_x_discrete(labels = c("no", "yes")) +
    scale_y_continuous(trans = asin.trans, breaks = c(0.9, 0.95, 0.98, 0.99)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/di/ais/gam_di_ais_surface.active_v3.png", height = 4.5, width = 13)



