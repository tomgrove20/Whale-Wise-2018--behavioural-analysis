### V3: GAM: SURFACE FEEDING, FOCAL

### 24/05/2023
### Tom Grove
### tomgrove20@yahoo.co.uk

# 1. transformations
# 2. Collinearity
# 4. Final GAM (GCV + backward selection)
# 5. Results + plotting

# note: PQL is used in gamm or gamm4, which performs poorly with binary data. So we're sticking with gam for surface activity and surface feeding


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
tot <- read.csv("intermediate-products/response-var-dfs/surface-active+feeding_focal.csv") %>%
  dplyr::select(-c("lat.whale", "lon.whale")) %>% drop_na() %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  filter(dist<1000) %>% # more likely to see surfacings that are breaches at longer distances
  filter(vess.accel.min.300>-10) %>% # temporary filtering until boat stuff properly processes
  filter(vess.accel.max.300<0.5) %>% # filtering out until boat stuff properly processes
  # we also need to convert m/s to mph!
  mutate(across(contains("speed"), ~.*conv), # speed and sd speed
         across(contains("accel"), ~.*conv.sq)) # acceleration

# quick look at the data set
colSums(is.na(tot)) 

# now we need to add encounter minute to this (rough!!)
follows.encmin <- read.csv("intermediate-products/follows/follows_start+end_consecutive_rough.csv") %>%
  mutate(folnum.unique = as.factor(folnum.unique))

tot <- tot %>% left_join(follows.encmin) %>%
  group_by(enc.num) %>%
  mutate(encounter.minute = as.numeric(difftime(datetime, start, unit = "mins"))) %>% ungroup()


#---------------- TRANSFORMATIONS --------------------

# define explanatory variables, then view distribution
vars <- c("dist", "seastate", "group", "day.of.year", "vess.di.60", "vess.accel.60", "vess.speed.60", "vess.speed.var.60", "encounter.minute")
for (i in vars) {
  if(class(tot[,i]) == "numeric"){ # if numeric, provide a histogram
    print(ggplot(data = tot, aes(x = get(i))) + geom_histogram() + xlab(i))
  } else { # if not numeric, provide a bar plot
    print(ggplot(data = tot, aes(x = get(i))) + geom_bar() + xlab(i))
  }
  readline(prompt = "next plot")
}

# accel  max not ideal but difficult to transform since crossing 0. DI is left-skewed and speed var is right skewed

# speed.var: log transformation
sort(tot$vess.speed.var.60); ggplot(data = tot, aes(x = log(vess.speed.var.60))) + geom_histogram()
# di: arcsin transformation
ggplot(data = tot, aes(x = asin(vess.di.60))) + geom_histogram()
# dist: log transformation
ggplot(data = tot, aes(x = log(dist))) + geom_histogram()

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("speed.var", "dist")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # log speed var
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"), # arcsin DI
         across(contains(c("encounter.minute")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         log.dist = log(dist))


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group", "year", "day.of.year", "log.dist", "arcsin.vess.di.60",  "vess.speed.60", "log.vess.speed.var.60","sqrt.encounter.minute")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/surface-feeding/foc/gam_surface-feeding_foc_collinear.png", scale = 1.2)

# as expected, correlated accel and accel.max
# going for accel max for now


#---------------- FULL MODEL + AUTOCORRELATION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year") # factor
var.num = c("day.of.year", "log.dist", "vess.speed.60",  "log.vess.speed.var.60", "arcsin.vess.di.60", "sqrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "REML", family =binomial))
summary(gam) # dev = 65.1%, ubre = -0.78

# View ACF plots
acf.pacf.plot(acf.structured(gam, 10, "gam"), 10)
ggsave("intermediate-products/gam-v3/surface-feeding/foc/gam_surface.feeding_foc_acf.png", width = 11, height = 3.6)

# trying with a bam
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
gam1 <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "fREML", family = binomial,
                            rho = -0.05, AR.start = ar.start, discrete = TRUE))
# HOW TO GET THE CORRECT RESIDUALS FROM A BAM? bam$rsd or residuals(bam)????
acf.pacf.plot(acf.structured(gam1, 10, "gam"), 10)


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year") # factor
var.num = c("day.of.year", "log.dist", "vess.speed.60",  "log.vess.speed.var.60", "arcsin.vess.di.60", "sqrt.encounter.minute") # numeric 
var.rand = c("folnum.unique") # random 

# formula
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "fREML", family = binomial,
                            rho = -0.05, AR.start = ar.start, discrete = TRUE))
summary(gam); AIC(gam) # dev = 54.9%, fREML = 880.1, AIC = 193.1

# final model (manual selection)
var.num = c("arcsin.vess.di.60", "sqrt.encounter.minute")
var.fac = c("group", "seastate", "year") 
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("bam", list(as.formula(f), data=as.name("tot"), method = "fREML", family = binomial,
                           rho = -0.1, AR.start = ar.start, discrete = TRUE))
summary(gam); AIC(gam) # dev = 54.4%, fREML = 881.0, AIC = 192.4

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam-v3/surface-feeding/foc/gam_surface-feeding_foc_final_v3.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam-v3/surface-feeding/foc/gam_surface-feeding_foc_final_v3.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam-v3/surface-feeding/foc/gam_surface-feeding_foc_results_v3.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
print(plot(b, allTerms = T), pages = 1)
var.num = c("arcsin.vess.di.60", "sqrt.encounter.minute")
var.fac = c("group", "seastate", "year") 
var.rand = c("folnum.unique")

# vessel DI
(sm <- plot(sm(b, 1)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "arcsin(vessel DI)", y = "s(arcsin(vessel DI))") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "surface.feeding", var = "arcsin.vess.di.60", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sin(arcsin.vess.di.60), y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Vessel DI", y = "Rate of surface feeding") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = asin.trans, breaks = c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/surface-feeding/foc/gam_surface-feeding_foc_vess.di.60_v3.png", height = 4.5, width = 13)


# encounter.minute
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(encounter minute)", y = "s(sqrt(encounter minute))") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "surface.feeding", var = "sqrt.encounter.minute", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = sqrt.encounter.minute^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Encounter minute", y = "Rate of surface feeding") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(0, 1,2,5,10,20,30)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/surface-feeding/foc/gam_surface-feeding_foc_enc.min_v3.png", height = 4.5, width = 13)


# group
v <- "group"
plotdata <- plot(b, allTerms = T, select = 4)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "surface.feeding", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Rate of surface feeding") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/surface-feeding/foc/gam_surface-feeding_foc_group_v3.png", height = 4.5, width = 13)


# group
v <- "seastate"
plotdata <- plot(b, allTerms = T, select = 5)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "surface.feeding", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Rate of surface feeding") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/surface-feeding/foc/gam_surface-feeding_foc_seastate_v3.png", height = 4.5, width = 13)


# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "surface.feeding", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Rate of surface feeding") + partial.theme + partial.grid.margin)
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam-v3/surface-feeding/foc/gam_surface-feeding_foc_year_v3.png", height = 4.5, width = 13)
