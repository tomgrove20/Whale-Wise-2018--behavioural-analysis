##################################
### GAM: SURFACE FEEDING, AIS ###
##################################
### 14/02/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

# 1. transformations
# 2. Collinearity
# 3. time frame selection
# 4. Final GAM (GCV + backward selection)
# 5. Results + plotting

# note: PQL is used in gamm or gamm4, which performs poorly with binary data. So we're sticking with gam for surface activity and surface feeding



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
tot <- read.csv("intermediate-products/response-var-dfs/surface-active+feeding_ais.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique"), as.factor) %>%
  mutate(dummy = 1) # important dummy variable to factor out random effect when predicting response

# quick look at the data set
colSums(is.na(tot)) 


#---------------- TRANSFORMATIONS --------------------

# define explanatory variables, then view distribution
vars <- c("seastate", "group", "day.of.year", "oak.meandist.30.1500", "rib.meandist.30.1500", "oak.dist.max.30.1500", "rib.dist.max.30.1500", "oak.num.30.1500", "rib.num.30.1500")
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
ggplot(data = tot, aes(x = rib.meandist.30.1500)) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(oak.meandist.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(oak.dist.max.30.1500))) + geom_histogram()
ggplot(data = tot, aes(x = sqrt(rib.dist.max.30.1500))) + geom_histogram()

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains("meandist"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt mean dist
         across(contains("dist.max"), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}")) # sqrt dist max


#---------------- COLLINEARITY --------------------

# define explanatory variables
vars <- c("seastate", "group", "day.of.year", "oak.num.30.1500", "rib.num.30.1500","sqrt.oak.meandist.30.1500", "sqrt.rib.meandist.30.1500", "sqrt.oak.dist.max.30.1500", "sqrt.rib.dist.max.30.1500",  "year")

# correlation matrix
x <- tot[,gsub("`","",vars)] %>% drop_na()
model.matrix(~0+., data=x) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size = 3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1))
ggsave("intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_collinear.png")

# as expected, several correlations.
# keeping meandist instead of dist.max for now. rib num has lower correlation with mean than max
# there is a worrying correlation for rib.num vs rib.meandist, but I'm keeping them both for now (marginal, 0.73)
# note: that worrying correlation is 0.61 if not transformed to sqrt, but has to happen unfortunately


#---------------- TIME/DISTANCE FRAME --------------------

# first defining explanatory variables (not running folnum.unique as random in this)
var.fac = c("group", "seastate", "year") # factor
var.rand = c("folnum.unique") # random 
var.num = c("day.of.year") # numeric 
# specifying time/distance frame
ts = c("30.1500", "10.500", "10.1000", "10.1500") # time/distance frame
var.t = c("oak.num", "rib.num", "sqrt.oak.meandist", "sqrt.rib.meandist") # varying vars!

# Combining together into a lagged data frame
t.df <- data.frame(var = rep(var.t, each = length(ts)), t = rep(ts,length(var.t)), t.orig = "10.1000") %>%
  mutate(var.orig = paste0(var,".",t.orig), var.target = paste0(var,".",t))

(var.dy <- c(var.num, unique(t.df$var.orig)))

# now running the loop!
for (i in 1:nrow(t.df)) {
  
  # variable list, different time/distance frame each time
  v <- var.dy %>% recode(!!t.df[i,"var.orig"] := t.df[i,"var.target"])
  
  f1 <- as.formula(paste("surface.feeding ~", 
                         paste("s(",v,", bs = 'cs', k=5)", collapse= "+"),"+", # smooths
                         paste(var.fac, collapse = "+")))
  
  # run binomial GAM
  gam <- do.call("gam", list(as.formula(f1), data=as.name("tot"), method = "GCV.Cp", family = binomial))
  
  # and add AIC to df!
  t.df[i,"AIC"] = AIC(gam)
  
  print (t.df[i,"var.target"]) # progress checker
}

t.df # oak.num: 30.1500. oak.dist: 10.1000. rib.num: 10.1500. rib.dist: 10.1000
ggplot(data = t.df, aes(x = t, y = AIC, group = var, color = var)) +
  geom_point(size = 2, alpha = 0.5) + geom_line() + 
  scale_color_brewer(palette = "Dark2") +
  scale_x_discrete(labels=c("10 min, 500 m", "10 min, 1000 m", "10 min, 1500 m", "30 min, 1500 m")) +
  labs(x = "", color = "Variable") + plot.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_time-dist_aic.png", dpi = 600, height = 5, width = 8)

# going for 10.1000 for oak and 10.1500 for rib


#---------------- BACKWARD SELECTION --------------------

# full model
# defining variable combinations
var.fac = c("group", "seastate", "year")
var.int = c("oak.num.10.1000", "rib.num.10.1500")
var.num = c("day.of.year", "sqrt.oak.meandist.10.1000", "sqrt.rib.meandist.10.1500")
var.rand = c("folnum.unique")

# formula
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.int,", bs = 'cs', k = 5)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# GAM!
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 57.4%, ubre = -0.737. Remove oak.num.10.1000

# Remove oak.num.10.1500 and rib.num.10.1500
# simply removing the int term!!
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 57.4%, ubre = -0.737. no change, so remove

# Remove sqrt.rib.meandist.10.1500
var.num = c("day.of.year", "sqrt.oak.meandist.10.1000")
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 57.2%, ubre = -0.737. bad and neutral, don't remove

# Remove sqrt.oak.meandist.10.1000
var.num = c("day.of.year", "sqrt.rib.meandist.10.1500")
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 57.2%, ubre = -0.737. bad and neutral, don't remove

# Remove day.of.year
var.num = c("sqrt.rib.meandist.10.1500", "sqrt.oak.meandist.10.1000")
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 57.1%, ubre = -0.731. bad changes, don't remove

# Remove day.of.year
var.num = c("day.of.year")
f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 57%, ubre = -0.737. bad and neutral, don't remove


# final model
var.fac = c("group", "seastate", "year")
var.num = c("day.of.year", "sqrt.oak.meandist.10.1000", "sqrt.rib.meandist.10.1500")
var.rand = c("folnum.unique")

f <- as.formula(paste("surface.feeding ~", 
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))
gam <- do.call("gam", list(as.formula(f), data=as.name("tot"), method = "GCV.Cp", family = binomial))
summary(gam) # dev = 57.4%, ubre = -0.738

# that is the final model! Let's check for concurvity
concurvity(gam); concurvity(gam, full = FALSE) # only high concurvity between folnum.unique and other vars. Not uncommon with random variables, ignoring for now

# saving the final model!
gam.final = gam; saveRDS(gam.final, "intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_final.rds")


#---------------- MODEL RESULTS --------------------
# gam <- readRDS("intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_final.rds")

# first saving estimates
res <- bind_rows(as.data.frame(summary(gam)$p.table), as.data.frame(summary(gam)$s.table)) %>%
  rownames_to_column("var")
write.csv(res, "intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_results.csv", row.names = FALSE)

# get ready for plotting!
b <- getViz(gam)

# checking normality assumptions
check(b, a.qq = list(method = "simul2", a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

## PARTIALS + PREDICTED
print(plot(b, allTerms = T), pages = 1)
var.fac = c("group", "seastate", "year")
var.num = c("day.of.year", "sqrt.oak.meandist.10.1000", "sqrt.rib.meandist.10.1500")
var.rand = c("folnum.unique")

# day.of.year
(sm <- plot( sm(b, 1) ) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "Julian day", y = "s(Julian day)") + partial.theme + partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "surface.feeding", var = "day.of.year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = day.of.year, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Julian day", y = "Rate of surface feeding") + partial.theme + partial.grid.margin +
    coord_cartesian(ylim = c(0,0.05)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_julian.png", height = 4.5, width = 13)

# sqrt.oak.meandist.10.1000
(sm <- plot(sm(b, 2)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(oak boat mean distance (m in 10 min, 1000 m))", y = "s(sqrt(oak boat mean distance))") + partial.theme +
    partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "surface.feeding", var = "sqrt.oak.meandist.10.1000", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = (sqrt.oak.meandist.10.1000)^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    scale_x_continuous(trans = "sqrt", breaks = c(50,200,500,1000,2000,3000)) +
    labs(x = "Oak boat mean distance (m in 10 min, 1000 m)", y = "Rate of surface feeding") + 
    partial.theme + partial.grid.margin + coord_cartesian(ylim = c(0,0.04)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_oak.meandist.10.1000.png", height=4.5, width=13)

# sqrt.rib.meandist.10.1500
(sm <- plot(sm(b, 3)) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(RIB mean distance (m in 10 min, 1500 m))", y = "s(sqrt(RIB mean distance))") + partial.theme +
    partial.grid.margin)
(pr <- pred.var(gam, tot, resp = "surface.active", var = "sqrt.rib.meandist.10.1500", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = (sqrt.rib.meandist.10.1500)^2, y = predict)) + geom_line() +
    scale_x_continuous(trans = "sqrt", breaks = c(10, 100, 500, 1000,2000, 3000)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "RIB mean distance (m in 10 min, 1500 m)", y = "Rate of surface feeding") + partial.theme + partial.grid.margin +
    coord_cartesian(ylim = c(0,0.04)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_rib.meandist.10.1500.png", height=4.5, width=13)

# group
v <- "group"
plotdata <- plot(b, allTerms = T, select = 5)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Group type", y = "Partial effect of group type"))
(pr <- pred.var.fac(gam, tot, resp = "surface.feeding", var = "group", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = group, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Group type", y = "Rate of surface feeding") + partial.theme + partial.grid.margin + 
    coord_cartesian(ylim = c(0,0.1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_group.png", height = 4.5, width = 13)

# sea state
v <- "seastate"
plotdata <- plot(b, allTerms = T, select = 6)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Sea state", y = "Partial effect of sea state"))
(pr <- pred.var.fac(gam, tot, resp = "surface.feeding", var = "seastate", var.fac, var.num, var.rand) %>% drop_na() %>%
    ggplot(aes(x = seastate, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Sea state", y = "Rate of surface feeding") + partial.theme + partial.grid.margin + 
    coord_cartesian(ylim = c(0,0.1)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_seastate.png", height = 4.5, width = 13)

# year
v <- "year"
plotdata <- plot(b, allTerms = T, select = 7)$plots[[1]]$data$fit %>% mutate(ci = 2*se)
jitter.pos <- jit.pos(plotdata,y,ci); jitter.height = jit.height(plotdata,y,ci)
(sm <- eval(var.fac.plot) +
    labs(x = "Year", y = "Partial effect of year"))
(pr <- pred.var.fac(gam, tot, resp = "surface.feeding", var = "year", var.fac, var.num, var.rand) %>% 
    ggplot(aes(x = year, y = predict)) + geom_boxplot(width = 0.4) +
    labs(x = "Year", y = "Rate of surface feeding") + partial.theme + partial.grid.margin + 
    coord_cartesian(ylim = c(0,0.15)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/surface-feeding/ais/gam_surface-feeding_ais_year.png", height = 4.5, width = 13)

