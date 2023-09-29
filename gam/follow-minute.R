###############################
### ANALYSIS: FOLLOW MINUTE ###
###############################
### 20/02/2022
### Tom Grove
### tomgrove20@yahoo.co.uk


### PACKAGES
packages <- c("tidyverse", "ppcor","RColorBrewer", "scales","survMisc", "Metrics", "corrplot", "lme4", "MASS", "mgcv", "tidymv", "mgcViz", "gridExtra", "gratia", "ggcorrplot", "cowplot", "ggpubr")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)


#---------------- FUNCTIONS + THEME --------------------

source("code/functions.R")
source("code/themes.R")


#---------------- DATA --------------------

tot <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv") %>%
  mutate(trip.id = paste0(as.Date(datetime),"_",vessel)) %>%
  group_by(trip.id, id.focal.catalogue) %>%
  # we only want the first follow for a trip for each animal
  filter(folnum.unique == first(folnum.unique)) %>%
  filter(!is.na(encounter.minute)) %>%
  group_by(folnum.unique) %>%
  filter(first(boat.ids.30.1500) == vessel | first(boat.ids.30.1500) == "") %>%
  mutate(sqrt.encounter.minute = sqrt(encounter.minute),
         folnum.unique = as.factor(folnum.unique)) %>%
  mutate(dummy = 1)

dive <- read.csv("intermediate-products/response-var-dfs/surface-interval_ais.csv") %>%
  rename(datetime = dt.end) %>%
  dplyr::select(datetime, ibi.mean, breath.num)

speed <- read.csv("intermediate-products/response-var-dfs/speed_ais.csv") %>%
  dplyr::select(datetime, speed, move.id.tot)

di <- read.csv("intermediate-products/response-var-dfs/di_ais.csv") %>%
  dplyr::select(datetime, DI)

tot <- tot %>% 
  left_join(dive) %>%
  left_join(speed) %>% 
  left_join(di)

ggplot(data = tot, aes(x = encounter.minute, y = dive.time)) + geom_point() + xlim(10,35) + geom_smooth(method = lm)

# speed + DI bootstrap
speed.boot <- read.csv("intermediate-products/bootstrap/bootstrap-speed.csv")
di.boot <- read.csv("intermediate-products/bootstrap/bootstrap-DI.csv")


#---------------- SURFACE ACTIVITY --------------------

df <- tot
summary(mod.sab <- gam(surface.active ~ s(sqrt.encounter.minute, bs = "cs", k = 10) +
                 s(folnum.unique, bs = "re", by = dummy), data = df, family = binomial))
plot(mod.sab, select = 1)

# plotting!
b <- getViz(mod.sab)
(sm <- plot( sm(b, 1) ) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(encounter minute)", y = "s(sqrt(encounter minute))") + partial.theme + partial.grid.margin)
(pr <- plot.folmin(mod.sab, df, resp = "surface.active", var = "sqrt.encounter.minute", var.rand = "folnum.unique") %>%
    ggplot(aes(x = sqrt.encounter.minute^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Encounter minute", y = "Rate of surface activity") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(0.1,1,2,5,10,20,30)) +
    coord_cartesian(ylim = c(0,0.15)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/follow-minute/gam_follow-minute_surface.active.png", height = 4.5, width = 13)


#---------------- SURFACE FEEDING --------------------

df <- tot
summary(mod.sfe <- gam(surface.feeding ~ s(sqrt.encounter.minute, bs = "cs", k = 10) +
                 s(folnum.unique, bs = "re", by = dummy), data = df, family = binomial))
plot(mod.sfe, select=1)


#---------------- DIVE TIME --------------------

# Need to filter for minute > 10
df <- tot %>% filter(encounter.minute>=10)
summary(mod.divetime <- gam(dive.time ~ s(sqrt.encounter.minute, bs = "cs", k = 10) +
                      s(folnum.unique, bs = "re", by = dummy), data = df, family = gaussian(link = "log")))

# plotting!
b <- getViz(mod.divetime)
(sm <- plot( sm(b, 1) ) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(encounter minute)", y = "s(sqrt(encounter minute))") + partial.theme + partial.grid.margin)
(pr <- plot.folmin(mod.divetime, df, resp = "dive.time", var = "sqrt.encounter.minute", var.rand = "folnum.unique") %>%
    ggplot(aes(x = sqrt.encounter.minute^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Encounter minute", y = "Dive time (sec)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(0.1,1,2,5,10,20,30)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/follow-minute/gam_follow-minute_dive.time.png", height = 4.5, width = 13)


#---------------- BREATH.NUM + MEAN IBI --------------------

# filtering for minute >1
df <- tot %>% filter(encounter.minute >= 1)

summary(mod.breath.num <- gam(breath.num ~ s(sqrt.encounter.minute, bs = "cs", k = 10) +
                        s(folnum.unique, bs = "re", by = dummy), data = df, family = poisson))
plot(mod.breath.num, select = 1)

summary(mod.ibi  <- gam(ibi.mean ~ s(sqrt.encounter.minute, bs = "cs", k = 10) +
                  s(folnum.unique, bs = "re", by = dummy), data = df, family = gaussian))
plot(mod.ibi, select = 1)


#---------------- SPEED --------------------

# GAM of speed against encounter minute
df <- tot
summary(mod.speed  <- gam(speed ~ s(sqrt.encounter.minute, bs = "cs", k = 10) +
                  s(folnum.unique, bs = "re", by = dummy), data = df, family = gaussian(link = "sqrt")))

# plotting!
b <- getViz(mod.speed)
(sm <- plot( sm(b, 1) ) + l_fitLine() + l_ciLine() + l_rug() + labs(x = "sqrt(encounter minute", y = "s(sqrt(encounter minute))") + partial.theme + partial.grid.margin)
(pr <- plot.folmin(mod.speed, df, resp = "speed", var = "sqrt.encounter.minute", var.rand = "folnum.unique") %>%
    ggplot(aes(x = sqrt.encounter.minute^2, y = predict)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
    labs(x = "Encounter minute", y = "Speed (km/h)") + partial.theme + partial.grid.margin +
    scale_x_continuous(trans = "sqrt", breaks = c(0.1,1,2,5,10,20,30)))
(comb <- as_ggplot(gridPrint(sm, pr, ncol = 2)) + 
    draw_plot_label(label = c("A", "B"), x = c(0, 0.5), y = c(1, 1), size = 20))
ggsave("intermediate-products/gam/follow-minute/gam_follow-minute_speed.png", height = 4.5, width = 13)


#---------------- SPEED BOOTSTRAP --------------------

n = 500 # number of iterations
boot.res <- data.frame()
for(i in 1:n) {
  
  # then replace speed with iteration
  sample <- tot %>% left_join(speed.boot %>% dplyr::select(datetime, move.id.tot, paste0("speed",i))) %>%
    dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  
  # now run the model
  summary <- summary(gam(speed ~ s(sqrt.encounter.minute, bs = "cs", k = 10) + s(folnum.unique, bs = "re"), 
                         data = sample, family = gaussian(link = "sqrt")))
  
  # now extracting important bits!
  smooth.names <- gsub("s\\(|\\)","",names(summary$chi.sq))
  smooth.edf <- t(as.data.frame(summary$edf)); rownames(smooth.edf) = NULL; colnames(smooth.edf) = paste0(smooth.names,"_estimate")
  smooth.p <- t(as.data.frame(summary$s.pv)); rownames(smooth.p) = NULL; colnames(smooth.p) = paste0(smooth.names,"_p")
  
  row <- cbind(smooth.edf, smooth.p)
  boot.res <- rbind(boot.res,row)
  
  if((i %% 50) == 0){print(paste0("iteration ",i))}
  
}

head(boot.res)


# plotting this out
# get effect size from original model
edf <- summary(mod.speed)$edf[1]

# what percentage of runs had p<0.05?
p <- (boot.res %>% summarise(p = paste0(sum(sqrt.encounter.minute_p<=0.05)/n()*100,"% p<0.05   ")))[1,1]

# then smooth vars
boot.res %>% ggplot() +
    geom_histogram(aes(x = sqrt.encounter.minute_estimate), alpha = 0.4, bins = 60, color = "black", fill = "lightgrey") +
    geom_text(aes(x = Inf, y = Inf, label = p),  
              vjust = 1.5, hjust = 1, size = 3.5) +
    geom_vline(aes(xintercept = edf), linetype = "dashed", alpha = 0.5, color = "red", size = 1) +
    labs(x = "EDF", y = "Count") +
    plot.theme + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
ggsave("intermediate-products/gam/follow-minute/gam_follow-minute_speed-boot_v2.png", width = 5, height = 3.5, dpi = 900)


#---------------- DI --------------------

df <- tot %>%
  mutate(arcsin.DI = asin(DI))

summary(mod.di <- gam(arcsin.DI ~ s(sqrt.encounter.minute, bs = "cs", k = 10) +
                    s(folnum.unique, bs = "re", by = dummy), data = df, family = gaussian))
plot(mod.di, select = 1) # not important



