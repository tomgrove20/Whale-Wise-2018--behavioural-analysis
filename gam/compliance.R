##############################################
### ANALYSIS: COMPLIANCE VS NON COMPLIANCE ###
##############################################
### 21/02/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

# how does whale behaviour change in compliance vs non-compliance? a series of ANOVAs
# within 300 m, 5-6 mph
# don't spend more than 20-30 minutes with a cetacean
# never approach directly in front
# when possible, stop the propellor
# using no boats within 10 min and 1500 m
 

### PACKAGES
packages <- c("tidyverse", "ppcor","RColorBrewer", "scales","survMisc", "Metrics", "corrplot", "lme4", "MASS", "mgcv", "tidymv", "mgcViz", "gridExtra", "gratia", "ggcorrplot", "cowplot", "ggpubr", "ggsignif")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)


#---------------- FUNCTIONS + THEME --------------------

source("code/functions.R")
source("code/themes.R")


#---------------- DATA --------------------

# total!
tot <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv") %>%  
  filter(!is.na(focal.vess.speed.10)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>%
  dplyr::select(datetime, dist,focal.vess.speed.10) %>%
  mutate(focal.vess.speed.10 = focal.vess.speed.10*2.23694) %>% # m/s to mph
  mutate(compliant.300 = ifelse(focal.vess.speed.10>6,"no", "yes"),
         compliant.50 = ifelse(focal.vess.speed.10>3,"no", "yes"))
tot.join <- tot %>% dplyr::select(-dist, -focal.vess.speed.10)

# surface feeding + activity
sab.sfe <- read.csv("intermediate-products/response-var-dfs/surface-active+feeding_ais.csv") %>%  
  filter(!is.na(focal.vess.speed.10)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>%  left_join(tot.join)

# dive stuff
dive <- read.csv("intermediate-products/response-var-dfs/dive-time_ais.csv") %>% 
  filter(!is.na(focal.vess.speed.10)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>%  left_join(tot.join)

# surfacing stuff
surfacing <- read.csv("intermediate-products/response-var-dfs/surface-interval_ais.csv") %>%
  mutate(datetime = dt.end) %>% 
  filter(!is.na(focal.vess.speed.10)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>%  left_join(tot.join)

# speed
speed <- read.csv("intermediate-products/response-var-dfs/speed_ais.csv") %>%
  filter(!is.na(focal.vess.speed.10)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>%  left_join(tot.join)

# di
di <- read.csv("intermediate-products/response-var-dfs/di_ais.csv") %>% 
  filter(!is.na(focal.vess.speed.10)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>%  left_join(tot.join)

# speed + DI bootstrap
speed.boot <- read.csv("intermediate-products/bootstrap/bootstrap-speed.csv")
di.boot <- read.csv("intermediate-products/bootstrap/bootstrap-DI.csv")


#---------------- PLOT SPEEDS --------------------

supp.labs <- c("Approaching zone: 50-300 m", "Encounter zone: <50 m")
names(supp.labs) <- c("within.300", "within.50")

tot %>% ungroup() %>%
  mutate(within.300 = ifelse(dist>300 | dist<50, NA, 1),
         within.50 = ifelse(dist>50, NA, 1)) %>%
  gather(key = "within", value = "val", within.300, within.50) %>%
  filter(!is.na(val)) %>%
  ggplot(aes(x = focal.vess.speed.10)) + 
  geom_histogram(color = "black", fill = "lightgrey", binwidth = 0.25, center = 0.125) +
  geom_vline(data = data.frame(within = c("within.300", "within.50"),
                               speed = c(6,3)), aes(xintercept = speed), 
             color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~within, ncol = 1, scales = "free_y",
               labeller = labeller(within = supp.labs)) + 
  plot.theme +
  labs(x = "Focal vessel speed (mph in previous 10 seconds)", y = "Number of surfacings")
ggsave("intermediate-products/compliance-effect/vess-speed_compliance_hist.png", dpi = 900)


#---------------- 5-6 mph within 50-300 m --------------------

# what proportion of surfacings were >6 mph?
df <- tot %>% filter(dist<=300 & dist>50)
nrow(df %>% filter(focal.vess.speed.10>6))/nrow(df) # 23.5% compliance

# number of follows?
f <- sab.sfe %>% filter(dist<=300 & dist>50 & !is.na(compliant.300)) %>% group_by(folnum.unique) %>% summarise(speed = max(focal.vess.speed.10))
nrow(f %>% filter(speed>(6/2.23)))/nrow(f)

# how about rates of surface feeding and activity
(df <- sab.sfe %>% filter(dist<=300 & dist>50 & !is.na(compliant.300)) %>%
  group_by(compliant.300, surface.active) %>%
  summarise(n = n()) %>% spread(key = "surface.active", value = "n") %>%
  column_to_rownames("compliant.300")) # no = 0.021, yes = 0.055
chisq.test(df) # X-squared = 4.626, df = 1, p-value = 0.03149

(df <- sab.sfe %>% filter(dist<=300 & dist>50 & !is.na(compliant.300)) %>%
    group_by(compliant.300, surface.feeding) %>%
    summarise(n = n()) %>% spread(key = "surface.feeding", value = "n") %>%
    column_to_rownames("compliant.300")) # no = 0.031, yes = 0.073
chisq.test(df) # X-squared = 6.8238, df = 1, p-value = 0.008995


# how about dive times?
df <- dive %>% filter(dist<= 300 & dist>50 & !is.na(compliant.300)) # 198 dives only, 38 = not compliant
summary(aov(log(dive.time)~compliant.300, data = df)) # F = 7.384, df = 1, p = 0.00695

(p.dive <- ggplot(data = df, aes(x = compliant.300, y = exp(log(dive.time)))) + geom_boxplot(width = 0.5)+
  scale_y_continuous(trans = "log", breaks = c(10,20,50,100,200,400,800)) +
  labs(x = "", y = "Dive time (sec)") + plot.theme)

# how about number of breaths?
df <- surfacing %>% filter(dist.end<= 300 & dist.end>50 & !is.na(compliant.300))# 401 surfacing intervals, 67 = not
kruskal.test(breath.num~compliant.300, data = df) # chi2 = 0.08007, df = 1, p = 0.78

(p.breath <- ggplot(data = df, aes(x = compliant.300, y = breath.num)) + geom_boxplot(width = 0.5)+
  labs(x = "", y = "Breaths per surfacing interval") + plot.theme)

# IBI
df <- df %>% filter(!is.na(ibi.mean))# 77 surfacing intervals, 60 = not
kruskal.test(ibi.mean~compliant.300, data = df)  # chi2 = 0.013, df = 1, p = 0.91

(p.ibi <- ggplot(data = df, aes(x = compliant.300, y = ibi.mean)) + geom_boxplot(width = 0.5)+
  labs(x = "", y = "Mean inter-breath interval (sec)") + plot.theme)

# speed
df <- speed %>% filter(dist<= 300 & dist>50 & !is.na(compliant.300)) # 931 rows, 293 not compliant
summary(aov(sqrt(speed)~compliant.300, data = df)) # F = 16.12, df = 1, p = 0.000064
df %>% group_by(compliant.300) %>% summarise(mean = mean(speed))

(p.speed <- ggplot(data = df, aes(x = compliant.300, y = speed)) + geom_boxplot(width = 0.5)+
  labs(x = "Speed compliance in approaching zone", y = "Speed (km/h)") + plot.theme)

# bootstrap speed
df <- df %>% left_join(speed.boot)
aov.boot <- data.frame()
n = 500
for (i in 1:n) {
  dat <- df %>% dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  test <- summary(aov(sqrt(speed)~compliant.300, data = dat))[[1]]
  row <- data.frame(test = "aov", var = "compliant.300", iteration = i, f = test[1,4], p = test[1,5])
  aov.boot <- bind_rows(aov.boot, row)
}
write.csv(aov.boot, "intermediate-products/compliance-effect/comp-300_speed_aov_boot.csv", row.names = FALSE)
nrow(aov.boot %>% filter(p<0.05))/nrow(aov.boot) # all iterations significant

# DI
df <- di %>% filter(dist<= 300 & dist>50 & !is.na(compliant.300)) # 768 rows, 184 not compliant
kruskal.test(asin(DI)~compliant.300, data = df) # chi2 = 0.26, df = 1, p = 0.61
df %>% group_by(compliant.300) %>% summarise(mean = mean(DI))

(p.di <- ggplot(data = df, aes(x = compliant.300, y = DI)) + geom_boxplot(width = 0.5)+
  scale_y_continuous(trans = asin.trans, breaks = c(0,0.2, 0.4, 0.6, 0.8, 0.9, 0.99)) +
  labs(x = "", y = "Directness index") + plot.theme)

# bootstrap DI
df <- df %>% left_join(di.boot)
aov.boot <- data.frame()
n = 500
for (i in 1:n) {
  dat <- df %>% dplyr::select(-DI) %>% rename(DI = paste0("DI",i))
  test <- kruskal.test(asin(DI)~compliant.300, data = dat) 
  row <- data.frame(test = "kruskal", var = "compliant.300", iteration = i, 
                    chi.sq = as.numeric(test$statistic), p = test$p.value)
  aov.boot <- bind_rows(aov.boot, row)
}
write.csv(aov.boot, "intermediate-products/compliance-effect/comp-300_DI_aov_boot.csv", row.names = FALSE)
nrow(aov.boot %>% filter(p<0.05))/nrow(aov.boot) # no iterations significant

# Bring it all together!
plot_grid(p.dive, p.breath, p.ibi, p.speed, p.di, align = "hv")
ggsave("intermediate-products/compliance-effect/approaching-zone-vars_no-yes.png", dpi = 600, height = 8, width = 14, bg = "white")

# only significant vars
bind_rows(speed %>% dplyr::select(datetime, speed, compliant.300), dive %>% dplyr::select(datetime, dive.time, compliant.300)) %>%
  gather(key = "var", value = "val", speed, dive.time) %>%
  mutate(var = factor(var, levels = c("dive.time", "speed"),
                      labels = c("Dive time (sec)", "Swim speed (km/h)"))) %>% 
  ggplot(aes(x = compliant.300, y = val)) +
  geom_boxplot(width = 0.5) +
  facet_wrap(~var, scales = "free_y", strip.position = "left") +
  geom_signif(data = data.frame(var = "Swim speed (km/h)"),
              aes(y_position=20, xmin="no", xmax="yes",
                  annotations="*"), textsize = 6, tip_length=0, manual = T) +
  labs(x = "Speed compliance in approaching zone", y = "") + 
  scale_y_continuous(expand = expansion(mult = c(0, .07)))+
  plot.theme + theme(strip.background = element_blank(), strip.placement = "outside",
                     panel.spacing = unit(1.5,"lines"))
ggsave("intermediate-products/compliance-effect/approaching-zone-sig-vars_v2.png", dpi = 600, height = 4, width = 9)

#---------------- 20-30 minutes with a cetacean --------------------

enc <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv") %>%
  group_by(folnum.unique) %>% summarise(follow.length = max(encounter.minute, na.rm = TRUE)) %>%
  mutate(compliant = ifelse(follow.length>30, "no", "yes")) %>%
  filter(follow.length > 0)
nrow(filter(enc, compliant == "no"))/nrow(enc) # 5% of follows were >30 minutes. At the moment, we don't have info on encounter duration


#---------------- Encounter zone, stop the propellor --------------------

# what speed = stop prop?
prop.stop = 3 # mph

df <- tot %>% filter(dist<=50) # 190 surfacings

# what proportion of surfacings were >3 mph?
nrow(df %>% filter(compliant.50 == "no"))/nrow(df) # 32.9% no compliance with other boats in vicinity

# how about rates of surface feeding and activity
(df <- sab.sfe %>% filter(dist<=50 & !is.na(compliant.50)) %>%
    group_by(compliant.50, surface.active) %>%
    summarise(n = n()) %>% spread(key = "surface.active", value = "n") %>%
    column_to_rownames("compliant.50")) # only 4 surface-active total anyway!
chisq.test(df) # X-squared = 2.3918e-30, df = 1, p-value = 1

(df <- sab.sfe %>% filter(dist<=50 & !is.na(compliant.50)) %>%
    group_by(compliant.50, surface.feeding) %>%
    summarise(n = n()) %>% spread(key = "surface.feeding", value = "n") %>%
    column_to_rownames("compliant.50")) # no = 0.031, yes = 0.073
chisq.test(df) # X-squared = 3.6316, df = 1, p-value = 0.05669


# how about dive times?
df <- dive %>% filter(dist<= 50 & !is.na(compliant.50)) # 40 dives only, 18 = not compliant
summary(aov(log(dive.time)~compliant.50, data = df)) # compliant.50  1  0.833  0.8331   1.253  0.269
df %>% group_by(compliant.50) %>% summarise(dive.time = mean(dive.time))

(p.dive <- ggplot(data = df, aes(x = compliant.50, y = exp(log(dive.time)))) + geom_boxplot(width = 0.5)+
  scale_y_continuous(trans = "log", breaks = c(10,20,50,100,200,400,800)) +
  labs(x = "", y = "Dive time (sec)") + plot.theme)


# how about number of breaths?
df <- surfacing %>% filter(dist.end<= 50 & !is.na(compliant.50))# 24 surfacing intervals
kruskal.test(breath.num~compliant.50, data = df) # chi-squared = 11.254, df = 1, p-value = 0.0007947
df %>% group_by(compliant.50) %>% summarise(mean = mean(breath.num)) # 1.82 = no, 5.92 = yes
(p.breath <- ggplot(data = df, aes(x = compliant.50, y = breath.num)) + geom_boxplot(width = 0.5)+
  labs(x = "", y = "Breaths per surfacing interval") + plot.theme)

# IBI
df <- df %>% filter(!is.na(ibi.mean))# 7 surfacing intervals, not worth it


# speed
df <- speed %>% filter(dist<= 50 & !is.na(compliant.50)) # 115 rows
summary(aov(sqrt(speed)~compliant.50, data = df)) # compliant.50   1  0.501  0.5010   2.619  0.108
df %>% group_by(compliant.50) %>% summarise(speed = mean(speed))
(p.speed <- ggplot(data = df, aes(x = compliant.50, y = speed)) + geom_boxplot(width = 0.5)+
  scale_y_continuous(trans = "sqrt", breaks = c(1,2,5,10,15,20)) +
  labs(x = "Speed <3 mph in encounter zone", y = "Speed (km/h)") + plot.theme)

# bootstrap speed
df <- df %>% left_join(speed.boot)
aov.boot <- data.frame()
n = 500
for (i in 1:n) {
  dat <- df %>% dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  test <- summary(aov(sqrt(speed)~compliant.50, data = dat))[[1]]
  row <- data.frame(test = "aov", var = "compliant.50", iteration = i, f = test[1,4], p = test[1,5])
  aov.boot <- bind_rows(aov.boot, row)
}
write.csv(aov.boot, "intermediate-products/compliance-effect/comp-50_speed_aov_boot.csv", row.names = FALSE)
nrow(aov.boot %>% filter(p<0.05))/nrow(aov.boot) # all iterations significant

# DI
df <- di %>% filter(dist<= 50 & !is.na(compliant.50)) # 100 rows
kruskal.test(asin(DI)~compliant.50, data = df) # chi-squared = 0.084706, df = 1, p-value = 0.771
df %>% group_by(compliant.50) %>% summarise(DI = mean(DI))
(p.di <- ggplot(data = df, aes(x = compliant.50, y = DI)) + geom_boxplot(width = 0.5)+
  scale_y_continuous(trans = asin.trans, breaks = c(0,0.2, 0.4, 0.6, 0.8, 0.9, 0.99)) +
  labs(x = "", y = "Directness index") + plot.theme)

# bootstrap DI
df <- df %>% left_join(di.boot)
aov.boot <- data.frame()
n = 500
for (i in 1:n) {
  dat <- df %>% dplyr::select(-DI) %>% rename(DI = paste0("DI",i))
  test <- kruskal.test(asin(DI)~compliant.50, data = dat) 
  row <- data.frame(test = "kruskal", var = "compliant.50", iteration = i, 
                    chi.sq = as.numeric(test$statistic), p = test$p.value)
  aov.boot <- bind_rows(aov.boot, row)
}
write.csv(aov.boot, "intermediate-products/compliance-effect/comp-50_DI_aov_boot.csv", row.names = FALSE)
nrow(aov.boot %>% filter(p<0.05))/nrow(aov.boot) # no iterations significant

# Bring it all together!
plot_grid(p.dive, p.breath, p.speed, p.di, align = "hv")
ggsave("intermediate-products/compliance-effect/encounter-zone-vars_no-yes.png", dpi = 600, height = 8, width = 9.2, bg = "white")
