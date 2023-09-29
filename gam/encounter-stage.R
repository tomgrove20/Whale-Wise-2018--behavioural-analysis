#################################
### ANALYSIS: ENCOUNTER STAGE ###
#################################
### 20/02/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

# a mixture of chi-squared tests and ANOVAs

### PACKAGES
packages <- c("tidyverse", "ppcor","RColorBrewer", "scales","survMisc", "Metrics", "corrplot", "lme4", "MASS", "mgcv", "tidymv", "mgcViz", "gridExtra", "gratia", "ggcorrplot", "cowplot", "ggpubr", "ggeffects")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)


#---------------- FUNCTIONS + THEME --------------------

source("code/functions.R")
source("code/themes.R")


#---------------- DATA --------------------

# total
tot <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv") %>%
  mutate_at(.vars = c("surface.active", "surface.feeding"), funs(recode(., `1` = "yes", `0` = "no"))) %>%
  mutate_at(.vars = c("surface.active", "surface.feeding"), funs(as.factor)) %>%
  group_by(folnum.unique) %>% mutate(dive.time = lead(dive.time)) %>% ungroup()
# decided to keep dive time attached to second surfacing (allows us to look at after!)

surfacing <- read.csv("intermediate-products/response-var-dfs/surface-interval_ais.csv") %>%
  rename(datetime = dt.end) %>% dplyr::select(datetime, ibi.mean, breath.num)

speed <- read.csv("intermediate-products/response-var-dfs/speed_ais.csv") %>%
  dplyr::select(datetime, speed)

di <- read.csv("intermediate-products/response-var-dfs/di_ais.csv") %>%
  dplyr::select(datetime, DI)

tot <- tot %>%
  left_join(surfacing, by = "datetime") %>% left_join(speed, by = "datetime") %>% left_join(di, by = "datetime")

# before/during/after
tot.stage <- tot %>% ungroup() %>%
  filter(!is.na(encounter.stage)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>% # going for 10 minutes, 1500 m
  mutate(encounter.stage = fct_relevel(encounter.stage, "before", "during", "after"))

# approach/encounter/depart
tot.enc <- tot %>% ungroup() %>%
  filter(!is.na(vess.encounter.stage)) %>%
  filter(boat.ids.10.1500 == vessel | boat.ids.10.1500 == "") %>% # going for 10 minutes, 1500 m
  mutate(vess.encounter.stage = recode(vess.encounter.stage, `-1` = "approach", `0` = "encounter", `1` = "departure"),
         vess.encounter.stage = fct_relevel(vess.encounter.stage, "approach", "encounter", "departure"))

# quick look at the data set
colSums(is.na(tot))

# speed + DI bootstrap
speed.boot <- read.csv("intermediate-products/bootstrap/bootstrap-speed.csv")
di.boot <- read.csv("intermediate-products/bootstrap/bootstrap-DI.csv")


#######################
### ENCOUNTER STAGE ###
#######################
# before 1st <300 m, during, after last 300 m

#---------------- Surface activity + feeding --------------------

# surface feeding = chi-squared test
(df = tot.stage %>% filter(dist<1000) %>%
   group_by(encounter.stage, surface.feeding) %>% summarise(n = n()) %>% ungroup() %>%
   spread(key = "surface.feeding", value = "n") %>% as.data.frame() %>% mutate(yes = replace_na(yes,0)))
(test <-  chisq.test(df %>% column_to_rownames("encounter.stage"))) # X-squared = 2.2853, df = 2, p-value = 0.319

# surface activity = chi-squared test
(df = tot.stage %>% filter(dist<1000) %>%
    group_by(encounter.stage, surface.active) %>% summarise(n = n()) %>% ungroup() %>%
    spread(key = "surface.active", value = "n") %>% as.data.frame() %>% mutate(yes = replace_na(yes,0)))
 # before  = 0.03, during = 0.06, after = 0.12
(test <-  chisq.test(df %>% column_to_rownames("encounter.stage"))) # X-squared = 3.8974, df = 2, p-value = 0.1425

# plotting!
supp.labs <- c("Surface active", "Surface.feeding"); names(supp.labs) <- c("surface.active", "surface.feeding")
tot.stage %>% filter(dist<1000) %>% 
  gather(key = "type", value = "value", surface.active, surface.feeding) %>%
  group_by(encounter.stage, type, value) %>% summarise(n = n()) %>% 
  group_by(encounter.stage, type) %>% mutate(prop = n/sum(n), n = sum(n)) %>% ungroup %>%
  ggplot(aes(x = encounter.stage, y = prop)) +
    geom_bar(aes(fill = value), stat = "identity") +
  geom_text(aes(label= n, y = 1.07)) + 
  facet_wrap(~type, ncol = 1, labeller = labeller(type = supp.labs)) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Encounter stage", y = "Proportion of surfacings") + plot.theme + theme(legend.position = "none")
ggsave("intermediate-products/encounter-stage/enc-stage_surface-active+feeding.png", width = 5, height = 7)


#---------------- Dive time --------------------

# log.dive.time = kruskal.wallis
df = tot.stage %>% filter(dive.time>30) # 437 dives
df %>% group_by(encounter.stage) %>% summarise(n = n(), mean = mean(dive.time)) # 41 before, 396 during
summary(test <- aov(log(dive.time) ~ encounter.stage, data = df)) # 1   1.18  1.1830   1.861  0.173

(p.dive <- df %>% bind_rows(data.frame(encounter.stage = "after")) %>%
    mutate(encounter.stage = fct_relevel(encounter.stage, "before", "during", "after")) %>%
    ggplot(aes(x = encounter.stage, y = dive.time)) +
  geom_point(position = position_jitter(width=0.23), size = 2, color = "grey", alpha = 0.3)+
  geom_boxplot(size = 0.8, width = 0.5, fill = "transparent") + labs(x = "", y = "Dive time (sec)") +
  scale_y_continuous(trans = "log", breaks = c(10,20,50,100,200,400,800)) + plot.theme)


#---------------- Surfacing interval --------------------

# number of breaths
df = tot.stage %>% filter(!is.na(breath.num) & !is.na(encounter.stage)) # only 528!
df %>% group_by(encounter.stage) %>% summarise(n = n(), mean = mean(breath.num)) # only 9 before!

# IBI 
df = tot.stage %>% filter(!is.na(ibi.mean) & !is.na(encounter.stage)) # only 100!
df %>% group_by(encounter.stage) %>% summarise(n = n(), mean = mean(ibi.mean)) # only 1 before!


#---------------- Speed and DI --------------------

# speed
df = tot.stage %>% filter(!is.na(encounter.stage) & !is.na(speed))# 1266
df %>% group_by(encounter.stage) %>% summarise(n = n(), mean = mean(speed)) # some clear differences ...
summary(test <- aov(sqrt(speed) ~ encounter.stage, data = df)) # 2    2.6  1.3049   4.503 0.0113 *

(p.speed <- ggplot(data = df, aes(x = encounter.stage, y = speed)) +
  geom_point(position = position_jitter(width=0.23), size = 2, color = "grey", alpha = 0.3)+
  geom_boxplot(size = 0.8, width = 0.5, fill = "transparent") + labs(x = "", y = "Speed (km/h)") +
  scale_y_continuous(trans = "sqrt", breaks = c(0.5, 2, 5, 10, 15, 20)) + plot.theme)

# bootstrap speed
df <- df %>% left_join(speed.boot)
aov.boot <- data.frame();n = 500
for (i in 1:n) {
  dat <- df %>% dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  test <- summary(aov(sqrt(speed)~encounter.stage, data = dat))[[1]]
  row <- data.frame(test = "aov", var = "encounter.stage", iteration = i, f = test[1,4], p = test[1,5])
  aov.boot <- bind_rows(aov.boot, row)
}
write.csv(aov.boot, "intermediate-products/encounter-stage/enc-stage_speed_aov_boot.csv", row.names = FALSE)
nrow(aov.boot %>% filter(p<0.05))/nrow(aov.boot) # all iterations significant

# DI
df = tot.stage %>% filter(!is.na(DI) & !is.na(encounter.stage)) # 1058
df %>% group_by(encounter.stage) %>% summarise(n = n(), mean = mean(DI)) # some clear differences ...
kruskal.test(asin(DI) ~ encounter.stage, data = df)

(p.di <- ggplot(data = df, aes(x = encounter.stage, y = DI)) +
  geom_point(position = position_jitter(width=0.23), size = 2, color = "grey", alpha = 0.3)+
  geom_boxplot(size = 0.8, width = 0.5, fill = "transparent") + labs(x = "Encounter stage", y = "Directness index") +
  scale_y_continuous(trans = "asn", breaks = c(0.05, 0.2, 0.5, 0.7, 0.9, 0.99)) + plot.theme +
    coord_cartesian(ylim = c(0.2,1)))

# bootstrap DI
df <- df %>% left_join(di.boot)
aov.boot <- data.frame(); n = 500
for (i in 1:n) {
  dat <- df %>% dplyr::select(-DI) %>% rename(DI = paste0("DI",i))
  test <- kruskal.test(asin(DI)~encounter.stage, data = dat) 
  row <- data.frame(test = "kruskal", var = "encounter.stage", iteration = i, 
                    chi.sq = as.numeric(test$statistic), p = test$p.value)
  aov.boot <- bind_rows(aov.boot, row)
}
write.csv(aov.boot, "intermediate-products/encounter-stage/enc-stage_DI_aov_boot.csv", row.names = FALSE)
nrow(aov.boot %>% filter(p<0.05))/nrow(aov.boot) # no iterations significant

# Bring it all together!
plot_grid(p.dive, p.speed, p.di, ncol = 1, nrow = 3, align = "hv")
ggsave("intermediate-products/encounter-stage/encounter-stage-vars.png", dpi = 600, height = 10, width = 7, bg = "white")


##############################
### VESSEL ENCOUNTER STAGE ###
#############################
# approach (to 200 m, >2 m/s), encounter (within 200 m), depart (>200 m, >2 m/s)

# we need to plot dist vs speed to justify having 
ggplot(data = tot, aes(x = dist, y = focal.vess.speed.10*2.23694)) +
  xlim(0,800) +
  geom_point(alpha = 0.05) + geom_smooth(color = "black") + 
  labs(x = "Distance (m)", y = "Vessel speed (mph, 10 sec)") + plot.theme
ggsave("intermediate-products/encounter-stage/dist-vs-speed.10.png", dpi = 900)

#---------------- Surface activity + feeding --------------------

# surface activity = chi-squared test
df = tot.enc %>% filter(dist<1000 & vess.encounter.stage != "departure") %>%
   group_by(vess.encounter.stage, surface.active) %>% summarise(n = n()) %>% ungroup() %>%
   spread(key = "surface.active", value = "n") %>% as.data.frame() %>% mutate(yes = replace_na(yes,0))
df %>% mutate(perc = yes/(yes+no))
(test <-  chisq.test(df %>% column_to_rownames("vess.encounter.stage")))

# surface feeding = chi-squared test
df = tot.enc %>% filter(dist<1000 & vess.encounter.stage != "departure") %>%
    group_by(vess.encounter.stage, surface.feeding) %>% summarise(n = n()) %>% ungroup() %>%
    spread(key = "surface.feeding", value = "n") %>% as.data.frame() %>% mutate(yes = replace_na(yes,0))
df %>% mutate(perc = yes/(yes+no))
(test <-  chisq.test(df %>% column_to_rownames("vess.encounter.stage")))

# plotting!
supp.labs <- c("Surface active", "Surface.feeding"); names(supp.labs) <- c("surface.active", "surface.feeding")
tot.enc %>% filter(dist<1000) %>% 
  gather(key = "type", value = "value", surface.active, surface.feeding) %>%
  group_by(vess.encounter.stage, type, value) %>% summarise(n = n()) %>% 
  group_by(vess.encounter.stage, type) %>% mutate(prop = n/sum(n), n = sum(n)) %>% ungroup %>%
  ggplot(aes(x = vess.encounter.stage, y = prop)) +
  geom_bar(aes(fill = value), stat = "identity") +
  geom_text(aes(label= n, y = 1.07)) + 
  facet_wrap(~type, ncol = 1, labeller = labeller(type = supp.labs)) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Vessel behaviour", y = "Proportion of surfacings") + plot.theme + theme(legend.position = "none")

ggsave("intermediate-products/encounter-stage/vess-enc-stage_surface-active+feeding.png", width = 5, height = 7)


#---------------- Dive time --------------------

# log.dive.time = kruskal.wallis
df = tot.enc %>% filter(dive.time>30 & !is.na(vess.encounter.stage)) # 420 dives
df %>% group_by(vess.encounter.stage) %>% summarise(n = n(), mean = mean(dive.time)) # 34 before, 159 during
df <- df %>% filter(vess.encounter.stage != "departure")
summary(test <- aov(log(dive.time) ~ vess.encounter.stage, data = df))

(p.dive <- ggplot(data = df, aes(x = vess.encounter.stage, y = dive.time)) +
  geom_point(position = position_jitter(width=0.23), size = 2, color = "grey", alpha = 0.3)+
  geom_boxplot(size = 0.8, width = 0.5, fill = "transparent") + labs(x = "", y = "Dive time (sec)") +
  scale_y_continuous(trans = "log", breaks = c(10,20,50,100,200,400,800)) + plot.theme)


#---------------- Surfacing interval --------------------

# number of breaths
df = tot.enc %>% filter(!is.na(breath.num) & !is.na(vess.encounter.stage)) # only 106
df %>% group_by(vess.encounter.stage) %>% summarise(n = n(), mean = mean(breath.num)) # 17 before and 102 after
(test <- kruskal.test(breath.num ~ vess.encounter.stage, data = df))

(p.breath <- ggplot(data = df, aes(x = vess.encounter.stage, y = breath.num)) +
    geom_point(position = position_jitter(width=0.23), size = 2, color = "grey", alpha = 0.3)+
    geom_boxplot(size = 0.8, width = 0.5, fill = "transparent") + labs(x = "", y = "Breaths per surfacing interval") +
    plot.theme)

# not enough for IBI


#---------------- Speed and DI --------------------

# speed
df = tot.enc %>% filter(!is.na(speed) & !is.na(vess.encounter.stage))# 488
df %>% group_by(vess.encounter.stage) %>% summarise(n = n(), mean = mean(speed)) # some clear differences ...
df <- df %>% filter(vess.encounter.stage != "departure")
summary(test <- aov(sqrt(speed) ~ vess.encounter.stage, data = df))

(p.speed <- ggplot(data = df, aes(x = vess.encounter.stage, y = speed)) +
  geom_point(position = position_jitter(width=0.23), size = 2, color = "grey", alpha = 0.3)+
  geom_boxplot(size = 0.8, width = 0.5, fill = "transparent") + labs(x = "Vessel behaviour", y = "Speed (km/h)") +
  scale_y_continuous(trans = "sqrt", breaks = c(0.5, 2, 5, 10, 15, 20)) + plot.theme)

# bootstrap speed
df <- df %>% left_join(speed.boot)
aov.boot <- data.frame();n = 500
for (i in 1:n) {
  dat <- df %>% dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  test <- summary(aov(sqrt(speed)~vess.encounter.stage, data = dat))[[1]]
  row <- data.frame(test = "aov", var = "vess.encounter.stage", iteration = i, f = test[1,4], p = test[1,5])
  aov.boot <- bind_rows(aov.boot, row)
}
write.csv(aov.boot, "intermediate-products/encounter-stage/vess-enc-stage_speed_aov_boot.csv", row.names = FALSE)
nrow(aov.boot %>% filter(p<0.05))/nrow(aov.boot) # all iterations significant


# DI
df = tot.enc %>% filter(!is.na(DI) & !is.na(vess.encounter.stage)) # 415
df %>% group_by(vess.encounter.stage) %>% summarise(n = n(), mean = mean(DI)) # some clear differences ...
df <- df %>% filter(vess.encounter.stage != "departure")
kruskal.test(asin(DI) ~ vess.encounter.stage, data = df)

(p.di <- ggplot(data = df, aes(x = vess.encounter.stage, y = DI)) +
  geom_point(position = position_jitter(width=0.23), size = 2, color = "grey", alpha = 0.3)+
  geom_boxplot(size = 0.8, width = 0.5, fill = "transparent") + labs(x = "", y = "Directness index") +
  scale_y_continuous(trans = "asn", breaks = c(0.05, 0.2, 0.5, 0.7, 0.9, 0.99)) + plot.theme)

# bootstrap DI
df <- df %>% left_join(di.boot)
aov.boot <- data.frame(); n = 500
for (i in 1:n) {
  dat <- df %>% dplyr::select(-DI) %>% rename(DI = paste0("DI",i))
  test <- kruskal.test(asin(DI)~vess.encounter.stage, data = dat) 
  row <- data.frame(test = "kruskal", var = "vess.encounter.stage", iteration = i, 
                    chi.sq = as.numeric(test$statistic), p = test$p.value)
  aov.boot <- bind_rows(aov.boot, row)
}
write.csv(aov.boot, "intermediate-products/encounter-stage/vess-enc-stage_DI_aov_boot.csv", row.names = FALSE)
nrow(aov.boot %>% filter(p<0.05))/nrow(aov.boot) # no iterations significant

# Bring it all together!
plot_grid(p.dive, p.breath, p.speed, p.di, ncol = 2, nrow = 2, align = "hv")
ggsave("intermediate-products/encounter-stage/vess-encounter-stage-vars.png", dpi = 600, height = 7, width = 11, bg = "white")



