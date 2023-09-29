#########################
### ANALYSIS: LOGGING ###
#########################
### 20/02/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

# a mixture of chi-squared tests and ANOVAs

### PACKAGES
packages <- c("tidyverse","RColorBrewer","nlme")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)


#---------------- FUNCTIONS + THEME --------------------

source("code/functions.R")
source("code/themes.R")
brewer.pal(n=9,"GnBu")
brewer.pal(n=9, "YlOrRd")


#---------------- DATA --------------------

# filtering:
# follows with logging and no surface feeding
# morning (= more likely to reflect actual resting behaviour instead of deep feeding or other behaviours that might involve logging)
# how does dive time and distance travelled vary with final distance and end speed (60)? 
# how does number of breaths (needs to be clear surfacing sequences) vary with start distance (can go +1) and end distance (can go -1) and end speed (can go -1)?
log.filter <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv") %>%
  group_by(folnum.unique) %>%
  # does a follow contain logging and not surface feeding?
  mutate(log = ifelse(grepl("log", behav.codes),1,0),
         log.follow = max(log, na.rm = TRUE),
         surface.feeding.follow = max(surface.feeding, na.rm = TRUE)) %>% 
  filter(log.follow == 1 & surface.feeding.follow == 0) %>%
  # needs to be before 1300
  filter(hour(datetime)<14) %>%
  dplyr::select(folnum.unique, log.follow) %>%
  distinct()# only 20 follows!

# now joining to dive and surface interval data sets
dive <- read.csv("intermediate-products/response-var-dfs/dive-time_ais.csv") %>%
  left_join(log.filter) %>% filter(log.follow == 1) # 34 clear dives

surface <- read.csv("intermediate-products/response-var-dfs/surface-interval_ais.csv") %>%
  left_join(log.filter) %>% filter(log.follow ==1) # 17 surfacing intervals

# how about speed?
speed <- read.csv("intermediate-products/response-var-dfs/speed_ais.csv") %>%
  left_join(log.filter) %>% filter(log.follow == 1, clear.dive == 1) # only 13 observations!
speed.unfiltered <- read.csv("intermediate-products/response-var-dfs/speed_ais.csv") %>% filter(clear.dive == 1)

# last thing, let's attach breath num to dives!
dive <- dive %>%
  left_join(surface %>% rename(datetime = dt.end) %>% dplyr::select(datetime, breath.num))


#---------------- DIVING PATTERNS --------------------

# how does number of breaths relate to start distance and dive.time.before?
mod <- lme(breath.num ~ dist.start + dive.time.before,random=~1|folnum.unique,data=surface)
anova(mod)
summary(mod)

ggplot(data = surface, aes(x = dist.start, y = breath.num)) +  
  geom_smooth(method = "lm", color = "black", alpha = 0.2) +
  geom_point(aes(color = as.factor(folnum.unique)), size = 3, alpha = 0.8) +
  scale_x_continuous(breaks = c(50,seq(100,700, by = 100))) + plot.theme + theme(legend.position = "none") +
  scale_y_continuous(breaks = c(1:14)) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Initial distance (m)", y = "Breaths per surfacing interval")
ggsave("intermediate-products/logging/log-follows_breath.num-vs-dist.start.png")

# how does distance travelled relate to final distance?
ggplot(data = speed, aes(x = dist, y = length)) + geom_point() + geom_smooth(method = "lm")


ggplot(data = dive, aes(x = dist, y = dive.time)) +
  geom_point() + xlim(0,200) + geom_smooth(method = "lm")





m2 <- lme(dive.time ~ dist, random=~1|folnum.unique, data = dive)
anova(m2)
plot(m2)

# conclusion: unlikely to use any of this!
