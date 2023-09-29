########################################
### ANALYSIS: INTERVAL APPROACH TYPE ###
########################################
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
brewer.pal(n=9,"GnBu")
brewer.pal(n=9, "YlOrRd")


#---------------- DATA --------------------

tot <- read.csv("intermediate-products/response-var-dfs/interval-approach_tot.csv")


#---------------- APPROACH ANGLE VS THINGS --------------------

# ANOVAs(Kruskal-Wallis) and chi-squared tests

# approach type vs fluke: chi-squared
df <- tot %>% 
  mutate(fluke.last = ifelse(fluke.last >=1, "Fluke", "No fluke"), 
         approach.type = recode(approach.type, back = "other", side = "other")) %>%
  group_by(approach.type, fluke.last) %>% summarise(n = n()) %>%
  spread(key = "fluke.last", value = "n") %>% ungroup() %>% as.data.frame() %>%
  column_to_rownames(var = "approach.type")
chisq.test(df) # not significant
# X-squared = 1.9294, df = 1, p-value = 0.1648
#       Fluke No fluke
# back     11        1
# front    70       16
# side     27        3

# approach type vs number of breaths
ggplot(data = tot, aes(x = dist.end, y = breath.num, group = approach.type, color = approach.type)) + geom_point()
kruskal.test(breath.num~approach.type, data = tot)
summary(lm(breath.num ~ approach.type + dist.end, data = tot))

# approach type vs DI
ggplot(data = tot, aes(x = approach.type, y = asin(DI.last))) + geom_boxplot()
summary(aov(asin(DI.last)~approach.type, data = tot))

# approach.type vs dive time after
ggplot(data = tot, aes(x = approach.type, y = dive.time.next)) + geom_boxplot()
summary(aov(log(dive.time.next)~approach.type, data = tot))
# Df Sum Sq Mean Sq F value Pr(>F)
# approach.type   2   0.19  0.0959   0.178  0.837