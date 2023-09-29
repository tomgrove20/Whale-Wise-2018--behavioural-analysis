###################################################
### BRT VARIABLE TRANSFORMATIONS & COLLINEARITY ###
###################################################
### 07/02/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

# BRTs will use the same explanatory variables, so we only have to explore transformations and collinearity once. BRTs are robust to both but transformations can still reduce skew and removing collinear variables can improve interpretability of the final model

### PACKAGES
packages <- c("tidyverse", "ppcor","RColorBrewer", "scales","Metrics", "corrplot")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

### THEME
source("code/themes.R")
brewer.pal(n=9,"GnBu")
brewer.pal(n=9, "YlOrRd")


############
### DATA ###
############

tot <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv") %>%
    mutate(tot.num.30.1500 = oak.num.30.1500 + rib.num.30.1500,
           tot.dist.30.1500 = oak.dist.30.1500 + rib.dist.30.1500,
           direction.time.a.num = as.numeric(gsub(":.*","",direction.time.a))*30+ as.numeric(gsub(".*:","",direction.time.a))*0.5)


################################
### VARIABLE TRANSFORMATIONS ###
################################
# this will be very manual!

## DIST
ggplot(data = tot, aes(x = dist)) + geom_histogram(binwidth = 100, center = 50) # good/bad/ugly
ggplot(data = tot, aes(x = log(dist))) + geom_histogram() # log transformation actually works super well!!

## CLOCK DIRECTION - not checking

## FOCAL SPEED
ggplot(data = tot, aes(x = focal.vess.speed.60)) + geom_histogram() # great

## FOCAL ACCELERATION
ggplot(data = tot, aes(x = focal.vess.accel.60)) + geom_histogram()

## FOCAL MIN/MAX ACCELERATION
ggplot(data = tot, aes(x = focal.vess.accel.max.60)) + geom_histogram()
ggplot(data = tot, aes(x = focal.vess.accel.max.60)) + geom_histogram() + xlim(0.5,2) # might want to filter out between 1 and 2?
ggplot(data = tot, aes(x = focal.vess.accel.min.60)) + geom_histogram()
ggplot(data = tot, aes(x = focal.vess.accel.min.60)) + geom_histogram() + xlim(-2,-0.5) # might want to filter out between -2 and -0.5?

## FOCAL DI
ggplot(data = tot, aes(x = focal.vess.di.60)) + geom_histogram() # not very good. various transformations didn't work (sqrt, cubert, log)

## AIS NUMBER OF BOATS
ggplot(data = tot, aes(x = tot.num.30.1500)) + geom_bar()

## AIS DISTANCE
ggplot(data = tot, aes(x = tot.dist.30.1500)) + geom_histogram()

## GROUP SIZE
ggplot(data = tot, aes(x = groupsize)) + geom_bar() # nothing we can do 

# for now, transforming slope (log+1), mxl (log), fronts_persist (log+1), sal (exp)
tot <- tot %>% mutate(`log(dist)` = log(dist)) 


####################
### COLLINEARITY ###
####################

## FOCAL VESSEL: using 30 seconds for now
vars <- c("groupsize",  "direction.time.a.num", "`log(dist)`", "focal.vess.speed.60", "focal.vess.accel.60", "focal.vess.di.60", "tot.num.30.1500", "tot.dist.30.1500", "beaufort", "day.of.year", "encounter.minute", "eng.off.on.1sur")

# Pearson's correlation
x <- tot[,gsub("`","",vars)] %>% drop_na()
pcor(x, method = "pearson")

# Correlation plots
corrplot.mixed(cor(x), order = "hclust", tl.pos = "lt", tl.col = "black", upper = "ellipse")
png("intermediate-products/brt/behaviours-vars-collinearity.png", width = 2000, height = 2000)
corrplot.mixed(cor(x), order = "hclust", tl.pos = "lt", tl.col = "black", upper = "ellipse", tl.cex = 3, cl.cex = 3, number.cex = 2.8)
dev.off()

# removing based on 0.8
# the only one with a significant correlation: tot.num.30.1500 and tot.dist.30.1500
