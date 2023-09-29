##########################
### AZIMUTH ERROR 2020 ###
##########################
# 13/09/2020
# Tom Grove (tom.grove@ed.ac.uk)

# Calculating the error of two azimuth measurement errors from whale watching vessels: using a rangefinder's electronic compass (Az function) or using a camera with reference landmarks 

# packages
packages <- c("tidyverse", "ggplot2", "sf", "geosphere", "lwgeom", "leaflet", "viridis")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

# ggplot theme
theme_set(theme_classic()) 
theme_update(axis.title.y = element_text(margin = margin(0,20,0,0)),
             axis.title.x = element_text(margin = margin(20,0,0,0)),
             text = element_text(size=15))

# Total raw data
az <- read.table("intermediate-products/error/error_raw-tot_202008_saeborg.txt", sep="\t") %>%
  mutate(datetime = as.POSIXct(datetime))

####################
### RANGE FINDER ###
####################
# already extracted (transcribed in follow)

# first, calculate errors from unchanged rf bearing (but with magnetic declination subtracted)
az <- az %>%
  mutate(err.az.rf = rangefinder.bearing - az.known) # simple error (difference)

# now, let's plot rangefinder.bearing vs az.known
ggplot(data = az) +
  geom_point(aes(x = az.known, y = rangefinder.bearing, color = dist.known/1000)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "GPS azimuth (deg.)", y = "RF azimuth (deg.)", color = "GPS distance (km)") +
  scale_color_viridis_c() # it seems we need a large correction factor (likely due to inaccuracy of rangefinder calibration, to be expected)
ggsave("intermediate-products/error/az-rf-unchanged-vs-gps_2020-stationary.png")

# we will also plot err.az.rf vs dist.known, to confirm what we see in the graph above
ggplot(data = az) +
  geom_point(aes(x = dist.known/1000, y = err.az.rf)) + 
  labs(x = "GPS distance (km)", y = "RF azimuth error (deg.)") +
  geom_hline(yintercept = 0)
# totally confirms what we thought above. Correction factor clearly needed
ggsave("intermediate-products/error/az-rf-unchanged-error-vs-dist_2020-stationary.png")

# to note, there seems to be something going on with short distances. Possibly because rangefinder distances were measured here, which takes a few seconds. At short distances, with a boat drifting, this can make a difference. This should be investigated later (see if we can get near to 15-16 degrees off?)

cor(az$az.known, az$rangefinder.bearing, method = "pearson", use="complete.obs") # 0.983 is pretty good (shows a clear relationship despite and correction needed)

# good correlation but clearly a poor intercept (likely because calibrating the rangefinder is pretty imprecise). To account for this, we will investigate the relationship between rf.az.err and dist.known. The intercept will be the correction factor (with opposite sign). 

lm.err.dist <- lm(err.az.rf ~ dist.known, data = az)
lm.err.dist

# slope (dist.known coeff) is close to zero, which is good.
# intercept is 15.8 = correction factor. We can now create a new field, az.rf.cor

cor = as.numeric(coef(lm.err.dist)["(Intercept)"]) # define correction factor

az <- az %>%
  mutate(az.rf.cor = rangefinder.bearing - cor, # correcting by cor
         err.az.rf.cor = az.rf.cor - az.known)

# now we can re-do the plots above
# az.rf.cor vs az.known
ggplot(data = az) +
  geom_point(aes(x = az.known, y = az.rf.cor, color = dist.known/1000)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "GPS azimuth (deg.)", y = "RF azimuth (deg.)", color = "GPS distance (km)") +
  scale_color_viridis_c() # it seems we need a large correction factor (likely due to inaccuracy of rangefinder calibration, to be expected)
ggsave("intermediate-products/error/az-rf-cor-vs-gps_2020-stationary.png")

# err.az.rf.cor vs dist.known, to confirm what we see in the graph above
ggplot(data = az) +
  geom_point(aes(x = dist.known/1000, y = err.az.rf.cor)) + 
  labs(x = "GPS distance (km)", y = "RF azimuth error (deg.)") +
  geom_hline(yintercept = 0) 
# smaller distances aren't very screwy!
ggsave("intermediate-products/error/az-rf-cor-error-vs-dist_2020-stationary.png")

# TOTAL STANDARD ERROR
# no RF delay
sd(az$err.az.rf.cor, na.rm=TRUE) /  
  sqrt(length(az$az.rf.cor[!is.na(az$az.rf.cor)])) # 0.280 degrees. Very good!

# CUMULATIVE STANDARD ERROR
# Finally, let's plot cumulative standard error vs distance

# cumulative SE, descending (starting with high)
az <- arrange(az, desc(dist.known))
for (i in 1:nrow(az)) {
  az$se.az.rf.desc[i] <- sd(az$err.az.rf.cor[1:i], na.rm = TRUE)/sqrt(length(az$err.az.rf.cor[1:i]))
}

# cumulative from small distance to large (ascending)
az <- arrange(az, dist.known)
for (i in 1:nrow(az)) {
  az$se.az.rf.asc[i] <- sd(az$err.az.rf.cor[1:i], na.rm = TRUE)/sqrt(length(az$err.az.rf.cor[1:i]))
}

# plot descending
ggplot(data = az, aes(x = dist.known, y = se.az.rf.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative rangefinder azimuth SE (desc.)")

# cumulative SE, ascending (starting with low)
ggplot(data = az, aes(x = dist.known, y = se.az.rf.asc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative rangefinder azimuth SE (asc.)")

# we will use descending cumulative SE. Let's plot up to 4000 m and save it
ggplot(data = az, aes(x = dist.known/1000, y = se.az.rf.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (km)", y = "RF azimuth SE (deg.)") +
    xlim(0,4.500) + ylim(0,1)  +
  theme(axis.ticks.length = unit(0.5, "cm"),
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
        axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")))
ggsave("intermediate-products/error/az-rf-cum-se-desc_2020-stationary.png")
# this includes 98/104 points

# this looks linear, so let's create a linear model (which we can use to extrapolate)
lm.az.rf <- lm(se.az.rf.desc ~ dist.known, data = az)
az$az.rf.se.pred <- predict(lm.az.rf, newdata = az) # predicted values based on known distances from the lm

# plot what this linear model looks like
ggplot(az, aes(x = dist.known, y=se.az.rf.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "RF azimuth SE (deg.)") +
  xlim(0,4500) + ylim(0,1) +
  geom_line(aes(y = az.rf.se.pred), color="red")  +
  theme(axis.ticks.length = unit(0.5, "cm"),
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
        axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")))
ggsave("intermediate-products/error/az-rf-cum-se-desc+model-line_2020-stationary.png")
# this includes 98/104 points

# CONCLUSION: pretty small errors. Some larger errors at small distances, which can likely be attributed to errors in the GPS positions themselves (3-4 m, makes a large difference < 100 m)

#############
### PHOTO ###
#############

# first, let's calculate Az between observer and cam.ref

for (i in 1:nrow(az)) {
  # azimuth
  az$az.ref[i] = (bearing(c(az$lon.obs[i], az$lat.obs[i]), c(az$lon.ref[i], az$lat.ref[i])) + 360) %% 360
}

az <- az %>%
  mutate(az.cam = az.ref + ((whale.pixel.x - ref.pixel.x)/img.x*fov),
         err.az.cam = az.cam - az.known)

# let's now plot az.cam vs az.known. They should be VERY similar!!
ggplot(data = az) +
  geom_point(aes(x = az.known, y = az.cam, color = dist.known/1000)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "GPS azimuth (deg.)", y = "Camera azimuth (deg.)", color = "GPS distance (km)") +
  scale_color_viridis_c() # they are pretty similar!
ggsave("intermediate-products/error/az-cam-vs-gps_2020.png")
# n = 72 points

cor(az$az.known, az$az.cam, method = "pearson", use = "complete.obs") # 0.999. Pretty good

# plotting error (rangefinder az - actual az)
ggplot(data = az) +
  geom_point(aes(x = dist.known/1000, y = err.az.cam)) + 
  labs(x = "GPS distance (km)", y = "Camera azimuth error (deg.)") +
  geom_hline(yintercept = 0) # tiny errors at greater distances. Much larger errors at short distances
ggsave("intermediate-products/error/az-cam-err-vs-dist_2020.png")
# very clearly, this is due to GPS errors at short distances. n = 72 points

# TOTAL SE
# calculating standard error
sd(az$err.az.cam, na.rm=TRUE) /  
  sqrt(length(az$az.cam[!is.na(az$az.cam)])) # 0.081 degrees. Very small!!

# CUMULATIVE SE
# cumulative from large distance to small
az <- arrange(az, desc(dist.known))
for (i in 1:nrow(az)) {
  az$se.az.cam.desc[i] <- sd(az$err.az.cam[1:i], na.rm = TRUE)/sqrt(length(az$err.az.cam[1:i]))
}

# cumulative from small distance to large
az <- arrange(az, dist.known)
for (i in 1:nrow(az)) {
  az$se.az.cam.asc[i] <- sd(az$err.az.cam[1:i], na.rm = TRUE)/sqrt(length(az$err.az.cam[1:i]))
}

# Finally, let's plot cumulative standard error vs distance
# cumulative SE, descending (starting with high)
ggplot(data = az, aes(x = dist.known, y = se.az.cam.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative camera azimuth SE (desc.)")
# very little change but it actually increases (as we would expect from previous plots)

# cumulative SE, ascending (starting with low)
ggplot(data = az, aes(x = dist.known, y = se.az.cam.asc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative camera azimuth SE (asc.)") # again, much higher at lower distances


# So camera azimuths match very well to GPS azimuths, as we would expect (we would be very concerned otherwise). Moreover, SE doesn't vary much with GPS distance. Therefore, we will use the total SE value for each distance. Whilst errors appear slightly greater at short distances, we will likely account for this via an error term for observer GPS position

# Finally, let's save the azimuth error data
az <- az %>%
  select(datetime, rangefinder.bearing, dist.known, az.known, as.numeric(ncol(az)-10):as.numeric(ncol(az)))

write.csv(az, "intermediate-products/error/error_az_2020_saeborg.csv")

###################
### CONCLUSIONS ###
###################
# Using both a range finder and a camera, azimuth measurement errors are generally small, with little difference between distances
# For a rangefinder, the SE for a specific angle will depend on the distance from the whale. This will be derived from the cumulative SE function (starting from high distances and cumulatively including shorter distances)
# For a camera, SE is tiny and does not vary with distance. Thus, a single value of 0.0068 degrees will be used
# The impact of these errors on movement pattern reconstruction will be investigated through resampling of angles from a distribution based on this random error.

