#####################
### AZIMUTH ERROR ###
#####################
# 11/07/2020
# Tom Grove (tom.grove@ed.ac.uk)

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
az <- read.table("intermediate-products/error/error_raw-tot_20190820_saeborg.txt", sep="\t") %>%
  mutate(datetime = as.POSIXct(datetime))


####################
### RANGE FINDER ###
####################
# already extracted (transcribed in follow)

# calculating errors (sqrt and squared to make all positive)
az <- az %>%
  mutate(err.az.rf = rangefinder.bearing - az.known) %>% # simple error (difference)
  mutate(err.az.rf.delay = rangefinder.bearing - az.delay.rf) # simple error with rf.delay (3 seconds) 

# to remove possible Constance instead of Hannah instances, add:
#  filter(notes != "on Constance not Hanna!") %>%

# plotting known azimuth against rangefinder bearing
# plotting with az from vessel location with no correction
ggplot(data = az) +
  geom_point(aes(x = az.known, y = rangefinder.bearing, color = dist.known/1000)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "GPS azimuth (deg.)", y = "RF azimuth (deg.)", color = "GPS distance (km)") +
  scale_color_viridis_c() # some big errors are obvious
ggsave("intermediate-products/error/az-rf-vs-gps_2019.png")

# plotting with vessel position 3 seconds later (to account for delay in recording rf point)
ggplot(data = az) +
  geom_point(aes(x = az.delay.rf, y = rangefinder.bearing, color = dist.known/1000)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "GPS azimuth (deg.)", y = "RF azimuth (deg.)", color = "GPS distance (km)") +
  scale_color_viridis_c() # still some pretty big errors
ggsave("intermediate-products/error/az-rf-vs-gps (+3 sec RF delay)_2019.png")

cor(az$az.known, az$rangefinder.bearing, method = "pearson", use="complete.obs") # 0.960 is pretty good

# plotting error vs dist.known
ggplot(data = az) +
    geom_point(aes(x = dist.known, y = err.az.rf)) + 
  labs(x = "GPS distance (m)", y = "RF azimuth error (deg.)") +
  geom_hline(yintercept = 0)
ggsave("intermediate-products/error/az-rf-error-vs-dist (+3 sec RF delay)_2019.png")
# some pretty big errors, larger at shorter distances. This is likely due to the time delay in getting a rangefinder point but we can't be sure

# Plot changes throughout the session
# normal
ggplot(data = az) +
  geom_point(aes(x = datetime, y = err.az.rf))

# let's try to create a lm and see what happens
lm.err.dist <- lm(err.az.rf ~ dist.known, data = az)
lm.err.dist # pretty small slope which is good! Large intercept though. let's try to remove small distances (perhaps GPS errors)

az.filter <- filter(az, dist.known > 400)

lm.err.dist.400 <- lm(err.az.rf ~ dist.known, data = az.filter)
lm.err.dist.400 # smaller slope (good) with an intercept of 8.5. Let's use this!

cor = as.numeric(coef(lm.err.dist.200)["(Intercept)"]) # define correction factor

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
ggsave("intermediate-products/error/az-rf-cor-vs-gps_2019-moving.png")

# err.az.rf.cor vs dist.known, to confirm what we see in the graph above
ggplot(data = az) +
  geom_point(aes(x = dist.known/1000, y = err.az.rf.cor)) + 
  labs(x = "GPS distance (km)", y = "RF azimuth error (deg.)") +
  geom_hline(yintercept = 0) 
# smaller distances aren't very screwy!
ggsave("intermediate-products/error/az-rf-cor-error-vs-dist_2019-moving.png")

# TOTAL STANDARD ERROR
# no RF delay
sd(az$az.rf.cor - az$az.known, na.rm=TRUE) /  
  sqrt(length(az$az.rf.cor[!is.na(az$az.rf.cor)])) # 1.48 degrees. Not too bad really

# CUMULATIVE SE
# we will also calculate cumulative standard error, in order of distance (both ascending and descending)
# cumulative from large distance to small
az <- arrange(az, desc(dist.known))
for (i in 1:nrow(az)) {
  az$se.az.rf.desc[i] <- sd(az$err.az.rf.cor[1:i], na.rm = TRUE)/sqrt(length(az$err.az.rf.cor[1:i]))
}

# cumulative from small distance to large
az <- arrange(az, dist.known)
for (i in 1:nrow(az)) {
  az$se.az.rf.asc[i] <- sd(az$err.az.rf.cor[1:i], na.rm = TRUE)/sqrt(length(az$err.az.rf.cor[1:i]))
}

# Finally, let's plot cumulative standard error vs distance
# cumulative SE, descending (starting with high)
ggplot(data = az, aes(x = dist.known, y = se.az.rf.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative rangefinder azimuth SE (desc.)")

# cumulative SE, ascending (starting with low)
ggplot(data = az, aes(x = dist.known, y = se.az.rf.asc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative rangefinder azimuth SE (asc.)")

# In this particular instance, makes more sense to use a single value for SE. We would not expect it to change across distances. It is likely larger at short distances due to GPS errors (although RF delay may have an impact)

# CONCLUSION: measured errors are very large for Az at short distances. We should measure errors when both are stationary to compare, since we likely don't have good enough GPS data. If relying solely on these data, use a single value for SE

#############
### PHOTO ###
#############

# first, let's calculate Az between observer and cam.ref

for (i in 1:nrow(az)) {
  # azimuth
  az$az.ref[i] = bearing(c(az$lon.obs[i], az$lat.obs[i]), c(az$lon.ref[i], az$lat.ref[i]))
}

az <- az %>%
  mutate(az.cam = az.ref + ((whale.pixel.x - ref.pixel.x)/img.x*fov),
         err.az.cam = az.cam - az.known)

# second, let's plot az.cam vs az.known. They should be VERY similar!!
ggplot(data = az) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = az.known, y = az.cam, color = dist.known/1000)) +
  labs(x = "GPS azimuth (deg.)", y = "Camera azimuth (deg.)", color = "GPS distance (km)") +
  scale_color_viridis_c() # they are pretty similar!
ggsave("intermediate-products/error/az-cam-vs-gps_2019-moving.png")

cor(az$az.known, az$az.cam, method = "pearson", use="complete.obs") # 0.997. Very good

# plotting error (rangefinder az - actual az)
ggplot(data = az) +
  geom_point(aes(x = dist.known, y = err.az.cam)) + # tiny errors at greater distances. Much larger errors at short distances
  labs(x = "GPS distance (m)", y = "Camera azimuth error (deg.)") +
  geom_hline(yintercept = 0)
ggsave("intermediate-products/error/az-cam-vs-dist_2019-moving.png")

ggplot(data = az) +
  geom_point(aes(x = datetime, y = err.az.cam)) +
  labs(x = "Time", y = "Camera azimuth error")

# TOTAL SE
sd(az$az.cam - az$az.known, na.rm=TRUE) /  
  sqrt(length(az$az.cam[!is.na(az$az.cam)])) # 0.29 degrees. Very small!!

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
  labs(x = "GPS distance (m)", y = "Camera azimuth SE (deg.)")
ggsave("intermediate-products/error/az-cam-se-desc_2019-moving.png") # I don't expect to use this error structure

# cumulative SE, ascending (starting with low)
ggplot(data = az, aes(x = dist.known, y = se.az.cam.asc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative camera azimuth SE (asc.)") 

# So camera azimuths match very well to GPS azimuths, but rangefinder az do not. I think this may be because of high moving speeds - hence errors are so consistently large at low distances. For now, we can safely use camera azimuths. To actually validate RF azimuths, we should carry out tests when both boats (ww and inflatable) are stationary, so that we have accurate locations

# Finally, let's save the azimuth error data
az <- az %>%
  select(datetime, rangefinder.bearing, dist.known, dist.delay.rf, az.known, az.delay.rf, as.numeric(ncol(az)-11):as.numeric(ncol(az)))

write.csv(az, "intermediate-products/error/error_az_20190820_saeborg.csv")

###################
### CONCLUSIONS ###
###################
# Rangefinder: Az errors quite large, variable and unpredictable. Very different to 2020. Unsure whether it's due to GPS errors/timing mismatches OR due to the movement of the boat increasing errors (could be either or both). To resolve this, this should be repeated, ensuring times are exactly the same on all GPS devices.
# Camera: small Az measurement errors
# For a rangefinder, we don't expect SE to vary with distance (although it does seem to be greater at shorter distances, unsure why yet). Thus, a single value of 1.48 degrees will be used as SE for all distances.
# For a camera, SE is tiny and does not vary with distance. Thus, a single value of 0.29 degrees will be used
# The impact of these errors on movement pattern reconstruction will be investigated through resampling of angles from a distribution based on this random error.

