######################
### DISTANCE ERROR ###
######################
# 11/07/2020
# Tom Grove (tom.grove@ed.ac.uk)

# packages
packages <- c("tidyverse", "ggplot2", "sf", "geosphere", "leaflet", "viridis", "rgeos")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

# ggplot theme
theme_set(theme_classic()) 
theme_update(axis.title.y = element_text(margin = margin(0,20,0,0)),
             axis.title.x = element_text(margin = margin(20,0,0,0)),
             text = element_text(size=15))

# Total raw data
dist <- read.table("intermediate-products/error/error_raw-tot_20190820_saeborg.txt", sep="\t") %>%
  mutate(datetime = as.POSIXct(datetime))

####################
### RANGE FINDER ### 
####################
# already extracted (transcribed in follow)

# calculating errors (sqrt and squared to make all positive)
# plotting known distance against rangefinder distance
rf.dist <- filter(dist, !is.na(rangefinder.distance))

rf.dist <-rf.dist %>%
  mutate(err.dist.rf = rangefinder.distance - dist.known) %>% # simple error (difference)
  mutate(err.dist.rf.delay = rangefinder.distance - dist.delay.rf) # simple error with rf.delay (3 seconds)

# to remove possible Constance instead of Hannah instances, add:
#  filter(notes != "on Constance not Hanna!") %>%

# plotting with dist from vessel location with no correction 
ggplot(data = rf.dist) +
  geom_point(aes(x = dist.known, y = rangefinder.distance)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "GPS distance (m)", y = "RF distance (m)") # pretty good!!
ggsave("intermediate-products/error/dist-rf-vs-gps_2019.png")

# +3 seconds RF correction
ggplot(data = rf.dist) +
  geom_point(aes(x = dist.delay.rf, y = rangefinder.distance)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "GPS distance", y = "Rangefinder distance") +
  scale_color_viridis_c() # also pretty good
ggsave("intermediate-products/error/dist-rf-vs-gps (+3 sec RF delay__2019.png")


# so, it seems that distance estimates are more accurate without rangefinder deployment correction of 3 seconds

# Look at Pearson linear correlation
cor(rf.dist$dist.known, rf.dist$rangefinder.distance, method = "pearson") # 0.998 is excellent!

# plotting error vs dist.known
ggplot(data = rf.dist) +
  geom_point(aes(x = dist.known, y = err.dist.rf)) + 
  labs(x = "GPS distance (m)", y = "Rangefinder distance error") +
  geom_hline(yintercept = 0)
# pretty small errors, evenly spread across distances. Comparatively larger at short distances
ggsave("intermediate-products/error/dist-rf-error-vs-dist_2019.png")

# plotting delayed error vs dist.known
ggplot(data = rf.dist) +
  geom_point(aes(x = dist.known, y = err.dist.rf.delay)) + 
  labs(x = "GPS distance (m)", y = "Rangefinder distance error") +
  geom_hline(yintercept = 0)
# pretty small errors, evenly spread across distances. Comparatively larger at short distances
ggsave("intermediate-products/error/dist-rf-error-vs-dist.png")

# Plot changes throughout the session
# normal
ggplot(data = rf.dist) +
  geom_point(aes(x = datetime, y = err.dist.rf)) # largest errors were when we were moving directly away from the boat (makes a lot of sense)

# + 3 second rf delay
ggplot(data = rf.dist) +
  geom_point(aes(x = datetime, y = err.dist.rf.delay)) # rf delay clearly doesn't work here. Far larger errors when moving directly away from the boat. RF delay doesn't seem to apply in the case of distance.

# TOTAL STANDARD ERROR
# no RF delay
sd(rf.dist$rangefinder.distance - rf.dist$dist.known, na.rm=TRUE) /  
  sqrt(length(rf.dist$rangefinder.distance[!is.na(rf.dist$rangefinder.distance)])) # 0.31 m. Pretty good!

# +3 seconds RF delay
sd(rf.dist$rangefinder.distance - rf.dist$dist.delay.rf, na.rm=TRUE) /  
  sqrt(length(rf.dist$rangefinder.distance[!is.na(rf.dist$rangefinder.distance)])) # 0.43, slightly higher!

# we will also calculate cumulative standard error, in order of distance (both ascending and descending)
# cumulative from large distance to small
rf.dist <- arrange(rf.dist, desc(dist.known))
for (i in 1:nrow(rf.dist)) {
  rf.dist$se.dist.rf.desc[i] <- sd(rf.dist$err.dist.rf[1:i], na.rm = TRUE)/sqrt(length(rf.dist$err.dist.rf[1:i]))
}

# cumulative from small distance to large
rf.dist <- arrange(rf.dist, dist.known)
for (i in 1:nrow(rf.dist)) {
  rf.dist$se.dist.rf.asc[i] <- sd(rf.dist$err.dist.rf[1:i], na.rm = TRUE)/sqrt(length(rf.dist$err.dist.rf[1:i]))
}

# Finally, let's plot cumulative standard error vs distance
# cumulative SE, descending (starting with high)
ggplot(data = rf.dist, aes(x = dist.known, y = se.dist.rf.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative rangefinder distance SE (desc.)")

# cumulative SE, ascending (starting with low)
ggplot(data = rf.dist, aes(x = dist.known, y = se.dist.rf.asc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative rangefinder distance SE (asc.)")
# In this particular instance, doesn't really make sense to use distance-specific errors. For error propagation, we should use total SE value. 

##############
### CAMERA ###
##############

# To calculate distance from camera, we need the following information
# - Photo dimension
# - Vertical pixel distance between whale and horizon
# - True/false horizon. If false, distance between observer and horizon
# - Height of observer above water

# First, we need to upload an Iceland shapefile (to calculate obs-false horizon distances)
land <- read_sf("intermediate-products/error/iceland_combined-features/iceland_combined-features.shp")
land
plot(land)

# now we need to calculate the distance between the observer and the land, at a specific bearing. To do this, we will draw a 100-km line in the direction of az.known

cam.dist.calc <- dist %>%
  filter(!is.na(horizon.distance.pixel)) %>% # we only want rows with horizon pixel distance
  arrange(datetime)
 
# creating end point of a line 100 km long in the desired bearing 
for(i in 1:nrow(cam.dist.calc)) {
  end <- destPoint(c(cam.dist.calc$lon.obs[i], cam.dist.calc$lat.obs[i]), b = cam.dist.calc$az.known[i], d = 100000)
  cam.dist.calc$lon.end[i] <- end[1,1]
  cam.dist.calc$lat.end[i] <- end[1,2]
}

# creating sf 'lines', containing full 100 km lines in specified bearing
lines <- cam.dist.calc %>%
  unite(geom.begin, lon.obs, lat.obs) %>% # collect start coords into one column for reshaping
  unite(geom.end, lon.end, lat.end) %>% # collect end coords into one column for reshaping
  gather(start_end, coords, geom.begin, geom.end) %>% # reshape to long
  separate(coords, c("LONG", "LAT"), sep = "_") %>% # convert our text coordinates back to individual numeric columns
  mutate_at(vars(LONG, LAT), as.numeric) %>%
  st_as_sf(coords = c("LONG", "LAT"), crs = 4326) %>% # create points
  st_transform(crs = 3057) %>% # using ISN93 CRS as it's safer to use planar coordinate systems when intersecting features
  group_by(datetime) %>% # make sure only datetime pairs are joined
  arrange(start_end, datetime) %>% # important that it's also arranged by datetime for later!
  summarise(do_union = FALSE) %>% # union points into lines. do_union = FALSE keeps points in the order you want
  st_cast("LINESTRING") %>% # cast as linestring (sf class)
  rename(geom.line = geometry) %>%
  ungroup
lines # still ISN93 (EPSG: 3057)
plot(lines)

# let's save this for plotting later
st_write(lines, "intermediate-products/error/intersect-lines-uncut/20190820_saeborg_error_intersect-lines_uncut.shp", append = FALSE)

# now let's look for the point of intersection between each line and the Icelandic land, and calculate the length of the resulting line

shortlines <- lines %>%
  as.data.frame %>%
  select(datetime) %>%
  arrange(datetime)


for (i in 1:nrow(shortlines)) {
  
  intersect <- st_cast(st_intersection(lines$geom.line[i], land), "POINT") # important to re cast the intersection product from MULTILINESTRING (line chopped up by intersections) to LINESTRING (separate line for each segment).
  
  if(length(intersect) == 0) { # if there is no intersection (true horizon!)
    shortlines$dist.hor.false[i] = NA
  } else {
    
    startpt <- st_cast(st_line_sample(lines$geom.line[i], sample = 0), "POINT") # recreate start point (first point on geom.line) 
    
    endpt <- intersect[which.min(st_distance(startpt, intersect))] # find which point from the intersecting lines corresponded with the minimum distance
    
    while(as.numeric(st_distance(startpt, endpt)) < 1500) { # horizon must be further than 1500 m away (arbitrary)
      intersect <- intersect[-1]
      endpt <- intersect[which.min(st_distance(startpt, intersect))] # find which point from the intersecting lines corresponded with the minimum distance
    }
    
    shortlines$dist.hor.false[i] = as.numeric(st_distance(startpt, endpt)) # minimum distance between intersected lines and start point. as.numeric removes the stupid units
    
    shortlines$geom.shortline[i] = st_cast(st_union(startpt, endpt), "LINESTRING") # creating truncated line between vessel location and intersected point
  }
}

shortlines <- st_as_sf(shortlines, crs = 3057)

# plotting this with Icelandic land
ggplot() +
  geom_sf(data = shortlines, color = "black") +
  geom_sf(data = land) +
  geom_sf(data = st_cast(st_combine(vess), "MULTILINESTRING"), color = "red") +
  coord_sf(xlim = c(min(st_coordinates(shortlines)[,1]), max(st_coordinates(shortlines)[,1])), ylim = c(min(st_coordinates(shortlines)[,2]), max(st_coordinates(shortlines)[,2]))) +
  theme_bw()  +
  theme(axis.text.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(10,0,0,0)),
        text = element_text(size=15))
ggsave("intermediate-products/error/20190820_error-hor-lines.png")
# looks like it worked!

st_write(shortlines, "intermediate-products/error/intersect-lines-cut/20190820_saeborg_error_intersect-lines-cut.shp", append = FALSE)

# OK, so this seems to have worked well. Now we need to calculate the distance

shortlines <- as.data.frame(shortlines) %>% # convert from sf back to normal df
  select(-geom.shortline)

cam.dist.calc <- cam.dist.calc %>% 
  left_join(shortlines, by = "datetime") # joining our distance-to-coast values (in m) to our full data set

# defining a few terms
sensor.y <- 14.9 # sensor height in mm for canon 77D
crop.factor <- 1 # The image isn't really being 'cropped', so no crop factor needed here IMO
re <- 6360 # radius of the earth (in km)

# calculating distance. Note: we use simpler TRUE horizon equations where possible. TRUE/FALSE based on whether we have distance to land value
cam.dist.df <- cam.dist.calc %>%
  mutate(
    # convert dist.hor.false from m to km
    dist.hor.false = dist.hor.false/1000,
    
    # total height(in km, so /1000)
    h = ifelse(stand.sit.a == "stand", ((height.plat.a/100) + stand.a)/1000, ((height.plat.a/100) + sit.a)/1000),
    
    # now look at dist.hor.true
    dist.hor.true = horizon(h*1000, r = re*1000)/1000,
    
    # theta = vertical angle between the whale and the horizon 
    theta = atan(((horizon.distance.pixel/img.y) * sensor.y)/(fl*crop.factor)),
    
    # alpha = angle above horizon to horizontal tangent. For now, only calculating for true horizon. use is.na(dist.hor.false) instead of horizontype == "t" for now. 
    alpha = ifelse(is.na(dist.hor.false)|dist.hor.false > dist.hor.true, atan(sqrt((2*re*h)+h^2)/re), NA),
    
    # gamma = angle between observer and false horizon. Remember to convert dist.hor to km (/1000)
    gamma = ifelse(is.na(alpha), dist.hor.false/re, NA),
    
    # l0 = straight-line distance between observer and horizon (not going around the surface of the earth), based on the cosine rule
    l0 = ifelse(is.na(alpha),sqrt((re^2) + ((re+h)^2) - (2*re*(re + h)*cos(gamma))), NA),
    
    # from here, calculate beta
    beta = ifelse(is.na(alpha), acos(((2*h*re) + (h^2) + (l0^2))/ (2*(re + h)*l0)) - theta, NA),
    
    # now calculate A, the angle between the whale and the horizontal (from the observer)
    A = ifelse(!is.na(alpha), alpha + theta, (pi/2) - beta),
    
    # cam.dist = camera distance, only for true horizon
    cam.dist = (((re + h)*sin(A)) - sqrt((re^2) - (((re + h)*cos(A))^2)))*1000)
                
# calculating errors (sqrt and squared to make all positive)
cam.dist.df <- cam.dist.df %>%
  mutate(err.dist.cam = cam.dist - dist.known) # simple error (difference)

# plotting with dist from vessel location with no correction 
ggplot(data = cam.dist.df) +
  geom_point(aes(x = dist.known, y = cam.dist)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "GPS distance (m)", y = "Camera distance (m)") +
  scale_color_viridis_c()  # pretty terrible!
  # xlim(0,1000) + ylim(0, 1000)
ggsave("intermediate-products/error/dist-cam-vs-gps_2019.png")

# Look at Pearson linear correlation
cor(test$dist.known, test$cam.dist, method = "pearson") # 0.9995 is great!

# plotting error vs dist.known
ggplot(data = cam.dist.df) +
  geom_point(aes(x = dist.known, y = err.dist.cam)) + 
  labs(x = "GPS distance (m)", y = "Camaera distance error") +
  geom_hline(yintercept = 0)
# Quite large errors at large distances
ggsave("intermediate-products/error/dist-cam-error-vs-dist_2019.png")

# Plot changes throughout the session
# normal
ggplot(data = cam.dist.df) +
  geom_point(aes(x = datetime, y = err.dist.cam)) # largest errors were when we were moving directly away from the boat (makes a lot of sense)

# TOTAL SE
sd(cam.dist.df$cam.dist - cam.dist.df$dist.known, na.rm=TRUE) /  
  sqrt(length(cam.dist.df$cam.dist[!is.na(cam.dist.df$cam.dist)])) # 2.5 m. Pretty good!

# we will also calculate cumulative standard error, in order of distance (both ascending and descending)
# cumulative from large distance to small
cam.dist.df <- arrange(cam.dist.df, desc(dist.known))
for (i in 1:nrow(cam.dist.df)) {
  cam.dist.df$se.dist.cam.desc[i] <- sd(cam.dist.df$err.dist.cam[1:i], na.rm = TRUE)/sqrt(length(cam.dist.df$err.dist.cam[1:i]))
}

# cumulative from small distance to large
cam.dist.df <- arrange(cam.dist.df, dist.known)
for (i in 1:nrow(cam.dist.df)) {
  cam.dist.df$se.dist.cam.asc[i] <- sd(cam.dist.df$err.dist.cam[1:i], na.rm = TRUE)/sqrt(length(cam.dist.df$err.dist.cam[1:i]))
}

# Finally, let's plot cumulative standard error vs distance
# cumulative SE, descending (starting with high)
ggplot(data = cam.dist.df, aes(x = dist.known, y = se.dist.cam.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Camera distance SE (m)") + xlim(0,1650) + ylim(0,10) +
  theme(axis.ticks.length = unit(0.5, "cm"),
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
        axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")))
ggsave("intermediate-products/error/dist-cam-cum-se-desc_2019-moving.png")

# cumulative SE, ascending (starting with low)
ggplot(data = cam.dist.df, aes(x = dist.known, y = se.dist.cam.asc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative rangefinder distance SE (asc.)")

# In this particular instance, descending makes far more sense than ascending due to larger errors at greater distances. Errors definitely vary with distance so we will create a linear model to extrapolate errors from this structure. No need to filter out large distances (unlike in 2020)

lm.dist.cam <- lm(se.dist.cam.desc ~ dist.known, data = cam.dist.df) # creating linear model
cam.dist.df$cam.dist.se.pred <- predict(lm.dist.cam, newdata = cam.dist.df) # predicted values based on known distances from the lm

ggplot(data = cam.dist.df, aes(x = dist.known, y = se.dist.cam.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Camera distance SE (m)") +
  xlim(0,1600) + ylim(0,10) +
  geom_line(aes(y = cam.dist.se.pred), color="red") +
  theme(axis.ticks.length = unit(0.5, "cm"),
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
        axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")))
ggsave("intermediate-products/error/dist-cam-cum-se-desc+model-line_2019-moving.png")

# Finally saving the dist data

cam.dist.df <- left_join(cam.dist.df, rf.dist)
write.csv(cam.dist.df, "intermediate-products/error/error_dist_20190820_saeborg.csv")

###################
### CONCLUSIONS ###
###################
# Overall, distance estimates produced by both rangefinder and photographic methods produced acceptably small errors (which now need to be propagated through the movement pattern models)
# Rangefinder: very small errors, don't vary much with distance. A single SE value of 0.31 m will be used
# Camera: larger errors but not too bad. Seem quite a bit larger (and more positive) than 2020. This makes sense as the inflatable was moving in 2019, so pixel distances will be smaller and calculated distances greater. Further, it was a rougher day and waves seemed to obscure boat sometimes 
# Errors vary greatly with distance, so a linear model based on cumulative SE (descending) will be used to interpolate values. Not recommending exclusion of values yet

