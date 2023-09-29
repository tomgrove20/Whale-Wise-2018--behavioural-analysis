###########################
### DISTANCE ERROR 2020 ###
###########################
# 13/09/2020
# Tom Grove (tom.grove@ed.ac.uk)

# Calculating the error of two distance measurement errors from whale watching vessels: using a laser rangefinder or a camera with horizon distances 

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
dist <- read.table("intermediate-products/error/error_raw-tot_202008_saeborg.txt", sep="\t") %>%
  mutate(datetime = as.POSIXct(datetime))

####################
### RANGE FINDER ### 
####################
# already extracted (transcribed in follow)
# note: for a few points in the 20200818 session, we added 2.5 m because a shot was likely taken of the bow and not the actual GPS position. You can compare with transcription

# calculating errors (sqrt and squared to make all positive)
dist <-dist %>%
  mutate(err.dist.rf = rangefinder.distance - dist.known) # simple error (difference)

# plotting with dist from vessel location
ggplot(data = filter(dist, !is.na(rangefinder.distance))) + # filter show only RF range of distances is shown
  geom_point(aes(x = dist.known, y = rangefinder.distance)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "GPS distance (m)", y = "RF distance (m)") # pretty good!!
ggsave("intermediate-products/error/dist-rf-vs-gps_2020-stationary.png")

# Look at Pearson linear correlation
cor(dist$dist.known, dist$rangefinder.distance, method = "pearson", use="complete.obs") # 0.999 is excellent!

# plotting error vs dist.known
ggplot(data = filter(dist, !is.na(rangefinder.distance))) +
  geom_point(aes(x = dist.known, y = err.dist.rf)) + 
  labs(x = "GPS distance (m)", y = "RF distance error (m)") +
  geom_hline(yintercept = 0)
# pretty small errors, evenly spread across distances. Comparatively larger at short distances
ggsave("intermediate-products/error/dist-rf-error-vs-dist_2020-stationary.png")

# TOTAL STANDARD ERROR
# no RF delay
sd(dist$rangefinder.distance - dist$dist.known, na.rm=TRUE) /  
  sqrt(length(dist$rangefinder.distance[!is.na(dist$rangefinder.distance)])) # 0.58 m. Pretty good!

# CUMULATIVE
# we will also calculate cumulative standard error, in order of distance (both ascending and descending)
# cumulative from large distance to small
dist <- arrange(dist, desc(dist.known))
for (i in 1:nrow(dist)) {
  dist$se.dist.rf.desc[i] <- sd(dist$err.dist.rf[1:i], na.rm = TRUE)/sqrt(length(dist$err.dist.rf[1:i]))
}

# cumulative from small distance to large
dist <- arrange(dist, dist.known)
for (i in 1:nrow(dist)) {
  dist$se.dist.rf.asc[i] <- sd(dist$err.dist.rf[1:i], na.rm = TRUE)/sqrt(length(dist$err.dist.rf[1:i]))
}

# Finally, let's plot cumulative standard error vs distance
# cumulative SE, descending (starting with high)
ggplot(data = filter(dist, !is.na(rangefinder.distance)), 
       aes(x = dist.known, y = se.dist.rf.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "RF distance SE (desc.)")

# cumulative SE, ascending (starting with low)
ggplot(data = filter(dist, !is.na(rangefinder.distance)), 
       aes(x = dist.known, y = se.dist.rf.asc)) +
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
land <- read_sf("data/iceland-isn93/is50v_strandlina_flakar_24122017.shp") %>% # keeping as ISN93 for now
  rename(geom.land = geometry) %>%
  select(LEGA, geom.land) %>%
  st_union() # an important step to prevent the creation of separate line features when finding intersections later
land
plot(land)

# this could come in handy in the future (all features combined), so let's save it!
st_write(land, "intermediate-products/error/iceland_combined-features/iceland_combined-features.shp", append = FALSE)

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
st_write(lines, "intermediate-products/error/intersect-lines-uncut/202008_saeborg_error_intersect-lines_uncut.shp", append = FALSE)

# now let's look for the point of intersection between each line and the Icelandic land, and calculate the length of the resulting line. 

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
  geom_sf(data = land) +
  geom_sf(data = st_cast(st_combine(vess), "MULTILINESTRING"), color = "red") +
  geom_sf(data = shortlines, color = "black") +
  coord_sf(xlim = c(min(st_coordinates(shortlines)[,1]), max(st_coordinates(shortlines)[,1])), ylim = c(min(st_coordinates(shortlines)[,2]), max(st_coordinates(shortlines)[,2]))) +
  theme_bw() +
  theme(axis.text.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(10,0,0,0)),
        text = element_text(size=15))
ggsave("intermediate-products/error/202008_error-hor-lines.png")
# looks like it worked!

st_write(shortlines, "intermediate-products/error/intersect-lines-cut_2020/intersect-lines-cut_2020.shp", append = FALSE)

# OK, so this seems to have worked well. Now we need to calculate the distance
# As part of this, for each row, we will also calculate the distance to the true horizon (dist.hor.true), which may often be lower than the apparent distance to the false horizon. If it lower than the false horizon, we adopt the true horizon method. 

shortlines <- as.data.frame(shortlines) %>% # convert from sf back to normal df
  select(-geom.shortline)

cam.dist.calc <- cam.dist.calc %>% 
  left_join(shortlines, by = "datetime") # joining our distance-to-coast values (in m) to our full data set

# defining a few terms
sensor.y <- 15 # sensor height in mm for canon 77D
crop.factor <- 1 # The image isn't really being 'cropped', so no crop factor needed here IMO
re <- 6360 # radius of the earth (in km) at 66 degrees N

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
  geom_point(aes(x = dist.known/1000, y = cam.dist/1000)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "GPS distance (km)", y = "Camera distance (km)") +
  scale_color_viridis_c()  # pretty terrible!
  # xlim(0,1000) + ylim(0, 1000)
ggsave("intermediate-products/error/dist-cam-vs-gps_2020-stationary.png")

# Look at Pearson linear correlation
cor(cam.dist.df$dist.known, cam.dist.df$cam.dist, method = "pearson", use="complete.obs") # 0.9995 is excellent!

# plotting error vs dist.known
ggplot(data = cam.dist.df) +
  geom_point(aes(x = dist.known/1000, y = err.dist.cam)) + 
  labs(x = "GPS distance (km)", y = "Camera distance error (m)") +
  geom_hline(yintercept = 0)
# Quite large errors at large distances
ggsave("intermediate-products/error/dist-cam-error-vs-dist_full_2020-stationary.png")

# plotting error vs dist.known
ggplot(data = cam.dist.df) +
  geom_point(aes(x = dist.known, y = err.dist.cam)) + 
  labs(x = "GPS distance (m)", y = "Camera distance error (m)") +
  geom_hline(yintercept = 0) +
  xlim(0,1000) +
  ylim(-50,50)
# Quite large errors at large distances
ggsave("intermediate-products/error/dist-cam-error-vs-dist_to-1000_2020-stationary.png")

# TOTAL STANDARD ERROR
# all distances
sd(cam.dist.df$cam.dist - cam.dist.df$dist.known, na.rm=TRUE) /  
  sqrt(length(cam.dist.df$cam.dist[!is.na(cam.dist.df$cam.dist)])) # 9.24 m. Not amazing but does include large distances but huge errors

# filter distances < 2600 m
cam.dist.filter <- cam.dist.df %>%
  select(cam.dist, dist.known, err.dist.cam) %>%
  filter(dist.known < 2600)

sd(cam.dist.filter$cam.dist - cam.dist.filter$dist.known, na.rm=TRUE) /  
  sqrt(length(cam.dist.filter$cam.dist[!is.na(cam.dist.filter$cam.dist)])) # 1.93 m. Pretty good!

# so from this, 2000 m seems a reasonable upper limit. Let's use cam.dist.filter with a distance of <2600 m to calculate cumulative standard error

# CUMULATIVE STANDARD ERROR
# we will also calculate cumulative standard error, in order of distance (both ascending and descending)
# cumulative from large distance to small
cam.dist.filter <- arrange(cam.dist.filter, desc(dist.known))
for (i in 1:nrow(cam.dist.filter)) {
  cam.dist.filter$se.dist.cam.desc[i] <- sd(cam.dist.filter$err.dist.cam[1:i], na.rm = TRUE)/sqrt(length(cam.dist.filter$err.dist.cam[1:i]))
}

# cumulative from small distance to large
cam.dist.filter <- arrange(cam.dist.filter, dist.known)
for (i in 1:nrow(cam.dist.filter)) {
  cam.dist.filter$se.dist.cam.asc[i] <- sd(cam.dist.filter$err.dist.cam[1:i], na.rm = TRUE)/sqrt(length(cam.dist.filter$err.dist.cam[1:i]))
}

# Finally, let's plot cumulative standard error vs distance
# cumulative SE, descending (starting with high)
ggplot(data = cam.dist.filter, aes(x = dist.known, y = se.dist.cam.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Camera distance SE (m)") +
  xlim(0,2000) + ylim(0,15) +
  theme(axis.ticks.length = unit(0.5, "cm"),
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
        axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")))
ggsave("intermediate-products/error/dist-cam-cum-se-desc_2020-stationary.png")


# cumulative SE, ascending (starting with low)
ggplot(data = cam.dist.df, aes(x = dist.known, y = se.dist.cam.asc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Cumulative rangefinder distance SE (asc.)")

# In this particular instance, descending makes far more sense than ascending due to larger errors at greater distances. Errors definitely vary with distance so we will create a linear model to extrapolate errors from this structure
 
lm.dist.cam <- lm(se.dist.cam.desc ~ dist.known, data = cam.dist.filter) # creating linear model
cam.dist.filter$cam.dist.se.pred <- predict(lm.dist.cam, newdata = cam.dist.filter) # predicted values based on known distances from the lm

ggplot(data = cam.dist.filter, aes(x = dist.known, y = se.dist.cam.desc)) +
  geom_line() +
  geom_rug(sides = "b") +
  labs(x = "GPS distance (m)", y = "Camera distance SE (m)") +
  xlim(0,2000) + ylim(0,15) +
  geom_line(aes(y = cam.dist.se.pred), color="red") +
  theme(axis.ticks.length = unit(0.5, "cm"),
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
        axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")))
ggsave("intermediate-products/error/dist-cam-cum-se-desc+model-line_2020-stationary.png")

# finally, we need to create a whole data set for distance measurements
# joining cam.dist.filter (SE info < 2500 m) and the master cam.dist.df
cam.dist.df <- left_join(cam.dist.df, cam.dist.filter)

# since the predecessor of cam.dist.df was dist (which contains all RF error data), this is perfect
write.csv(cam.dist.df, "intermediate-products/error/error_dist_2020_saeborg.csv")

###################
### CONCLUSIONS ###
###################
# Overall, distance estimates produced by both rangefinder and photographic methods produced acceptably small errors (which now need to be propagated through the movement pattern models)
# Rangefinder: very small errors, don't vary much with distance. A single SE value of 0.58 m will be used
# Camera: larger errors but still very small. Very large errors > 2500 m but this is unsurprising given the influence of waves on pixel distances (2-3 pixels can change value sby 100s of metres). Errors vary greatly with distance, so a linear model based on cumulative SE (<2600m, descending) will be used to interpolate values. Values with distances > 2600 m should be excluded.

