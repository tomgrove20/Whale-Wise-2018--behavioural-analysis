####################################################
### SHOULD WE PREFERENTIALLY USE CAMERA ANGLES?? ###
####################################################
### 26/01/2020
### Tom Grove
### tomgrove20@yahoo.co.uk

# How much can/should we rely on rangefinder azimuths when they are so poor?
# STEP 1: look at 100 random RF bearings paired with Cam bearings from each of 2018,19,20. Plot RF bearing error vs bearing

# packages
packages <- c("tidyverse", "ggplot2", "sf", "geosphere", "lwgeom", "exifr", "leaflet", "viridis", "data.table", "lubridate")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

######################
### DATA FILTERING ###
######################

# filtered movement data (from contiguous movement follows) from all three years 
loc <- st_read("final-products/follows_2018-20_Mn_positions.gpkg") %>%
  select(datetime, observer.a, observer.b, beaufort, declination, lon, lat, vessel, plat.a.code, plat.b.code, direction.time.a, direction.time.b,
         photo.name, img.y, img.x, fl, fov, ref.id, whale.pixel.x, ref.pixel.x,
         dist.rf, dist.cam, dist,
         az.rf, az.cam, az, move.id.tot) %>%
  mutate(sur = row_number(), # assign surfacing number before filtering (comparison with other R codes)
         year = year(datetime)) %>% 
  filter(!is.na(az.rf),
         !is.na(img.y)) 

# now we will split up into each year, randomise the row order and save
for (i in unique(loc$year)) {
  
  loc.filt = filter(loc, year == i) %>%
    slice(sample(1:n()))
  
  write.csv(loc.filt, paste0("intermediate-products/position-calc/az-rf_filtered_", i,".csv"))
}

# now go to Excel, we need 100 Az.cam comparisons, delete other rows, save, add to 'data' folder in this project and upload here!
loc.cam <- rbind(read.csv("data/position-calc/az-rf_filtered_2018.csv"),
                 read.csv("data/position-calc/az-rf_filtered_2019.csv"),
                 read.csv("data/position-calc/az-rf_filtered_2020.csv")) %>%
  select(-c(X, geom, X.1)) %>%
  mutate(az.rf = az.rf %% 360)
head(loc.cam)

# now upload camera refs
cam.ref <- read_sf("data/camera-angle-refs/camera-angle-refs.shp") %>% # currently in ISN93
  st_transform(crs = "EPSG:4326") %>% # transform to WGS84
  rename(geom.ref = geometry, ref.id = name) %>%
  mutate(lon.ref = unlist(map(geom.ref,1)),
         lat.ref = unlist(map(geom.ref,2))) %>%
  as.data.frame %>%
  select(-c("id", "XCOORD", "YCOORD", "geom.ref"))

# let's check that all camera angle references within the follows DF match to one in the cam ref shapefile
cam.check <- follows %>%
  select(ref.id) %>%
  distinct %>%
  left_join(cam.ref)
view(cam.check)
unique(cam.ref$ref.id) 

# if everything looks good, then we can join!
loc.cam <- loc.cam %>%
  left_join(cam.ref)

# calculate camera azimuth
loc.cam <- loc.cam %>% rowwise() %>% # ensures that each row is treated separately
  mutate(az.ref = (bearing(c(lon, lat), c(lon.ref, lat.ref)) + 360) %% 360,
         az.cam = (az.ref + ((whale.pixel.x - ref.pixel.x)/img.x)*fov) %% 360, # make sure always between 0 and 360
         az.rf.err = (az.rf-az.cam + 180) %% 360-180)
sum(!is.na(loc.cam$az.cam))

ggplot(data = loc.cam, aes(x = az.cam, y = az.rf.err, col = as.factor(year))) +
  geom_point()
sd(loc.cam$az.rf.err)

ggplot(data = loc.cam, aes(x = az.rf.err, col = as.factor(year))) +
  geom_histogram()

# what proportion are within 10 degrees?
nrow(loc.cam)
nrow(filter(loc.cam, abs(az.rf.err)<10))
nrow(filter(loc.cam, abs(az.rf.err)<20))

# see how error varies with distance
ggplot(data = loc.cam, aes(x = dist, y = az.rf.err, col = as.factor(year))) +
  geom_point()

# see how error varies with vessel
loc.cam %>% group_by(vessel) %>%
  summarise(mean = mean(az.rf.err), sd = sd(az.rf.err))


# nearly everything is within 40 degrees. At a few different distances away, what is the difference in position with a 40 degree error?
# we will start at (0,0), in metres. Then calculate position 

