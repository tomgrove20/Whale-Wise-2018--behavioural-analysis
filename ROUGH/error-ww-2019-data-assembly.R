### Whale measurement error: vessel
### 29/06/2020
### Tom Grove (tom.grove@ed.ac.uk)

# Here, we will look at errors for all angle/distance measurement methods from a whale-watching vessel
# Our reference values will be calculated from known GPS positions of the observer and inflatable (test whale)
# NOTE: remember to add compass corrections
# NOTE: remember to add Constance-Hannah corrections
# NOTE: start by calculating everything in EPSG 4326 (WGS84, lat/lon). This is important to calculate distances from photos (relies on earth's curvature)

# packages
packages <- c("tidyverse", "ggplot2", "sf", "geosphere", "lwgeom", "exifr", "leaflet", "viridis")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

# ggplot theme
theme_set(theme_classic()) 
theme_update(axis.title.y = element_text(margin = margin(0,10,0,0)),
             axis.title.x = element_text(margin = margin(10,0,0,0)))

## LOAD DATA

# error follow
follow <- read.csv("data/error/20190820_follows_error_vessel.csv") %>%
  mutate(datetime = as.POSIXct(paste0(date," ",timeutc), format = "%Y%m%d %H:%M:%S")) %>%
  mutate(time.rf = datetime - 2) %>% # takes about 3 seconds to fire the range finder
  mutate(photo = paste0(date,"_",tripstart,"_",vessel,"_error_",distance.photo,".JPG")) # photo file name

# if you wanted to subset. follow <- follow[96:114,]

# we need to correct rangefinder bearings by magnetic declination
# Location: (66.085024, -17.536986), centre of Skjalfandi Bay
# Date: 2019-08-20
# Magnetic declination = -11.55. Therefore, we add 11.55 on to every measurement
follow <- follow %>%
  mutate(rangefinder.bearing = rangefinder.bearing - 11.55) # to account for magnetic declination

# Vessel platform height data. This is crucial for photo distance calculation. Also important when looking for sources of error. Create separate tables for a and b (makes joining easier)
plat.height.a <- read.csv("data/platform-height_images.csv") %>%
  select(vessel, platform, height.plat, stand.sit) %>%
  rename(platform.a = platform, height.plat.a = height.plat, stand.sit.a = stand.sit)

plat.height.b <- read.csv("data/platform-height_images.csv") %>%
  select(vessel, platform, height.plat, stand.sit) %>%
  rename(platform.b = platform, height.plat.b = height.plat, stand.sit.b = stand.sit)

# observer height data. Needed to calculate precise height of camera. Again, create separate dfs for a and b
obs.a.height <- read.csv("data/observer-height.csv") %>%
  rename(observer.a = name, stand.a = stand, sit.a = sit)

obs.b.height <- read.csv("data/observer-height.csv") %>%
  rename(observer.b = name, stand.b = stand, sit.b = sit)

# before adding in vessel and observer height data, we need to ensure the platform codes are standardised
unique(filter(plat.height.a, vessel == "saeborg")$platform)
unique(follow$platform.b) # matches
unique(follow$platform.a) # matches
# all looks good, now we can join data

# adding height data
follow <- follow %>%
  left_join(plat.height.a, by = c("vessel", "platform.a")) %>% # adding vessel platform a heights
  left_join(plat.height.b, by = c("vessel", "platform.b")) %>% # adding vessel platform b heights
  left_join(obs.a.height, by = "observer.a") %>% # adding observer a heights
  left_join(obs.b.height, by = "observer.b") # adding observer b heights

head(follow)

files <- list.files(path = "data/error/20190820_1230_saeborg_error/", pattern = "*.JPG")

# photo data
photo.info <- read_exif(paste0("data/error/20190820_1230_saeborg_error/",files), tags=c("FileName", "DateTimeOriginal", "ExifImageWidth", "ExifImageHeight", "FocalLength", "FOV")) %>%
  rename(photo = FileName, date.time = DateTimeOriginal, img.x = ExifImageWidth, img.y = ExifImageHeight, fl = FocalLength, fov = FOV  ) %>% # renaming to make easier
  mutate(date.time = as.POSIXct(date.time, format = "%Y:%m:%d %H:%M:%S")) # convert date.time to proper format

# Might as well save these data for later (prevent need for photos to take up space)
write.csv(photo.info, "data/error/photo-data_error_20190820_saeborg_1230.csv")

# adding photo data
follow <- follow %>%
  left_join(photo.info, by = "photo")

# camera angle refs
cam.ref <- read_sf("data/camera-angle-refs/camera-angle-refs.shp") %>% # currently in ISN93
  st_transform(crs = "EPSG:4326") %>% # transform to WGS84
  rename(geom.ref = geometry, ref.id = name) %>%
  mutate(lon.ref = unlist(map(geom.ref,1)),
         lat.ref = unlist(map(geom.ref,2))) %>%
  select(-c("id", "XCOORD", "YCOORD"))
head(cam.ref)

cam.check <- follow %>%
  select(ref.id) %>%
  distinct(ref.id) %>%
  left_join(cam.ref)
cam.check # farm.past.tungen is empty
unique(cam.ref$ref.id) # tungen.past.farm

follow <- follow %>%
  mutate(ref.id = replace(ref.id, ref.id == "farm.past.tungen", "tungen.past.farm")) %>%
  left_join(cam.ref)

# vessel track
vess <- read.csv("data/error/20190820_1230_saeborg_gps.csv") %>%
  rename(east = Easting, north = Northing, date = ChangeDate, time = ChangeTime) %>%
  select(east, north, date, time) %>%
  mutate(datetime = as.POSIXct(paste0(date," ",time), format = "%Y-%m-%d %H:%M:%S")) %>%
  mutate(datetime = datetime - 3) %>% # we think vessel track is about 3 seconds ahead of inflatable (whale) track most of the way through, so taking 3 seconds off for the first part (deduced by comparing with camera angle)
  st_as_sf(coords = c("east", "north"), crs = "EPSG:3057") %>%
  st_transform(crs = "EPSG:4326") %>%
  rename(geom.vess = geometry) %>%
  mutate(lon.obs = unlist(map(geom.vess,1)),
         lat.obs = unlist(map(geom.vess,2))) %>%
  select(-c("date","time"))
head(vess)

# we also need to create a separate data frame for extracting locations of the vessel 3 seconds later, when the rangefinder delay might kick in
loc.rf <- vess %>% 
  rename(lon.rf = lon.obs, lat.rf = lat.obs, time.rf = datetime, geom.rf = geom.vess)

# test whale track (inflatable)
whale <- read.csv("data/error/20190820_1230_inflatable_loggerpro.csv", fileEncoding="utf-16", sep = "\t") %>%
  rename(datetime = time, lat = Latitude, lon = Longitude) %>%
  select(datetime, lat, lon) %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = "EPSG:4326", remove = FALSE) %>%
  rename(lon.whale = lon, lat.whale = lat, geom.whale.known = geometry) %>%
  as.data.frame # convert back to data frame
head(whale)

# Before getting into this, I recommend looking at both tracks in QGIS. Just to see they make sense
# I have done this and it seems to make sense

## ACTUAL DISTANCE/ANGLE
# First, we will calculate actual distance and bearing  
follow.tot <- follow %>%
  left_join(vess, by = "datetime") %>%
  left_join(whale, by = "datetime") %>%
  left_join(loc.rf, by = "time.rf") %>%
  mutate(dist.known = st_distance(geom.vess, geom.whale.known, by_element = TRUE)) %>%
  mutate(dist.known = as.numeric(str_remove(dist.known, "[m]"))) %>%
  mutate(dist.delay.rf = st_distance(geom.rf, geom.whale.known, by_element = TRUE)) %>%
  mutate(dist.delay.rf = as.numeric(str_remove(dist.delay.rf, "[m]")))

for (i in 1:nrow(follow.tot)) {
  # known azimuth when record started
  follow.tot$az.known[i] = (bearing(c(follow.tot$lon.obs[i], follow.tot$lat.obs[i]), c(follow.tot$lon.whale[i], follow.tot$lat.whale[i])) + 360) %% 360
  
  # known azimuth when vessel is 3 seconds later (accounting for range finder delay)
  follow.tot$az.delay.rf[i] = (bearing(c(follow.tot$lon.rf[i], follow.tot$lat.rf[i]), c(follow.tot$lon.whale[i], follow.tot$lat.whale[i])) + 360) %% 360
}

# we're then going to correct for points where we focused on Constance instead of Hanna (who was carrying the tablet). Therefore, we'll change az.known to account for this

follow.tot <- follow.tot %>%
  mutate(az.alter = ifelse(!is.na(hanna.px),((constance.px - hanna.px)/img.x)*fov, NA),
         rangefinder.bearing = ifelse(!is.na(hanna.px), rangefinder.bearing - az.alter, rangefinder.bearing))

head(follow.tot)

# saving as tab-separated txt file, not comma-separated csv, because geom columns would be a problem otherwise (commas would automatically separate into separate columns)
write.table(follow.tot, "intermediate-products/error/error_raw-tot_20190820_saeborg.txt", sep="\t")

# OK, so we have derived known distance and angles. Now let's extract each of the measurements

### PLOTS TO SAVE
# writing these at the end, since they're not part of the actual data prep or checking

# create sf for Iceland land
land.wgs <- st_transform(read_sf("data/iceland-isn93/is50v_strandlina_flakar_24122017.shp"), crs = 4326)

## Camera refs 
ggplot() +
  geom_sf(data = land.wgs, color = "lightgray") +
  geom_sf(data = cam.ref, color = "purple") +
  coord_sf(xlim = c(min(st_coordinates(cam.ref)[,1]), max(st_coordinates(cam.ref)[,1])), ylim = c(min(st_coordinates(cam.ref)[,2]), max(st_coordinates(cam.ref)[,2]))) +
  theme_bw() + # overriding theme_classic 
  theme(text = element_text(size=10))
ggsave("intermediate-products/error/skjalfandi-camera-refs.png")
