#################################################################
### Whale measurement error: vessel (both vessels stationary) ###
#################################################################
### 09/09/2020
### Tom Grove
### tomgrove20@yahoo.co.uk

# Here, we will look at errors for all angle/distance measurement methods from a whale-watching vessel. Both vessels were stationary (or approximately) when measurements were taking place
# On August 6th and 18th 2020, measurements were taken from Saeborg by Amelie and Beni. The test 'whale' was the research inflatable, with Alli as the captain and Aleksandra Lechawar as the GPS. On August 18th, Heimir Hardarson drove Saeborg out of the harbour
# Our reference values will be calculated from known GPS positions of the observer and inflatable (test whale)
# NOTE: remember to add compass corrections
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

file.1 <- "data/error/20200806_saeborg_error_follows.csv"
file.2 <- "data/error/20200818_saeborg_error_follows.csv"

# error follow (not adding time for rangefinder firing delay for now)
follow <- full_join(read.csv(file.1), read.csv(file.2)) %>%
  mutate(datetime = as.POSIXct(paste0(date," ",timeutc), format = "%Y%m%d %H:%M:%S")) %>%
  mutate(photo = paste0(date,"_",tripstart,"_",vessel,"_error_",distance.photo,".JPG")) # photo file name

# if you wanted to subset. follow <- follow[96:114,]

# we need to correct rangefinder bearings by magnetic declination
# Location: (66.041448, -17.353420), near Husavik Harbour
# Date: 2020-08-06
# Magnetic declination = -11.1. Therefore, we subtract 11.1 on to every measurement
follow <- follow %>%
  mutate(rangefinder.bearing = rangefinder.bearing - 11.1) # to account for magnetic declination

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

# photo data

files <- list.files(path = "data/error/2020_error_photos/", pattern = "*(?i).JPG", recursive = TRUE) # (?i) makes the search case-insensitive

photo.info <- read_exif(paste0("data/error/2020_error_photos/",files), tags=c("FileName", "DateTimeOriginal", "ExifImageWidth", "ExifImageHeight", "FocalLength", "FOV")) %>%
  rename(photo = FileName, date.time = DateTimeOriginal, img.x = ExifImageWidth, img.y = ExifImageHeight, fl = FocalLength, fov = FOV  ) %>% # renaming to make easier
  mutate(date.time = as.POSIXct(date.time, format = "%Y:%m:%d %H:%M:%S")) # convert date.time to proper format

# Might as well save these data for later (prevent need for photos to take up space)
write.csv(photo.info, "data/error/photo-data_error_2020_saeborg.csv")

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
cam.check
unique(cam.ref$ref.id) 

follow <- follow %>%
  left_join(cam.ref)


# vessel track

# for the 20200806 trip, we will just use a single location (since Saeborg was stationary the entire time)

# first session
vess1 <- as.data.frame(seq.POSIXt(min(filter(follow, date == 20200806)$datetime), max(filter(follow, date == 20200806)$datetime), by = "sec")) %>%
  rename(datetime = 1) %>% # creating list of every second from start to finish
  mutate(lat = 66.045785,
         lon = -17.345754) %>%
  mutate(lat = ifelse(datetime > "2020-08-06 21:13:30", 66.045735, lat),
         lon = ifelse(datetime > "2020-08-06 21:13:30", -17.345708, lon)) # a late addition. seemed their position changed in the second half (when they started range finder distances approximately)
head(vess1)

# second session
vess2 <- read.csv("data/error/20200818_1615_error_saeborg_gps.csv", fileEncoding="UTF-16LE", sep = "\t", header = TRUE) %>%
  rename(datetime = time, lat = Latitude, lon = Longitude) %>%
  select(datetime, lat, lon) %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"))
head(vess2)

# every four seconds, let's fill in missing seconds
secs <- as.data.frame(seq.POSIXt(min(vess2$datetime), max(vess2$datetime), by = "sec")) %>%
  rename(datetime = 1) # creating list of every second from start to finish

vess2 <- full_join(vess2, secs) # adding rows for missing seconds

# creating a function to fill missing values from nearest value 
fill.near <- function(dat) {
  N <- length(dat)
  na.pos <- which(is.na(dat))
  if (length(na.pos) %in% c(0, N)) {
    return(dat)
  }
  non.na.pos <- which(!is.na(dat))
  intervals  <- findInterval(na.pos, non.na.pos,
                             all.inside = TRUE)
  left.pos   <- non.na.pos[pmax(1, intervals)]
  right.pos  <- non.na.pos[pmin(N, intervals+1)]
  left.dist  <- na.pos - left.pos
  right.dist <- right.pos - na.pos
  
  dat[na.pos] <- ifelse(left.dist <= right.dist,
                        dat[left.pos], dat[right.pos])
  return(dat)
}

vess2 <- vess2 %>%
  arrange(datetime) %>%
  mutate(lat = fill.near(lat), lon = fill.near(lon))
head(vess2) # perfect!

# now we can join vessel tracks from the two sessions and convert to shapefile
vess <- full_join(vess1, vess2) %>%
  st_as_sf(coords = c("lon", "lat"), crs = "EPSG:4326", remove = FALSE) %>%
  rename(geom.vess = geometry, lon.obs = lon, lat.obs = lat)

# test whale track (inflatable)

file1 <- "data/error/20200806_2030_error_inflatable_gps.csv"
file2 <- "data/error/20200818_1615_error_inflatable_gps.csv"

whale <- full_join(read.table(file1, fileEncoding="utf-16", sep = "\t", header = TRUE),
                   read.table(file2, fileEncoding="utf-16", sep = "\t", header = TRUE))%>%
  rename(datetime = time, lat = Latitude, lon = Longitude) %>%
  select(datetime, lat, lon) %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = "EPSG:4326", remove = FALSE) %>%
  rename(lon.whale = lon, lat.whale = lat, geom.whale.known = geometry) %>%
  as.data.frame # convert back to data frame
head(whale)

# Before getting into this, I recommend looking at all tracks in QGIS. Just to see they make sense
# I have done this and it seems to make sense

## ACTUAL DISTANCE/ANGLE
# First, we will calculate actual distance and bearing  
follow.tot <- follow %>%
  left_join(vess, by = "datetime") %>%
  left_join(whale, by = "datetime") %>%
  mutate(dist.known = st_distance(geom.vess, geom.whale.known, by_element = TRUE)) %>%
  mutate(dist.known = as.numeric(str_remove(dist.known, "[m]")))

for (i in 1:nrow(follow.tot)) {
  
  follow.tot$az.known[i] = (bearing(c(follow.tot$lon.obs[i], follow.tot$lat.obs[i]), c(follow.tot$lon.whale[i], follow.tot$lat.whale[i])) + 360) %% 360

}

head(follow.tot)

# saving as tab-separated txt file, not comma-separated csv, because geom columns would be a problem otherwise (commas would automatically separate into separate columns)
write.table(follow.tot, "intermediate-products/error/error_raw-tot_202008_saeborg.txt", sep="\t")

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
  theme_bw() # overriding theme_classic
ggsave("intermediate-products/error/skjalfandi-camera-refs.png")
