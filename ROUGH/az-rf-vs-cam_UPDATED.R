###################################
### UPDATED::: COMPARING RF AND CAMERA AZs ###
###################################
### 21/12/2021
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

# we want good points to derive these errors. Let's remove all Gardar and all Nattfari ULF, UCF, UCC, UCB
loc <- filter(loc, vessel != "gardar", beaufort<3)
# for (i in 1:nrow(loc)) {
#  if (loc$vessel[i] == "nattfari" & c("ulf", "ucf", "ucc", "ucb") %in% loc$plat.b.code[i]) {
#    loc <- loc[-i,]
#  }
#}

# let's also remove distances <100 m
loc <- filter(loc, dist>100)

# now we only want to keep Tom after June for 2018 and Beverly after July for 2020. And we want BFSS <= 2
loc <- rbind(filter(loc, year == 2018 & month(datetime)>5 & observer.b == "tom"),
             filter(loc, year == 2019 & month(datetime)>7 & observer.b == "beverly"),
             filter(loc, year == 2020))
             
# finally, let's filter by points where the average speed (+- 3 seconds) was <5 km/h
# uploading GPS data
gps <- read.csv("intermediate-products/gps_obs-vessel_2018-20.csv") %>% select(-1) %>% mutate(datetime = as.POSIXct(datetime))
gps <- rbind(filter(gps, year(datetime) == 2018 & month(datetime)>5),
             filter(gps, year(datetime) == 2019 & month(datetime)>7),
             filter(gps, year(datetime) == 2020))

# calculate speed for each point
for (i in 1:nrow(loc)) {
  # define the time point
  t = loc$datetime[i] 
  # which row in the GPS thing contains that val?
  row = which(gps$datetime == t)[1]
  # do we have three valid seconds before and after
  if (isTRUE(t-3 == gps$datetime[row-3]) & isTRUE(t+3 == gps$datetime[row+3])) {
    coords = filter(gps, datetime>=t-3 & datetime<=t+3) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% summarise() %>% st_cast("LINESTRING")
    length = as.numeric(st_length(coords))
    loc$speed[i] = (length/1000)/(6/3600)
  }
  print(paste0(i, " out of ", nrow(loc)))
}

# let's say that speed must be <8
loc.filt <- filter(loc, speed<=10)

# now we will save
write.csv(loc.filt, "intermediate-products/position-calc/az-rf_filtered_more.csv", row.names = FALSE)


# now go to Excel, add Az.cam comparisons, delete other rows, save, add to 'data' folder in this project and upload here!
loc.cam <- read.csv("data/position-calc/az-rf_filtered_more.csv") %>%
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

ggplot(data = loc.cam, aes(x = az.cam, y = az.rf, col = as.factor(year))) +
  geom_point() + scale_color_brewer(palette = "Dark2") +
  labs(x = "Camera azimuth", y = "Rangefinder azimuth", col = "Year") +
  geom_abline()
ggsave("intermediate-products/error/az-cam-vs-rf_whale-points.png")
cor(loc.cam$az.cam, loc.cam$az.rf)
cor(loc.cam$az.rf.err, loc.cam$dist)
ggplot(data = loc.cam, aes(x = az.cam, y = az.rf.err, col = as.factor(year))) +
  geom_point()
sd(loc.cam$az.rf.err)
std.error <- function(x) sd(x)/sqrt(length(x))
std.error(loc.cam$az.rf.err)

ggplot(data = loc.cam, aes(x = az.rf.err, col = as.factor(year))) +
  geom_histogram()

# what proportion are within 10 degrees?
nrow(loc.cam)
nrow(filter(loc.cam, abs(az.rf.err)<10))
nrow(filter(loc.cam, abs(az.rf.err)<20))

# see how error varies with distance
ggplot(data = loc.cam, aes(x = dist, y = az.rf.err, col = as.factor(year))) +
  geom_point()

# nearly everything is within 40 degrees. At a few different distances away, what is the difference in position with a 40 degree error?
# we will start at (0,0), in metres. Then calculate position 

ggplot(data = loc.cam, aes(x = dist)) + geom_histogram()
ggplot(data = loc.cam, aes(x = speed, y = az.rf.err, col = as.factor(year))) +
  geom_point()
