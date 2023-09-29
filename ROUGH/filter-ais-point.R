#################################################
### CALCULATING BOAT NUMBER/DISTANCE FROM AIS ###
#################################################
### 09/01/2020
### Tom Grove
### tomgrove20@yahoo.co.uk

# Here, we aim to calculate the number of boats within a certain distance of a given point. Within this range, we will also calculate the distance of each vessel to the point. This can be used for behavioural and acoustic analyses

### PACKAGES
packages <- c("tidyverse", "ggplot2", "sf", "data.table")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

### PLOT THEME
theme_set(theme_bw()) 
map.theme <- theme(axis.text.y = element_text(margin = margin(0,10,0,0)),
                   axis.text.x = element_text(margin = margin(10,0,0,0)),
                   text = element_text(size=15))
 
# create cropped Skjalfandi shape file
skjalfandi <- st_read("data/iceland-isn93/is50v_strandlina_flakar_24122017.shp") %>%
  st_transform(crs = 4326) %>%
  st_crop(ymin = 65.94, ymax = 66.24, xmin = -18.04, xmax = -17.19)

### DATA 
# pre-matched 2020 AIS data. Easier to convert to a planar CRS (3057)
ais.20 <- read.csv("intermediate-products/ais/ais_20200611-20200907_ww_matched.csv") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% st_transform(crs = 3057) %>%
  mutate(posdate = as.POSIXct(posdate)) # create datetime format, easier for analyses
  
### CHECKING 
# quick filter of AIS data and plot
ais.20 %>%
  filter(posdate > as.POSIXct("2020-08-15 11:20:00") & posdate < as.POSIXct("2020-08-15 11:40:00")) %>%
  ggplot() +
  geom_sf(aes(color = vessel)) + geom_sf(data = skjalfandi) + 
  coord_sf(crs = 4326) +
  map.theme + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))


################
### ONE TIME ###
################

# defining reference point in 4326, then convert to 3057 for planar calculation
ref.x = -17.5
ref.y = 66.15
ref <- st_sfc(st_point(c(ref.x, ref.y)), crs = 4326) %>% st_transform(crs = 3057)

# defining reference distance and time
ref.time <- as.POSIXct("2020-08-15 11:30") # we will search 10 minutes either side of this time
ref.dist = 2000 # 2km for now

ggplot() +
  geom_sf(data = ref) +
  geom_sf(data = st_buffer(ref, ref.dist), fill = "transparent", size = 1) +
  geom_sf(data = skjalfandi) + coord_sf(crs = 4326) + map.theme +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

# now we can filter! 
filter.time <- ais.20 %>%
  # search for the 10 minutes either side!
  filter(posdate > ref.time - 10*60 & posdate < ref.time + 10*60)

filter <- st_sf(st_intersection(st_buffer(ref, ref.dist), filter.time)) %>% 
  rename(geometry = 1) %>% st_join(filter.time) %>%
  # calculate distance for each point
  mutate(dist = as.numeric(st_distance(geometry, ref, by_element = TRUE)))
head(filter)

ggplot() +
  geom_sf(data = st_buffer(ref, ref.dist), fill = "transparent", size = 1) +
  geom_sf(data = filter, aes(col = vessel), size = 2, alpha = 0.5) +
  geom_sf(data = skjalfandi) +
  coord_sf(crs = 4326) + map.theme + 
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  labs(col = "Vessel")

write.csv(filter, "ais-point-test_1-time.csv")

# Just for fun, let's look at the minimum distance and calculable min/mean/max speed of each vessel

# minimum distance
filter %>% group_by(vessel) %>% summarise(min = min(dist)) 

# min, mean, max speed
filter %>% arrange(posdate) %>% group_by(vessel) %>%
  mutate(
    # distance between vessel at time i and time i+1
    length = as.numeric(st_distance(geometry, geometry[row_number() + 1], by_element = T)),
    # time difference between vessel at time i and time i + 1
    time.diff = as.numeric(difftime(posdate[row_number() + 1], posdate, units = "secs")),
    # speed = distance/time
    speed = (length/time.diff) * (60*60/1000)) %>%
  # min, mean and max speed for each vessel
  summarise(mean = mean(speed, na.rm = TRUE),
            min = min(speed, na.rm = TRUE),
            max = max(speed, na.rm = TRUE))


##################
### MANY TIMES ###
##################

# times could be a written vector or column from a data frame
times <- c(as.POSIXct("2020-08-15 11:40:00"), as.POSIXct("2020-08-16 11:40:00"), as.POSIXct("2020-08-17 11:40:00"))

# now we can filter! 
filter.time <- ais.20 %>%
  # search for the 10 minutes either side!
  filter(inrange(posdate, times - 10*60, times + 10*60))

# plot how filtered times look
ggplot() +
  geom_sf(data = filter.time, aes(col = as.factor(as.Date(posdate)))) +
  geom_sf(data = skjalfandi) +
  geom_sf(data = st_buffer(ref, ref.dist), col = "black", fill = "transparent", size = 1) + 
  scale_color_brewer(palette = "Dark2") +
  coord_sf(crs = 4326) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  labs(col = "Date") + map.theme

# now filter by space
filter <- st_sf(st_intersection(st_buffer(ref, ref.dist), filter.time)) %>% rename(geometry = 1) %>%
  st_join(filter.time) %>%
  # calculate distance for each point
  mutate(dist = as.numeric(st_distance(geometry, ref, by_element = TRUE)))

# and view full filter!
ggplot() +
  geom_sf(data = st_buffer(ref, ref.dist), col = "black", fill = "transparent", size = 1) +
  geom_sf(data = filter, aes(col = as.factor(as.Date(posdate))), size = 1.5) +
  geom_sf(data = skjalfandi) +
  coord_sf(crs = 4326) + map.theme + 
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  labs(col = "Date")

write.csv(filter, "ais-point-test_many-times.csv")



