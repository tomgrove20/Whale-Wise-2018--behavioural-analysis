####################################################
### Number of whale watching trips from AIS data ###
####################################################

# We have drawn a line from across the entrance the Husavik Harbour. The idea is that, each trip, the boat should cross this line twice. Thus, by counting the number of crosses and dividing by 2 for each whale watching boat, we can calculate the daily number of trips!

# packages
packages <- c("tidyverse", "ggplot2", "viridis", "sf", "fuzzyjoin", "chron", "mapview")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

# upload AIS data (WW boats only)
ais <- read.csv("intermediate-products/ais/ais_20200611-20200907_ww_matched.csv") %>%
  mutate(posdate = as.POSIXct(posdate))

# upload harbour entrance line shapefile
harbour <- read_sf("data/harbour-entrance-line.gpkg")

# let's first erase anything between times 0000 and 0750 (will make sense later, also saves processing time)
ais <- ais %>%
  mutate(time = times(format(posdate, "%H:%M:%S")),
         date = as.Date(posdate)) %>%
  filter(time > "07:30:00" & time < "23:00:59")
head(ais)

# if you only want specific dates, use this code!
# date.start <- as.Date("2020-01-01") # change!
# date.end <- as.Date("2020-12-01") # change!
# ais <- filter(ais, date >= as.Date(date.start) & date <= as.Date(date.end))

# now, we can calculate the number of intersections per boat with the harbour entrance line. Generally, number of trips = 0.5*number of intersections. For some midnight trips, only one intersection (0000-0750 removed), so we round up!

ais.trips <- ais %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% # convert to sf for intersection calculations
  group_by(date, vessel) %>%
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING") %>%
  mutate(
    # count number of points using mapview package, round up to nearest number (in case of late trips that return after midnight)
    n.intersect = npts(st_intersection(harbour)),
    # now convert to number of trips by dividing by 2 and rounding up to nearest integer (in case of late trips that return after midnight)
    n.trip = round(n.intersect/2, digits = 0)) %>%
    select(-geometry) %>% as.data.frame()

head(ais.trips)

# finally, let's summarise this to number of trips per day!
ais.trips.short <- ais.trips %>%
  group_by(date) %>%
  summarise(n.trip = sum(n.trip)) %>%
  complete(date = seq.Date(min(date), max(date), by="day")) %>% # fill in missing dates 
  mutate(n.trip = ifelse(is.na(n.trip),0,n.trip)) # for those missing dates, number of trips = 0
head(ais.trips.short)
summary(ais.trips.short) # for min, max, mean




