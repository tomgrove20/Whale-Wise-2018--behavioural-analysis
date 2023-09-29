### Joining PocketGIS tracks
### 27/06/2020
### Tom Grove, tomgrove20@yahoo.co.uk

# We used PocketGIS to track location every second while conducting behavioural observation on whale-watching boats. We need to check the overlap between PocketGIS tracks and AIS locations. If AIS locations do not cover full range of Pocket, we will have to ask Marianne Rasmussen (mhr@hi.is) to ask the Icelandic coastguard for more AIS data.
# Covering 2018 and 2019
# We need to remove all rangefinder points from this data set, leaving only vessel (observer) tracks

# packages
packages <- c("tidyverse", "purrr")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

getwd() # make sure it corresponds to Rproj location

## MERGING

# upload pocketGIS data
files <- list.files(path = "data/pocketgis/whale-watching", pattern ="*.csv", recursive = TRUE)
files # worked!

# now let's join into one
points.list <- lapply(paste0("data/pocketgis/whale-watching/",files), read.csv)

points <- points.list %>% reduce(full_join)

## FILTERING

# look at unique values for each thing
unique(points$type) # let's remove anything that says 'whale' (rangefinder points)
unique(points$Height) # let's remove all zeros and anything with more than 1 dp (rangefinder points)
unique(points$WHALE.BOAT)
unique(points$whale.boat.) # let's remove anything that says whale (rangefinder points)
unique(points$X9)
unique(points$X) # let#s remove anything that says whale
unique(points$X.1)
unique(points$X.2)
unique(points$X.3)

# When filtering, we will create a new column called 'height.10', = height * 10. We will then filter out any records where height.10 is not a whole number
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

points.filter <- points %>%
  mutate(height.10 = Height * 10) %>% # creating column = height * 10. If this isn't an integer, it shouldn't be present in the data frame
  filter(type != "whale" | is.na(type)) %>% # removing anything that is whale
  filter(whale.boat. != "whale" | is.na(whale.boat.)) %>% # removing anything that is whale
  filter(Height != 0) %>% # remove height = 0
  filter(is.wholenumber(height.10)) %>% # filter out any measured point
  filter(X != "whale" | is.na(X)) %>% # removing anything that is whale
  rename(east = Easting, north = Northing, date = ChangeDate, time = ChangeTime, height =  Height) %>%
  select(date, time, east, north, height)

head(points.filter) # looking good!

# save
write.csv(points.filter, "intermediate-products/pocketgis_vessel_2018-19.csv")
