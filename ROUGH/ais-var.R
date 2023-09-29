###################################
### Vessel information from AIS ###
###################################
### 23/11/2021
### Tom Grove (tomgrove20@yahoo.co.uk), Whale Wise

# Aim: for each follow position, calculate vessel-related variables from AIS data. Calculate the number of vessels, min/max/mean/SD; average oak boat speed; average RIB speed; presence and absence of RIB and oak boats; within different distances. 

# packages
packages <- c("tidyverse", "sf", "lubridate")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

############
### DATA ###
############
# AIS data, total
ais <- read.csv("intermediate-products/ais/ais_20180505-20180831_ww_matched.csv") %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  mutate(posdate = as.POSIXct(posdate),
         daynum = as.numeric(as.Date(posdate))) %>% # create day number field for quicker filtering
  group_by(vessel) %>%
  mutate(n = row_number()) %>% # to later detect whether filtered positions represent consecutive reported points
  ungroup()

# we need to assign each vessel to oak or rib
vess.type <- data.frame(vessel = c("salka", "nattfari", "saeborg", "andvari", "fanney", "sylvia", "bjossi", "gardar", "a-helga", "a-kibba",  "kjoi", "kria", "faldur", "a-sigga", "haukur", "a-johanna", "opal", "hildur", "a-sigga-2"),
                        type = c("oak", "oak", "oak", "oak", "oak", "oak", "oak", "oak", "rib", "rib", "rib", "rib", "oak", "rib", "oak", "rib", "oak", "oak", "rib"))

ais <- left_join(ais, vess.type)

# Humpback data, positions calculated
mn <- st_read("final-products/follows_2018-20_Mn_positions.gpkg") %>%
  mutate(daynum = as.numeric(as.Date(as.POSIXct(datetime)))) %>% # create day number field for quicker filtering
  filter(datetime > as.POSIXct("2018-07-10 00:00:00") & datetime < as.POSIXct("2018-07-15 23:59:59")) # just to keep quick for now!!

##########################
### FILTER/CALCULATION ###
##########################
# defining the lag time you want, to create a time window for variable calculations 
lag <- 20 # minutes

# we also want short and long distance thresholds, to allow calculation of variables at two different spatial scales
d.short <- 1500 # metres
d.long <- 5000 # metres
d.list <- c(d.short, d.long)

mn.calc <- mn
for (i in 1:nrow(mn.calc)) {
  # filter by times
  filter.time <- ais %>% filter(daynum == mn.calc$daynum[i]) %>%
    filter(posdate > mn.calc$datetime[i]-(lag*60) & posdate < mn.calc$datetime[i])
  
  # filter by location
  filter <- st_sf(st_intersection(st_buffer(select(mn.calc[i,],geom), d.short), filter.time))
  
  # presence of oak and rib in last 30 minutes
  mn.calc$oak.pres[i] = "oak" %in% filter$type; rib.pres = "rib" %in% filter$type
  
  # number of oaks and ribs
  mn.calc$num.oak[i] = length(unique(filter(filter, type == "oak")$vessel))
  mn.calc$num.rib[i] = length(unique(filter(filter, type == "rib")$vessel))
  
  # max/average/SD time, dist, speed of oaks and RIBs. This is a bit complex. 
  speeds <- filter %>% group_by(vessel) %>%
    mutate(posdate.1 = posdate[row_number()+1],
           geom.1 = geom[row_number()+1],
           valid = n[row_number()+1]) %>%
    filter(!is.na(posdate.1)) %>% rowwise() %>%
    mutate(time = as.numeric(posdate.1 - posdate)*60,
           dist = as.numeric(st_distance(geom, geom.1)), # distance in metres (remove units with as.numeric)
           speed = dist/time) 
  
  mn.calc$time.tot[i] = sum(speeds$time) # total time vessels occupy the area, seconds
  mn.calc$time.oak[i] = sum(filter(speeds, type == "oak")$time)
  mn.calc$time.rib[i] = sum(filter(speeds, type == "rib")$time)
  
  mn.calc$dist.tot[i] = sum(speeds$dist) # total distance covered by all boats, metres
  mn.calc$dist.rib[i] = sum(filter(speeds, type == "rib")$dist)
 
  mn.calc$speed.oak[i] = mn.calc$dist.oak[i]/mn.calc$time.oak[i] # total and max speeds for oak and rib
  mn.calc$speed.rib[i] = mn.calc$dist.rib[i]/mn.calc$time.rib[i]
  mn.calc$max.speed.oak[i] = max(filter(speeds, type == "oak")$speed, na.rm = TRUE)
  mn.calc$max.speed.rib[i] = max(filter(speeds, type == "rib")$speed, na.rm = TRUE)
  
  print(paste0(100*i/nrow(mn.calc)," %")) # progress checker
}

# a couple of plots for fun!
ggplot(data = mn.calc) +
  geom_point(aes(x = datetime, y = dist.oak))
ggplot(data = mn.calc) +
  geom_point(aes(x = datetime, y = speed.oak))

# how many follows had only one boat for the entire follow? AN IMPORTANT FILTER!
test <- mn.calc %>% rowwise() %>% mutate(num.tot = sum(num.oak + num.rib)) %>% group_by(folnum.unique) %>% summarise(num.max = max(num.tot))
ggplot(data = test) +
    geom_bar(stat = "identity", aes(x = folnum.unique, y = num.max))



