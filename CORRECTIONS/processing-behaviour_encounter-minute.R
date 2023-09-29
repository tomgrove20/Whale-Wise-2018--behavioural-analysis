#########################################
### playing with behavioural data processing ###
#########################################
### 23/05/2023
### Tom Grove
### tomgrove20@yahoo.co.uk

# I want to see whether I can construct a new encounter minute variable that can be incorporated into focal data sets
# packages
packages <- c("tidyverse", "data.table", "sf", "nngeo", "lwgeom", "geosphere")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)


follows.encmin <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv") %>%
  mutate(trip.id = paste0(as.Date(datetime),"_",vessel)) %>%
  dplyr::select(trip.id, datetime, folnum.unique, id.focal.catalogue) %>%
  group_by(trip.id, id.focal.catalogue, folnum.unique) %>%
  summarise(start = first(datetime), end = last(datetime)) %>%
  group_by(trip.id, id.focal.catalogue) %>%
  mutate(timediff = ifelse(row_number()>1, difftime(start,lag(end), units = "mins"), 0)) %>%
  # now defining an enc.num variable
  mutate(enc.num = ifelse(is.na(timediff), 0, as.integer(factor(cumsum(timediff >= 20))))) %>%
  group_by(trip.id, id.focal.catalogue, enc.num) %>%
  mutate(enc.num = cur_group_id())
write.csv(follows.encmin, "intermediate-products/follows/follows_start+end_consecutive_rough.csv", row.names = FALSE)


# what I need to do next: for the swim speed foc model, left_join follows.encmin to the tot data frame. You can then use the enc.num grouping and start time of each follow (from follows.encmin) to recalculate encounter minute. Then sqrt it and add to the model to check results