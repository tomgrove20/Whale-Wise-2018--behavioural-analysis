### Filtering AIS data for a specific, short time (could be any time at all)
### 01/07/2020
### Tom Grove

# e.g. can use to validate track times
# e.g. can use to decipher AIS vessel codes

# packages
packages <- c("tidyverse", "sf")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

getwd()

# with the newest ais files
ais <- ais.18
# read.csv("data/skjalfandi_ais_20180505-20180831.csv") %>%
#  mutate(posdate = as.POSIXct(posdate, format = "%Y-%m-%d %H:%M:%S"))

ais.filter <- ais %>%
  filter(posdate > "2018-07-04 11:28:00" & posdate < "2018-07-06 11:40:00")

write.csv(ais.filter, "intermediate-products/ais-check.csv")

# with an older ais file
ais.old <- read.csv("data/skjalfandi_ais_201805-08.csv") %>%
  mutate(posdate = as.POSIXct(posdate, format = "%Y-%m-%d %H:%M:%S"))

ais.filter.old <- ais.old %>%
  filter(posdate > "2018-07-04 11:28:00" & posdate < "2018-07-06 11:40:00")

write.csv(ais.filter.old, "intermediate-products/ais-check_old.csv")



### if we want to filter a processed data set for ww boats only
ais.ww <- read.csv("intermediate-products/ais/ais_2018-20_ww_matched.csv") %>%
  mutate(posdate = as.POSIXct(posdate)) %>%
  filter(posdate > "2018-07-04 08:00:00" & posdate < "2018-07-04 15:00:00")

ais.ww <- ais.ww %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  group_by(vessel) %>% arrange(posdate) %>%
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING")
plot(ais.ww)

ggplot() +
  geom_sf(data = land, color = "black", fill = "black" ) +
  geom_sf(data = ais.ww, aes(color = vessel)) + # filter to may 2018
  coord_sf(crs = 4326, xlim = c(-17.90391, -17.24308), ylim = c(65.98223, 66.18384)) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(axis.text.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(10,0,0,0)),
        text = element_text(size=15),
        panel.grid.major = element_line(colour = "transparent"),
        legend.position = "none")

ggsave("intermediate-products/ais-ww-check.png")

st_write(ais.ww, "intermediate-products/ais-ww-check.shp")


  