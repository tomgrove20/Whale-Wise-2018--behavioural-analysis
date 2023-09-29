

tom <- lines %>%
  dplyr::select(datetime, vessel, boat.ids.30.1500) %>%
  mutate(other.boat.check = ifelse(vessel == boat.ids.30.1500, 0, 1))


t <- follows.tot[3298,] %>% ungroup() %>% dplyr::select(datetime, lon, lat, lon.whale, lat.whale, boat.ids.10.1500) %>% as.data.frame()

h <- hey[3392,] %>% ungroup() %>% dplyr::select(datetime, lon, lat, lon.whale, lat.whale) %>% as.data.frame()

a <- ais[2979,] %>% ungroup() %>% dplyr::select(datetime, lon, lat, lon.whale, lat.whale, boat.ids.10.1500) %>% as.data.frame()


sprintf("%.10f",a$lat.whale)
sprintf("%.10f",t$lat.whale)
sprintf("%.10f",h$lat.whale)

hey <- read.csv("intermediate-products/follows_2018-20_cleaned_calc_full.csv")
