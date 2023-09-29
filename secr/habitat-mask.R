### Creating a habitat mask in secr
### Tom Grove
### 13.02.2020

library(dplyr)
library(secr)
library(ggplot2)
library(sf)
library(raster)
library(sp)

# More information on habitat masks: https://www.otago.ac.nz/density/pdfs/secr-habitatmasks.pdf

# Skjalfandi Bay study area:
x_coord <- c(537000, 582000, 582000, 537000)
y_coord <- c(655000, 655000, 607000, 607000)
xym <- cbind(x_coord, y_coord)
xym
p = Polygon(xym)
ps = Polygons(list(p),1)
skjal.poly.sq = SpatialPolygons(list(ps))

crs(skjal.poly.sq) <- "+proj=lcc +lat_1=64.25 +lat_2=65.75 +lat_0=65 +lon_0=-19 +x_0=500000 +y_0=500000 +ellps=GRS80 +units=m +no_defs" # crs for isn93
data <- data.frame(f = 99.9)
skjal.poly.sq <- SpatialPolygonsDataFrame(skjal.poly.sq, data)

# Upload polygon of Iceland outline. 
iceland <- st_read("data/iceland-isn93/is50v_strandlina_flakar_24122017.shp")

ggplot() +
  geom_sf(data=iceland, color = "black", fill = "gray")

# Now let's crop to skjalfandi
iceland <- as_Spatial(iceland)
crs(iceland)

skjal <- crop(iceland, extent(skjal.poly.sq))

skjal.plot <- st_as_sf(skjal)

ggplot() +
  geom_sf(data=skjal.plot, color = "black", fill = "gray") 

# Looks good! The polygon skjal can now be used to create a habitat mask

# Next, create the traps (necessary for spatially explicit cr)

n.cols <- (582000 - 537000)/1000
n.rows <- (655000 - 607000)/1000

skjal.traps <- make.grid(nx = n.cols, ny = n.rows, spacex = 1000, spacey = 1000, originxy = c(537000, 607000))

# Now, create the habitat mask

skjal.mask <- make.mask(skjal.traps, buffer = 12000, nx = n.cols, ny = n.rows, poly.habitat = FALSE, poly = skjal, keep.poly = TRUE)
# keep.poly = TRUE: retains polygon in this layer
# poly.habitat = FALSE: habitat lies outside the poylgon
