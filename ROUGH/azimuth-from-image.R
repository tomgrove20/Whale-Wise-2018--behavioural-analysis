### Measuring whale angle using a reference object
### Tom Grove
### 17/04/2020

# Often, when observing whales, we fail to measure azimuth with an electronic range finder. In this case, it may still be possible to calculate azimuth from images containing a whale and a stable reference object whose position can be deduced (e.g. landform, building). 

# We need the following in a table to conduct this:
# - Observer location
# - Horizontal pixel distance between horizontal midpoint and whale midpoint
# - Horizontal pixel distance between horizontal midpoint and reference object
# - Lens/camera frame of view horizontal angle (requires knowledge of sensor size and focal length)
# - Image width (pixels)

# Equation taken from last page of the following document: http://128.148.32.110/courses/cs148/tutorials/range_bearing_estimation.pdf
# - I(x) = horizontal distance between midpoint and object
# - V(x) = image width (pixels) / 2. right > 0, left < 0
# - I(theta) = horizontal angle (azimuth) of object relative to camera orientation
# - V(theta) = frame of view angle / 2

# Packages
library(dplyr)
library(exifr)

# Let's use an example: observer at a known position in Skjalfandi Bay takes an photo, containing two reference objects: each 2D edge of Lundey.

# Image: 20180601_1300_andvari_77.JPG (in 'data')
# Image time: 14:33:27
# Observer location (from pocketGIS): 566583.613 E,623836.771 N (ISN93)
# Location 1 (left):
# Location 2 (right):
# Angle between two lines coming from the two points (from QGIS): 2.998 degrees

# Now let's compare with image measurements (in ImageJ):
(img.dat <- read_exif("data/20180601_1300_andvari_77.JPG", tags = c("ExifImageWidth", "FocalLength", "FOV" ))) # Note, unsure if FOV (field of view) is the same as AOV (angle of view)
image.width <- 5472
midpoint <- 5472/2
sensor.x <- 22.5 # for Canon 70D
p1.x <- 1830
p2.x <- 5459
dist.x <- p2.x - p1.x
aov <- (2*atan(sensor.x/(2*img.dat$FocalLength)))*(180/pi) # angle of view

# Calculating angle:
p1.p2.az <- (dist.x/image.width)*aov
# 2.998

# Clearly, these angles match up - perfect!

# CONCLUSION: the above approach will be suitable to calculate whale azimuth.