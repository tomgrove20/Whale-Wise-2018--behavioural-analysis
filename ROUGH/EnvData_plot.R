#####################################
### plotting environmental data #####
### Amelie Laute ####################
### 22.03.2020 ######################
#####################################


### explaining the plot coding

#layout
layout(matrix(1:10, ncol = 1),                 #how many plots in one
       widths = 1,                             #width of each plot
       heights = c(1.35,1,1,1,1,1,1,1,1,1.7),  #height of each plot (1st and last different due to margins)
       respect = FALSE)                        #?

#each plot
par(mar = c(0, 4.1, 0, 4.1))                                   #margin of this plot
plot(EnvData_2018$Time, EnvData_2018$`Gust [m/s]`,             #data source (x,y)
     type="h",                                                 #Histogram
     xaxt='n', yaxt='n',                                       #no x and y-Axis (numbers!)
     ylab="")                                                  #no y-Axis label
axis(side = 4, at = pretty(range(EnvData_2018$`Gust [m/s]`)))  #add y-Axis on the right instead (1=top, 2=left, 3=bottom, 4=right)
mtext("Gust [m/s]", side = 4, line = 2, cex=0.7)               #add y-Axis label on the right (line=distance to plot, cex=size of text)



############
### 2018 ###
############

library(readxl)
EnvData_2018 <- read_excel("Environmental_data/data_sheets/EnvData_2018.xlsx")
View(EnvData_2018)

# single plots

plot(EnvData_2018$Time, EnvData_2018$`Wind speed [m/s]`, type="h")
plot(EnvData_2018$Time, EnvData_2018$`Gust [m/s]`, type="h")
plot(EnvData_2018$Time, EnvData_2018$`Temperature [?C]`, type="h")
plot(EnvData_2018$Time, EnvData_2018$`Air pressure [hPa]`, type="h")
plot(EnvData_2018$Time, EnvData_2018$`Measured sea height [m]`, type="h")
plot(EnvData_2018$Time, EnvData_2018$`Precipitation [mm]`, type="h")
plot(EnvData_2018$Time, EnvData_2018$`Tide [m]`, type="h")
plot(EnvData_2018$Time, EnvData_2018$`Current direction [degrees]`, type="h")
plot(EnvData_2018$Time, EnvData_2018$`Current speed [m/s]`, type="h")
plot(EnvData_2018$Time, EnvData_2018$`Surge [m]`, type="h")
plot(EnvData_2018$Time, EnvData_2018$`Sea height [m]`, type="h")


#large plot including 10 rows

par(mar=c(6,6,4,4))   #not used

layout(matrix(1:10, ncol = 1), widths = 1, heights = c(1.35,1,1,1,1,1,1,1,1,1.7), respect = FALSE)

par(mar = c(0, 4.1, 2.1, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Wind speed [m/s]`, type="h", xaxt='n', ylab="")
title(main="Environmental data 2018")
axis(side = 2, at = pretty(range(EnvData_2018$`Wind speed [m/s]`))) 
mtext("Wind speed [m/s]", side = 2, line = 2.5, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Gust [m/s]`, type="h", xaxt='n', yaxt='n', ylab="")
axis(side = 4, at = pretty(range(EnvData_2018$`Gust [m/s]`))) 
mtext("Gust [m/s]", side = 4, line = 2, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Temperature [°C]`, type="h", xaxt='n', ylab="")
axis(side = 2, at = pretty(range(EnvData_2018$`Temperature [°C]`))) 
mtext("Temperature [°C]", side = 2, line = 2.5, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Air pressure [hPa]`, type="h", xaxt='n',  yaxt='n', ylab="")
axis(side = 4, at = pretty(range(EnvData_2018$`Air pressure [hPa]`))) 
mtext("Air pressure [hPa]", side = 4, line = 2, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Measured sea height [m]`, type="h", xaxt='n', ylab="")
axis(side = 2, at = pretty(range(EnvData_2018$`Measured sea height [m]`))) 
mtext("Measured sea height [m]", side = 2, line = 2.5, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Precipitation [mm]`, type="h", xaxt='n',  yaxt='n', ylab="")
axis(side = 4, at = pretty(range(EnvData_2018$`Precipitation [mm]`))) 
mtext("Precipitation [mm]", side = 4, line = 2, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Tide [m]`, type="h", xaxt='n', ylab="")
axis(side = 2, at = pretty(range(EnvData_2018$`Tide [m]`))) 
mtext("Tide [m]", side = 2, line = 2.5, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Current direction [degrees]`, type="h", xaxt='n', yaxt='n', ylab="")
axis(side = 4, at = pretty(range(EnvData_2018$`Current direction [degrees]`))) 
mtext("Current direction [degrees]", side = 4, line = 2, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Surge [m]`, type="h", xaxt='n', ylab="")
axis(side = 2, at = pretty(range(EnvData_2018$`Surge [m]`))) 
mtext("Surge [m]", side = 2, line = 2.5, cex=0.7)

par(mar = c(4.1, 4.1, 0, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Current speed [m/s]`, type="h", xlab="", yaxt='n', ylab="")
axis(side = 4, c(0, 0.005, 0.01)) 
mtext("Current speed [m/s]", side = 4, line = 2, cex=0.7)
mtext("Date", side = 1, line = 2, cex=0.9)


# not included (would be 11th)

par(mar = c(4.1, 4.1, 0, 4.1))
plot(EnvData_2018$Time, EnvData_2018$`Sea height [m]`, type="h", ylab="")
axis(side = 2, at = pretty(range(EnvData_2018$`Sea height [m]`))) 
mtext("Sea height [m]", side = 2, line = 2.5, cex=0.7)




############
### 2020 ###
############

EnvData_2020 <- read_excel("Environmental_data/data_sheets/EnvData_2020.xlsx")
View(EnvData_2020)

# single plots

plot(EnvData_2020$Time, EnvData_2020$`Wind speed [m/s]`, type="h")
plot(EnvData_2020$Time, EnvData_2020$`Gust [m/s]`, type="h")
plot(EnvData_2020$Time, EnvData_2020$`Temperature [?C]`, type="h")
plot(EnvData_2020$Time, EnvData_2020$`Air pressure [hPa]`, type="h")
plot(EnvData_2020$Time, EnvData_2020$`Measured sea height [m]`, type="h")
plot(EnvData_2020$Time, EnvData_2020$`Precipitation [mm]`, type="h")
plot(EnvData_2020$Time, EnvData_2020$`Tide [m]`, type="h")
plot(EnvData_2020$Time, EnvData_2020$`Current direction [degrees]`, type="h")
plot(EnvData_2020$Time, EnvData_2020$`Current speed [m/s]`, type="h")
plot(EnvData_2020$Time, EnvData_2020$`Surge [m]`, type="h")
plot(EnvData_2020$Time, EnvData_2020$`Sea height [m]`, type="h")

#large plot including 10 rows

par(mar=c(6,6,4,4))

layout(matrix(1:10, ncol = 1), widths = 1, heights = c(1.35,1,1,1,1,1,1,1,1,1.7), respect = FALSE)

par(mar = c(0, 4.1, 2.1, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Wind speed [m/s]`, type="h", xaxt='n', ylab="")
title(main="Environmental data 2020")
axis(side = 2, at = pretty(range(EnvData_2020$`Wind speed [m/s]`))) 
mtext("Wind speed [m/s]", side = 2, line = 2.5, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Gust [m/s]`, type="h", xaxt='n', yaxt='n', ylab="")
axis(side = 4, at = pretty(range(EnvData_2020$`Gust [m/s]`))) 
mtext("Gust [m/s]", side = 4, line = 2, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Temperature [°C]`, type="h", xaxt='n', ylab="")
axis(side = 2, at = pretty(range(EnvData_2020$`Temperature [°C]`))) 
mtext("Temperature [°C]", side = 2, line = 2.5, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Air pressure [hPa]`, type="h", xaxt='n',  yaxt='n', ylab="")
axis(side = 4, at = pretty(range(EnvData_2020$`Air pressure [hPa]`))) 
mtext("Air pressure [hPa]", side = 4, line = 2, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Measured sea height [m]`, type="h", xaxt='n', ylab="")
axis(side = 2, at = pretty(range(EnvData_2020$`Measured sea height [m]`))) 
mtext("Measured sea height [m]", side = 2, line = 2.5, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Precipitation [mm]`, type="h", xaxt='n',  yaxt='n', ylab="")
axis(side = 4, at = pretty(range(EnvData_2020$`Precipitation [mm]`))) 
mtext("Precipitation [mm]", side = 4, line = 2, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Tide [m]`, type="h", xaxt='n', ylab="")
axis(side = 2, at = pretty(range(EnvData_2020$`Tide [m]`))) 
mtext("Tide [m]", side = 2, line = 2.5, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Current direction [degrees]`, type="h", xaxt='n', yaxt='n', ylab="")
axis(side = 4, c(0, 150, 300)) 
mtext("Current direction [degrees]", side = 4, line = 2, cex=0.7)

par(mar = c(0, 4.1, 0, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Surge [m]`, type="h", xaxt='n', ylab="")
axis(side = 2, at = pretty(range(EnvData_2020$`Surge [m]`))) 
mtext("Surge [m]", side = 2, line = 2.5, cex=0.7)

par(mar = c(4.1, 4.1, 0, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Current speed [m/s]`, type="h", xlab="", yaxt='n', ylab="")
axis(side = 4, c(0, 0.01, 0.02)) 
mtext("Current speed [m/s]", side = 4, line = 2, cex=0.7)
mtext("Date", side = 1, line = 2, cex=0.9)



# not included (would be 11th)

par(mar = c(4.1, 4.1, 0, 4.1))
plot(EnvData_2020$Time, EnvData_2020$`Sea height [m]`, type="h", ylab="")
axis(side = 2, at = pretty(range(EnvData_2020$`Sea height [m]`))) 
mtext("Sea height [m]", side = 2, line = 2.5, cex=0.7)

