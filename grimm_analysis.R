## TF5 Grimm data analysis #
# Load packages ####
require(ncdf4)
require(openair)
require(ggplot2)
require(rworldmap)
require(ggmap)
require(mapproj)

# Load the data and create date fields ####
# Where are the data
data.path <- "/home/gustavo/data/TF5_JapanNZvoyage/2013_Feb/"
grimm.file <- "GrimmDump.txt.parsed"
gps.file <- "gps_cn.csv"
# Grimm data
grimm.data <- read.csv(paste0(data.path,grimm.file),sep = "")
grimm.data$date <- ISOdatetime(year = grimm.data$Year+2000,
                              month = grimm.data$Month, 
                              day = grimm.data$Day,
                              hour = grimm.data$Hour,
                              min = grimm.data$Minute, 
                              sec = grimm.data$Second, tz = "UTC")
dp=c(265,290,324,374,424,474,539,614,675,748,894,1140,1442,1789,2236,2739,3240,3742,4472,5701,6982,7984,9220,11180,13693,16202,18708,22361,27386,30984,35000)
dlogDp=c(0.0492180227, 0.0299632234, 0.0669467896, 0.057991947, 0.0511525224, 0.0457574906, 0.0644579892, 0.0494853631, 0.0321846834,  0.057991947, 0.096910013, 0.1139433523,  0.0901766303,  0.096910013, 0.096910013, 0.079181246, 0.0669467896,  0.057991947, 0.096910013, 0.1139433523,  0.0621479067,  0.0543576623,  0.0705810743,  0.096910013, 0.079181246, 0.0669467896,  0.057991947, 0.096910013, 0.079181246, 0.0280287236,  1)

grimm.dndlog <- grimm.data[,names(grimm.data)[17:48]]
grimm.dndlog[,1:31] <- (grimm.data[,17:47] - c(grimm.data[,18:47],0))/dlogDp

# GPS data
gps.data <- read.csv(paste0(data.path,gps.file), stringsAsFactors = FALSE)
gps.data$date <- as.POSIXct(paste(gps.data$Date,gps.data$Time),format = "%d/%m/%Y %H:%M:%S",tz = 'UTC')

# Merge the data ####
march_data <- merge(gps.data,grimm.dndlog,by = 'date', all = TRUE)
march_data.1min <- timeAverage(march_data,avg.time = '1 min')
march_data.5min <- timeAverage(march_data,avg.time = '5 min')

# Summary plots ####
# Position of the ship
newmap <- getMap(resolution = "low")
plot(newmap,xlim = c(160,175), ylim = c(-45,0))
points(march_data.5min$Longitude,march_data.5min$Latitude,col = 'red', cex = .6)
centreLat <- mean(march_data.1min$Latitude,na.rm = TRUE)
centreLon <- mean(march_data.1min$Longitude,na.rm = TRUE)

map <- get_map(location = c(centreLat,centreLon),zoom  = 3, maptype = "terrain")
ggmap(map) + 
  geom_point(aes(x=Longitude,y=Latitude,colour = log10(n265)),size = 3,data = march_data.1min, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")

# Analysis:
# For each datapoint, calculate the 48hr backtrajectories (Hysplit) and add the "land covered" to get a total of "land influenced airmass" for each datapoint ... explore other metrics.

