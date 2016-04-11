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
gps.file <- "gps_data_dat.txt"
# Grimm data
grimm.data <- read.csv(paste0(data.path,grimm.file),sep = "")
grimm.data$date <- ISOdatetime(year = grimm.data$Year+2000,
                              month = grimm.data$Month, 
                              day = grimm.data$Day,
                              hour = grimm.data$Hour,
                              min = grimm.data$Minute, 
                              sec = grimm.data$Second, tz = "UTC")

# GPS data
gps.data <- read.csv(paste0(data.path,gps.file), sep = "", stringsAsFactors = FALSE)
gps.data$date <- as.POSIXct(paste(gps.data$Date,gps.data$Time),format = "%d/%m/%Y %H:%M:%S",tz = 'UTC')

# Merge the data ####
march_data <- merge(gps.data,grimm.data,by = 'date', all = TRUE)
march_data.1min <- timeAverage(march_data,avg.time = '1 min')
march_data.5min <- timeAverage(march_data,avg.time = '5 min')

# Summary plots ####
# Position of the ship
newmap <- getMap(resolution = "low")
plot(newmap,xlim = c(160,175), ylim = c(-45,0))
points(march_data.5min$Longitude,march_data.5min$Latitude,col = 'red', cex = .6)

map <- get_map(location = c(170,0),zoom  = 3, maptype = "terrain")
ggmap(map) + 
  geom_point(aes(x=Longitude,y=Latitude,colour = log10(n265)),size = 3,data = march_data.1min, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")
