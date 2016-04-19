## TF5 Grimm data analysis #
# Load packages ####
require(ncdf4)
require(openair)
require(ggplot2)
require(rworldmap)
require(ggmap)
require(mapproj)
require(corrgram)
require(corrplot)

# Load the data and create date fields ####
# Where are the data
data.path <- "/home/gustavo/data/TF5_JapanNZvoyage/2013_Feb/"
grimm.file <- "GrimmDump.txt.parsed"
gps.file <- "gps_cpc.txt"
# Grimm data
grimm.data <- read.csv(paste0(data.path,grimm.file),sep = "")
grimm.data$date <- ISOdatetime(year = grimm.data$Year+2000,
                              month = grimm.data$Month, 
                              day = grimm.data$Day,
                              hour = grimm.data$Hour,
                              min = grimm.data$Minute, 
                              sec = grimm.data$Second, tz = "UTC")
grimm.data$Year<-NULL
grimm.data$Month<-NULL
grimm.data$Day<-NULL
grimm.data$Hour<-NULL
grimm.data$Minute<-NULL
grimm.data$Second<-NULL
grimm.data$Error<-NULL
dp=c(265,290,324,374,424,474,539,614,675,748,894,1140,1442,1789,2236,2739,3240,3742,4472,5701,6982,7984,9220,11180,13693,16202,18708,22361,27386,30984,35000)
sp = 4*pi*(dp/2)^2
vp = (4/3)*pi*(dp/2)^3
dlogDp=c(0.0492180227, 0.0299632234, 0.0669467896, 0.057991947, 0.0511525224, 0.0457574906, 0.0644579892, 0.0494853631, 0.0321846834,  0.057991947, 0.096910013, 0.1139433523,  0.0901766303,  0.096910013, 0.096910013, 0.079181246, 0.0669467896,  0.057991947, 0.096910013, 0.1139433523,  0.0621479067,  0.0543576623,  0.0705810743,  0.096910013, 0.079181246, 0.0669467896,  0.057991947, 0.096910013, 0.079181246, 0.0280287236,  1)
grimm.data$Ngrimm <- rowSums(grimm.data[,1:31]*dlogDp)*10/(60*1.5)

# GPS data
gps.data <- read.csv(paste0(data.path,gps.file),sep = "", stringsAsFactors = FALSE)
gps.data$date <- as.POSIXct(paste(gps.data$Date,gps.data$Time),format = "%d/%m/%Y %H:%M:%S",tz = 'UTC')
gps.data$Date <- NULL
gps.data$Time <- NULL

# Merge the data ####
march_data <- merge(gps.data,grimm.data,by = 'date', all = TRUE)
march_data.1min <- timeAverage(march_data,avg.time = '1 min')
march_data.10min <- timeAverage(march_data,avg.time = '10 min')
march_data.30min <- timeAverage(march_data,avg.time = '30 min')

# Move to dN/dLogDp
dn10.1min <- (march_data.1min$N10 - march_data.1min$Ngrimm)/log10(265/10)
dn10.10min <- (march_data.10min$N10 - march_data.10min$Ngrimm)/log10(265/10)
dn10.30min <- (march_data.30min$N10 - march_data.30min$Ngrimm)/log10(265/10)
march_data.1min$dN10 <- dn10.1min
march_data.10min$dN10 <- dn10.10min
march_data.30min$dN10 <- dn10.30min

# Summary plots ####
# Position of the ship
newmap <- getMap(resolution = "low")
plot(newmap,xlim = c(160,175), ylim = c(-45,0))
points(march_data.30min$Longitude,march_data.30min$Latitude,col = 'red', cex = .6)
centreLat <- mean(march_data.10min$Latitude,na.rm = TRUE)
centreLon <- mean(march_data.10min$Longitude,na.rm = TRUE)
centreLat <- -40.025
centreLon <- 172.964
map <- get_map(location = c(centreLon,centreLat),zoom  = 7, maptype = "terrain")
ggmap(map) + 
  geom_point(aes(x=Longitude,y=Latitude,colour = log10(n265)),size = 3,data = march_data.10min, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")

ggmap(map) + 
  geom_point(aes(x=Longitude,y=Latitude,colour = log10(N10)),size = 3,data = march_data.10min, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")

# Size distribution

# Before
size_dist_data <- subset(march_data.10min, subset = (date>as.POSIXct("2013-02-27 20:00:00", tz = "UTC") & date<as.POSIXct("2013-02-27 21:00:00", tz = "UTC")))[,c(5:35)]
ave_size_dist_before <- data.frame(dp = dp, dndlogdp = colMeans(size_dist_data,na.rm = TRUE)*sp)

size_dist_data <- subset(march_data.10min, subset = (date>as.POSIXct("2013-02-27 21:14:00", tz = "UTC") & date<as.POSIXct("2013-02-27 22:38:00", tz = "UTC")))[,c(5:35)]
ave_size_dist_during <- data.frame(dp = dp, dndlogdp = colMeans(size_dist_data,na.rm = TRUE)*sp)

size_dist_data <- subset(march_data.10min, subset = (date>as.POSIXct("2013-02-27 22:40:00", tz = "UTC") & date<as.POSIXct("2013-02-27 23:52:00", tz = "UTC")))[,c(5:35)]
ave_size_dist_after <- data.frame(dp = dp, dndlogdp = colMeans(size_dist_data,na.rm = TRUE)*sp)

#size_dist_data <- subset(march_data.10min, subset = (date>as.POSIXct("2013-02-27 04:53:00", tz = "UTC") & date<as.POSIXct("2013-02-27 05:52:00", tz = "UTC")))[,c(5:35)]
#ave_size_dist_after2 <- data.frame(dp = dp, dndlogdp = colMeans(size_dist_data,na.rm = TRUE)*vp)

ggplot(ave_size_dist_before,aes(x=dp,y=dndlogdp))+
  geom_line(aes(colour = 'before'))+
  geom_line(aes(y=ave_size_dist_during$dndlogdp,colour = 'during'))+
  geom_line(aes(y=ave_size_dist_after$dndlogdp,colour = 'after'))+
#  geom_line(aes(y=ave_size_dist_after2$dndlogdp,colour = 'after2'))+
  #scale_color_discrete()+
  scale_x_log10()+
  scale_y_log10()

# During
size_dist_data <- timeAverage(march_data.10min, start.date = "", end.date = "")
ave_size_dist <- data.frame(dp = c(51,dp), dndlogdp = colMeans(size_dist_data,na.rm = TRUE))
ggplot(ave_size_dist,aes(x=dp,y=dndlogdp))+
  geom_point(colour = 'red')+
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 6))+
  scale_x_log10()+
  scale_y_log10()

# After
size_dist_data <- timeAverage(march_data.10min, start.date = "", end.date = "")
ave_size_dist <- data.frame(dp = c(51,dp), dndlogdp = colMeans(size_dist_data,na.rm = TRUE))
ggplot(ave_size_dist,aes(x=dp,y=dndlogdp))+
  geom_point(colour = 'red')+
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 6))+
  scale_x_log10()+
  scale_y_log10()



# Correlation to identify aerosol modes

for_correl <- march_data.10min[,c(37,5:35)]
find_channels <- cor(for_correl,use = 'pairwise')

corrgram(find_channels)
corrplot(find_channels,method = 'circle')


# Analysis:
# For each datapoint, calculate the 48hr backtrajectories (Hysplit) and add the "land covered" to get a total of "land influenced airmass" for each datapoint ... explore other metrics.

