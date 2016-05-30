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
require(opentraj)
require(doParallel)
require(rgdal)
require(sp)

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

# Analysis:
# For each datapoint, calculate the 48hr backtrajectories (Hysplit) and add the "land covered" to get a total of "land influenced airmass" for each datapoint ... explore other metrics.

# Load Land use data from LUCAS
land_use <- readOGR("/home/gustavo/data/TF5_JapanNZvoyage/gis_layers/lucas-nz-land-use-map-1990-2008-2012-v016","lucas-nz-land-use-map-1990-2008-2012-v016")

# Restrict analysis to "near NZ" (26 to 28 Feb)
data_for_backtrajectories <- subset(march_data.10min, subset = (date>as.POSIXct("2013-02-26 04:30:00", tz = "UTC") & date<as.POSIXct("2013-03-01 00:00:00", tz = "UTC")))
data_for_backtrajectories$forest <- NA*data_for_backtrajectories$N10
data_for_backtrajectories$grassland <- data_for_backtrajectories$forest
data_for_backtrajectories$cropland <- data_for_backtrajectories$forest
data_for_backtrajectories$wetland <- data_for_backtrajectories$forest
data_for_backtrajectories$settlements <- data_for_backtrajectories$forest
data_for_backtrajectories$otherland <- data_for_backtrajectories$forest

# Summary plots ####
# Position of the ship
centreLat <- mean(data_for_backtrajectories$Latitude,na.rm = TRUE)
centreLon <- mean(data_for_backtrajectories$Longitude,na.rm = TRUE)
map <- get_map(location = c(centreLon,centreLat),zoom  = 4, maptype = "terrain")
ggmap(map) + 
  geom_point(aes(x=Longitude,y=Latitude,colour = log10(n265)),size = 3,data = data_for_backtrajectories, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")

n_points <- length(data_for_backtrajectories$date)
for (point_nr in (1:n_points)) {
  sample_traj <- ProcTraj(data_for_backtrajectories$Latitude[point_nr],
                          data_for_backtrajectories$Longitude[point_nr],
                          hour.interval = 1,
                          name = "SOAP1",
                          start.hour = format(data_for_backtrajectories$date[point_nr],format="%H:00"),
                          end.hour = format(data_for_backtrajectories$date[point_nr],format="%H:00"),
                          '~/data/hysplit/trunk/working/',
                          '~/repositories/TF5_analysis/',
                          hours = -48, height = 50,
                          '~/data/hysplit/trunk',
                          dates = format(data_for_backtrajectories$date[point_nr],format="%Y-%m-%d"))
  # extract lon-lat pairs from trajectory
  xy = sample_traj[,c(8,7)]
  # convert trajectory dataframe into "spatial point data frame"
  traj_spdf <- SpatialPointsDataFrame(coords = xy,
                                      data = sample_traj,
                                      proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  land_use_for_trajectories <- over(traj_spdf,land_use)
  data_for_backtrajectories$forest[point_nr] <- sum((land_use_for_trajectories$LUC_ID <=73 &  land_use_for_trajectories$LUC_ID >=71),na.rm = TRUE)
  data_for_backtrajectories$grassland[point_nr] <- sum((land_use_for_trajectories$LUC_ID <=76 &  land_use_for_trajectories$LUC_ID >=74),na.rm = TRUE)
  data_for_backtrajectories$cropland[point_nr] <- sum((land_use_for_trajectories$LUC_ID <=78 &  land_use_for_trajectories$LUC_ID >=77),na.rm = TRUE)
  data_for_backtrajectories$wetland[point_nr] <- sum((land_use_for_trajectories$LUC_ID <=80 &  land_use_for_trajectories$LUC_ID >=79),na.rm = TRUE)
  data_for_backtrajectories$settlements[point_nr] <- sum((land_use_for_trajectories$LUC_ID ==81),na.rm = TRUE)
  data_for_backtrajectories$otherland[point_nr] <- sum((land_use_for_trajectories$LUC_ID ==82),na.rm = TRUE)
}

data_for_backtrajectories$any_land <- with(data_for_backtrajectories,forest+grassland+cropland+wetland+settlements+otherland)

ggmap(map) + 
  geom_point(aes(x=Longitude,y=Latitude,colour = otherland),size = 3,data = data_for_backtrajectories, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")
