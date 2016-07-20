## TF5 Grimm data analysis #
# Load packages ####
library('ncdf4')
library('openair')
library('ggplot2')
library('rworldmap')
library('ggmap')
library('mapproj')
library('corrgram')
library('corrplot')
library('openair')
library('opentraj')
library('doParallel')
library('rgdal')
library('sp')
library('gridExtra')
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
all_march_data <- merge(gps.data,grimm.data,by = 'date', all = TRUE)
march_data <- subset(all_march_data,subset = (!is.na(all_march_data$Latitude)&(!is.na(all_march_data$Longitude))))
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
map <- get_map(location = c(centreLon,centreLat),zoom  = 7, maptype = "terrain")
ggmap(map) + 
  geom_point(aes(x=Longitude,y=Latitude,colour = log10(n265)),size = 3,data = march_data.10min, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")

# Analysis:
# For each datapoint, calculate the 48hr backtrajectories (Hysplit) and add the "land covered" to get a total of "land influenced airmass" for each datapoint ... explore other metrics.

# Load Land use data from LUCAS
#land_use <- readOGR("/home/gustavo/data/TF5_JapanNZvoyage/gis_layers/lucas-nz-land-use-map-1990-2008-2012-v016","lucas-nz-land-use-map-1990-2008-2012-v016")
land_use_data <- readGDAL("/home/gustavo/data/TF5_JapanNZvoyage/gis_layers/land_use_global_MOD12C1_T3.tif")

# Restrict analysis to 26 to 28 Feb ... when the ship was "near NZ"
data_for_backtrajectories <- subset(march_data.10min, subset = (date>as.POSIXct("2013-02-26 04:00:00", tz = "UTC") & date<as.POSIXct("2013-02-28 23:59:00", tz = "UTC")))
# Do not restrict analysis ... the whole trip!
data_for_backtrajectories <- march_data.10min
data_for_backtrajectories$land_use <- NA*data_for_backtrajectories$N10
data_for_backtrajectories$forest <- data_for_backtrajectories$land_use
data_for_backtrajectories$settlements <- data_for_backtrajectories$land_use
data_for_backtrajectories$other_veg <- data_for_backtrajectories$land_use
n_points <- length(data_for_backtrajectories$date)

# get the number of phisical cores availables
cores <- detectCores()
#
cl <- makeCluster(cores)

registerDoParallel(cl)
all_trajectories <- foreach(point_nr=1:n_points,.packages=c("opentraj"),.combine=rbind) %dopar%
  {
 #output.file.name<-""
 #output.file.name<-paste0("traj", "_", as.character(point_nr), "_")
  ProcTraj(data_for_backtrajectories$Latitude[point_nr],
           data_for_backtrajectories$Longitude[point_nr],
           hour.interval = 1,
           name = 'output.file.name',
           start.hour = format(data_for_backtrajectories$date[point_nr],format="%H:00"),
           end.hour = format(data_for_backtrajectories$date[point_nr],format="%H:00"),
           '~/data/hysplit/trunk/working/',
           '~/repositories/TF5_analysis/',
           hours = -96, height = 30,
           '~/data/hysplit/trunk',
           ID = point_nr,
           dates = format(data_for_backtrajectories$date[point_nr],format="%Y-%m-%d"),
           clean.files = TRUE)
  }
dump_output <- foreach(point_nr = 1:n_points,.packages = c("sp")) %dopar%  {
  all_trajectories[(((point_nr-1)*49+1):((point_nr)*49)),1] <- point_nr
  sample_traj <- all_trajectories[(((point_nr-1)*49+1):((point_nr)*49)),]
  xy <- sample_traj[,c(8,7)]
  traj_spdf <- SpatialPointsDataFrame(coords = xy,
                                      data = sample_traj,
                                      proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  land_use_for_trajectories <- over(traj_spdf,land_use_data)
  land_use <- sum(land_use_for_trajectories$band1 >0,na.rm = TRUE)
  forest <- sum(land_use_for_trajectories$band1 >=5 &
                                                      land_use_for_trajectories$band1 <= 8 ,na.rm = TRUE)
  other_veg <- sum(land_use_for_trajectories$band1 >=1 &
                                                        land_use_for_trajectories$band1 <= 4 ,na.rm = TRUE)
  settlements <- sum(land_use_for_trajectories$band1 ==10 ,na.rm = TRUE)
  c(land_use, forest, other_veg, settlements)
}
for (point_nr in (1:n_points)) {
  data_for_backtrajectories$land_use[point_nr] <- dump_output[[point_nr]][1]
  data_for_backtrajectories$forest[point_nr] <- dump_output[[point_nr]][2]
  data_for_backtrajectories$other_veg[point_nr] <- dump_output[[point_nr]][3]
  data_for_backtrajectories$settlements[point_nr] <- dump_output[[point_nr]][4]
}
centreLat <- mean(data_for_backtrajectories$Latitude,na.rm = TRUE)
centreLon <- mean(data_for_backtrajectories$Longitude,na.rm = TRUE)
map <- get_map(location = c(centreLon,centreLat),zoom  = 3, maptype = "terrain")

p1 <- ggmap(map) + 
  ggtitle("log10(n265)")+
  geom_point(aes(x=Longitude,y=Latitude,colour = log10(n265)),size = 3,data = data_for_backtrajectories, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")
# ggsave("./map_01.pdf",width=20,height = 30, units = "cm")

p2 <- ggmap(map) + 
  ggtitle("Land")+
  geom_point(aes(x=Longitude,y=Latitude,colour = land_use),size = 3,data = data_for_backtrajectories, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")
# ggsave("./map_02.pdf",width=20,height = 30, units = "cm")

p3 <- ggmap(map) + 
  ggtitle("Forest")+
  geom_point(aes(x=Longitude,y=Latitude,colour = forest),size = 3,data = data_for_backtrajectories, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")
# ggsave("./map_03.pdf",width=20,height = 30, units = "cm")

p4 <- ggmap(map) + 
  ggtitle("Settlements")+
  geom_point(aes(x=Longitude,y=Latitude,colour = settlements),size = 3,data = data_for_backtrajectories, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")
# ggsave("./map_04.pdf",width=20,height = 30, units = "cm")

p5 <- ggmap(map) + 
  ggtitle("Other Vegetation")+
  geom_point(aes(x=Longitude,y=Latitude,colour = other_veg),size = 3,data = data_for_backtrajectories, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")
# ggsave("./map_05.pdf",width=20,height = 30, units = "cm")

centreLat <- mean(all_trajectories$lat,na.rm = TRUE)
centreLon <- mean(all_trajectories$lon,na.rm = TRUE)
map2 <- get_map(location = c(centreLon,centreLat),zoom  = 3, maptype = "terrain")
p6 <- ggmap(map2) + 
  ggtitle("All Trajectories")+
  geom_point(aes(x=lon,y=lat,colour = height),size = 3,data = all_trajectories, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")
# ggsave("./map_06.pdf",width=20,height = 30, units = "cm")


ggplot_list <- Filter(function(x) is(x, "ggplot"), mget(ls()))
ggsave("./All_plots_30m.pdf", marrangeGrob(grobs = ggplot_list, nrow=1, ncol=1))


timePlot(data_for_backtrajectories,pollutant = c('N10','n265','n3240','land_use','forest','other_veg','settlements'))

# Explore certain trajectories ####
datapoint <- n_points - 5
sample_traj <- all_trajectories[(((point_nr-1)*49+1):((point_nr)*49)),]
sample_traj <- all_trajectories
centreLat <- mean(sample_traj$lat,na.rm = TRUE)
centreLon <- mean(sample_traj$lon,na.rm = TRUE)
map2 <- get_map(location = c(centreLon,centreLat),zoom  = 3, maptype = "terrain")
ggmap(map2) + 
  ggtitle("All Trajectories")+
  geom_point(aes(x=lon,y=lat,colour = receptor),size = 3,data = sample_traj, alpha = .3) +
  scale_colour_gradient(low = "white",high = "red")



