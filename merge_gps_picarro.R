# Picarro-GPS merging ####
# Author: Gustavo Olivares
# Date: 2016-07-20
# License: GPLv3
##########################

# Load Libraries ####
library('ggplot2')
library('openair')
# Set data filepaths ####
base_path_2012 <- '/home/gustavo/data/GPS_PICARRO/Jun2012/'
base_path_2013 <- '/home/gustavo/data/GPS_PICARRO/May2013/'
# 2012 data processing ####
# Load data ####
gps <- read.csv(paste0(base_path_2012,"gps_201206.csv"), sep=",")
names(gps) <- c('Day','Time','lat','lon')
gps$date <- as.POSIXct(paste(gps$Day,gps$Time),format = '%d/%m/%Y %H:%M:%OS', tz='UTC')
picarro <- read.csv(paste0(base_path_2012,"picarro_201206.txt"), sep="")
picarro$date <- as.POSIXct(paste(picarro$DATE,picarro$TIME),format = '%d/%m/%y %H:%M:%OS', tz='UTC')

# Merge data ####
merged_data <- merge(gps,picarro,by = 'date', all = TRUE)
# Patch GPS data ####
# The beginning
merged_data$lat[1] <- gps$lat[1]
merged_data$lon[1] <- gps$lon[1]
# The end
n_merged <- length(merged_data$date)
n_gps <- length(gps$date)
merged_data$lat[n_merged] <- gps$lat[n_gps]
merged_data$lon[n_merged] <- gps$lon[n_gps]

# Resample to 1 minutes ####
june2012.1minute <- timeAverage(merged_data,avg.time = '1 min')
write.table(june2012.1minute,file = paste0(base_path_2012,'june2012_1min.csv'),sep = '\t',row.names = FALSE)

# 2013 data processing ####
# Load data ####
gps <- read.csv(paste0(base_path_2013,"gps_201305.csv"), sep=",")
names(gps) <- c('Day','Time','lat','lon')
gps$date <- as.POSIXct(paste(gps$Day,gps$Time),format = '%d/%m/%Y %H:%M:%OS', tz='UTC')
picarro <- read.csv(paste0(base_path_2013,"picarro_201305.txt"), sep="")
picarro$date <- as.POSIXct(paste(picarro$DATE,picarro$TIME),format = '%d/%m/%y %H:%M:%OS', tz='UTC')

# Merge data ####
merged_data <- merge(gps,picarro,by = 'date', all = TRUE)
# Patch GPS data ####
# The beginning
merged_data$lat[1] <- gps$lat[1]
merged_data$lon[1] <- gps$lon[1]
# The end
n_merged <- length(merged_data$date)
n_gps <- length(gps$date)
merged_data$lat[n_merged] <- gps$lat[n_gps]
merged_data$lon[n_merged] <- gps$lon[n_gps]

# Resample to 1 minutes ####
may2013.1minute <- timeAverage(merged_data,avg.time = '1 min')
write.table(may2013.1minute,file = paste0(base_path_2013,'may2013_1min.csv'),sep = '\t',row.names = FALSE)
