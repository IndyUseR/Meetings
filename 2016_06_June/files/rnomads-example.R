library("rNOMADS")
library("XLConnect")
library("dplyr")
library("sendmailR")

GFS.list <- GetDODSDates(abbrev="gfs_0p50", archive=F, request.sleep=1) #returns list of available dates with GFS data in a numerical list

current.GFS.url <- tail(GFS.list$url,1)
GFS.model.runs <- GetDODSModelRuns(current.GFS.url)

GFS.12z <- tail(GFS.model.runs$model.run,1) # pulls the actual 12z GFS model instead of the analysis 

#run when needed
#GFS.info.variables <- GetDODSModelRunInfo(current.GFS.url, GFS.12z) # view all weather variables the model has to offer 

lon <- c(480, 620) # USA Longitude
lat <- c(240, 280) # USA Latitude

variables <- c("tmax2m","tmin2m", "weasdsfc", "ugrd10m", "vgrd10m")

GFS.Model.Data <- DODSGrab(current.GFS.url, GFS.12z, variables, time = c(6,14), lon, lat)


###################### read in station data with lat-lon #################################

stationlist <- loadWorkbook(filename = "StationList.xls") #read in station data
stationlist <- readWorksheet(stationlist, sheet = 1, startCol = 1, startRow = 1, endCol = 4, endRow = 232)#make it into a dataframe 

#GFS.test <- matrix(data=NA, nrow= length(stationlist$Station.List)*9,ncol = length(variables))  #nrows is number of stations mutliplied by the .Rdata files(9), ncols in # of weather variables

Lat <- as.numeric(stationlist$Lat)
Lon <- as.numeric(stationlist$Lon)


#### New Buildprofile function completely circumvents having to use ModelGrid
GFS.test <- BuildProfile(GFS.Model.Data, Lon, Lat, spatial.average = TRUE, points = 4)


# write code to unpack data from the big list

#GFS.test <-  matrix(unlist(GFS.test[[1]]$profile.data),ncol=5,byrow=TRUE) # take the 1 and turn it into an i into a for loop 

GFS.data <- NULL

#matrix(data=NA, nrow= length(stationlist$Station.List)*9,ncol = length(variables))

for( i in 1:length(GFS.test)) {
  
  GFS.data1 <-  matrix(unlist(GFS.test[[i]]$profile.data),ncol=5,byrow=TRUE)
  GFS.data <- rbind(GFS.data1, GFS.data)
}


GFS.data <- data.frame(GFS.data) #change the data from a matrix format to a data frame

colnames(GFS.data) <- variables #add column names to the data frame with all the model data 

###################### DPLYR Section ##############################

GFS.data <- tbl_df(GFS.data)# attach data frame to local environment so row/column names are easier to type

GFS.data <- GFS.data %>%  ## add srf.wind to the data frame by using pythagorean theorem
  select(tmax2m, tmin2m, weasdsfc, ugrd10m, vgrd10m) %>%
  mutate(sfc.wind= sqrt((ugrd10m^2) +(vgrd10m^2))*2.24)

GFS.data <- GFS.data[,-(4:5)]# remove u and v wind columns 

GFS.data  <- GFS.data %>% #round surface snowfall to look 3 decimal number to look cleaner
  select(tmax2m, tmin2m, weasdsfc, sfc.wind) %>%
  mutate(snowfall=(round(weasdsfc, digits = 3)))%>%
  mutate(snowfall= ifelse(snowfall<0,0,snowfall))%>%   # take negative values away
  mutate(snowfall= (snowfall*.039370))# convert to inches 

GFS.data <- GFS.data[,-3]# remove old snowfall column

GFS.data<- GFS.data %>% #convert each temp row to fahrenheit 
  select(tmax2m, tmin2m, snowfall, sfc.wind) %>%
  mutate(Tmax= 1.8*(tmax2m-273)+32) %>%
  mutate(Tmin = 1.8*(tmin2m-273)+32)

GFS.data <- GFS.data[,-(1:2)] #remove old Kelvin columns 

#GFS.data <- GFS.data %>% 
#select(Tmax, Tmin, snowfall, sfc.wind) %>%
#mutate(Tavg= (Tmax+Tmin)/2)

#GFS.data <- GFS.data[,-(1:2)] #remove high and low temp columns

stationlist9 <- as.character(rep(stationlist$Station.List, each= 9)) #repeat each station name for each of the 9 forecast invervals in the data frame

GFS.data <- cbind(stationlist9, GFS.data)# bind the row names and NWP data frame 

GFS.data<- GFS.data %>% group_by(stationlist9) %>% # use the dplyr package to summarize the whole data frame by temp, snowfall, and wind
  summarize(snowfall = ((max(snowfall)-min(snowfall))*10), wind = max(sfc.wind), maxT =max(Tmax), minT= min(Tmin), avgT= (colMeans(rbind(maxT,minT)))) %>%
  mutate(HDD = (65-avgT)) %>% 
  mutate(Fantasy.Score= (snowfall*10) + HDD + wind)

GFS.data <- GFS.data[,-(4:6)]# remove old temp columns. 

dataframe.order <- order(-xtfrm(GFS.data$Fantasy.Score))

GFS.data <- GFS.data[dataframe.order,]

write.csv(GFS.data, "fantasyblizzardforecastDA.csv")
