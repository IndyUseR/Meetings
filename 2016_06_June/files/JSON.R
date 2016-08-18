pullWundergroundJSON <- function(CityCode,CityName){
  library("jsonlite")
  
  # CityCode can be a 4-digit airport ID, zipcode, State/City combo
  # Wunderground no. of calls per day: http://www.wunderground.com/weather/api/d/832f2cc267c8aa47/edit.html?error=norecurly
  
  # Concatenates the pieces of the Wunderground URL with user input
  vectorOfText <- c("http://api.wunderground.com/api/832f2cc267c8aa47/hourly10day/q/",CityCode,".json")
  text <- paste(vectorOfText, collapse="")
  forecast.df <- as.data.frame(fromJSON(txt = text))# grab wunderground data from web and assigns it
  
  #Extracts just the datasets needed
  output <- cbind(forecast.df$hourly_forecast.FCTTIME,forecast.df$hourly_forecast.temp, forecast.df$hourly_forecast.dewpoint, forecast.df$hourly_forecast.condition, forecast.df$hourly_forecast.pop)
  
  # Exclude unneeded columns
  output <- output[, -c(2:5,8,11:18,20:26,27,29,31)]

  
  filename <- paste(CityName,"_JSON.csv")# this pastes the city in each csv 
  
  write.csv(output, file = filename)# creates a csv file in the specified working directory. 

}