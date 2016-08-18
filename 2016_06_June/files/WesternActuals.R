#westernActuals

  
  date <- as.Date(Sys.time())
  date.minus1 <- date - 1       # to get yesterday's date
  month <- format(date.minus1, "%m")
  day <- format(date.minus1, "%d") 
  #day <- format(date - 1, "%d") 
  year <- format(date.minus1,"%Y")    # 4-digit year

  
  file.date <- c(year, month, day)
  url.date <- paste(file.date, collapse = "-")
  
  
  #calls weatherData package to active so script can run, will not run if package is not activated
  library("weatherData")
  
  #station ID's are defined in a list 
  ID <- c("TUS","PHX","LAS","LAX","SAC","PSP")
  
  #temperature variables are defined before start of the loop as numeric vectors with the same length 
  #as the ID vector. aka however many stations we have called, in this case 6. 
  tempMin <- numeric(length(ID))
  tempMax <- numeric(length(ID))
  
  # i loops over the length of ID which in this case is 6. So the script loops 6 times through the 
  #getDailyMinMaxTemp function. 
  for(i in 1:length(ID)){
    temp <- getSummarizedWeather(ID[i], url.date)
    tempMin[i] = temp[4] #allows 16 min temps to be produced, min temp is in the second column of each row
    tempMax[i] = temp[2] #allows 16 max temps to be produced, max temp is in the fourth column of each row
   }
  
  #cbind combinds three seperate data columns created above into one file which is named output 
  output <- cbind(ID,tempMin,tempMax)
  filename <- paste("Western Actuals",url.date,".csv")# this pastes the date in each csv so we don't lose track of data
  
  write.csv(output, file = filename)# creates a csv file in the specified working directory. 

