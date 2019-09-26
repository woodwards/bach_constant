# base flow analysis following "lfstat" method
# Gustard, A. & Demuth, S. (2009) (Eds) Manual on Low-flow Estimation and Prediction. 
# Operational Hydrology Report No. 50, WMO-No. 1029, 136p.

library(tidyverse)
library(lfstat)
library(lubridate)
library(FlowScreen)

out_path <- "run - eckhardt_priors_narrow/"
source("read_data8.r")

startcalib <- min(options$startcalib)
endcalib <- max(options$endcalib)
nyears <- 15
catchments <- unique(data$file)
i <- 1
for (i in 1:8){
  catchment <- catchments[[i]]
  print(catchment)
  cdata <- data %>%
    filter(file==catchment) %>%
    mutate(row=date-min(date)+1) %>%
    filter(row>=startcalib & row<=endcalib) %>%
    mutate(
      date2=as.POSIXct(date*(60*60*24), origin="1899-12-30", tz="GMT"),
      date3=as.Date(date, origin="1899-12-30", tz="GMT"),
      day=day(date3),
      month=month(date3),
      year=year(date3)
    )
  ldata <- createlfobj(cdata, hyearstart=1)
   # bfi <- BFI(ldata)   
   # print(paste("BFI Fixed    =", bfi))
  a <- recession(ldata, method="MRC", seglen = 7, threshold = 70)
  a <- exp(-1/a)
  print(paste("Recession =", a))
  # cdata$baseflow2 <- bf_eckhardt(cdata$flow, 0.98, 0.8)
  cdata$baseflow2 <- bf_eckhardt(cdata$flow, a, 0.8)
  bfi2 <- sum(cdata$baseflow2) / sum(cdata$flow)
  print(paste("BFI Eckhardt 0.8 =", bfi2))
  # rain <- sum(cdata$rain)/15
   # print(paste("rain =", rain))
   
}


