rm(list = ls())
cat("\014")

library(raster)
library(ncdf4)
"%>%" = magrittr::'%>%'

setwd("C:/NDG/duque-gardeazabal/pages")

## download current monthly value ----
ersst_url <- "https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/netcdf/"

list.files(pattern = "\\.nc") %>% file.remove()

today <- Sys.time() %>% as.POSIXlt()
today$mon <- today$mon -1
f.today <- substr(today,1,7) %>% gsub("-","",.) %>% paste0("ersst.v5.",.,".nc")

download.file(paste0(ersst_url,f.today), 
              "SST_today.nc",
              mode = "wb") # important command in windows

SSTa <- brick("SST_today.nc", varname="ssta")


## calculate AMM & Atl3 ----

lat_Weight.Mean <- function(R){
  aux <- rasterToPoints(R)
  R.weight <- ( aux[,3] *sqrt(cos(3.14159 * aux[,2]/180)) )%>%
    mean(., na.rm=TRUE )
  return( R.weight)
} # Special function for calculated the mean with latitude weighting, also for the two datasets

E.sst.AMM <- lat_Weight.Mean( crop(SSTa, extent(c(360 -70 , 360 -15 ,5,25)))) - 
  lat_Weight.Mean( crop(SSTa, extent(c(360 -40, 360 -0,-25,-5))))

E.sst.AEN <- lat_Weight.Mean(crop(SSTa, extent(c(360-20, 360-0 ,-3,3))))

write.csv(rbind(E.sst.AMM, E.sst.AEN), "Current_AMM_Atl3.csv", col.names = F, row.names = F)

## update value on webpage ----
ch.file.Atl <- function(x){
  a <- readLines(x)
  
  find <- grep("Meridional mode:", a, fixed = T) %>% as.numeric()
  b <- read.csv("Current_AMM_Atl3.csv",header=T) %>% round(.,2)

  substr(a[ find -2 ], 7, 16) <- as.character(today)

  if (b[1,]>0){
    substr(a[ find ], 18, 22) <- format(b[1,], nsmall=2) %>% paste0(" ", .)
  } else { substr(a[ find ], 18, 22) <- format(b[1,], nsmall=2)}
  
  if (b[2,]>0){
    substr(a[ find +2 ], 18, 22) <- format(b[2,], nsmall=2) %>% paste0(" ", .)
  } else { substr(a[ find +2 ], 18, 22) <- format(b[2,], nsmall=2)}
  
  write.table(a, x, row.names = FALSE,col.names = FALSE,quote = FALSE)
}

ch.file.Atl("./part01/Atl_impacts.md")
