library(data.table)
library(httr)
library(rgdal)
library(jsonlite)

testit <- function(x)
{
  p1 <- proc.time()
  Sys.sleep(x)
  proc.time() - p1 # The cpu usage should be negligible
}

API_KEY = ''

setwd("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/")
catchprop=fread("data/catchprop.csv")
uniquecatch=unique(catchprop$stat_id)
names(catchprop)

catchprop_api=list()
for(j in seq_along(uniquecatch)){
  print(j)
  myurl<-paste0('https://hydapi.nve.no/api/v1/Stations?StationId=',uniquecatch[j],'&Parameter=Stations&Active=0')
  httpResponse <- GET(myurl, add_headers('X-API-Key' = API_KEY, 'Accept' = 'application/json'))
  results = fromJSON(content(httpResponse, "text"))

  catchprop_api[[j]]=data.table(stat_id=results$data$stationId,
                       latitude=results$data$latitude,
                       longitude=results$data$longitude,
                       utm_east_z33=results$data$utmEast_Z33,
                       utm_north_z33=results$data$utmNorth_Z33,
                       station_name=results$data$stationName,
                       drainage_basin_key=results$data$drainageBasinKey,
                       regine_area=results$data$regineNo,
                       gradient_1085=results$data$gradient1085,
                       gradient_basin=results$data$gradientBasin,
                       gradient_river=results$data$gradientRiver,
                       length_km_basin=results$data$lengthKmBasin,
                       length_km_river=results$data$lengthKmRiver,
                       perc_agricul=results$data$percentAgricul,
                       perc_bog=results$data$percentBog,
                       perc_eff_bog=results$data$percentEffBog,
                       perc_eff_lake=results$data$percentEffLake,
                       perc_lake=results$data$percentLake,
                       perc_glacier=results$data$percentGlacier,
                       perc_forest=results$data$percentForest,
                       perc_mountain=results$data$percentMountain,
                       perc_urban=results$data$percentUrban,
                       area_total=results$data$drainageBasinArea,
                       area_norway=results$data$drainageBasinAreaNorway,
                       height_minimum=results$data$heightMinimum,
                       height_maximum=results$data$heightMaximum,
                       height_hypso_10=results$data$heightHypso10,
                       height_hypso_20=results$data$heightHypso20,
                       height_hypso_30=results$data$heightHypso30,
                       height_hypso_40=results$data$heightHypso40,
                       height_hypso_50=results$data$heightHypso50,
                       height_hypso_60=results$data$heightHypso60,
                       height_hypso_70=results$data$heightHypso70,
                       height_hypso_80=results$data$heightHypso80,
                       height_hypso_90=results$data$heightHypso90,
                       specific_runoff=results$data$specificDischarge,
                       annual_runoff=results$data$annualRunoff)

  testit(1)
}

catchprop_api=rbindlist(catchprop_api)
fwrite(catchprop_api,"data/catchprop_nveapi.csv")
