#Here is documentation of hydapi:
# https://hydapi.nve.no/UserDocumentation/#gettingstarted
#Here you can build your own url.
#https://hydapi.nve.no/swagger/index.html?urls.primaryName=V1
library(data.table)
library(httr)
library(rgdal)
library(jsonlite)
API_KEY = ''
#---------------------------------------------------------------------------------------------------------------------------------------------------------#

catchnames=substr(names(read.table("/nr/samba/user/roksvag/Nyttig/nve_api/runoffdata_hydrological_test.txt")),2,100) #from old nve project.


#------------------------------------------------------------------------------------------
catch_prop=fread("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Data/catchment_properties_thea.csv") #from old nve project "metadataMappe": H:\NVE_August2020\RawData\Data_2020
catch_prop[,stat_id:=paste0(catch_prop$regine,".",catch_prop$main,".0")]

#catch_prop=catch_prop[is.na(first_year_regulation)==TRUE & regulation_part_reservoirs==0,] #figure out what to remove here.
catch_prop=catch_prop[stat_id%in%catchnames,]
catch_prop=catch_prop[is.na(first_year_regulation)==TRUE & regulation_part_reservoirs<=0.2,]


#------------------------------------------------------------------------------------------
#Trenger kanskje ikke disse. Nok med utm east og north? Bruke kun 1 seNorge-punkt + høyde?
setwd('/nr/samba/user/roksvag/Nyttig/nve_api/Feltgrenser/')
catchments <- readOGR(dsn = '.', layer = 'Felt_Avreninngkart')
polygons<-spTransform(catchments,CRS("+proj=utm +zone=33 +datum=WGS84"))
polygons$stID=polygons$stasjonNr #Make sure that the names of the polygons are stored in polygons$stID.
polygons$area_km2=polygons$Shape_Area/10^6 #Make sure that the areas of the polygons are stored in polygons$area_km2. #This is only important if we are going to use the plotAreas()-funksjon.
#------------------------------------------------------------------------------------------

savelist=list();k=1
statID_d=unique(catch_prop$stat_id)

for(stationID in statID_d){
  start_time="1980-01-01"
  end_time="2022-11-20"

  print(stationID)

  # URL for reading daily data (regular)
  myurl<-paste0('https://hydapi.nve.no/api/v1/Observations?StationId=',stationID,'&Parameter=1001&ResolutionTime=day&VersionNumber=1&ReferenceTime=',start_time,'%2F',end_time,'&QualityTypes=3&TimeOffset=PT0H')

  httpResponse <- GET(myurl, add_headers('X-API-Key' = API_KEY, 'Accept' = 'application/json'))

  results = fromJSON(content(httpResponse, "text"))
  if(is.null(results$data$observations)==FALSE){
    if(length(results$data$observations[[1]])!=0){
      dates<-as.POSIXct(strptime(sub('T',' ',substr(as.character(results$data$observations[[1]][,1]),1,19)),"%Y-%m-%d %H:%M:%S"))

      data<-results$data$observations[[1]][,2]
      tosave=data.table(results$data$observations[[1]])
      tosave[,stat_id:=stationID]
      tosave[,stationName:=results$data$stationName]

      savelist[[k]]=tosave
      k=k+1
    }
  }
}

savelist=rbindlist(savelist)
fwrite(savelist,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Data/NVEdata/catchday.csv")

catch_prop=catch_prop[catch_prop$stat_id%in%savelist$stat_id,]
catch_prop[,N:=.N,stat_id]
catch_prop=catch_prop[N==1,]
fwrite(catch_prop,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Data/NVEdata/catchprop.csv")


savelist=savelist[savelist$stat_id%in%catch_prop$stat_id,]
fwrite(savelist,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Data/NVEdata/catchday.csv")

#We have data from all these catchments:
plot(catch_prop[,.(utm_east_z33,utm_north_z33)])


#-------#-------#-------#-------#-------#-------#-------#-------#-------
plot_id=which(polygons$stID %in%catch_prop$stat_id)

plot(NA,xlim=c(-1*10^5,11*10^5),ylim=c(6400000,8000000));
#axis(1);axis(2)
for(j in plot_id){
  lines(polygons@polygons[[j]]@Polygons[[1]]@coords,col="black")
  this_stat=polygons$stID[j]
  points(catch_prop[stat_id==this_stat,.(utm_east_z33,utm_north_z33)])

}

polygondata=data.table(poly_area=polygons$area_km2[plot_id],stat_id=polygons$stasjonNr[plot_id],poly_name=polygons$stasjonNav[plot_id])

catch_prop_subset=catch_prop[,.(area_total,stat_id,station_name,utm_east_z33,utm_north_z33,longitude,latitude)]
catch_prop_subset[,N:=.N,stat_id]
catch_prop_subset=catch_prop_subset[N==1,]

polygondata_check=merge(polygondata,catch_prop_subset,by="stat_id")

#Stemmer på en prikk: Polygonene ser ut til å samsvare med "catch info".
plot(polygondata_check[,.(poly_area,area_total)]);lines(c(0,15000),c(0,15000),col="red",lwd=2)


#-----------Sjekke at navna er de samme i catch_prop og i fila som er lasta ned fra api-----------------#
savelist2=merge(savelist,catch_prop_subset,by="stat_id")
sammenligning=unique(savelist2[,.(stationName,station_name)])
#ser likt ut!#
#----------------------------------------------------------------------------------------#


if(0){
  # URL for reading instataneous data (irregular)
  myurl<-'https://hydapi.nve.no/api/v1/Observations?StationId=12.171.0&Parameter=1001&ResolutionTime=inst&VersionNumber=1&ReferenceTime=2010-02-02%2F2015-03-03&QualityTypes=2&TimeOffset=PT0H'
  httpResponse <- GET(myurl, add_headers('X-API-Key' = API_KEY, 'Accept' = 'application/json'))

  results = fromJSON(content(httpResponse, "text"))
  dates<-as.POSIXct(strptime(sub('T',' ',substr(as.character(results$data$observations[[1]][,1]),1,19)),"%Y-%m-%d %H:%M:%S"))
  data<-results$data$observations[[1]][,2]
  plot(dates,data,type='l')

  # URL for reading hourly data (regular)
  myurl<-'https://hydapi.nve.no/api/v1/Observations?StationId=12.171.0&Parameter=1001&ResolutionTime=hour&VersionNumber=1&ReferenceTime=2010-02-02%2F2012-03-03&QualityTypes=2&TimeOffset=PT0H'
  httpResponse <- GET(myurl, add_headers('X-API-Key' = API_KEY, 'Accept' = 'application/json'))
  results = fromJSON(content(httpResponse, "text"))
  dates<-as.POSIXct(strptime(sub('T',' ',substr(as.character(results$data$observations[[1]][,1]),1,19)),"%Y-%m-%d %H:%M:%S"))
  data<-results$data$observations[[1]][,2]
  plot(dates,data,type='l')
}
