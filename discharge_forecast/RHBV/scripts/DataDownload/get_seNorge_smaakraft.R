library(data.table)
library(geoR)
#library(devtools)
library(rgdal)
library(proj4)
library(seNorge)

#install.packages("rgdal")

#To run:
#source("/nr/project/stat/RenewableEnergy/Smaakraft/HBVpackage/RHBV/scripts/DataDownload/get_seNorge_smaakraft.R")

#------------------------------------------------------------#
#usethis::use_git_config(user.name = "tjroksva", user.email = "thea.roksvag@gmail.com")
#usethis::create_github_token()
#credentials::set_github_pat("senorgetok2")
#remotes::install_github("NorskRegnesentral/seNorge",force=TRUE)

shapelist=list.files("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Data/Smaakraftdata/Feltgrenser/")

sN_locs=get_coord_key()
coords=sN_locs


catchweather=c()
lonlat=c()
k=1
for(j in shapelist){
  print(j)
  setwd(paste0('/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Data/Smaakraftdata/Feltgrenser/',j,'/'))
  catchments <- readOGR(dsn = '.', layer = 'NedbfeltF_v3')

  polygons=spTransform(catchments, CRS("+proj=longlat +datum=WGS84"))
  polygons$stID=polygons$vassdragNr #Make sure that the names of the polygons are stored in polygons$stID.
  polygons$area_km2=polygons$areal_km2 #Make sure that the areas of the polygons are stored in polygons$area_km2. #This is only important if we are going to use the plotAreas()-funksjon.
  polygons$name=j


  currpoly=polygons@polygons[[1]]@Polygons[[1]]@coords
  pointsinpoly=point.in.polygon(coords$lon, coords$lat, currpoly[,1], currpoly[,2],mode.checked=FALSE)


  theselocs=coords[pointsinpoly==1,]

  if(dim(theselocs)[1]==0){
    theselocs=coords[which.min(sqrt((coords$lon-currpoly[1,1])^2+(coords$lat-currpoly[1,2])^2)),]
    print("This was a tiny catchment.")
  }

  print(dim(theselocs))

  cleanname=polygons$name
  cleanname=tolower(cleanname)
  cleanname=gsub('ø','oe',cleanname)
  cleanname=gsub('æ','ae',cleanname)
  cleanname=gsub('å','aa',cleanname)
  cleanname=gsub(' ','_',cleanname)

  theselocs[,catchid:=polygons$vassdragNr]
  theselocs[,catchname:=cleanname]

  weatherdat=load_sN(years =1960:1989,
                     locs = theselocs[,.(lon,lat)],
                     vars = c('prec','st',"st_min",'st_max'),
                     dir = sN_dir(),
                     add_elevation = FALSE,
                     mc_cores = 4)

  weatherdat=merge(weatherdat,theselocs[,.(lon,lat,catchname,elev)],by=c("lon","lat"))

  catchweather=rbind(catchweather,
                     weatherdat[,.("prec"=round(mean(prec),2),
                                   "st"=round(mean(st),2),
                                   "st_min"=round(mean(st_min),2),
                                   "st_max"=round(mean(st_max),2),
                                   "elev"=round(mean(elev,na.rm=TRUE),2)),.(date,catchname)])

  currlonlat=unique(weatherdat[,.(lon,lat)])
  currlonlat= currlonlat[,.("lon"=round(mean(lon),2),
                            "lat"=round(mean(lat),2))]

  lonlat=rbind(lonlat,currlonlat)



  if(k %in% c(1,10,20,30,40,50,60,70,80,90,100,110,120,130)){
    fwrite(catchweather,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/seNorge_daily_smaakraft_1960.csv")
    fwrite(lonlat,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/seNorge_lonlat_smaakraft_1960.csv")
  }
  k=k+1

}

fwrite(catchweather,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/seNorge_daily_smaakraft_1960.csv")
fwrite(lonlat,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/seNorge_lonlat_smaakraft_1960.csv")


#cacthweather=fread("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/seNorge_daily_smaakraft.csv")



#---------------------Download sN coordinates used for averaging-----------------------------------------#
shapelist=list.files("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Data/Smaakraftdata/Feltgrenser/")
coordlist=list()
k=1
for(j in shapelist){
  print(j)
  setwd(paste0('/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Data/Smaakraftdata/Feltgrenser/',j,'/'))
  catchments <- readOGR(dsn = '.', layer = 'NedbfeltF_v3')

  polygons=spTransform(catchments, CRS("+proj=longlat +datum=WGS84"))
  polygons$stID=polygons$vassdragNr #Make sure that the names of the polygons are stored in polygons$stID.
  polygons$area_km2=polygons$areal_km2 #Make sure that the areas of the polygons are stored in polygons$area_km2. #This is only important if we are going to use the plotAreas()-funksjon.
  polygons$name=j

  cleanname=polygons$name
  cleanname=tolower(cleanname)
  cleanname=gsub('ø','oe',cleanname)
  cleanname=gsub('æ','ae',cleanname)
  cleanname=gsub('å','aa',cleanname)
  cleanname=gsub(' ','_',cleanname)


  currpoly=polygons@polygons[[1]]@Polygons[[1]]@coords
  pointsinpoly=point.in.polygon(coords$lon, coords$lat, currpoly[,1], currpoly[,2],mode.checked=FALSE)


  theselocs=coords[pointsinpoly==1,]

  if(dim(theselocs)[1]==0){
    theselocs=coords[which.min(sqrt((coords$lon-currpoly[1,1])^2+(coords$lat-currpoly[1,2])^2)),]
    print("This was a tiny catchment.")
  }

  theselocs$stat_id=cleanname
  coordlist[[k]]=theselocs
  k=k+1

}
coordlist=rbindlist(coordlist)

fwrite(coordlist,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/sN_allcoords_smaakraft.csv")

coordlist=fread("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/sN_allcoords_smaakraft.csv")


xy <- SpatialPoints(cbind(Longitude=coordlist$lon,Latitude=coordlist$lat),proj4string = CRS("+proj=longlat"))
res=spTransform(xy, "+proj=utm +zone=33 +datum=WGS84")

coordlist[,utmx:=res@coords[,1]]
coordlist[,utmy:=res@coords[,2]]
fwrite(coordlist,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/sN_allcoords_smaakraft.csv")

#---------------------------------------------------------------------------------------------------------------------#

