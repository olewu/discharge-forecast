library(devtools)
#usethis::use_git_config(user.name = "tjroksva", user.email = "thea.roksvag@gmail.com")
#credentials::set_github_pat("senorgetok")
#remotes::install_github("NorskRegnesentral/seNorge",force=TRUE)
library(seNorge)
library(data.table)
library(geoR)
library(rgdal)
#------------------------------------------------------------#
catchprop=fread("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/catchprop.csv")
catchids=catchprop$stat_id


#------------------------------------------------------------#
setwd('/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Data/NVEdata/Feltgrenser/')
catchments <- readOGR(dsn = '.', layer = 'Felt_Avreninngkart')
#polygons<-spTransform(catchments,CRS("+proj=utm +zone=33 +datum=WGS84"))
polygons=spTransform(catchments, CRS("+proj=longlat +datum=WGS84"))
polygons$stID=polygons$stasjonNr #Make sure that the names of the polygons are stored in polygons$stID.
polygons$area_km2=polygons$Shape_Area/10^6 #Make sure that the areas of the polygons are stored in polygons$area_km2. #This is only important if we are going to use the plotAreas()-funksjon.


poly_ind=which(polygons$stID%in%catchids)
#-----------------------------------------------------------#
sN_locs=get_coord_key()
#----------------------------------------------------------------------------------------------#

catchweather=c()

for(i in 1:length(poly_ind)){
  print(i)
  currpoly=polygons@polygons[[poly_ind[i]]]@Polygons[[1]]@coords
  coords=sN_locs
  pointsinpoly=point.in.polygon(coords$lon, coords$lat, currpoly[,1], currpoly[,2],mode.checked=FALSE)
  theselocs=coords[pointsinpoly==1,]
  theselocs[,catchid:=polygons$stID[poly_ind[i]]]


  weatherdat=load_sN(years =1990:2022,
                     locs = theselocs[,.(lon,lat)],
                     vars = c('prec','st',"st_min",'st_max'),
                     dir = sN_dir(),
                     add_elevation = FALSE,
                     mc_cores = 8)


  weatherdat=merge(weatherdat,theselocs[,.(lon,lat,catchid,elev)],by=c("lon","lat"))

  catchweather=rbind(catchweather,weatherdat[,.("prec"=round(mean(prec),2),"st"=round(mean(st),2),"st_min"=round(mean(st_min),2),"st_max"=round(mean(st_max),2),"elev"=round(mean(elev,na.rm=TRUE),2)),.(date,catchid)])

  if(i%in%c(10,20,30,40,50,60,70,80,90,100,120,150,180,200,240)){
    fwrite(catchweather,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/seNorge_daily.csv")
    #fwrite(catchweather,"/nr/samba/user/roksvag/Prosjekter2023/seNorge_daily.csv")

  }
}

fwrite(catchweather,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/seNorge_daily242.csv")
#fwrite(catchweather,"/nr/samba/user/roksvag/Prosjekter2023/seNorge_daily242.csv")



#-----------------------------------------------------------------------------#

#Square version:
#lonbb=c(min(currpoly[,1]),max(currpoly[,1]))
#latbb=c(min(currpoly[,2]),max(currpoly[,2]))
#theselocs=sN_locs[lon<=lonbb[2] & lon>=lonbb[1] & lat>= latbb[1] & lat<=latbb[2],]
#theselocs[,catchid:=polygons$stID[poly_ind[i]]]

coordlist=list()
k=1
for(i in 1:length(poly_ind)){
  print(i)
  currpoly=polygons@polygons[[poly_ind[i]]]@Polygons[[1]]@coords
  coords=sN_locs
  pointsinpoly=point.in.polygon(coords$lon, coords$lat, currpoly[,1], currpoly[,2],mode.checked=FALSE)
  theselocs=coords[pointsinpoly==1,]
  theselocs[,catchid:=polygons$stID[poly_ind[i]]]

  coordlist[[k]]=theselocs
  k=k+1

}

coordlist=rbindlist(coordlist)

fwrite(coordlist,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/sN_allcoords_nve.csv")

library(data.table)
coordlist=fread("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/sN_allcoords_nve.csv")

xy <- SpatialPoints(cbind(Longitude=coordlist$lon,Latitude=coordlist$lat),proj4string = CRS("+proj=longlat"))
res=spTransform(xy, "+proj=utm +zone=33 +datum=WGS84")

coordlist[,utmx:=res@coords[,1]]
coordlist[,utmy:=res@coords[,2]]
fwrite(coordlist,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/sN_allcoords_nve.csv")

