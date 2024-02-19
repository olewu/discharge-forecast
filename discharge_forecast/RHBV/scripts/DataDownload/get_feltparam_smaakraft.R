library(data.table)
library(imputeTS)
library(sp)
# library(rgdal)

setwd("/projects/NS9873K/owul/projects/discharge_forecast/")

LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')) ## for example
  res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
  return(as.data.frame(res))
}#


catchprop=fread(file="RHBV/data/catchprop_nveapi.csv")
#hbv_params=rbind(fread(file="HBVpackage/RHBV/data/catchparams_MSE_1_100.csv"),
                 #fread(file="HBVpackage/RHBV/data/catchparams_MSE_101_242.csv"))



file="pdf_catchment_props/feltparam_smaakraft.csv"
feltparam_sk=fread(file)
#-----feltparam_smaakraft is made from Ole's skript, catchment_data_from_pdf.
#It is both in RHBV and in "Smaakraftdata/catchment_data".

#----------Data cleaning--------------------------#
feltparam_sk=transpose(feltparam_sk, fill=NA, ignore.empty=FALSE, keep.names=NULL, make.names=1)

names(catchprop)
names(feltparam_sk)

#Make sure "feltparam_sk" has the same column names as "catchprop" for the variables we are using:
setnames(feltparam_sk,c("V1","Vassdragsnr.","Beregn.punkt","Fylke.","Areal (A)","HoydeMIN","Hoyde10","Hoyde20","Hoyde30","Hoyde40","Hoyde50","Hoyde60","Hoyde70","Hoyde80","Hoyde90","HoydeMAX","Elvegradent1085 (EG,1085)",
                        "Elvleengde (EL)","Skog (ASKOG)","Sjo (ASJO)","Snaufjell (ASF)","Avrenning 1961-90 (QN)"),
         c("stat_id","catch_id","position","fylke","area_total","height_minimum","height_hypso_10","height_hypso_20","height_hypso_30","height_hypso_40","height_hypso_50","height_hypso_60","height_hypso_70",
           "height_hypso_80","height_hypso_90","height_maximum","gradient_1085","length_km_river","perc_forest","perc_lake","perc_mountain","specific_runoff"))

#Take out the most important variables:
feltparam=feltparam_sk[,.(stat_id,catch_id,position,fylke,area_total,height_minimum,height_hypso_10,height_hypso_20,height_hypso_30,height_hypso_40,height_hypso_50,
                height_hypso_60,height_hypso_70,height_hypso_80,height_hypso_90,height_maximum,gradient_1085,length_km_river,perc_forest,perc_lake,perc_mountain,specific_runoff)]

#Fix position information from "position" string:
feltparam[,is_east:=unlist(gregexpr('E', position))!=-1]
feltparam[is_east==TRUE,utm_east_z33:=as.numeric(substr(position,1,unlist(gregexpr('E', position))-1))]
feltparam[is_east==FALSE,utm_east_z33:=-1*as.numeric(substr(position,1,unlist(gregexpr('W', position))-1))]
feltparam[is_east==TRUE,utm_north_z33:=as.numeric(substr(position,unlist(gregexpr('E', position))+1,unlist(gregexpr('N', position))-1))]
feltparam[is_east==FALSE,utm_north_z33:=as.numeric(substr(position,unlist(gregexpr('W', position))+1,unlist(gregexpr('N', position))-1))]

#Plot the coordinates with the NVE catch coordinates:
# plot(feltparam[,.(utm_east_z33,utm_north_z33)],col="green")
# points(catchprop[,.(utm_east_z33,utm_north_z33)],col="red")
#the projection looks reasonable!

#remove nordic letters and only use small letters:
feltparam[,stat_id:=gsub('ø','oe',stat_id)]
feltparam[,stat_id:=gsub('å','aa',stat_id)]
feltparam[,stat_id:=gsub('æ','ae',stat_id)]
feltparam[,stat_id:=gsub('_i','_1',stat_id)]
feltparam[,stat_id:=gsub('_ii','_2',stat_id)]

#also remove spaces and replace by "_":
feltparam[,stat_id:=gsub(' ','_',stat_id)]

#---------Fix missing height hypso-----------#
uniquecatch=unique(feltparam$stat_id)

feltparam2=list()
for(j in 1:length(uniquecatch)){
  curr=feltparam[stat_id==uniquecatch[j]]
  if(any(curr=="-m")==FALSE){
    feltparam2=rbind(feltparam2,curr)
    next
  }

  subcurr=copy(curr[,.(height_minimum,height_hypso_10,height_hypso_20,height_hypso_30,height_hypso_40,height_hypso_50,
                  height_hypso_60,height_hypso_70,height_hypso_80,height_hypso_90,height_maximum)])
  hypso=as.vector(unlist(subcurr))

  ind=which(hypso=="-m")
  hypso[ind]=NA
  hypso=as.numeric(hypso)
  hypso=na_interpolation(hypso)

  toreplace=colnames(subcurr)[ind]
  toreplace_ind=which(colnames(curr)%in%toreplace)

  curr=data.frame(curr)
  curr[,toreplace]=hypso[ind]

  feltparam2=rbind(feltparam2,curr)
}

feltparam=data.table(feltparam2)


#-----Need to add information about climatology for the Smaakraft catchments------#
sN=fread("RHBV/data/seNorge_daily_smaakraft.csv")
sN_lonlat=fread("RHBV/data/seNorge_lonlat_smaakraft.csv")

sN[,month:=month(date)]
sN[,year:=year(date)]

sN[,annual_prec:=sum(prec),.(year,catchname)]
sN[,wet_days:=sum(prec>=1)/365,.(year,catchname)]

sN[month %in%c(6,7,8),summer_prec:=sum(prec),.(year,catchname)]
sN[month %in%c(9,10,11),fall_prec:=sum(prec),.(year,catchname)]


sN[month %in%c(6,7,8),summer_temp:=mean(st,na.rm=TRUE),.(year,catchname)]
sN[month %in%c(12,1,2),winter_temp:=mean(st,na.rm=TRUE),.(year,catchname)]


sN[,mean_annual_prec:=mean(prec),.(catchname)]
sN[,mean_summer_prec:=mean(summer_prec,na.rm=TRUE),.(catchname)]
sN[,mean_fall_prec:=mean(fall_prec,na.rm=TRUE),.(catchname)]
sN[,mean_wet_days:=mean(wet_days),.(catchname)]

sN[,mean_summer_temp:=mean(summer_temp,na.rm=TRUE),.(catchname)]
sN[,mean_winter_temp:=mean(winter_temp,na.rm=TRUE),.(catchname)]

sN_summary=unique(sN[,.(mean_summer_prec,mean_fall_prec,mean_annual_prec,mean_wet_days,mean_summer_temp,mean_winter_temp,catchname)])
setnames(sN_summary,"catchname","stat_id")

propandclim=merge(sN_summary,feltparam,by="stat_id")
#------------------------------------------------------------------#

#From utm to lon/lat:
utm=propandclim[,.(utm_east_z33,utm_north_z33)]

LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')) ## for example
  res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
  return(as.data.frame(res))
}#


lonlat=LongLatToUTM(utm$utm_east_z33,utm$utm_north_z33,33)

propandclim[,longitude:=lonlat$coords.x1]
propandclim[,latitude:=lonlat$coords.x2]


fwrite(propandclim,"/projects/NS9873K/owul/smaakraft/data_for_ole/smaakraft_catchments/smaakraft_prop_and_clim.csv")

