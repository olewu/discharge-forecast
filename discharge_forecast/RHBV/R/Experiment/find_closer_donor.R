library(data.table)
library(RHBV)
library(imputeTS)
setwd("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/")

#---------------------------------------------------------------------------------------------------#
#-----Upload data files-----#
smaakraft_prop=fread("Data/Smaakraftdata/smaakraft_prop_and_clim.csv")
smaakraft_prop[,catch_category:="smaakraft"]
nve_prop=fread(file="HBVpackage/RHBV/data/nve_prop_and_clim.csv")
nve_prop[,catch_category:="nve"]
nve_prop=nve_prop[perc_glacier<=0.1]

sN=fread("HBVpackage/RHBV/data/seNorge_daily_smaakraft.csv")
sN_old=fread("HBVpackage/RHBV/data/seNorge_daily_smaakraft_1960.csv")

sN=rbind(sN_old,sN) #if we want to include 1960-1989.
sN=sN[order(catchname,date)]

hbv_params=rbind(fread(file="HBVpackage/RHBV/data/catchparams_MSE_1_100.csv"),
                 fread(file="HBVpackage/RHBV/data/catchparams_MSE_101_242.csv"))

startyear=1960 #1960 for the case where we use seNorge from 1960 ->. Otherwise 1990.
#--------------------------------------------------------------------------------------------------#

all_prop=rbind(smaakraft_prop,nve_prop,fill=TRUE)
save_var=c("stat_id","catch_category","utm_east_z33","utm_north_z33","longitude","latitude","area_total","height_minimum","height_hypso_10","height_hypso_20","height_hypso_30","height_hypso_40",
           "height_hypso_50","height_hypso_60","height_hypso_70","height_hypso_80","height_hypso_90","height_maximum","gradient_1085","length_km_river","mean_summer_prec","mean_fall_prec",
           "perc_forest","perc_lake","perc_mountain","mean_summer_temp","mean_winter_temp","specific_runoff")
all_prop=all_prop[,.SD,.SDcols=save_var]
#--------------------------------------------------------------------------------------------------#



par(mfrow=c(3,3))
set.seed(2)

#Donor selection properties:
#-----------------------------------#

n_donors=5
properties_to_include=c("utm_east_z33","utm_north_z33","area_total","height_maximum","mean_summer_prec","mean_fall_prec","specific_runoff")

#------------------------------------#

allcatchments=unique(all_prop[catch_category=="smaakraft",stat_id])
res=list();donorsave=list();k=1
for(j in 1:length(allcatchments)){
  print(j)
  currcatch=allcatchments[j] #shuold be a smaakraft catchment.

  #Catch properties for nve catchments and target smaakraft catchment.
  our_properties=copy(all_prop[stat_id%in%unique(c(hbv_params$stat_id,currcatch))])
  target_area=our_properties[stat_id==currcatch]$area_total

  target_latitude=our_properties[stat_id==currcatch]$latitude
  target_longitude=our_properties[stat_id==currcatch]$longitude

  utm_x=our_properties[stat_id==currcatch]$utm_east_z33
  utm_y=our_properties[stat_id==currcatch]$utm_north_z33

  donortab=find_my_donors(our_properties,my_id=currcatch,properties_to_include,to_scale=TRUE)
  mydonors=donortab[2:(n_donors+1),]
  donor_id=unique(mydonors$stat_id)


  hypsovec=unlist(our_properties[stat_id==currcatch,.(height_minimum,height_hypso_10,height_hypso_20,height_hypso_30,
                                                      height_hypso_40,height_hypso_50,height_hypso_60,height_hypso_70,
                                                      height_hypso_80,height_hypso_90,height_maximum)])
  target_weather=sN[catchname==currcatch]
  target_elev_bands=make_elevation_bands(target_area,hypsovec)
  target_basin_data=make_basin_data(target_elev_bands,refh=target_weather$elev[1])

  target_weather[,N:=.N,date]
  target_weather=target_weather[N==1,]

  alldates=data.table(date=seq(as_date("1960-01-01"),as_date("2022-12-30"),by=1))
  target_weather=merge(alldates,target_weather,all.x=TRUE)
  target_weather[,catchname:=currcatch]



  if(any(is.na(c(target_weather$st,target_weather$prec))==TRUE)){
    target_weather$st=na_interpolation(target_weather$st,maxgap=5)
    target_weather$prec=na_interpolation(target_weather$prec,maxgap=5)
  }

  weatherlist=distribute_weather(target_weather$st,target_weather$prec,datevec=target_weather$date,
                                 basin_data=target_basin_data,latitude=round(target_latitude,2))
  weatherlist$datevec=as_date(target_weather$date)

  donorsim_mm=c()
  for(i in 1:n_donors){
    curr_params=hbv_params[stat_id %in% donor_id[i]]
    vf_donor=sim_HBV(mypars=curr_params$calibparams,weatherlist)
    donorsim_mm=rbind(donorsim_mm,vf_donor$tf_module)
  }
  to_mm_const=1/(sum(target_basin_data$`area(km2)`)*10^6)*60*60*24*10^3
  donormedian_mm=apply(donorsim_mm,2,median)
  donormedian_cumecs=apply(donorsim_mm*to_mm_const^{-1},2,median)

  res=rbind(res,data.table(stat_id=currcatch,date=weatherlist$datevec,q_sim_mm=donormedian_mm,q_sim_cumecs=donormedian_cumecs,utm_x=utm_x,utm_y=utm_y,lon=target_longitude,lat=target_latitude,area=target_area))

  donorsave_tmp=donortab[1:6]
  donorsave_tmp[,target_catchment:=donortab[1,stat_id]]
  donorsave[[k]]=donorsave_tmp
  k=k+1

  #save donor information.
}

donorsave=rbindlist(donorsave)
fwrite(donorsave,"/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Results/historical_simulations/experimental_vannmerker/donordata_vannmerker_experimental.csv")
