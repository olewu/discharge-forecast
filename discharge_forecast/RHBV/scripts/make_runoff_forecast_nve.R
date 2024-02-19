library(RHBV)
setwd("/projects/NS9873K/owul/projects/discharge_forecast/")
#----------------------------------------------------------------------------------#
catchprop_nve=fread(file="results/catchment_properties/nve/catchprop_nveapi.csv")

type="NSE"
hbv_params=rbind(fread(file=paste0("results/catchment_HBV_parameters/nve/catchparams_",type,"_1_100.csv")),
                 fread(file=paste0("results/catchment_HBV_parameters/nve/catchparams_",type,"_101_242.csv")))

DATE = "2024-02-18" #Sys.Date() # "2024-01-13"

# startyear=as.integer(format(Sys.Date(), "%Y")) - 3

sN=unique(fread(file=paste("results/forecast_input/nve/fc_init_",DATE,"T06:00:00Z_merge_sn.csv",sep=""))) #replace by smaakraft weather data.
vfdat_nve=fread(file="data/historical_data/sildre_nve/catchday.csv")

catchprop_nve=catchprop_nve[order(area_total),]
allcatchments=unique(catchprop_nve$stat_id)
#----------------------------------------------------------------------------------#
set.seed(2)
catchs=1:90

res=list();donorsave=list();k=1
res_vf=list()
for(jj in catchs){

  currcatch=allcatchments[jj] #shuold be a smaakraft catchment.
  catchprop_nve[stat_id==currcatch,area_total]

  print(jj)
 
  n_donors=5
  properties_to_include=c("utm_east_z33","utm_north_z33","area_total","height_hypso_50","height_maximum","gradient_1085","length_km_river","mean_summer_prec","mean_fall_prec",
                                                    "perc_forest","perc_lake","perc_mountain","mean_summer_temp","mean_winter_temp","specific_runoff")
  #----------------------------------------------------------------------------------#

  our_properties=copy(catchprop_nve[stat_id%in%hbv_params$stat_id])
  
  target_area=our_properties[stat_id==currcatch]$area_total

  target_latitude=our_properties[stat_id==currcatch]$latitude
  target_longitude=our_properties[stat_id==currcatch]$longitude

  utm_x=our_properties[stat_id==currcatch]$utm_east_z33
  utm_y=our_properties[stat_id==currcatch]$utm_north_z33

  donortab=find_my_donors(our_properties,my_id=currcatch,properties_to_include,to_scale=TRUE)
  
  if(is.null(donortab)==TRUE){next}
  
  print(currcatch)

  mydonors=donortab[2:(n_donors+1),]
  donor_id=unique(mydonors$stat_id)

  #------------------------------------------------------------#
  #Prepare target weather data (make sure it is in "our properties")
  target_latitude=our_properties[stat_id==currcatch]$latitude
  target_area=our_properties[stat_id==currcatch]$area_total
  target_longitude=our_properties[stat_id==currcatch]$longitude

  hypsovec=unlist(our_properties[stat_id==currcatch,.(height_minimum,height_hypso_10,height_hypso_20,height_hypso_30,
                                                      height_hypso_40,height_hypso_50,height_hypso_60,height_hypso_70,
                                                      height_hypso_80,height_hypso_90,height_maximum)])



  if(length(hypsovec)==0){next}
  target_weather=sN[catchname==currcatch]
  # target_elev_bands=make_elevation_bands(target_area,hypsovec)
  target_elev_bands=data.table(from=hypsovec[1:(length(hypsovec)-1)],to=hypsovec[2:(length(hypsovec))],areal=target_area/(length(hypsovec)-1))
  target_basin_data=make_basin_data(target_elev_bands,refh=target_weather$elev[1])

  if(any(is.na(c(target_weather$st,target_weather$prec))==TRUE)){
    target_weather$st=na_interpolation(target_weather$st,maxgap=5)
    target_weather$prec=na_interpolation(target_weather$prec,maxgap=5)
  }


  weatherlist=distribute_weather(target_weather$st,target_weather$prec,datevec=target_weather$date,
                                 basin_data=target_basin_data,latitude=round(target_latitude,2))
  weatherlist$datevec=as_date(target_weather$date)
  #-------------------------------------------------------------#
  true_params=hbv_params[stat_id ==currcatch]
  vf_calib=sim_HBV(mypars=true_params$calibparams,weatherlist)

  donorsim_mm=c()
  for(j in 1:n_donors){
    curr_params=hbv_params[stat_id %in% donor_id[j]]
    vf_donor=sim_HBV(mypars=curr_params$calibparams,weatherlist)
    donorsim_mm=rbind(donorsim_mm,vf_donor$tf_module)
  }

  to_mm_const=1/(target_area*10^6)*60*60*24*10^3
  donormedian_mm=apply(donorsim_mm,2,median)
  donormedian_cumecs=apply(donorsim_mm*to_mm_const^{-1},2,median)

  vf_calib_cumecs=vf_calib$tf_module*to_mm_const^{-1}

  res=rbind(res,data.table(stat_id=currcatch,date=weatherlist$datevec,q_sim_mm=donormedian_mm,q_sim_cumecs=donormedian_cumecs,utm_x=utm_x,utm_y=utm_y,lon=target_longitude,lat=target_latitude,area=target_area))

  res_vf=rbind(res_vf,data.table(stat_id=currcatch,date=weatherlist$datevec,q_sim_mm=vf_calib$tf_module,q_sim_cumecs=vf_calib_cumecs,utm_x=utm_x,utm_y=utm_y,lon=target_longitude,lat=target_latitude,area=target_area))

  donorsave_tmp=donortab[1:6]
  donorsave_tmp[,target_catchment:=donortab[1,stat_id]]
  donorsave[[k]]=donorsave_tmp
  k=k+1

}

donorsave=rbindlist(donorsave)

res[,month:=month(date)]
res[,year:=year(date)]
res[,day:=day(date)]
res=res[-(1:365)] #remove "burnin" year.

res_vf[,month:=month(date)]
res_vf[,year:=year(date)]
res_vf[,day:=day(date)]
res_vf=res_vf[-(1:365)] #remove "burnin" year. year>startyear

print('saving')

#-----------------------------------------------------------------------------#
fwrite(res,paste("results/discharge_forecast/nve/daily_",DATE,"T06:00:00Z.csv",sep=""))
fwrite(res_vf,paste("results/discharge_forecast/nve/daily_local_",DATE,"T06:00:00Z.csv",sep=""))
fwrite(donorsave,paste("results/discharge_forecast/nve/donordata_",DATE,"T06:00:00Z.csv",sep=""))
#-----------------------------------------------------------------------------#

