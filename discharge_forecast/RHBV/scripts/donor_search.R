#Script for finding suitable donor catchments.
library(RHBV)
library(data.table)

setwd("/projects/NS9873K/owul/projects/discharge_forecast/")

catchprop_nve=fread(file="results/catchment_properties/nve/catchprop_nveapi.csv")
catchday=fread(file="data/historical_data/sildre_nve/catchday.csv") #1961-1990

hbv_params=rbind(fread(file="results/catchment_HBV_parameters/nve/catchparams_MSE_1_100.csv"),
                 fread(file="results/catchment_HBV_parameters/nve/catchparams_MSE_101_242.csv"))

#catchprop=catchprop[perc_glacier<0.1]

hbv_params=hbv_params[stat_id%in%catchprop$stat_id]
possible_donors=unique(hbv_params$stat_id)



#------------------------------------------------------------------------------#

catchday=merge(catchday,catchprop[,.(area_total,stat_id)],by="stat_id")
catchday=catchday[is.na(value)==FALSE,]
catchday[,to_mm_const:=1/(area_total*10^6)*60*60*24*10^3]
catchday[,value:=value*to_mm_const]

catchday[,month:=month(time)]
catchday[,year:=year(time)]

catchday[,N_days:=.N,.(year,stat_id)]
catchday=catchday[N_days>=350]
catchday[,annual_runoff:=sum(value),.(year,stat_id)]
catchyear=unique(catchday[,.(annual_runoff,year,stat_id,area_total)])
catchyear[,N_years:=.N,.(stat_id)]
catchyear=catchyear[N_years>=10,]

catchmeanrunoff=catchyear[,.("mean_annual"=mean(annual_runoff),"area_total"=area_total),.(stat_id)]
catchmeanrunoff=unique(catchmeanrunoff)
catchmeanrunoff[,to_cubic_const:=((area_total*10^6/1000)/(60*60*24*365))]
catchmeanrunoff[,cumecs:=mean_annual*to_cubic_const]
catchmeanrunoff[,to_specific_const:=1000/(area_total)]
catchmeanrunoff[,specific_runoff:=cumecs*to_specific_const]

catchday[,N_days:=.N,.(year,month,stat_id)]
catchday=catchday[N_days>=25]
catchmonth=catchday[,.("monthly_runoff"=sum(value)),.(year,month,stat_id,area_total)]
catchmonth[,N_years:=.N,.(month,stat_id)]
catchmonth=catchmonth[N_years>=10,]
catchmonth=unique(catchmonth[,.("mean_monthly"=mean(monthly_runoff),"area_total"=area_total),.(stat_id,month)])

#we dont use this for finding donors. However, we can use it for evaluations:
#fwrite(catchmeanrunoff,file="HBVpackage/RHBV/data/specific_runoff.csv")
#fwrite(catchmonth,file="HBVpackage/RHBV/data/specific_runoff_month.csv")

######catchprop=merge(catchprop,catchmeanrunoff[,.(stat_id,specific_runoff)],by="stat_id")


#1 l, 0.001 m^3
# 1000 l, 1 m^3

#-------------------------------------------------------------------------------#
sN=fread(file="data/historical_data/senorge/nve/seNorge_daily242.csv")
sN[,month:=month(date)]
sN[,year:=year(date)]

sN[,annual_prec:=sum(prec),.(year,catchid)]
sN[,wet_days:=sum(prec>=1)/365,.(year,catchid)]

sN[month %in%c(6,7,8),summer_prec:=sum(prec),.(year,catchid)]
sN[month %in%c(9,10,11),fall_prec:=sum(prec),.(year,catchid)]


sN[month %in%c(6,7,8),summer_temp:=mean(st,na.rm=TRUE),.(year,catchid)]
sN[month %in%c(12,1,2),winter_temp:=mean(st,na.rm=TRUE),.(year,catchid)]


sN[,mean_annual_prec:=mean(prec),.(catchid)]
sN[,mean_summer_prec:=mean(summer_prec,na.rm=TRUE),.(catchid)]
sN[,mean_fall_prec:=mean(fall_prec,na.rm=TRUE),.(catchid)]
sN[,mean_wet_days:=mean(wet_days),.(catchid)]

sN[,mean_summer_temp:=mean(summer_temp,na.rm=TRUE),.(catchid)]
sN[,mean_winter_temp:=mean(winter_temp,na.rm=TRUE),.(catchid)]

sN_summary=unique(sN[,.(mean_summer_prec,mean_fall_prec,mean_annual_prec,mean_wet_days,mean_summer_temp,mean_winter_temp,catchid)])
setnames(sN_summary,"catchid","stat_id")
#------------------------------------------------------------------------------#

propandclim=merge(sN_summary,catchprop,by="stat_id")
propandclim=propandclim[stat_id%in%possible_donors,]

#The smallest catchments available
propandclim=propandclim[order(area_total)]
target_catchments=propandclim[1:80,stat_id]

fwrite(propandclim,"results/catchment_properties/nve/nve_prop_and_clim.csv")

#------------------------------------------------------------------------------#
n_donors=5;
k=1

our_properties=propandclim
donorlist=list()
for(strategy in 1:6){
  if(strategy==1){
    properties_to_include=c("utm_east_z33","utm_north_z33")
    nametag="utm"
  }

  if(strategy==2){
    properties_to_include=c("utm_east_z33","utm_north_z33","area_total","specific_runoff")
    nametag="utm_area_specific"
  }

  if(strategy==3){
    properties_to_include=c("utm_east_z33","utm_north_z33","area_total","height_hypso_50","height_maximum","gradient_1085","length_km_river","mean_summer_prec","mean_fall_prec",
                            "perc_forest","perc_lake","perc_mountain","mean_summer_temp","mean_winter_temp","specific_runoff")
    nametag="utm_area_specific_hypso_perc_clim"
  }

  if(strategy==4){
    properties_to_include=c("utm_east_z33","utm_north_z33","specific_runoff")
    nametag="utm_specific"
  }

  if(strategy==5){
    properties_to_include=c("utm_east_z33","utm_north_z33","specific_runoff","area_total","height_hypso_50","height_maximum","gradient_1085",
                            "length_km_river","mean_summer_prec","mean_fall_prec",
                            "mean_summer_temp","mean_winter_temp","specific_runoff")
    nametag="utm_area_specific_clim"
  }

  if(strategy==6){
    properties_to_include=c("utm_east_z33","utm_north_z33","specific_runoff","area_total")
    nametag="utm_area_specific"
  }


  for(targetnum in 1:length(target_catchments)){
    currcatch=target_catchments[targetnum]

    donortab=find_my_donors(our_properties,my_id=currcatch,properties_to_include,to_scale=TRUE)

    mydonors=donortab[1:(n_donors+1),]

    donorlist[[k]]=data.table(target_id=currcatch,donor_id= mydonors$stat_id[2:(n_donors+1)],distance=mydonors$distance[2:(n_donors+1)],nametag=nametag,strategy=strategy)
    k=k+1

  }
}

donorlist=rbindlist(donorlist)
fwrite(donorlist,"results/donors/donorlist.csv")

#--------To test----------------------------------------
properties_to_include=c("utm_east_z33","utm_north_z33","specific_runoff","area_total","height_hypso_50","height_maximum","mean_summer_prec","mean_fall_prec",
                        "mean_summer_temp","mean_winter_temp","specific_runoff", "perc_forest","perc_lake","perc_mountain")
donortab=find_my_donors(our_properties,my_id="19.96.0",properties_to_include=properties_to_include,to_scale=TRUE)

plot(our_properties[,.(utm_east_z33,utm_north_z33)])
points(our_properties[stat_id %in%target_catchments,.(utm_east_z33,utm_north_z33)],col="green")
points(donortab[1:1,.(utm_east_z33,utm_north_z33)],col="red",pch=19)
points(donortab[2:6,.(utm_east_z33,utm_north_z33)],col="blue",pch=3)
#-----------------------------------------

