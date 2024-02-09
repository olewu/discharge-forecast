library(RHBV)
setwd("/projects/NS9001K/owul/projects/discharge_forecast/")
#----------------------------------------------------------------------------------#
catchprop_nve=fread(file="results/catchment_properties/nve/catchprop_nveapi.csv")

type="NSE"
hbv_params=rbind(fread(file=paste0("results/catchment_HBV_parameters/nve/catchparams_",type,"_1_100.csv")),
                 fread(file=paste0("results/catchment_HBV_parameters/nve/catchparams_",type,"_101_242.csv")))

sN=unique(fread(file="data/historical_data/senorge/nve/seNorge_daily242.csv"))
sN_old=unique(fread("data/historical_data/senorge/nve/seNorge_daily_nve_1960-1989.csv"))

sN=rbind(sN_old,sN) #if we want to include 1960-1989.

vfdat_nve=fread(file="data/historical_data/sildre_nve/catchday.csv")

catchprop_nve=catchprop_nve[order(area_total),]
allcatchments=unique(catchprop_nve$stat_id)
#----------------------------------------------------------------------------------#
par(mfrow=c(3,3))
set.seed(2)
catchs=1:90

for(jj in catchs){

  currcatch=allcatchments[jj]
  catchprop_nve[stat_id==currcatch,area_total]

  n_donors=5
  properties_to_include=c("utm_east_z33","utm_north_z33","area_total","height_hypso_50","height_maximum","gradient_1085","length_km_river","mean_summer_prec","mean_fall_prec",
                                                    "perc_forest","perc_lake","perc_mountain","mean_summer_temp","mean_winter_temp","specific_runoff")
  #----------------------------------------------------------------------------------#

  our_properties=copy(catchprop_nve[stat_id%in%hbv_params$stat_id])

  donortab=find_my_donors(our_properties,my_id=currcatch,properties_to_include,to_scale=TRUE)
  if(is.null(donortab)==TRUE){next}

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
  target_weather=sN[catchid==currcatch]
  target_elev_bands=make_elevation_bands(target_area,hypsovec)
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


  to_mm_const=1/(sum(target_basin_data$`area(km2)`)*10^6)*60*60*24*10^3
  donormedian_mm=apply(donorsim_mm,2,median)
  donormedian_cumecs=apply(donorsim_mm*to_mm_const^{-1},2,median)

  vf_calib_mm=vf_calib$tf_module
  vf_calib_cumecs=vf_calib_mm*to_mm_const^{-1}

  res=data.table(stat_id=currcatch,date=weatherlist$datevec,q_sim_mm=donormedian_mm,q_sim_cumecs=donormedian_cumecs,q_calib_mm=vf_calib_mm,q_calib_cumecs=vf_calib_cumecs)
  res=res[year(date)>1990,]
  sN_check=sN[year(date)>1990 & catchid == currcatch]

  #-------------------------------------------------------------------------------------------------------------------#
  #Check specific runoff:
  #ind=1:365+365*5
  #plot(res$date[ind],res$q_sim_mm[ind],type="l",xlab="Date",ylab="Streamflow m^3/s",lwd=2,ylim=c(-20,100));title(currcatch);grid()
  #lines(sN_check$date[ind],sN_check$prec[ind],col=rgb(0,0,1,alpha=0.5))
  #lines(sN_check$date[ind],sN_check$st[ind],col=rgb(1,0,0,alpha=0.5))

  true_specific_runoff=fread(file="data/historical_data/sildre_nve/specific_runoff.csv")[stat_id==currcatch] #ok, beregnet fra data i den aktuelle perioden.
  true_monthly_runoff=fread(file="data/historical_data/sildre_nve/specific_runoff_month.csv")[stat_id==currcatch] #ok, beregnet fra data i den aktuelle perioden.


  res[,month:=month(date)]
  res[,year:=year(date)]

  res=res[year>=1990]#remove burn-in
  res[,annual_runoff_mm:=sum(q_sim_mm),.(year)]
  res[, monthly_runoff_mm:=sum(q_sim_mm),.(year,month)]
  res[, mean_monthly_runoff_mm:=mean(monthly_runoff_mm),.(month)]
  res[, mean_annual_runoff_mm:=mean(annual_runoff_mm),]

  res[,annual_runoff_calib_mm:=sum(q_calib_mm),.(year)]
  res[, monthly_runoff_calib_mm:=sum(q_calib_mm),.(year,month)]
  res[, mean_monthly_runoff_calib_mm:=mean(monthly_runoff_calib_mm),.(month)]
  res[, mean_annual_runoff_calib_mm:=mean(annual_runoff_calib_mm),]

  MAF_true_mm=true_specific_runoff$mean_annual
  MAF_pred_mm=unique(res$mean_annual_runoff_mm)
  MAF_calib_mm=unique(res$mean_annual_runoff_calib_mm)

  specific_true=true_specific_runoff$specific_runoff
  specific_pred=MAF_pred_mm*true_specific_runoff$to_cubic_const * true_specific_runoff$to_specific_const

  MAF_month_mm=unique(res[order(month),mean_monthly_runoff_mm])
  MAF_month_calib_mm=unique(res[order(month),mean_monthly_runoff_calib_mm])

  if(length(true_monthly_runoff$mean_monthly)==0){next}
  plot(1:12,MAF_month_mm,type="o",pch=21,bg="skyblue",cex=1.5,xlab="Month",ylim=quantile(c(MAF_month_mm,true_monthly_runoff$mean_monthly,MAF_month_calib_mm),c(0,1))+c(-5,5),ylab="mm/month",cex.axis=1.5,cex.lab=1.5);grid()
  lines(1:12,true_monthly_runoff$mean_monthly,type="o",cex=1.5,pch=21,bg="orange")
  lines(1:12,MAF_month_calib_mm,type="o",cex=1.5,pch=21,bg="yellow")

  abline(h=MAF_true_mm/12,col="orange",lty=2)
  abline(h=MAF_pred_mm/12,col="skyblue",lty=2)
  abline(h=MAF_calib_mm/12,col="yellow",lty=3)

  #legend("topleft",c("Predicted","True","Calib","Annual/12"),pt.bg=c("skyblue","orange","yellow","skyblue"),pch=c(21,21,21,NA),lty=c(NA,NA,NA,2),cex=2)
  title(paste0("ID ",currcatch,", lon ",round(target_longitude,1),", lat ",round(target_latitude,1)," area ",round(target_area,2)))
  # Legg også til hva som skjer med de sanne parameterne til feltet for å se
  # om man gjør det mye dårligere med donorfelt.
  #-------------------------------------------------------------------------------------------------------------------#

}

if(0){
  #catch jj=67
  wd="data/historical_data/sildre_nve/"
  catchday=fread(paste0(wd,"catchday.csv"))[,.(time,value,stat_id)]

  toplot=res
  catchday=catchday[stat_id==toplot$stat_id[1]]
  catchday[,date:=as_date(paste0(year(time),"-",month(time),"-",day(time)))]
  catchday=merge(catchday[,.(date,value)],res,by="date")
  catchday=catchday[year(date)==2020 & month(date)%in%c(1:6)]

  par(mar=c(5,5,5,5))
  plot(catchday$date, catchday$q_calib_cumecs,type="l",col="orange",lwd=2,ylim=c(0,20),ylab="m^3/s",xlab="Day",cex.axis=1.5,cex.lab=1.5);grid()
  lines(catchday$date,catchday$value,type="l",col="darkblue",lwd=2)
  #legend("topright",col=c("darkblue","orange"),c("Observed","Simulated"),cex=2,lty=1,lwd=2)

}
