library(data.table)
library(RHBV)
library(imputeTS)
setwd("/projects/NS9001K/owul/projects/discharge_forecast/")

#---------------------------------------------------------------------------------------------------#
#-----Upload data files-----#
smaakraft_prop=fread("results/catchment_properties/smaakraft/smaakraft_prop_and_clim.csv")
smaakraft_prop[,catch_category:="smaakraft"]
nve_prop=fread(file="results/catchment_properties/nve/nve_prop_and_clim.csv")
nve_prop[,catch_category:="nve"]
nve_prop=nve_prop[perc_glacier<=0.1]

sN=fread("data/historical_data/senorge/smaakraft/seNorge_daily_smaakraft.csv")
# sN_old=fread("RHBV/data/seNorge_daily_smaakraft_1960.csv")

# sN=rbind(sN_old,sN) #if we want to include 1960-1989.
sN=sN[order(catchname,date)]

hbv_params=rbind(fread(file="results/catchment_HBV_parameters/nve/catchparams_MSE_1_100.csv"),
                 fread(file="results/catchment_HBV_parameters/nve/catchparams_MSE_101_242.csv"))

startyear=1990 #1960 for the case where we use seNorge from 1960 ->. Otherwise 1990.
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
properties_to_include=c("utm_east_z33","utm_north_z33","area_total","height_hypso_50","height_maximum","gradient_1085","length_km_river","mean_summer_prec","mean_fall_prec",
                        "perc_forest","perc_lake","perc_mountain","mean_summer_temp","mean_winter_temp","specific_runoff")

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

  alldates=data.table(date=seq(as_date("1990-01-01"),as_date("2022-12-30"),by=1))
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

#-----------------------------------------------------------------------------------------------------------------------#

res[,month:=month(date)]
res[,year:=year(date)]
res[,day:=day(date)]
res=res[year>startyear] #remove "burnin" year.
res[,monthly_sum_mm:=sum(q_sim_mm),.(month,year,stat_id,utm_x,utm_y,lon,lat,area)]
res[,annual_sum_mm:=sum(q_sim_mm),.(year,stat_id,utm_x,utm_y,lon,lat,area)]

#------------------------------------------------------------------------------------------------------------------------#

monthly_means=res[,.("monthly_mean_cumecs"=mean(q_sim_cumecs),"monthly_mean_mm"=mean(monthly_sum_mm)),.(month,stat_id,utm_x,utm_y,lon,lat,area)]
annual_means=res[,.("annual_mean_cumecs"=mean(q_sim_cumecs),"annual_mean_mm"=mean(annual_sum_mm)),.(stat_id,utm_x,utm_y,lon,lat,area)]

#------------------------------------------------------------------------------------------------------------------------#
fwrite(annual_means,"results/discharge_historical_simulations/smaakraft/annual_means.csv")
fwrite(monthly_means,"results/discharge_historical_simulations/smaakraft/monthly_means.csv")
fwrite(res,"results/discharge_historical_simulations/smaakraft/daily_sim.csv")
fwrite(donorsave,"results/discharge_historical_simulations/smaakraft/donordata.csv")
#------------------------------------------------------------------------------------------------------------------------#


#---Lag noen Norges-kart------#
# library(sp);library(ggplot2);library(rgdal);library(viridis)

# setwd("/nr/project/stat/ClimateFutures/FoodProduction/Blending/Norgeomriss/")
# norge <- readOGR(dsn = ".", layer = "norge")
# norge=spTransform(norge, CRS("+proj=longlat +datum=WGS84"))

# #projnorge=proj4string(norge)
# #norge=spTransform(norge, CRS(projnorge))
# #norgecoords=norge@polygons[[1]]@Polygons[[1]]@coords
# #norgecoords=data.frame(x=norgecoords[,1]/10000,y=norgecoords[,2]/10000,id=1)
# themeinfo=theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),legend.title = element_text(size=20),legend.text = element_text(size=20),strip.text=element_text(size=20))

# pdf(file=paste0("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Results/historical_simulations/07072023/monthly_runoff.pdf"),width=14,height=8)
# ggplot()+
#   geom_polygon(aes(x=long, y=lat, group = id), data=norge, fill="white",color="gray",lwd=0.3)+
#   geom_point(data=monthly_means,
#              mapping=aes(x=lon,y=lat,color=monthly_mean_mm),pch=19,cex=2)+facet_wrap(~as.factor(month))+
#   xlab("Lon")+ylab("Lat")+theme_classic()+scale_color_gradientn(colours = rev(viridis(20)),name="Runoff (mm)")+themeinfo
# dev.off()
# #xlim(c(3,32))+ylim(c(55,72))


# pdf(file=paste0("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Results/historical_simulations/07072023/annual_runoff.pdf"),width=6,height=4)
# ggplot()+
#   geom_polygon(aes(x=long, y=lat, group = id), data=norge, fill="white",color="gray",lwd=0.3)+
#   geom_point(data=annual_means,
#              mapping=aes(x=lon,y=lat,color=annual_mean_mm),pch=19,cex=2)+
#   xlab("Lon")+ylab("Lat")+theme_classic()+scale_color_gradientn(colours = rev(viridis(20)),name="Runoff (mm)")+themeinfo
# dev.off()



# for(j in 1:length(allcatchments)){
#   png(file=paste0("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Results/historical_simulations/07072023/catchfig/cumecs/",allcatchments[j],"_monthly.png"),width=600,height=400)
#   plot(monthly_means[stat_id==allcatchments[j],.(month,monthly_mean_cumecs)],
#        type="o",xlab="Måned",ylab="Gjennomsnittlig vannføring (m^3/s)",pch=21,bg="skyblue",cex=2,
#        cex.lab=1.5,cex.axis=1.5,cex.main=1.5);
#   title(paste0("Catch ", allcatchments[j],", lon ",
#                round(donorsave[target_catchment==allcatchments[j],longitude][1],1),", lat ",
#                round(donorsave[target_catchment==allcatchments[j],latitude][1],1),", area ",
#                round(donorsave[target_catchment==allcatchments[j],area_total][1],1) ))
#   grid()
#   points(monthly_means[stat_id==allcatchments[j],.(month,monthly_mean_cumecs)],type="o",pch=21,bg="skyblue",cex=2)
#   dev.off()

# }


# for(j in 1:length(allcatchments)){
#   png(file=paste0("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Results/historical_simulations/07072023/catchfig/mm/",allcatchments[j],"_monthly.png"),width=600,height=400)
#   plot(monthly_means[stat_id==allcatchments[j],.(month,monthly_mean_mm)],
#        type="o",xlab="Måned",ylab="Vannføring (mm/month)",pch=21,bg="skyblue",cex=2,
#        cex.lab=1.5,cex.axis=1.5,cex.main=1.5,ylim=c(0,800));
#   title(paste0("Catch ", allcatchments[j],", lon ",
#                round(donorsave[target_catchment==allcatchments[j],longitude][1],1),", lat ",
#                round(donorsave[target_catchment==allcatchments[j],latitude][1],1),", area ",
#                round(donorsave[target_catchment==allcatchments[j],area_total][1],1) ))
#   grid()
#   points(monthly_means[stat_id==allcatchments[j],.(month,monthly_mean_mm)],type="o",pch=21,bg="skyblue",cex=2)
#   dev.off()

# }


# #---------------------------------------------------------------------------------------------------------------------------------#
# setwd("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/")
# #Compare with sweco, seN and true gWh
# history=fread("Data/Smaakraftdata/Reference_sim/history.csv",header=TRUE)

# history[,stat_id:=tolower(stat_id)]
# history[,stat_id:=gsub('ø','oe',stat_id)]
# history[,stat_id:=gsub('å','aa',stat_id)]
# history[,stat_id:=gsub('æ','ae',stat_id)]

# obs=matrix(unlist(t(history[,2:13])),ncol=1)
# place=rep(unique(history$stat_id),each=12)
# month=rep(1:12,length(unique(place)))
# history=data.table(stat_id=place,month=month,prod_history=as.vector(obs))
# compare=merge(monthly_means,history,c("stat_id","month"))


# sweco=fread("Data/Smaakraftdata/Reference_sim/sweco_sim.csv",header=TRUE)
# sweco[,stat_id:=tolower(stat_id)]
# sweco[,stat_id:=gsub('ø','oe',stat_id)]
# sweco[,stat_id:=gsub('å','aa',stat_id)]
# sweco[,stat_id:=gsub('æ','ae',stat_id)]

# obs=matrix(unlist(t(sweco[,2:13])),ncol=1)
# place=rep(unique(sweco$stat_id),each=12)
# month=rep(1:12,length(unique(place)))
# sweco=data.table(stat_id=place,month=month,sweco_history=as.vector(obs))
# compare=merge(compare,sweco,c("stat_id","month"))


# sN_original=fread("Data/Smaakraftdata/Reference_sim/seNorge_sim.csv",header=TRUE)
# sN_original[,stat_id:=tolower(stat_id)]
# sN_original[,stat_id:=gsub('ø','oe',stat_id)]
# sN_original[,stat_id:=gsub('å','aa',stat_id)]

# sN_original[,stat_id:=gsub('æ','ae',stat_id)]

# obs=matrix(unlist(t(sN_original[,2:13])),ncol=1)
# place=rep(unique(sN_original$stat_id),each=12)
# month=rep(1:12,length(unique(place)))
# sN_original=data.table(stat_id=place,month=month,sNoriginal_history=as.vector(obs))
# compare=merge(compare,sN_original,c("stat_id","month"))


# tocheck=unique(compare$stat_id)
# for(j in 1:length(tocheck)){
#   png(file=paste0("/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/Results/historical_simulations/07072023/catchfig/prod/",tocheck[j],"_prodcompare.png"),width=1200,height=1000)
#   par(mfrow=c(2,1))
#   plot(compare[stat_id==tocheck[j],.(month,monthly_mean_cumecs)],type="o",xlab="Month",ylab="Monthly mean (m^3/s)",pch=21,bg="skyblue",cex=1.5);grid();title(tocheck[j])
#   legend("topright",cex=1.5,pch=19,col=c("skyblue","yellow","red","darkgreen"),c("CF","Swe","sN","True"),horiz=TRUE)
#   ymin=min(compare[stat_id==tocheck[j],.(prod_history,sweco_history,sNoriginal_history)])
#   ymax=max(compare[stat_id==tocheck[j],.(prod_history,sweco_history,sNoriginal_history)])
#   plot(compare[stat_id==tocheck[j],.(month,sweco_history)],type="o",pch=21,bg="yellow",cex=1.5,ylim=c(ymin,ymax),ylab="Production");grid();title(tocheck[j])
#   points(compare[stat_id==tocheck[j],.(month,sNoriginal_history)],type="o",pch=21,bg="red",cex=1.5)
#   points(compare[stat_id==tocheck[j],.(month,prod_history)],type="o",pch=21,bg="darkgreen",cex=1.5)
#   dev.off()
# }

