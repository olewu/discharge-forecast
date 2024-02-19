library(RHBV)
#-----------------------------------------------------------------------------------------------------------#
#Prepare data:
setwd("/projects/NS9873K/owul/projects/discharge_forecast/")

catchday=fread("data/historical_data/sildre_nve/catchday.csv")
catchday=catchday[,.(time,value,stat_id,stationName)]
catchprop=fread("results/catchment_properties/nve/catchprop.csv")
weather=fread("data/historical_data/senorge/nve/seNorge_daily242.csv")
setnames(weather,"catchid","stat_id")

catchprop=catchprop[,.(regine_area,main_no,point_no,station_name,longitude,latitude,height_minimum,height_hypso_10,height_hypso_20,height_hypso_30,height_hypso_40,height_hypso_50,
                       height_hypso_60,height_hypso_70,height_hypso_80,height_hypso_90,height_maximum,perc_agricul,perc_bog,perc_eff_bog,perc_eff_lake,
                       perc_forest,perc_glacier,perc_lake,perc_mountain,perc_urban, gradient_basin,gradient_river,area_total)]

catchprop[,stat_id:=paste0(regine_area,".",main_no,".0")]

catchday=merge(catchday,catchprop,"stat_id")
catchday[,date:=as_date(substr(time,1,10))]
alldat=merge(catchday,weather,c("date","stat_id"),all.x=FALSE)
numdata=alldat[,.N,stat_id]
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


#----------------------------------
library(foreach)
#https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
parallel::detectCores()
n.cores <- 8 #parallel::detectCores()-1

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)#PSOCK or FORK
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
#----------------------------------


unique_ids=unique(alldat[,stat_id])

maxyears=10
minyears=5

#length(unique_ids)
reslist=foreach(j = 1:100) %dopar%{
  this_stat_id=unique_ids[j]
  test_station=alldat[stat_id==this_stat_id]

  test_station=prepare_data(test_station)

  if(dim(test_station)[1]< 365*minyears){
    return(list())
  }

  if(dim(test_station)[1]> 365*maxyears){#max number of years,set to 10 here.
    test_station=test_station[(dim(test_station)[1]-365*maxyears+1):(dim(test_station)[1]),]
  }
  mindate=min(test_station$date)
  maxdate=max(test_station$date)

  subarea=test_station$area_total[1]/10
  from=c(test_station$height_minimum[1],test_station$height_hypso_10[1],test_station$height_hypso_20[1],test_station$height_hypso_30[1],
         test_station$height_hypso_40[1],test_station$height_hypso_50[1],test_station$height_hypso_60[1],test_station$height_hypso_70[1],
         test_station$height_hypso_80[1],test_station$height_hypso_90[1])
  to=c(test_station$height_hypso_10[1],test_station$height_hypso_20[1],test_station$height_hypso_30[1],
       test_station$height_hypso_40[1],test_station$height_hypso_50[1],test_station$height_hypso_60[1],test_station$height_hypso_70[1],
       test_station$height_hypso_80[1],test_station$height_hypso_90[1],test_station$height_maximum[1])

  if(any(is.na(c(from,to,subarea))==TRUE)){
    return(list())
  }

  basin=data.table(from=from,to=to,areal=rep(subarea,10))
  basin=make_basin_data(basin,refh=test_station$elev[1])


  to_mm_const=1/(sum(basin$`area(km2)`)*10^6)*60*60*24*10^3
  test_station[,qt:=value*to_mm_const]
  test_station=unique(test_station)

  if(any(is.na(c(test_station$st,test_station$prec))==TRUE)){
    test_station$st=na_interpolation(test_station$st,maxgap=5)
    test_station$prec=na_interpolation(test_station$prec,maxgap=5)
  }

  if(any(is.na(c(test_station$st,test_station$prec,test_station$qt,test_station$elev[1]))==TRUE)){
    return(list(j=j))
  }

  weatherlist=distribute_weather(test_station$st,test_station$prec,test_station$date,
                                 basin_data=basin,latitude=round(test_station$latitude[1],0))

  HBVpars=set_hbv_params()
  HBVfit=calibrate_HBV(weatherlist,test_station$qt,HBVpars,score_function="NS")

  optimparams=HBVfit$optim$bestmem
  reslist=data.table(calibparams=optimparams,paramnum=1:length(optimparams),stat_id=this_stat_id,mindate=mindate,maxdate=maxdate)

  return(reslist)
}


parallel::stopCluster(cl = my.cluster)

reslist=rbindlist(reslist)
fwrite(reslist,file="results/catchment_HBV_parameters/nve/catchparams_NSE_1_100.csv")


#How to simulate:
#sim_HBV(mypars,weatherlist,init_snow=0,init_soil=0,init_routing=c(0,0,0))
