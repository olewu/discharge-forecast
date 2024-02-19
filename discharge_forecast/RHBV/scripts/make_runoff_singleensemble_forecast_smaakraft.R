library(data.table)
library(RHBV)
library(imputeTS)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (ensemble number).", call.=FALSE)
} else if (length(args) > 1) {
  print("Only one argument required. Additional arguments will be ignored.")
} else {
  ensn = as.integer(args[1])
}

setwd("/projects/NS9873K/owul/projects/discharge_forecast/")

#-----------------------------------------------------------------------------#
#-----Upload data files-----#
smaakraft_prop=fread("results/catchment_properties/smaakraft/smaakraft_prop_and_clim.csv")
smaakraft_prop[,catch_category:="smaakraft"]
nve_prop=fread(file="results/catchment_properties/nve/nve_prop_and_clim.csv")
nve_prop[,catch_category:="nve"]
nve_prop=nve_prop[perc_glacier<=0.1]

DATE = Sys.Date() # as.Date("2024-02-04")
startyear=as.integer(format(DATE, "%Y")) - 3

FC=fread(paste("results/ens_forecast_input/smaakraft/fc_init_",DATE,".csv",sep=""))

FC=FC[order(catchname,date)]

# get senorge data
sn_path = 'data/regular_downloads/senorge/smaakraft'
file_list <- list.files(path = sn_path, pattern = "seNorge_.*\\.csv")

sN_merge <- data.frame()
for (file in file_list) {
  # Extract the date from the file name using regular expressions
  date_str <- gsub("^[^0-9]+\\_|\\.[A-Za-z]+$", "", file)
  file_date <- as.Date(date_str, format = "%Y%m%d")

  if (file_date <= DATE) {
    data <- read.csv(paste(sn_path,file,sep='/'), header = TRUE)
    data <- data[,c('catchname', 'date', 'prec', 'st', 'elev')]
    sN_merge <- rbind(sN_merge, data)
  }
}

# get archive senorge data from startyear to 2022:
sN_archive = 'data/historical_data/senorge/smaakraft/seNorge_daily_smaakraft.csv'
sN_arch = read.csv(sN_archive, header = TRUE)
sN_arch = sN_arch[,c('catchname', 'date', 'prec', 'st', 'elev')]
sN_arch_sel <- subset(sN_arch, date >= as.Date(sprintf("%s-01-01",startyear)))

sn_full <- unique(rbind(sN_arch_sel,sN_merge))
sn_full = sn_full[order(sn_full$catchname,sn_full$date),]

hbv_params=rbind(fread(file="results/catchment_HBV_parameters/nve/catchparams_MSE_1_100.csv"),
                 fread(file="results/catchment_HBV_parameters/nve/catchparams_MSE_101_242.csv"))


#-----------------------------------------------------------------------------#

all_prop=rbind(smaakraft_prop,nve_prop,fill=TRUE)
save_var=c("stat_id","catch_category","utm_east_z33","utm_north_z33","longitude","latitude","area_total","height_minimum","height_hypso_10","height_hypso_20","height_hypso_30","height_hypso_40",
           "height_hypso_50","height_hypso_60","height_hypso_70","height_hypso_80","height_hypso_90","height_maximum","gradient_1085","length_km_river","mean_summer_prec","mean_fall_prec",
           "perc_forest","perc_lake","perc_mountain","mean_summer_temp","mean_winter_temp","specific_runoff")
all_prop=all_prop[,.SD,.SDcols=save_var]
#-----------------------------------------------------------------------------#

set.seed(2)

#Donor selection properties:
#-----------------------------------#

n_donors=5
properties_to_include=c("utm_east_z33","utm_north_z33","area_total","height_hypso_50","height_maximum","gradient_1085","length_km_river","mean_summer_prec","mean_fall_prec",
                        "perc_forest","perc_lake","perc_mountain","mean_summer_temp","mean_winter_temp","specific_runoff")

#------------------------------------#

allcatchments=unique(all_prop[catch_category=="smaakraft",stat_id])
donorsave=list(); k=1

print(ensn)

# select one ensemble member:
pcol = sprintf('prec_%s',ensn-1)
tcol = sprintf('st_%s',ensn-1)
columnms = c("catchname","date",pcol,tcol)
fc_subs <- FC[, .SD, .SDcols = columnms] # -1 since indexing in python begins on 0!
names(fc_subs)[names(fc_subs) %in% c(pcol,tcol)] = c("prec","st")

# merge the ensemble member with seNorge:
sn_tf = copy(sn_full) # [, -which(names(sn_full) == "elev")]
# names(sn_tf)[names(sn_tf) %in% c("prec","st")] = c(pcol,tcol)
sn_tf$date = as.Date(sn_tf$date, format = "%Y-%m-%d")
# fc_sn_merge = merge(sn_tf, fc_subs, by = 'catchname', all = TRUE)
fc_subs$elev = NA
fc_subs$date = as.Date(fc_subs$date, format = "%Y-%m-%d")
fc_sn_merge = as.data.table(rbind(sn_tf,fc_subs))

res_ens = list()
new_colnames = c(sprintf("q_sim_mm_%s",ensn-1),sprintf("q_sim_cumecs_%s",ensn-1))

for(j in 1:length(allcatchments)){
  print(j)
  currcatch=allcatchments[j] #should be a smaakraft catchment.

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
  target_weather=fc_sn_merge[catchname==currcatch]
  target_elev_bands=make_elevation_bands(target_area,hypsovec)
  target_basin_data=make_basin_data(target_elev_bands,refh=target_weather$elev[1])

  target_weather[,N:=.N,date]
  target_weather=target_weather[N==1,]

  alldates=data.table(date=seq(as.Date(paste(startyear,"-01-01",sep='')),as_date(DATE+25),by=1))
  target_weather=merge(alldates,target_weather,all.x=TRUE)
  target_weather[,catchname:=currcatch]


  if(any(is.na(c(target_weather$st,target_weather$prec))==TRUE)){
    target_weather$st=na_interpolation(target_weather$st,maxgap=5)
    target_weather$prec=na_interpolation(target_weather$prec,maxgap=5)
  }

  weatherlist=distribute_weather(target_weather$st,target_weather$prec,datevec=target_weather$date,
                                basin_data=target_basin_data,latitude=round(target_latitude,2))
  weatherlist$datevec=as_date(target_weather$date)

  # Run HBV model 'n_donor' times with the catchment weather and the parameters
  # of the 'n_donor' most similar catchments
  donorsim_mm=c()
  for(i in 1:n_donors){
    curr_params=hbv_params[stat_id %in% donor_id[i]]
    vf_donor=sim_HBV(mypars=curr_params$calibparams,weatherlist)
    donorsim_mm=rbind(donorsim_mm,vf_donor$tf_module)
  }
  to_mm_const=1/(target_area*10^6)*60*60*24*10^3
  donormedian_mm=apply(donorsim_mm,2,median)
  donormedian_cumecs=apply(donorsim_mm*to_mm_const^{-1},2,median)

  new = data.table(
    stat_id=currcatch,
    date=weatherlist$datevec,
    q_sim_mm=donormedian_mm,
    q_sim_cumecs=donormedian_cumecs,
    utm_x=utm_x,utm_y=utm_y,
    lon=target_longitude,lat=target_latitude,
    area=target_area
  )

  setnames(new,c("q_sim_mm","q_sim_cumecs"),new_colnames)

  res_ens = rbind(res_ens, new)

  #save donor information.
  if(ensn == 1){
    donorsave_tmp=donortab[1:6]
    donorsave_tmp[,target_catchment:=donortab[1,stat_id]]
    donorsave[[k]]=donorsave_tmp
    k=k+1
  }
}

if(ensn == 1){
  donorsave=rbindlist(donorsave)
  fwrite(donorsave,paste("results/discharge_forecast/smaakraft/donordata_",DATE,".csv",sep=""))
}

res_ens[,month:=month(date)]
res_ens[,year:=year(date)]
res_ens[,day:=day(date)]
res_ens=res_ens[-(1:365)] # remove first year (spin-up)

savefc_directory <- paste("results/discharge_forecast/smaakraft/daily21d_",DATE,sep="")

# Check if the directory already exists
if(!file.exists(savefc_directory)){
  # If it doesn't exist, create the directory
  dir.create(savefc_directory, recursive = TRUE)
  cat("Directory created:", savefc_directory, "\n")
}

ENSN = sprintf("%02d", ensn)
fwrite(res_ens,paste(savefc_directory,"/daily21d_",DATE,"_ens",ENSN,".csv",sep=""))

print(DATE)