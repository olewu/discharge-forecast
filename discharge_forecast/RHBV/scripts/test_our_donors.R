library(RHBV)

#use donor_search.R to find donors with different strategies.
#In this script we test which strategy is better.
#Donorlist.csv is made from donor_search.R. See this script for strategies.

setwd("/projects/NS9001K/owul/projects/discharge_forecast/")
catchprop=fread(file="results/catchment_properties/nve/catchprop_nveapi.csv")

#change from MSE to NSE:
type="NSE"
hbv_params=rbind(fread(file=paste0("results/catchment_HBV_parameters/nve/catchparams_",type,"_1_100.csv")),
                 fread(file=paste0("results/catchment_HBV_parameters/nve/catchparams_",type,"_101_242.csv")))


donorlist=fread("results/donors/donorlist.csv")
sN=unique(fread(file="data/historical_data/senorge/nve/seNorge_daily242.csv"))

vfdat=fread("data/historical_data/sildre_nve/catchday.csv")


#------------------------------------------------------------------------------------------------------------------------#
target_ids=unique(donorlist$target_id)

index=1
MSE_all=list()
cor_all=list()
NSq_all=list()

for(curr_strategy in 1:6){
    print(paste0("Neste strategi er nummer: ",curr_strategy))
  for(j in 1:length(target_ids)){
    print(j)
    curr_target=target_ids[j]
    if(dim(hbv_params[stat_id==curr_target])[1]==0){
      next
    }

    curr_donors=donorlist[target_id==curr_target & strategy==curr_strategy,donor_id]
    true_params=hbv_params[stat_id==curr_target]

    test_station=vfdat[stat_id==curr_target]
    test_station[,date:=time]
    vf_true=data.table(qt=test_station$value,date=test_station$date)
    vf_true$date=substr(vf_true$date,1,10)
    vf_true$date=as_date(vf_true$date)
    vf_true=vf_true[!(is.na(vf_true$qt)==TRUE),]



    alldonors=c()
    for(k in 1:length(curr_donors)){
      curr_params=hbv_params[stat_id== (curr_donors[k])]
      target_latitude=catchprop[stat_id==curr_target]$latitude
      target_area=catchprop[stat_id==curr_target]$area_total
      hypsovec=unlist(catchprop[stat_id==curr_target,.(height_minimum,height_hypso_10,height_hypso_20,height_hypso_30,
                                                    height_hypso_40,height_hypso_50,height_hypso_60,height_hypso_70,
                                                    height_hypso_80,height_hypso_90,height_maximum)])

      target_weather=sN[catchid==curr_target]
      target_elev_bands=make_elevation_bands(target_area,hypsovec)
      target_basin_data=make_basin_data(target_elev_bands,refh=target_weather$elev[1])


      if(any(is.na(c(target_weather$st,target_weather$prec))==TRUE)){
        target_weather$st=na_interpolation(target_weather$st,maxgap=5)
        target_weather$prec=na_interpolation(target_weather$prec,maxgap=5)
      }

      weatherlist=distribute_weather(target_weather$st,target_weather$prec,datevec=target_weather$date,
                         basin_data=target_basin_data,latitude=round(target_latitude,2))
      weatherlist$datevec=as_date(target_weather$date)


      indtokeep=which(weatherlist$datevec%in%vf_true$date)
      datestokeep=weatherlist$datevec[indtokeep]
      datestokeep=vf_true$date[which(vf_true$date %in% datestokeep)]


      #Compare simulations with donor parameters with the simulations from true calib parameters#
      vf_donor=sim_HBV(mypars=curr_params$calibparams,weatherlist)
      vf_target=sim_HBV(mypars=true_params$calibparams,weatherlist)

      #Only evaluate for dates where we have streamflow obs.
      to_mm_const=1/(sum(target_basin_data$`area(km2)`)*10^6)*60*60*24*10^3

      compare_dates=vf_true$date[which(vf_true$date %in% datestokeep)]
      qt_obs=vf_true$qt[which(vf_true$date %in% compare_dates)]
      qt_donor=(vf_donor$tf_module[which(weatherlist$datevec%in%compare_dates)])*to_mm_const^{-1}
      qt_target=vf_target$tf_module[which(weatherlist$datevec%in%compare_dates)]*to_mm_const^{-1}

      alldonors=rbind(alldonors,qt_donor)
    }

    if(any(is.na(qt_obs)==TRUE)){
      qt_obs=na_interpolation(qt_obs,maxgap=5)
    }

    qt_donor_mean=colMeans(alldonors)
    qt_donor_median=apply(alldonors,2,median)

    alldonors=rbind(alldonors,qt_donor_mean,qt_donor_median)

    #-----How close the donor is to the streamflow simulated by the "true" params----#
    #MSE_2sim=apply(alldonors,1,MSEq,qobs=qt_target)
    #cor_2sim=apply(alldonors,1,cor,qt_target)

    #-----How close the donor hbv params model is to the actual observed streamflow----#
    MSE_obs_simdonor=apply(alldonors,1,MSEq,qobs=qt_obs)
    NSq_obs_simdonor=apply(alldonors,1,NSq,qobs=qt_obs)
    cor_obs_simdonor=apply(alldonors,1,cor,qt_obs)


    #-----How close the target hbv params is to the actual observed streamflow----#
    MSE_obs_simtarget=MSEq(qt_obs,qt_target)
    NSq_obs_simtarget=NSq(qt_obs,qt_target)
    cor_obs_simtarget=cor(qt_obs,qt_target)
    #------------------------------------------------------------------------------------------#

    MSE_res=c(MSE_obs_simdonor,MSE_obs_simtarget)
    cor_res=c(cor_obs_simdonor,cor_obs_simtarget)
    NSq_res=c(NSq_obs_simdonor,NSq_obs_simtarget)


    names(MSE_res)=c("donor1","donor2","donor3","donor4","donor5","donor_mean","donor_median","target_reference")
    names(cor_res)=c("donor1","donor2","donor3","donor4","donor5","donor_mean","donor_median","target_reference")
    names(NSq_res)=c("donor1","donor2","donor3","donor4","donor5","donor_mean","donor_median","target_reference")

    MSE_res_dt=data.table(t(MSE_res))
    cor_res_dt=data.table(t(cor_res))
    NSq_res_dt=data.table(t(NSq_res))
    MSE_res_dt[,stat_id:=curr_target]
    cor_res_dt[,stat_id:=curr_target]
    NSq_res_dt[,stat_id:=curr_target]
    MSE_res_dt[,strategy:=curr_strategy]
    cor_res_dt[,strategy:=curr_strategy]
    NSq_res_dt[,strategy:=curr_strategy]

    MSE_all[[index]]=MSE_res_dt
    cor_all[[index]]=cor_res_dt
    NSq_all[[index]]=NSq_res_dt

    index=index+1
  }
}

MSE_all=rbindlist(MSE_all)
cor_all=rbindlist(cor_all)
NSq_all=rbindlist(NSq_all)

#Change from /MSE/ to /NSE/
#fwrite(MSE_all,paste0("results/donors/",type,"/MSE_donors.csv"))
#fwrite(cor_all,paste0("results/donors/",type,"/cor_donors.csv"))
#fwrite(NSq_all,paste0("results/donors/",type,"/NSq_donors.csv"))

#------------------------------------------------------------------------------------------------------------------------------------#
#Change from /MSE/ to /NSE/
MSE_all=fread(paste0("results/donors/",type,"/MSE_donors.csv"))
cor_all=fread(paste0("results/donors/",type,"/cor_donors.csv"))
NSq_all=fread(paste0("results/donors/",type,"/NSq_donors.csv"))


par(mfrow=c(2,3),cex.axis=1.5,cex.lab=1.5,mar=c(5,5,5,5));#par(pty="s")
plot(y=cor_all[strategy==1,donor_median],cor_all[strategy==1,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 1")
plot(y=cor_all[strategy==2,donor_median],cor_all[strategy==2,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 2")
plot(y=cor_all[strategy==3,donor_median],cor_all[strategy==3,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 3")
plot(y=cor_all[strategy==4,donor_median],cor_all[strategy==4,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 4")
plot(y=cor_all[strategy==5,donor_median],cor_all[strategy==5,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 5")
plot(y=cor_all[strategy==6,donor_median],cor_all[strategy==6,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 6")


par(mfrow=c(2,3),cex.axis=1.5,cex.lab=1.5,mar=c(5,5,5,5));#par(pty="s")
plot(y=cor_all[strategy==1,donor1],cor_all[strategy==1,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 1")
plot(y=cor_all[strategy==2,donor1],cor_all[strategy==2,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 2")
plot(y=cor_all[strategy==3,donor1],cor_all[strategy==3,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 3")
plot(y=cor_all[strategy==4,donor1],cor_all[strategy==4,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 4")
plot(y=cor_all[strategy==5,donor1],cor_all[strategy==5,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 5")
plot(y=cor_all[strategy==6,donor1],cor_all[strategy==6,target_reference],cex=1.5,xlim=c(0,1),ylim=c(0,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 6")


par(mfrow=c(2,3),cex.axis=1.5,cex.lab=1.5,mar=c(5,5,5,5));#par(pty="s")
plot(y=NSq_all[strategy==1,donor_median],NSq_all[strategy==1,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 1")
plot(y=NSq_all[strategy==2,donor_median],NSq_all[strategy==2,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 2")
plot(y=NSq_all[strategy==3,donor_median],NSq_all[strategy==3,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 3")
plot(y=NSq_all[strategy==4,donor_median],NSq_all[strategy==4,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 4")
plot(y=NSq_all[strategy==5,donor_median],NSq_all[strategy==5,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 5")
plot(y=NSq_all[strategy==6,donor_median],NSq_all[strategy==6,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),col=rgb(red=1,green=0.5,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 6")


par(mfrow=c(2,3),cex.axis=1.5,cex.lab=1.5,mar=c(5,5,5,5));#par(pty="s")
plot(y=NSq_all[strategy==1,donor1],NSq_all[strategy==1,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 1")
plot(y=NSq_all[strategy==2,donor1],NSq_all[strategy==2,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 2")
plot(y=NSq_all[strategy==3,donor1],NSq_all[strategy==3,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 3")
plot(y=NSq_all[strategy==4,donor1],NSq_all[strategy==4,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 4")
plot(y=NSq_all[strategy==5,donor1],NSq_all[strategy==5,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 5")
plot(y=NSq_all[strategy==6,donor1],NSq_all[strategy==6,target_reference],cex=1.5,xlim=c(-0.5,1),ylim=c(-0.5,1),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest dono, strategy 6")

par(mfrow=c(2,3),cex.axis=1.5,cex.lab=1.5,mar=c(5,5,5,5));#par(pty="s")
plot(y=MSE_all[strategy==1,donor_median],MSE_all[strategy==1,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 1")
plot(y=MSE_all[strategy==2,donor_median],MSE_all[strategy==2,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 2")
plot(y=MSE_all[strategy==3,donor_median],MSE_all[strategy==3,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 3")
plot(y=MSE_all[strategy==4,donor_median],MSE_all[strategy==4,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 4")
plot(y=MSE_all[strategy==5,donor_median],MSE_all[strategy==5,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),xlab="True target params",ylab="Donor params",pch=19,col=rgb(red=1,green=0.5,blue=0,alpha=0.5));lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 5")
plot(y=MSE_all[strategy==6,donor_median],MSE_all[strategy==6,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),col=rgb(red=1,green=0.5,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Donor median, strategy 6")


par(mfrow=c(2,3),cex.axis=1.5,cex.lab=1.5,mar=c(5,5,5,5));#par(pty="s")
plot(y=MSE_all[strategy==1,donor1],MSE_all[strategy==1,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 1")
plot(y=MSE_all[strategy==2,donor1],MSE_all[strategy==2,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 2")
plot(y=MSE_all[strategy==3,donor1],MSE_all[strategy==3,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 3")
plot(y=MSE_all[strategy==4,donor1],MSE_all[strategy==4,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 4")
plot(y=MSE_all[strategy==5,donor1],MSE_all[strategy==5,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest donor, strategy 5")
plot(y=MSE_all[strategy==6,donor1],MSE_all[strategy==6,target_reference],cex=1.5,xlim=c(0,3),ylim=c(0,3),col=rgb(red=0.5,green=0.8,blue=0,alpha=0.5),pch=19,xlab="True target params",ylab="Donor params");grid();lines(c(-100,100),c(-100,100));grid(); title("Closest dono, strategy 6")


#------------------------------------------------------------------------------------------------------------------------------------#
library(ggplot2)
themeinfo=theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),legend.title = element_text(size=20),legend.text = element_text(size=20),plot.title = element_text(size = 20))+ theme(legend.position = "none")

ggplot(cor_all,aes(x=as.factor(strategy),y=donor_median,col=as.factor(strategy)))+geom_boxplot()+ylim(c(0.4,1))+xlab("Strategy")+ylab("Correlation")+themeinfo
ggplot(cor_all,aes(x=as.factor(strategy),y=donor1,col=as.factor(strategy)))+geom_boxplot()+ylim(c(0.4,1))+xlab("Strategy")+ylab("Correlation")+themeinfo
ggplot(cor_all,aes(x=as.factor(strategy),y=target_reference,col=as.factor(strategy)))+geom_boxplot()+ylim(c(0.4,1))+xlab("Strategy")+ylab("Correlation")+themeinfo

ggplot(NSq_all,aes(x=as.factor(strategy),y=donor_median,col=as.factor(strategy)))+geom_boxplot()+ylim(c(-0.5,1))+xlab("Strategy")+ylab("NSE")+themeinfo
ggplot(NSq_all,aes(x=as.factor(strategy),y=donor1,col=as.factor(strategy)))+geom_boxplot()+ylim(c(-0.5,1))+xlab("Strategy")+ylab("NSE")+themeinfo
ggplot(NSq_all,aes(x=as.factor(strategy),y=target_reference,col=as.factor(strategy)))+geom_boxplot()+ylim(c(-0.5,1))+xlab("Strategy")+ylab("NSE")+themeinfo

ggplot(MSE_all,aes(x=as.factor(strategy),y=donor_median,col=as.factor(strategy)))+geom_boxplot()+xlab("Strategy")+ylab("MSE")+themeinfo+ylim(c(0,3))
ggplot(MSE_all,aes(x=as.factor(strategy),y=donor1,col=as.factor(strategy)))+geom_boxplot()+xlab("Strategy")+ylab("MSE")+themeinfo+ylim(c(0,3))
ggplot(MSE_all,aes(x=as.factor(strategy),y=target_reference,col=as.factor(strategy)))+geom_boxplot()+xlab("Strategy")+ylab("MSE")+themeinfo+ylim(c(0,3))


#Sjekk hvor gode resultater vi får på måneds-snitt.
