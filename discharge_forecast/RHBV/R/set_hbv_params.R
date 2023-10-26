set_hbv_params=function(){
  #assuming daily data
  #Snow
  snowp<-c(1.2,1,0,2.5)
  snownames<-c('SFCF','Tr','Tt','fm')
  snowdesc<-c('snowfall correction factor [-]','solid and liquid precipitation threshold temperature [ºC]','melt temperature [ºC]','snowmelt factor [mm/°C.delta t]')
  snowparams<-data.frame(p=snowp,type='snow',desc=snowdesc)
  dimnames(snowparams)[[1]]<-snownames
  snowparams


  #Soil
  soilp<-c(200,0.8,1.15)
  soilnames<-c('FC','LP','beta')
  soildesc<-c('fictitious soil field capacity [mm]','parameter to get actual ET [-]; 0<LP<1','parameter that allows non-linear relations between input-output runoff [-]')
  soilparams<-data.frame(p=soilp,type='soil',desc=soildesc)
  dimnames(soilparams)[[1]]<-soilnames
  soilparams

  #Route
  routep<-c(0.4,0.15,0.002,0.9,0.1)
  routenames<-c('K0','K1','K2','UZL','PERC')
  routedesc<-c('STZ storage constant [1/delta t]; 1>K0','SUZ storage constant [1/delta t]; K0>K1','SLZ storage constant [1/delta t]; K1>K2','maximum rate flux between STZ and SUZ; [mm/delta t]','maximum rate flux between STZ and SUZ [mm/delta t]; UZL>PERC')

  #UZL er vannstand i karet?
  routeparams<-data.frame(p=routep,type='route',desc=routedesc)
  dimnames(routeparams)[[1]]<-routenames
  routeparams

  #UH
  UHp<-1.5
  UHnames<-'Bmax'
  UHdesc<-'base of the UH triangle [timestep]'
  UHparams<-data.frame(p=UHp,type='UH',desc=UHdesc)
  dimnames(UHparams)[[1]]<-UHnames
  UHparams

  PCp<-1.0
  PCnames<-'Pcorr'
  PCdesc<-'precipitation correction [factor]'
  PCparams<-data.frame(p=PCp,type='PC',desc=PCdesc)
  dimnames(PCparams)[[1]]<-PCnames
  PCparams

  #-----------------------------------------------#
  #Insert lower and upper bound of parameters
  #-----------------------------------------------#
  snowparams$l<-c(0.8,-2,-3,1)
  snowparams$u<-c(3,3,3,6)
  snowparams<-snowparams[,c('p','l','u','type','desc')]

  delta<-1e-4
  soilparams$l<-c(20,delta,1)
  soilparams$u<-c(600,1,4)
  soilparams<-soilparams[,c('p','l','u','type','desc')]
  soilparams

  routeparams$l<-c(0.3,0.1,0.001,5.0,0.1)
  routeparams$u<-c((1-delta),0.3,0.1,100,2.5)
  routeparams<-routeparams[,c('p','l','u','type','desc')]
  routeparams

  UHparams$l<-1
  UHparams$u<-7
  UHparams<-UHparams[,c('p','l','u','type','desc')]
  UHparams


  PCparams$l<-0.5
  PCparams$u<-3.0
  PCparams<-PCparams[,c('p','l','u','type','desc')]
  PCparams


  HBVpars<-rbind(snowparams,soilparams,routeparams,UHparams, PCparams)
  #------------------------------------------------------------------------------#

  HBVpars$l[8:10]=rep(0.001,3)
  HBVpars$u[8:10]=c(0.999,0.999,0.999)

  HBVpars$l[11:12]=c(0.01,2)
  HBVpars$u[11:12]=c(100,100)

  return(HBVpars)

}
