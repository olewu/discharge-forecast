#' make_basin_data
#'
#' @param areas
#' @param refh The height of location of the precip/temp gauge. If you use gridded data from e.g. seNorge, use mean elevation in the area.
#'
#' @return A table with catchment information on a suitable format for further use.
#' @export
#'
#' @examples
make_basin_data=function(areas,refh){
  areas$h<-(areas$from+areas$to)/2

  #----------------------------------------------------------------------------------------#
  # distribute ther data in 10 elevation bands
  # set columns
  nareas<-nrow(areas)
  col_eb    <- paste0('eb_', 1:nareas) # elevation band names
  col_area  <- areas$areal
  col_rel   <- round( col_area / sum(col_area), 4)
  col_heigh <- areas$h

  basin_data <- data.frame('elev_band' = col_eb, 'area(km2)' = col_area,
                           'rel_area' = col_rel, 'h(masl)' = col_heigh,'refh'=refh,
                           check.names = F)
  return(basin_data)
}

make_elevation_bands=function(total_area,hypsovec){
  #use the returned dt as input in basin_data (areas).

  from=hypsovec[1:(length(hypsovec)-1)]
  to=hypsovec[2:(length(hypsovec))]
  subarea=total_area/(length(hypsovec)-1)

  basin=data.table(from=from,to=to,areal=subarea)

  return(basin)
}


#' distribute_weather()
#'
#' @param tempvec
#' @param precipvec
#' @param datevec
#' @param basin_data
#' @param latitude
#'
#' @return
#' @export
#'
#' @examples
distribute_weather=function(tempvec,precipvec,datevec,basin_data,latitude){
  time_step="daily"
  latitude #for evaporation calculations

  if(mean(tempvec,na.rm=TRUE)>100){
    tempvec=tempvec-273.15
    print("Convert Kelvin temperatures to Celsius.")
  }


  nareas=dim(basin_data)[1]
  refh=basin_data$refh[1]
  col_heigh=basin_data$`h(masl)`

  # meteorological data distribution
  mat_tair <- matrix(0, nrow = length(tempvec), ncol = nareas)
  mat_prec <- matrix(0, nrow = length(precipvec), ncol = nareas)
  mat_pet  <- matrix(0, nrow = length(tempvec), ncol = nareas)

  #uses functions Temp_model and Precip_model to manipulate temperature and precipitation by elevation band
  #Elevation band by elevation band...
  for(i in 1:nareas){
    #Temperature decreases with 0.6 degree for each 100m increase if no precipitation. References to height of gauge refh
    t_noprecip<-Temp_model(model = 1, inputData = tempvec,zmeteo = refh, ztopo = col_heigh[i], param = c(-0.6) )
    #Temperature decreases with 0.5 degree for each 100m increase if precipitation
    t_precip<-Temp_model(model = 1, inputData = tempvec,zmeteo = refh, ztopo = col_heigh[i], param = c(-0.5) )
    #Create a temporar variable holding the no-precipitation temps
    tt<-t_noprecip
    #For days with precipitation the temps are replaced with noprecip temps
    tt[precipvec>0]<-t_precip[precipvec>0]
    mat_tair[ , i]<-tt

    #Precipitation increases with 8% for each 100 m increase. Reference to height of gauge refh.
    mat_prec[ , i] <- Precip_model(model = 1, inputData = precipvec,zmeteo = refh, ztopo = col_heigh[i], param = c(8) )

  }

  mat_pet=mat_prec*NA
  for(i in 1 : nareas){
     mat_pet[,i]<-PE_Oudin(yday(datevec), Temp=mat_tair[,i], Lat=latitude, LatUnit = "deg", TimeStepIn = time_step, TimeStepOut = time_step)
  }

  weatherlist=list(basin=basin_data,mat_pet=mat_pet,mat_tair=mat_tair,mat_prec=mat_prec,datevec=datevec)

  return(weatherlist)
}


#' calibrate_HBV()
#'
#' @param weatherlist
#' @param streamflow
#' @param param_limits
#' @param score_function
#'
#' @return
#' @export
#'
#' @examples
calibrate_HBV=function(weatherlist,streamflow,param_limits,score_function="MSE",burnin=365){
  HBVpars=param_limits

  system.time(f11<-DEoptim(fn=obj_fun_hbv,
                             lower=HBVpars$l,
                             upper=HBVpars$u,
                             control=DEoptim.control(NP=NA,itermax=600,reltol=1e-4,steptol=50,trace=10,parallelType = 0),
                           weatherlist=weatherlist,qq=streamflow,score_function=score_function,rout_model=3,rm=-c(1:burnin)))

  f11$optim$bestmem[8:10]=rev(sort(f11$optim$bestmem[8:10]))
  f11$optim$bestmem[11:12]=rev(sort(f11$optim$bestmem[11:12]))


  if(f11$optim$bestmem[8]==f11$optim$bestmem[9] | f11$optim$bestmem[8]==f11$optim$bestmem[10] | f11$optim$bestmem[9]==f11$optim$bestmem[10]){
    f11$optim$bestmem[8]=f11$optim$bestmem[8]+rnorm(1,mean=0,sd=0.0001)
    f11$optim$bestmem[9]=f11$optim$bestmem[9]+rnorm(1,mean=0,sd=0.0001)
    f11$optim$bestmem[10]=f11$optim$bestmem[10]+rnorm(1,mean=0,sd=0.0001)
    f11$optim$bestmem[f11$optim$bestmem[8:10]>1]=1
    f11$optim$bestmem[f11$optim$bestmem[8:10]<0]=0
    f11$optim$bestmem[8:10]=rev(sort(f11$optim$bestmem[8:10]))
  }

  if(f11$optim$bestmem[11]==f11$optim$bestmem[12]){
    f11$optim$bestmem[11]=f11$optim$bestmem[11]+rnorm(1,mean=0,sd=0.0001)
    f11$optim$bestmem[12]=f11$optim$bestmem[12]+rnorm(1,mean=0,sd=0.0001)
    f11$optim$bestmem[11:12]=rev(sort(f11$optim$bestmem[11:12]))
  }



  return(f11)
}



#' hydrological_hbv()
#'
#' @param basin
#' @param tair
#' @param precip
#' @param pet
#' @param param_snow
#' @param param_soil
#' @param param_route
#' @param param_tf
#' @param init_snow
#' @param init_soil
#' @param init_routing
#' @param returnDetails
#' @param rout_model
#'
#' @return
#' @export
#'
#' @examples
hydrological_hbv <- function(basin,tair,precip,pet,
                             param_snow,param_soil,param_route,param_tf,
                             init_snow = 0,init_soil = 0,init_routing = c(0, 0, 0) ,
                             returnDetails=FALSE,rout_model=3){
  ## brief arguments description
  # basin: data frame with the same structure of the data("semi_distributed_hbv) (colnames included).
  # tair: numeric matrix with air temperature inputs.
  # precip: numeric matrix with precipitation inputs.
  # pet: numeric matrix with potential eavapotranspiration inputs.
  # param_snow: numeric vector with snow module parameters.
  # param_soil: numeric vector with soil moisture parameters.
  # param_routing: numeric vector with the routing parameters.
  # param_tf: numeric vector with the transfer function parameter.
  # init_snow: numeric value with initial snow water equivalent. Default value being 20 mm.
  # init_soil: numeric value with initial soil moisture content. Default value being 100 mm.
  # init_routing: numeric vector with bucket water initial values. Default values are 0 mm.
  ## output
  # simulated streamflow series.

  #Sort the variables:
  param_route=c(rev(sort(param_route[1:3])),rev(sort(param_route[4:5])))

  if(param_route[1]==param_route[2] | param_route[1]==param_route[3] | param_route[3]==param_route[2]){
    param_route[1]=param_route[1]+rnorm(1,mean=0,sd=0.0001)
    param_route[2]=param_route[2]+rnorm(1,mean=0,sd=0.0001)
    param_route[3]=param_route[3]+rnorm(1,mean=0,sd=0.0001)
    param_route[param_route[1:3]>1]=1
    param_route[param_route[1:3]<0]=0
    param_route[1:3]=rev(sort(param_route[1:3]))
  }

  if(param_route[4]==param_route[5]){
    param_route[4]=param_route[4]+rnorm(1,mean=0,sd=0.0001)
    param_route[5]=param_route[5]+rnorm(1,mean=0,sd=0.0001)
    param_route[4:5]=rev(sort(param_route[4:5]))
  }

  n_it <- nrow(basin)
  # create output lists
  snow_module  <- list()
  soil_module  <- list()
  route_module <- list()
  tf_module    <- list()
  if(rout_model==3) init_routing = init_routing[1:2]
  # snow and soil module in every elevation band

  for(i in 1:n_it){
    snow_module[[ i ]] <-
      SnowGlacier_HBV(model = 1, inputData = cbind(tair[ , i], precip[ , i]),
                      initCond =  c(init_snow, 2), param = param_snow)
    soil_module[[ i ]] <-
      Soil_HBV(model = 1, inputData = cbind(snow_module[[i]][ , 5] , pet[ , i]),
               initCond = c(init_soil, basin[i, 'rel_area']), param = param_soil )

  } # end for

  # get total soil discharge
  soil_disch <- lapply(X = 1:n_it, FUN = function(x){
    out <- soil_module[[x]][ , 1]
  })
  soil_disch <- Reduce(f = `+`, x = soil_disch)

  # route module
  route_module <- Routing_HBV(model = rout_model, lake = F, inputData = as.matrix(soil_disch),
                              initCond = init_routing, param = param_route )

  # transfer function
  tf_module <- round(
    UH(model = 1, Qg = route_module[ , 1], param = param_tf), 4
  )
  if(returnDetails){
    retval<-list(snow_module=snow_module,
                 soil_module=soil_module,
                 route_module=route_module,
                 tf_module=tf_module)
  }
  else{
    retval<-tf_module
  }
  return(retval)
}






#' obj_fun_hbv()
#'
#' @param params
#' @param weatherlist
#' @param qq
#' @param rm
#' @param score_function
#' @param nsw
#' @param rout_model
#'
#' @return
#' @export
#'
#' @examples
obj_fun_hbv<-function(params,weatherlist,qq,rm=-c(1:365),score_function='MSE',nsw=0.3,rout_model=3){

  #-----Add something that makes it possible to handle "overlapping" parameter limits-----#


 # params[8:10]=rev(sort(params[8:10]))
  #params[11:12]=rev(sort(params[11:12]))



  hbv_sim<-hydrological_hbv(basin=weatherlist$basin,
                            tair=weatherlist$mat_tair,
                            precip=weatherlist$mat_prec,
                            pet=weatherlist$mat_pet,
                            param_snow=params[1:4],
                            param_soil=params[5:7],
                            param_route=params[8:12],
                            param_tf=params[13],
                            init_snow = 20,
                            init_soil = 0,
                            init_routing = c(0, 0, 0),
                            rout_model=rout_model)

  qsim<-hbv_sim
  qsim<-qsim[rm] #remove the beginning of the dataset affected by initial conditions.  "Burn in"

  qq<-qq[rm]

  #ANE:
  if(score_function=='ANE')retval<-ANEq(qq,qsim)
  #Mean squared error
  if(score_function=='MSE')retval<-MSEq(qq,qsim)
  #Nash and Sutcliffe
  if(score_function=='NS'){
    retval<--NSq(qq,qsim)
  }
  if(score_function=='LNS'){
    retval<--LNNSq(qq,qsim)
  }
  if(score_function=='VE'){
    retval<-abs(VEq(qq,qsim))
  }
  #Relativ bias
  if(score_function=='RB'){
    retval<-RBq(qq,qsim)
  }
  if(score_function=='COMB'){
    retval<--(1-nsw)*NSq(qq,qsim)+abs(VEq(qq,qsim))
  }
  #Correlation
  if(score_function=='Cor'){
    if(!is.na(sd(qsim)) & sd(qsim)>1e-6){
      retval<--cor(qq,qsim)
    }
    else retval<-0
    if(is.na(retval))retval<-0
  }

  return(retval)
}


#' sim_HBV()
#'
#' @param mypars
#' @param weatherlist
#' @param init_snow
#' @param init_soil
#' @param init_routing
#'
#' @return
#' @export
#'
#' @examples
sim_HBV=function(mypars,weatherlist,init_snow=0,init_soil=0,init_routing=c(0,0,0)){
  #mypars comes from calibrate_HBV, f11$optim$bestmem
  basin_data=weatherlist$basin
  mat_tair=weatherlist$mat_tair
  mat_pet=weatherlist$mat_pet
  mat_prec=weatherlist$mat_prec

  mydistfit<-hydrological_hbv(basin=basin_data,
                              tair=mat_tair,
                              precip=mat_prec,
                              pet=mat_pet,
                              param_snow=mypars[1:4],
                              param_soil=mypars[5:7],
                              param_route=mypars[8:12],
                              param_tf=mypars[13],
                              init_snow = init_snow,
                              init_soil = init_soil,
                              init_routing =init_routing,
                              returnDetails=TRUE)
  return(mydistfit)
}
