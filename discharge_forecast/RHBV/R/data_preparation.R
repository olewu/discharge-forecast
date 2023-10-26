#' prepare_data
#' @description Imputes missing values and/or return the longest sequence of non-missing data. Imputation method: Linear interpolation through the imputeTS package.
#' @param dt A data table with columns date (daily resolution) and value. The value includes data, e.g.  streamflow, precip or temp, that should be imputed/prepared. Other columns than value and date can also be present.
#' @param maxgap Maximum day gap allowed for imputation. If the gap is larger than maxgap, missing values are not predicted/imputed. Instead, a data table with the longest data sequence is returned.
#' @return The data table with missing values imputed values. If there are too many missing values in a row, the longest sequence of data without missing values is returned.
#' @export
#'
#' @examples
prepare_data=function(dt,maxgap=5){
  dt=copy(dt)
  dt=dt[order(date),]
  dt=dt[is.na(value)==FALSE,]
  dt[,daydiff:=c(1,diff(date)),]

  #date og value må være med#
  timestart=min(dt$date)
  timeend=max(dt$date)


  uniquediffs=unique((dt$daydiff)[is.na(dt$daydiff)==FALSE])

  if(length(uniquediffs)==1){
    dt[,daydiff:=NULL]
    return(dt)
  }

  alldates=seq(timestart,timeend,by=1)
  alldates=data.table(date=alldates)
  dt=merge(alldates,dt,by="date",all.x=TRUE)
  dt=dt[order(date),]
  dt[,value:=na_interpolation(value,maxgap=maxgap)]


  dt[,has_obs:=1-is.na(value)]



  #--------Find longest sequence of data-------------#
  run_lengths = rle(dt$has_obs == 1)
  (result = max(run_lengths$lengths[run_lengths$values]))
  only_a = ifelse(run_lengths$values, run_lengths$lengths, 0)
  longest_run_index = which.max(only_a)
  index = sum(run_lengths$lengths[seq_len(longest_run_index - 1)]) + 1
  dt=dt[index : (index + result - 1),]


  #dt[,daydiff:=c(1,diff(date))]
  dt[,daydiff:=NULL]
  dt[,has_obs:=NULL]



  return(dt)
}

