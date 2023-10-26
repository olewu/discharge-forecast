#' Find donor catchments for your ungauged target catchment based on catchment properties
#'
#' @param our_properties A data table with (standardized) catchment attributes of the target catchment (typically ungauged) and donor catchments (gauged). One of the columns should be "stat_id" and include catchment id numbers.
#' @param my_id The id of the target catchment (ungauged).
#' @param properties_to_include Vector with the name of the attributes/covariates that should be included in the donor search. Should be a subset of the column names in "theirproperties" and "myproperties". No categorical variables allowed.
#' @param to_scale True if the properties should be standardized (recommended).
#'
#' @return A data table similar to "our_properties" with the catchments' (L2-)distance to the target catchment added in a new column. The distance is here taken across all catchments attributes mentioned in "properties_to_include". The resulting table is sorted according to distance. The target catchment is included at the top with distance 0.
#' @export
#'
#' @examples
#' catchprop=fread(file="/nr/project/stat/ClimateFutures/RenewableEnergy/Smaakraft/HBVpackage/RHBV/data/catchprop.csv")
#' donortab=find_my_donors(our_properties,my_id="30.3.0",properties_to_include=c("utm_east_z33","utm_north_z33"),to_scale=TRUE)
find_my_donors=function(our_properties,my_id,properties_to_include,to_scale=TRUE){
  stat_ids=our_properties$stat_id
  our_important_properties=our_properties[,.SD,.SDcols=colnames(our_properties)%in%properties_to_include]
  our_important_properties=data.table(apply(our_important_properties,2,as.numeric))

  if(to_scale==TRUE){
    our_important_properties=data.table(scale(our_important_properties))
  }


  my_important_properties=our_important_properties[which(stat_ids==my_id),]

  if(dim(my_important_properties)[1]==0){return(NULL)}


  their_important_properties=our_important_properties

  distvec=rep(0,dim(their_important_properties)[1])
  for(jj in 1:(dim(my_important_properties)[2])){
    number=(unlist(my_important_properties[,.SD,.SDcols=jj])-unlist(their_important_properties[,.SD,.SDcols=jj]))^2

    if(is.na(mean(number,na.rm=TRUE))==FALSE){
      distvec=distvec+number
    }
  }

  distvec=sqrt(distvec)
  distancetable=cbind(distance=distvec,our_properties)
  distancetable=distancetable[order(distance)]

  return(distancetable)

}
