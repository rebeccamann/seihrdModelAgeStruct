
.loadParameters <- function(){
  parametersPath <- system.file("parameters.json", package = "simpleSeihrdModel", mustWork = TRUE)
  params <- jsonlite::fromJSON(parametersPath)
  
  youngAgeGrp = params$'0-19'
  adultAgeGrp = params$'20-64'
  elderlyAgeGrp = params$'65+'
  
  allAgeGrpParams = list(youngAgeGrp, adultAgeGrp, elderlyAgeGrp)
  return(allAgeGrpParams)
}
