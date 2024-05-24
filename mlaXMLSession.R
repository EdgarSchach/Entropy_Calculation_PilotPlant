# require("XML")
#### not any more in use: require("DescTools")


#' Checker of the status of the current grouping scheme
#'
#' Checks whether the current grouping system is active. Gropuings are defined through an XML
#' file like those used in MLA DataViewer
#'
#' @return a boolean: TRUE if active, FALSE if inacive.
#' @export
#'
#' @examples
#'  mlaGetGroupingStatus()
mlaGetGroupingStatus = function(){
  return(options()$mlaGrouping$active)
}



#' Setter of the status of the current grouping scheme
#'
#' Sets the current grouping system as active (TRUE) or inactive (FALSE). Gropuings are defined through an XML
#' file like those used in MLA DataViewer
#'
#' @param status boolean, is the status active?
#'
#' @return a boolean answering the question: was the status correcty set? if FALSE, then a warning is generated
#' explaining what went wrong.
#' @export
#'
#' @examples
#'  mlaSetGroupingStatus(TRUE)
mlaSetGroupingStatus = function(status){
  if(is.logical(status) & length(status)==1){
    o = options()$mlaGrouping
    if(length(o)==0){
      warning("no grouping yet defined")
      return(FALSE)
    }
    o$active = status
    options(mlaGrouping=o)
    return(TRUE)
  }else{
    warning("status must be a single logical value: nothing done")
    return(FALSE)
  }
}


#' Setter of the grouping scheme from an XML file
#'
#' Sets a new grouping system as default, and activates it. Gropuings are defined through an XML
#' file like those used in MLA DataViewer.
#'
#' @param file name of the XML fil defining the grouping
#'
#' @return invisibly, the grouping structure, namely a list of three elements:
#' \begin{itemize}
#'  \item \code{active}: a boolean, is this grouping active?
#'  \item \code{groupoingMatrix}: a matrix giving the contribution of each original mineral (in rows) to each mineral group (in columns)
#'  \item \code{colorTable}: a vector of hexadecimal colors for each mineral group
#' \end{itemize}
#' @export
#'
#' @examples
#' print(mlaSetXMLGroupsAndPalette("./R/grouping.xml"))
mlaSetXMLGroupsAndPalette = function(file="grouping.xml"){
 xmlGroups = XML::xmlParse(file)
 nodeList = XML::xmlToList(xmlGroups)

 # set of group names
 nameSet = sapply(nodeList[1:(length(nodeList)-2)], function(x){
   x[[".attrs"]]
 })

 # red coordinate for each group
 redSet = sapply(1:(length(nodeList)-2), function(i){
   as.integer(nodeList[[i]][["Group"]][[1]][3])
 })
 # green coordinate for each group
 greenSet = sapply(1:(length(nodeList)-2), function(i){
   as.integer(nodeList[[i]][["Group"]][[2]][3])
 })
 # blue coordinate for each group
 blueSet = sapply(1:(length(nodeList)-2), function(i){
   as.integer(nodeList[[i]][["Group"]][[3]][3])
 })
 # mineral names belonging to each group
 mineralSets = lapply(1:(length(nodeList)-2), function(i){
  xx = nodeList[[i]]
  xn = which(names(xx)=="Entry")
  sapply(xx[xn], function(x) x["Name"])
 })

 # create matrix and vector representations of the grouping table
 names(mineralSets) = nameSet
 groupingTable = sapply(1:length(mineralSets), function(i) rep(names(mineralSets)[i], length(mineralSets[[i]])))
 groupingTable = unlist(groupingTable)
 names(groupingTable) = unlist(mineralSets)

 groupingMatrix = matrix(0, ncol=length(mineralSets), nrow=length(groupingTable))
 rownames(groupingMatrix) = names(groupingTable)
 colnames(groupingMatrix) = names(mineralSets)

 for(i in 1:length(mineralSets)){
  groupingMatrix[mineralSets[[i]],names(mineralSets)[i]]=1
 }

 # create color vector
 colorTable = rgb(red=redSet, green=greenSet, blue=blueSet, maxColorValue=255, names=names(mineralSets))

 mlaGrouping = list(active=TRUE, groupingMatrix= groupingMatrix, colorTable= colorTable)
 options(mlaGrouping = mlaGrouping)
 invisible(mlaGrouping)
}





#' Reads a grouping scheme from an XML file
#'
#' Reads a grouping system from a file. Gropuings are defined through an XML
#' file like those used in MLA DataViewer.
#'
#' @param file name of the XML fil defining the grouping
#'
#' @return the \code{groupoingMatrix}, a matrix giving the contribution of each original mineral
#' (in rows) to each mineral group (in columns). This is NEITHER stored in any sense in R NOR it
#' conditions further computations, and is solely provided as a helping means for reading XML
#' grouping files. Use mlaSetXMLGroupsAndPalette(file) for reading, storing an activating a grouping
#' structure.
#'
#' @export
#'
#' @examples
#' mlaGetXMLGroups("./R/grouping.xml")
mlaGetXMLGroups = function(file="grouping.xml"){
  xmlGroups = XML::xmlParse(file)
  nodeList = XML::xmlToList(xmlGroups)

  # set of group names
  nameSet = sapply(nodeList[1:(length(nodeList)-2)], function(x){
    x[[".attrs"]]
  })

  # red coordinate for each group
  redSet = sapply(1:(length(nodeList)-2), function(i){
    as.integer(nodeList[[i]][["Group"]][[1]][3])
  })
  # green coordinate for each group
  greenSet = sapply(1:(length(nodeList)-2), function(i){
    as.integer(nodeList[[i]][["Group"]][[2]][3])
  })
  # blue coordinate for each group
  blueSet = sapply(1:(length(nodeList)-2), function(i){
    as.integer(nodeList[[i]][["Group"]][[3]][3])
  })
  # mineral names belonging to each group
  mineralSets = lapply(1:(length(nodeList)-2), function(i){
    xx = nodeList[[i]]
    xn = which(names(xx)=="Entry")
    sapply(xx[xn], function(x) x["Name"])
  })

  # create matrix and vector representations of the grouping table
  names(mineralSets) = nameSet
  groupingTable = sapply(1:length(mineralSets), function(i) rep(names(mineralSets)[i], length(mineralSets[[i]])))
  groupingTable = unlist(groupingTable)
  names(groupingTable) = unlist(mineralSets)

  groupingMatrix = matrix(0, ncol=length(mineralSets), nrow=length(groupingTable))
  rownames(groupingMatrix) = names(groupingTable)
  colnames(groupingMatrix) = names(mineralSets)

  for(i in 1:length(mineralSets)){
    groupingMatrix[mineralSets[[i]],names(mineralSets)[i]]=1
  }

  return(groupingMatrix)
}

