edgelist2Mat <- function(EdgeList, ColWeightName = NULL){
  # edgelist2Mat() returns a matrix with lower taxa in rows and upper taxa in columns
  # if ColWeightName is not null, non-zero elements of the matrix are filled with the corresponding weight
  if (!is.data.frame(EdgeList)){
    stop("EdgeList must be a dataframe with at least two columns ('LowerTaxon' and 'UpperTaxon').")
  }
  if ((is.null(EdgeList$LowerTaxon)) | (is.null(EdgeList$UpperTaxon))){
    stop("EdgeList must be a dataframe with at least two columns ('LowerTaxon' and 'UpperTaxon').")
  }
  if ((!is.null(ColWeightName)) & (!is.character(ColWeightName))){
    stop("ColWeightName must be a character string.")
  }
  
  LTList <- sort(unique(EdgeList$LowerTaxon))
  UTList <- sort(unique(EdgeList$UpperTaxon))
  
  FWMat <- table(subset(EdgeList, select = c(LowerTaxon, UpperTaxon)))
  if (!is.null(ColWeightName)){
    NPred <- ncol(FWMat)
    for (Sp in 1:NPred){
      IndexIntSp <- which(EdgeList$UpperTaxon == colnames(FWMat)[Sp])
      IndexIntSp <- IndexIntSp[order(EdgeList$LowerTaxon[IndexIntSp])]
      PreyOfSp <- EdgeList$LowerTaxon[IndexIntSp]
      FWMat[PreyOfSp, Sp] <- EdgeList[IndexIntSp, ColWeightName]
    }
  }
  return(FWMat)
}