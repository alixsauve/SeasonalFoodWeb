shannonSeason <- function(EdgeList, SummerWeightName, WinterWeightName){
  if (!is.data.frame(EdgeList)){
    stop("EdgeList must be a dataframe with at least two columns ('LowerTaxon' and 'UpperTaxon').")
  }
  if ((is.null(EdgeList$LowerTaxon)) | (is.null(EdgeList$UpperTaxon))){
    stop("EdgeList must be a dataframe with at least two columns ('LowerTaxon' and 'UpperTaxon').")
  }
  if ((!is.character(SummerWeightName)) | (!is.character(WinterWeightName))){
    stop("SummerWeightName and WinterWeightName must be character strings.")
  }
  LTList <- sort(unique(EdgeList$LowerTaxon)); LGSize <- length(LTList)
  UTList <- sort(unique(EdgeList$UpperTaxon)); UGSize <- length(UTList)
  
  HSeason <- vector("list", 2); names(HSeason) <- c("LowerGuild", "UpperGuild")
  # lower guild
  HSeason_LG <- c()
  for (Sp in LTList){
    IntSp <- which(EdgeList$LowerTaxon == Sp)
    TotWeights <- sum(unlist(subset(EdgeList, select = c(SummerWeightName, WinterWeightName), LowerTaxon == Sp)), na.rm = TRUE)
    SummerWeight <- sum(EdgeList[IntSp, SummerWeightName], na.rm = TRUE)
    WinterWeight <- sum(EdgeList[IntSp, WinterWeightName], na.rm = TRUE)
    SummerRat <- SummerWeight/TotWeights; WinterRat <- WinterWeight/TotWeights
    if ((SummerRat != 0) & (WinterRat != 0)){
      HSeason_Sp <- - SummerRat * log(SummerRat) - WinterRat * log(WinterRat)
    }
    else {
      HSeason_Sp <- 0
    }
    HSeason_LG <- c(HSeason_LG, HSeason_Sp)
  }
  HSeason$LowerGuild <- data.frame(Taxon = LTList, HYear = HSeason_LG)
  # upper guild
  HSeason_UG <- c()
  for (Sp in UTList){
    IntSp <- which(EdgeList$UpperTaxon == Sp)
    TotWeights <- sum(unlist(subset(EdgeList, select = c(SummerWeightName, WinterWeightName), UpperTaxon == Sp)), na.rm = TRUE)
    SummerWeight <- sum(EdgeList[IntSp, SummerWeightName], na.rm = TRUE)
    WinterWeight <- sum(EdgeList[IntSp, WinterWeightName], na.rm = TRUE)
    SummerRat <- SummerWeight/TotWeights; WinterRat <- WinterWeight/TotWeights
    if ((SummerRat != 0) & (WinterRat != 0)){
      HSeason_Sp <- - SummerRat * log(SummerRat) - WinterRat * log(WinterRat)
    }
    else {
      HSeason_Sp <- 0
    }
    HSeason_UG <- c(HSeason_UG, HSeason_Sp)
  }
  HSeason$UpperGuild <- data.frame(Taxon = UTList, HYear = HSeason_UG)
  
  return(HSeason)
}