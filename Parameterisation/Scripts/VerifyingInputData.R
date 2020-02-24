# This script reads the tables describing the food web parameterisation and verify that the modelled species and the interaction lists match.
rm(list=ls(all=TRUE))

library(RColorBrewer)

E <- 0.1 # conversion efficiency

# define working directories
PAR_DIR <- "../OutputTables"

setwd(PAR_DIR)
PopParam <- read.csv("YearPopParam.csv", header = TRUE, stringsAsFactors = FALSE)
IntParamI <- read.csv("YearIntParam_TypeI.csv", header = TRUE, stringsAsFactors = FALSE)
IntParamII <- read.csv("YearIntParam_TypeII.csv", header = TRUE, stringsAsFactors = FALSE)

PreyList_FRI <- sort(unique(IntParamI$LowerTaxon))
PreyList_FRII <- sort(unique(IntParamII$LowerTaxon))

PredList <- sort(unique(c(IntParamI$UpperTaxon, IntParamII$UpperTaxon)))

# are prey lists the same?
match(PreyList_FRI, PreyList_FRII)
match(PreyList_FRI, PreyList_FRII)-seq(1:length(PreyList_FRI))
PreyList <- sort(unique(c(PreyList_FRI, PreyList_FRII)))
# are interaction list the same?
setdiff(IntParamI$LowerTaxon, IntParamII$LowerTaxon)
setdiff(IntParamI$UpperTaxon, IntParamII$UpperTaxon)

# is the prey list in PopParam the same as in the interaction list?
PreyListPopParam <- PopParam$Taxon[PopParam$TrophicLevel == "Prey"]
match(PreyListPopParam, PreyList_FRI)
any(is.na(match(PreyListPopParam, PreyList_FRI)))
PopParam <- PopParam[match(c(PredList, PreyList), PopParam$Taxon), ]

write.csv(PopParam, file = "YearPopParam.csv", row.names = FALSE, quote = FALSE)
