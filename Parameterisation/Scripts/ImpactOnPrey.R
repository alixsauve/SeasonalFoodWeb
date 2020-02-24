rm(list = ls())

DAT_DIR <- "../../Data/BPF_Data"
WK_DIR <- "../../Data/ProcessedData"

# load data
setwd(WK_DIR)
PredPref <- read.csv("PredPref.csv", header = TRUE, stringsAsFactors = FALSE)
SpQuantities <- read.csv("SpDensitySeasons_2.csv", header = TRUE, stringsAsFactors = FALSE)

setwd(DAT_DIR)
PredInfo <- read.csv("Table5-1.csv", header = TRUE, dec = ",", stringsAsFactors = FALSE)
PredInfo$NdayPresent_W[is.na(PredInfo$NdayPresent_W)] <- 0
PredInfo$NdayPresent_S[is.na(PredInfo$NdayPresent_S)] <- 0

PreyList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Prey"])
PredList <- sort(PredInfo$Taxon)
NPred <- nrow(PredInfo)

# verify order of species
SpQuantities <- SpQuantities[match(c(PredList, PreyList), SpQuantities$Taxon),]

ImpactOnPrey <- data.frame(LowerClade = character(), LowerTaxon = character(), UpperClade = character(), UpperTaxon = character(),
                           IntakePrey_W = numeric(), IntakePrey_S = numeric())

for (Sp in 1:NPred){
  ImpactOnPrey_Sp <- subset(PredPref, UpperTaxon == PredList[Sp])
  ImpactOnPrey_Sp <- ImpactOnPrey_Sp[order(ImpactOnPrey_Sp$LowerTaxon),]
  
  NDaysSSp <- PredInfo$NdayPresent_S[Sp]
  NDaysWSp <- PredInfo$NdayPresent_W[Sp]
  
  # calculate the seasonal intake of prey
  if ((!"IntakePrey_W" %in% colnames(ImpactOnPrey_Sp)) & (!"IntakePrey_S" %in% colnames(ImpactOnPrey_Sp))){
    ImpactOnPrey_Sp$IntakePrey_W <- NA
    ImpactOnPrey_Sp$IntakePrey_S <- NA
  }
  
  DFISp <- PredInfo$DailyFoodIntake[Sp]

  # making sure that the sum of the diet fractions is 1 (might not be the case when not all interactions are reported)
  ImpactOnPrey_Sp$FracBioInDiet_W <- ImpactOnPrey_Sp$FracBioInDiet_W/sum(ImpactOnPrey_Sp$FracBioInDiet_W, na.rm = TRUE) * 100
  ImpactOnPrey_Sp$FracBioInDiet_S <- ImpactOnPrey_Sp$FracBioInDiet_S/sum(ImpactOnPrey_Sp$FracBioInDiet_S, na.rm = TRUE) * 100
  
  ImpactOnPrey_Sp$IntakePrey_W <- DFISp * ImpactOnPrey_Sp$FracBioInDiet_W / 100 * NDaysWSp
  ImpactOnPrey_Sp$IntakePrey_W[is.na(ImpactOnPrey_Sp$IntakePrey_W)] <- 0
  ImpactOnPrey_Sp$IntakePrey_S <- DFISp * ImpactOnPrey_Sp$FracBioInDiet_S / 100 * NDaysSSp
  ImpactOnPrey_Sp$IntakePrey_S[is.na(ImpactOnPrey_Sp$IntakePrey_S)] <- 0
  
  ImpactOnPrey_Sp <- ImpactOnPrey_Sp[order(ImpactOnPrey_Sp$LowerTaxon),]
  ImpactOnPrey_Sp$UpperClade <- PredInfo$Clade[Sp]; ImpactOnPrey_Sp$LowerClade <- SpQuantities$Clade[match(ImpactOnPrey_Sp$LowerTaxon, SpQuantities$Taxon)]

  ImpactOnPrey_Sp <- subset(ImpactOnPrey_Sp, select = colnames(ImpactOnPrey))
  ImpactOnPrey <- rbind(ImpactOnPrey, ImpactOnPrey_Sp)
  
  if (sum(ImpactOnPrey_Sp$IntakePrey_W, na.rm = TRUE) == 0){
    print(Sp)
  }
}

ImpactOnPrey$IntakePrey_W <- round(ImpactOnPrey$IntakePrey_W, 4)
ImpactOnPrey$IntakePrey_S <- round(ImpactOnPrey$IntakePrey_S, 4)
ImpactOnPrey$IntakePrey_Unit <- "g/y"

setwd(WK_DIR)
ImpactOnPrey$IntakePrey_W[SpQuantities$IsMigratoryBird[match(ImpactOnPrey$LowerTaxon, SpQuantities$Taxon)]] <- 0
ImpactOnPrey <- ImpactOnPrey[(ImpactOnPrey$IntakePrey_W != 0) | (ImpactOnPrey$IntakePrey_S != 0), ]
write.csv(ImpactOnPrey, file = "PerCapitaIntakePreyBiomass.csv", quote = FALSE, row.names = FALSE)

PreyList <- sort(unique(ImpactOnPrey$LowerTaxon))
PredList <- sort(unique(ImpactOnPrey$UpperTaxon))
SpQuantities <- SpQuantities[!is.na(match(SpQuantities$Taxon, c(PredList, PreyList))), ]
write.csv(SpQuantities, file = "SpDensitySeasonsFinal.csv", row.names = FALSE, quote = FALSE)