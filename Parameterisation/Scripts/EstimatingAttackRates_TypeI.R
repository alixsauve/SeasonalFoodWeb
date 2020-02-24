# In this script, we infer the attack rates as defined in a functional response of type I.

rm(list=ls(all=TRUE))

library(RColorBrewer)

E <- 0.1 # conversion efficiency

# define working directories
DAT_DIR <- "../../Data/BPF_Data"
WK_DIR <- "../../Data/ProcessedData"
PAR_DIR <- "../../Parameterisation/OutputTables"

# load data
setwd(DAT_DIR)
PredInfo <- read.csv("Table5-1.csv", header = TRUE, dec = ",", stringsAsFactors = FALSE)
SubsetFW <- read.csv("IntListSubsets.csv", header = TRUE, stringsAsFactors = FALSE)

setwd(WK_DIR)
IntakePrey <- read.csv("PerCapitaIntakePreyBiomass.csv", header = TRUE, stringsAsFactors = FALSE)
SpDensSeason <- read.csv("SpDensitySeasonsFinal.csv", header = TRUE, stringsAsFactors = FALSE)

# species list
PreyList <- sort(SpDensSeason$Taxon[SpDensSeason$Trophic_level == "Prey"])
PredList <- sort(PredInfo$Taxon)
NPred <- nrow(PredInfo)

SpDensSeason <- SpDensSeason[match(c(PredList, PreyList), SpDensSeason$Taxon),]

# here we force predators' impact to be 0 if the bird migrates during winter
IntakePrey$IntakePrey_W[SpDensSeason$IsMigratoryBird[match(IntakePrey$LowerTaxon, SpDensSeason$Taxon)]] <- 0

# some numbers on community size
IntakePrey$IntNo <- 1:nrow(IntakePrey)
# summer
IntakePreySummer <- subset(IntakePrey, IntakePrey_S > 0)
PreySummer <- unique(IntakePreySummer$LowerTaxon)
PredSummer <- unique(IntakePreySummer$UpperTaxon)
# winter
IntakePreyWinter <- subset(IntakePrey, IntakePrey_W > 0)
PreyWinter <- unique(IntakePreyWinter$LowerTaxon)
PredWinter <- unique(IntakePreyWinter$UpperTaxon)
# shared species and interactions
length(intersect(PredWinter, PredSummer))
length(intersect(PreyWinter, PreySummer))
length(intersect(IntakePreySummer$IntNo, IntakePreyWinter$IntNo))
# season specific species and interactions
length(which(!PreySummer %in% PreyWinter))
length(which(!PredSummer %in% PredWinter))
length(which(!(IntakePreySummer$IntNo %in% IntakePreyWinter$IntNo)))
# specific to winter
length(which(!PreyWinter %in% PreySummer))
length(which(!PredWinter %in% PredSummer))
length(which(!(IntakePreyWinter$IntNo %in% IntakePreySummer$IntNo)))


FracS <- 168/365
FracW <- 197/365

# make sure that SubsetFW contains the same prey than SpDensSeason
IndexInt <- which((!is.na(match(SubsetFW$LowerTaxon, IntakePrey$LowerTaxon))) &
                    (!is.na(match(SubsetFW$UpperTaxon, IntakePrey$UpperTaxon))))
SubsetFW <- SubsetFW[IndexInt,]

# calculate attack rate for each interaction
# create and initiate two columns for attack rates
IntakePrey$AttackRate_S <- NA; IntakePrey$AttackRate_W <- NA
for (Sp in PredList){
  # list of prey
  PreyOfSp <- sort(unique(IntakePrey$LowerTaxon[IntakePrey$UpperTaxon == Sp]))
  # their densities
  DensPrey_S <- SpDensSeason$Summer_density[match(PreyOfSp, SpDensSeason$Taxon)]
  DensPrey_W <- SpDensSeason$Winter_density[match(PreyOfSp, SpDensSeason$Taxon)]
  
  # prey for which seasonal densities are missing
  IndexMissingDens <- which((is.na(DensPrey_S)) & (is.na(DensPrey_W)))# & (is.na(DensPrey_Sp)))
  DensPrey_S[IndexMissingDens] <- SpDensSeason$Mean_density[match(PreyOfSp[IndexMissingDens], SpDensSeason$Taxon)]
  DensPrey_W[IndexMissingDens] <- DensPrey_S[IndexMissingDens]
  
  # their body masses
  BodyMassPrey <- SpDensSeason$Body_mass[match(PreyOfSp, SpDensSeason$Taxon)]
  
  # their biomass
  BiomassPrey_S <- DensPrey_S*BodyMassPrey
  BiomassPrey_W <- DensPrey_W*BodyMassPrey
  
  # for taxa with fixed densities
  IndexFixedDens <- which((PreyOfSp == "Carcass") | (PreyOfSp == "Hymenoptera") | (PreyOfSp == "Lumbricidae_sp."))
  BiomassPrey_S[IndexFixedDens] <- 1; BiomassPrey_W[IndexFixedDens] <- 1; #BiomassPrey_Sp[IndexFixedDens] <- 1
  
  # the attack rates
  IntakePrey$AttackRate_S[IntakePrey$UpperTaxon == Sp] <- IntakePrey$IntakePrey_S[IntakePrey$UpperTaxon == Sp] / (FracS*BiomassPrey_S)
  IntakePrey$AttackRate_W[IntakePrey$UpperTaxon == Sp] <- IntakePrey$IntakePrey_W[IntakePrey$UpperTaxon == Sp] / (FracW*BiomassPrey_W)
}

which((IntakePrey$AttackRate_W == 0) & (IntakePrey$IntakePrey_W != 0))

# some species are predated specifically during winter, but we do not know their densities during winter
# we use then their spring densities
IndexWSpecificInt <- which((IntakePrey$IntakePrey_S == 0) & (IntakePrey$IntakePrey_W != 0) & (is.na(IntakePrey$AttackRate_W)))
DensSpOfWSpecificInt <- SpDensSeason$Spring_density[match(IntakePrey$LowerTaxon[IndexWSpecificInt], SpDensSeason$Taxon)]
BodyMassSpOfWSpecificInt <- SpDensSeason$Body_mass[match(IntakePrey$LowerTaxon[IndexWSpecificInt], SpDensSeason$Taxon)]
IntakePrey$AttackRate_W[IndexWSpecificInt] <- IntakePrey$IntakePrey_W[IndexWSpecificInt] / (FracW * BodyMassSpOfWSpecificInt * DensSpOfWSpecificInt)

IntakePrey$AttackRate_S[is.na(IntakePrey$AttackRate_S)] <- 0
IntakePrey$AttackRate_W[is.na(IntakePrey$AttackRate_W)] <- 0

# first complete data frame with species densities to ease data handling while drawing the figures
IntakePrey$LowerTaxonBodyMass <- NA; IntakePrey$UpperTaxonBodyMass <- NA
for (Sp in PredList){
  BodyMassSp <- SpDensSeason$Body_mass[SpDensSeason$Taxon == Sp]
  IntakePrey$UpperTaxonBodyMass[IntakePrey$UpperTaxon == Sp] <- BodyMassSp
  PredInfo$BodyMass[PredInfo$Taxon == Sp] <- BodyMassSp
}
for (Sp in PreyList){
  BodyMassSp <- SpDensSeason$Body_mass[SpDensSeason$Taxon == Sp]
  IntakePrey$LowerTaxonBodyMass[IntakePrey$LowerTaxon == Sp] <- BodyMassSp
}
IntakePrey$PredToPreyBodyMass <- IntakePrey$UpperTaxonBodyMass / IntakePrey$LowerTaxonBodyMass

# now calculate mean killing rates and magnitude of fluctuations
# convert NA into 0
IntakePrey$IntakePrey_S[is.na(IntakePrey$IntakePrey_S)] <- 0
IntakePrey$IntakePrey_W[is.na(IntakePrey$IntakePrey_W)] <- 0

# create a parameter dataframe for each predation frequencies
FreqPredList <- vector("list", 5)
FreqPredList[2:5] <- c("Season", "Monthly", "Weekly", "Daily")
for (FreqPred in FreqPredList){
  FreqPred <- unlist(FreqPred)
  IntakePreyFreq <- IntakePrey
  if (!is.null(FreqPred)){
    IntakePreyFreq$IntakePrey_W[!SubsetFW[[paste(FreqPred, "W", sep = "")]]] <- 0
    IntakePreyFreq$IntakePrey_S[!SubsetFW[[paste(FreqPred, "S", sep = "")]]] <- 0
    IntakePreyFreq <- IntakePreyFreq[(IntakePreyFreq$IntakePrey_S != 0) | (IntakePreyFreq$IntakePrey_W != 0), ]
  }
  
  # create a new dataframe for readability
  PredParam <- subset(IntakePrey, select = -c(IntakePrey_S, IntakePrey_W))
  PredParam$MeanG <- NA; PredParam$Epsilon <- NA; PredParam$SeasonDom <- NA
  
  # calculate mean killing rate
  PredParam$MeanG <- (IntakePrey$AttackRate_W + IntakePrey$AttackRate_S) / 2
  # calculate magnitude of the seasonal forcing
  PredParam$Epsilon <- abs(IntakePrey$AttackRate_W - IntakePrey$AttackRate_S)/(IntakePrey$AttackRate_W + IntakePrey$AttackRate_S)
  
  # summer dominant interactions
  IndexIntS <- which(IntakePrey$AttackRate_S >= IntakePrey$AttackRate_W)
  PredParam$SeasonDom[IndexIntS] <- "S"
  
  # for winter dominant interactions
  IndexIntW<- which(IntakePrey$AttackRate_W > IntakePrey$AttackRate_S)
  PredParam$SeasonDom[IndexIntW] <- "W"
  
  # format the data frame to match our needs in the simulation script
  PredParam <- subset(PredParam, select = c(LowerTaxon, UpperTaxon, SeasonDom, MeanG, Epsilon))
  colnames(PredParam) <- c("LowerTaxon", "UpperTaxon", "SeasonDom", "G", "E")
  # PredParam <- subset(PredParam, select = c(LowerTaxon, UpperTaxon, Season, SeasonDom, MeanG, Epsilon))
  # colnames(PredParam) <- c("LowerTaxon", "UpperTaxon", "Season", "SeasonDom", "G", "E")
  PredParam$E[is.na(PredParam$E)] <- 0
  
  setwd(PAR_DIR)
  write.csv(PredParam, file = paste("YearIntParam_TypeI", FreqPred, ".csv", sep = ""), row.names = FALSE)
}

# represent the parameter values ruling predation
UTBodyMassRange <- c(1e1, 1e6)
BodyMassRatioRange <- c(1e-2, 2e4)
AttackRateRange <- c(1e-3, 1e6)
HandlingTimeRange <- c(1e-7, 1e-3)
IntCol <- brewer.pal(2, "Set1")

PredParam$IntCol <- IntCol[2];
PredParam$IntCol[IntakePrey$UpperClade == "Raptor"] <- IntCol[1]

PredInfo$SpCol <- IntCol[1]; PredInfo$SpCol[PredInfo$Clade != "Raptor"] <- IntCol[2]

plot(1, 1, log = "xy", type = "n", xlim = BodyMassRatioRange, ylim = AttackRateRange,
     xlab = "Body mass ratio", ylab = "Summer discovery rate (ha/g/y)")
points(IntakePrey$PredToPreyBodyMass, PredParam$G, pch = 20,
       col = adjustcolor(PredParam$IntCol, alpha.f = 0.5))
legend("topleft", legend = c("Raptor", "Mammal predator"), bty = "n", pch = 20, col = IntCol)

SpDensInitSimu <- subset(SpDensSeason, select = c(CladeBis, Clade, Taxon, Trophic_level, Body_mass, Spring_density))
colnames(SpDensInitSimu) <- c("CladeSimpler", "Clade", "Taxon", "TrophicLevel", "BodyMass_g", "InitDensity_Nha")
IndexMissingDens <- which((is.na(SpDensInitSimu$InitDensity_Nha)) & (!is.na(SpDensInitSimu$BodyMass_g)))
SpDensInitSimu$InitDensity_Nha[IndexMissingDens] <- SpDensSeason$Mean_density[IndexMissingDens]
# for ungulates, take their winter densities
SpDensInitSimu$InitDensity_Nha[SpDensInitSimu$Clade == "Ungulate"] <- SpDensSeason$Winter_density[SpDensInitSimu$Clade == "Ungulate"]

write.csv(SpDensInitSimu, file = "SpDensBiomass.csv", row.names = FALSE, quote = FALSE)
