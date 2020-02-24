# In this script, I infer the attack rates as defined in a functional response of type II.

rm(list=ls(all=TRUE))

library(RColorBrewer)

E <- 0.1 # conversion efficiency

# define working directories
DAT_DIR <- "../../Data/BPF_Data"
WK_DIR <- "../../Data/ProcessedData"
PAR_DIR <- "../../Parameterisation/OutputTables"

# load data
setwd(WK_DIR)
IntakePrey <- read.csv("PerCapitaIntakePreyBiomass.csv", header = TRUE, stringsAsFactors = FALSE)
BirthRatesPred <- read.csv("PredatorsBirthRates.csv", header = TRUE, stringsAsFactors = FALSE)
PrefDiet <- read.csv("PredPref.csv", header = TRUE, stringsAsFactors = FALSE)
SpDensSeason <- read.csv("SpDensitySeasonsFinal.csv", header = TRUE, stringsAsFactors = FALSE)

setwd(DAT_DIR)
PredInfo <- read.csv("Table5-1.csv", header = TRUE, dec = ",", stringsAsFactors = FALSE)

# simplify clade list for the upper taxa
IntakePrey$UpperClade[IntakePrey$UpperClade == "Raptor"] <- "Bird"
IntakePrey$UpperClade[IntakePrey$UpperClade != "Bird"] <- "Mammal"

# here we force predators' impact to be 0 if the bird migrates during winter
IntakePrey$IntakePrey_W[SpDensSeason$IsMigratoryBird[match(IntakePrey$LowerTaxon, SpDensSeason$Taxon)]] <- 0
IntakePrey <- IntakePrey[(IntakePrey$IntakePrey_W != 0) | (IntakePrey$IntakePrey_S != 0), ]

# # WE NEED TO SUBSET PREFDIET
# PrefDiet$FracBioInDiet_W[SpDensSeason$IsMigratoryBird[match(PrefDiet$LowerTaxon, SpDensSeason$Taxon)]] <- 0
# PrefDiet <- PrefDiet[(PrefDiet$FracBioInDiet_S != 0) | (PrefDiet$FracBioInDiet_W != 0), ]

# species list
PreyList <- sort(SpDensSeason$Taxon[SpDensSeason$Trophic_level == "Prey"])
PredList <- sort(PredInfo$Taxon)
NPred <- nrow(PredInfo)

# calculate handling time
PredInfo$HandlingTime <- NA; PredInfo$HandlingTime_Unit <- "y"
# for (Sp in PredList){
for (Sp in PredList){
  DFI_Sp <- PredInfo$DailyFoodIntake[PredInfo$Taxon == Sp]
  DFI_JuvSp <- PredInfo$DailyFoodIntake_Juv[PredInfo$Taxon == Sp]
  BR_Sp <- BirthRatesPred$B[BirthRatesPred$Taxon == Sp]
  PredInfo$HandlingTime[PredInfo$Taxon == Sp] <- 1/((DFI_Sp + BR_Sp*DFI_JuvSp)*365)
}

FracS <- 168/365
FracW <- 197/365

# calculate total discovery rates
PredInfo$TotDiscoveryRate_S <- NA; PredInfo$TotDiscoveryRate_W <- NA
for (Sp in PredList){
  IntakePrey_Sp <- subset(IntakePrey, UpperTaxon == Sp)
  
  TotIntake_S <- sum(IntakePrey_Sp$IntakePrey_S, na.rm = TRUE)
  TotIntake_W <- sum(IntakePrey_Sp$IntakePrey_W, na.rm = TRUE)
  
  TotDiscoveryRate_Sp_S <- 1/FracS*TotIntake_S / (1 - 1/FracS * TotIntake_S * PredInfo$HandlingTime[PredInfo$Taxon == Sp])
  TotDiscoveryRate_Sp_S <- ifelse(is.na(TotDiscoveryRate_Sp_S), 0, TotDiscoveryRate_Sp_S)
  PredInfo$TotDiscoveryRate_S[PredInfo$Taxon == Sp] <- TotDiscoveryRate_Sp_S
  
  TotDiscoveryRate_Sp_W <- 1/FracW*TotIntake_W / (1 - 1/FracW * TotIntake_W * PredInfo$HandlingTime[PredInfo$Taxon == Sp])
  TotDiscoveryRate_Sp_W <- ifelse(is.na(TotDiscoveryRate_Sp_W), 0, TotDiscoveryRate_Sp_W)
  PredInfo$TotDiscoveryRate_W[PredInfo$Taxon == Sp] <- TotDiscoveryRate_Sp_W
  
  if ((PredInfo$TotDiscoveryRate_S[PredInfo$Taxon == Sp] < 0) | (PredInfo$TotDiscoveryRate_W[PredInfo$Taxon == Sp] < 0)){
    print(paste(Sp, " -> Intake must be below ", 1 / PredInfo$HandlingTime[PredInfo$Taxon == Sp], "; estimates {S, W}: {", TotIntake_S, ", ", TotIntake_W, "}", sep = ""))
  }
}

# calculate attack rate for each interaction
# create and initiate two columns for attack rates
IntakePrey$AttackRate_S <- NA; IntakePrey$AttackRate_W <- NA
for (Sp in PredList){
  # list of prey
  PreyOfSp <- sort(unique(IntakePrey$LowerTaxon[IntakePrey$UpperTaxon == Sp]))
  # their densities
  DensPreyOfSp_S <- SpDensSeason$Summer_density[match(PreyOfSp, SpDensSeason$Taxon)]
  DensPreyOfSp_W <- SpDensSeason$Winter_density[match(PreyOfSp, SpDensSeason$Taxon)]
  # prey for which seasonal densities are missing
  IndexMissingDens <- which((is.na(DensPreyOfSp_S)) & (is.na(DensPreyOfSp_W)))
  DensPreyOfSp_S[IndexMissingDens] <- SpDensSeason$Mean_density[match(PreyOfSp[IndexMissingDens], SpDensSeason$Taxon)]
  DensPreyOfSp_W[IndexMissingDens] <- DensPreyOfSp_S[IndexMissingDens]
  
  # their body masses
  BodyMassPreyOfSp <- SpDensSeason$Body_mass[match(PreyOfSp, SpDensSeason$Taxon)]
  # their biomasses
  BiomassPreyOfSp_S <- DensPreyOfSp_S*BodyMassPreyOfSp
  BiomassPreyOfSp_W <- DensPreyOfSp_W*BodyMassPreyOfSp
  
  # for taxa with fixed densities
  IndexFixedDens <- which((PreyOfSp == "Carcass") | (PreyOfSp == "Hymenoptera") | (PreyOfSp == "Lumbricidae_sp."))
  BiomassPreyOfSp_S[IndexFixedDens] <- 1; BiomassPreyOfSp_W[IndexFixedDens] <- 1
  
  # the total prey intakes by the predator
  TotIntake_S <- PredInfo$TotDiscoveryRate_S[PredInfo$Taxon == Sp]
  TotIntake_W <- PredInfo$TotDiscoveryRate_W[PredInfo$Taxon == Sp]
  
  # dietary preferences of the predator Sp
  PrefDiet_Sp <- subset(PrefDiet, UpperTaxon == Sp)
  PrefDiet_Sp <- subset(PrefDiet_Sp, LowerTaxon != Sp)
  PrefDiet_Sp <- PrefDiet_Sp[order(PrefDiet_Sp$LowerTaxon),]
  PrefDiet_Sp <- PrefDiet_Sp[match(PreyOfSp, PrefDiet_Sp$LowerTaxon), ] # keep those under focus
  Pref_S <- PrefDiet_Sp$FracBioInDiet_S; Pref_W <- PrefDiet_Sp$FracBioInDiet_W
  
  # the attack rates
  IntakePrey$AttackRate_S[IntakePrey$UpperTaxon == Sp] <- Pref_S * TotIntake_S / (100 * BiomassPreyOfSp_S)
  IntakePrey$AttackRate_W[IntakePrey$UpperTaxon == Sp] <- Pref_W * TotIntake_W / (100 * BiomassPreyOfSp_W)
  if (length(IntakePrey$AttackRate_S[IntakePrey$UpperTaxon == Sp]) != length(Pref_S * TotIntake_S / (100 * BiomassPreyOfSp_S))){
    stop()
  }
}

IntakePrey$AttackRate_S[is.na(IntakePrey$AttackRate_S)] <- 0
IntakePrey$AttackRate_W[is.na(IntakePrey$AttackRate_W)] <- 0

# first complete data frame with species densities to ease data handling while drawing the figures
IntakePrey$LowerTaxonBodyMass <- NA; IntakePrey$UpperTaxonBodyMass <- NA
for (Sp in PredList){
  BodyMassSp <- SpDensSeason$Body_mass[SpDensSeason$Taxon == Sp]
  IntakePrey$UpperTaxonBodyMass[IntakePrey$UpperTaxon == Sp] <- BodyMassSp
  IntakePrey$HandlingTime[IntakePrey$UpperTaxon == Sp] <- PredInfo$HandlingTime[PredInfo$Taxon == Sp]
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

# create a new dataframe for readability
PredParam <- subset(IntakePrey, select = -c(IntakePrey_S, IntakePrey_W))
PredParam$MeanG <- NA; PredParam$Epsilon <- NA; PredParam$SeasonDom <- NA
for (Sp in PredList){
  IndexIntSp <- which(PredParam$UpperTaxon == Sp)
  
  # calculate mean killing rate
  PredParam$MeanG[IndexIntSp] <- (IntakePrey$AttackRate_W[IndexIntSp] + IntakePrey$AttackRate_S[IndexIntSp]) / 2
  PredParam$Epsilon[IndexIntSp] <- abs(IntakePrey$AttackRate_S[IndexIntSp]-IntakePrey$AttackRate_W[IndexIntSp])/(IntakePrey$AttackRate_S[IndexIntSp]+IntakePrey$AttackRate_W[IndexIntSp])
  
  # summer dominant interactions
  IndexIntSSp <- which((PredParam$UpperTaxon == Sp) & (IntakePrey$AttackRate_S >= IntakePrey$AttackRate_W))
  PredParam$SeasonDom[IndexIntSSp] <- "S"
  # for winter dominant interactions
  IndexIntWSp <- which((PredParam$UpperTaxon == Sp) & (IntakePrey$AttackRate_W > IntakePrey$AttackRate_S))
  PredParam$SeasonDom[IndexIntWSp] <- "W"
}

PredParam <- subset(PredParam, select = c(LowerTaxon, UpperTaxon, SeasonDom, MeanG, Epsilon, HandlingTime))
colnames(PredParam) <- c("LowerTaxon", "UpperTaxon", "SeasonDom", "G", "E", "H")
PredParam$E[is.na(PredParam$E)] <- 0

setwd(PAR_DIR)
write.csv(PredParam, file = "YearIntParam_TypeII.csv", row.names = FALSE)

# represent the parameter values ruling predation
UTBodyMassRange <- c(1e1, 1e6)
BodyMassRatioRange <- c(1e-2, 2e4)
AttackRateRange <- c(1e-2, 1e7)
HandlingTimeRange <- c(1e-7, 1e-3)
IntCol <- brewer.pal(2, "Set1")

PredParam$IntCol <- IntCol[2];
PredParam$IntCol[IntakePrey$UpperClade == "Bird"] <- IntCol[1]

PredInfo$SpCol <- IntCol[1]; PredInfo$SpCol[PredInfo$Clade != "Raptor"] <- IntCol[2]

par(mfrow = c(1, 2)); c <- 1
plot(1, 1, log = "xy", type = "n", xlim = UTBodyMassRange, ylim = HandlingTimeRange,
     xlab = "Predator body mass", ylab = "Handling time (y/g)")
points(PredInfo$BodyMass, PredInfo$HandlingTime, pch = 20,
       col = adjustcolor(PredInfo$SpCol, alpha.f = 0.5))
legend("topleft", legend = toupper(letters[c]), bty = "n"); c <- c + 1
legend("topright", legend = c("Raptor", "Mammal predator"), bty = "n", pch = 20, col = IntCol)
plot(1, 1, log = "xy", type = "n", xlim = BodyMassRatioRange, ylim = AttackRateRange,
     xlab = "Body mass ratio", ylab = "Summer discovery rate (ha/g/y)")
points(IntakePrey$PredToPreyBodyMass, PredParam$G, pch = 20,
       col = adjustcolor(PredParam$IntCol, alpha.f = 0.5))
legend("topleft", legend = toupper(letters[c]), bty = "n"); c <- c + 1