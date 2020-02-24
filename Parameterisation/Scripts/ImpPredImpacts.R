# In this script, we examine the per-capita impact predators have on their prey, as well as the impact prey have on each individuals.
# The aim is to simplify the food web so as to account mostly for the most important interactions, "important" in terms of strength on the population growth.

rm(list = ls())

library(bipartite)

E <- 0.1 # conversion efficiency

DAT_DIR <- "../../Data/BPF_Data"
WK_DIR <- "../../Data/ProcessedData"

# load data
setwd(WK_DIR)
IntakePrey <- read.csv("PerCapitaIntakePreyBiomass.csv", header = TRUE, stringsAsFactors = FALSE)
SpQuantities <- read.csv("SpDensitySeasonsFinal.csv", header = TRUE, stringsAsFactors = FALSE)

# here we force predators' impact to be 0 if the bird migrates during winter
IntakePrey$IntakePrey_W[SpQuantities$IsMigratoryBird[match(IntakePrey$LowerTaxon, SpQuantities$Taxon)]] <- 0
IntakePrey <- IntakePrey[(IntakePrey$IntakePrey_W != 0) | (IntakePrey$IntakePrey_S != 0), ]

# simplify clade list for the upper taxa
IntakePrey$UpperClade[IntakePrey$UpperClade == "Raptor"] <- "Bird"
IntakePrey$UpperClade[IntakePrey$UpperClade != "Bird"] <- "Mammal"
# simplify clade list for the lower taxa
IntakePrey$LowerClade[grep("forme", IntakePrey$LowerClade)] <- "Bird"
IntakePrey$LowerClade[IntakePrey$LowerClade == "Passerine"] <- "Bird"
IntakePrey$LowerClade[!is.na(match(IntakePrey$LowerClade, c("Rodent", "Lagomorpha", "Bat", "Shrew/Mole", "Ungulate")))] <- "Mammal"
IntakePrey$LowerClade[is.na(match(IntakePrey$LowerClade, c("Bird", "Mammal", "Reptile", "Amphibian", "Fish")))] <- "Other"

PreyList <- SpQuantities$Taxon[SpQuantities$Trophic_level == "Prey"]
PredList <- SpQuantities$Taxon[SpQuantities$Trophic_level == "Predator"]

IntTypes <- unique(subset(IntakePrey, select = c(LowerClade, UpperClade)))
# set colour code
IntTypes$IntCol[IntTypes$UpperClade == "Bird"] <- "#E41A1C"
IntTypes$IntCol[IntTypes$UpperClade != "Bird"] <- "#377EB8"
# set symbol code
IntTypes$IntPch[IntTypes$LowerClade == "Bird"] <- 15
IntTypes$IntPch[IntTypes$LowerClade == "Mammal"] <- 16
IntTypes$IntPch[IntTypes$LowerClade == "Reptile"] <- 17
IntTypes$IntPch[IntTypes$LowerClade == "Amphibian"] <- 18
IntTypes$IntPch[IntTypes$LowerClade == "Fish"] <- 1
IntTypes$IntPch[IntTypes$LowerClade == "Other"] <- 8 # this should not display because R_k are unknown

# set colour code for each interaction
IntakePrey$IntCol <- NA; IntakePrey$IntPch <- NA
for (Type in 1:nrow(IntTypes)){
  Index_Type <- which((IntakePrey$LowerClade == IntTypes$LowerClade[Type]) & (IntakePrey$UpperClade == IntTypes$UpperClade[Type]))
  IntakePrey$IntCol[Index_Type] <- IntTypes$IntCol[Type]
  IntakePrey$IntPch[Index_Type] <- IntTypes$IntPch[Type]
}

IntakePrey$LowerTaxonDens <- NA; IntakePrey$UpperTaxonDens <- NA
IntakePrey$LowerTaxonBodyMass <- NA; IntakePrey$UpperTaxonBodyMass <- NA
SpQuantities$Degree <- NA
for (Sp in PreyList){
  # DensSp <- SpQuantities$MeanDensity_Nha[SpQuantities$Taxon == Sp]
  BodyMassSp <- SpQuantities$Body_mass[SpQuantities$Taxon == Sp]
  Index_IntSp <- which(IntakePrey$LowerTaxon == Sp)
  # IntakePrey$LowerTaxonDens[Index_IntSp] <- DensSp
  IntakePrey$LowerTaxonBodyMass[Index_IntSp] <- BodyMassSp
  SpQuantities$Degree[SpQuantities$Taxon == Sp] <- length(Index_IntSp)
}
for (Sp in PredList){
  # DensSp <- SpQuantities$MeanDensity_Nha[SpQuantities$Taxon == Sp]
  BodyMassSp <- SpQuantities$Body_mass[SpQuantities$Taxon == Sp]
  Index_IntSp <- which(IntakePrey$UpperTaxon == Sp)
  # IntakePrey$UpperTaxonDens[Index_IntSp] <- DensSp
  IntakePrey$UpperTaxonBodyMass[Index_IntSp] <- BodyMassSp
  SpQuantities$Degree[SpQuantities$Taxon == Sp] <- length(Index_IntSp)
}

setwd(WK_DIR)
pdf("PerCapitaIntake_BasedOnBiomasses.pdf", width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(5, 4, 1, 1) + 0.1)
# summer
# mammals
plot(1, 1, type = "n", xlim = c(0, max(SpQuantities$Degree, na.rm = TRUE)), ylim = c(1e-4, 1e5),
     pch = 20, log = "y", xlab = "Interactions", ylab = "Per capita intake (N/ha)", main = "Summer")
abline(h = c(1, 6, 24, 168), lty = 2, lwd = 0.5)
text(40, exp(0.5) * c(1, 6, 24, 168), labels = c("1/season", "1/month", "1/week", "1/day"), cex = 0.8, adj = 0)
PredMamList <- unique(IntakePrey$UpperTaxon[IntakePrey$UpperClade == "Mammal"])
for (Sp in PredMamList){
  IntakePreyBySp <- subset(IntakePrey, UpperTaxon == Sp)
  points(sort(IntakePreyBySp$IntakePrey_S/IntakePreyBySp$LowerTaxonBodyMass, decreasing = TRUE), pch = IntakePreyBySp$IntPch, col = IntakePreyBySp$IntCol)
  lines(sort(IntakePreyBySp$IntakePrey_S/IntakePreyBySp$LowerTaxonBodyMass, decreasing = TRUE), col = IntakePreyBySp$IntCol)
  print(min(IntakePreyBySp$IntakePrey_S/IntakePreyBySp$LowerTaxonBodyMass, na.rm = TRUE))
}
legend("topleft", legend = "A", bty = "n") #; legend("topright", legend = "Summer", bty = "n", text.font = 2)
legend("top", legend = c("Bird", "Mammal", "Reptile", "Amphibian", "Fish", "Other"),
       pch = c(15, 16, 17, 18, 1, 8), bty = "n", ncol = 2)

# raptors
plot(1, 1, type = "n", xlim = c(0, max(SpQuantities$Degree, na.rm = TRUE)), ylim = c(1e-4, 1e5),
     pch = 20, log = "y", xlab = "Interactions", ylab = "Per capita intake (N/ha)", main = "Summer")
abline(h = c(1, 6, 24, 168), lty = 2, lwd = 0.5)
text(40, exp(0.5) * c(1, 6, 24, 168), labels = c("1/season", "1/month", "1/week", "1/day"), cex = 0.8, adj = 0)
PredBirList <- unique(IntakePrey$UpperTaxon[IntakePrey$UpperClade == "Bird"])
for (Sp in PredBirList){
  IntakePreyBySp <- subset(IntakePrey, UpperTaxon == Sp)
  points(sort(IntakePreyBySp$IntakePrey_S/IntakePreyBySp$LowerTaxonBodyMass, decreasing = TRUE), pch = IntakePreyBySp$IntPch, col = IntakePreyBySp$IntCol)
  lines(sort(IntakePreyBySp$IntakePrey_S/IntakePreyBySp$LowerTaxonBodyMass, decreasing = TRUE), col = IntakePreyBySp$IntCol)
}
legend("topleft", legend = "B", bty = "n") #; legend("topright", legend = "Summer", bty = "n", text.font = 2)
# winter
# mammals
plot(1, 1, type = "n", xlim = c(0, max(SpQuantities$Degree, na.rm = TRUE)), ylim = c(1e-4, 1e5),
     pch = 20, log = "y", xlab = "Interactions", ylab = "Per capita intake (N/ha)", main = "Winter")
abline(h = c(1, 6, 28, 197), lty = 2, lwd = 0.5)
text(40, exp(0.5) * c(1, 6, 24, 168), labels = c("1/season", "1/month", "1/week", "1/day"), cex = 0.8, adj = 0)
PredMamList <- unique(IntakePrey$UpperTaxon[IntakePrey$UpperClade == "Mammal"])
for (Sp in PredMamList){
  IntakePreyBySp <- subset(IntakePrey, UpperTaxon == Sp)
  points(sort(IntakePreyBySp$IntakePrey_W/IntakePreyBySp$LowerTaxonBodyMass, decreasing = TRUE), pch = IntakePreyBySp$IntPch, col = IntakePreyBySp$IntCol)
  lines(sort(IntakePreyBySp$IntakePrey_W/IntakePreyBySp$LowerTaxonBodyMass, decreasing = TRUE), col = IntakePreyBySp$IntCol)
}
legend("topleft", legend = "C", bty = "n") #; legend("topright", legend = "Winter", bty = "n", text.font = 2)
# raptors
plot(1, 1, type = "n", xlim = c(0, max(SpQuantities$Degree, na.rm = TRUE)), ylim = c(1e-4, 1e5),
     pch = 20, log = "y", xlab = "Interactions", ylab = "Per capita intake (N/ha)", main = "Winter")
abline(h = c(1, 6, 28, 197), lty = 2, lwd = 0.5)
text(40, exp(0.5) * c(1, 6, 24, 168), labels = c("1/season", "1/month", "1/week", "1/day"), cex = 0.8, adj = 0)
PredBirList <- unique(IntakePrey$UpperTaxon[IntakePrey$UpperClade == "Bird"])
for (Sp in PredBirList){
  IntakePreyBySp <- subset(IntakePrey, UpperTaxon == Sp)
  points(sort(IntakePreyBySp$IntakePrey_W/IntakePreyBySp$LowerTaxonBodyMass, decreasing = TRUE), pch = IntakePreyBySp$IntPch, col = IntakePreyBySp$IntCol)
  lines(sort(IntakePreyBySp$IntakePrey_W/IntakePreyBySp$LowerTaxonBodyMass, decreasing = TRUE), col = IntakePreyBySp$IntCol)
}
legend("topleft", legend = "D", bty = "n") #; legend("topright", legend = "Winter", bty = "n", text.font = 2)
dev.off()

# Now, let's try to subset the network
# first, we set a dataframe with the per capita intake of prey
PerCapitaIntakePrey <- subset(IntakePrey, select = c(LowerClade, LowerTaxon, UpperClade, UpperTaxon))
# winter
PerCapitaIntakePrey$PCIntakePrey_W <- IntakePrey$IntakePrey_W/IntakePrey$LowerTaxonBodyMass
ThrsListW <- c(1, 6, 28, 197); names(ThrsListW) <- c("SeasonW", "MonthlyW", "WeeklyW", "DailyW")
for (Thrsd in ThrsListW){
  PerCapitaIntakePrey[[names(ThrsListW)[ThrsListW == Thrsd]]] <- TRUE
  PerCapitaIntakePrey[[names(ThrsListW)[ThrsListW == Thrsd]]][PerCapitaIntakePrey$PCIntakePrey_W < Thrsd] <- FALSE
}
# summer
PerCapitaIntakePrey$PCIntakePrey_S <- IntakePrey$IntakePrey_S/IntakePrey$LowerTaxonBodyMass
ThrsListS <- c(1, 6, 24, 168); names(ThrsListS) <- c("SeasonS", "MonthlyS", "WeeklyS", "DailyS")
for (Thrsd in ThrsListS){
  PerCapitaIntakePrey[[names(ThrsListS)[ThrsListS == Thrsd]]] <- TRUE
  PerCapitaIntakePrey[[names(ThrsListS)[ThrsListS == Thrsd]]][PerCapitaIntakePrey$PCIntakePrey_S < Thrsd] <- FALSE
}

setwd(DAT_DIR)
SubsetFW <- subset(PerCapitaIntakePrey, select = -c(PCIntakePrey_S, PCIntakePrey_W))
write.csv(SubsetFW, file = "IntListSubsets.csv", quote = FALSE, row.names = FALSE)

# represent the corresponding food webs
setwd(WK_DIR)
pdf("SubsetFoodWeb.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
for (ThrsdIndex in 1:4){
  # subset food web
  PerCapitaIntakePrey_Thrsd <- subset(PerCapitaIntakePrey, (PCIntakePrey_S >= ThrsListS[ThrsdIndex]) | (PCIntakePrey_W >= ThrsListW[ThrsdIndex]))
  PerCapitaIntakePrey_Thrsd$PCIntakePrey_S[PerCapitaIntakePrey_Thrsd$PCIntakePrey_S < ThrsListS[ThrsdIndex]] <- 0
  PerCapitaIntakePrey_Thrsd$PCIntakePrey_W[PerCapitaIntakePrey_Thrsd$PCIntakePrey_W < ThrsListW[ThrsdIndex]] <- 0
  PredList_Thrsd <- sort(unique(PerCapitaIntakePrey_Thrsd$UpperTaxon)); NPred_Thrsd <- length(PredList_Thrsd)
  PreyList_Thrsd <- sort(unique(PerCapitaIntakePrey_Thrsd$LowerTaxon)); NPrey_Thrsd <- length(PreyList_Thrsd)
  NLinks_Thrsd <- nrow(PerCapitaIntakePrey_Thrsd)
  
  # create an array for each season
  SummerFW <- subset(PerCapitaIntakePrey_Thrsd, select = c(LowerTaxon, UpperTaxon, PCIntakePrey_S))
  MatSummerFW <- matrix(0, NPrey_Thrsd, NPred_Thrsd, dimnames = list(PreyList_Thrsd, PredList_Thrsd))
  WinterFW <- subset(PerCapitaIntakePrey_Thrsd, select = c(LowerTaxon, UpperTaxon, PCIntakePrey_W))
  MatWinterFW <- matrix(0, NPrey_Thrsd, NPred_Thrsd, dimnames = list(PreyList_Thrsd, PredList_Thrsd))
  
  # fill the matrices
  for (Int in 1:NLinks_Thrsd){
    LT <- PerCapitaIntakePrey_Thrsd$LowerTaxon[Int]; UT <- PerCapitaIntakePrey_Thrsd$UpperTaxon[Int]
    BodyMassLT <- SpQuantities$Body_mass[SpQuantities$Taxon == LT]
    MatSummerFW[LT, UT] <- SummerFW$PCIntakePrey_S[Int] * BodyMassLT
    MatWinterFW[LT, UT] <- WinterFW$PCIntakePrey_W[Int] * BodyMassLT
  }
  
  PredIndex <- which(!is.na(match(PredList, PredList_Thrsd)))
  colnames(MatSummerFW) <- PredIndex
  plotweb(MatSummerFW, bor.col.interaction = adjustcolor("grey", alpha.f = 0.5), text.low.col = "white",
          method = "normal", col.interaction = adjustcolor("grey80", alpha.f = 0.5),
          y.lim = c(-2, 2), x.lim = c(-0.1, 2), arrow = "up.center", high.y = 1.75, low.y = 0.4)
  text(-0.1, 2.05, labels = toupper(letters[ThrsdIndex]))
  colnames(MatWinterFW) <- PredIndex
  plotweb(MatWinterFW, bor.col.interaction = adjustcolor("grey", 0.5), text.low.col = "white",
          method = "normal", col.interaction = adjustcolor("grey80", alpha.f = 0.5),
          arrow = "up.center", add = TRUE, high.y = -0.4, low.y = -1.75)
  
  print(c(NLinks_Thrsd, NPred_Thrsd, NPrey_Thrsd))
}
dev.off()

