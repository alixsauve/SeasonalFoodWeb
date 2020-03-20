# This scripts draws a three-panel figure with
# A) Predation seasonality
# B) Trophic similarity between seasons
# C) A comparison of proportional generality between seasons

rm(list=ls(all=TRUE))

library(RColorBrewer)
library(bipartite)

E <- 0.1

# define working directories
DAT1_DIR <- "../Data/BPF_Data"
DAT2_DIR <- "../ProcessedData"
FIG_DIR <- "../../Figures/FigureOutputs"

# load functions
source("shannonSeason.R")
source("edgelist2Mat.R")
source("dietBrayCurtis.R")
source("jaccardDistDiet.R")

# load data
setwd(DAT1_DIR)
PredDensPB <- read.csv("PredPostBreedingBiomass.csv", header = TRUE, stringsAsFactors = FALSE)
setwd(DAT2_DIR)
IntakePrey <- read.csv("PerCapitaIntakePreyBiomass.csv", header = TRUE, stringsAsFactors = FALSE)
SpDensSeason <- read.csv("SpDensitySeasonsFinal.csv", header = TRUE, stringsAsFactors = FALSE)
SpQuantities <- read.csv("SpDensBiomass.csv", header = TRUE, stringsAsFactors = FALSE)

# convert predators' densities in N/ha (currently in N/10 km2)
# NB: 1km2 = 100ha -> 10km2 = 1e3ha
PredDensPB$PostBreedingDensity_Nha <- PredDensPB$PostBreedingDensity_N10km2 / 1e3
PredDensPB$MeanDensity_Nha <- PredDensPB$MeanDensity_N10km2 / 1e3

# here we force predators' impact to be 0 if the bird migrates during winter
IntakePrey$IntakePrey_W[SpDensSeason$IsMigratoryBird[match(IntakePrey$LowerTaxon, SpDensSeason$Taxon)]] <- 0
IntakePrey <- IntakePrey[(IntakePrey$IntakePrey_W != 0) | (IntakePrey$IntakePrey_S != 0), ]

NLinks <- nrow(IntakePrey)
# convert per capita biomass intake into population biomass impact times predators' densities -> total biomass harvested
for (Int in 1:NLinks){
  UT <- IntakePrey$UpperTaxon[Int]
  DensUT_S <- PredDensPB$PostBreedingDensity_Nha[PredDensPB$Taxon == UT]
  DensUT_W <- PredDensPB$MeanDensity_Nha[PredDensPB$Taxon == UT]
  IntakePrey$IntakePrey_S[Int] <- IntakePrey$IntakePrey_S[Int] * DensUT_S
  IntakePrey$IntakePrey_W[Int] <- IntakePrey$IntakePrey_W[Int] * DensUT_W
}

# non constant compartments
IndexSp <- which(!is.na(SpQuantities$InitDensity_Nha[SpQuantities$TrophicLevel == "Prey"]))

# build two bipartite matrices for each season
SummerIntakePrey <- subset(IntakePrey, select = c(LowerTaxon, UpperTaxon, IntakePrey_S), (!is.na(IntakePrey_S)) & (IntakePrey_S != 0))
SummerMat <- edgelist2Mat(SummerIntakePrey, "IntakePrey_S")
WinterIntakePrey <- subset(IntakePrey, select = c(LowerTaxon, UpperTaxon, IntakePrey_W), (!is.na(IntakePrey_W)) & (IntakePrey_W != 0))
WinterMat <- edgelist2Mat(WinterIntakePrey, "IntakePrey_W")

# build two bipartite matrices for each season of same size
SummerIntakePrey_Full <- subset(IntakePrey, select = c(LowerTaxon, UpperTaxon, IntakePrey_S))
SummerIntakePrey_Full$IntakePrey_S[is.na(SummerIntakePrey_Full$IntakePrey_S)] <- 0
SummerMat_Full <- edgelist2Mat(SummerIntakePrey_Full, "IntakePrey_S")
WinterIntakePrey_Full <- subset(IntakePrey, select = c(LowerTaxon, UpperTaxon, IntakePrey_W))
WinterIntakePrey_Full$IntakePrey_W[is.na(WinterIntakePrey_Full$IntakePrey_W)] <- 0
WinterMat_Full <- edgelist2Mat(WinterIntakePrey_Full, "IntakePrey_W")

# species biomass densities for each season
PreyBiomassSummer <- (SpDensSeason$Summer_density*SpDensSeason$Body_mass)[match(rownames(SummerMat), SpQuantities$Taxon)]
PreyBiomassWinter <- (SpDensSeason$Winter_density*SpDensSeason$Body_mass)[match(rownames(WinterMat), SpQuantities$Taxon)]
# prey for which seasonal densities are missing
PreyBiomassSummer[is.na(PreyBiomassSummer)] <- (SpQuantities$InitDensity_Nha*SpQuantities$BodyMass_g)[match(rownames(SummerMat)[is.na(PreyBiomassSummer)], SpDensSeason$Taxon)]
PreyBiomassWinter[is.na(PreyBiomassWinter)] <- (SpQuantities$InitDensity_Nha*SpQuantities$BodyMass_g)[match(rownames(WinterMat)[is.na(PreyBiomassWinter)], SpDensSeason$Taxon)]

# predators' densities
PredDensSummer <- PredDensPB$PostBreedingDensity_Nha[match(colnames(SummerMat), PredDensPB$Taxon)]
PredDensWinter <- PredDensPB$MeanDensity_Nha[match(colnames(SummerMat), PredDensPB$Taxon)]
# predators' biomass densities
PredBiomassSummer <- PredDensPB$PostBreedingBiomass_gha[match(colnames(SummerMat), PredDensPB$Taxon)]
PredBiomassWinter <- (SpQuantities$InitDensity_Nha*SpQuantities$BodyMass_g)[match(colnames(WinterMat), SpQuantities$Taxon)]

# network metrics
# species level
SummerMetricsSp <- specieslevel(SummerMat, index = c("partner diversity", "effective partners", "proportional generality"),
                                high.abun = PredBiomassSummer, low.abun = PreyBiomassSummer)
WinterMetricsSp <- specieslevel(WinterMat, index = c("partner diversity", "effective partners", "proportional generality"),
                                high.abun = PredBiomassWinter, low.abun = PreyBiomassSummer)
SpHSeason <- shannonSeason(IntakePrey, "IntakePrey_S", "IntakePrey_W")

IntakeByPred <- subset(SpQuantities, select = c(Taxon, Clade), TrophicLevel == "Predator")
IntakeByPred$Intake_W <- NA; IntakeByPred$Intake_S <- NA
for (Sp in IntakeByPred$Taxon){
  IntakeByPred$Intake_W[IntakeByPred$Taxon == Sp] <- sum(IntakePrey$IntakePrey_W[IntakePrey$UpperTaxon == Sp], na.rm = TRUE)
  IntakeByPred$Intake_S[IntakeByPred$Taxon == Sp] <- sum(IntakePrey$IntakePrey_S[IntakePrey$UpperTaxon == Sp], na.rm = TRUE)
}

# dissimilarity of predators
DissimPartners <- c(); DissimPartners_Bin <- c()
for (Sp in colnames(SummerMat_Full)){
  DissimPartnersSp <- dietBrayCurtis(SummerMat_Full[, Sp], WinterMat_Full[, Sp])
  DissimPartners <- c(DissimPartners, DissimPartnersSp)
  DietSSp <- SummerMat_Full[, Sp]; DietSSp[DietSSp != 0] <- 1
  DietWSp <- WinterMat_Full[, Sp]; DietWSp[DietWSp != 0] <- 1
  DietDissim <- 1
  if ((sum(DietSSp) != 0) & (sum(DietWSp) != 0)){
    DietDissim <- jaccardDistDiet(DietSSp, DietWSp)
  }
  DissimPartners_Bin <- c(DissimPartners_Bin, DietDissim)
}
for (Sp in rownames(SummerMat_Full)){
  DissimPartnersSp <- dietBrayCurtis(SummerMat_Full[Sp, ], WinterMat_Full[Sp, ])
  DissimPartners <- c(DissimPartners, DissimPartnersSp)
  PredSOfSp <- SummerMat_Full[Sp, ]; PredSOfSp[PredSOfSp != 0] <- 1
  PredWOfSp <- WinterMat_Full[Sp, ]; PredWOfSp[PredWOfSp != 0] <- 1
  PredDissim <- 1
  if ((sum(PredSOfSp) != 0) & (sum(PredWOfSp) != 0)){
    PredDissim <- jaccardDistDiet(PredSOfSp, PredWOfSp)
  }
  DissimPartners_Bin <- c(DissimPartners_Bin, PredDissim)
}

# group level
SummerMetricsGr <- grouplevel(SummerMat, weighted = FALSE,
                              index = c("number of species", "mean number of links", "partner diversity", "mean number of shared partners", "niche overlap"))
WinterMetricsGr <- grouplevel(WinterMat, weighted = FALSE,
                              index = c("number of species", "mean number of links", "partner diversity", "mean number of shared partners", "niche overlap"))
RatMetricsGr <- data.frame(S2WRat = SummerMetricsGr/WinterMetricsGr, Metric = names(SummerMetricsGr))
RatMetricsGr$Level <- "Group"
# network level
SummerMetricsNet <- networklevel(SummerMat, index = c("connectance", "links per species"))#, "interaction evenness"))
WinterMetricsNet <- networklevel(WinterMat, index = c("connectance", "links per species"))#, "interaction evenness"))
RatMetricsNet <- data.frame(S2WRat = SummerMetricsNet/WinterMetricsNet, Metric = names(SummerMetricsNet))
RatMetricsNet$Level <- "Network"

SharedPrey <- intersect(rownames(SummerMetricsSp$`lower level`), rownames(WinterMetricsSp$`lower level`))
SharedPred <- intersect(rownames(SummerMetricsSp$`higher level`), rownames(WinterMetricsSp$`higher level`))

setwd(FIG_DIR)
pdf("Figure_2.pdf", width = 7, height = 2.65)

par(mfrow = c(1, 3))
# PLOT 1 (A)
# Season generalism
par(mar = c(5, 5, 1, 1) + 0.1)
plot(sort(1-SpHSeason$LowerGuild$HYear, decreasing = TRUE),
     ylab = "Predation seasonality", xlab = "Species", pch = 20, col = "grey", ylim = c(0, 1))
points(sort(1-SpHSeason$UpperGuild$HYear, decreasing = TRUE), pch = 20)
legend("topleft", legend = toupper(letters[1]), bty = "n")

# PLOT 2 (B)
par(mar = c(5, 5, 1, 1) + 0.1)
plot(sort(1 - DissimPartners[SpQuantities$TrophicLevel == "Prey"], decreasing = TRUE),
     ylab = "Trophic similarity between seasons", xlab = "Species", pch = 20, col = "grey", ylim = c(0, 1))
points(sort(1 - DissimPartners[SpQuantities$TrophicLevel == "Predator"], decreasing = TRUE), pch = 20)
points(sort(1 - DissimPartners_Bin[SpQuantities$TrophicLevel == "Prey"], decreasing = TRUE), col = "grey")
points(sort(1 - DissimPartners_Bin[SpQuantities$TrophicLevel == "Predator"], decreasing = TRUE))
legend("topleft", legend = toupper(letters[2]), bty = "n")

# PLOT 3 (C)
par(mar = c(5, 5, 1, 1) + 0.1)
plot(SummerMetricsSp$`lower level`[SharedPrey, "proportional.generality"],
     WinterMetricsSp$`lower level`[SharedPrey, "proportional.generality"],
     ylab = expression(e^(H^W)[i] / e^(H[max]^W)), xlab = expression(e^(H^S)[i] / e^(H[max]^S)), pch = 20,
     col = adjustcolor("grey", alpha.f = 0.5), log = "xy",
     xlim = 10^c(-2, 2), ylim = 10^c(-2, 2))
points(SummerMetricsSp$`higher level`[SharedPred, "proportional.generality"],
       WinterMetricsSp$`higher level`[SharedPred, "proportional.generality"],
       col = adjustcolor("black", alpha.f = 0.5), pch = 20); abline(0, 1)
MeanHiSPrey <- mean(SummerMetricsSp$`lower level`[SharedPrey, "proportional.generality"])
MeanHiSPred <- mean(SummerMetricsSp$`higher level`[SharedPred, "proportional.generality"])
MeanHiWPrey <- mean(WinterMetricsSp$`lower level`[SharedPrey, "proportional.generality"])
MeanHiWPred <- mean(WinterMetricsSp$`higher level`[SharedPred, "proportional.generality"])
points(c(MeanHiSPred, MeanHiSPrey), c(MeanHiWPred, MeanHiWPrey), pch = 13, col = c("black", "grey"), cex = 1.5)
legend("topleft", legend = toupper(letters[3]), bty = "n")
dev.off()

# testing differences in eHi between summer and winter
EffPartSum <- c(SummerMetricsSp$`lower level`[SharedPrey, "proportional.generality"], SummerMetricsSp$`higher level`[SharedPred, "proportional.generality"])
EffPartWint <- c(WinterMetricsSp$`lower level`[SharedPrey, "proportional.generality"], WinterMetricsSp$`higher level`[SharedPred, "proportional.generality"])
shapiro.test(EffPartSum-EffPartWint) # non normal distribution of differences
wilcox.test(EffPartSum, EffPartWint, paired = TRUE, alternative = "greater") # paired non-parametric test

# predators
EffPartSum_Pred <- SummerMetricsSp$`higher level`[SharedPred, "proportional.generality"]
EffPartWint_Pred <- WinterMetricsSp$`higher level`[SharedPred, "proportional.generality"]
shapiro.test(EffPartSum_Pred-EffPartWint_Pred) # non normal distribution of differences
wilcox.test(EffPartSum_Pred, EffPartWint_Pred, paired = TRUE, alternative = "greater") # paired non-parametric test

# prey species
EffPartSum_Prey <- SummerMetricsSp$`lower level`[SharedPrey, "proportional.generality"]
EffPartWint_Prey <- WinterMetricsSp$`lower level`[SharedPrey, "proportional.generality"]
shapiro.test(EffPartSum_Prey-EffPartWint_Prey) # non normal distribution of differences
wilcox.test(EffPartSum_Prey, EffPartWint_Prey, paired = TRUE, alternative = "greater") # paired non-parametric test


cbind(c(rownames(SummerMetricsSp$`lower level`[SharedPrey, ]), rownames(SummerMetricsSp$`higher level`[SharedPred, ])),
      c(SummerMetricsSp$`lower level`[SharedPrey, "proportional.generality"]-WinterMetricsSp$`lower level`[SharedPrey, "proportional.generality"],
        SummerMetricsSp$`higher level`[SharedPred, "proportional.generality"]-WinterMetricsSp$`higher level`[SharedPred, "proportional.generality"]))

# saving species level metrics
SpLevelMetrics <- subset(SpQuantities, select = c(Clade, Taxon, TrophicLevel))
SpLevelMetrics$PredSeas <- c(SpHSeason$UpperGuild$HYear, SpHSeason$LowerGuild$HYear)
write.csv(SpLevelMetrics, "SpLevelMetrics_20190910.csv", quote = FALSE, row.names = FALSE)
