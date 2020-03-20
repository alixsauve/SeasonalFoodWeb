# This script draws species SSE and SADs.

rm(list=ls(all=TRUE))

EXTINCT_THRS <- 1e-6

# define directories
PAR_DIR <- "../../Parameterisation/OutputTables"
SIM_DIR <- "../../Simulations/Outputs"
DAT2_DIR <- "../../Data/BPF_Data"
FIG_DIR <- "../../Figures/FigureOutputs"
DAT1_DIR <- "../../Data/ProcessedData"

# load observed densities
setwd(DAT1_DIR)
SpDensities <- read.csv("SpDensitySeasonsFinal.csv", header = TRUE, stringsAsFactors = FALSE)
SpDensities$Spring_density[SpDensities$Trophic_level == "Predator"] <- SpDensities$Mean_density[SpDensities$Trophic_level == "Predator"]
SpDensities$Unit_Spring_density[SpDensities$Trophic_level == "Predator"] <- SpDensities$Unit_Mean_density[SpDensities$Trophic_level == "Predator"]
setwd(DAT2_DIR)
PredDensPB <- read.csv("PredPostBreedingBiomass.csv", header = TRUE, stringsAsFactors = FALSE)

# load parameterisation
setwd(PAR_DIR)
PopParam <- read.csv("YearPopParam.csv", header = TRUE, stringsAsFactors = FALSE)
IntParamI <- read.csv("YearIntParam_TypeI.csv", header = TRUE, stringsAsFactors = FALSE)
IntParamII <- read.csv("YearIntParam_TypeII.csv", header = TRUE, stringsAsFactors = FALSE)

# load time series
setwd(SIM_DIR)
TSmI_GiModel <- read.csv("TS_SmI_DDModel.csv", header = FALSE, stringsAsFactors = FALSE)
TSmII_GiModel <- read.csv("TS_SmII_DDModel.csv", header = FALSE, stringsAsFactors = FALSE)
TSmI_GiData <- read.csv("TS_SmI_DDData.csv", header = FALSE, stringsAsFactors = FALSE)
TSmII_GiData <- read.csv("TS_SmII_DDData.csv", header = FALSE, stringsAsFactors = FALSE)

# prey list
PreyIndices <- which(PopParam$TrophicLevel == "Prey"); PreyList <- PopParam$Taxon[PreyIndices]
NPrey <- length(PreyIndices)-sum(PopParam$IsCst[PreyIndices])
# predator list
PredIndices <- which(PopParam$TrophicLevel == "Predator"); PredList <- PopParam$Taxon[PredIndices]
NPred <- length(PredIndices)-sum(PopParam$IsCst[PredIndices])

# subset SpDensities with the species actually modelled
SpDensities <- SpDensities[match(PopParam$Taxon, SpDensities$Taxon), ]

# types of models
GiTypes <- c("GiModel", "GiData")
FRTypes <- c("TSmI", "TSmII")

# time series treatment
for (FR in FRTypes){
  for (Gi in GiTypes){
    TS <- get(paste(FR, Gi, sep = "_"))
    TimeVect <- TS[, 1]; TS <- as.matrix(TS[, -1])
    TS[TS < EXTINCT_THRS] <- 0
    TS[, is.na(TS[1, ])] <- NA
    assign(paste(FR, Gi, sep = "_"), TS)
  }
}

# convert densities to biomasses
SpBiomass <- SpDensities
Seasons <- c("Spring", "Summer", "Autumn", "Winter")
for (S in Seasons){
  DensColName <- paste(S, "density", sep = "_"); UnitColName <- paste("Unit", S, "density", sep = "_")
  SeasonDens2Convert <- which(SpDensities[, UnitColName] == "N/ha")
  SpBiomass[SeasonDens2Convert, DensColName] <- SpBiomass[SeasonDens2Convert, DensColName] * SpBiomass$Body_mass[SeasonDens2Convert]
  SpBiomass[SeasonDens2Convert, UnitColName] <- "g/ha"
}
SpBiomass$Autumn_density[match(PredDensPB$Taxon, SpBiomass$Taxon)] <- PredDensPB$PostBreedingBiomass_gha
SpBiomass$Unit_Autumn_density[match(PredDensPB$Taxon, SpBiomass$Taxon)] <- "g/ha"

# for some species, we do not have seasonal densities
# we compare then with mean densities
IndicesNonSeasDens <- which((is.na(SpBiomass$Spring_density)) & (is.na(SpBiomass$Summer_density)))
SpBiomass$Spring_density[IndicesNonSeasDens] <- SpBiomass$Mean_density[IndicesNonSeasDens]
SpBiomass$Summer_density[IndicesNonSeasDens] <- SpBiomass$Mean_density[IndicesNonSeasDens]

# species symbol
SpBiomass$SpPch <- NA
PchListPred <- c(17, 16); names(PchListPred) <- c("Bird", "Mammal")
PchListPrey <- c(3, 2, 22, 1, 8, 5); names(PchListPrey) <- unique(SpBiomass$CladeBis[PreyIndices])
SpBiomass$SpPch <- PchListPrey[match(SpBiomass$CladeBis, names(PchListPrey))]
SpBiomass$SpPch[SpBiomass$Trophic_level == "Predator"] <- SpBiomass$SpPch[SpBiomass$Trophic_level == "Predator"] + 15

# species symbol
SpBiomass$SpPch <- NA
PchListPred <- c(17, 16); names(PchListPred) <- c("Bird", "Mammal")
PchListPrey <- c(3, 2, 22, 1, 8, 5); names(PchListPrey) <- unique(SpBiomass$CladeBis[PreyIndices])
SpBiomass$SpPch <- PchListPrey[match(SpBiomass$CladeBis, names(PchListPrey))]
SpBiomass$SpPch[SpBiomass$Trophic_level == "Predator"] <- SpBiomass$SpPch[SpBiomass$Trophic_level == "Predator"] + 15

IndexSpring10Y <- match(89:99, TimeVect); IndexSummer10Y <- match(89.25:99.25, TimeVect)
IndexAutumn10Y <- match(89.5:99.5, TimeVect)
PBSpCol <- rep("black", nrow(SpBiomass)); PBSpCol[SpBiomass$Trophic_level == "Predator"] <- "chocolate3"

setwd(FIG_DIR)
pdf("Figure_7.pdf", width = 7, height = 5.25)
layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE), heights = c(2, 2, 0.8))
par(mar = c(4, 6.5, 2, 1), cex.lab = 1.2); SubPlot <- 1
for (FR in FRTypes){
  SummerDens_PerGi <- c(); SpringDens_PerGi <- c(); AutumnDens_PerGi <- c()
  SpSSE_Summer_PerGi <- c(); SpSSE_Spring_PerGi <- c(); SpSSE_Autumn_PerGi <- c() 
  Extinctions_PerGi <- c()
  for (Gi in GiTypes){
    TS <- get(paste(FR, Gi, sep = "_"))

    # mean spring density
    MeanSpringDens <- apply(TS, 2, function(X){return(mean(X[IndexSpring10Y]))})
    MeanSpringDens[tail(TS, 1) == 0] <- NA; SpringDens_PerGi <- cbind(SpringDens_PerGi, MeanSpringDens)
    SpSSE_Spring <- apply(cbind(SpBiomass$Spring_density, MeanSpringDens), 1, function(X){return(diff(X)^2)})
    SpSSE_Spring_PerGi <- cbind(SpSSE_Spring_PerGi, SpSSE_Spring)
    
    # mean Summer density
    MeanSummerDens <- apply(TS, 2, function(X){return(mean(X[IndexSummer10Y]))})
    MeanSummerDens[tail(TS, 1) == 0] <- NA; SummerDens_PerGi <- cbind(SummerDens_PerGi, MeanSummerDens)
    SpSSE_Summer <- apply(cbind(SpBiomass$Summer_density, MeanSummerDens), 1, function(X){return(diff(X)^2)})
    SpSSE_Summer_PerGi <- cbind(SpSSE_Summer_PerGi, SpSSE_Summer)
    
    # mean Autumn density
    MeanAutumnDens <- apply(TS, 2, function(X){return(mean(X[IndexAutumn10Y]))})
    MeanAutumnDens[tail(TS, 1) == 0] <- NA; AutumnDens_PerGi <- cbind(AutumnDens_PerGi, MeanAutumnDens)
    SpSSE_Autumn <- apply(cbind(SpBiomass$Autumn_density, MeanAutumnDens), 1, function(X){return(diff(X)^2)})
    SpSSE_Autumn_PerGi <- cbind(SpSSE_Autumn_PerGi, SpSSE_Autumn)
  }

  # species SSE
  plot(NA, NA, log = "xy",
       xlab = bquote("Squared error"[(g[i])[Model]]), ylab = bquote("Squared error"[(g[i])[Data]]),
       xlim = 10^c(-6, 9), ylim = 10^c(-6, 9)); abline(0, 1)
  SpPch <- ifelse(PopParam$TrophicLevel == "Predator", 1, 20)
  points(SpSSE_Summer_PerGi[PreyIndices, 1], SpSSE_Summer_PerGi[PreyIndices, 2], pch = SpBiomass$SpPch[PreyIndices], col = adjustcolor("black", alpha.f = 0.5), cex = 1.2)
  points(SpSSE_Spring_PerGi[, 1], SpSSE_Spring_PerGi[, 2], pch = SpBiomass$SpPch, col = adjustcolor("chartreuse4", alpha.f = 0.5), cex = 1.2)
  points(SpSSE_Autumn_PerGi[PredIndices, 1], SpSSE_Autumn_PerGi[PredIndices, 2], pch = SpBiomass$SpPch[PredIndices], col = adjustcolor("chocolate3", alpha.f = 0.5), cex = 1.2)
  
  legend("topleft", legend = toupper(letters[SubPlot]), bty = "n", cex = 1.2); SubPlot <- SubPlot + 1
  mtext(ifelse(FR == "TSmI", "Type I", "Type II"), side = 2, font = 2, line = 5)
  
  ObsRankedSpringPredDens <- sort(SpBiomass$Spring_density[PredIndices], decreasing = TRUE)
  ObsRankedAutumnPredDens <- sort(SpBiomass$Autumn_density[PredIndices], decreasing = TRUE)
  # SAD Predators
  plot(seq(1, NPrey, length.out = NPred), ObsRankedAutumnPredDens*1e10, type = "l", log = "y",
       ylim = 10^c(-6, 14), xlim = c(1, NPrey), yaxt = "n",
       xlab = "Species rank", ylab = "Biomass density (g/ha)", lwd = 10, col = adjustcolor("chocolate3", alpha.f = 0.5))
  lines(seq(1, NPrey, length.out = NPred), ObsRankedSpringPredDens*1e10, lwd = 10, col = adjustcolor("chartreuse4", alpha.f = 0.5))#; points(ObsRankedSpringDens, pch = 22, col = "chartreuse4", bg = "cornflowerblue"); 
  lines(seq(1, NPrey, length.out = NPred), sort(AutumnDens_PerGi[PredIndices, 1]*1e10, decreasing = TRUE), col = "chocolate3")
  points(seq(1, NPrey, length.out = NPred), sort(AutumnDens_PerGi[PredIndices, 1]*1e10, decreasing = TRUE), pch = 22, bg = "white", col = "chocolate3", cex = 1)
  lines(seq(1, NPrey, length.out = NPred), sort(SpringDens_PerGi[PredIndices, 1]*1e10, decreasing = TRUE), col = "chartreuse4")
  points(seq(1, NPrey, length.out = NPred), sort(SpringDens_PerGi[PredIndices, 1]*1e10, decreasing = TRUE), pch = 22, bg = "white", col = "chartreuse4", cex = 1)
  lines(seq(1, NPrey, length.out = NPred), sort(AutumnDens_PerGi[PredIndices, 2]*1e10, decreasing = TRUE), col = "chocolate3")
  points(seq(1, NPrey, length.out = NPred), sort(AutumnDens_PerGi[PredIndices, 2]*1e10, decreasing = TRUE), pch = 22, bg = "grey", col = "chocolate3", cex = 1)
  lines(seq(1, NPrey, length.out = NPred), sort(SpringDens_PerGi[PredIndices, 2]*1e10, decreasing = TRUE), col = "chartreuse4")
  points(seq(1, NPrey, length.out = NPred), sort(SpringDens_PerGi[PredIndices, 2]*1e10, decreasing = TRUE), pch = 22, bg = "grey", col = "chartreuse4", cex = 1)
  
  abline(h = 1e-2*1e10)
  axis(3, at = seq(1, NPrey, length.out = NPred), labels = seq(1, NPred))
  axis(2, at = 10^c(seq(-6, 4, 2)), labels = 10^seq(-6, 4, 2))
  axis(2, at = 10^seq(8, 14, 2), labels = 10^seq(-2, 4, 2))
  
  # SAD prey
  ObsRankedSummerPreyDens <- sort(SpBiomass$Summer_density[PreyIndices], decreasing = TRUE)
  ObsRankedSpringPreyDens <- sort(SpBiomass$Spring_density[PreyIndices], decreasing = TRUE)
  lines(seq(1, NPrey), ObsRankedSummerPreyDens, lwd = 10, col = adjustcolor("black", alpha.f = 0.5))
  lines(seq(1, NPrey), ObsRankedSpringPreyDens, lwd = 10, col = adjustcolor("chartreuse4", alpha.f = 0.5))#; points(ObsRankedSpringDens, pch = 22, col = "chartreuse4", bg = "cornflowerblue"); 
  
  lines(sort(SummerDens_PerGi[PreyIndices, 1], decreasing = TRUE)); points(sort(SummerDens_PerGi[PreyIndices, 1], decreasing = TRUE), pch = 22, bg = "white", cex = 1)
  lines(sort(SummerDens_PerGi[PreyIndices, 2], decreasing = TRUE)); points(sort(SummerDens_PerGi[PreyIndices, 2], decreasing = TRUE), pch = 22, bg = "grey", cex = 1)
  lines(sort(SpringDens_PerGi[PreyIndices, 1], decreasing = TRUE), col = "chartreuse4"); points(sort(SpringDens_PerGi[PreyIndices, 1], decreasing = TRUE), pch = 22, bg = "white", col = "chartreuse4", cex = 1)
  lines(sort(SpringDens_PerGi[PreyIndices, 2], decreasing = TRUE), col = "chartreuse4"); points(sort(SpringDens_PerGi[PreyIndices, 2], decreasing = TRUE), pch = 22, bg = "grey", col = "chartreuse4", cex = 1)

  legend("topleft", legend = toupper(letters[SubPlot]), bty = "n"); SubPlot <- SubPlot + 1
}
par(mar = c(0, 4, 0, 0))
plot.new()
legend("left", legend = c("Spring", "Summer", "Autumn"),
       col = c("chartreuse4", "black", "chocolate3"), pch = 20, bty = "n",
       title = expression(bold("Season")), cex = 1.2)
legend("right", legend = sort(names(PchListPrey)), pch = PchListPrey[order(names(PchListPrey))],
       bty = "n", ncol = 2, cex = 1.2, title = expression(bold("Species type")))
plot.new()
legend("center", legend = c("Observed", expression((hat(g)[i])[Model]), expression((hat(g)[i])[Data]),
                            "Spring", "Summer", "Autumn"),
       pch = c(NA, 22, 22, 22, 22, 22), pt.bg = c("cornflowerblue", "white", "grey", "white", "white", "white"),
       col = c("black", "black", "black", "chartreuse4", "black", "chocolate3"), lty = 1, lwd = c(10, 1, 1, 1, 1, 1), ncol = 2, bty = "n", cex = 1.2)
dev.off()
