# This script draws final versus inital densities of predators and prey for each type of model.
rm(list=ls(all=TRUE))

# loading libraries
library(TeachingDemos)

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
SpBiomass$Summer_density[match(PredDensPB$Taxon, SpBiomass$Taxon)] <- PredDensPB$PostBreedingBiomass_gha
SpBiomass$Unit_Summer_density[match(PredDensPB$Taxon, SpBiomass$Taxon)] <- "g/ha"

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


IndexSpring10Y <- match(89:99, TimeVect); IndexSummer10Y <- match(89.25:99.25, TimeVect)
PBSpCol <- rep("black", nrow(SpBiomass)); PBSpCol[SpBiomass$Trophic_level == "Predator"] <- "chocolate3"
setwd(FIG_DIR)
pdf("Figure_5.pdf", width = 7, height = 5.5)
layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE), heights = c(2, 2, 0.8))
par(mar = c(4, 6, 2.5, 1), cex.lab = 1.2); SubPlot <- 1
for (FR in FRTypes){
  for (Gi in GiTypes){
    TS <- get(paste(FR, Gi, sep = "_"))
    # mean spring density
    MeanSpringDens <- apply(TS, 2, function(X){return(mean(X[IndexSpring10Y]))})
    MeanSpringDens[tail(TS, 1) == 0] <- NA
    # mean Summer density
    MeanSummerDens <- apply(TS, 2, function(X){return(mean(X[IndexSummer10Y]))})
    MeanSummerDens[tail(TS, 1) == 0] <- NA
    
    plot(SpBiomass$Summer_density, MeanSummerDens, log = "xy", col = adjustcolor(PBSpCol, alpha.f = 0.5),
         pch = SpBiomass$SpPch, cex = 1.2, xlim = 10^c(-3, 5), ylim = 10^c(-3, 5),
         xlab = "Observed density (g/ha)", ylab = "Simulated density (g/ha)")
    points(SpBiomass$Spring_density, MeanSpringDens, cex = 1.2,
           pch = SpBiomass$SpPch, col = adjustcolor("chartreuse4", alpha.f = 0.5))
    abline(0, 1); abline(1, 1, col = "grey"); abline(-1, 1, col = "grey")
    legend("topleft", legend = toupper(letters[SubPlot]), bty = "n", cex = 1.2)
    
    if (Gi == "GiModel"){
      mtext(ifelse(FR == "TSmI", "Type I", "Type II"), side = 2, font = 2, line = 4.5)
    }
    if (FR == "TSmI"){
      mtext(ifelse(Gi == "GiModel", expression(bold((hat(g)[i])[Model])), expression(bold((hat(g)[i])[Data]))), side = 3, line = 0.5)
    }
    par(mar = c(0, 0, 0, 0), cex.axis = 0.6, tcl = -0.1, mgp = c(3, 0.1, 0))
    subplot(hist(log10(MeanSummerDens/SpBiomass$Summer_density), main = NULL, yaxt = "n", xlab = NA, ylab = NA, xlim = c(-7, 5), lwd = 0.5,
                 breaks = seq(-6.5, 5, by = 0.5), border = "white", col = "black"),
            x = 1e3-100, y = 1e-1, size = c(1.5, 0.5))
    print(range(log10(abs(MeanSummerDens-SpBiomass$Summer_density)), na.rm = TRUE))
    par(mar = c(4, 6, 2.5, 1), cex.axis = 1, tcl = -0.5, mgp = c(3, 1, 0)); SubPlot <- SubPlot + 1
  }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("Spring", "Summer", "Autumn"),
       col = c("chartreuse4", "black", "chocolate3"), pch = 20, bty = "n",
       title = expression(bold("Season")), cex = 1.2)
plot.new()
legend("center", legend = sort(names(PchListPrey)), pch = PchListPrey[order(names(PchListPrey))],
       bty = "n", ncol = 3, title = expression(bold("Species type")), cex = 1.2)
dev.off()
