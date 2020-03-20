# This script draws phase portraits for the different models.

# This script draws time series of the simulations with inserts of phase portraits (or should I do the other way around?).

rm(list=ls(all=TRUE))

EXTINCT_THRS <- 1e-6

# define directories
PAR_DIR <- "../../Parameterisation/OutputTables"
SIM_DIR <- "../../Simulations/Outputs"
FIG_DIR <- "../../Figures/FigureOutputs"

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

# constant resources
CstIndices <- which(PopParam$IsCst == 1)
NCst <- sum(PopParam$IsCst)

# prey list
PreyIndices <- which(PopParam$TrophicLevel == "Prey"); PreyList <- PopParam$Taxon[PreyIndices]
NPrey <- length(PreyIndices)-sum(PopParam$IsCst[PreyIndices])
# predator list
PredIndices <- which(PopParam$TrophicLevel == "Predator"); PredList <- PopParam$Taxon[PredIndices]
NPred <- length(PredIndices)-sum(PopParam$IsCst[PredIndices])

# types of models
GiTypes <- c("GiModel", "GiData")
FRTypes <- c("TSmI", "TSmII")

# time series treatment
Persistence <- matrix(0, 2, 2, dimnames = list(FRTypes, GiTypes))
for (FR in FRTypes){
  for (Gi in GiTypes){
    TS <- get(paste(FR, Gi, sep = "_"))
    TimeVect <- TS[, 1]; TS <- TS[, -1]
    TS[TS < EXTINCT_THRS] <- 0
    TS[, is.na(TS[1, ])] <- NA
    assign(paste(FR, Gi, sep = "_"), TS)
    
    # sum densities within the same trophic level
    TS_Prey <- apply(TS[, PreyIndices], 1, sum, na.rm = TRUE); TS_Pred <- apply(TS[, PredIndices], 1, sum, na.rm = TRUE)
    TS_PredPrey <- cbind(TS_Pred, TS_Prey)
    assign(paste(FR, Gi, "PredPrey", sep = "_"), TS_PredPrey)
    
    # extinction times
    ExtinctTimeIndices <- apply(TS, 2, function(X){ExtTime <- which(X < 1e-6); return(ifelse(length(ExtTime) == 0, NA, min(ExtTime)))})
    assign(paste(FR, Gi, "ExtTimeInd", sep = "_"), ExtinctTimeIndices)
    
    # calculate persistence
    Persistence[FR, Gi] <- length(which((tail(TS, 1) > 0) & (!is.na(tail(TS, 1)))))/(NPrey + NPred - NCst) 
  }
}

# last 10 years
par(mfrow = c(1, 2)); SubPlot <- 1
Index10Y <- which(TimeVect >= 90)
for (FR in FRTypes){
  plot(NA, NA, log = "xy",
       ylim = 10^c(1, 3), xlim = c(1e4, 10^5),
       xlab = "Prey density", ylab = "Predator density", main = paste("Type ", ifelse(FR == FRTypes[1], "I", "II"), sep = ""))
  legend("topleft", legend = toupper(letters[SubPlot]), bty = "n"); SubPlot <- SubPlot + 1
  for (Gi in GiTypes){
    TS_PredPrey <- get(paste(FR, Gi, "PredPrey", sep = "_"))
    lines(TS_PredPrey[Index10Y, 2], TS_PredPrey[Index10Y, 1], col = ifelse(Gi == "DDModel", "black", "grey"))
  }
}
legend("topright", legend = c(expression((g[i])[Model]), expression((g[i])[Data])), col = c("black", "grey"), lty = 1)

# with transients
par(mfrow = c(2, 2)); SubPlot <- 1
for (FR in FRTypes){
  for (Gi in GiTypes){
    TS <- get(paste(FR, Gi, sep = "_"))
    plot(TS_PredPrey[, 2], TS_PredPrey[, 1], type = "l", col = "grey",
         ylim = c(20, 500), xlim = c(1.5e4, 10^5), log = "xy",
         xlab = "Prey density", ylab = "Predator density")
    lines(TS_PredPrey[Index10Y, 2], TS_PredPrey[Index10Y, 1])
    ExtTimeInd <- get(paste(FR, Gi, "ExtTimeInd", sep = "_"))
    points(TS_PredPrey[ExtTimeInd, 2], TS_PredPrey[ExtTimeInd, 1], pch = 3, cex = 0.5)
    legend("topleft", legend = toupper(letters[SubPlot]), bty = "n"); SubPlot <- SubPlot + 1
  }
}
legend("topright", legend = c(expression('T'<'90 years'), expression('T'>='90 years'), "Species extinction"),
       bty = "n", col = c("grey", "black", "black"), lty = c(1, 1, NA), pch = c(NA, NA, 3))

setwd(FIG_DIR)
pdf("Figure_4.pdf", width = 7, height = 5.25)
# phase portraits and time series
layout(matrix(c(1, 2, 2, 3, 4, 4,
                1, 1, 1, 3, 3, 3,
                5, 6, 6, 7, 8, 8,
                5, 5, 5, 7, 7, 7), 4, 6, byrow = TRUE))
SubPlot <- 1
for (FR in FRTypes){
  for (Gi in GiTypes){
    TS_PredPrey <- get(paste(FR, Gi, "PredPrey", sep = "_"))
    # phase portrait
    par(mar = c(5, 6, 0, 1) + 0.1, cex.lab = 1.2)
    plot(TS_PredPrey[, 2], TS_PredPrey[, 1], type = "l", col = "grey",
         ylim = c(20, 5e4), xlim = c(1.8e4, 10^5), log = "xy",
         xlab = expression(sum(R[k], k)), ylab = "",
         frame.plot = FALSE, yaxt = "n")
    axis(2, at = c(20, 50, 100, 200, 500), labels = c(20, 50, 100, 200, 500))
    mtext(expression(sum(C[i], i)), side = 2, line = 2, cex = 0.8, at = 100)
    lines(TS_PredPrey[Index10Y, 2], TS_PredPrey[Index10Y, 1])
    ExtTimeInd <- get(paste(FR, Gi, "ExtTimeInd", sep = "_"))
    points(TS_PredPrey[ExtTimeInd, 2], TS_PredPrey[ExtTimeInd, 1], pch = 3)
    legend("topleft", legend = toupper(letters[SubPlot]), bty = "n", cex = 1.2); SubPlot <- SubPlot + 1
    legend("bottomright", legend = paste("Per. = ", round(Persistence[FR, Gi], 2), sep = ""), cex = 1.2, bty = "n")
    if (Gi == GiTypes[1]){
      mtext(ifelse(FR == "TSmI", "Type I", "Type II"), side = 2, font = 2, cex = 1.2, line = 4.5)
    }
    
    # time series
    TS <- get(paste(FR, Gi, sep = "_"))
    par(mar = c(4, 4, 2.5, 1))
    plot(NA, NA, xlim = c(90, 100), ylim = c(1e-6, 1e6), xlab = "Time", ylab = "Species density", log = "y", frame = FALSE)
    for (Sp in (NPrey+NPred):1){
      ColSp <- ifelse(PopParam$TrophicLevel[Sp] == "Prey", "grey", "black")
      lines(TimeVect, TS[, Sp], col = adjustcolor(ColSp, alpha.f = 0.3))
    }
    if (FR == FRTypes[1]){
      mtext(ifelse(Gi == "GiModel", expression(bold((hat(g)[i])[Model])), expression(bold((hat(g)[i])[Data]))),
            side = 3, cex = 1.2, line = 0.5)
    }
  }
}
dev.off()
