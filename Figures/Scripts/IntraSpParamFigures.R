# In this script, I group 4 figures describing parameters driving intra-specific processes
# in our model of BPF community dynamics.

# I aim at a 2 by 2 figure, with A) intrinsic growth rate for prey species, B) intra-specific competition rates for prey species,
# C) baseline mortality rates for predator species, D) density-dependent mortality rates for predator species.


rm(list=ls(all=TRUE))

library(RColorBrewer)

ALPHA <- 0.75

# define working directory
DAT1_DIR <- "../../Data/BPF_Data"
DAT2_DIR <- "../../Data/ProcessedData"
PAR_DIR <- "../../Parameterisation/OutputTables"
FIG_DIR <- "../../Figures/FigureOutputs"  

# load data
setwd(DAT1_DIR)
IntrinsicGR_Clades <- read.csv("HypAverageIntrinsicGR.csv", header = TRUE, stringsAsFactors = FALSE)
IntrinsicGR_Clades$MeanR[IntrinsicGR_Clades$Clade == "Fish"] <- IntrinsicGR_Clades$MeanR[IntrinsicGR_Clades$Clade == "Fish"]
setwd(DAT2_DIR)
SpDensities <- read.csv("SpDensitySeasonsFinal.csv", header = TRUE, stringsAsFactors = FALSE)

setwd(PAR_DIR)
SpQuantities <- read.csv("SpDensBiomass.csv", header = TRUE, stringsAsFactors = FALSE)
PopParam_GiModel <- read.csv("YearPopParam.csv", header = TRUE, stringsAsFactors = FALSE)
PopParam_GiModel <- subset(PopParam_GiModel, select = -DDData)
PopParam_GiData <- read.csv("YearPopParam.csv", header = TRUE, stringsAsFactors = FALSE)
PopParam_GiData <- subset(PopParam_GiData, select = -DDModel)
PopParam <- subset(PopParam_GiData, select = c(Taxon, Clade, TrophicLevel, BodyMass, R, BETA, M, IsCst))
PopParam$DDModel <- PopParam_GiModel$DD; PopParam$DDData <- PopParam_GiData$DD
PopParam$CladeBis <- SpDensities$CladeBis[match(PopParam$Taxon, SpDensities$Taxon)]

PopParam$R[PopParam$Clade == "Fish"] <- PopParam$R[PopParam$Clade == "Fish"] * 15

PredIndices <- which(PopParam$TrophicLevel == "Predator"); PreyIndices <- which(PopParam$TrophicLevel == "Prey")

# a list of clades, for each we provide a data frame with taxa and index in PopParam
Clades <- vector("list", nrow(IntrinsicGR_Clades)); names(Clades) <- IntrinsicGR_Clades$Clade
# Voles/Mice
VolesMice <- which((PopParam$Clade == "Rodent") & (PopParam$TrophicLevel == "Prey"))
VolesMice <- data.frame(Taxon = PopParam$Taxon[VolesMice], Index = VolesMice, stringsAsFactors = FALSE)
Clades$VolesMice <- VolesMice
# Amphibians
Amphibians <- which((PopParam$Clade == "Amphibian") & (PopParam$TrophicLevel == "Prey"))
Amphibians <- data.frame(Taxon = PopParam$Taxon[Amphibians], Index = Amphibians, stringsAsFactors = FALSE)
Clades$Amphibians <- Amphibians
# Shrews
Shrews <- which((PopParam$Clade == "Shrew/Mole") & (PopParam$Taxon != "Talpa_europaea") & (PopParam$TrophicLevel == "Prey"))
Shrews <- data.frame(Taxon = PopParam$Taxon[Shrews], Index = Shrews, stringsAsFactors = FALSE)
Clades$Shrews <- Shrews
# Birds
Birds <- c(grep("forme", PopParam$Clade), 
           which(((PopParam$Clade == "Bird") | (PopParam$Clade == "Passerine")) & (PopParam$TrophicLevel == "Prey")))
Birds <- data.frame(Taxon = PopParam$Taxon[Birds], Index = Birds, stringsAsFactors = FALSE)
Clades$Birds <- Birds
# Reptiles
Reptiles <- which((PopParam$Clade == "Reptile") & (PopParam$TrophicLevel == "Prey"))
Reptiles <- data.frame(Taxon = PopParam$Taxon[Reptiles], Index = Reptiles, stringsAsFactors = FALSE)
Clades$Reptiles <- Reptiles
# Fish
Fish <- which((PopParam$Clade == "Fish") & (PopParam$TrophicLevel == "Prey"))
Fish <- data.frame(Taxon = PopParam$Taxon[Fish], Index = Fish, stringsAsFactors = FALSE)
Clades$Fish <- Fish

# pick a colour palette for the different clades
CladeCol <- brewer.pal(nrow(IntrinsicGR_Clades) + 3, "Set1"); names(CladeCol) <- c(as.character(IntrinsicGR_Clades$Clade), "Other mammals", "Crustaceans", "Ungulates")
names(CladeCol)[1] <- "Voles & mice"
# assign a colour to each species
PopParam$SpCol <- NA
for (SpGroup in 1:length(Clades)){
  PopParam$SpCol[Clades[[SpGroup]]$Index] <- CladeCol[SpGroup]
}
PopParam$SpCol[(PopParam$Taxon == "Lepus_europaeus") | (PopParam$Taxon == "Undet._bats") | (PopParam$Taxon == "Talpa_europaea")] <- CladeCol[SpGroup + 1]
PopParam$SpCol[PopParam$Taxon == "Astacus_sp."] <- CladeCol[SpGroup + 2]
PopParam$SpCol[PopParam$Clade == "Ungulate"] <- CladeCol[SpGroup + 3]

PopParam$SpCol[(PopParam$TrophicLevel == "Predator") & (PopParam$Clade == "Bird")] <- "black"
PopParam$SpCol[(PopParam$TrophicLevel == "Predator") & (PopParam$Clade != "Bird")] <- "white"

# species symbol
PopParam$SpPch <- NA
PchListPred <- c(17, 16); names(PchListPred) <- c("Bird", "Mammal")
PchListPrey <- c(3, 2, 22, 1, 8, 5); names(PchListPrey) <- unique(PopParam$CladeBis[PreyIndices])
PopParam$SpPch <- PchListPrey[match(PopParam$CladeBis, names(PchListPrey))]
PopParam$SpPch[PopParam$TrophicLevel == "Predator"] <- PopParam$SpPch[PopParam$TrophicLevel == "Predator"] + 15

setwd(FIG_DIR)
pdf("Figure_3.pdf", width = 7, height = 7.875)
layout(matrix(c(1, 2, 3, 3, 4, 5, 6, 7), nrow = 4, byrow = TRUE), heights = c(2, 2, 2, 0.8))
par(mar = c(4, 6, 2, 1), cex.lab = 1.2); SubPlot <- 1

# prey intrinsic growth rate
plot(PopParam$BodyMass, PopParam$R, pch = PopParam$SpPch,# col = PopParam$SpCol,
     ylim = 10^c(-1, 2), xlim = 10^c(0, 6), log = "xy",
     xlab = "Body mass (g)", ylab = "Intrinsic growth rate (/y)")
legend("topleft", legend = toupper(letters[1]), bty = "n", cex = 1.2)

# now add segments to highlight body mass scaling and reference mean species
for (SpGroup in names(Clades)){
  # we draw AB segments
  RangeBodyMass <- range(PopParam$BodyMass[Clades[[SpGroup]]$Index], na.rm = TRUE)
  XCoord <- RangeBodyMass + c(-1, 1) * 0.5 # coordinates of point A
  RangeIGR <- range(PopParam$R[Clades[[SpGroup]]$Index], na.rm = TRUE)
  YCoord <- rev(RangeIGR) + c(-1, 1) * 0.5 * (ALPHA-1) # coordinates of point B
  segments(XCoord[1], YCoord[1], XCoord[2], YCoord[2], lty = 2, col = CladeCol[which(names(Clades) == SpGroup)])
  points(IntrinsicGR_Clades$MeanBodyMass[IntrinsicGR_Clades$Clade == SpGroup],
         IntrinsicGR_Clades$MeanR[IntrinsicGR_Clades$Clade == SpGroup],
         pch = 24, col = CladeCol[which(names(Clades) == SpGroup)], bg = "black")
}
legend("topright", legend = c("Voles and mice", "Amphibians", "Shrews", "Birds", "Reptiles", "Fish"),
       lty = 2, col = CladeCol, pch = 24, pt.bg = "black", bty = "n")

# time-varying r_k(t)
TimeVect <- seq(0, 1, by = 1e-3)
plot(TimeVect, sin(2*pi*TimeVect) + 1, type = "l", xlab = "Time (y)", ylab = expression(r[k](t)), yaxt = "n", ylim = c(0, 2.5))
axis(2, at = c(0, 1, 2), labels = c(0, expression(bar(r)[k]), expression(2*bar(r)[k])))
rect(0, 0.5, 0.5, 2, border = NA, col = adjustcolor("grey", alpha.f = 0.5)); text(0.25, 0.75, "Breeding period")
rect(0, 0, 0.5, 0.25, border = NA, col = adjustcolor("goldenrod", alpha.f = 0.5)); text(0.25, 0.125, "Summer")
rect(0.5, 0, 1, 0.25, border = NA, col = adjustcolor("darkblue", alpha.f = 0.5)); text(0.75, 0.125, "Winter")
arrows(0, 1.25, 0, 1, angle = 15, length = 0.1, col = "chartreuse4"); text(0, 1.4, "Spring", col = "chartreuse4", font = 3, adj = c(0.1, 0.5))
arrows(1, 1.25, 1, 1, angle = 15, length = 0.1, col = "chartreuse4"); text(1, 1.4, "Spring", col = "chartreuse4", font = 3, adj = c(0.9, 0.5))
arrows(0.25, 2.25, 0.25, 2, angle = 15, length = 0.1); text(0.25, 2.4, "Summer peak", font = 3)
arrows(0.5, 1.25, 0.5, 1, angle = 15, length = 0.1, col = "chocolate3"); text(0.5, 1.4, "Autumn", font = 3, col = "chocolate3")
legend("topleft", legend = toupper(letters[2]), bty = "n", cex = 1.2)

# prey intra-specific competition
IntraSpCompPrey <- subset(PopParam, (TrophicLevel == "Prey") & (Clade != "Fish"))
IntraSpCompPrey <- subset(IntraSpCompPrey,  (!is.na(BodyMass)) & ((!is.na(BETA))))
LogBodyMass <- log10(IntraSpCompPrey$BodyMass); LogBETA <- log10(IntraSpCompPrey$BETA)
LMBetaBM <- lm(LogBETA ~ LogBodyMass)
CoefLMBetaBM<- coef(LMBetaBM)
NewXValues <- data.frame(XValues = sort(LogBodyMass))
PreyBeta <- predict(LMBetaBM, NewXValues, interval = "confidence")

plot(PopParam$BodyMass, PopParam$BETA,
     ylim = 10^c(-6, 5), xlim = 10^c(0, 6), log = "xy", type = "n",
     xlab = "Body mass (g)", ylab = "Intra-specific competition rate\n(ha/g/y)")
# draw the prediction interval
polygon(10^c(sort(LogBodyMass), sort(LogBodyMass, decreasing = TRUE)),
        10^c(PreyBeta[order(LogBodyMass), "lwr"], rev(PreyBeta[order(LogBodyMass), "upr"])),
        col = adjustcolor("grey", alpha.f = 0.7), border = NA)
lines(10^range(unique(LogBodyMass)), (10^coef(LMBetaBM)[1])*(range(unique(10^LogBodyMass)))^coef(LMBetaBM)[2], lwd = 1)
points(PopParam$BodyMass, PopParam$BETA, xaxt = "n", yaxt = "n", pch = PopParam$SpPch)
legend("topleft", legend = toupper(letters[3]), bty = "n", cex = 1.2)

# predators' baseline mortality
PopParamPred <- subset(PopParam, TrophicLevel == "Predator")
LogBodyMass_Bird <- log10(PopParamPred$BodyMass[PopParamPred$Clade == "Bird"])
LogM_Bird <- log10(PopParamPred$M[PopParamPred$Clade == "Bird"])
LMMBM_Bird <- lm(LogM_Bird ~ LogBodyMass_Bird)
CoefLMMBM_Bird <- coef(LMMBM_Bird)
NewXValues_Bird <- data.frame(XValues = sort(LogBodyMass_Bird))
PredM_Bird <- predict(LMMBM_Bird, NewXValues_Bird, interval = "confidence")

LogBodyMass_Mam <- log10(PopParamPred$BodyMass[PopParamPred$Clade == "Mammal"])
LogM_Mam <- log10(PopParamPred$M[PopParamPred$Clade == "Mammal"])
LMMBM_Mam <- lm(LogM_Mam ~ LogBodyMass_Mam)
CoefLMMBM_Mam <- coef(LMMBM_Mam)
NewXValues_Mam <- data.frame(XValues = sort(LogBodyMass_Mam))
PredM_Mam <- predict(LMMBM_Mam, NewXValues_Mam, interval = "confidence")

plot(PopParam$BodyMass, PopParam$M, log = "xy",
     ylim = 10^c(-2, 0), xlim = 10^c(1, 5), type = "n",
     xlab = "Body mass (g)", ylab = "Baseline mortality")
# draw the prediction interval
polygon(10^c(sort(LogBodyMass_Bird), sort(LogBodyMass_Bird, decreasing = TRUE)),
        10^c(PredM_Bird[order(LogBodyMass_Bird), "lwr"], rev(PredM_Bird[order(LogBodyMass_Bird), "upr"])),
        col = adjustcolor("grey", alpha.f = 0.7), border = NA)
lines(10^range(unique(LogBodyMass_Bird)), (10^coef(LMMBM_Bird)[1])*(range(unique(10^LogBodyMass_Bird)))^coef(LMMBM_Bird)[2])
polygon(10^c(sort(LogBodyMass_Mam), sort(LogBodyMass_Mam, decreasing = TRUE)),
        10^c(PredM_Mam[order(LogBodyMass_Mam), "lwr"], rev(PredM_Mam[order(LogBodyMass_Mam), "upr"])),
        col = adjustcolor("grey", alpha.f = 0.7), border = NA)
lines(10^range(unique(LogBodyMass_Mam)), (10^coef(LMMBM_Mam)[1])*(range(unique(10^LogBodyMass_Mam)))^coef(LMMBM_Mam)[2], lty = 2)

points(PopParam$BodyMass, PopParam$M, pch = PopParam$SpPch)
legend("topleft", legend = toupper(letters[4]), bty = "n", cex = 1.2)

# predators' density-dependent mortality

plot(PopParam$BodyMass, PopParam$DDModel,
     ylim = 10^c(-2, 4), xlim = 10^c(1, 5), log = "xy", type = "n",
     ylab = "Density dependent mortality\n(ha/g/y)", xlab = "Body mass (g)")

LogBodyMass_Bird <- log10(PopParamPred$BodyMass[PopParamPred$Clade == "Bird"])
LogDD_Bird <- log10(PopParamPred$DDData[PopParamPred$Clade == "Bird"])
LMDDBM_Bird <- lm(LogDD_Bird ~ LogBodyMass_Bird)
CoefLMDDBM_Bird <- coef(LMDDBM_Bird)
NewXValues_Bird <- data.frame(XValues = sort(LogBodyMass_Bird))
PredDD_Bird <- predict(LMDDBM_Bird, NewXValues_Bird, interval = "confidence")

LogBodyMass_Mam <- log10(PopParamPred$BodyMass[PopParamPred$Clade == "Mammal"])
LogDD_Mam <- log10(PopParamPred$DDData[PopParamPred$Clade == "Mammal"])
LMDDBM_Mam <- lm(LogDD_Mam ~ LogBodyMass_Mam)
CoefLMDDBM_Mam <- coef(LMDDBM_Mam)
NewXValues_Mam <- data.frame(XValues = sort(LogBodyMass_Mam))
PredDD_Mam <- predict(LMDDBM_Mam, NewXValues_Mam, interval = "confidence")

# draw the prediction interval
polygon(10^c(sort(LogBodyMass_Bird), sort(LogBodyMass_Bird, decreasing = TRUE)),
        10^c(PredDD_Bird[order(LogBodyMass_Bird), "lwr"], rev(PredDD_Bird[order(LogBodyMass_Bird), "upr"])),
        col = adjustcolor("grey47", alpha.f = 0.5), border = NA)
lines(10^range(unique(LogBodyMass_Bird)), (10^coef(LMDDBM_Bird)[1])*range(unique(10^LogBodyMass_Bird))^coef(LMDDBM_Bird)[2], lwd = 2)
polygon(10^c(sort(LogBodyMass_Mam), sort(LogBodyMass_Mam, decreasing = TRUE)),
        10^c(PredDD_Mam[order(LogBodyMass_Mam), "lwr"], rev(PredDD_Mam[order(LogBodyMass_Mam), "upr"])),
        col = adjustcolor("grey47", alpha.f = 0.5), border = NA)
lines(10^range(LogBodyMass_Mam), (10^coef(LMDDBM_Mam)[1])*range(10^LogBodyMass_Mam)^coef(LMDDBM_Mam)[2], lty = 2)


LogBodyMass_Bird <- log10(PopParamPred$BodyMass[PopParamPred$Clade == "Bird"])
LogDD_Bird <- log10(PopParamPred$DDModel[PopParamPred$Clade == "Bird"])
LMDDBM_Bird <- lm(LogDD_Bird ~ LogBodyMass_Bird)
CoefLMDDBM_Bird <- coef(LMDDBM_Bird)
NewXValues_Bird <- data.frame(XValues = sort(LogBodyMass_Bird))
PredDD_Bird <- predict(LMDDBM_Bird, NewXValues_Bird, interval = "confidence")

LogBodyMass_Mam <- log10(PopParamPred$BodyMass[PopParamPred$Clade == "Mammal"])
LogDD_Mam <- log10(PopParamPred$DDModel[PopParamPred$Clade == "Mammal"])
LMDDBM_Mam <- lm(LogDD_Mam ~ LogBodyMass_Mam)
CoefLMDDBM_Mam <- coef(LMDDBM_Mam)
NewXValues_Mam <- data.frame(XValues = sort(LogBodyMass_Mam))
PredDD_Mam <- predict(LMDDBM_Mam, NewXValues_Mam, interval = "confidence")

# draw the prediction interval
polygon(10^c(sort(LogBodyMass_Bird), sort(LogBodyMass_Bird, decreasing = TRUE)),
        10^c(PredDD_Bird[order(LogBodyMass_Bird), "lwr"], rev(PredDD_Bird[order(LogBodyMass_Bird), "upr"])),
        col = adjustcolor("grey87", alpha.f = 0.5), border = NA)
lines(range(10^LogBodyMass_Bird), (10^coef(LMDDBM_Bird)[1])*range(10^LogBodyMass_Bird)^coef(LMDDBM_Bird)[2], col = "white")
polygon(10^c(sort(LogBodyMass_Mam), sort(LogBodyMass_Mam, decreasing = TRUE)),
        10^c(PredDD_Mam[order(LogBodyMass_Mam), "lwr"], rev(PredDD_Mam[order(LogBodyMass_Mam), "upr"])),
        col = adjustcolor("grey87", alpha.f = 0.5), border = NA)
lines(range(10^LogBodyMass_Mam), (10^coef(LMDDBM_Mam)[1])*range(10^LogBodyMass_Mam)^coef(LMDDBM_Mam)[2], lty = 2, col = "white")
# for (Pred in PredIndices){
#   segments(log(PopParam$BodyMass)[Pred], log(PopParam$DDModel)[Pred], log(PopParam$BodyMass)[Pred], log(PopParam$DDData)[Pred])
# }
points(PopParam$BodyMass, PopParam$DDModel, bg = "grey", col = "black", pch = ifelse(PopParam$Clade == "Mammal", 21, 24))
points(PopParam$BodyMass, PopParam$DDData, pch = PopParam$SpPch)
legend("topleft", legend = toupper(letters[5]), bty = "n", cex = 1.2)
legend("topright", legend = c(expression((hat(g)[i])[Model]), expression((hat(g)[i])[Data])), bty = "n",
       pch = c(21, 20), col = "black", pt.bg = c("grey", NA))

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("Predator", "Prey"), pch = c(20, 1), bty = "n",
       title = expression(bold("Trophic position")), cex = 1.2, ncol = 2)
plot.new()
legend("center", legend = sort(names(PchListPrey)), pch = PchListPrey[order(names(PchListPrey))],
       bty = "n", ncol = 3, title = expression(bold("Species type")), cex = 1.2)
dev.off()
