# In this script, we draw the food web of the forest of Bialowieza, with the estimated intakes of prey

rm(list = ls())

# load library
library(bipartite)
library(RColorBrewer)

# define working directories
DAT_DIR <- "../BPF_Data"
FIG_DIR <- "../../Figures/FigureOutputs"
WK_DIR <- "../../Data/ProcessedData"

# load data
setwd(WK_DIR)
ImpactOnPrey <- read.csv("PerCapitaIntakePreyBiomass.csv", header = TRUE, stringsAsFactors = FALSE)

# remove interactions with taxa classified as "Other"
OtherRes <- sort(unique(ImpactOnPrey$LowerTaxon[ImpactOnPrey$LowerClade == "Other"]))
ImpactOnPrey <- ImpactOnPrey[-which((ImpactOnPrey$LowerClade == "Other") | (ImpactOnPrey$LowerTaxon == "Astacus_sp.")), ]

SpQuantities <- read.csv("SpDensBiomass.csv", header = TRUE, stringsAsFactors = FALSE)
SpQuantities <- SpQuantities[-which(!is.na(match(SpQuantities$Taxon, OtherRes))),] # remove taxa classified as "Other"
setwd(DAT_DIR)
PredDensPB <- read.csv("PredPostBreedingBiomass.csv", header = TRUE, stringsAsFactors = FALSE)

# convert predators' densities in N/ha (currently in N/10 km2)
# NB: 1km2 = 100ha -> 10km2 = 1e3ha
PredDensPB$PostBreedingDensity_Nha <- PredDensPB$PostBreedingDensity_N10km2 / 1e3
PredDensPB$MeanDensity_Nha <- PredDensPB$MeanDensity_N10km2 / 1e3

PreyList <- sort(unique(ImpactOnPrey$LowerTaxon)); NPrey <- length(PreyList)
PreyDens <- SpQuantities$MeanDensity[SpQuantities$TrophicLevel == "prey"]
PredList <- sort(unique(ImpactOnPrey$UpperTaxon)); NPred <- length(PredList)
PredDens <- SpQuantities$MeanDensity[SpQuantities$TrophicLevel == "predator"]
NLinks <- nrow(ImpactOnPrey)
NSp <- NPrey + NPred

# rename species groups
SpGroup <- unique(data.frame(Taxon = c(ImpactOnPrey$UpperTaxon, ImpactOnPrey$LowerTaxon),
                             Group = c(ImpactOnPrey$UpperClade, ImpactOnPrey$LowerClade), stringsAsFactors = FALSE))
SpGroup <- SpGroup[match(SpQuantities$Taxon, SpGroup$Taxon),]; SpGroup$TrophicLevel <- SpQuantities$TrophicLevel
SpGroup$Group[c(grep("forme", SpGroup$Group), which((SpGroup$Group == "Passerine") | (SpGroup$Group == "Raptor")))] <- "Bird"
SpGroup$Group[(SpGroup$Group != "Bird") & (SpGroup$Group != "Fish") & (SpGroup$Group != "Reptile") & (SpGroup$Group != "Amphibian")] <- "Mammal"


Groups <- unique(SpGroup$Group[!is.na(SpGroup$Group)])
GroupsCol <- brewer.pal(length(Groups), "Dark2"); names(GroupsCol) <- Groups
SpGroup$SpCol <- GroupsCol[as.character(SpGroup$Group)]

# set seasonal interactions to 0 for the season during which they do not happen
ImpactOnPrey$IntakePrey_W[is.na(ImpactOnPrey$IntakePrey_W)] <- 0
ImpactOnPrey$IntakePrey_S[is.na(ImpactOnPrey$IntakePrey_S)] <- 0

# subset food web, create an array for each season
SummerFW <- subset(ImpactOnPrey, select = c(LowerTaxon, UpperTaxon, IntakePrey_S))
MatSummerFW <- matrix(0, NPrey, NPred, dimnames = list(PreyList, PredList))
WinterFW <- subset(ImpactOnPrey, select = c(LowerTaxon, UpperTaxon, IntakePrey_W))
MatWinterFW <- matrix(0, NPrey, NPred, dimnames = list(PreyList, PredList))

# convert per capita biomass intake into population biomass impact times predators' densities -> total biomass harvested
for (Int in 1:NLinks){
  LT <- ImpactOnPrey$LowerTaxon[Int]; UT <- ImpactOnPrey$UpperTaxon[Int]
  DensUT_S <- PredDensPB$PostBreedingDensity_Nha[PredDensPB$Taxon == UT]
  DensUT_W <- PredDensPB$MeanDensity_Nha[PredDensPB$Taxon == UT]
  # MatSummerFW[LT, UT] <- SummerFW$IntakePrey_S[Int] * DensUT_S
  # MatWinterFW[LT, UT] <- WinterFW$IntakePrey_W[Int] * DensUT_W
  MatSummerFW[LT, UT] <- ImpactOnPrey$IntakePrey_S[Int] * DensUT_S
  MatWinterFW[LT, UT] <- ImpactOnPrey$IntakePrey_W[Int] * DensUT_W
}

rownames(MatSummerFW) <- seq(1, NPrey); rownames(MatWinterFW) <- seq(1, NPrey)
colnames(MatSummerFW) <- seq(1, NPred); colnames(MatWinterFW) <- seq(1, NPred)

PredList <- unlist(lapply(strsplit(PredList, "_"), paste, collapse = " "))
PreyList <- unlist(lapply(strsplit(PreyList, "_"), paste, collapse = " "))
# PreyCol <- SpGroup$SpCol[SpQuantities$TrophicLevel == "prey"]
  
ImpactOnPrey_S <- rowSums(MatSummerFW); ImpactOnPrey_W <- rowSums(MatWinterFW)

setwd(FIG_DIR)
pdf("SeasonalFW_BiomassIntake_BPF.pdf", width = 10, height = 8)
plotweb(MatSummerFW, bor.col.interaction = adjustcolor("grey", 0.5), text.low.col = "white", # ifelse(ImpactOnPrey_S < 1e4, "white", "black"),
        method = "normal", low.abun = PreyDens, high.abun = PredDens, empty = FALSE, col.interaction = adjustcolor("grey80", alpha.f = 0.5),
        col.high = SpGroup$SpCol[SpGroup$TrophicLevel == "Predator"], col.low = SpGroup$SpCol[SpGroup$TrophicLevel == "Prey"],
        bor.col.high = SpGroup$SpCol[SpQuantities$TrophicLevel == "Predator"], bor.col.low = SpGroup$SpCol[SpGroup$TrophicLevel == "Prey"],
        y.lim = c(-3, 2), x.lim = c(-0.1, 2.25), arrow = "up.center", high.y = 1.75, low.y = 0.4)
text(-0.1, 2, labels = "Summer")
plotweb(MatWinterFW, bor.col.interaction = adjustcolor("grey", 0.5), text.low.col = "white", # ifelse(ImpactOnPrey_W < 1e4, "white", "black"),
        method = "normal", low.abun = PreyDens, high.abun = PredDens, empty = FALSE, col.interaction = adjustcolor("grey80", alpha.f = 0.5),
        col.high = SpGroup$SpCol[SpQuantities$TrophicLevel == "Predator"], col.low = SpGroup$SpCol[SpQuantities$TrophicLevel == "Prey"],
        bor.col.high = SpGroup$SpCol[SpQuantities$TrophicLevel == "Predator"], bor.col.low = SpGroup$SpCol[SpQuantities$TrophicLevel == "Prey"],
        arrow = "up.center", add = TRUE, high.y = -0.4, low.y = -1.75) 
text(-0.1, 0, labels = "Winter")
legend("bottomleft", legend = Groups, col = GroupsCol, bty = "n", pch = 15, cex = 0.8)
legend("bottom", legend = paste(seq(1:NPred), PredList, sep = " - "), bty = "n", ncol = 3, cex = 0.8)
dev.off()

PredTotIntake_W <- colSums(MatWinterFW); PredTotIntake_S <- colSums(MatSummerFW)
PredTotIntake <- colSums(MatSummerFW + MatWinterFW)
