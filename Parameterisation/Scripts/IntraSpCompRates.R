# In this script, we estimate the intra-specific competition rate of prey species from various groups (voles and mices, amphibians, shrews, moles, birds, reptiles, ungulates).

rm(list=ls(all=TRUE))

library(RColorBrewer)

ALPHA <- 0.75 # scaling parameter (3/4 power law)

# define working directory
DAT_DIR <- "../../Data/BPF_Data"
WK_DIR <- "../../Data/ProcessedData"
PAR_DIR <- "../../Parameterisation/OutputTables"

# load data
setwd(WK_DIR)
SpQuantities <- read.csv("SpDensitySeasonsFinal.csv", header = TRUE, stringsAsFactors = FALSE)
PreyList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Prey"])
PredList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Predator"])
# verify order of species
SpQuantities <- SpQuantities[match(c(PredList, PreyList), SpQuantities$Taxon),]

setwd(PAR_DIR)
PopParam <- read.csv("YearPopParam.csv", header = TRUE, stringsAsFactors = FALSE)

# species lists
PreyList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Prey"])
PredList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Predator"])

# verify order of species
SpQuantities <- SpQuantities[match(c(PredList, PreyList), SpQuantities$Taxon),]
# now subset SpQuantities to the species list in PopParam
SpQuantities <- SpQuantities[match(PopParam$Taxon, SpQuantities$Taxon),]

# missing maximal densities
# pick the maximal densities observed within a year if no other value is provided
IndexMissingMaxDens <- which((is.na(SpQuantities$Max_density)) & (SpQuantities$Trophic_level == "Prey"))
MaxSeasonDens <- apply(subset(SpQuantities, select = c(Spring_density, Summer_density, Autumn_density, Winter_density)), 1, max, na.rm = TRUE)
MaxSeasonDens[is.infinite(MaxSeasonDens)] <- NA
SpQuantities$Max_density[IndexMissingMaxDens] <- MaxSeasonDens[IndexMissingMaxDens]

# rodents
PerIncrease <- (SpQuantities$Max_density - SpQuantities$Mean_density)[SpQuantities$Clade == "Rodent"] / SpQuantities$Mean_density[SpQuantities$Clade == "Rodent"]
PerIncrease <- mean(PerIncrease, na.rm = TRUE)
SpQuantities$Max_density[(SpQuantities$Clade == "Rodent") & (is.na(SpQuantities$Max_density))] <- SpQuantities$Mean_density[(SpQuantities$Clade == "Rodent") & (is.na(SpQuantities$Max_density))]*(1 + PerIncrease)
# passerines
PerIncrease <- (SpQuantities$Max_density - SpQuantities$Mean_density)[SpQuantities$Clade == "Passerine"] / SpQuantities$Mean_density[SpQuantities$Clade == "Passerine"]
PerIncrease <- mean(PerIncrease, na.rm = TRUE)
SpQuantities$Max_density[(SpQuantities$Clade == "Passerine") & (is.na(SpQuantities$Max_density))] <- SpQuantities$Mean_density[(SpQuantities$Clade == "Passerine") & (is.na(SpQuantities$Max_density))]*(1 + PerIncrease)
# anseriformes
PerIncrease <- (SpQuantities$Max_density - SpQuantities$Mean_density)[SpQuantities$Clade == "Anseriforme"] / SpQuantities$Mean_density[SpQuantities$Clade == "Anseriforme"]
PerIncrease <- mean(PerIncrease, na.rm = TRUE)
SpQuantities$Max_density[(SpQuantities$Clade == "Anseriforme") & (is.na(SpQuantities$Max_density))] <- SpQuantities$Mean_density[(SpQuantities$Clade == "Anseriforme") & (is.na(SpQuantities$Max_density))]*(1 + PerIncrease)
# charadriiformes
PerIncrease <- (SpQuantities$Max_density - SpQuantities$Mean_density)[SpQuantities$Clade == "Charadriiforme"] / SpQuantities$Mean_density[SpQuantities$Clade == "Charadriiforme"]
PerIncrease <- mean(PerIncrease, na.rm = TRUE)
SpQuantities$Max_density[(SpQuantities$Clade == "Charadriiforme") & (is.na(SpQuantities$Max_density))] <- SpQuantities$Mean_density[(SpQuantities$Clade == "Charadriiforme") & (is.na(SpQuantities$Max_density))]*(1 + PerIncrease)
# piciformes
PerIncrease <- (SpQuantities$Max_density - SpQuantities$Mean_density)[SpQuantities$Clade == "Piciforme"] / SpQuantities$Mean_density[SpQuantities$Clade == "Piciforme"]
PerIncrease <- mean(PerIncrease, na.rm = TRUE)
SpQuantities$Max_density[(SpQuantities$Clade == "Piciforme") & (is.na(SpQuantities$Max_density))] <- SpQuantities$Mean_density[(SpQuantities$Clade == "Piciforme") & (is.na(SpQuantities$Max_density))]*(1 + PerIncrease)
# birds
BirdsIndex <- sort(c(grep("forme", SpQuantities$Clade), which(SpQuantities$Clade == "Passerine"), which(SpQuantities$Clade == "Bird")))
PerIncrease <- (SpQuantities$Max_density - SpQuantities$Mean_density)[BirdsIndex] / SpQuantities$Mean_density[BirdsIndex]
PerIncrease <- mean(PerIncrease, na.rm = TRUE)
MissingMaxBirdsIndex <- intersect(which((is.na(SpQuantities$Max_density)) & (SpQuantities$Trophic_level == "Prey")), BirdsIndex)
SpQuantities$Max_density[MissingMaxBirdsIndex] <- SpQuantities$Mean_density[MissingMaxBirdsIndex]*(1 + PerIncrease)
# other taxa
MissingMaxTaxaIndex <- which((is.na(SpQuantities$Max_density)) & (SpQuantities$Trophic_level == "Prey") &
                               ((SpQuantities$Clade != "Canidae") & (SpQuantities$Clade != "Felidae")))
PerIncrease <- (SpQuantities$Max_density - SpQuantities$Mean_density)[SpQuantities$Trophic_level == "Prey"] / SpQuantities$Mean_density[SpQuantities$Trophic_level == "Prey"]
PerIncrease <- mean(PerIncrease, na.rm = TRUE)
SpQuantities$Max_density[MissingMaxTaxaIndex] <- SpQuantities$Mean_density[MissingMaxTaxaIndex]*(1 + PerIncrease)

# calculate maximum biomasses based on maximum densities
SpQuantities$Max_biomass <- SpQuantities$Max_density * SpQuantities$Body_mass
SpQuantities$Unit_Max_biomass <- "g/ha"

# set colour code
# pick a colour palette for the different clades
Clades <- data.frame(CladeIndex = 1:9, CladeNames = c("Voles & mice", "Amphibians", "Shrews", "Birds", "Reptiles", "Fish", "Other mammals", "Crustaceans", "Ungulates"))
Clades$CladeCol <- brewer.pal(nrow(Clades), "Set1")

# correspondences between Clades and clade names in PopParam
CladeCorres <- vector("list", nrow(Clades))
CladeCorres[[1]] <- "Rodent"; CladeCorres[[2]] <- "Amphibian"; CladeCorres[[3]] <- "Shrew/Mole"; CladeCorres[[5]] <- "Reptile"; CladeCorres[[6]] <- "Fish"
CladeCorres[[4]] <- c("Raptor", "Passerine", "Bird", unique(PopParam$Clade[grep("forme", PopParam$Clade)]))
CladeCorres[[7]] <- c("Lagomorpha", "Bat")
CladeCorres[[8]] <- "Crustacean"; CladeCorres[[9]] <- "Ungulate"

# assign a colour to each species
PopParam$SpCol <- NA
for (SpGroup in 1:nrow(Clades)){
  PopParam$SpCol[!is.na(match(PopParam$Clade, CladeCorres[[SpGroup]]))] <- Clades$CladeCol[SpGroup]
}
PopParam$SpCol[PopParam$Taxon == "Talpa_europaea"] <- Clades$CladeCol[Clades$CladeNames == "Other mammals"]

# calculate the intra-specific competition rate of prey as r_k / max(R_k)
PopParam$BETA <- NA
PopParam$BETA[PopParam$TrophicLevel == "Prey"] <- PopParam$R[PopParam$TrophicLevel == "Prey"] / SpQuantities$Max_biomass[SpQuantities$Trophic_level == "Prey"]
PopParam$BETA <- round(PopParam$BETA, 5)

# draw a plot
layout(matrix(c(1, 1, 1, 1, 2), 1, 5))
par(mar = c(5, 4, 1, 1) + 0.1)
plot(log(SpQuantities$Body_mass), log(PopParam$BETA), bg = PopParam$SpCol, pch = 21,
     xlab = "Body mass (g)", ylab = "Intra-specific competition rate (ha/year/g)",
     xaxt = "n", yaxt = "n")
axis(1, at = log(10^seq(0, 6)), labels = 10^seq(0, 6))
axis(2, at = log(10^seq(-4, 10, 2)), labels = 10^seq(-4, 10, 2))
# the legend
par(mar = c(5, 1, 1, 1) + 0.1)
plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA)
legend("topleft", legend = Clades$CladeNames, pch = 21, pt.bg = Clades$CladeCol, bty = "n")

PopParam <- subset(PopParam, select = -SpCol)
write.csv(PopParam, file = "YearPopParam.csv", quote = FALSE, row.names = FALSE)
