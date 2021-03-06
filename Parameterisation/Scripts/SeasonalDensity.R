rm(list=ls(all=TRUE))

# define directories
DAT_DIR <- "../../Data/ProcessedData"

# load densities
setwd(DAT_DIR)
SpDensities <- read.csv("SpDensitySeasons_1.csv", header = TRUE, stringsAsFactors = FALSE)
PreyList <- sort(SpDensities$Taxon[SpDensities$Trophic_level == "Prey"])
PredList <- sort(SpDensities$Taxon[SpDensities$Trophic_level == "Predator"])
# verify order of species
SpDensities <- SpDensities[match(c(PredList, PreyList), SpDensities$Taxon),]

# convert all densities to N/ha
# spring
IndexUnit2Change <- ((SpDensities$Unit_Spring_density == "N/km2") & (!is.na(SpDensities$Unit_Spring_density)))
if (any(IndexUnit2Change)){
  SpDensities$Spring_density[IndexUnit2Change] <- SpDensities$Spring_density[IndexUnit2Change]/100
  SpDensities$Unit_Spring_density[IndexUnit2Change] <- "N/ha"
}
# summer
IndexUnit2Change <- ((SpDensities$Unit_Summer_density == "N/km2") & (!is.na(SpDensities$Unit_Summer_density)))
if (any(IndexUnit2Change)){
  SpDensities$Summer_density[IndexUnit2Change] <- SpDensities$Summer_density[IndexUnit2Change]/100
  SpDensities$Unit_Summer_density[IndexUnit2Change] <- "N/ha"
}
# autumn
IndexUnit2Change <- ((SpDensities$Unit_Autumn_density == "N/km2") & (!is.na(SpDensities$Unit_Autumn_density)))
if (any(IndexUnit2Change)){
  SpDensities$Autumn_density[IndexUnit2Change] <- SpDensities$Autumn_density[IndexUnit2Change]/100
  SpDensities$Unit_Autumn_density[IndexUnit2Change] <- "N/ha"
}
# winter
IndexUnit2Change <- ((SpDensities$Unit_Winter_density == "N/km2") & (!is.na(SpDensities$Unit_Winter_density)))
if (any(IndexUnit2Change)){
  SpDensities$Winter_density[IndexUnit2Change] <- SpDensities$Winter_density[IndexUnit2Change]/100
  SpDensities$Unit_Winter_density[IndexUnit2Change] <- "N/ha"
}
# year
IndexUnit2Change <- ((SpDensities$Unit_Mean_density == "N/km2") & (!is.na(SpDensities$Unit_Mean_density)))
if (any(IndexUnit2Change)){
  SpDensities$Mean_density[IndexUnit2Change] <- SpDensities$Mean_density[IndexUnit2Change]/100
  SpDensities$Unit_Mean_density[IndexUnit2Change] <- "N/ha"
}

# for our estimates of attack rates, we need winter densities
# here we make a couple of hypotheses to complete this dataset for prey abundances during winter
# for which clades are we missing information?
Clades_MissingWinterAbund <- unique(SpDensities$Clade[(is.na(SpDensities$Winter_density)) & (SpDensities$Trophic_level == "Prey")])
# for amphibians, rodents, shrews, we assume spring densities are pretty much the same as winter's
SpDensities$Winter_density[(is.na(SpDensities$Winter_density)) & (SpDensities$Trophic_level == "Prey")] <- SpDensities$Spring_density[(is.na(SpDensities$Winter_density)) & (SpDensities$Trophic_level == "Prey")]
# remove winter densities for migratory birds
SpDensities$Winter_density[SpDensities$IsMigratoryBird] <- NA
SpDensities$Unit_Winter_density[SpDensities$IsMigratoryBird] <- NA
write.csv(SpDensities, file = "SpDensitySeasons_2.csv", quote = FALSE, row.names = FALSE)

PreyDensities <- subset(SpDensities, Trophic_level == "Prey")

# for birds, mammals (rodents, shrews/moles) and amphibians, we have seasonal densities
DensRange <- 5*10^c(-6, 2)
SeasonPch <- c(8, 21, 25, 16); SeasonBg <- c("black", "white", "white", "black")
SeasonCol <- c("grey", "black", "black", "black")

pdf("SeasonPreyDensity.pdf", width = 5, height = 4)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1) + 0.1, cex.axis = 0.8)
# birds
plot(NA, type = "n", xlim = c(1, 12), ylim = DensRange, log = "y", xaxt = "n",
     xlab = "Month", ylab = "Density (N/ha)")
axis(1, at = seq(1, 12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
IndexBirds <- which(PreyDensities$CladeBis == "Bird")
for (Sp in PreyDensities$Taxon[IndexBirds]){
  IndexSp <- which(PreyDensities$Taxon == Sp)
  DensSp <- c(PreyDensities$Winter_density[IndexSp],
              PreyDensities$Spring_density[IndexSp],
              PreyDensities$Summer_density[IndexSp],
              PreyDensities$Autumn_density[IndexSp])
  lines(c(1, 4, 7, 10), DensSp)
  points(c(1, 4, 7, 10), DensSp, pch = SeasonPch, bg = SeasonBg)
}
mtext("A) Birds", side = 3, at = 4, line = 0.5, cex = 0.6)
# amphibians
plot(NA, type = "n", xlim = c(1, 12), ylim = DensRange, log = "y", xaxt = "n",
     xlab = "Month", ylab = "Density (N/ha)")
axis(1, at = seq(1, 12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
IndexAmphs <- which(PreyDensities$CladeBis == "Amphibian")
for (Sp in PreyDensities$Taxon[IndexAmphs]){
  IndexSp <- which(PreyDensities$Taxon == Sp)
  DensSp <- c(PreyDensities$Winter_density[IndexSp],
              PreyDensities$Spring_density[IndexSp],
              PreyDensities$Summer_density[IndexSp],
              PreyDensities$Autumn_density[IndexSp])
  lines(c(1, 4, 7, 10), DensSp)
  points(c(1, 4, 7, 10), DensSp, pch = SeasonPch, bg = SeasonBg, col = SeasonCol)
}
mtext("B) Amphibians", side = 3, at = 4, line = 0.5, cex = 0.6)
# rodents
plot(NA, type = "n", xlim = c(1, 12), ylim = DensRange, log = "y", xaxt = "n",
     xlab = "Month", ylab = "Density (N/ha)")
axis(1, at = seq(1, 12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
IndexRodents <- which(PreyDensities$Clade == "Rodent")
for (Sp in PreyDensities$Taxon[IndexRodents]){
  IndexSp <- which(PreyDensities$Taxon == Sp)
  DensSp <- c(PreyDensities$Winter_density[IndexSp],
              PreyDensities$Spring_density[IndexSp],
              PreyDensities$Summer_density[IndexSp],
              PreyDensities$Autumn_density[IndexSp])
  lines(c(1, 4, 7, 10), DensSp)
  points(c(1, 4, 7, 10), DensSp, pch = SeasonPch, bg = SeasonBg, col = SeasonCol)
}
mtext("C) Rodents", side = 3, at = 4, line = 0.5, cex = 0.6)
# shrews and moles
plot(NA, type = "n", xlim = c(1, 12), ylim = DensRange, log = "y", xaxt = "n",
     xlab = "Month", ylab = "Density (N/ha)")
axis(1, at = seq(1, 12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
IndexShrews <- which(PreyDensities$Clade == "Shrew/Mole")
for (Sp in PreyDensities$Taxon[IndexShrews]){
  IndexSp <- which(PreyDensities$Taxon == Sp)
  DensSp <- c(PreyDensities$Winter_density[IndexSp],
              PreyDensities$Spring_density[IndexSp],
              PreyDensities$Summer_density[IndexSp],
              PreyDensities$Autumn_density[IndexSp])
  lines(c(1, 4, 7, 10), DensSp)
  points(c(1, 4, 7, 10), DensSp, pch = SeasonPch, bg = SeasonBg, col = SeasonCol)
}
mtext("D) Shrews and moles", side = 3, at = 4, line = 0.5, cex = 0.6)
dev.off()
