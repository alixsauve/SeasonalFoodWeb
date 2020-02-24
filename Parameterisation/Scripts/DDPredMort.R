# In this script, we draw densitiy-dependent mortality rates of predatory species in BPF.

rm(list=ls(all=TRUE))

library(RColorBrewer)

ALPHA <- 0.75 # scaling parameter (3/4 power law)

# define working directory
DAT_DIR <- "../../Data/BPF_Data"
WK_DIR <- "../../Data/ProcessedData"
PAR_DIR <- "../../Parameterisation/OutputTables"

# load data
setwd(DAT_DIR)
# information on predators
PredInfo <- read.csv("Table5-1.csv", header = TRUE, stringsAsFactors = FALSE, dec = ",")
# species densities
PostBreedingDens <- read.csv("PredPostBreedingBiomass.csv", header = TRUE, stringsAsFactors = FALSE)

setwd(WK_DIR)
# highest longevities
LifeExpectPred <- read.csv("PredMortality.csv", header = TRUE, stringsAsFactors = FALSE, dec = ",")
LifeExpectPred <- subset(LifeExpectPred, TrophicLevel == "predator")
LifeExpectPred$BaseMort <- 1/LifeExpectPred$LifeExpectCapti
LifeExpectPred$BaseMort[is.na(LifeExpectPred$BaseMort)] <- 1/LifeExpectPred$LifeExpectWild[is.na(LifeExpectPred$BaseMort)]
# species densities
SpQuantities <- read.csv("SpDensitySeasonsFinal.csv", header = TRUE, stringsAsFactors = FALSE)
PreyList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Prey"])
PredList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Predator"])
# verify order of species
SpQuantities <- SpQuantities[match(c(PredList, PreyList), SpQuantities$Taxon),]
PredQuantities <- subset(SpQuantities, Trophic_level == "Predator")
PredQuantities$Mean_biomass <- PredQuantities$Mean_density*PredQuantities$Body_mass
PredQuantities$unit_Mean_biomass <- "g/ha"
# birth rates
BirthRatesPred <- read.csv("PredatorsBirthRates.csv", header = TRUE, stringsAsFactors = FALSE)

# load population parameter table
setwd(PAR_DIR)
PopParam <- read.csv("YearPopParam.csv", header = TRUE, stringsAsFactors = FALSE)

# maximum intake per individual predator
BirthRatesPred$MaxPerCapIntake <- 365*(PredInfo$DailyFoodIntake + BirthRatesPred$B*PredInfo$DailyFoodIntake_Juv)

# indicate whether species is considered a constant resource
PopParam$IsCst <- 0
PopParam$IsCst[(PopParam$Clade == "Insect/Arachnid") | (PopParam$Clade == "Plant") |
               (PopParam$Clade == "Mollusk/Lumbricidae") |
               (PopParam$Taxon == "Canis_familiaris") | (PopParam$Taxon == "Felis_catus") | (PopParam$Taxon == "Gallus_gallus") |
               (PopParam$Clade == "Other")] <- 1

# simplify clades
LifeExpectPred$Clade[LifeExpectPred$Clade != "Raptor"] <- "Mammal"
NPred <- nrow(LifeExpectPred)

SpGroup <- unique(LifeExpectPred$Clade)
SpGroupCol <- brewer.pal(length(SpGroup), "Set1"); names(SpGroupCol) <- SpGroup
SpCol <- rep(NA, NPred)
for (Sp in 1:NPred){
  SpCol[Sp] <- SpGroupCol[LifeExpectPred$Clade[Sp]]
}
SpPch <- rep(21, NPred); SpPch[PopParam$Clade[PopParam$TrophicLevel == "Predator"] == "Mustelidae"] <- 23

# DD = (max(b_i) - m_i)/(max(C_i) * 1.5)
BirthRatesPred$M <- PopParam$M[match(BirthRatesPred$Taxon, PopParam$Taxon)]
BirthRatesPred$BodyMass <- PopParam$BodyMass[match(BirthRatesPred$Taxon, PopParam$Taxon)]
BirthRatesPred$DD <- (0.1*BirthRatesPred$MaxPerCapIntake/BirthRatesPred$BodyMass - BirthRatesPred$M)/(PostBreedingDens$PostBreedingBiomass_gha[match(BirthRatesPred$Taxon, PostBreedingDens$Taxon)] * 1.5)
BirthRatesPred$DD <- round(BirthRatesPred$DD, 2)

Index_Mam <- which(LifeExpectPred$Clade != "Raptor")
Index_Rapt <- which(LifeExpectPred$Clade == "Raptor")

LM_DDvsBM_Mam <- lm(log(BirthRatesPred$DD[Index_Mam]) ~ log(BirthRatesPred$BodyMass[Index_Mam]))
CoefLM_Mam <- coef(LM_DDvsBM_Mam)
LM_DDvsBM_Rapt <- lm(log(BirthRatesPred$DD[Index_Rapt]) ~ log(BirthRatesPred$BodyMass[Index_Rapt]))
CoefLM_Rapt <- coef(LM_DDvsBM_Rapt)

par(mar = c(5, 4, 1, 1) + 0.1)
plot(log(BirthRatesPred$BodyMass), log(BirthRatesPred$DD),
     xaxt = "n", yaxt = "n", ylim = log(10^c(-2, 3)),
     bg = SpCol, pch = SpPch, ylab = "Density dependent mortality", xlab = "Body mass (in g)")
abline(CoefLM_Mam[1], CoefLM_Mam[2], lty = 2, col = SpGroupCol["Mammal"])
abline(CoefLM_Rapt[1], CoefLM_Rapt[2], lty = 2, col = SpGroupCol["Raptor"])
axis(1, at = log(c(50, 100, 200, 500, 1000, 2000, 5000, 1e4, 1e5)), labels = c(50, 100, 200, 500, 1000, 2000, 5000, 1e4, 1e5))
axis(2, at = log(10^seq(-2, 2)), labels = 10^seq(-2, 2))
legend("topright", legend = c("Bird", "Mammal"), pt.bg = SpGroupCol, pch = 21, bty = "n")

PopParam$DDModel[PopParam$TrophicLevel == "Predator"] <- BirthRatesPred$DD

# an alternative would be to infer the intrinsic growth rate based on population increase during the reproduction season
# dC_i/dt = C_i * (b_i - m_i - g_i * C_i) = C_i * (r_i - g_i * C_i)
# g_i = (b_i - m_i)/mean(C_i) = r_i / mean(C_i)
BirthRatesPred$DDData <- BirthRatesPred$R/(SpQuantities$Mean_density[match(BirthRatesPred$Taxon, SpQuantities$Taxon)] * BirthRatesPred$BodyMass)

PopParam$DDData[PopParam$TrophicLevel == "Predator"] <- BirthRatesPred$DDData
PopParam <- subset(PopParam, select = -DD)
write.csv(PopParam, file = "YearPopParam.csv", row.names = FALSE, quote = FALSE)

plot(PopParam$DDModel, PopParam$DDData, log = "xy",
     xlab = expression((g[i])[model]), ylab = expression((g[i])[data]),
     xlim = 10^c(-1, 3), ylim = 10^c(-1, 3), pch = 20)
text(PopParam$DDModel, PopParam$DDData, labels = PopParam$Taxon, cex = 0.8)
abline(0, 1); abline(1, 1, col = "grey")
