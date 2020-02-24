# In this script, I draw a figure representing baseline mortality rates of predators against their

rm(list = ls())

library(RColorBrewer)

# define working directory
WK_DIR <- "../../Data/ProcessedData"
PAR_DIR <- "../../Parameterisation/OutputTables"

setwd(WK_DIR)
LifeExpectPred <- read.csv("PredMortality.csv", header = TRUE, stringsAsFactors = FALSE, dec = ",")
LifeExpectPred <- subset(LifeExpectPred, TrophicLevel == "predator")
LifeExpectPred$BaseMort <- 1/LifeExpectPred$LifeExpectCapti
LifeExpectPred$BaseMort[is.na(LifeExpectPred$BaseMort)] <- 1/LifeExpectPred$LifeExpectWild[is.na(LifeExpectPred$BaseMort)]

# simplify clades
LifeExpectPred$Clade[LifeExpectPred$Clade != "Raptor"] <- "Mammal"
NPred <- nrow(LifeExpectPred)

SpGroup <- unique(LifeExpectPred$Clade)
SpGroupCol <- brewer.pal(length(SpGroup), "Set1"); names(SpGroupCol) <- SpGroup
SpCol <- rep(NA, NPred)
for (Sp in 1:NPred){
  SpCol[Sp] <- SpGroupCol[LifeExpectPred$Clade[Sp]]
}

Index_Mam <- which(LifeExpectPred$Clade != "Raptor")
Index_Rapt <- which(LifeExpectPred$Clade == "Raptor")

LifeExpectPred$Clade <- as.factor(LifeExpectPred$Clade)
AOVMortBMClade <- aov(log(BaseMort) ~ log(BodyMass) * Clade, data = LifeExpectPred)
CoefAOVMortBMClade <- coef(AOVMortBMClade)

par(mar = c(5, 4, 1, 1) + 0.1)
plot(log(LifeExpectPred$BodyMass), log(LifeExpectPred$BaseMort), xlim = log(c(0.05, 35)),
     ylim = log(c(0.025, 0.17)), xaxt = "n", yaxt = "n",
     xlab = "Body mass (in kg)", ylab = "Baseline mortality (/year)", pch = 21, bg = SpCol)
axis(1, at = log(c(0.05, 0.1, 1, seq(0, 40, 5))), labels = c(0.05, 0.1, 1, seq(0, 40, 5)))
axis(2, at = log(seq(0.05, 0.15, 0.05)), labels = seq(0.05, 0.15, 0.05))
abline(CoefAOVMortBMClade[1], CoefAOVMortBMClade[2], lty = 2, col = SpGroupCol["Mammal"])
abline(CoefAOVMortBMClade[1] + CoefAOVMortBMClade[3], CoefAOVMortBMClade[2] + CoefAOVMortBMClade[4], lty = 2, col = SpGroupCol["Raptor"])
legend("topright", legend = c("Bird", "Mammal"), pt.bg = SpGroupCol, pch = 21, bty = "n")

setwd(PAR_DIR)
PopParam <- read.csv("YearPopParam.csv", header  = TRUE, stringsAsFactors = FALSE)
PopParam$M[PopParam$TrophicLevel == "Predator"] <- round(LifeExpectPred$BaseMort, 5)
write.csv(PopParam, file = "YearPopParam.csv", row.names = FALSE, quote = FALSE)