# This script reads tables of predators' preferences and take the seasonal preferences
# for the interactions specified by Saavedra et al. 2016, Ecology.

# sum(as.numeric(RawDietSp$FracBioInDiet_W[RawDietSp$FracBioInDiet_W != "+"]))
# sum(as.numeric(RawDietSp$FracBioInDiet_S[RawDietSp$FracBioInDiet_S != "+"]))

rm(list = ls())

DAT_DIR <- "../../Data/BPF_Data"
BIP_DIR <- "../../Data/SaavedraEtAl_2016_Supplements"
WK_DIR <- "../../Data/ProcessedData"

# load information on each predator
setwd(DAT_DIR)
PredInfo <- read.csv("Table5-1.csv", header = TRUE, dec = ",", stringsAsFactors = FALSE)

# load species list and densities across seasons
setwd(WK_DIR)
SpQuantities <- read.csv("SpDensitySeasons_2.csv", header = TRUE, stringsAsFactors = FALSE)

# retrieve Saavedra et al.'s bipartite networks
setwd(BIP_DIR)
SummerMat <- as.matrix(read.csv("Supplement_A_summer_data.csv", header = TRUE, row.names = 1))
WinterMat <- as.matrix(read.csv("Supplement_B_winter_data.csv", header = TRUE, row.names = 1))

# species list
PredList <- sort(unique(c(colnames(SummerMat), colnames(WinterMat))))
PreyList <- sort(unique(c(rownames(SummerMat), rownames(WinterMat))))
# summer species list
PredSummerList <- colnames(SummerMat); PreySummerList <- rownames(SummerMat)

# winter species list
PredWinterList <- colnames(WinterMat); PreyWinterList <- rownames(WinterMat)
# correct prey winter list with information on migration
PreyWinterList2 <- SpQuantities$Taxon[(SpQuantities$Trophic_level == "Prey") & (!SpQuantities$IsMigratoryBird)]

GenusLevelPrey <- PreyList[grep("_sp.", PreyList)] # whole genera predated by predators
GenusLevelPrey <- GenusLevelPrey[(GenusLevelPrey != "Arachnida_sp.") & (GenusLevelPrey != "Carduelis_spinus") &
                                   (GenusLevelPrey != "Gastropoda_sp.") & (GenusLevelPrey != "Lumbricidae_sp.") &
                                   (GenusLevelPrey != "Sylvia_sp.")]

# Annual food web
Mat <- matrix(0, length(PreyList), length(PredList), dimnames = list(PreyList, PredList))
Mat[PreySummerList, PredSummerList] <- SummerMat
Mat[PreyWinterList, PredWinterList] <- Mat[PreyWinterList, PredWinterList] + WinterMat
Mat[Mat != 0] <- 1

# create and initialise a data frame on interaction preferences
IntListPref <- data.frame(LowerTaxon = character(), UpperTaxon = character(), FracBioInDiet_S = numeric(), FracBioInDiet_W = numeric()) 

# for each predator
setwd(DAT_DIR)
for (IndexSp in 1:length(PredList)){
  Sp <- PredInfo$Taxon[IndexSp] # its name
  PreySp <- names(Mat[Mat[, Sp] != 0, Sp]) # its diet in Saavedra et al.'s diet
  DietSp <- data.frame(LowerTaxon = PreySp, UpperTaxon = Sp, FracBioInDiet_S = NA, FracBioInDiet_W = NA, stringsAsFactors = FALSE) # create and initialise a data frame on this predator's diet
  
  # specify if prey are available during summer and winter
  # summer prey
  IndexPreyS <- match(PreySummerList, DietSp$LowerTaxon); IndexPreyS <- IndexPreyS[!is.na(IndexPreyS)]
  if (length(intersect(Sp, PredWinterList)) == 1){
    DietSp$FracBioInDiet_S[IndexPreyS] <- 0
  }
  # winter prey
  IndexPreyW <- match(PreyWinterList, DietSp$LowerTaxon); IndexPreyW <- IndexPreyW[!is.na(IndexPreyW)]
  if (length(intersect(Sp, PredWinterList)) == 1){
    DietSp$FracBioInDiet_W[IndexPreyW] <- 0
  }
  
  TableName <- paste(PredInfo$TableDiet[IndexSp], ".csv", sep = "") # the name of the table with diet preferences
  print(TableName)
  RawDietSp <- read.csv(TableName, header = TRUE, stringsAsFactors = FALSE, dec = ".") # its observed diet
  
  # add a column to the RawDietSp table to identify which observed prey is used in Saavedra et al.'s bipartite web
  RawDietSp$IsUsed <- FALSE
  
  # get the information ready to be identified
  IndexSharedID <- match(DietSp$LowerTaxon, RawDietSp$Scientific_name)
  if (!is.null(RawDietSp$FracBioInDiet_S)){
    DietSp$FracBioInDiet_S[!is.na(IndexSharedID)] <- RawDietSp$FracBioInDiet_S[IndexSharedID[!is.na(IndexSharedID)]]
  }
  if (!is.null(RawDietSp$FracBioInDiet_W)){
    DietSp$FracBioInDiet_W[!is.na(IndexSharedID)] <- RawDietSp$FracBioInDiet_W[IndexSharedID[!is.na(IndexSharedID)]]
  }
  
  RawDietSp$IsUsed[!is.na(match(RawDietSp$Scientific_name, DietSp$LowerTaxon))] <- TRUE
  
  # are there prey with alternative IDs?
  OtherPrey <- unique(RawDietSp$Comment[!is.na(RawDietSp$Comment)])
  OtherPrey <- intersect(OtherPrey, DietSp$LowerTaxon)
  for (Prey in OtherPrey){
    if (!is.null(RawDietSp$FracBioInDiet_S)){
      FracInS <- RawDietSp$FracBioInDiet_S[grep(Prey, RawDietSp$Comment)]
      FracInS[FracInS == "+"] <- 0.01; FracInS <- as.numeric(FracInS)
      FracInS <- sum(FracInS, na.rm = TRUE)
      if (!is.na(DietSp$FracBioInDiet_S[DietSp$LowerTaxon == Prey])){
        DietSp$FracBioInDiet_S[DietSp$LowerTaxon == Prey] <- as.numeric(DietSp$FracBioInDiet_S[DietSp$LowerTaxon == Prey]) + FracInS
      }
      else {
        DietSp$FracBioInDiet_S[DietSp$LowerTaxon == Prey] <- FracInS
      }
    }
    if (!is.null(RawDietSp$FracBioInDiet_W)){
      FracInW <- RawDietSp$FracBioInDiet_W[grep(Prey, RawDietSp$Comment)]
      FracInW[FracInW == "+"] <- 0.01; FracInW <- as.numeric(FracInW)
      FracInW <- sum(FracInW, na.rm = TRUE)
      if (!is.na(DietSp$FracBioInDiet_W[DietSp$LowerTaxon == Prey])){
        DietSp$FracBioInDiet_W[DietSp$LowerTaxon == Prey] <- as.numeric(DietSp$FracBioInDiet_W[DietSp$LowerTaxon == Prey]) + FracInW
      }
      else {
        DietSp$FracBioInDiet_W[DietSp$LowerTaxon == Prey] <- FracInW
      }
    }
  }
  RawDietSp$IsUsed[!is.na(match(RawDietSp$Comment, OtherPrey))] <- TRUE
  
  # get the names of prey which are not identified in the raw diet
  MissingPrey <- DietSp$LowerTaxon[(is.na(DietSp$FracBioInDiet_S)) & (is.na(DietSp$FracBioInDiet_W))]
  if (length(MissingPrey) > 0){
    GenusMissingPrey <- sapply(strsplit(MissingPrey, split = "_"), "[[", 1)
    MissingPrey <- paste(GenusMissingPrey, "_sp.", sep = "")
    IndexMissingPrey <- which((is.na(DietSp$FracBioInDiet_S)) & (is.na(DietSp$FracBioInDiet_W)))
    for (Prey in MissingPrey){
      if (length(grep(Prey, RawDietSp$Scientific_name)) == 1){
        if (!is.null(RawDietSp$FracBioInDiet_S)){
          FracInS <- RawDietSp$FracBioInDiet_S[grep(Prey, RawDietSp$Scientific_name)]
          FracInS[FracInS == "+"] <- 0.01
          DietSp$FracBioInDiet_S[IndexMissingPrey[MissingPrey == Prey]] <- as.numeric(FracInS)
        }
        if (!is.null(RawDietSp$FracBioInDiet_W)){
          FracInW <- RawDietSp$FracBioInDiet_W[grep(Prey, RawDietSp$Scientific_name)]
          FracInW[FracInW == "+"] <- 0.01
          DietSp$FracBioInDiet_W[IndexMissingPrey[MissingPrey == Prey]] <- as.numeric(FracInW)
        } 
      }
    }
  }
  RawDietSp$IsUsed[!is.na(match(RawDietSp$Scientific_name, MissingPrey))] <- TRUE
  
  DietSp$FracBioInDiet_S[DietSp$FracBioInDiet_S == "+"] <- 0.01; DietSp$FracBioInDiet_S <- as.numeric(DietSp$FracBioInDiet_S)
  DietSp$FracBioInDiet_W[DietSp$FracBioInDiet_W == "+"] <- 0.01; DietSp$FracBioInDiet_W <- as.numeric(DietSp$FracBioInDiet_W)
  
  # now, distribute other diet fractions within the shortened prey list based on species genus
  PreyGenus <- sapply(strsplit(DietSp$LowerTaxon, split = "_"), "[[", 1)
  PreyGenus[PreyGenus == "Undet."] <- NA
  UnusedTaxa <- RawDietSp$Scientific_name[!RawDietSp$IsUsed]
  GenusRaw <- UnusedTaxa[grep("_sp.", UnusedTaxa)]
  for (Genus in GenusRaw){
    IndexCorrPrey <- which(paste(PreyGenus, "_sp.", sep = "") == Genus)
    if (length(IndexCorrPrey) != 0){
      RawDietSp$IsUsed[grep(Genus, RawDietSp$Scientific_name)] <- TRUE
      if (!is.null(RawDietSp$FracBioInDiet_S)){
        IndexCorrPreyS <- intersect(IndexCorrPrey, IndexPreyS)
        FracInS <- RawDietSp$FracBioInDiet_S[grep(Genus, RawDietSp$Scientific_name)]
        FracInS[FracInS == "+"] <- 0.01
        DietSp$FracBioInDiet_S[IndexCorrPreyS] <- as.numeric(DietSp$FracBioInDiet_S[IndexCorrPreyS]) + as.numeric(FracInS)/length(IndexCorrPreyS)
      }
      if (!is.null(RawDietSp$FracBioInDiet_W)){
        IndexCorrPreyW <- intersect(IndexCorrPrey, IndexPreyW)
        FracInW <- RawDietSp$FracBioInDiet_W[grep(Genus, RawDietSp$Scientific_name)]
        FracInW[FracInW == "+"] <- 0.01
        DietSp$FracBioInDiet_W[IndexCorrPreyW] <- as.numeric(DietSp$FracBioInDiet_W[IndexCorrPreyW]) + as.numeric(FracInW)/length(IndexCorrPreyW)
      }
    }
  }
  
  # get genus-level prey in the modelled diet
  GenusLevelPrey_Sp <- intersect(DietSp$LowerTaxon, GenusLevelPrey)
  for (Prey in GenusLevelPrey_Sp){
    Genus <- strsplit(Prey, "_")[[1]][1]
    IndexPreyInRawDiet <- grep(Genus, RawDietSp$Scientific_name)
    RawDietSp$IsUsed[IndexPreyInRawDiet] <- TRUE
    if (!is.null(RawDietSp$FracBioInDiet_S)){
      FracInS <- RawDietSp$FracBioInDiet_S[IndexPreyInRawDiet]
      FracInS[FracInS == "+"] <- 0.01
      DietSp$FracBioInDiet_S[grep(Prey, DietSp$LowerTaxon)] <- sum(as.numeric(FracInS))
    }
    if (!is.null(RawDietSp$FracBioInDiet_W)){
      FracInW <- RawDietSp$FracBioInDiet_W[IndexPreyInRawDiet]
      FracInW[FracInW == "+"] <- 0.01
      DietSp$FracBioInDiet_W[grep(Prey, DietSp$LowerTaxon)] <- sum(as.numeric(FracInW))
    }
  }
  
  # undetermined deer (Cervidae)
  if (length(grep("Undet._deer", RawDietSp$Comment)) != 0){
    CervidaeIndex <- c(grep("Capreolus", DietSp$LowerTaxon), grep("Cervus", DietSp$LowerTaxon))
    CervidaeIndexRawData <- c(grep("Capreolus", RawDietSp$Scientific_name), grep("Cervus", RawDietSp$Scientific_name))
    if (length(CervidaeIndex) != 0){
      if (!is.null(RawDietSp$FracBioInDiet_S)){
        FracInS <- as.numeric(RawDietSp$FracBioInDiet_S[grep("Undet._deer", RawDietSp$Comment)])
        DietSp$FracBioInDiet_S[CervidaeIndex] <- as.numeric(DietSp$FracBioInDiet_S[CervidaeIndex]) + FracInS/length(CervidaeIndexRawData)
      }
      if (!is.null(RawDietSp$FracBioInDiet_W)){
        FracInW <- as.numeric(RawDietSp$FracBioInDiet_W[grep("Undet._deer", RawDietSp$Comment)])
        DietSp$FracBioInDiet_W[CervidaeIndex] <- as.numeric(DietSp$FracBioInDiet_W[CervidaeIndex]) + FracInW/length(CervidaeIndexRawData)
      }
      RawDietSp$IsUsed[RawDietSp$Comment == "Undet._deer"] <- TRUE
    }
  }
  # undetermined shrew
  if (length(grep("Undet._shrews", RawDietSp$Comment)) != 0){
    ShrewIndex <- c(grep("Sorex", DietSp$LowerTaxon), grep("Neomys", DietSp$LowerTaxon))
    ShrewIndexRawData <- c(grep("Sorex", RawDietSp$Scientific_name), grep("Neomys", RawDietSp$Scientific_name))
    ShrewIndexRawData <- ShrewIndexRawData[-grep("_sp.", RawDietSp$Scientific_name[ShrewIndexRawData])]
    if (length(ShrewIndex) != 0){
      if (!is.null(RawDietSp$FracBioInDiet_S)){
        FracInS <- as.numeric(RawDietSp$FracBioInDiet_S[grep("Undet._shrews", RawDietSp$Comment)])
        DietSp$FracBioInDiet_S[ShrewIndex] <- as.numeric(DietSp$FracBioInDiet_S[ShrewIndex]) + FracInS/length(ShrewIndexRawData)
      }
      if (!is.null(RawDietSp$FracBioInDiet_W)){
        FracInW <- as.numeric(RawDietSp$FracBioInDiet_W[grep("Undet._shrews", RawDietSp$Comment)])
        DietSp$FracBioInDiet_W[ShrewIndex] <- as.numeric(DietSp$FracBioInDiet_W[ShrewIndex]) + FracInW/length(ShrewIndexRawData)
      }
      RawDietSp$IsUsed[RawDietSp$Comment == "Undet._shrews"] <- TRUE
    }
  }
  # undetermined fish
  if (length(grep("Undet._fish", RawDietSp$Comment)) != 0){
    FishList <- c("Abramis_brama", "Blicca_bjoerkna", "Carassius_carassius", "Esox_lucius", "Gasterosteus_aculeatus", "Gobio_gobio", "Gymnocephalus_cernuus", "Leuciscus_idus", "Lota_lota", "Misgurnus_fossilis", "Perca_fluviatilis", "Rutilus_rutilus", "Scardinius_erythrophthalmus", "Tinca_tinca")
    FishIndex <- match(FishList, DietSp$LowerTaxon); FishIndex <- FishIndex[!is.na(FishIndex)]
    if (length(FishIndex) != 0){
      if (!is.null(RawDietSp$FracBioInDiet_S)){
        FishIndexS <- intersect(FishIndex, IndexPreyS)
        FracInS <- as.numeric(RawDietSp$FracBioInDiet_S[grep("Undet._fish", RawDietSp$Common_name)])
        DietSp$FracBioInDiet_S[FishIndexS] <- DietSp$FracBioInDiet_S[FishIndexS] + FracInS/length(FishIndexS)
      }
      if (!is.null(RawDietSp$FracBioInDiet_W)){
        FishIndexW <- intersect(FishIndex, IndexPreyW)
        FracInW <- as.numeric(RawDietSp$FracBioInDiet_W[grep("Undet._fish", RawDietSp$Common_name)])
        DietSp$FracBioInDiet_W[FishIndex] <- DietSp$FracBioInDiet_W[FishIndex] + FracInW/length(FishIndex)
      }
    }
    RawDietSp$IsUsed[grep("Undet._fish", RawDietSp$Common_name)] <- TRUE
  }
  
  # predator-specific treatments
  if (TableName == "Table4-8.csv"){
    if (!is.null(RawDietSp$FracBioInDiet_S)){
      DietSp$FracBioInDiet_S[DietSp$LowerTaxon == "Natrix_natrix"] <- RawDietSp$FracBioInDiet_S[RawDietSp$Common_name == "Undet._reptile"]
    }
    if (!is.null(RawDietSp$FracBioInDiet_W)){
      DietSp$FracBioInDiet_W[DietSp$LowerTaxon == "Natrix_natrix"] <- RawDietSp$FracBioInDiet_W[RawDietSp$Common_name == "Undet._reptile"]
    }
    RawDietSp$IsUsed[RawDietSp$Common_name == "Undet._reptile"] <- TRUE
  }
  if (TableName == "Table4-40.csv"){
    if (!is.null(RawDietSp$FracBioInDiet_S)){
      # fish
      DietSp$FracBioInDiet_S[DietSp$LowerTaxon == "Perca_fluviatilis"] <- RawDietSp$FracBioInDiet_S[RawDietSp$Common_name == "Fish"]
      # insects
      IndexInsects <- which((DietSp$LowerTaxon == "Heteroptera") | (DietSp$LowerTaxon == "Hymenoptera") | (DietSp$LowerTaxon == "Lepidoptera"))
      FracInS <- RawDietSp$FracBioInDiet_S[RawDietSp$Common_name == "Other_and_Undet._insects"];
      FracInS[FracInS == "+"] <- 0.01; FracInS <- as.numeric(FracInS)
      DietSp$FracBioInDiet_S[IndexInsects] <- FracInS / length(IndexInsects)
    }
    if (!is.null(RawDietSp$FracBioInDiet_W)){
      # fish
      DietSp$FracBioInDiet_W[DietSp$LowerTaxon == "Perca_fluviatilis"] <- RawDietSp$FracBioInDiet_W[RawDietSp$Common_name == "Fish"]
      # insects
      IndexInsects <- which((DietSp$LowerTaxon == "Heteroptera") | (DietSp$LowerTaxon == "Hymenoptera") | (DietSp$LowerTaxon == "Lepidoptera"))
      FracInW <- RawDietSp$FracBioInDiet_W[RawDietSp$Common_name == "Other_and_Undet._insects"];
      FracInW[FracInW == "+"] <- 0.01; FracInW <- as.numeric(FracInW)
      DietSp$FracBioInDiet_W[IndexInsects] <- FracInW / length(IndexInsects)
    }
    RawDietSp$IsUsed[RawDietSp$Common_name == "Fish"] <- TRUE
    RawDietSp$IsUsed[RawDietSp$Common_name == "Other_and_Undet._insects"] <- TRUE
  }
  if (TableName == "Table4-25.csv"){
    # ide or gudgeon
    IndexTaxa <- which((DietSp$LowerTaxon == "Leuciscus_idus") | (DietSp$LowerTaxon == "Gobio_gobio"))
    IndexTaxaS <- intersect(IndexTaxa, IndexPreyS)
    DietSp$FracBioInDiet_S[IndexTaxaS] <- DietSp$FracBioInDiet_S[IndexTaxaS] + as.numeric(RawDietSp$FracBioInDiet_S[grep("Ide_or_gudgeon", RawDietSp$Common_name)])/length(IndexTaxaS)
    IndexTaxaW <- intersect(IndexTaxa, IndexPreyW)
    DietSp$FracBioInDiet_W[IndexTaxa] <- DietSp$FracBioInDiet_W[IndexTaxa] + as.numeric(RawDietSp$FracBioInDiet_W[grep("Ide_or_gudgeon", RawDietSp$Common_name)])/length(IndexTaxaW)
    RawDietSp$IsUsed[grep("Ide_or_gudgeon", RawDietSp$Common_name)] <- TRUE
    # Undet. Cyprinidae
    Cyprinidae <- c("Abramis_brama", "Blicca_bjoerkna", "Carassius_carassius", "Gobio_gobio", "Leuciscus_idus", "Rutilus_rutilus", "Scardinius_erythrophthalmus", "Tinca_tinca")
    IndexTaxa <- match(Cyprinidae, DietSp$LowerTaxon)
    IndexTaxaS <- intersect(IndexTaxa, IndexPreyS)
    DietSp$FracBioInDiet_S[IndexTaxaS] <- DietSp$FracBioInDiet_S[IndexTaxaS] + as.numeric(RawDietSp$FracBioInDiet_S[grep("Undet._Cyprinidae", RawDietSp$Common_name)])/length(IndexTaxaS)
    IndexTaxaW <- intersect(IndexTaxa, IndexPreyW)
    DietSp$FracBioInDiet_W[IndexTaxaW] <- DietSp$FracBioInDiet_W[IndexTaxaW] + as.numeric(RawDietSp$FracBioInDiet_W[grep("Undet._Cyprinidae", RawDietSp$Common_name)])/length(IndexTaxaW)
    RawDietSp$IsUsed[grep("Undet._Cyprinidae", RawDietSp$Common_name)] <- TRUE
    # Undet. Cobitidae
    IndexTaxa <- which(DietSp$LowerTaxon == "Misgurnus_fossilis")
    IndexTaxaS <- intersect(IndexTaxa, IndexPreyS)
    DietSp$FracBioInDiet_S[IndexTaxaS] <- DietSp$FracBioInDiet_S[IndexTaxaS] + as.numeric(RawDietSp$FracBioInDiet_S[grep("Undet._Cobitidae", RawDietSp$Common_name)])
    IndexTaxaW <- intersect(IndexTaxa, IndexPreyW)
    DietSp$FracBioInDiet_W[IndexTaxaW] <- DietSp$FracBioInDiet_W[IndexTaxaW] + as.numeric(RawDietSp$FracBioInDiet_W[grep("Undet._Cobitidae", RawDietSp$Common_name)])
    RawDietSp$IsUsed[grep("Undet._Cobitidae", RawDietSp$Common_name)] <- TRUE
  }
  if (TableName == "p268.csv"){
    # Undet. birds
    PreyedBirds <- c("Columba_palumbus", "Erithacus_rubecula", "Turdus_merula", "Turdus_philomelos")
    IndexTaxa <- match(PreyedBirds, DietSp$LowerTaxon)
    DietSp$FracBioInDiet_S[IndexTaxa] <- as.numeric(RawDietSp$FracBioInDiet_S[grep("Undet._birds", RawDietSp$Common_name)])/length(PreyedBirds)
    DietSp$FracBioInDiet_W[IndexTaxa] <- as.numeric(RawDietSp$FracBioInDiet_W[grep("Undet._birds", RawDietSp$Common_name)])/length(PreyedBirds)
    RawDietSp$IsUsed[grep("Undet._birds", RawDietSp$Common_name)] <- TRUE
    # Undet. small mammals
    PreyedMamm <- c("Sorex_araneus", "Talpa_europaea")
    IndexTaxa <- match(PreyedMamm, DietSp$LowerTaxon)
    DietSp$FracBioInDiet_S[IndexTaxa] <- as.numeric(RawDietSp$FracBioInDiet_S[grep("Undet._small_mammals", RawDietSp$Common_name)])/length(PreyedMamm)
    DietSp$FracBioInDiet_W[IndexTaxa] <- as.numeric(RawDietSp$FracBioInDiet_W[grep("Undet._small_mammals", RawDietSp$Common_name)])/length(PreyedMamm)
    RawDietSp$IsUsed[grep("Undet._small_mammals", RawDietSp$Common_name)] <- TRUE
  }
  
  DietSp$FracBioInDiet_S[DietSp$FracBioInDiet_S == "+"] <- 0.01; DietSp$FracBioInDiet_S <- as.numeric(DietSp$FracBioInDiet_S)
  DietSp$FracBioInDiet_W[DietSp$FracBioInDiet_W == "+"] <- 0.01; DietSp$FracBioInDiet_W <- as.numeric(DietSp$FracBioInDiet_W)
  
  # copy-paste summer's preferences if winter's are not specified, and vice versa
  # and make sure that sum of preferences are equal 100
  if (length(intersect(Sp, PredWinterList)) == 1){ # if the predator is active on the study area during winter
    IndexPreyW <- which(!is.na(match(DietSp$LowerTaxon, PreyWinterList2)))
    # simplify dietary preferences / preferences are set to zero for species which are expected to migrate out of the study area during winter
    DietSp$FracBioInDiet_W[-IndexPreyW] <- 0
    if (sum(DietSp$FracBioInDiet_W, na.rm = TRUE) == 0){
      DietSp$FracBioInDiet_W[IndexPreyW] <- DietSp$FracBioInDiet_S[IndexPreyW]
    }
    DietSp$FracBioInDiet_W[IndexPreyW] <- DietSp$FracBioInDiet_W[IndexPreyW]/sum(DietSp$FracBioInDiet_W[IndexPreyW])*100
  }
  if (length(intersect(Sp, PredSummerList)) == 1){ # if the predator is active on the study area during summer
    IndexPreyS <- which(!is.na(match(DietSp$LowerTaxon, PreySummerList)))
    DietSp$FracBioInDiet_S[-IndexPreyS] <- 0
    if (sum(DietSp$FracBioInDiet_S, na.rm = TRUE) == 0){
      DietSp$FracBioInDiet_S[IndexPreyS] <- DietSp$FracBioInDiet_W[IndexPreyS]
    }
    DietSp$FracBioInDiet_S[IndexPreyS] <- DietSp$FracBioInDiet_S[IndexPreyS]/sum(DietSp$FracBioInDiet_S[IndexPreyS])*100
  }
  
  print(sum(DietSp$FracBioInDiet_S, na.rm = TRUE)); print(sum(DietSp$FracBioInDiet_W, na.rm = TRUE))
  IntListPref <- rbind(IntListPref, DietSp)
}

CstRes <- SpQuantities$Taxon[SpQuantities$Clade == "Other"]
PreyListModel <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Prey"])
IntListPref <- IntListPref[!is.na(match(IntListPref$LowerTaxon, PreyListModel)), ] # keep only prey which we model
# remove interactions with constant resources which are below 10% biomass within the predators' diets
IntListPref <- subset(IntListPref, (is.na(match(LowerTaxon, CstRes))) | ((FracBioInDiet_S >= 10) | (FracBioInDiet_W >= 10)))

for (IndexSp in 1:length(PredList)){
  Sp <- PredInfo$Taxon[IndexSp] # its name
  if ((sum(IntListPref$FracBioInDiet_S[IntListPref$UpperTaxon == Sp]) != 0) | (!is.na(sum(IntListPref$FracBioInDiet_S[IntListPref$UpperTaxon == Sp])))){
    IntListPref$FracBioInDiet_S[IntListPref$UpperTaxon == Sp] <- IntListPref$FracBioInDiet_S[IntListPref$UpperTaxon == Sp]/sum(IntListPref$FracBioInDiet_S[IntListPref$UpperTaxon == Sp], na.rm = TRUE)
  }
  if ((sum(IntListPref$FracBioInDiet_W[IntListPref$UpperTaxon == Sp]) != 0) & (!is.na(sum(IntListPref$FracBioInDiet_W[IntListPref$UpperTaxon == Sp])))){
    IntListPref$FracBioInDiet_W[IntListPref$UpperTaxon == Sp] <- IntListPref$FracBioInDiet_W[IntListPref$UpperTaxon == Sp]/sum(IntListPref$FracBioInDiet_W[IntListPref$UpperTaxon == Sp], na.rm = TRUE)
  }
}
IntListPref$FracBioInDiet_S <- IntListPref$FracBioInDiet_S*100
IntListPref$FracBioInDiet_W <- IntListPref$FracBioInDiet_W*100

setwd(WK_DIR)
write.csv(IntListPref, file = "PredPref.csv", quote = FALSE, row.names = FALSE)
