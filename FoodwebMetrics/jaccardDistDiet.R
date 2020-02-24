# jaccardDistDiet() returns the distance between two diets.

jaccardDistDiet <- function(diet1, diet2){
  if ((!is.vector(diet1)) | (!is.vector(diet1))){
    stop("diet1 and diet2 are two vectors of same length describing two diets.")
  }
  if (length(diet1) != length(diet2)){
    stop("diet1 and diet2 must have the same length.")
  }
  
  diet1 <- which(diet1 != 0)
  diet2 <- which(diet2 != 0)
  dietInter <- length(intersect(diet1, diet2)) # size of the intersection between the two diets
  dietUnion <- length(union(diet1, diet2)) # size of the union of the two diets
  jaccIndex <- dietInter/dietUnion
  jaccDist <- 1-jaccIndex
  return(jaccDist)
}