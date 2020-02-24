# dietBrayCurtis() returns the Bray-Curtis dissimilarity between two diets (vectors)

dietBrayCurtis <- function(diet1, diet2){
  # diet1 and diet2 are two vectors of same length describing two diets
  
  if ((!is.vector(diet1)) | (!is.vector(diet1))){
    stop("diet1 and diet2 are two vectors of same length describing two diets.")
  }
  if (length(diet1) != length(diet2)){
    stop("diet1 and diet2 must have the same length.")
  }
  BC_min <- 0
  BC_sum <- 0
  
  N <- length(diet1) # number of potentiel resources in the diet
  for (i in 1:N){
    BC_min <- BC_min+min(diet1[i], diet2[i])
    BC_sum <- BC_sum+(diet1[i]+diet2[i])
  }
  BC <- 1 - 2*BC_min/BC_sum
  return(BC)
}