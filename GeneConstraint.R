#Function for a constraint score
constraint.function <- function(obs,exp,adjustment)
  {
  numerator <- obs-exp
  denominator <- sqrt(exp)
  
  constraint <- numerator/denominator
  constraint <- constraint*adjustment
  
  return(constraint)
  }

#Function to normalise a observed or expected frequency
normalise.function <- function(value,range,whole)
  {
  normalisation <- whole/range
  newvalue <- value*normalisation
  return(newvalue)
  }

#Import the list of gene sizes
genesizes <- read.table("GeneSizes.txt", header=TRUE)

#Import the list of variants
variants <- read.table("Variants.txt", header=TRUE)

#Import the results of the gene scans
brca1scan <- read.table("BCRA1Scan.txt", header=TRUE)
nf1scan <- read.table("NF1Scan.txt", header=TRUE)


