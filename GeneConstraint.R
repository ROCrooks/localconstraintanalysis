#Function for a constraint score
constraint.function <- function(obs,exp,adjustment)
  {
  numerator <- as.numeric(obs)-as.numeric(exp)
  denominator <- sqrt(as.numeric(exp))
  
  constraint <- numerator/denominator
  constraint <- constraint*as.numeric(adjustment)
  
  return(constraint)
  }

#Function to normalise a observed or expected frequency
normalise.function <- function(value,range,whole)
  {
  normalisation <- as.numeric(whole)/as.numeric(range)
  newvalue <- as.numeric(value)*normalisation
  return(newvalue)
  }

#Function that processes variant constraint scores 
variants.processing.function <- function(x,genedetails)
  {
  #Retrieve the details about the gene
  wholegenesize <- genedetails[x['Gene'],'Size']
  adjustmentfactor <- [gene,'Adjustment_Factor']
  
  #Calculate the unnormalised constraint metrics
  constraint.gene <- constraint.function(genedetails[x['Gene'],'Observed'],genedetails[x['Gene'],'Expected'],adjustmentfactor)
  constraint.X15 <- constraint.function(x['X15_Observed'],x['X15_Expected'],adjustmentfactor)
  constraint.X30 <- constraint.function(x['X30_Observed'],x['X30_Expected'],adjustmentfactor)
  constraint.X60 <- constraint.function(x['X60_Observed'],x['X60_Expected'],adjustmentfactor)
  constraint.X90 <- constraint.function(x['X90_Observed'],x['X90_Expected'],adjustmentfactor)
  
  #Calculate normalised constraint metrics
  constraint.X15.normalised <- constraint.function(normalise.function(x['X15_Observed'],15,wholegenesize),normalise.function(x['X15_Expected'],15,wholegenesize),adjustmentfactor)
  constraint.X30.normalised <- constraint.function(normalise.function(x['X30_Observed'],30,wholegenesize),normalise.function(x['X30_Expected'],30,wholegenesize),adjustmentfactor)
  constraint.X60.normalised <- constraint.function(normalise.function(x['X60_Observed'],60,wholegenesize),normalise.function(x['X60_Expected'],60,wholegenesize),adjustmentfactor)
  constraint.X90.normalised <- constraint.function(normalise.function(x['X90_Observed'],90,wholegenesize),normalise.function(x['X90_Expected'],90,wholegenesize),adjustmentfactor)
  
  #Compile metrics into vector
  output <- c("Gene_Constraint"=constraint.gene,"X15_Constraint"=constraint.X15,"X30_Constraint"=constraint.X30,"X60_Constraint"=constraint.X60,"X90_Constraint"=constraint.X90,"X15_Constraint_Normalised"=constraint.X15.normalised,"X30_Constraint_Normalised"=constraint.X30.normalised,"X60_Constraint_Normalised"=constraint.X60.normalised,"X90_Constraint_Normalised"=constraint.X90.normalised)
  return(output)
  }

#Gene scan processing function
genescan.processing.function <- function(x,genedetails,gene)
  {
  #Retrieve the details about the gene
  wholegenesize <- genedetails[gene,'Size']
  adjustmentfactor <- [gene,'Adjustment_Factor']
  
  #Calculate the unnormalised constraint metrics
  constraint.X15 <- constraint.function(x['X15_Observed'],x['X15_Expected'],adjustmentfactor)
  constraint.X30 <- constraint.function(x['X30_Observed'],x['X30_Expected'],adjustmentfactor)
  constraint.X60 <- constraint.function(x['X60_Observed'],x['X60_Expected'],adjustmentfactor)
  constraint.X90 <- constraint.function(x['X90_Observed'],x['X90_Expected'],adjustmentfactor)
  
  #Calculate normalised constraint metrics
  constraint.X15.normalised <- constraint.function(normalise.function(x['X15_Observed'],15,wholegenesize),normalise.function(x['X15_Expected'],15,wholegenesize),adjustmentfactor)
  constraint.X30.normalised <- constraint.function(normalise.function(x['X30_Observed'],30,wholegenesize),normalise.function(x['X30_Expected'],30,wholegenesize),adjustmentfactor)
  constraint.X60.normalised <- constraint.function(normalise.function(x['X60_Observed'],60,wholegenesize),normalise.function(x['X60_Expected'],60,wholegenesize),adjustmentfactor)
  constraint.X90.normalised <- constraint.function(normalise.function(x['X90_Observed'],90,wholegenesize),normalise.function(x['X90_Expected'],90,wholegenesize),adjustmentfactor)
  
  #Compile metrics into vector
  output <- c("X15_Constraint"=constraint.X15,"X30_Constraint"=constraint.X30,"X60_Constraint"=constraint.X60,"X90_Constraint"=constraint.X90,"X15_Constraint_Normalised"=constraint.X15.normalised,"X30_Constraint_Normalised"=constraint.X30.normalised,"X60_Constraint_Normalised"=constraint.X60.normalised,"X90_Constraint_Normalised"=constraint.X90.normalised)
  return(output)
  }

#Attach constraints to the selected array
bind.constraints.function <- function(array,constraints)
  {
  array$X15_Constraint <- constraints['X15_Constraint',]
  array$X30_Constraint <- constraints['X30_Constraint',]
  array$X60_Constraint <- constraints['X60_Constraint',]
  array$X90_Constraint <- constraints['X90_Constraint',]
  array$X15_Constraint_Normalised <- constraints['X15_Constraint_Normalised',]
  array$X30_Constraint_Normalised <- constraints['X30_Constraint_Normalised',]
  array$X60_Constraint_Normalised <- constraints['X60_Constraint_Normalised',]
  array$X90_Constraint_Normalised <- constraints['X90_Constraint_Normalised',]
  
  return(array)
  }

#Import the list of gene sizes
genedetails <- read.table("GeneDetails.txt", header=TRUE)

#Import the list of variants
variants.data.frame <- read.table("Variants.txt", header=TRUE)

#Import the results of the gene scans
brca1scan.data.frame <- read.table("BCRA1Scan.txt", header=TRUE)
nf1scan.data.frame <- read.table("NF1Scan.txt", header=TRUE)

#Calculate constraint metrics both raw and normalised
variants.constraints <- apply(variants.data.frame,1,variants.processing.function,genedetails)
variants.data.frame$Gene_Constraint <- variants.constraints['Gene_Constraint',]
variants.data.frame <- bind.constraints.function(variants.data.frame,variants.constraints)

#Calculate constraint metrics across entire gene for BRCA1 and NF1
brca1scan.constraints <- apply(brca1scan.data.frame,1,genescan.processing.function,genedetails,"BRCA1")
brca1scan.data.frame <- bind.constraints.function(brca1scan.data.frame,brca1scan.constraints)
nf1scan.constraints <- apply(nf1scan.data.frame,1,genescan.processing.function,genedetails,"NF1")
nf1scan.data.frame <- bind.constraints.function(nf1scan.data.frame,nf1scan.constraints)

#Calculate correlations between gene constraint and regional constraint
window.15bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'],variants.data.frame[,'X15_Constraint_Normalised'])
window.30bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'],variants.data.frame[,'X30_Constraint_Normalised'])
window.60bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'],variants.data.frame[,'X60_Constraint_Normalised'])
window.90bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'],variants.data.frame[,'X90_Constraint_Normalised'])

#Draw plot of local constraints at different window sizes
png("variantslocalconstraints.png",width=720,height=720)
par(mfrow=c(2,2))
par(oma=c(4,4,0,0)) 
par(mar=c(2,2,1,1))
plot(variants.data.frame[,'Gene_Constraint'],variants.data.frame[,'X15_Constraint_Normalised'],xlab='',ylab='',main="+/- 15bp")
plot(variants.data.frame[,'Gene_Constraint'],variants.data.frame[,'X30_Constraint_Normalised'],xlab='',ylab='',main="+/- 30bp")
plot(variants.data.frame[,'Gene_Constraint'],variants.data.frame[,'X60_Constraint_Normalised'],xlab='',ylab='',main="+/- 60bp")
plot(variants.data.frame[,'Gene_Constraint'],variants.data.frame[,'X90_Constraint_Normalised'],xlab='',ylab='',main="+/- 90bp")
mtext('Gene Constraint Score',side = 1, outer = TRUE, line = 2)
mtext('Regional Constraint Score',side = 2, outer = TRUE, line = 2)
dev.off()
