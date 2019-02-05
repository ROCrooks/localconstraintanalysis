#Function that calculates the constraint score
constraint.function <- function(obs,exp,adjustment)
  {
  #Calculate Z score
  numerator <- as.numeric(obs)-as.numeric(exp)
  denominator <- sqrt(as.numeric(exp))
  constraint <- numerator/denominator
  
  #Adjust using gene wide adjustment factor
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
variants.processing.function <- function(variant,genedetails)
  {
  #Retrieve the details about the gene
  wholegenesize <- genedetails[variant['Gene'],'Size']
  adjustment <- genedetails[variant['Gene'],'Adjustment']
  constraint.gene <- genedetails[variant['Gene'],'Constraint']
  
  constraint.region.raw <- constraint.function(variant['Region_Observed'],variant['Region_Expected'],adjustment)
  
  constraint.region.normalised <- constraint.function(normalise.function(variant['Region_Observed'],variant['Size'],wholegenesize),normalise.function(variant['Region_Expected'],variant['Size'],wholegenesize),adjustment)
  
  #Compile metrics into vector
  output <- c("Gene_Constraint"=constraint.gene,"Raw_Constraint"=constraint.region.raw,"Normalised_Constraint"=constraint.region.normalised)
  return(output)
  }

#Generate an example of a Z score
examplezscores <- function(value)
  {
  double = value*2
  zscore = value-double/sqrt(double)
  return(zscore)
  }

#Import the list of gene sizes
genedetails <- read.table("GeneDetails.txt", header=TRUE)

#Import the list of variants
variants.data.frame <- read.table("Variants.txt", header=TRUE)

#Calculate constraint metrics both raw and normalised
variants.constraints <- apply(variants.data.frame,1,variants.processing.function,genedetails)
variants.data.frame$Gene_Constraint <- variants.constraints['Gene_Constraint',]
variants.data.frame$Raw_Constraint <- variants.constraints['Raw_Constraint',]
variants.data.frame$Normalised_Constraint <- variants.constraints['Normalised_Constraint',]

#Calculate correlations between gene constraint and regional constraint
window.15bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V15bp'],variants.data.frame[,'Raw_Constraint'][variants.data.frame['Type'] == 'V15bp'])
window.30bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V30bp'],variants.data.frame[,'Raw_Constraint'][variants.data.frame['Type'] == 'V30bp'])
window.60bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V60bp'],variants.data.frame[,'Raw_Constraint'][variants.data.frame['Type'] == 'V60bp'])
window.90bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V90bp'],variants.data.frame[,'Raw_Constraint'][variants.data.frame['Type'] == 'V90bp'])

#Draw plot of local constraints at different window sizes
png("variantslocalconstraints.png",width=720,height=720)
par(mfrow=c(2,2))
par(oma=c(4,4,0,0)) 
par(mar=c(2,2,1,1))
plot(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V15bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V15bp'],xlab='',ylab='',main="+/- 15bp")
plot(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V30bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V30bp'],xlab='',ylab='',main="+/- 30bp")
plot(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V60bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V60bp'],xlab='',ylab='',main="+/- 60bp")
plot(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V90bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V90bp'],xlab='',ylab='',main="+/- 90bp")
mtext('Gene Constraint Score',side = 1, outer = TRUE, line = 2)
mtext('Regional Constraint Score',side = 2, outer = TRUE, line = 2)
dev.off()

#Draw gene scan plot of BRCA1
png("brca1scan.png",width=1080,height=720)
plot(
  variants.data.frame[,'Position']
  [variants.data.frame['Type'] == 'S15bp' & variants.data.frame['Gene'] == 'BRCA1']
  ,variants.data.frame[,'Normalised_Constraint']
  [variants.data.frame['Type'] == 'S15bp' & variants.data.frame['Gene'] == 'BRCA1']
  ,type='l'
  ,xlab='Position'
  ,ylab='Local Constraint'
  ,main='BRCA1 Regional Constraint Scan'
  ,lwd=2,col="yellow"
  )
lines(
  variants.data.frame[,'Position']
  [variants.data.frame['Type'] == 'S30bp' & variants.data.frame['Gene'] == 'BRCA1']
  ,variants.data.frame[,'Normalised_Constraint']
  [variants.data.frame['Type'] == 'S30bp' & variants.data.frame['Gene'] == 'BRCA1']
  ,lwd=2,col="blue"
)
lines(
  variants.data.frame[,'Position']
  [variants.data.frame['Type'] == 'S60bp' & variants.data.frame['Gene'] == 'BRCA1']
  ,variants.data.frame[,'Normalised_Constraint']
  [variants.data.frame['Type'] == 'S60bp' & variants.data.frame['Gene'] == 'BRCA1']
  ,lwd=2,col="green"
)
lines(
  variants.data.frame[,'Position']
  [variants.data.frame['Type'] == 'S90bp' & variants.data.frame['Gene'] == 'BRCA1']
  ,variants.data.frame[,'Normalised_Constraint']
  [variants.data.frame['Type'] == 'S90bp' & variants.data.frame['Gene'] == 'BRCA1']
  ,lwd=2,col="red"
)
legend(legend=c("15bp","30bp","60bp","90bp"),col=c("yellow","blue","green","red"),lwd=c("2","2","2","2"),lty=c(1,1,1,1),"topright")
dev.off()

#Make the plot of Z score examples
size = 100
zexamples <- data.frame(row.names=c(1:size))
zexamples$values <- c(1:size)
zscores <- apply(zexamples,1,examplezscores)
zexamples$zscore <- zscores
zexamples

png("zexamples.png",width=720,height=720)
plot(
  zexamples[,'values']
  ,zexamples[,'zscore']
  ,type='l'
  ,xlab='Observed Value'
  ,ylab='Z Score'
  ,main='Change in Z Score as Observed Increases'
)
dev.off()
