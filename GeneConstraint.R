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
  #Constraints in the gene details table are not adjusted
  constraint.gene <- genedetails[variant['Gene'],'Constraint']*adjustment
  
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

#Generate a gene vs regional variant graph
variantgraph <- function(x,y,name)
  {
  plot(x,y,xlab='',ylab='',main=name)
  }

#Draw gene scan plot
genescangraph <- function(gene)
  {
  #Declare file name
  filename = paste(c(gene,"scan.png"),sep="",collapse="")
  title = paste(c(gene," Regional Constraint Scan"),sep="",collapse="")
  
  #Draw gene scan plot of BRCA1
  png(filename,width=1080,height=720)
  plot(
    variants.data.frame[,'Position']
    [variants.data.frame['Type'] == 'S15bp' & variants.data.frame['Gene'] == gene]
    ,variants.data.frame[,'Normalised_Constraint']
    [variants.data.frame['Type'] == 'S15bp' & variants.data.frame['Gene'] == gene]
    ,type='l'
    ,xlab='Position'
    ,ylab='Local Constraint'
    ,main=title
    ,lwd=2,col="yellow"
  )
  lines(
    variants.data.frame[,'Position']
    [variants.data.frame['Type'] == 'S30bp' & variants.data.frame['Gene'] == gene]
    ,variants.data.frame[,'Normalised_Constraint']
    [variants.data.frame['Type'] == 'S30bp' & variants.data.frame['Gene'] == gene]
    ,lwd=2,col="blue"
  )
  lines(
    variants.data.frame[,'Position']
    [variants.data.frame['Type'] == 'S60bp' & variants.data.frame['Gene'] == gene]
    ,variants.data.frame[,'Normalised_Constraint']
    [variants.data.frame['Type'] == 'S60bp' & variants.data.frame['Gene'] == gene]
    ,lwd=2,col="green"
  )
  lines(
    variants.data.frame[,'Position']
    [variants.data.frame['Type'] == 'S90bp' & variants.data.frame['Gene'] == gene]
    ,variants.data.frame[,'Normalised_Constraint']
    [variants.data.frame['Type'] == 'S90bp' & variants.data.frame['Gene'] == gene]
    ,lwd=2,col="red"
  )
  #Draw 0 line
  abline(h=0)
  #Draw horizontal line at the 3SD line for statistical significance
  abline(h=3.09)
  legend(legend=c("15bp","30bp","60bp","90bp"),col=c("yellow","blue","green","red"),lwd=c("2","2","2","2"),lty=c(1,1,1,1),"topright")
  dev.off()
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
window.15bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V15bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V15bp'])
window.30bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V30bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V30bp'])
window.60bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V60bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V60bp'])
window.90bp.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V90bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V90bp'])
window.exon.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'Exon'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'Exon'])
window.domain.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'Domain'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'Domain'])

#Draw plot of local constraints at different window sizes
png("variantslocalconstraints.png",width=1080,height=720)
par(mfrow=c(2,3))
par(oma=c(4,4,0,0)) 
par(mar=c(2,2,1,1))
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V15bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V15bp'],"+/- 15bp")
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V30bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V30bp'],"+/- 30bp")
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'Exon'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'Exon'],"Exon")
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V60bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V60bp'],"+/- 60bp")
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V90bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V90bp'],"+/- 90bp")
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'Domain'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'Domain'],"Domain")
mtext('Gene Constraint Score',side = 1, outer = TRUE, line = 2)
mtext('Regional Constraint Score',side = 2, outer = TRUE, line = 2)
dev.off()

#Make th gene scan plots for BRCA1 and NF1
genescangraph("BRCA1")
genescangraph("NF1")

#Make the plot of Z score examples
size = 100
zexamples <- data.frame(row.names=c(1:size))
zexamples$values <- c(1:size)
zscores <- apply(zexamples,1,examplezscores)
zexamples$zscore <- zscores

#Plot graph showing that Z score increases as Observed increases
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

#Count the numbers of variants and each window size
v15.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V15bp',])
v30.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V30bp',])
v60.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V60bp',])
v90.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V90bp',])
exon.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'Exon',])
domain.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'Domain',])

#Count the number of variants where gene constraint is already >3.09
v15.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V15bp' & variants.data.frame$Gene_Constraint >= 3.09,])
v30.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V30bp' & variants.data.frame$Gene_Constraint >= 3.09,])
v60.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V60bp' & variants.data.frame$Gene_Constraint >= 3.09,])
v90.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V90bp' & variants.data.frame$Gene_Constraint >= 3.09,])
exon.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'Exon' & variants.data.frame$Gene_Constraint >= 3.09,])
domain.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'Domain' & variants.data.frame$Gene_Constraint >= 3.09,])

#Count the number of instances where a local constraint exceeds 3.09, but gene constraint doesn't
v15.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V15bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
v30.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V30bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
v60.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V60bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
v90.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V90bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
exon.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'Exon' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
domain.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'Domain' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
                                           
v15.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V15bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])
v30.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V30bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])
v60.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V60bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])
v90.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V90bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])
exon.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'Exon' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])
domain.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'Domain' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])

#Calculate percentages
v15.relevent.percent.normalised <- round((v15.relevantnormalised.variants/v15.total.variants)*100, digits=3)
v30.relevent.percent.normalised <- round((v30.relevantnormalised.variants/v30.total.variants)*100, digits=3)
v60.relevent.percent.normalised <- round((v60.relevantnormalised.variants/v60.total.variants)*100, digits=3)
v90.relevent.percent.normalised <- round((v90.relevantraw.variants/v90.total.variants)*100, digits=3)
exon.relevent.percent.normalised <- round((exon.relevantraw.variants/exon.total.variants)*100, digits=3)
domain.relevent.percent.normalised <- round((domain.relevantraw.variants/domain.total.variants)*100, digits=3)

v15.relevent.percent.raw <- round((v15.relevantraw.variants/v15.total.variants)*100, digits=3)
v30.relevent.percent.raw <- round((v30.relevantraw.variants/v30.total.variants)*100, digits=3)
v60.relevent.percent.raw <- round((v60.relevantraw.variants/v60.total.variants)*100, digits=3)
v90.relevent.percent.raw <- round((v90.relevantraw.variants/v90.total.variants)*100, digits=3)
exon.relevent.percent.raw <- round((exon.relevantraw.variants/exon.total.variants)*100, digits=3)
domain.relevent.percent.raw <- round((domain.relevantraw.variants/domain.total.variants)*100, digits=3)

#Calculate uplifts
v15.relevent.uplift.normalised <- round((v15.relevantnormalised.variants/v15.constrainedgenes.variants)*100, digits=3)
v30.relevent.uplift.normalised <- round((v30.relevantnormalised.variants/v30.constrainedgenes.variants)*100, digits=3)
v60.relevent.uplift.normalised <- round((v60.relevantnormalised.variants/v60.constrainedgenes.variants)*100, digits=3)
v90.relevent.uplift.normalised <- round((v90.relevantraw.variants/v90.constrainedgenes.variants)*100, digits=3)
exon.relevent.uplift.normalised <- round((exon.relevantraw.variants/exon.constrainedgenes.variants)*100, digits=3)
domain.relevent.uplift.normalised <- round((domain.relevantraw.variants/domain.constrainedgenes.variants)*100, digits=3)

v15.relevent.uplift.raw <- round((v15.relevantraw.variants/v15.constrainedgenes.variants)*100, digits=3)
v30.relevent.uplift.raw <- round((v30.relevantraw.variants/v30.constrainedgenes.variants)*100, digits=3)
v60.relevent.uplift.raw <- round((v60.relevantraw.variants/v60.constrainedgenes.variants)*100, digits=3)
v90.relevent.uplift.raw <- round((v90.relevantraw.variants/v90.constrainedgenes.variants)*100, digits=3)
exon.relevent.uplift.raw <- round((exon.relevantraw.variants/exon.constrainedgenes.variants)*100, digits=3)
domain.relevent.uplift.raw <- round((domain.relevantraw.variants/domain.constrainedgenes.variants)*100, digits=3)

#Write data to output table text file
write(paste(c("Window Type","Constrained Genes","Constrained Normalised","Constrained Unnormalised","Uplift Normalised","Uplift Unnormalised"),collapse="\t"),"outputtable.txt")
output.table.header.v15 <- "+/- 15bp"
output.table.normalised.v15 <- paste(c(v15.relevent.percent.normalised,"% (",v15.relevantnormalised.variants,"/",v15.total.variants,")"),collapse="")
output.table.raw.v15 <- paste(c(v15.relevent.percent.raw,"% (",v15.relevantraw.variants,"/",v15.total.variants,")"),collapse="")
write(paste(c(output.table.header.v15,output.table.normalised.v15,output.table.raw.v15),collapse="\t"),"outputtable.txt",append=TRUE)
output.table.header.v30 <- "+/- 30bp"
output.table.normalised.v30 <- paste(c(v30.relevent.percent.normalised,"% (",v30.relevantnormalised.variants,"/",v30.total.variants,")"),collapse="")
output.table.raw.v30 <- paste(c(v30.relevent.percent.raw,"% (",v30.relevantraw.variants,"/",v30.total.variants,")"),collapse="")
write(paste(c(output.table.header.v30,output.table.normalised.v30,output.table.raw.v30),collapse="\t"),"outputtable.txt",append=TRUE)
output.table.header.v60 <- "+/- 60bp"
output.table.normalised.v60 <- paste(c(v60.relevent.percent.normalised,"% (",v60.relevantnormalised.variants,"/",v60.total.variants,")"),collapse="")
output.table.raw.v60 <- paste(c(v60.relevent.percent.raw,"% (",v60.relevantraw.variants,"/",v60.total.variants,")"),collapse="")
write(paste(c(output.table.header.v60,output.table.normalised.v60,output.table.raw.v60),collapse="\t"),"outputtable.txt",append=TRUE)
output.table.header.v90 <- "+/- 90bp"
output.table.normalised.v90 <- paste(c(v90.relevent.percent.normalised,"% (",v90.relevantnormalised.variants,"/",v90.total.variants,")"),collapse="")
output.table.raw.v90 <- paste(c(v90.relevent.percent.raw,"% (",v90.relevantraw.variants,"/",v90.total.variants,")"),collapse="")
write(paste(c(output.table.header.v90,output.table.normalised.v90,output.table.raw.v90),collapse="\t"),"outputtable.txt",append=TRUE)
output.table.header.exon <- "Exon"
output.table.normalised.exon <- paste(c(exon.relevent.percent.normalised,"% (",exon.relevantnormalised.variants,"/",exon.total.variants,")"),collapse="")
output.table.raw.exon <- paste(c(exon.relevent.percent.raw,"% (",exon.relevantraw.variants,"/",exon.total.variants,")"),collapse="")
write(paste(c(output.table.header.exon,output.table.normalised.exon,output.table.raw.exon),collapse="\t"),"outputtable.txt",append=TRUE)
output.table.header.domain <- "Domain"
output.table.normalised.domain <- paste(c(domain.relevent.percent.normalised,"% (",domain.relevantnormalised.variants,"/",domain.total.variants,")"),collapse="")
output.table.raw.domain <- paste(c(domain.relevent.percent.raw,"% (",domain.relevantraw.variants,"/",domain.total.variants,")"),collapse="")
write(paste(c(output.table.header.domain,output.table.normalised.domain,output.table.raw.domain),collapse="\t"),"outputtable.txt",append=TRUE)