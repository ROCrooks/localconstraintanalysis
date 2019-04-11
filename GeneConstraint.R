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

#Generate a gene vs regional variant graph
variantgraph <- function(x,y,name)
  {
  plot(x,y,xlab='',ylab='',main=name,cex.main=1.5)
  }

#Draw gene scan plot
genescangraph <- function(gene)
  {
  #Declare file name
  filename <- paste(c(gene,"scan.png"),sep="",collapse="")
  title <- paste(c(gene," Regional Constraint Scan"),sep="",collapse="")
  
  #Define domain boundaries according to gene
  if (gene == "BRCA1")
    {
    domain.starts <- c(4924,5266)
    domain.ends <- c(5208,5565)
    }
  else if (gene == "NF1")
    {
    domain.starts <- c(3703,4738)
    domain.ends <- c(4353,5214)
    }
  #Draw gene scan plot of BRCA1
  png(filename,width=1080,height=720)
  plot(
    #Highlight areas within domains
    panel.first = rect(domain.starts, -1e6, domain.ends, 1e6, col='grey70', border=NA),
    xaxs = "i",
    yaxs = "i",
    
    variants.data.frame[,'Position']
    [variants.data.frame['Type'] == 'S15bp' & variants.data.frame['Gene'] == gene]
    ,variants.data.frame[,'Normalised_Constraint']
    [variants.data.frame['Type'] == 'S15bp' & variants.data.frame['Gene'] == gene]
    ,type='l'
    ,xlab='Position'
    ,ylab='Local Constraint'
    ,cex.lab = 1.5
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
  lines(
    variants.data.frame[,'Position']
    [variants.data.frame['Type'] == 'SExon' & variants.data.frame['Gene'] == gene]
    ,variants.data.frame[,'Normalised_Constraint']
    [variants.data.frame['Type'] == 'SExon' & variants.data.frame['Gene'] == gene]
    ,lwd=2,col="purple"
  )
  #Draw 0 line
  abline(h=0)
  #Draw horizontal line at the 3SD line for statistical significance
  abline(h=3.09)
  legend(legend=c("15bp","30bp","60bp","90bp","Exon"),col=c("yellow","blue","green","red","purple"),lwd=c("2","2","2","2"),lty=c(1,1,1,1),"topright")
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
window.exon.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'VExon'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'VExon'])
window.domain.correlation <- cor(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'VDomain'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'VDomain'])
write(paste(c("Correlations (R2) were: 15bp = ",round(window.15bp.correlation, digits=3),", 30bp = ",round(window.30bp.correlation, digits=3),", 60bp = ",round(window.60bp.correlation, digits=3),", 90bp = ",round(window.90bp.correlation, digits=3),", within exon = ",round(window.exon.correlation, digits=3),", within domain = ",round(window.domain.correlation, digits=3)),collapse=""),"outputcorrelations.txt")

#Draw plot of local constraints at different window sizes
png("variantslocalconstraints.png",width=1080,height=720)
par(mfrow=c(2,3))
par(oma=c(4,4,0,0)) 
par(mar=c(2,2,1,1))
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V15bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V15bp'],"+/- 15bp")
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V30bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V30bp'],"+/- 30bp")
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'VExon'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'VExon'],"Exon")
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V60bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V60bp'],"+/- 60bp")
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'V90bp'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'V90bp'],"+/- 90bp")
variantgraph(variants.data.frame[,'Gene_Constraint'][variants.data.frame['Type'] == 'VDomain'],variants.data.frame[,'Normalised_Constraint'][variants.data.frame['Type'] == 'VDomain'],"Domain")
mtext('Gene Constraint Score',side = 1, outer = TRUE, line = 2, cex=1.5)
mtext('Regional Constraint Score',side = 2, outer = TRUE, line = 2, cex=1.5)
dev.off()

#Make th gene scan plots for BRCA1 and NF1
genescangraph("BRCA1")
genescangraph("NF1")

#Make the plot of Z score examples
z.expected <- c(0:100)
z.score <- ((z.expected*2)-z.expected)/(sqrt(z.expected))
zexamples <- data.frame(Expected=z.expected,ZScore=z.score)

#Plot graph showing that Z score increases as Observed increases
png("zexamples.png",width=720,height=720)
plot(
  zexamples[,'Expected'],
  zexamples[,'ZScore'],
  type='l',
  xlab='Observed Value',
  ylab='Z Score',
  xaxs = "i",
  yaxs = "i",
  cex.lab = 1.5,
  xlim = c(0,100),
  ylim = c(0,10)
)
dev.off()

#Count the numbers of variants and each window size
v15.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V15bp',])
v30.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V30bp',])
v60.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V60bp',])
v90.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V90bp',])
exon.total.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'VExon',])
domain.total.variants <- nrow(unique(subset(variants.data.frame, Type == "VDomain", select=c(Key, VariantName))))

#Count the number of variants where gene constraint is already >3.09
v15.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V15bp' & variants.data.frame$Gene_Constraint >= 3.09,])
v30.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V30bp' & variants.data.frame$Gene_Constraint >= 3.09,])
v60.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V60bp' & variants.data.frame$Gene_Constraint >= 3.09,])
v90.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V90bp' & variants.data.frame$Gene_Constraint >= 3.09,])
exon.constrainedgenes.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'VExon' & variants.data.frame$Gene_Constraint >= 3.09,])
domain.constrainedgenes.variants <- nrow(unique(subset(variants.data.frame, Type == "VDomain" & Gene_Constraint >= 3.09, select=c(Key, VariantName))))

#Calculate the percentage of variants where gene constraint is already >3.09
v15.constrainedgenes.percentage <- round((v15.constrainedgenes.variants/v15.total.variants)*100, digits=3)
v30.constrainedgenes.percentage <- round((v30.constrainedgenes.variants/v30.total.variants)*100, digits=3)
v60.constrainedgenes.percentage <- round((v60.constrainedgenes.variants/v60.total.variants)*100, digits=3)
v90.constrainedgenes.percentage <- round((v90.constrainedgenes.variants/v90.total.variants)*100, digits=3)
exon.constrainedgenes.percentage <- round((exon.constrainedgenes.variants/exon.total.variants)*100, digits=3)
domain.constrainedgenes.percentage <- round((domain.constrainedgenes.variants/domain.total.variants)*100, digits=3)

#Count the number of instances where a local constraint exceeds 3.09, but gene constraint doesn't
v15.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V15bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
v30.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V30bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
v60.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V60bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
v90.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V90bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
exon.relevantnormalised.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'VExon' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Normalised_Constraint >= 3.09,])
domain.relevantnormalised.variants <- nrow(unique(subset(variants.data.frame, Type == "VDomain" & Gene_Constraint < 3.09 & Normalised_Constraint >= 3.09, select=c(Key, VariantName))))

v15.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V15bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])
v30.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V30bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])
v60.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V60bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])
v90.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'V90bp' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])
exon.relevantraw.variants <- nrow(variants.data.frame[variants.data.frame$Type == 'VExon' & variants.data.frame$Gene_Constraint < 3.09 & variants.data.frame$Raw_Constraint >= 3.09,])
domain.relevantraw.variants <- nrow(unique(subset(variants.data.frame, Type == "VDomain" & Gene_Constraint < 3.09 & Raw_Constraint >= 3.09, select=c(Key, VariantName))))

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

#Write data for uplift to output table text file
write(paste(c("Window Type","Constrained Genes","Constrained Normalised","Constrained Unnormalised","Uplift Normalised","Uplift Unnormalised"),collapse="\t"),"outputtable.txt")
output.table.header.v15 <- "+/- 15bp"
output.table.constrained.gene.v15 <- paste(c(v15.constrainedgenes.percentage,"% (",v15.constrainedgenes.variants,"/",v15.total.variants,")"),collapse="")
output.table.normalised.v15 <- paste(c(v15.relevent.percent.normalised,"% (",v15.relevantnormalised.variants,"/",v15.total.variants,")"),collapse="")
output.table.raw.v15 <- paste(c(v15.relevent.percent.raw,"% (",v15.relevantraw.variants,"/",v15.total.variants,")"),collapse="")
output.table.uplift.normalised.v15 <- paste(c(v15.relevent.uplift.normalised,"%"),collapse="")
output.table.uplift.raw.v15 <- paste(c(v15.relevent.uplift.raw,"%"),collapse="")
write(paste(c(output.table.header.v15,output.table.constrained.gene.v15,output.table.normalised.v15,output.table.raw.v15,output.table.uplift.normalised.v15,output.table.uplift.raw.v15),collapse="\t"),"outputtable.txt",append=TRUE)
output.table.header.v30 <- "+/- 30bp"
output.table.constrained.gene.v30 <- paste(c(v30.constrainedgenes.percentage,"% (",v30.constrainedgenes.variants,"/",v30.total.variants,")"),collapse="")
output.table.normalised.v30 <- paste(c(v30.relevent.percent.normalised,"% (",v30.relevantnormalised.variants,"/",v30.total.variants,")"),collapse="")
output.table.raw.v30 <- paste(c(v30.relevent.percent.raw,"% (",v30.relevantraw.variants,"/",v30.total.variants,")"),collapse="")
output.table.uplift.normalised.v30 <- paste(c(v30.relevent.uplift.normalised,"%"),collapse="")
output.table.uplift.raw.v30 <- paste(c(v30.relevent.uplift.raw,"%"),collapse="")
write(paste(c(output.table.header.v30,output.table.constrained.gene.v30,output.table.normalised.v30,output.table.raw.v30,output.table.uplift.normalised.v30,output.table.uplift.raw.v30),collapse="\t"),"outputtable.txt",append=TRUE)
output.table.header.v60 <- "+/- 60bp"
output.table.constrained.gene.v60 <- paste(c(v60.constrainedgenes.percentage,"% (",v60.constrainedgenes.variants,"/",v60.total.variants,")"),collapse="")
output.table.normalised.v60 <- paste(c(v60.relevent.percent.normalised,"% (",v60.relevantnormalised.variants,"/",v60.total.variants,")"),collapse="")
output.table.raw.v60 <- paste(c(v60.relevent.percent.raw,"% (",v60.relevantraw.variants,"/",v60.total.variants,")"),collapse="")
output.table.uplift.normalised.v60 <- paste(c(v60.relevent.uplift.normalised,"%"),collapse="")
output.table.uplift.raw.v60 <- paste(c(v60.relevent.uplift.raw,"%"),collapse="")
write(paste(c(output.table.header.v60,output.table.constrained.gene.v60,output.table.normalised.v60,output.table.raw.v60,output.table.uplift.normalised.v60,output.table.uplift.raw.v60),collapse="\t"),"outputtable.txt",append=TRUE)
output.table.header.v90 <- "+/- 90bp"
output.table.constrained.gene.v90 <- paste(c(v90.constrainedgenes.percentage,"% (",v90.constrainedgenes.variants,"/",v90.total.variants,")"),collapse="")
output.table.normalised.v90 <- paste(c(v90.relevent.percent.normalised,"% (",v90.relevantnormalised.variants,"/",v90.total.variants,")"),collapse="")
output.table.raw.v90 <- paste(c(v90.relevent.percent.raw,"% (",v90.relevantraw.variants,"/",v90.total.variants,")"),collapse="")
output.table.uplift.normalised.v90 <- paste(c(v90.relevent.uplift.normalised,"%"),collapse="")
output.table.uplift.raw.v90 <- paste(c(v90.relevent.uplift.raw,"%"),collapse="")
write(paste(c(output.table.header.v90,output.table.constrained.gene.v90,output.table.normalised.v90,output.table.raw.v90,output.table.uplift.normalised.v90,output.table.uplift.raw.v90),collapse="\t"),"outputtable.txt",append=TRUE)
output.table.header.exon <- "Exon"
output.table.constrained.gene.exon <- paste(c(exon.constrainedgenes.percentage,"% (",exon.constrainedgenes.variants,"/",exon.total.variants,")"),collapse="")
output.table.normalised.exon <- paste(c(exon.relevent.percent.normalised,"% (",exon.relevantnormalised.variants,"/",exon.total.variants,")"),collapse="")
output.table.raw.exon <- paste(c(exon.relevent.percent.raw,"% (",exon.relevantraw.variants,"/",exon.total.variants,")"),collapse="")
output.table.uplift.normalised.exon <- paste(c(exon.relevent.uplift.normalised,"%"),collapse="")
output.table.uplift.raw.exon <- paste(c(exon.relevent.uplift.raw,"%"),collapse="")
write(paste(c(output.table.header.exon,output.table.constrained.gene.exon,output.table.normalised.exon,output.table.raw.exon,output.table.uplift.normalised.exon,output.table.uplift.raw.exon),collapse="\t"),"outputtable.txt",append=TRUE)
output.table.header.domain <- "Domain"
output.table.constrained.gene.domain <- paste(c(domain.constrainedgenes.percentage,"% (",domain.constrainedgenes.variants,"/",domain.total.variants,")"),collapse="")
output.table.normalised.domain <- paste(c(domain.relevent.percent.normalised,"% (",domain.relevantnormalised.variants,"/",domain.total.variants,")"),collapse="")
output.table.raw.domain <- paste(c(domain.relevent.percent.raw,"% (",domain.relevantraw.variants,"/",domain.total.variants,")"),collapse="")
output.table.uplift.normalised.domain <- paste(c(domain.relevent.uplift.normalised,"%"),collapse="")
output.table.uplift.raw.domain <- paste(c(domain.relevent.uplift.raw,"%"),collapse="")
write(paste(c(output.table.header.domain,output.table.constrained.gene.domain,output.table.normalised.domain,output.table.raw.domain,output.table.uplift.normalised.domain,output.table.uplift.raw.domain),collapse="\t"),"outputtable.txt",append=TRUE)

#Extract variant names where the local constraint made a difference to classification
difference.v15.Normalised <- subset(variants.data.frame, Type == "V15bp" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09, select=c(Key, VariantName))
difference.v30.Normalised <- subset(variants.data.frame, Type == "V30bp" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09, select=c(Key, VariantName))
difference.v60.Normalised <- subset(variants.data.frame, Type == "V60bp" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09, select=c(Key, VariantName))
difference.v90.Normalised <- subset(variants.data.frame, Type == "V90bp" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09, select=c(Key, VariantName))
difference.Exon.Normalised <- subset(variants.data.frame, Type == "VExon" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09, select=c(Key, VariantName))
difference.Domain.Normalised <- subset(variants.data.frame, Type == "VDomain" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09, select=c(Key, VariantName))
difference.Domain.Normalised <- unique(difference.Domain.Normalised)

difference.v15.Unnormalised <- subset(variants.data.frame, Type == "V15bp" & Gene_Constraint < 3.09 & Raw_Constraint > 3.09, select=c(Key, VariantName))
difference.v30.Unnormalised <- subset(variants.data.frame, Type == "V30bp" & Gene_Constraint < 3.09 & Raw_Constraint > 3.09, select=c(Key, VariantName))
difference.v60.Unnormalised <- subset(variants.data.frame, Type == "V60bp" & Gene_Constraint < 3.09 & Raw_Constraint > 3.09, select=c(Key, VariantName))
difference.v90.Unnormalised <- subset(variants.data.frame, Type == "V90bp" & Gene_Constraint < 3.09 & Raw_Constraint > 3.09, select=c(Key, VariantName))
difference.Exon.Unnormalised <- subset(variants.data.frame, Type == "VExon" & Gene_Constraint < 3.09 & Raw_Constraint > 3.09, select=c(Key, VariantName))
difference.Domain.Unnormalised <- subset(variants.data.frame, Type == "VDomain" & Gene_Constraint < 3.09 & Raw_Constraint > 3.09, select=c(Key, VariantName))
difference.Domain.Unnormalised <- unique(difference.Domain.Unnormalised)

#Make plot of the number of variants that are locally constrained
typecounts <- c(
nrow(difference.v15.Normalised),
nrow(difference.v30.Normalised),
nrow(difference.v60.Normalised),
nrow(difference.v90.Normalised),
nrow(difference.Exon.Normalised),
nrow(difference.Domain.Normalised),
nrow(difference.v15.Unnormalised),
nrow(difference.v30.Unnormalised),
nrow(difference.v60.Unnormalised),
nrow(difference.v90.Unnormalised),
nrow(difference.Exon.Unnormalised),
nrow(difference.Domain.Unnormalised))

variant.frequency.data <- matrix(typecounts, nrow = 2, byrow=TRUE, dimnames = list(c("Unnormalised","Normalised"), c("+/- 15bp","+/- 30bp","+/- 60bp","+/- 90bp","Exon","Domain")))
png("variant-counts.png",width=1080,height=720)
barplot(variant.frequency.data,
        xlab = "Type",
        ylab = "Variants",
        col = c("red","green"),
        beside = TRUE,
        cex.lab = 1.5
)
legend("topright",
       c("Normalised","Unnormalised"),
       fill = c("red","green")
)
dev.off()

#Find BRCA1 variants where normalised local constraint makes a difference
difference.brca1.v15.Normalised <- subset(variants.data.frame, Type == "V15bp" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09 & Gene == "BRCA1", select=c(Key, VariantName))
difference.brca1.v30.Normalised <- subset(variants.data.frame, Type == "V30bp" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09 & Gene == "BRCA1", select=c(Key, VariantName))
difference.brca1.v60.Normalised <- subset(variants.data.frame, Type == "V60bp" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09 & Gene == "BRCA1", select=c(Key, VariantName))
difference.brca1.v90.Normalised <- subset(variants.data.frame, Type == "V90bp" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09 & Gene == "BRCA1", select=c(Key, VariantName))
difference.brca1.Exon.Normalised <- subset(variants.data.frame, Type == "VExon" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09 & Gene == "BRCA1", select=c(Key, VariantName))
difference.brca1.Domain.Normalised <- subset(variants.data.frame, Type == "VDomain" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09 & Gene == "BRCA1", select=c(Key, VariantName))
difference.brca1.Domain.Normalised <- unique(difference.brca1.Domain.Normalised)

#Calculate the percentage of variants where local constraint makes a difference
get.percent.make.difference.gene <- function(gene,variants)
  {
  all.variants.gene <- nrow(unique(subset(variants, Gene == gene & Type != "S15bp" & Type != "S30bp" & Type != "S60bp" & Type != "S90bp" & Type != "SExon" & Type != "SDomain", select=c(Key, VariantName))))
  difference.variants.gene <- nrow(unique(subset(variants, Gene == gene & Type != "S15bp" & Type != "S30bp" & Type != "S60bp" & Type != "S90bp" & Type != "SExon" & Type != "SDomain" & Gene_Constraint < 3.09 & Normalised_Constraint > 3.09, select=c(Key, VariantName))))
  percent.difference.gene <- (difference.variants.gene/all.variants.gene)*100
  #Compile metrics into vector
  output <- percent.difference.gene
  return(output)
  }
unique.genes <- unique(subset(variants.data.frame,, select=c(Gene)))
percentage.differences <- apply(unique.genes,1,get.percent.make.difference.gene,variants.data.frame)
genes.variant.difference.percentage <- data.frame(Gene=unique.genes,PercentDiff=percentage.differences)

#Data to display on bar plot
display.genes.variant.difference.percentage <- subset(genes.variant.difference.percentage, PercentDiff > 0, select=c(Gene, PercentDiff))
#Colors of bars based on the panel
display.genes.variant.difference.percentage$Color <- c(
  "red",
  "blue",
  "blue",
  "blue",
  "yellow",
  "red",
  "yellow",
  "blue",
  "blue",
  "yellow",
  "yellow",
  "yellow",
  "yellow",
  "yellow",
  "yellow",
  "blue",
  "yellow",
  "yellow",
  "yellow",
  "yellow",
  "green",
  "green",
  "orange",
  "orange",
  "yellow",
  "yellow",
  "yellow",
  "purple",
  "green",
  "blue",
  "green",
  "yellow",
  "yellow",
  "yellow",
  "yellow")
display.genes.variant.difference.percentage <- display.genes.variant.difference.percentage[order(display.genes.variant.difference.percentage$PercentDiff),]

#Display chart of % of variants in each gene where local constraint makes a difference
png("gene-make-difference.png",width=1080,height=720)
par(mar=c(4,7.5,2,1))
barplot(
  display.genes.variant.difference.percentage$PercentDiff, 
  names.arg = display.genes.variant.difference.percentage$Gene,
  col = display.genes.variant.difference.percentage$Color,
  horiz = TRUE,
  xlim = c(0,100),
  ann = FALSE,
  las = 1,
  xaxs = "i",
  yaxs = "i")
legend("right",
       c("Breast Cancer","Aortopathy","DSD","Noonan","Lynch","Ovarian Cancer"),
       fill = c("red","blue","yellow","green","orange","purple")
)
mtext(side = 1, text = "Percentage of Variants", line = 2, cex = 1.5)
mtext(side = 2, text = "Gene", line = 6, cex = 1.5)
dev.off()

#Count the number not displayed for the figure legend
not.displayed <- nrow(genes.variant.difference.percentage)-nrow(display.genes.variant.difference.percentage)
not.displayed.percent <- (not.displayed/nrow(genes.variant.difference.percentage))*100
text.for.not.display <- c(round(not.displayed.percent, digits = 1),"% (",not.displayed,"/",nrow(genes.variant.difference.percentage),")"," of the genes are not shown here, as they had no variants that were locally constrained, but not constrained gene wide.")
write(paste(text.for.not.display,collapse=""),"genes-no-difference.txt")

#Make U scores graph
uscores.data.frame <- read.table("uscores.txt", header=TRUE)
subset.uscores.data.frame <- subset(uscores.data.frame,Nucleotide >= 600 & Nucleotide <= 800, select=c(Nucleotide, UScore))
png("uscores.png",width=1080,height=720)
par(mar=c(4,4,4,4))
plot(uscores.data.frame$Nucleotide,uscores.data.frame$UScore,type="h",col="blue",
     ann = FALSE,
     xaxs = "i",
     yaxs = "i")
lines(subset.uscores.data.frame$Nucleotide,subset.uscores.data.frame$UScore,type="h",col="red")
mtext(side = 1, text = "Nucleotide in BRCA1", line = 3, cex = 1.5)
mtext(side = 2, text = "U Score", line = 2, cex = 1.5)
dev.off()

#Calculate number of variants expected in the subset and return it to text file
total.uscore <- sum(uscores.data.frame$UScore)
subset.uscore <- sum(subset.uscores.data.frame$UScore)
percentage.uscore <- (subset.uscore/total.uscore)*100
total.expected <- 508.9
subset.expected <- (total.expected/100)*percentage.uscore
write(paste(c("as the percentage of the area under the curve found between nucleotides 600 and 800 is ",round(percentage.uscore,digits=2),"%, and the total expected number of variants in ExAc is ",total.expected,", the number of variants expected to fall in the region between 600 and 800 is ",round(subset.expected,digits=1)),collapse=""),"uscoreexplain.txt")

#Boxplot of local constraint variability
boxplot.data.frame <- subset(variants.data.frame,Type == "V15bp" | Type == "V30bp" | Type == "V60bp" | Type == "V90bp" | Type == "VExon" | Type == "VDomain",select=c(Gene,Raw_Constraint,Normalised_Constraint))
genes.constraint <- unique(subset(variants.data.frame,,c(Gene,Gene_Constraint)))
boxplot.lims = c(-10,10)
png("boxwhisker-normalised.png",width=1080,height=720)
par(mar=c(8,4,1,4))
boxplot(boxplot.data.frame$Normalised_Constraint~boxplot.data.frame$Gene,
        outline=FALSE,
        ylim = boxplot.lims,
        xaxs = "i",
        yaxs = "i",
        ann = FALSE,
        las = 2,
        col="yellow",
        medcol="red"
        )
mtext(side = 1, text = "Gene", line = 6, cex = 1.5)
mtext(side = 2, text = "Local Constraint (Normalised)", line = 2, cex = 1.5)
abline(h=0)
abline(h=3.09)
par(new=T)
plot(genes.constraint$Gene,genes.constraint$Gene_Constraint,
     ylim = boxplot.lims,
     xaxs = "i",
     yaxs = "i",
     xaxt = "n",
     yaxt = "n"
)
dev.off()

boxplot.lims = c(-5,5)
png("boxwhisker-raw.png",width=1080,height=720)
par(mar=c(8,4,1,4))
boxplot(boxplot.data.frame$Raw_Constraint~boxplot.data.frame$Gene,
        outline=FALSE,
        ylim = boxplot.lims,
        xaxs = "i",
        yaxs = "i",
        #xlab = "Gene",
        #ylab = "Local Constraint (Normalized)",
        ann = FALSE,
        las = 2,
        col="yellow",
        medcol="red"
)
mtext(side = 1, text = "Gene", line = 6, cex = 1.5)
mtext(side = 2, text = "Local Constraint (Normalised)", line = 2, cex = 1.5)
abline(h=0)
abline(h=3.09)
par(new=T)
plot(genes.constraint$Gene,genes.constraint$Gene_Constraint,
     ylim = boxplot.lims,
     xaxs = "i",
     yaxs = "i",
     xaxt = "n",
     yaxt = "n"
)
dev.off()
