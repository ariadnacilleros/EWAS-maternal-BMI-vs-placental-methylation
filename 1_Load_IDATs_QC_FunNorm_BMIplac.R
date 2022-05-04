##############################################################################################
#### (Pre) Load relevant packages and functions ##############################################
##############################################################################################

## Load relevant packages:
library(minfi)
library(RColorBrewer)	

## Set working directory, must be where iDATs and Illumina sample sheet are located:
setwd("C://Users//TEVERSO//Desktop//D_NHBCS//S_Placenta_450KDNAM_Array//Preprocessed with BMIQ - Copy")

##############################################################################################
#### Read in iDATs, perform basic QC (detection p-values), and functional normalization ######
##############################################################################################

#### Read in sheet with batch and ID info:
targets 	= read.metharray.sheet(getwd())
RGSet	= read.metharray.exp(targets = targets)
#Features  Samples 
#  622399      346 
save(RGSet, file = "RGChannelSet_IDAT_Read_in.RData")

#### Generate batch file
pd <- pData(RGSet)

#### Get raw beta values:
Mset = preprocessRaw(RGSet)
RatioSet = ratioConvert(Mset, what = "both", keepCN = TRUE)
RawBetas = getBeta(RatioSet)
qc = getQC(Mset)
png(file="Sample_QC_Plot.png",width=1000,height=1000,pointsize=20)
plotQC(qc)
dev.off()

## Functional Normalization (This step takes a while and requires lots of memory):
Mset.norm = preprocessFunnorm(RGSet, sex = NULL, bgCorr = T, dyeCorr = T, nPCs = 2, verbose = TRUE)
FunNormBetas = getBeta(Mset.norm)

## Check whether reported sex matches inferred sex ()
pred.sex = data.frame(getSex(Mset.norm, cutoff = -2))
plotSex(getSex(Mset.norm, cutoff = -2))
pred.sex$ID = as.character(pd$Sample_Name)
## Can use below file to check for agreement with reported sex:
write.table(pred.sex, file = "PredictedSex.csv",quote = TRUE, sep = ",", row.names = TRUE, col.names = TRUE, qmethod = "double")

## Check distributions before and after normalization
library(RColorBrewer)	
png(file="BetaValue_Distributions_BeforeQC.png",width=1400,height=700,pointsize=12)
par(mfrow=c(1,2))
densityPlot(RawBetas, sampGroups = pd$Slide, legend=FALSE, main = "Raw Betas", xlab = "Beta")
densityPlot(FunNormBetas, sampGroups = pd$Slide, legend=FALSE, main = "FunNorm Adjusted Betas", xlab = "Beta")
dev.off()

## Use detection p-values to exclude samples:
detect.p = detectionP(RGSet, type = "m+u")
binning <- detect.p
binning[binning < 0.01] = 0
binning[binning > 0.01] = 1

## Identify potential problem samples (poor det p-values for > 1% of probes)
prop.table = NULL
  for (i in 1:length(colnames(detect.p))){
      table = table(binning[,i])
      prop.table.new = prop.table(table)
      prop.table = rbind(prop.table, prop.table.new)
      }
rownames(prop.table) = colnames(binning)
prop.table = as.data.frame(prop.table)
colnames(prop.table) = c("Pass", "Fail")
ptable = prop.table[order(prop.table[,2],decreasing = FALSE),]

## Review prop.table to determine which samples should be dropped
failed = rownames(subset(ptable, Fail > 0.01)) 
length(failed)	## If no samples fail (length(fail) = 0), then can skip these next two lines:
RGSet = RGSet[,-c(which(rownames(pd) %in% failed))] ## Remove samples with poor detection, then re-run FunNorm
pd = pData(RGSet)

## Identify and drop individual probes with at least one probe with poor detection:
detect.p = detectionP(RGSet, type = "m+u")  	
detect.p[detect.p > 0.01] = NA
filtered.p = na.omit(detect.p) 
probeQC = data.frame(SamplesExcluded = ifelse(length(failed) > 0,length(failed),"None"),
				ProbesExcluded = (dim(detect.p)[1]-dim(filtered.p)[1]),
				RemainingProbes = length(rownames(filtered.p)))
write.table(probeQC, file = "ProbeQC.csv",quote = TRUE, sep = ",", row.names = TRUE, col.names = TRUE, qmethod = "double")

## Re-run FunNorm if samples were dropped:
Mset.norm = preprocessFunnorm(RGSet, sex = NULL, bgCorr = T, dyeCorr = T, nPCs = 2, verbose = TRUE)
FunNormBetas = getBeta(Mset.norm)
## Exclude bad probes from the RawBeta and FunNormBeta data:
intersect = intersect(rownames(Mset.norm), rownames(filtered.p))
length(intersect)
Mset.fn.sub = Mset.norm[intersect,]
FunNormBetas.sub = getBeta(Mset.fn.sub)

## Check density plots after excluding the poorly-detected probes:
png(file="BetaValue_Distributions_AfterQC_FunNorm.png",width=700,height=700,pointsize=12)
densityPlot(FunNormBetas.sub, sampGroups = pd$Slide, legend=FALSE, main = "PostQC - Normalized Beta", xlab = "Beta")
dev.off()

## Get annotations probes:
annot = getAnnotation(RGSet) 

## The 'FunNormBetas.sub' object includes data for samples and probes with good detection p-values, and that have undergone functional normalization:
save(annot,pd,FunNormBetas.sub,Mset.fn.sub, file = "1_BasicQC_and_FunctionalNormalization.RData")

