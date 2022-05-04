##############################################################################################
#### (Pre) Load relevant packages and functions ##############################################
##############################################################################################
setwd("C://Users//TEVERSO//Desktop//D_NHBCS//S_Placenta_450KDNAM_Array//Preprocessed with BMIQ - Copy")
## Load Libraries and Data
library(minfi)
library(sva)
load("2_BMIQ_Adjsuted.RData")

##############################################################################################
#### (3) Use ComBat to remove batch effect ###################################################
##############################################################################################

## Read in a function that will calculate the variance of each row
rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE,
                    SumSquares=FALSE, twopass=FALSE) {
   if (SumSquares) return(rowSums(x^2, na.rm, dims))
   N <- rowSums(!is.na(x), FALSE, dims)
   Nm1 <- if (unbiased) N-1 else N
   if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
                      sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
   (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
 }

mval <- apply(FunNorm.BMIQ, 2, function(x) log2((x)/(1-x)))

## Calculate the variance of each probe and remove any with a variance of 0 prior to Combat.
vars = as.matrix(rowVars(mval))

## Replace all probes with no variance with NA and remove them from the FunNorm set
vars[vars == 0] = NA
vars = na.omit(vars)
intersect = intersect(rownames(vars), rownames(mval))
print(length(intersect))
fn.sub = FunNorm.BMIQ[intersect,]
mval = mval[intersect,]

## Ensure Objects are aligned
table(ifelse(rownames(pd) == colnames(mval),"Match","Off")) ## All should Match

## Group small batches with others:
pd$Batch <- NA
pd$Batch <- pd$Sample_Plate
batch <- pd$Batch

###########################################################################
## Check variation in array data associated with batch (ie. Slide/plate/box)
###########################################################################

## Run a principle component analysis to determine if there are any remaining
## batch effects following data normalization.

PCobj = prcomp(t(mval), retx = T, center = T, scale. = T)
boxplot(PCobj$x,col="grey",frame=F)

# Can use Scree plot to determine number of PCs to keep
plot(PCobj,type="line",cex.lab=1.5, cex.main=1.5) 
abline(v=5,lty=3, col="red")	# Variance being to flatten at 5th PC

# Extract the PCs from the PCobj object
PCs = PCobj$x

# Extract the proportion of variability and cummulative proportion of 
# varibility explained by the top R PCs.
	
R = 5
propvar = summary(PCobj)$importance["Proportion of Variance", 1:R]
cummvar = summary(PCobj)$importance["Cumulative Proportion", 1:R]
		
## Generate plots of the resulting PCs
	
# Plot of the proportion of variability explained by the top R PCs
# Plot of the cummulative proportion of variability explained by the top R PCs

par(mfrow=c(1,2))	
par(mar = c(5,5,4,2))
barplot(propvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Variation Explained (%)", cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
par(mar = c(5,5,4,2))
barplot(cummvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Cummulative Variation Explained (%)",cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
	
# Plot of PCX and PCY; by Batch
PCs = PCobj$x
PCs =PCs[,1:5]
Prin.comp<-merge(PCs,pd, by = "row.names",all=T) 

png(file="PC_Variation_by_Plate.png",width=900,height=900,pointsize=12)
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Batch), xlab = "PC1", ylab = "PC2")
legend("topright", legend=levels(as.factor(Prin.comp$Batch)), text.col=seq_along(levels(as.factor(Prin.comp$Batch))))
dev.off()
## Can further test via ANOVA, to look for strong variations in PCs by Batch:
anova(lm(Prin.comp$PC1~Prin.comp$Batch))
anova(lm(Prin.comp$PC2~Prin.comp$Batch))
anova(lm(Prin.comp$PC3~Prin.comp$Batch))
anova(lm(Prin.comp$PC4~Prin.comp$Batch))
anova(lm(Prin.comp$PC5~Prin.comp$Batch))

## Model matrix for batch-corrections (May need to adjust model matrix to 'protect' coefficients (study specific)):
mod <- model.matrix(~1, data=pd)

## Run ComBat to remove slide # batch effects

combat.adj = ComBat(mval,batch = batch, mod = mod)

## Check to see if batch effect was succesfully removed
PCobj = prcomp(t(combat.adj), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs =PCs[,1:5]
Prin.comp<-merge(PCs,pd, by = "row.names",all=T) 

#### Check whether batches are still distinguished by first and second PC:
png(file="PC_Variation_by_Batch_AfterComBat.png",width=900,height=900,pointsize=12)
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Batch), xlab = "PC1(12.2%)", ylab = "PC2(9.3%)")
legend("topright", legend=levels(as.factor(Prin.comp$Batch)), text.col=seq_along(levels(as.factor(Prin.comp$Batch))))
dev.off()


###########################################################################
## Convert the batch-adjusted M-values back into betas:

expit2 = function(x) 2^x/(1+2^x)
fnbetas.adj = expit2(combat.adj)

## Save normalized and batch-adjusted beta values
save(fnbetas.adj,pd,annot, file = "3_ComBat_Adjusted.RData")


