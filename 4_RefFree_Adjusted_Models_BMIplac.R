#################################################################################
# Genome-wide DNA methylation scan using RefFreeEWAS to adjust for cellular heterogeneity
# RefFreeEWAS Houseman et al. 2016 (PMID: 27358049)
# Code provided by Todd Everson - adapted by N. Fernandez-Jimenez
################################################################################

setwd("C://Users//TEVERSO//Desktop//P_PACE Projects//Maternal Smoking Placental DNAM//Analysis//Todd Test")
#################################################################################
# Step 0: Load Packages, Functions and Data (Data need to be normalized, and standardized across probe-types [BMIQ + FunNorm is ideal])

## For RefFreeEWAS analyses:
library(readr)
library(RefFreeEWAS)

## For Categorical Splitting
library(Hmisc)

## For RLM Linear models:
library(data.table)	# to process results
library(MASS) 		# rlm function for robust linear regression
library(sandwich) 	# Huber?s estimation of the standard error
library(lmtest) 		# to use coeftest
library(parallel) 	# to use multicore approach - part of base R

# load methylation data as fully processed beta values
load(file="3_ComBat_Adjusted.RData")

# load Chen Probe annotation file and exclude SNP and Cross Hybridizing probes:
load(file="ChenProbeIDs.RData") #### Data from Chen et al (2013) PMID: 23314698
annot2$SNPs = annot2[,"Global"] #### Can change Global to "AFR", "AMR", "ASN", or "EUR" to match content specific allelic frequency; Use "Global" if population is very diverse
kept.probes = subset(annot2, as.character(Name) %in% rownames(fnbetas.adj) & sex %in% "Keep" & CRSH %in% "Keep" & SNPs %in% "Keep")

# Drop problematic probes:
beta_matrix = fnbetas.adj[which(rownames(fnbetas.adj) %in% as.character(kept.probes$Name)),]

# load dataset of exposure variables and covariates, this object should be named 'pheno_data':
load(file="NHBCS_SMK_EWAS.RData")	## Rownames of pheno_data need be same IDs as Colnames of beta_matrix ##

# make sure that exposure/covariate data are being treated in the correct way by R 
# Factors: Maternal Smoking; Parity; Maternal Education
# Continuous: Maternal Age

# Ensure columns from methylation beta-value matrix align with phenotype data
table(ifelse(rownames(pheno_data)==colnames(beta_matrix),"Matched","--NOT MATCHED--"))
# if any "NOT MATCHED", check pheno_data file and re-order to match IDs for beta_matrix
PHENO = pheno_data

#### !!! for betas object, ensure that rows are probes and columns are samples

#################################################################################
# Step 1: Subset methylation matrix by variance, estimate K (number of cell types)

## Calculated variance for each CpG site, only using most variable 10,000 sites:
v 			= apply(beta_matrix,1,var,na.rm=TRUE) 
y.tops 		= beta_matrix[names(v[order(v,decreasing=TRUE)][1:10000]),]

# Get an array of solutions (for different values of K), using 50 iterations per solution:
cellmixArray	= RefFreeCellMixArray(y.tops,Klist=1:12,iters=50)

# Do the bootstrap procedure (Houseman Supplemental S3)
cellmixArrayBoot	= RefFreeCellMixArrayDevianceBoots(
      			cellmixArray,            # Array object
      			Y=y.tops,                # Data on which array was based
      			R=100,                   # Number of bootstraps
      			bootstrapIterations=10)  # Number of iterations per boot

# Show winsorized mean deviance per K
wnsrMeanDev		= apply(cellmixArrayBoot[-1,], 2, mean, trim=0.25)
wnsrMeanDev

# Choose K based on minimum deviance
Kchoose 		= which.min(wnsrMeanDev)
Kchoose 		# estimated number of cell-types

#################################################################################
# Step 2: Estimated cell-type proportions per sample (Omega) and methylation profiles of estimated cell-types (Mu)

# PART A: Use above K with minimum deviance to generate estimated cell-mix (Omega)
Omega 		= cellmixArray[[ Kchoose ]]$Omega # The cell mixture matrix

pdf("PACE_STUDYNAME_Omega_heatmap.pdf")
heatmap(Omega)	# Visualize cell-mix across samples
dev.off()

pdf("PACE_STUDYNAME_Omega_OutlierScreening.pdf")
par(mfrow=c(2,round(Kchoose/2)))
for(cell in 2:Kchoose){
plot(Omega[,1],Omega[,cell],xlab="Omega-1",ylab=paste("Omega-",cell,sep=""))
}
dev.off()

#### !!!! Examine OutlierScreening plots, if obvious extreme outliers are generated in these plots, repeate Step 2, Part A, by
#### !!!! Replacing 	Omega = cellmixArray[[ Kchoose ]]$Omega 		with 		Omega = cellmixArray[[ (Kchoose-1) ]]$Omega 
#### !!!! Can repeat this process (Kchoose-2), (Kchoose-3), etc. until no more outliers are being generated
#### !!!! Once Omega doesn't contain extreme outliers, proceed with Part B

# Part B: Mu for entire CpG array (contribution of individual CpGs to the cell-mixture estimation) for the entire array:
Mu 			= RefFreeCellMix(
      			Y=y.tops,  					# Data on which the solution was based 
      			mu0=cellmixArray[[ Kchoose ]]$Mu, 	# Mu from the solution
      			Yfinal=beta_matrix,      		# The whole array
      			iters=1)$Mu  				# One iteration is enough here
dim(Mu) 

# Row variances to identify the most highly variable CpGs by Cell-Mix:
Mu.Row.Var 			= as.data.frame(as.matrix(apply(Mu,1,var,na.rm=TRUE)))
Mu.Row.Var$MVar		= as.numeric(as.factor(cut2(Mu.Row.Var[,1],g=20,levels.mean=TRUE)))
Mu.Row.Var			= Mu.Row.Var[order(Mu.Row.Var[,1]*(-1)),]
colnames(Mu.Row.Var) 	= c("Var","RankVar")

write.table(Mu.Row.Var, "PACE_STUDYNAME_Variability_By_CellMix_CpGs.txt", na="NA")

#################################################################################
# Step 3: See if cell-mixtures are associated with phenotypes/covariates of interest:

## Create indicator for which type of table to create:
table.type = ifelse(lapply(PHENO, class) %in% c("numeric", "integer"),"Continuous","Categorical")

## Create Table of rank-sum associations with categorical variables:
Cat.Table.p = data.frame(matrix(NA,nrow=table(table.type)[1],ncol=(dim(Omega)[2])))
for(tc in seq(1:dim(Omega)[2])){
Cat.Table.p[,c(tc)] =	t(as.data.frame(
				apply(as.data.frame(PHENO[,which(table.type %in% "Categorical")]),2,function(kt) 
				c(kruskal.test(Omega[,tc],as.factor(kt))[c("p.value")]))
		)	)
rownames(Cat.Table.p) = names(PHENO[,which(table.type %in% "Categorical")])
}
names(Cat.Table.p) = c(paste("KW_Pval_",colnames(Omega),sep=""))

## Create Table of rank-sum associations with continuous variables (Kendall for duplciate values):
Cont.Table.p = data.frame(matrix(NA,nrow=table(table.type)[2],ncol=(dim(Omega)[2])))
for(tc in seq(1:dim(Omega)[2])){
Cont.Table.p[,c(tc)] =	t(as.data.frame(
				apply(as.data.frame(PHENO[,which(table.type %in% "Continuous")]),2,function(kt) 
				c(cor.test(Omega[,tc],c(as.numeric(kt)),method="kendall")[c("p.value")]))
		)	)
rownames(Cont.Table.p) = names(PHENO[,which(table.type %in% "Continuous")])
}
names(Cont.Table.p) = c(paste("Tau_Pval_",colnames(Omega),sep=""))

write.csv(Cat.Table.p, "PACE_STUDYNAME_CatCovariate_CellMix_Associations.csv", na="NA")
write.csv(Cont.Table.p, "PACE_STUDYNAME_ContCovariate_CellMix_Associations.csv", na="NA")

# Remove 1 cell-type to reduce multi-colinearity
Omega.s = Omega[ ,1:(length(colnames(Omega))-1)]

# Ensure CellMix matrix matches methylation beta-value matrix and other phenotype data
table(ifelse(rownames(Omega.s)==rownames(PHENO),"MATCHED","Not Matched"))
table(ifelse(rownames(Omega.s)==colnames(beta_matrix),"MATCHED","Not Matched"))

#################################################################################

# Step 4: R-code for detecting outliers using gaphunter
# prepared by: Sinjini Sikdar and adapted by: N Fernandez-Jimenez - nora.fernandez@ehu.eus (16/03/2018)
# Line 74 - modified by MKL(mikyeong.lee@nih.gov) to keep the digits after the decimal ## Dec 12, 2017 
## A final outlier-removed beta file ("Beta4ewas.Rda")

# Required input & parameters
N <- ncol(beta_matrix) ## No. of samples in your data.
gapsize <- 0.3    ## Gap size between the outliers and all other values. We have used 0.3 in our analysis.
cutoff <- ifelse(5>(0.25/100)*N,5/N,0.0025)   ## This cutoff is chosen for detecting probes with outliers. We have chosen this
## cutoff such that a probe can have a maximum of 5 or 0.0025 of the total number                                                
## of samples (whichever is larger) as outliers. Can change it if required.

## Note: If gap hunter returns error as it may not detect any outlier, please reduce the gap size or increase the cutoff.

# Program for gap-hunter
## log files for identifying missing values (NAs) before running gaphunter
probes_log_before <- data.frame(probes=rownames(beta_matrix),no_missing_before=apply(beta_matrix, 1, function(x) sum(is.na(x))))
samples_log_before <- data.frame(samples=colnames(beta_matrix),no_missing_before=apply(beta_matrix, 2, function(x) sum(is.na(x))))

## imputation (as gaphunter cannot handle missing values)
beta1 <- t(apply(beta_matrix,1,function(x) ifelse(is.na(x),median(x,na.rm=T),x)))          

## gaphunter
gapres <- gaphunter(beta1,keepOutliers=TRUE,threshold=gapsize,outCutoff=cutoff)
gapres1 <- gaphunter(beta1,keepOutliers=FALSE,threshold=gapsize,outCutoff=cutoff)
all_signals_probes <- gapres$proberesults
all_signals_samples <- gapres$sampleresults
all_signals_all <- merge(all_signals_probes,all_signals_samples,by.x="row.names",by.y="row.names")
without_outlier_signals_probes <- gapres1$proberesults
outlier <- setdiff(rownames(all_signals_probes),rownames(without_outlier_signals_probes))  ## probes with outliers
outlier_signals_all <- all_signals_all[which(all_signals_all[,1]%in%outlier),]          
no_incl <- max(outlier_signals_all$Groups)+2
colnames(outlier_signals_all) <- c("probes",colnames(outlier_signals_all)[2:no_incl],colnames(beta1))

new_beta <- NULL
for(p in 1:nrow(outlier_signals_all)){
  group_number <- as.numeric(strsplit(names(which.max(outlier_signals_all[p,3:no_incl])),"")[[1]][6])
  df2 <- outlier_signals_all[p,c(rep(TRUE,no_incl),outlier_signals_all[p,-c(1:no_incl)]!=group_number)]
  sub_int <- colnames(df2)[-c(1:no_incl)]
  beta_out <- beta_matrix[which(rownames(beta_matrix)%in%df2$probes),]
  beta_out[names(beta_out)%in%sub_int] <- NA
  new_beta <- rbind(new_beta,c(probes=df2$probes,beta_out))
}
rownames(new_beta) <- new_beta[,1]
new_beta_1 <- data.frame(new_beta[,-1],check.names=FALSE)
beta2 <- as.data.frame(beta_matrix)
beta2[match(rownames(new_beta_1),rownames(beta2)), ] <- new_beta_1

# Outputs
## Final beta matrix (saved as beta2 in .Rda format) with all outliers detected as NAs. 
## We have the same beta matrix as before (input matrix) with additional NAs for extreme outliers.
beta_matrix[is.na(new_beta_1)] <- NA; save(beta_matrix, file = "Beta4ewas.Rda") ## Use this for further analysis
load("Beta4ewas.Rda")
## log files in .csv formats. These log files will have columns with:
## probe/sample names: "probes"/"samples", 
## number of missing values (NAs) before running gaphunter: "no_missing_before", 
## number of missing values (NAs including outliers) after running gaphunter: "no_missing after",
## number of outliers identified by gaphunter: "outlier_gaphunter", and 
## the percentage of (non-missing) samples/probes identified as outliers by gaphunter: "percent_outlier_gaphunter".

probes_log_after <- data.frame(probes=rownames(beta2),no_missing_after=apply(beta2, 1, function(x) sum(is.na(x))))
log_probes_combined <- merge(probes_log_before,probes_log_after,by.x="probes",by.y="probes",all.x=TRUE,all.y=TRUE)
log_probes_combined$outlier_gaphunter <- log_probes_combined$no_missing_after-log_probes_combined$no_missing_before
log_probes_combined$percent_outlier_gaphunter <- (log_probes_combined$outlier_gaphunter/(ncol(beta2)-log_probes_combined$no_missing_before))*100
write.csv(log_probes_combined,"Gaphunter_CpGs_STUDY_DATE_ANALYST.csv")

samples_log_after <- data.frame(samples=colnames(beta2),no_missing_after=apply(beta2, 2, function(x) sum(is.na(x))))
log_samples_combined <- merge(samples_log_before,samples_log_after,by.x="samples",by.y="samples",all.x=TRUE,all.y=TRUE)
log_samples_combined$outlier_gaphunter <- log_samples_combined$no_missing_after-log_samples_combined$no_missing_before
log_samples_combined$percent_outlier_gaphunter <- (log_samples_combined$outlier_gaphunter/(nrow(beta2)-log_samples_combined$no_missing_before))*100
write.csv(log_samples_combined,"Gaphunter_samples_STUDY_DATE_ANALYST.csv")   

# transpose betas so that rows are samples and columns are probes
beta_matrix = t(beta_matrix)

# Ensure CellMix matrix matches methylation beta-value matrix and other phenotype data
table(ifelse(rownames(Omega.s)==rownames(PHENO),"MATCHED","Not Matched"))
table(ifelse(rownames(Omega.s)==rownames(beta_matrix),"MATCHED","Not Matched"))

#################################################################################
# Step 5: Conduct EWAS Models
# Standard code from PACE

#### Modeling with covariate adjustment
#Add function for running the model (4 covariates in addition to cell-mix)
RLMtest = function(meth_matrix,methcol,exposure, X1, X2, X3, X4) {
  mod = try(rlm(meth_matrix[, methcol]~exposure+X1+X2+X3+X4,maxit=200))
  if(class(mod) == "try-error"){
    print(paste("error thrown by column", methcol))
    invisible(rep(NA, 3))
  }else cf = coeftest(mod, vcov=vcovHC(mod, type="HC0"))
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
}

#Run adjusted MWAS
system.time(ind.res1 <- mclapply(setNames(seq_len(ncol(beta_matrix)), dimnames(beta_matrix)[[2]]), RLMtest, meth_matrix=beta_matrix, 
                                 exposure=PHENO$mbmi
                                 , X1=PHENO$mage	# Maternal Age
                                 , X2=PHENO$parity2c	# Parity
                                 , X3=PHENO$meducation3c	# Maternal Education
                                 , X4=PHENO$smokpreg	# Smoking in pregnancy
))
Sys.time()

#### Modeling with covariate and cell-mix adjustment
#Add function for running the model (4 covariates in addition to cell-mix)
RLMtest = function(meth_matrix,methcol,exposure, X1, X2, X3, X4, Cells) {
  mod = try(rlm(meth_matrix[, methcol]~exposure+X1+X2+X3+X4+Cells,maxit=200))
  if(class(mod) == "try-error"){
    print(paste("error thrown by column", methcol))
    invisible(rep(NA, 3))
  }else cf = coeftest(mod, vcov=vcovHC(mod, type="HC0"))
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
}

#Run adjusted MWAS
system.time(ind.res2 <- mclapply(setNames(seq_len(ncol(beta_matrix)), dimnames(beta_matrix)[[2]]), RLMtest, meth_matrix=beta_matrix, 
                                 exposure=PHENO$mbmi
                                 , X1=PHENO$mage	# Maternal Age
                                 , X2=PHENO$parity2c	# Parity
                                 , X3=PHENO$meducation3c	# Maternal Education
                                 , X4=PHENO$smokpreg	# Smoking in Pregnancy
                                 , Cells = Omega.s))
Sys.time()


#Process covariate adjusted results
setattr(ind.res1, 'class', 'data.frame')
setattr(ind.res1, "row.names", c(NA_integer_,4))
setattr(ind.res1, "names", make.names(names(ind.res1), unique=TRUE))
probelistnames1 = names(ind.res1)
all.results1 = t(data.table(ind.res1))
all.results1 = data.table(all.results1)
all.results1[, probeID := probelistnames1]
setnames(all.results1, c("BETA","SE", "P_VAL", "probeID")) # rename columns
setcolorder(all.results1, c("probeID","BETA","SE", "P_VAL"))
rm(probelistnames1, ind.res1)

#Lambda calculation
lambda1 = qchisq(median(all.results1$P_VAL,na.rm=T), df = 1, lower.tail = F)/
  qchisq(0.5, 1)
## Save analysis sample sizes & lambda:
N1 = nrow(beta_matrix)

#Q-Q plots
library(qqman)
pdf("PACE_STUDYNAME_Model1_QQ_Plot.pdf")
qq(all.results1$P_VAL,main="PACE_STUDYNAME_QQ-plot")
dev.off()

## For zipping files:
library(R.utils)

# export table of results
write.table(all.results1, "PACE_STUDYNAME_Model1_date.txt",na="NA")
gzip("PACE_STUDYNAME_Model1_date.txt")

#Process covariate and reffree adjsuted results
setattr(ind.res2, 'class', 'data.frame')
setattr(ind.res2, "row.names", c(NA_integer_,4))
setattr(ind.res2, "names", make.names(names(ind.res2), unique=TRUE))
probelistnames2 = names(ind.res2)
all.results2 = t(data.table(ind.res2))
all.results2 = data.table(all.results2)
all.results2[, probeID := probelistnames2]
setnames(all.results2, c("BETA","SE", "P_VAL", "probeID")) # rename columns
setcolorder(all.results2, c("probeID","BETA","SE", "P_VAL"))
rm(probelistnames2, ind.res2)

#Lambda calculation
lambda2 = qchisq(median(all.results2$P_VAL,na.rm=T), df = 1, lower.tail = F)/
  qchisq(0.5, 1)
## Save analysis sample sizes & lambda:
N2 = nrow(beta_matrix)

#Q-Q plots
pdf("PACE_STUDYNAME_Model2_QQ_Plot.pdf")
qq(all.results2$P_VAL,main="PACE_STUDYNAME_QQ-plot")
dev.off()

# export table of results
write.table(all.results2, "PACE_STUDYNAME_Model2_date.txt",na="NA")
gzip("PACE_STUDYNAME_Model2_date.txt")

# Export lambda for each
write.table(data.frame(Model = c("Covariate Adjusted","Covariate and RefFree Adjusted"), Lambda = c(lambda1,lambda2), N=c(N1,N2)), "PACE_STUDYNAME_Lambdas.txt", na="NA")



