## DISCOVERING THE DIAGNOSTIC AND PROGNOSTIC LncRNA SIGNATURES IN CRC PATIENTS USING MACHINE LEARNING APPROACHES. ##

############################## PART 1 DIFFERENTIAL GENE EXPRESSION ANALYSIS ############################################################
# NOTE:For differential gene expression analysis we use the UCSC xena	toil recomputed data. The original TCGA-COAD project 
# only contain 41 nontumor samples, and thus less powerfull to discriminate the differentially expressed genes.
########################################################################################################################################

# 1.Download the required data files.

# url="https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz"
# url="https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_fpkm.gz"
# url="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/gencode.v23.long_noncoding_RNAs.gtf.gz"
# url="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.long_noncoding_RNAs.gtf.gz"

# 2.Processing the long noncoding RNAs annotation file, and checking the sample and lncRNA information. 
## cd /mnt/c/Users/woodhaha/Desktop/CRC_data_mining
## echo|awk 'BEGIN{FS==" " ; OFS==" "}{if($3=="gene" ){print $10,$12,$14}}' gencode.v23.long_noncoding_RNAs.gtf> "lnc_RNAs_annotation_v22.txt"


## NOTE: R will automatically change the character "-" in gene symbol to character ".", using stringr::str_replace(X$Symbol, "\\.{1}", "-" )

# 3.Loading packages and parameter settings.
##
setwd("E:/CRC_data_mining")
#load("E:/CRC_data_mining/Results/Documents.rdata")
par(mfrow = c(1, 1), omi = c(0.5, 0.5, 0.5, 0.5), las = 1, family = "Arial")
par(cex.lab=1,cex.axis=1,cex.main=1)##default omi = c(0.5, 0.5, 0.5, 0.5), mgp = c(3, 1, 0
options(digits = 4,stringsAsFactors=F,scipen=0)

###windows(width = 7, height = 7, pointsize = 12)

LoadRpak <- function(pkg){new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
                          if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
	                      sapply(pkg, require, character.only = TRUE)}

.libPaths("C:/R Working Directory/win-library/3.6")# Using R 3.6 package library 
packages <- c("ggplot2", "dplyr", "reshape2", "RColorBrewer", "scales", "grid","regclass","caret","glmnet",
"pROC","AppliedPredictiveModeling","dplyr","ggpubr","ranger","ggplot2","ggfortify","tidyr","survminer","survival", "qdap", "table1", "xtable", "glmnet", "penalizedSVM", "parallel", "randomForestSRC",
"ggRandomForests","xtable","peperr","tgp","mlegp","pamr","lattice","c060","limma")

LoadRpak(packages)

##
Phenodata <- read.table("TcgaTargetGTEX_phenotype.txt", sep = "\t",
						header = T, stringsAsFactors = F, col.names = c("sample", "category",
						"tissue", "site", "type", "gender", "study"))

# 4.Check the mapping annotation files 
Annotation23 = data.table::fread("gencode_v23_RNAs_annotation.txt", sep = ";", header = F, data.table = F)[, 1:3]
colnames(Annotation23) = c("ENSG_id", "Type", "Symbol")
##
lncRNA_v22 <- read.table("lnc_RNAs_annotation_v22.txt", strip.white = T, blank.lines.skip = T, header = F,
					 stringsAsFactors = F, sep = ";")	[,c(1:3)]

lncRNA_v22 <- cbind(do.call(rbind, stringr::str_split(lncRNA_v22$V1, " ")), lncRNA_v22)
lncRNA_v22 <- lncRNA_v22[,-6]
colnames(lncRNA_v22) <- c("Chr","start","end","strand","ENSG_id", "Type", "Symbol")
lncRNA_v22$Position <- paste(lncRNA_v22$Chr, ":", lncRNA_v22$start, "-", lncRNA_v22$end, " ", lncRNA_v22$strand)
lncRNA_v22 <- lncRNA_v22 %>% dplyr::select(ENSG_id,Symbol, Position, Type)

##
lncRNA_v23 <- read.table("lnc_RNAs_annotation_v23.txt", strip.white = T, blank.lines.skip = T, header = F,
					 stringsAsFactors = F, sep = ";")[, c(1:3)] ## processed from (gencode.v23.long_noncoding_RNAs.gtf) file#

lncRNA_v23 <- cbind(do.call(rbind, stringr::str_split(lncRNA_v23$V1, " ")), lncRNA_v23)
lncRNA_v23 <- lncRNA_v23[, -6]
colnames(lncRNA_v23) = c("Chr", "start", "end", "strand", "ENSG_id", "Type", "Symbol")
lncRNA_v23$Position <- paste(lncRNA_v23$Chr, ":", lncRNA_v23$start, "-", lncRNA_v23$end, " ", lncRNA_v23$strand)
lncRNA_v23 <- lncRNA_v23 %>% dplyr::select(ENSG_id, Symbol, Position, Type)

##
table(duplicated(lncRNA_v23$Symbol)) ##27 lncRNA have same gene symbol
table(duplicated(lncRNA_v22$Symbol)) ##22 lncRNA have same gene symbol
table(lncRNA_v23$ENSG_id %in% lncRNA_v22$ENSG_id) ##The v22 and v23 have some discrepancies
##Ckeck the duplicated genes symbols.

(lncRNA_v23[lncRNA_v23$Symbol %in% lncRNA_v23$Symbol[duplicated(lncRNA_v23$Symbol)],])
(lncRNA_v22[lncRNA_v22$Symbol %in% lncRNA_v22$Symbol[duplicated(lncRNA_v22$Symbol)],])


# 5.Extract the tcga-gtex project colon samples.
colon.sample <- dplyr::filter(Phenodata, site == "Colon");
table(colon.sample$study) # GTEX(308) and TCGA(331)

expr <- data.table::fread("TcgaTargetGtex_rsem_gene_fpkm", sep = "\t", header = T, data.table = F) ##log2(FPKM+0.001)
#
table(expr$sample %in% Annotation23$ENSG_id) ## All genes were mapped to  gencode_v23_RNAs_annotation.
rownames(expr) <- expr$sample;
## Extracting matched TCGA-GTEX colon samples
expr <- expr[, colnames(expr) %in% colon.sample$sample] 
dim(expr) ## 60498 genes (mapped to gencode_v23_RNAs_annotation) and  639 samples.

# 6.Filter genes that with very low expression.
##Remove genes with certain percentage of value lower than a cutoff (value less than 0.1 in more than 60% of the samples).
##Because such low-expressed features tend to reflect noise and correlations based on values that are mostly zero aren't really meaningful. 
filter1 <- function(x) { length(x[x < 0.1]) / length(x) > 0.6 }
table(apply(expr, 1, filter1)) ## 15255 remained.
##
expr <- expr[apply(expr, 1, filter1) == FALSE,] #45243 genes removed.  
rownames(expr) <- gsub("-", ".", rownames(expr))
colnames(expr) <- gsub("-", ".", colnames(expr))

# 7.Count the numbers of samples in each groups.
cat(length(grep("TCGA", colnames(expr))), "TCGA samples\n");
cat(length(grep("GTEX", colnames(expr))), "GTEX samples")
group=ifelse(grepl(".11$|GTEX", colnames(expr)), "NonTumor", "Tumor")
table(group)  #NonTumor(349) and   Tumor(290) 
    
# 8.Differetial gene expression analysis in tcga-gtex using limma package.
## The data is already log2(FPKM+0.001)	transformed according to the UCSC xena	toil protocol.

design <- model.matrix(~0 + factor(group))
colnames(design) <- levels(factor(group));
head(design)
contrast_matrix <- makeContrasts(Tumor - NonTumor, levels = design) ## Tumor minus NonTumor.

fit = lmFit(expr, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit,trend = T, robust =T) ## ebayes: Empirical Bayes Statistics for Differential Expression.

DEA <- topTable(fit, adjust.method = "BH", number = Inf, sort.by = "p")	#FDR adjusted p val

## Ensemble gene id to gene symbol. 

DEA$Symbol <- lookup(rownames(DEA),Annotation23[,-2])
write.table(DEA, "limma_DEA.txt", quote = F)

##  Matched Differentially expressed long non coding RNA. 
lncRNA_DEA <- DEA[rownames(DEA) %in% lncRNA_v23$ENSG_id,]
lncRNA_DEA$ENSG_id <- rownames(lncRNA_DEA)
table(duplicated(lncRNA_DEA$Symbol))
lncRNA_DEA[duplicated(lncRNA_DEA$Symbol),]

write.table(lncRNA_DEA, "TCGA_GTEX_lncRNA_DEA.txt", quote = F)

# 9.CLINICAL FLLOWUP INFORMATIONS IN THE TCGA PROJECT
clinical = read.csv("tcga_coad_clinical_followup_data.csv",header = T, stringsAsFactors = F)

## Variable type tansformation
clinical = clinical %>% dplyr::select( - lymphnodes.positive.by.HE) %>% mutate(
  stage=factor(stage, levels = c("i", "ii", "iii", "iv"), labels = c("i", "ii", "iii", "iv")),
  vital.status = factor(vital.status, levels = c("alive", "dead", "P-value"),
  labels = c("alive", "dead","P-value")),
  venous.invasion =factor(venous.invasion,levels = c("no", "yes"), labels = c("no", "yes")),
  perineural.invasion=factor(perineural.invasion,levels = c("no", "yes"), labels = c("no", "yes")),
  person.neoplasm.status =factor(person.neoplasm.status, levels = c("tumor free", "with tumor"), labels = c("tumor free", "with tumor")),
  lymphatic.invasion = factor(lymphatic.invasion, levels = c("no", "yes"), labels = c("no", "yes")),
  gender= factor(gender, levels = c("male", "female"), labels = c("male", "female")) )
##
str(clinical)
summary(clinical)

# 10.CREATE A TABLE FOR SUMMARY INFORMATIONS AND STATISTICAL TEST .

## Reference:url="https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html"
table1::label(clinical$gender) <- "Gender"
table1::label(clinical$age) <- "Age"
table1::label(clinical$vital.status) <- "Status"
table1::label(clinical$OS) <- "Overall surrivall"
table1::label(clinical$DFS) <- "Disease free survival"
table1::label(clinical$stage) = "Stages"
table1::units(clinical$age) <- "years"
table1::units(clinical$OS) <- "days"
table1::units(clinical$DFS) <- "days"
table1::units(clinical$preoperative.CEA) <- "ng/ml"

##Defind a fuction  for statistical test

source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/rndr.R")
##CLINICAL FLLOWUP INFORMATIONS 
table.1 = table1(~OS + DFS + age + gender + preoperative.CEA + stage +
	 person.neoplasm.status + lymphatic.invasion + perineural.invasion +
	 venous.invasion | vital.status, topclass = "Rtable1-zebra",
	 footnote = "Note: Categorical variable using Chi-Square tests and T-test for continous variable.", droplevels = F, render = rndr,
	 render.strat=rndr.strat, overall=F,data=clinical) #+ location 

# 11.Basic survival analysis with clinical informations

## Different stages,genders and age group (above60,low60)
## 
clinicdata <- clinical%>% dplyr::select(id,stage,OS,status,gender,age)
clinicdata$Age.group=ifelse(clinicdata$age>=60,"above60","low60")
clinicdata=clinicdata%>%mutate(Age.group=factor(Age.group))

#
fit <- survfit(Surv(OS/365, status) ~ stage, data = clinicdata)
survdiff(formula = Surv(OS/356, status) ~ stage, data = clinicdata)
table(clinicdata$stage)

ggsurvplot(fit, data = clinicdata, risk.table = F, pval = T, ggtheme = theme_survminer(),
		   title = "Overall survival probability of stage group in TCGA-COAD",
		   xlab = "Time(years)", legend = "top", surv.median.line = "hv")



######################################### PART2 SURVIVAL ANALYSIS. #####################################################################
#NOTE: Unfortunately, the UCSC xena	toil recomputed data don't contain all the clinical follow up information.So we use the TCGA-COAD orignal data.
# gene exprssion data (HTSEQ FPKM) for surviaval analysis. All genes are annotated with gencode-v22 RNAs annotation gtf file. 
# The gene expression data downloaded from the https://portal.gdc.cancer.gov/ databaes.
#(TCGAbiolink,GDCRNATools,can be used for downloading the data)
########################################################################################################################################

# 1.Extract lncrna exprssion data from coad htseq fpkm dataset.
Exprs <- data.table::fread("COAD_HTSeq_FPKM.txt", sep = "\t", header = T, data.table = F,stringsAsFactors=F)
colnames(Exprs)<-gsub("-",".",colnames(Exprs))
##
group.2 = ifelse(grepl(".11$", colnames(Exprs)), "NonTumor", "Tumor")
table(group.2) ## 41 NonTumor and 479 Tumor with technical replications  
 
Exprs <- Exprs[, - grep("\\.11", colnames(Exprs))] ##Filter out nontumor samples for survival analysis.
##
Annotation22 <- read.table("gencode_v22_RNAs_annotation.txt", strip.white = T, 
                            blank.lines.skip = T, header = F, stringsAsFactors = F, sep = ";")[, c(1:3)]
colnames(Annotation22) = c("ENSG_id", "Type", "Symbol")

table(Exprs$ID %in% Annotation22$ENSG_id)## 60483 genes mapped to gencode_v22_RNAs_annotation
##
rownames(Exprs) <- Exprs$ID
Exprs <- Exprs[, -1]

###Some of the tumor have technical replicate samples (22 samples),and use the average value to represent.
###Note: R does not allow same rownames,so we transpose the expression matrix and collapsed it using avearrays()function. 

table(duplicated(substr(colnames(Exprs), 1, 12)))
colnames(Exprs) <- substr(colnames(Exprs), 1, 12)
Exprs <- as.data.frame(limma::avearrays(Exprs))
##
table(duplicated(colnames(Exprs)))## 456 unique sample left.


lncRNA_Exprs <- Exprs[rownames(Exprs) %in% lncRNA_v22$ENSG_id,]
## 15900 lncRNA been Extracted using lncRNA_v22

lncRNA_Exprs <- as.data.frame(log2(as.matrix(lncRNA_Exprs) + 0.001)) #Normaliztion
lncRNA_Exprs <- as.data.frame(t(lncRNA_Exprs))

# 2.Filter lncRNAs that have  nearly zero (<0.1) expression in 60% samples.
filter2 <- function(x) { length(x[x <0.1]) / length(x) > 0.6}

table(apply(lncRNA_Exprs, 2, filter2))

lncRNA_Exprs <- lncRNA_Exprs[, apply(lncRNA_Exprs, 2, filter2) == F]

lncRNA_Exprs <- data.frame(id = rownames(lncRNA_Exprs), lncRNA_Exprs)

write.table(lncRNA_Exprs, "colon_lncRNA_Exprs.txt", quote = F,sep=",",row.names=F) 
##
##lncRNA_Exprs <- data.table::fread("colon_lncRNA_Exprs.txt", sep = ",", header = T, data.table = F)

##
#3.Merge the clinical and lncrna expression matrix
clinicdata <- clinical%>% dplyr::select(id,stage,OS,status)
table(is.na(clinicdata$stage))
clinicdata <- na.omit(clinicdata) ## with complete stage informations
#
Clinic_lncRNA_Exprs <- merge(clinicdata, lncRNA_Exprs,by="id")
#
write.table(Clinic_lncRNA_Exprs, "Clinic_lncRNA_Exprs.txt", quote = F)

## checking the not mereged patients
clinicdata[clinicdata$id %in% Clinic_lncRNA_Exprs$id == F,]
## Note: Ten patients missing clincal stage information and three patients don't match the epxression data,finally 428 patients left.

##4.Univariate kaplan meier  and univariate cox proportional hazards  survival analysis.
## kaplan meier non-parametric analysis for each lncrna by gene expression group, high and low group are separated by median gene expression.

kaplan_dat <-Clinic_lncRNA_Exprs[, !colnames(Clinic_lncRNA_Exprs)%in%c("id","stage","status","OS")] 
### Tansformming the exprssion data star with the 5th colunm
OS <- Clinic_lncRNA_Exprs$OS
Status <- as.numeric(as.vector(Clinic_lncRNA_Exprs$status))
#

survplot <- function(x) { 
	OS <- Clinic_lncRNA_Exprs$OS
	Status <- as.numeric(as.vector(Clinic_lncRNA_Exprs$status))
	group = ifelse(x > median(x), 'high', 'low')
	dat <- data.frame(group = group, time = OS, status = Status)
        fit <- survfit(Surv(OS/365, status) ~ group, data = dat)
        ggsurvplot(fit, data=dat,risk.table = F, pval = T, ggtheme = theme_survminer(),
			 title = paste("Overall survival probability of   in TCGA-COAD"),
			 xlab = "Time(years)", legend = "top", surv.median.line = "hv")
	}
survplot(kaplan_dat$ENSG00000258122.1) #	choose a random gene
## Batch logrank survival analysis.
batch_log_rank_os <- apply(kaplan_dat, 2, function(values) {
	group = ifelse(values > median(values), 'high', 'low') 
	dat <- data.frame(group = group, time = OS, status = Status, stringsAsFactors = F)
	fit = surv_fit(Surv(time, status)~ group, data = dat)
	p.val = surv_pvalue(fit)})

##combine the result list to dataframe
log_rank_os <- do.call(rbind, batch_log_rank_os)
log_rank_os$ENSG_id <- rownames(log_rank_os)
log_rank_os$Symbol = lookup(log_rank_os$ENSG_id, lncRNA_v22[, -c(3,4)])
##
log_rank_os <- dplyr::select(log_rank_os, ENSG_id, pval, Symbol) %>% mutate(logrank.pvalue = pval) %>% dplyr::select(logrank.pvalue, ENSG_id,Symbol)
##
log_rank_os <-na.omit(log_rank_os)
table(duplicated(log_rank_os$Symbol))
log_rank_os[duplicated(log_rank_os$Symbol),]
##
write.table(log_rank_os, "log_rank_os.txt", quote = F)
#sig_log_rank_os <- log_rank_os %>% dplyr::filter(logrank.pvalue < 0.05)

## Apply the univariate coxph function to multiple covariates at once
#Since the other clinical feautre contains too much missing value, only stage information is included. 
covariates <- colnames(Clinic_lncRNA_Exprs)[!colnames(Clinic_lncRNA_Exprs) %in% c("id", "status", "OS")]
coxdata = Clinic_lncRNA_Exprs %>% dplyr::select(-id, - status)
#
cox_os_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS, Status)~', x))) # 
##
cox_os_models <- lapply(cox_os_formulas, function(x) { coxph(x, data =coxdata ) }) ## 

#Extract data (HR.pvalue,wald.test,beta,CI,)
#Note: The stage is a four levle categrical variable, doesn't comparable to the other result list when converting to dataframe. 
# A tricky way is first extract it beforehand.
cox_stage=cox_os_models[[1]] ## stage 
summary(cox_stage)
#
library(broom)
tidy_cox_stage=tidy(cox_stage,exponentiate = T, conf.int = T,conf.level = 0.95)
tidy_cox_stage$beta=cox_stage[["coefficients"]]
tidy_cox_stage = tidy_cox_stage %>% dplyr::select(-std.error)%>%mutate(
CI=paste0("(",round(conf.low,3),"-",round(conf.high,3),")"))%>%transmute(term=term,beta=beta,HR = estimate, CI=CI,wald.test=statistic,HR.pvalue=p.value)

cox_os_models[[1]]	=NULL
##
coxph_os <- lapply(cox_os_models,function(x) {
						x = summary(x)
						HR.pvalue = signif(x$wald["pvalue"], digits = 3)
						wald.test = signif(x$wald["test"], digits = 3)
						beta = signif(x$coef[1], digits = 2);
						HR = signif(x$coef[2], digits = 2);
						HR.confint.lower = signif(x$conf.int[, "lower .95"], 3)
						HR.confint.upper = signif(x$conf.int[, "upper .95"], 3)
						HR.confint = paste0(" (",
		                                HR.confint.lower, "-", HR.confint.upper, ")")
						res = c(beta, HR, HR.confint, wald.test, HR.pvalue)
						names(res) = c("beta", "HR", "CI", "wald.test","HR.pvalue")
						return(res)})
coxph_os_res <- as.data.frame(do.call(rbind, coxph_os),stringsAsFactors = F)
coxph_os_res = data.frame(ENSG_id = rownames(coxph_os_res), coxph_os_res)
#
coxph_os_res <- coxph_os_res %>% mutate(beta = as.numeric(beta), HR = as.numeric(HR), wald.test = as.numeric(wald.test),  HR.pvalue = as.numeric(HR.pvalue) )

coxph_os_res$Symbol <- lookup(coxph_os_res$ENSG_id, lncRNA_v22)
coxph_os_res <- na.omit(coxph_os_res)

#coxph_os_res<-rbind(tidy_cox_stage,coxph_os_res)
write.csv(coxph_os_res, "coxph_os.csv",row.names=F)

#
lncRNA_DEA <- lncRNA_DEA %>% dplyr::select(logFC, adj.P.Val, Symbol, ENSG_id)
#
DEA_OS_analysis <- merge(merge(lncRNA_DEA, log_rank_os, by = "ENSG_id"), coxph_os_res, by = "ENSG_id")
DEA_OS_analysis <- DEA_OS_analysis %>% dplyr::select(-Symbol.x, -Symbol.y)
#
write.csv(DEA_OS_analysis, "DEA_OS_analysis.csv")
#

## Set criterion to extract significant lncRNAs
# Note: If the 	criterion is too strigent there might be less candidate been left.
#	significant survival analysis results(abs(logFC)>=1,HR>1,

#sig_DEA_OS <- dplyr::filter(DEA_OS_analysis, logrank.pvalue <= 0.05 & HR.pvalue <= 0.05 & adj.P.Val <= 0.05&abs(logFC) >= 0.5)

#dim(sig_DEA_OS)

# Inference with Regularized Cox and Generalized Linear Models, and select the nozero coffecients variables(features).
# With 10 fold cross validation. The selected features are at the optimal Log lamda value.

survdata = Clinic_lncRNA_Exprs[, colnames(Clinic_lncRNA_Exprs) %in% c("id", "OS", "stage", "status", DEA_OS_analysis$ENSG_id )]	 
survdata <- survdata %>% dplyr::mutate(status = as.numeric(as.vector(status)))

# glmfit <- glmnet(y = Surv(survdata$OS, survdata$status), x = as.matrix(survdata[, 5:ncol(survdata)]),
	                 #family = "cox", nlambda = 50, alpha = 0.5, standardize = T)
# plot(glmfit)

#####	Parameter tuning for the elastic net penalized Cox model( Ridge Plus Lasso with 10-fold cross validation.)

for (i in 0:10) {
	  set.seed(12345)
	  require(doParallel)
	  cl=makeCluster(3)
	  registerDoParallel(cl)
	  assign(paste("cvfit", i, sep = ""),
	  cv.glmnet(y = Surv(survdata$OS, survdata$status), x = as.matrix(survdata[, 5:ncol(survdata)]),
	  family = "cox", standardize = T, alpha = i/10, nlambda = 50, nfolds = 10, parallel = T))
}
stopCluster(cl)
# Plot the Solution Paths

par(mfrow = c(1, 4) ,mar = c(5, 2, 5, 2), las = 1, family = "Arial")
plot(cvfit10, main = "LASSO")
plot(cvfit0, main = "Ridge")
plot(cvfit5, main = paste("Elastic Net","alpha=0.5"))
plot(cvfit1, main = paste("Elastic Net", "alpha = 0.1"))
###

###Searching for best tunning Elastic Net parameters  (alpha and lamda)  with epsgo	function.

bounds <- t(data.frame(alpha = c(0, 1)))
colnames(bounds) <- c("lower", "upper")
set.seed(12345)
foldid <- balancedFolds(class.column.factor = cbind(time = survdata$OS, status = survdata$status)[, 2], cross.outer = 10)
table(foldid, cbind(time = survdata$OS, status = survdata$status)[, 2])
#(Efficient Parameter Selection via Global Optimization)
fit <- epsgo(Q.func = "tune.glmnet.interval", bounds = bounds, parms.coding = "none", seed = 1234,
			 fminlower = -100, x = as.matrix(survdata[, 5:ncol(survdata)]), y = cbind(time = survdata$OS, status = survdata$status), family = "cox", foldid = foldid, type.min = "lambda.min",
			 type.measure = "deviance", verbose = T)

sum_para <- summary(fit, verbose = T) #summary_parameter
plot(sum_para)

deviance = as.matrix(as.data.frame.matrix(xtabs(deviance ~ alpha + lambda, data = as.data.frame(sum_para$info))))
library(plotly)
deviance_plot = plotly::plot_ly(x = ~sum_para$info$alpha, y = ~sum_para$info$lambda, z = ~deviance) %>% add_surface(color = c("#FC4E07", "#E7B800", "#00AFBB"))
htmlwidgets::saveWidget(deviance_plot, "deviance_plot.html")
## 
require(doParallel)
cl = makeCluster(3)
registerDoParallel(cl)
set.seed(1234)	###
cv.glmfit <- cv.glmnet(y = Surv(survdata$OS, survdata$status), x = as.matrix(survdata[, 5:ncol(survdata)]),
			 family = "cox", parallel = T, standardize = T, alpha = sum_para$opt.alpha, nlambda = 50)
#
stopCluster(cl)
#
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), las = 1, family = "Arial")

plot(cv.glmfit, main="Elastic Net Regularized CoxPH Regression",xlab=expression(log(lambda)))
abline(h = sum_para$opt.error, lty = 3)

res.cv <-cv.glmfit$glmnet.fit

cof <- coef(res.cv, s = cv.glmfit$lambda.min)
cof = as.data.frame(as.matrix(cof))##transform dgCMatrix to dataframe
cof$ENSG_id = row.names(cof)
cof$Symbol = lookup(cof$ENSG_id, lncRNA_v22)

colnames(cof) = c("beta", "ENSG_id","Symbol")
cv.glmfit$lambda.min
cv.glmfit$lambda.1se
(nonozero.coef=dplyr::filter(cof, beta != 0) )

nonozero.coef <- nonozero.coef %>% arrange(desc(beta))

c060::Plot.coef.glmnet(cvfit = cv.glmfit, betas = nonozero.coef$ENSG_id)

(res.cv <- cv.glmfit$glmnet.fit)

## Another way to plot the results
cof <- coef(res.cv, s = cv.glmfit$lambda.min)
cofn <- cof[which(cof != 0)]
names.cof <- rownames(cof)
cofn <- cof[which(cof != 0)]
names(cofn) <- names.cof[which(cof != 0)]
bet <- res.cv$beta[match(names(cofn), rownames(res.cv$beta)),]
glmnet:::plotCoef(bet, lambda = res.cv$lambda, df = res.cv$df, dev = res.cv$dev.ratio,
				 xvar = "lambda", add = F, col = "red", label = T,main="Coefficients shrankage paths", xlab = expression(log(lambda)))
abline(v = log(cv.glmfit$lambda.min), lty = 3)
abline(v = log(cv.glmfit$lambda.1se), lty = 3)
####
glmnet:::plotCoef(bet, lambda = res.cv$lambda, df = res.cv$df, dev = res.cv$dev.ratio,
				 xvar = "dev", add = F, col = "red", label = T, main = "Coefficients shrankage paths of Elastic Net Regularized Cox PH Regression model", xlab = "Fraction Deviance Explained")
abline(v = 0.05070, lty = 3);abline(v = 0, lty = 3)

#A summary table of the features of the selected lncRNAs.
ggtextdata = DEA_OS_analysis[DEA_OS_analysis$ENSG_id %in% nonozero.coef$ENSG_id,]
ggtextdata$Type = lookup(ggtextdata$ENSG_id, lncRNA_v22[, - c(2, 3)])
ggtextdata$Position = lookup(ggtextdata$ENSG_id, lncRNA_v22[, c(1, 3)])

ggtextdata <- ggtextdata %>% dplyr::select(Symbol, ENSG_id, logFC, adj.P.Val,
										   HR, CI, HR.pvalue, logrank.pvalue, Type, Position) %>% mutate(HR.pvalue = round(HR.pvalue, 4),logrank.pvalue = round(logrank.pvalue, 4), adj.P.Val = round(adj.P.Val, 4)) %>% arrange(desc(HR)) #Position, Type

table.2 = ggtexttable(ggtextdata, theme = ttheme("mBlue", base_size = 8), rows = NULL)

#
ggdotchart1 = ggdotchart(ggtextdata, x = "Symbol", y = "logFC", dot.size = "HR", color = "HR.pvalue",
			add = "segments", sorting = "descending",
			ggtheme = theme_pubr()) +
			geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
			scale_color_gradientn(colors = c("#FC4E07", "#E7B800", "#00AFBB")) +
			theme(legend.position = "top") +
			labs(title = "Features of candidate lncRNAs") +
			theme(plot.title = element_text(size = 12, hjust = 0.5, lineheight = 0.2),
			axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 0.1)) +
			geom_rug(alpha = 0.8, size = 0.5, color = "grey")


##5.Data are randomly partitioned into train (70%) and test (30%) according to vital status.
survdata$group <- ifelse(survdata$status == "0", "alive", "dead")
survdata <- survdata %>% mutate(group = factor(group, levels = c("alive", "dead"), labels = c("alive", "dead")))
table(survdata$group)
#
set.seed(12345) ##Set random seeds for reproducible results.
trainIndex <- createDataPartition(survdata$group, p = .7, list = FALSE, times = 1)

Train <- survdata[trainIndex,]
Test <- survdata[-trainIndex,]

table(Train$group)
table(Test$group)

Train <- Train %>% dplyr::select(-group)
Test <- Test %>% dplyr::select(-group)

table(Train$status)
table(Test$status)
Train <- Train[, colnames(Train) %in% c("id", "OS", "stage", "status", nonozero.coef$ENSG_id)]
Test <- Test[, colnames(Test) %in% c("id", "OS", "stage", "status", nonozero.coef$ENSG_id)]

##5.Check the information between train and validating sets.create a html style summary table.

clinical = clinical[clinical$id %in% Clinic_lncRNA_Exprs$id,]
clinical$group = ifelse(clinical$id %in% Train$id, "Training", "Validating")
clinical = clinical %>% dplyr::mutate(group = factor(group, levels = c("Training", "Validating", "P-value"), labels = c("Training", "Validating", "P-value")))
clinical = clinical %>% dplyr::mutate(status = factor(status))
#
table1::label(clinical$group) = "Group"
#
source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/rndr2.R")
table.3 = table1(~OS  + age + status + gender + preoperative.CEA + stage  + lymphatic.invasion + perineural.invasion + venous.invasion | group,
				 topclass = "Rtable1-zebra", droplevels = F, render = rndr2, render.strat = rndr.strat, overall = F, data = clinical,
	       footnote = "Note:Chi-sq test for categorical variables,continuity correction is only used in the 2-by-2 case, and t-test for continous variables.")
###


##6.Multivariable cox proportional hazards model  analysis and feature selection
##Note: _t for Train data,_v for Validating data.
cox_data_train <- Train %>%dplyr:: select(-id) %>% mutate(status = as.numeric(as.vector(status)))
#
cox_data_test <- Test %>% dplyr::select(-id) %>% mutate(status = as.numeric(as.vector(status)))
#

cox <- coxph(Surv(OS, status) ~ ., data = cox_data_train)
summary(cox)

stepcox = step(cox, trace = 1, scale = 0,
		 direction = c("both", "backward", "forward"),
		 steps = 1000, k = 2) ##trace =0 supressing the 
summary(stepcox)
#coxph.detail(stepcox)												  
coef(stepcox)
##
stepcox_res <- as.data.frame(summary(stepcox)[["coefficients"]])

lookup(rownames(stepcox_res), lncRNA_v23)
stepcox_res$Symbol <- lookup(rownames(stepcox_res), lncRNA_v22[-c(3,4)])
stepcox_res$Items = rownames(stepcox_res)


##Summary table														
table.4= ggtexttable(stepcox_res, theme = ttheme("mBlue", base_size = 8), rows = NULL) #

cox_vip_list = names(stepcox$means) ##8 lncRNAs selected by step( backward?Cforward ) coxph analysis
#####

##Forestplot
ggforest(stepcox, main = "Impotance of candidates in Overall survival of TCGA-COAD")
ggtexttable(DEA_OS_analysis[DEA_OS_analysis$ENSG_id %in% cox_vip_list,], theme = ttheme("mBlue", base_size =10), rows = NULL)
#####
print(DEA_OS_analysis[DEA_OS_analysis$ENSG_id %in% cox_vip_list,]) ###

##7.Multivariable cox proportional hazards model risk prediction on training and validating data sets

## CoxPH risk prediction model evaluated by ROC curve and AUC
library(survivalROC)
source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/PlotsurvROC.R")
par(mfrow = c(1, 1))
PlotsurvROC(1,Train,Test)
PlotsurvROC(3,Train,Test)
PlotsurvROC(5, Train, Test)
##
source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/CALAUC.R")
AUC = data.frame(time = rep(seq(1, 15, 0.05), 2),
			group = rep(c("Training", "Testing"), each = 0.5*length(rep(seq(1, 15, 0.05), 2))),
	    AUC=c(sapply(seq(1, 15, 0.05), CALAUC_train),sapply(seq(1,15,0.05),CALAUC_test)))

(AUC %>% dplyr::filter(group == "Training") %>% arrange(desc(AUC)))
(AUC %>% dplyr::filter(group == "Testing") %>% arrange(desc(AUC)))
ggplot(AUC, aes(x = time, y = AUC, colour = group)) +
	geom_line(size = 1.2) +scale_color_manual(values = c("#FC4E07","#0073C2FF" ))	+
	geom_rug(sides = "l") + labs(title = "Multivariable CoxPH model Risk prediction Performance", y = "Area Under Curve (AUC)", x = "Observation Time (years)") +
	geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey", size = 0.5) +
	geom_vline(xintercept = 3, linetype = "dashed", color = "grey", size = 0.5) +
	geom_vline(xintercept = 5, linetype = "dashed", color = "grey", size = 0.5)	+
    theme(plot.title = element_text(size = 20, hjust = 0.4, lineheight = 0.2), axis.text.x = element_text(size = 6)) +
	theme_pubr() 


##8.Visualize with survminer patients in risk groups
source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/GGsurvplot.R")
GGsurvplot(cox_vip_list, Train, stepcox)
GGsurvplot(cox_vip_list, Test, stepcox)


##Risk score distribution among patients
risk = function(dataset) {
	newdata = dataset[, colnames(dataset) %in% c(cox_vip_list, "status", "stage", "OS", "id")]
	riskScore <- predict(stepcox, type = "risk", newdata = newdata, reference = "strata")
	risk <- as.factor(ifelse(riskScore > median(riskScore), "high", "low"))
	risk_df <- cbind.data.frame(newdata[, c("id", "OS", "status")], riskScore = riskScore, risk = risk)
	risk_df <- risk_df %>% dplyr::mutate(OS = as.numeric(OS), riskScore = as.numeric(riskScore)) %>% arrange(desc(riskScore)) %>% mutate(log10.riskScore = log10(riskScore), zscore = (riskScore - mean(riskScore))/sd(riskScore))
        risk_df$survival.status = ifelse(risk_df$status == 0, "Censored", "Dead")
        risk_df$survival.status <- as.factor(risk_df$survival.status)
	return(risk_df)
	}


risk_df <-rbind( data.frame(risk(Train), group = rep("Training", nrow(risk(Train)))),
          data.frame(risk(Test), group = rep("Testing", nrow(risk(Test)))))

ggboxplot(risk_df, x = "group", y = "log10.riskScore", color = "survival.status", palette = c("#00AFBB", "#FC4E07"), add = "jitter") +
	stat_compare_means(aes(group = survival.status), method = "t.test", label = "p.signif" + stat_compare_means(label.y = 14))

#
#plot(risk_df$zscore,main = "Risk score distribution of TCGA-COAD patients")
#show_point_shapes()
par(mfrow = c(1, 1), mar = c(1, 1,1,1))
ggdotchart2 = ggdotchart(risk_df, x = "id", y = "log10.riskScore", shape = "group",
	 color = "survival.status",
	 title = "Risk score distribution of TCGA-COAD patients",
	 repel = T, ylab = F, xlab = expression(paste(Log[10], " ( RiskScore)")), size = 1, palette = c('#0073C2FF', '#FC4E07'),
	 ggtheme = theme_pubr()) + coord_flip() +
	 theme(axis.text.y = element_text(size = 5,, vjust = 1), plot.title = element_text(size = 15, hjust = 0.5, lineheight = 0.2)) +
	 geom_hline(yintercept = median(risk_df$log10.riskScore), linetype = "dashed", color = "darkgrey") +
 	 geom_rug(aes(color = survival.status), data = risk_df, size = 0.1, alpha = 1, position = "jitter") +
	 scale_shape_manual(values = c(18, 20)) +
	 theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) 
         #geom_jitter(position = position_jitter(height = .05)) #,axis.title.y = element_blank())

library(dplyr)

summary_dat <- Clinic_lncRNA_Exprs[, colnames(Clinic_lncRNA_Exprs) %in% c("id", cox_vip_list)]

row.names(summary_dat) = summary_dat$id

meregedata = merge(risk_df, summary_dat, by = "id")
meregedata = meregedata %>% arrange(desc(log10.riskScore))
meregedata =meregedata %>% dplyr::select(-id,-status) %>% mutate(
group = factor(group, levels = c("Testing", "Training"), labels = c("Testing", "Training")),
survival.status = factor(survival.status, levels = c("Censored", "Dead", "P-value"),
labels = c("Censored", "Dead", "P-value")))
#
str(meregedata)
colnames(meregedata)[8:15] =lookup(colnames(meregedata)[8:15], Annotation23[, -2])

#CREATE A TABLE FOR SUMMARY INFORMATIONS AND STATISTICAL TEST .
library(table1)
table1::label(meregedata$OS) <- "Overall surrivall"
table1::label(meregedata$survival.status) <- "Status"
table1::label(meregedata$group) <- "Group"
table1::units(meregedata$OS) <- "days"

source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/rndr3.R")
table1(~OS + riskScore + risk + log10.riskScore + group + `RP11-440D17.3` + `EIF3J-AS1` + `TNRC6C-AS1` + MIR210HG + `AC113189.5` + LINC00261 + `CTB-25B13.12` + AC114730.3 | survival.status, group, topclass = "Rtable1-zebra", data = meregedata, droplevels = F, render = rndr3, overall = F, render.strat = rndr.strat) #


##9.Randomforestsrc survival analysis (RF-SRC)).
rfdata_train <- Train %>% dplyr::select(-id) %>% mutate(OS = OS/365)
rfdata_train <- rfdata_train[, colnames(rfdata_train) %in% c("stage", "OS", "status", cox_vip_list)]

rfdata_test <-	Test%>%dplyr:: select( - id) %>% dplyr::mutate(OS = OS/365) ###convert days to years
rfdata_test <- rfdata_test[, colnames(rfdata_test) %in% c("stage", "OS", "status", cox_vip_list)]

set.seed(12345)
options(rf.cores = 4)
## Tunning hyperparamters to obtain miminal OOB error
(rf_tune= tune(Surv(OS, status) ~., data = rfdata_train, mtryStart = floor(sqrt(ncol(rfdata_train))),
          nodesizeTry = c(1:9, seq(10, 100, by = 5)), ntreeTry = 500,
          stepFactor = 1.25, improve = 1e-3, strikeout = 3, maxIter = 50,
          trace = T, doBest = T))
#
rftun_res=as.data.frame(rf_tune[["results"]])

mtry = rftun_res$mtry
nodesize= rftun_res$nodesize
OOBerror=as.matrix(as.data.frame.matrix(xtabs(err~nodesize+mtry,data=rftun_res))) 
library(plotly)
OOBerror_plot = plotly::plot_ly(x = ~mtry, y = ~nodesize, z = ~OOBerror) %>% add_surface(color = c("#FC4E07", "#E7B800", "#00AFBB")) %>%
	layout(autosize = T, scene = list(xaxis = list(range = c(0, 5)),
	yaxis = list(range = c(0, 10)),
	zaxis = list(range = c(0, 0.3))))
#	
htmlwidgets::saveWidget(OOBerror_plot, "OOBerror_plot.html")

rf_tune$optimal		   
#
source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/PlotRFtune.R")
PlotRFtune(rf_tune)
##
options(rf.cores = 4)
set.seed(1234)#
(rfsrc <- rfsrc(Surv(OS, status) ~ ., data = rfdata_train, nsplit = 3,
	   mtry = rf_tune$optimal["mtry"],
	   nodesize=rf_tune$optimal["nodesize"],
	   bootstrap = "by.root",
	   samptype = "swor",
           na.action = "na.impute", nimpute = 1,
	   tree.err = T, split.depth = "by.tree",
           importance = "permute", ntree = 500, do.trace = T, statistics=T))

plot(rfsrc, main = "RandomForest OOB error and Vimp")

###
##prediction on valdation data
(rf.pred <- predict(rfsrc,newdata =rfdata_test,na.action = "na.impute"))

##10.Survival information

gg_data_train = as.data.frame(gg_survival(interval = "OS", censor = "status", data = rfdata_train,type = c("kaplan") ))
gg_data_train$Group=rep("training",nrow(gg_data_train ))
gg_data_test=as.data.frame(gg_survival(interval = "OS", censor = "status", data = rfdata_test,type = c("kaplan") ))
gg_data_test$Group=rep("testing",nrow(gg_data_test))
ggdata=rbind(gg_data_test,gg_data_train)
#SURVIVAL PROBABILITY OF PATIENTS IN TCGA-COAD
ggplot(ggdata, aes(x = time, y = surv, colour = Group)) +
	geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.05) +
	geom_step(size = 1) +
	geom_rug() + labs(title = "Survival Probability of patients in TCGA-COAD", y = "Survival Probability", x = "Observation Time (years)") +
	geom_vline(xintercept = 5, linetype = "dashed", color = "grey", size = 0.5) +
	geom_hline(yintercept = 0.6, linetype = "dashed", color = "grey", size = 0.5)+
        theme(plot.title = element_text(size = 20, hjust = 1, lineheight = 0.2), axis.text.x = element_text(size = 6)) +theme_pubr() 
#CUMULATIVE HAZARD OF PATIENTS IN TCGA-COAD	 geom_step ,geom_line
ggplot(ggdata, aes(x = time, y = cum_haz, colour = Group)) +
	geom_step(size = 1) +
	geom_rug() + labs(title = "Cumulative Hazard of patients in TCGA-COAD", y = "Cumulative Hazard", x = "Observation Time (years)") +
	geom_vline(xintercept = 5, linetype = "dashed", color = "grey", size = 0.5) +
	geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey", size = 0.5)+
        theme(plot.title = element_text(size = 20, hjust = 1, lineheight = 0.2), axis.text.x = element_text(size = 6)) +theme_pubr() 


#RANDOM FOREST PREDICTED SURVIVAL OF PATIENTS IN TRAINING DATA
#Random forest predicted survival. Blue lines correspond to censored observations,
#red lines correspond to patients who experienced the event(death) .

plot(gg_rfsrc(rfsrc)) + theme(legend.position = c(0.1, 0.2)) +
	labs(title = "Random forest predicted survival of patients in training data", y = "Survival Probability", x = "Observation Time (years)") +
	theme(plot.title = element_text(size = 20, hjust = 0, lineheight = 0.2), axis.text.x = element_text(size = 6)) + geom_vline(xintercept = c(1, 3), linetype = "dashed") +
	coord_cartesian(x = c(0, 4)) + theme_pubr() + geom_rug() +
 	annotate(geom = "text", x = 1, y = 0.3,
	label = paste("OOB error",round(mean(na.omit(rfsrc$err.rate))*100,3),"%"),
	color = "black") +scale_color_manual(values = c("#0073C2FF", "#FC4E07"))

#RANDOM FOREST PREDICTED SURVIVAL OF PATIENTS IN VALIDATING DATA
plot(gg_rfsrc(rf.pred)) + theme(legend.position = c(0.1, 0.2)) +
	labs(title = "Random forest predicted survival of patients in testing data", y = "Survival Probability", x = "Observation Time (years)") +
	theme(plot.title = element_text(size = 20, hjust = 0, lineheight = 0.2), axis.text.x = element_text(size = 6)) +
	geom_vline(xintercept = c(1, 3), linetype = "dashed") +
	coord_cartesian(x = c(0, 4)) + theme_pubr() + theme_pubr() + geom_rug() +
        annotate(geom = "text", x = 1, y = 0.3,
	label = paste("OOB error", round(mean(na.omit(rf.pred$err.rate))*100, 3), "%"),
	color = "black") +scale_color_manual(values = c("#0073C2FF", "#FC4E07"))
#Variable importance (VIMP) was originally defined in CART using a measure involving surrogate
#variablesThe most popular VIMP method uses a prediction error approach involving \noising-up" each variable in turn. VIMP for a variable values xv is the di erence between prediction error when xv is randomly permuted, compared to
#prediction error under the observed values. In VIMP, prognostic risk factors are determined by testing the forest prediction under alternative data settings, ranking the most important variables according to their impact on predictive ability of the forest.

vimp = gg_vimp(rfsrc)
vimp$vars= stringr::str_replace(vimp$vars, "\\.{1}", "-" )
vimp <- vimp %>% arrange(desc(vimp))

plot(gg_vimp(rfsrc)) +
	theme(legend.position ="top") +
	labs(fill = "VIMP > 0", title = "Impotance of candidates in OS of TCGA-COAD") +	
	theme(plot.title = element_text(size = 15, hjust = 0.8, lineheight = 0.2), axis.text.y = element_text(size =10))
	
vimplist = dplyr::filter(vimp, positive == "TRUE") %>% arrange(desc(vimp)) %>% dplyr::select(vars)
vimplist = as.vector(vimplist$vars)

varsel <- var.select(rfsrc)
gg_md <- gg_minimal_depth(varsel)
print(gg_md)

##
plot(gg_md) + labs( title = "Minimal depth of candidates in OS of TCGA-COAD ") + theme(plot.title = element_text(size = 15, hjust = 0.8, lineheight = 0.2), axis.text.y = element_text(size = 10))
#
plot(gg_minimal_vimp(gg_md)) +
	theme(legend.position = "top") + labs(title = "Impotance of candidate on Overall Survival of of TCGA-COAD") + theme(plot.title = element_text(size = 15, hjust = 0.8, lineheight = 0.2), axis.text.y = element_text(size = 10))
#
gg_v <- gg_variable(rfsrc, time = c(1,3,5), time.labels = c("1 Years","3 Years","5 Years"))
#
plot(gg_v, xvar = "stage", alpha = 0.6) +
		 labs(title = "Survival Probility of stages") +
		 theme(legend.position = "top") +
		 labs(y = "Survival", x = "Tumor Stage") + scale_color_manual(values = c("#0073C2FF", "#FC4E07"))	+
		 theme(plot.title = element_text(size = 15, hjust = 0.5, lineheight = 0.2)) 

#
plot(gg_v, xvar = colnames(gg_v)[2:9],panel = T, alpha = 0.1) + theme(legend.position = "top") +
		 labs(title = "LncRNA Expression dependent Overall Survival Probability", y = "Survival Probability", x = "lncRNA Expression log2(FPKM+0.001)") +
		 scale_color_manual(values = c("#0073C2FF", "#FC4E07")) + coord_cartesian(ylim = c(0.4, 1))
 	         theme(plot.title = element_text(size = 12, hjust = 0.5, lineheight = 0.2))

#COMPUTE THE PREDICTION ERROR AND C INDEX FIT CANDIDATE COX AND RANDOMFORESTSRC MODELS AND COMPUTE THE KAPLAN - MEIER ESTIMATE. set.seed(12345)
#estimation of the expected Brier score in general survival models with right - censored event times
# prodlim Functions For Estimating Probabilities From Right Censored Data

source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/Bier.R")
#The Brier Score measures the accuracy of probabilistic predictions. with 100 bootstrap cross-validation  permution
#
Bier_train <- Bier(Train)
Bier_test <- Bier(Test)

par(mfrow=c(2,2))
plot(Bier_train, smooth = T, lwd = 1.5, legend.cex = 0.8, legend.lty = 1, type = "s", add.refline = T, xlim = c(0, 3650), ylim = c(0, 0.5), xlab = "Time(days)",
			 ylab = "Prediction error (Brier score)", col = c("black","green", "#FC4E07", "#0073C2FF")) 
lines(Bier_train$time, unlist(Bier_t$AppErr["Reference"]), "s", lwd = 1.5)
#
plot(Bier_test, smooth = T, lwd = 1.5, legend.cex = 0.8, legend.lty = 1, type = "s", add.refline = T, xlim = c(0, 3650), ylim = c(0, 0.5), xlab = "Time(days)",
			 ylab = "Prediction error (Brier score)")
lines(Bier_test$time, unlist(Bier_v$AppErr["Reference"]), "s", lwd = 1.5)

#The concordance probability(C-index) is the frequency of concordant pairs among all pairs of subjects. It can be usedto measure and compare the discriminative power of a risk prediction models. 

source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/Cindex.R")

Cindex_train <- Cindex(Train)
Cindex_test <- Cindex(Test)
#
plot(Cindex_train, smooth = T, lwd=1.5,legend.cex = 0.8, legend.lty = 1, type = "s", add.refline = T, xlim = c(0, 3650), ylim = c(0, 0.5), xlab = "Time(days)",
			 ylab = "Concordance Index (C-index)", col = c("green", "#FC4E07", "#0073C2FF"))
#lines(Cindex_train$time, unlist(Cindex_train$AppCindex["RandomForest"]), "s", lwd = 1.5,col="#0073C2FF")
#
plot(Cindex_test, smooth = T, lwd = 1.5, legend.cex = 0.8, legend.lty = 1, type = "s", add.refline = T, xlim = c(0, 3650), ylim = c(0, 0.5), xlab = "Time(days)",
			 ylab = "Concordance Index (C-index)", col = c("green", "#FC4E07", "#0073C2FF"))
#lines(Cindex_test$time, unlist(Cindex_test$AppCindex["RandomForest"]), "s", lwd = 1.5, col = "#0073C2FF")


##################################### PART 3 CLASSIFICATION MODELS (Machine learning classification)#####################################

#preprocessing data using caret preProcess fuction with BoxCox method to remove skewness and leading to normality.
#Ref:https://topepo.github.io/caret/index.html
##train data preprocessing
#########################################################################################################################################

classification_data <- expr[rownames(expr) %in% cox_vip_list,]
classification_data <- as.data.frame(t(classification_data))
colnames(classification_data) = lookup(colnames(classification_data), lncRNA_v23[, c(1, 2)])
#
classification_data<- data.frame(group=ifelse(grepl(".11$|GTEX", rownames(classification_data)), "NonTumor", "Tumor"),classification_data)
#
normalization <- preProcess(classification_data, method = "BoxCox") #"center", "scale",
classification_data <- predict(normalization, newdata = classification_data)
classification_data$group = factor(classification_data$group ,levels = c("NonTumor", "Tumor"), labels = c("NonTumor", "Tumor"))

table(classification_data$group)
#
write.csv(classification_data, "Machine_learning_data.csv")
#classification_data = read.csv("Machine_learning_data.csv", row.names = 1)

## check near Zero Variables and high correlated variables
nz <- nearZeroVar(classification_data[, -1], saveMetrics = TRUE)
descrCor <- cor(classification_data[, -1])
highCorr <- sum(abs(descrCor[upper.tri(descrCor)]) > .8)
summary(descrCor[upper.tri(descrCor)])
#
library(AppliedPredictiveModeling)
#
par(mar = c(6, 6, 10, 6) + 0.1)
## plot lncRNA expression features (distribution and boxplot)
par(mfrow = c(2, 1))
transparentTheme(trans = .9)
featurePlot(x =classification_data[,-1], center = F, scale = F, main = "Feature density distribution Plot",
			y = classification_data$group, plot = "density", scales = list(x = list(relation = "free"), y = list(relation = "free")),
			adjust = 0.5, pch = "|", layout = c(4, 2), auto.key = list(columns = 5))


featurePlot(x = classification_data[,-1], center = F, scale = F,
			y = classification_data$group,plot = "box", main = "Feature Box-Whisker Plot",
			scales = list(y = list(relation = "free"), x = list(rot = 90)),
			layout = c(4,2), auto.key = list(columns = 2))

## Data partitioned into 70% for training and 30% for testing.

set.seed(1234)
ind = createDataPartition(classification_data$group, p = 0.7, list = FALSE)
train <- classification_data[ind,]
test <- classification_data[-ind,]
## Custom Control Pararmeters with Cross-Validated (10 fold, repeated 5 times)
options(digits = 4)

custom <- trainControl(method = "repeatedcv", number = 10,
					   repeats = 5, verboseIter = T,
					   classProbs = T, savePredictions = 'final',
					   summaryFunction = twoClassSummary) #,search = 'random'


logisreg <- train(group ~ ., data = train, method = "glmStepAIC", family = "binomial",
				 tuneLength = 10,
				 trControl = custom) # method = 'glmStepAIC',glm
logisreg
summary(logisreg)

#par(mar = c(5, 5, 8, 6))
ggplot(varImp(logisreg, scale = T), main = "Importance of features in Logistic Regression for Training data")
varImp(logisreg, scale = T)

# Regularized Elastic Net Regression model.
# the elastic net is a regularized regression method that linearly combines the L1 and L2 penalties of the lasso and ridge methods.

options(digits = 4)
set.seed(12345)
#alpha (Mixing Percentage) alpha=0 ridge, alpha=1 lasso,lambda (Regularization Parameter)
elastnet <- train(group ~ ., data = train, method = "glmnet",
				         tuneLength = 10,trControl = custom, metric = "ROC",
					  tuneGrid = expand.grid(alpha = seq(0, 1, length = 10),
					  lambda = seq(0.0001, 0.1, length = 5)))

elastnet$results
elastnet$bestTune ## best tunning parameters
(best = elastnet$finalModel)

summary(elastnet)


#options(max.print= 5000, width =1000)
coef(best, s = elastnet$bestTune$lambda)
Beta=elastnet[["finalModel"]][["beta"]]
Lamda=elastnet[["finalModel"]][["lambda"]]

## plot Resampling Profile of Regularized Elastic Net Regression for Training data.

p1 <- ggplot(elastnet) + geom_vline(xintercept = elastnet$bestTune$alpha, col = "darkgrey", lty = 2) +
	theme_bw() + theme(legend.position = "top") +
	labs(title = "Resampling Profile of Elastic Net Regularized  Regression for Training data",
	     x = paste("Mixing Percentage", "(", expression(alpha), ")")) +
	theme(plot.title = element_text(size = 14, hjust = 1, lineheight = 0.2)) 
	#annotate("text", x =  elastnet$bestTune$alpha, y = 0.9985, label = "bestTune")

par(mar = c(5, 5, 8, 5))
plot(elastnet$finalModel, xvar = "lambda", label = T,
	     main = "Coefficents shrinkage of Elastic Net Regularized Logistic Regression for Training data",
	     las=1, xlab = expression(log(lambda)),xlim=c(-11,10))
       abline(v =log( elastnet$bestTune$lambda), col = "darkgrey", lty = 2)

plot(elastnet$finalModel, xvar = "dev", label = T, main = "Coefficents shrinkage of Elastic Net regression for Training data")

ggplot(varImp(elastnet, scale = T), main = "Importantance of features in Elastic Net regression for Training data")  
varImp(elastnet, scale = T)

##
sparsmatrix = as.data.frame(t(as.matrix(elastnet[["finalModel"]][["beta"]])))
#
Beta <- sparsmatrix
lambda = best$lambda

# Prepare data for plot
library(tidyr)
Beta = cbind(lambda, Beta)
Beta = Beta %>% gather(symbol, value, EIF3J.AS1:AC114730.3) %>% mutate(lambda = log(lambda)) ##wide to length
# Plot
ggplot(Beta, aes(lambda, value)) +
	geom_line(aes(colour = symbol)) +
	ggtitle("Coefficents shrinkage of Elastic Net Regularized Logistic Regression") +
	xlab(expression(log(lambda))) + theme_bw() +
	ylab("Coefficient Value") +theme(legend.position = "top")+
  geom_vline(xintercept = log(elastnet$bestTune$lambda), col = "darkgrey", lty = 2)

summary(elastnet)

##
## RandomForest model
set.seed(12345)
library(doParallel)
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
#
RF = train(group ~ ., data = train, method = "ranger",, tuneLength = 10,
	trControl = custom, metric = "ROC", tuneGrid = expand.grid(mtry = (1:ncol(train)-1),
	splitrule = "gini", min.node.size = seq(1, 10)))

stopCluster(cl)

p2 <-ggplot(RF) + geom_vline(xintercept = RF$bestTune$mtry, col = "darkgrey", lty = 2) +
			theme_bw() + theme(legend.position = "top") +
			labs(title = "Resampling Profile of randomForest model for Training data",
			x = "Randomly Selected Predictors(mtry)") +
	                theme(plot.title = element_text(size = 14, hjust = 1, lineheight = 0.2))

RF$bestTune$mtry
#ggplot(varImp(RF, scale = T))
#varImp(RF, scale = T)

#Support Vector Machines with Radial Basis Function Kernel

set.seed(12345)
library(doParallel)
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
custom <- trainControl(method = "repeatedcv", number = 10,repeats = 5, verboseIter = T,classProbs = T, savePredictions = 'final',summaryFunction = twoClassSummary, allowParallel=T)

SVM <- train(group ~ ., data = train, method = "svmRadial",
			 tuneLength = 10,
			 trControl = custom, metric = "ROC",
			 tuneGrid = expand.grid(sigma = seq(0,1, length = 10), C = seq(0.1, 2, length = 10)))

stopCluster(cl)
SVM$bestTune$C

#plot(SVM, main = "Resampling Profile of Support Vector Machine for Training data")
p3 <- ggplot(SVM) + theme_bw() + theme(legend.position = "top") + labs(title = "Resampling Profile of Support Vector Machine for Training data", x = "Optimal Selected Predictors(sigma)") +
      geom_vline(xintercept = SVM$bestTune$sigma, col = "darkgrey", lty = 2) +
      theme(plot.title = element_text(size = 14, hjust = 1, lineheight = 0.2))

ggplot(varImp(SVM, scale = T)) 
SVM$bestTune

#Neural Network
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
set.seed(12345)

NNET <- train(group ~ ., data = train, method = "nnet",tuneLength = 10,
						 trControl = custom, metric = "ROC", allowParallel = T ,
						 tuneGrid = expand.grid(size = c(1:15), decay = seq(0.1, 1, length = 5)))

NNET
plot(varImp(NNET, scale = T), main = "Importantance of features in Neural Network model for Training data")

p4 <- ggplot(NNET) + theme_bw() + theme(legend.position = "top") + labs(title = "Resampling Profile of Neural Network for Training data", x = "Hidden Units(Size)") +
	    geom_vline(xintercept = NNET$bestTune$size, col = "darkgrey", lty = 2) +
      theme(plot.title = element_text(size = 14, hjust = 1, lineheight = 0.2))

stopCluster(cl)
NNET$bestTune$size

###
## Model evaluation using Accuracy(ACC),Area under curve (AUC) 

source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/Model evaluation.R")

Modelperformance = rbind(model_evaluation(logisreg),
						             model_evaluation(RF),
						             model_evaluation(SVM),
		                             model_evaluation(elastnet),
		                             model_evaluation(NNET))

Modelperformance = data.frame(Models = rep(c("LogisReg", "RandomForest", "SVM", "ElastNetReg", "NeuralNet"), each = 2), Modelperformance)

ggtexttable(Modelperformance, theme = ttheme("mBlue", base_size = 8), rows = NULL) 

source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/multiplot.R")
multiplot(p1,  p4,p3, p2, cols = 2)
## 
source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/Plot_model_ROC.R")
Plot_model_ROC(NNET, train) 
Plot_model_ROC(NNET, test, "#E7B800", add = T)

source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/Plot multiple ROC Curve.R")
Plot_multi_ROC(train)
Plot_multi_ROC(test)

#####   Lift Curves
source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/Plot_Lift_Curve.R")

LiftPlot(train)
LiftPlot(test)
## Calibration Curves

source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/CaliPlot.R")
CaliPlot(train)
CaliPlot(test)

## Compare Models
model.list <- list(Logistic_Regression = logisreg, Random_Forest = RF, Elastnet = elastnet, SVM = SVM,Neural_Network=NNET)
res <- resamples(model.list)
summary(res)
par(mfrow = c(1, 4))
bwplot(res, layout = c(4, 1), main = "Model performance Comparision on Training data")
xyplot(res, metric = "ROC")#,"Model performance Comparision matrix on Training data")
splom(res)
#densityplot(res);
difValues <- diff(res)
difValues
SU = summary(difValues)
trellis.par.set(theme)
bwplot(difValues, layout = c(4, 1))
#####
