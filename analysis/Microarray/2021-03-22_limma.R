
# 2021-03-22   LIMMA

library(limma)


## Set cohort data as in Tao-Script:


# A label "cohort" is defined which will be used in the script below to tag any output files

# AffyDataset, 936 samples finding cohort with low/int RS:
cohort <- "s936_"
# clin info
clin <- t.header.find.lowint
clin$Xcombi_id <- rownames(clin)
# mat from AffyDataset, 936 samples finding cohort with low/int RS, cols=Xcombi_id  rows=AffyIDs:
mat <- log(t(tcdm.find.lowint))




### TK-EDIT:   SELECT SAMPLE SUBGROUP   ###

### TK-EDIT:   Generate "type" as factor age>=55:old or age<=50:young
clin$type <- NA
clin$type[clin$age>=55] <- "old"
clin$type[clin$age<=50] <- "young"
clin$type <- as.factor(clin$type)

### TK-EDIT:   EXCLUDE SAMPLES with age between 51 and 54 years  (clin$type=NA)
clin <- clin[!is.na(clin$type),]
mat <- mat[,colnames(mat) %in% rownames(clin)]



## Screen for Diff Expr Genes using limma()
#    include type and batch as covariates


type <- clin$type
batch <- factor(clin$datas_new09...TransBIG.40.)

table(type, batch)
#           batch
# type      1   2   3   4   6  12  16  17  20  21  25  33  40
# old      86  14 125  13  21  22  78  66   9  21  89  29  11
# young     4  32  20  14  16  10  11  49   7  25  30   2  61

 

# Construct model matrix including both  batch and type
X <- model.matrix(~batch+type)

# Use limma to fit gene expression in mat to model matrix
fit <- lmFit(mat,X)

# Empirical Bayes Statistics for Differential Expression
fit <- eBayes(fit)
results <- decideTests(fit)
summary(results)


# volcano plot validation-cohort-A

fn=paste0("./out/",cohort,"volcano-plot.pdf")
pdf(fn,w=8,h=8)
volcanoplot(fit, coef="typeyoung")
plotMD(fit, coef="typeyoung", status=results)
dev.off()


# Select 1000 top probesets as diff expressed genes

deg <- topTable(fit, coef="typeyoung", n=1000)

fn=paste0("./out/",cohort,"DEG-1000list.txt")
write.table(deg,file=fn,row.names=TRUE, col.names = NA, quote=FALSE, sep="\t")




######################################################################################################



####  Repeat complete DEG-Analysis with validation-2 cohort

# First save list of 1000 DEG from finding cohort:
deg1 <- deg

# Variable deg will be overwritten at the end of the subsequent script:


# A label "cohort" is defined which will be used in the script below to tag any output files

# AffyDataset, 669 samples VALIDATION cohort with low/int RS:
cohort <- "s669val_"
# clin info
clin <- t.header.valid.lowint
clin$Xcombi_id <- rownames(clin)
# mat from AffyDataset, 669 samples validation cohort with low/int RS, cols=Xcombi_id  rows=AffyIDs:
mat <- log(t(tcdm.valid.lowint))



### TK-EDIT:   SELECT SAMPLE SUBGROUP   ###

### TK-EDIT:   Generate "type" as factor age>=55:old or age<=50:young
clin$type <- NA
clin$type[clin$age>=55] <- "old"
clin$type[clin$age<=50] <- "young"
clin$type <- as.factor(clin$type)

### TK-EDIT:   EXCLUDE SAMPLES with age between 51 and 54 years  (clin$type=NA)
clin <- clin[!is.na(clin$type),]
mat <- mat[,colnames(mat) %in% rownames(clin)]


## Screen for Diff Expr Genes using limma()
#    include type and batch as covariates


type <- clin$type
batch <- factor(clin$datas_new09...TransBIG.40.)

table(type, batch)
#          batch
# type     9 14 23 29 31 37 42 44 47 50 52 53 60 62
# old     94 67 53 21 10  2 40 34  9 24 56  0  1 36
# young   31  3  5  6 22  0 28  2  8 20 34  1  1  1


# Construct model matrix including both  batch and type
X <- model.matrix(~batch+type)

# Use limma to fit gene expression in mat to model matrix
fit <- lmFit(mat,X)

# Empirical Bayes Statistics for Differential Expression
fit <- eBayes(fit)
results <- decideTests(fit)
summary(results)



# volcano plot validation-cohort-A

fn=paste0("./out/",cohort,"volcano-plot.pdf")
pdf(fn,w=8,h=8)
volcanoplot(fit, coef="typeyoung")
plotMD(fit, coef="typeyoung", status=results)
dev.off()


# Select 1000 top probesets as diff expressed genes


deg <- topTable(fit, coef="typeyoung", n=1000)

fn=paste0("./out/",cohort,"DEG-1000list.txt")
write.table(deg,file=fn,row.names=TRUE, col.names = NA, quote=FALSE, sep="\t")

deg2 <- deg   # deg list for validation-cohort


######################################################################################################


# Now compare deg1 and deg2 for list of genes found in both analyses from the different cohorts:
sum(rownames(deg1) %in% rownames(deg2))
# 380

# Combine the 380 genes found in both lists of 1000 deg in one list
deg12 <- deg1[rownames(deg1) %in% rownames(deg2), ]
deg12 <- cbind(deg12, deg2[rownames(deg2) %in% rownames(deg12), ])
colnames(deg12) <- c(paste0(colnames(deg1), "_find"), paste0(colnames(deg2), "_valid"))

# check for same direction of change in both cohorts
deg12$samedir_find_val <- (sign(deg12$logFC_find)==sign(deg12$logFC_valid))

# percentage of gene with same direction
sum(deg12$samedir_find_val) / nrow(deg12)
# 64.7 %

# Save the list of differentially expressed genes
fn=paste0("./out/","DEG-383list_find_and_valid.txt")
write.table(deg12,file=fn,row.names=TRUE, col.names = NA, quote=FALSE, sep="\t")


