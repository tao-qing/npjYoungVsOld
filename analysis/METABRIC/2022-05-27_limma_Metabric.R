


# 2022-05-27   LIMMA Analysis of Metabric cohort

# Based on the EDITED VERSION from 2021-11-11: log() changed to log2()  !!!


library(limma)


###############

### The following sample and expression data have to be avaiable:
#  clin.2.lowint
#  mb.expr
#
#  These are provided by script:
#     "2022-05-23_Metabric-Analyses_AgeERpos.R"

###############

### First use the subset of 956 samples with lowint RS (from the ERpos cohort based both on mRNA AND IHC)
###    for the subsequent analysis here:

clin <- clin.2.lowint
rownames(clin) <- clin$SAMPLE_ID2
expr <- mb.expr[,colnames(mb.expr)%in%clin$SAMPLE_ID2] # select expr data for selected samples
expr <- cbind(mb.expr[,c(1:2)], expr)  # re-add gene names and EntrezID, dataframe
mat <- expr
rownames(mat) <- mat[,1]
mat <- as.matrix(mat[,-(1:2)]) # only expression data as matrix with genenames as rownames


################

# A label "cohort" is defined which will be used in the script below to tag any output files

cohort <- "MB-n956-lowintRS_"

################


### TK-EDIT:   SELECT SAMPLE SUBGROUP   ###

### Generate "type" as factor age>=55:old or age<=50:young
clin$type <- NA
clin$type[clin$AGE_AT_DIAGNOSIS>=55] <- "old"
clin$type[clin$AGE_AT_DIAGNOSIS<=50] <- "young"
clin$type <- as.factor(clin$type)

### TK-EDIT:   EXCLUDE SAMPLES with age between 51 and 54 years  (clin$type=NA)
clin <- clin[!is.na(clin$type),]
mat <- mat[,colnames(mat) %in% rownames(clin)]


## Screen for Diff Expr Genes using limma()
#    include type as covariate

type <- clin$type

table(type)
# type
# old young 
# 710   157 

 

# Construct model matrix including type
X <- model.matrix(~type)

# Use limma to fit gene expression in mat to model matrix
fit <- lmFit(mat,X)

# Add Gene names
fit$genes <- row.names(fit)


# Empirical Bayes Statistics for Differential Expression
fit <- eBayes(fit)
results <- decideTests(fit)
summary(results)


# volcano plot 

fn=paste0("./out/",cohort,"volcano-plot.pdf")
pdf(fn,w=8,h=8)
volcanoplot(fit, coef="typeyoung", highlight = 15, names=fit$genes)
plotMD(fit, coef="typeyoung", status=results)
dev.off()

# Summary for typeyoung: 1281 Down, 1913 Up (19089 NotSig) --> total 3194

# Select XXXXXXXX top probesets as diff expressed genes

# deg <- topTable(fit, coef="typeyoung", n=XXXXXXXXXXX)

# fn=paste0("./out/",cohort,"DEG-list.txt")
# write.table(deg,file=fn,row.names=TRUE, col.names = NA, quote=FALSE, sep="\t")

#   Save list of significant DEG from the finding cohort analysis
#   since the variable deg will be overwritten in the subsequent script for cohort-B:
# deg1 <- deg



#################################






######################################################################################################




