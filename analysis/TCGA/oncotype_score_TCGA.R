setwd("/Users/taoqing/Box Sync/ERPosExpression/")

rm(list=ls())
gc()

load("./data/TCGA_BRCA.RData")

source("./scripts/oncotypedx1")
source("./scripts/oncotypedx2_withoutBAG1")

library(genefu)

rs_mb <- oncotypedx1(data= t(mat), annot= anno, do.mapping=FALSE, verbose=TRUE)
save(rs_mb,file="./out/TCGA_breast_cancer_OncotypeDX_score.RData")

#without BAG1
rs_mb_noBAG1 <- oncotypedx2(data= t(mat), annot= anno, do.mapping=FALSE, verbose=TRUE)

plot(rs_mb$score,rs_mb_noBAG1$score)





