#setwd("/Users/taoqing/Box Sync/ERPosExpression/2022-05-10_Metabric/")
setwd("/Users/taoqing/Box Sync/ERPosExpression/2022-05-10_Metabric/")

rm(list=ls())
gc()

load("./data/Metabric_data_for_oncotype.RData")

source("./oncotypedx2_withoutBAG1")

require(genefu)

anno=as.data.frame(anno)

sig.oncotypedx = sig.oncotypedx[sig.oncotypedx$symbol!="BAG1",]

rs.mb <- oncotypedx2(data= t(mat1), annot= anno, do.mapping=FALSE, verbose=TRUE)

save(rs.mb,file="./out/Metabric_OncotypeDX_score.RData")



