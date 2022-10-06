setwd("/Users/taoqing/Box Sync/ERPosExpression/")
library(data.table)
library(dplyr)

#clin=fread("./data/TCGA_BRCA_Clinical.txt",data.table=F,h=T)
load("./data/TCGA_BRCA.RData")

#clin=clin %>%
#     mutate(type=ifelse(Age<=50,"Young",ifelse(Age>=55,"old","other"))) %>%
#     filter(type%in%c("Young","old") & HER2<15.17 & ERstatus=="Positive")
     

table(clin$type)
#age 
tapply(clin$Age,clin$type,summary)

#race
table(clin$Race,clin$type)
chisq.test(clin$Race,clin$type)

#T stage 
table(clin$pathologic_stage,clin$type)
chisq.test(clin$pathologic_stage,clin$type)

#Lymph node 
table(clin$lymphnodes_status,clin$type)
chisq.test(clin$lymphnodes_status,clin$type)

#subtype 
table(clin$Subtype_mRNA,clin$type)
chisq.test(clin$Subtype_mRNA,clin$type)


# menopause
table(clin$menopause_status_fin,clin$type)
chisq.test(clin$menopause_status_fin,clin$type)



clin1=fread("/Users/taoqing/Box Sync/GermlineSomatic (tao.qing@yale.edu)/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/PanCan_ClinicalData_V4_wAIM_filtered10389.txt",stringsAsFactors = F,data.table=F)

clin2=clin1[which(clin1$bcr_patient_barcode%in%clin$bcr_patient_barcode),]



chisq.test(cbind(c(200,66,15),c(426,107,51)))
chisq.test(cbind(c(99,47,135),c(243,111,230)))
chisq.test(cbind(c(168,18,62,33),c(267,238,46,33)))
chisq.test(cbind(c(17,13,55,196),c(30,184,35,335)))

