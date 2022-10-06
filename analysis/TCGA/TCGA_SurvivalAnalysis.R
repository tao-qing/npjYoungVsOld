setwd("/Users/qingtao/Box Sync/ERPosExpression/")
library(maftools)
library(data.table)
library(dplyr)
require('ggfortify')
library(survminer)
library(survival)

source("./TCGASurvival/run_survival_function.R")
source("./TCGASurvival/samstein_functions.R")
source("./TCGASurvival/createOncoMatrix.R")
source("./TCGASurvival/ggkmTable.R")

#prepare data
fn="/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/Clinical_Liu_Cell2018/mmc1.xlsx"
survMat = data.frame(readxl::read_xlsx(fn))


clin=fread("./data/TCGA_BRCA_Clinical.txt",data.table=F,h=T)
clin=clin%>%mutate(type=ifelse(Age<=50,"Young",ifelse(Age>=55,"Old","other")))%>%filter(type%in%c("Young","Old") & HER2<15.17 & ERstatus=="Positive")
clin=clin%>%select("bcr_patient_barcode","ERstatus","HER2status","type")%>%mutate(Tumor_Sample_Barcode=bcr_patient_barcode)


survMat1=survMat%>%filter(bcr_patient_barcode%in%clin$bcr_patient_barcode)%>%mutate(Type=clin$type[pmatch(bcr_patient_barcode,clin$bcr_patient_barcode)])

#survMat1$OS=as.numeric(survMat1$OS)
#survMat1$Type=factor(survMat1$Type)
#survMat1$OS.time=as.numeric(survMat1$OS.time)
#if(any(is.na(survMat1$OS.time))){
#  survMat1=survMat1[-which(is.na(survMat1$OS.time)),]
#}


#tiff("TCGA_youngVsold_OS.tiff",h=2000,w=2300,res=300)
#pdf("TCGA_youngVsold_OS.pdf",h=10,w=10)
p1=plotReadableKM(data=survMat1, predictor="Type", surv_time = "OS.time",status="OS", legend_lbls=NULL, ptitle=NULL)
ggsave("TCGA_youngVsold_OS.pdf",plot=p1$plot,w = 7, h = 7, useDingbat=F,limitsize = FALSE)
#dev.off()

p2=plotReadableKM(data=survMat1, predictor="Type", surv_time = "PFI.time",status="PFI", legend_lbls=NULL, ptitle=NULL)
ggsave("TCGA_youngVsold_PFI.pdf",plot=p2$plot,w = 7, h = 7, useDingbat=F,limitsize = FALSE)

p3=plotReadableKM(data=survMat1, predictor="Type", surv_time = "DFI.time",status="DFI", legend_lbls=NULL, ptitle=NULL)
ggsave("TCGA_youngVsold_DFI.pdf",plot=p3$plot,w = 7, h = 7, useDingbat=F,limitsize = FALSE)

