#germline variant clinical score association
setwd("/Users/Administrator/Box/Huang_lab/manuscripts/CNV_immune/analysis/")
source("./global_aes_out.R")
source("./dependency_files_tq.R")
source("./stat_functions.R")
source("./load_somatic.R")

### MAIN ###
#immune response features for consideration
#response=read.table("/Users/Administrator/Box/Huang_lab/manuscripts/GermlineSomatic/analysis/germline_immune_cov/out/allAssoc_iterms_plot 2",sep="\t",h=F)
response=read.table("/Users/g1w2x3f2012/Box/Huang_lab/manuscripts/GermlineSomatic/analysis/germline_immune_cov/out/allAssoc_iterms_plot 2",sep="\t",h=F)

#immune score
# immuneProfile=read.table("/Users/Administrator/Box/Huang_lab/manuscripts/GermlineSomatic/analysis/germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile=read.table("/Users/g1w2x3f2012/Box/Huang_lab/manuscripts/GermlineSomatic/analysis/germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile=immuneProfile[-which(is.na(immuneProfile$value)),]
immuneProfile=immuneProfile[which(immuneProfile$variable%in%response[,2]),]


# CNV files
  # Taylor data
  library(tidyr)
  library(dplyr)

  cnv_taylor <- read.table(gzfile("C:/Users/g1w2x3f2012/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/CNV_Taylor_CancerCell2018/all_thresholded.by_genes_whitelisted.tsv.gz"), h=T,sep="\t",stringsAsFactors=FALSE)
  # cnvfile <- read.table(gzfile("C:/Users/g1w2x3f2012/box/all_thresholded.by_genes_whitelisted.tsv.gz"), h=T,sep="\t",stringsAsFactors=FALSE)
  cnv_long <- pivot_longer(cnv_taylor, cols = 4:10716, names_to = "bcr_patient_barcode", values_to = "Value") #%>%
  cnv_genes_long = cnv_long %>%
       mutate(bcr_patient_barcode = substr(bcr_patient_barcode,1,12))
  # all the samples from tumor of different patients, no duplicates
  sample_pat_info <- data.frame(sample_id = colnames(cnv_taylor), pat_id = substr(colnames(cnv_taylor),1,12))[4:10716,]
  sample_dup_info <- sample_pat_info[duplicated(sample_pat_info$pat_id),]
  
  cnv_genes_long$bcr_patient_barcode <- gsub("\\.", "-", cnv_genes_long$bcr_patient_barcode)
  cnv_genes_long <- cnv_genes_long[,c(1,4,5)]
  colnames(cnv_genes_long)[1] = "geneSymbol"
  # save(cnv_genes_long, file = "cnv_genes_taylor.rda")
  # load  cnv_gene file
  # load("cnv_genes_taylor.rda")
  
  
  # fn = "../data/cnv_genes_long.txt"
  # write.table(cnv_genes_long, quote=F, sep="\t", file = fn, row.names = F)
  # cnv_genes_long <- read.table(file = fn, h=T,sep="\t",stringsAsFactors=FALSE)
  
  # Taylor arm_level data
  cnv_arm <- read.table("/Users/g1w2x3f2012/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/CNV_Taylor_CancerCell2018/PANCAN_ArmCallsAndAneuploidyScore_092817.txt",h=T,sep="\t",stringsAsFactors=FALSE)
  cnv_arm = data.frame(bcr_patient_barcode = substr(cnv_arm$Sample,1,12), cnv_arm [,3:42],stringsAsFactors = FALSE)
  # cnv_test = data.frame(bcr_patient_barcode = substr(cnvfile$Sample,1,12), cnvfile [,2:42])
  cnv_dup_barcode <- cnv_arm[duplicated(cnv_arm$bcr_patient_barcode),]$bcr_patient_barcode
  cnv_armfile <- cnv_arm [-which(cnv_arm$bcr_patient_barcode%in%cnv_dup_barcode), -2]
  
###### Adj arm level alt#####
  # cytoband annotation by Taylor data
  fn = "../data/cytoband_gene.txt"
  # write.table(cyto_anno, quote=F, sep="\t", file = fn, row.names = F)
  cyto_anno <- read.table(file = fn, h=T,sep="\t",stringsAsFactors=FALSE)
  
  # edit colnames of  armfile & turn into longer form
  colnames(cnv_armfile)[c(26:28,39,40)] <- paste0("X", c(13:15,21,22),"q") 
  cnv_armfile_long <- cnv_armfile %>%
    pivot_longer(cols = 2:40, names_to = "arm", values_to = "arm_value")
  
# Covariates gender & TMB
  # ##### clinical files #####
  # clin_f = "../../../Huang_lab_data/TCGA_PanCanAtlas_2018/tcga_combined_study_clinical_data.tsv"
  clin_f = "../../../Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/clinical_PANCAN_patient_with_followup.tsv"
  # clin_f = "../../../Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/PanCan_ClinicalData_V4_wAIM_filtered10389.txt"
  clin = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_f, stringsAsFactors=FALSE)
  
  tmbMat=immuneProfile[immuneProfile$variable=="TMB",]
  
# filter overlapped samples
  length(unique(immuneProfile$bcr_patient_barcode))
  length(unique(cnv_genes_long$bcr_patient_barcode))
  length(unique(cnv_armfile_long$bcr_patient_barcode))
  length(unique(clin$bcr_patient_barcode))
  length(unique(tmbMat$bcr_patient_barcode))
  
  samples=intersect(immuneProfile$bcr_patient_barcode, cnv_genes_long$bcr_patient_barcode)
  samples3 = intersect(samples, cnv_armfile_long$bcr_patient_barcode)
  samples4 = intersect(samples3, tmbMat$bcr_patient_barcode)
  samples5 = intersect(samples4, clin$bcr_patient_barcode)
  
  label="overlaped_n9336"
  
  clin_selected <- clin[which(clin$bcr_patient_barcode %in% samples5),]
  sexMat <- data.frame(bcr_patient_barcode = clin_selected$bcr_patient_barcode, sex = clin_selected$gender)
  immuneProfile=immuneProfile[which(immuneProfile$bcr_patient_barcode%in%samples5),]
  cnv_long=cnv_genes_long[which(cnv_genes_long$bcr_patient_barcode%in%samples5),]
  
  
# read input files
gene_sample = data.frame(table(cnv_long$bcr_patient_barcode, cnv_long$Gene))
colnames(gene_sample) = c("bcr_patient_barcode","Gene","Freq")
gene_sample=gene_sample[gene_sample$Freq!=0,]

cnv_genes_pos <- cnv_long [abs(cnv_long$Value) > 1, ] %>%
                mutate(AltType = ifelse(Value == 2, "AMP", "DEL")) %>%
                mutate(Gene =  paste(geneSymbol, AltType, sep = "_")) %>%
                mutate(Freq = 1 )
gene_sample_postive <- cnv_genes_pos [, c(1,2,5,6)]


##### individual cancer type analysis #####
clin <- clin_selected
cancers = unique(immuneProfile$TCGA_Study)
genes_cnv = unique(gene_sample_postive$Gene)

rm(clin_selected, cnv_arm, cnv_armfile)
rm(cnv_genes_long, cnv_long)

save.image(paste0("1_CNV_immune_association_fullgenes",".RData"))


#this is a time consuming step, so I move the analysis to a super computer cluster for parallel computating

#setwd("/gpfs/ysm/project/tq37/GermlineSomatic/analysis/germline_immune_cov")
setwd("/home/tq37/project/GermlineSomatic/analysis/germline_immune_cov")
load("1_CNV_immune_association_fullgenes.RData")

library(data.table)
core=32
require(foreach)
require(doParallel)
registerDoParallel(core)

finalMatrix <- foreach(gene=genes_cnv, .combine=rbind) %dopar% {
  tt=NULL
  gene_sample_g = gene_sample_postive [gene_sample_postive$Gene==gene,]
  var_exp_g = merge(immuneProfile,gene_sample_g,by="bcr_patient_barcode",all.x=T)
  var_exp_g$Freq[is.na(var_exp_g$Freq)] = 0
  var_exp_g$Freq[var_exp_g$Freq != 0 ] = 1
  geneSymbol = na.omit(unique(var_exp_g$geneSymbol))
  arm = cyto_anno[which(cyto_anno$geneSymbol == geneSymbol),]$arm
  gene_armvalue = cnv_armfile_long [cnv_armfile_long$arm == arm,]
  for (cancer in cancers){
    var_exp_g_c = var_exp_g[var_exp_g$TCGA_Study %in% cancer,]
    for (clinicPheno in unique(var_exp_g_c$variable)){
      var_exp_g_c_s = var_exp_g_c[var_exp_g_c$variable == clinicPheno,]
      gene_path_count=sum(var_exp_g_c_s$Freq)
      if((dim(var_exp_g_c_s)[1]>10) && sum(var_exp_g_c_s$Freq) > 2){
        pheno=as.character(as.matrix(response[response[,2]==clinicPheno,]))
        df=var_exp_g_c_s
        df$TMB=tmbMat$value[pmatch(df$bcr_patient_barcode,tmbMat$bcr_patient_barcode)]
        df$sex=sexMat$sex[pmatch(df$bcr_patient_barcode,sexMat$bcr_patient_barcode)]
        df=merge(df, gene_armvalue, by = "bcr_patient_barcode")
        if(all(is.na(df$value)) | all(is.na(var_exp_g_c_s$value[var_exp_g_c_s$Freq!=0])) | all(is.na(var_exp_g_c_s$value[var_exp_g_c_s$Freq==0]))){
          next
        }
        
        # run GLM
        w = wilcox.test(var_exp_g_c_s$value[var_exp_g_c_s$Freq==0],var_exp_g_c_s$value[var_exp_g_c_s$Freq!=0])
        wP = w$p.value
        wWstat = w$statistic
        #cancer_gene_stat = run_glm(var_exp_g_c_s,yi="value",xi="Freq",covi=c("Age","Sex","Race"),ytype="Continuous")
        
        if(pheno[4]=="log2"){
          df$value=log2(df$value+as.numeric(pheno[5]))
        }
        
        
        if(clinicPheno=="TMB"){
          #pheno[3]="Poisson"
          if(length(which(df$value==0))>0){
            df=df[df$value!=0,]
          }
          df$value=log2(as.numeric(as.matrix(df$value)))
          #df$value=log2(as.numeric(as.matrix(df$value))+0.001)
        }
        
        #if(){
        #cancer_gene_stat=run_glm(df,yi="value",xi="Freq",covi=c("Age","PC1","PC2"),ytype=pheno[3],gene=gene,cancer=cancer)
        #}else{
          cancer_gene_stat=run_glm(df,yi="value",xi="Freq",covi=c("Age","PC1","PC2","sex","TMB", "arm_value"),ytype=pheno[3],gene=gene,cancer=cancer)
        #}
        #,"TCGA_Subtype"
        if(pheno[3]=="Poisson" | pheno[3]=="Binary"){
          cancer_gene_stat=cancer_gene_stat[-which(colnames(cancer_gene_stat)%in%c("F","z.value","xi_lvl1"))]
        }else{
          cancer_gene_stat=cancer_gene_stat[-which(colnames(cancer_gene_stat)%in%c("F","t.value","xi_lvl1"))]
        }
        colnames(cancer_gene_stat) = c("yi","ytype","xi","df","Deviance","Resid. Df","Resid. Dev","Pr(>F)","Estimate","Std..Error","Pr...t..","covars");
        # compile results
        #full_cancer_gene_stat = cbind(cancer,gene,clinicPheno,gene_path_count,wP,wWstat,cancer_gene_stat)
        full_cancer_gene_stat = cbind(cancer,gene,clinicPheno,gene_path_count,wP,wWstat,cancer_gene_stat)
        tt = rbind(tt, full_cancer_gene_stat)
      }
    }
  }
  print(grep(gene,genes_cnv))
  return(tt)
}


#"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
#tt=finalMatrix
colnames(finalMatrix) = c("cancer","gene","clinicPhenotype","gene_path_count","wilcoxP","W_stat","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance","p-value(Anova)","coefficient","StdError","p-value","covariants");
finalMatrix_overlap9035 <- finalMatrix %>%
                          separate(col = gene, into = c("AlterType","Gene")) %>%
                          mutate(FDR = p.adjust(finalMatrix[,"p-value"], method="BH"))%>%
                          arrange(FDR)
  
  
  
# finalMatrix $ FDR = p.adjust(finalMatrix[,"p-value"], method="BH")
# finalMatrix=finalMatrix[order(finalMatrix$FDR, decreasing=FALSE),]
save(finalMatrix_overlap9035, file = "CNV_ImmuneAssoc_103gene_overlap9035.rda")

#colnames(tt) = c("cancer","gene","clinicPhenotype","gene_path_count","wilcoxP","W_stat","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance","F_statistic","p-value(Anova)","coefficient","StdError","p-value","covariants");
#tt$FDR = p.adjust(tt[,"p-value"], method="fdr") # MAW new, calculates FDR based on the method from,
# Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57, 289â€?300.
#tt$wilcoxFDR = p.adjust(tt[,"wilcoxP"], method="fdr")



  
tn = paste0("./out/CNV_ImmuneAssoc_genelevel","byCancerType_cov.txt")
write.table(finalMatrix, quote=F, sep="\t", file = tn, row.names = F)

tn =paste0("./out/CNV_ImmuneAssoc_gene103_",label, "_byCancerType_cov.txt")
write.table(finalMatrix_overlap9035, quote=F, sep="\t", file = tn, row.names = F)



#### annotate cytoband
# Read data
cnv_genome <- read.table(gzfile("./out/CNV_ImmuneAssoc_fullgenes_overlaped_n9336_allcancerType.txt.gz"),sep="\t",h=T, stringsAsFactors = FALSE)
cnv_taylor <- read.table(gzfile("C:/Users/g1w2x3f2012/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/CNV_Taylor_CancerCell2018/all_thresholded.by_genes_whitelisted.tsv.gz"), h=T,sep="\t",stringsAsFactors=FALSE)
