setwd("/Users/taoqing/Box Sync/ERPosExpression/")

library(ggplot2)
library(ggsignif)
library(DescTools)
#library(kSamples)
library(dplyr)
library(data.table)
library(tibble)


SCANB_clin = fread("/Users/taoqing/Box Sync/ERPosExpression/data/SCANB/GSE96058_clinicalInformation.txt",data.table = F) %>%
       mutate(Age = gsub("age at diagnosis: ","", Age),
              Tumor_size = gsub("tumor size: ","",Tumor_size),
              LymphNode_Status = gsub("lymph node status: |Node","",LymphNode_Status),
              ER_Status = gsub("er status: ","",ER_Status),
              PR_Status = gsub("pgr status: ","",PR_Status),
              HER2 = gsub("her2 status: ","",HER2),
              KI67 = gsub("ki67 status: ","",KI67),
              PAM50 = gsub("pam50 subtype: ","",PAM50),
              OS = gsub("overall survival days: ","",OS),
              OS_Event = gsub("overall survival event: ","",OS_Event),
              Endocrine_Treated = gsub("endocrine treated: ","",Endocrine_Treated),
              Chemo_Treated = gsub("chemo treated: ","",Chemo_Treated),
              ) %>% 
       filter(!grepl("repl",Samples)) %>%
       filter(ER_Status == 1 & HER2 == 0) %>%
       filter(Age <=50 | Age >= 55) %>%
       mutate(type=ifelse(Age<=50,"Young",ifelse(Age>=55,"Old","other")))
       

       

SCANB_mat = fread("/Users/taoqing/Box Sync/ERPosExpression/data/SCANB/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",data.table = F) %>%
      column_to_rownames(., var = "V1") 
SCANB_mat = SCANB_mat[,pmatch(SCANB_clin$Samples,colnames(SCANB_mat))]


#OncotypeDX score
source("./scripts/oncotypedx3_withoutCTSL2")
library(genefu)
library(biomaRt)

hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")

ensembl <- useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')

#tmp = listAttributes(mart = ensembl)

mapping <- getBM(
  attributes = c('entrezgene_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = rownames(SCANB_mat),
  mart = ensembl
)

#sig.oncotypedx$symbol[which(!sig.oncotypedx$symbol %in% mapping$hgnc_symbol)]
SCANB_mat_sub = SCANB_mat[rownames(SCANB_mat)%in%sig.oncotypedx$symbol,]
mapping = mapping[pmatch(rownames(SCANB_mat_sub),mapping$hgnc_symbol),]
colnames(mapping)=c("EntrezGene.ID","symbol")
rownames(mapping)=mapping$symbol

rs_mb <- oncotypedx3(data= t(SCANB_mat_sub), annot = mapping, do.mapping=FALSE, verbose=TRUE)

clin = SCANB_clin %>%
       mutate(OncotypeDx_score = rs_mb$score[SCANB_clin$Samples]) %>% 
       slice_min(OncotypeDx_score, prop = 0.8) # select bottom 80%

mat = SCANB_mat[,clin$Samples]

save(mat,clin,file="./data/SCANB_breastcancer_data.RData")



load("./data/SCANB_breastcancer_data.RData")

table(clin$type)
#age 
tapply(as.numeric(clin$Age),clin$type,summary)

#race
table(clin$Race,clin$type)
chisq.test(clin$Race,clin$type)

#T stage 
clin$pathologic_stage = ifelse(as.numeric(clin$Tumor_size)>=5, "T3&4",ifelse(as.numeric(clin$Tumor_size)<5,"T1&2","NA"))
table(clin$pathologic_stage,clin$type)
chisq.test(clin$pathologic_stage,clin$type)

table(clin$pathologic_stage,clin$type)[,2]/sum(table(clin$pathologic_stage,clin$type)[,2])
table(clin$pathologic_stage,clin$type)[,1]/sum(table(clin$pathologic_stage,clin$type)[,1])

#Lymph node 
table(clin$LymphNode_Status,clin$type)
chisq.test(clin$LymphNode_Status,clin$type)

table(clin$LymphNode_Status,clin$type)[,2]/sum(table(clin$LymphNode_Status,clin$type)[,2])
table(clin$LymphNode_Status,clin$type)[,1]/sum(table(clin$LymphNode_Status,clin$type)[,1])


#Histological grade, n (%)
clin$Histological_grade = ifelse(clin$NHG%in%c("nhg: G1","nhg: G2"),"G1&2",ifelse(clin$NHG%in%c("nhg: G3"),"G3","NA"))
table(clin$Histological_grade,clin$type)
chisq.test(clin$Histological_grade,clin$type)

table(clin$Histological_grade,clin$type)[,2]/sum(table(clin$Histological_grade,clin$type)[,2])
table(clin$Histological_grade,clin$type)[,1]/sum(table(clin$Histological_grade,clin$type)[,1])

#subtype 
table(clin$PAM50,clin$type)
chisq.test(clin$PAM50,clin$type)

table(clin$PAM50,clin$type)[,2]/sum(table(clin$PAM50,clin$type)[,2])
table(clin$PAM50,clin$type)[,1]/sum(table(clin$PAM50,clin$type)[,1])




#treatment
clin$treatment = ifelse(clin$Chemo_Treated == 1,"Chemotherapy",
                        ifelse(clin$Endocrine_Treated == 1,
                               "Endocrine treatment only",
                               ifelse(clin$Chemo_Treated == 0 & clin$Endocrine_Treated==0, 
                                      "No adjuvant treatment",
                                      "NA")
                        )
                        )

  
  

table(clin$treatment,clin$type)
chisq.test(clin$treatment,clin$type)

table(clin$treatment,clin$type)[,2]/sum(table(clin$treatment,clin$type)[,2])
table(clin$treatment,clin$type)[,1]/sum(table(clin$treatment,clin$type)[,1])

# menopause
#table(clin$menopause_status_fin,clin$type)
#chisq.test(clin$menopause_status_fin,clin$type)


