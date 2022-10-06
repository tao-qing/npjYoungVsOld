
##############################################################
###
###  2022-05-10
###
###  Validation in METABRIC cohort
###
##############################################################


library(data.table)
library(dplyr)

# library(xlsx)
######## TK-EDIT: library(xlsx) ERROR Java-dependency !
######## read.xlsx() will be replaced by tidyverse read_excel()

library(ggplot2)
th = theme_bw() + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black")
) + theme(
  legend.position = 'none',
  axis.text.x = element_text(colour = "black", size = 12),
  axis.text.y = element_text(
    colour = "black",
    size = 12,
    hjust = 0.95
  ),
  axis.ticks = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "bold"),
  axis.title = element_text(size = 12, face = "bold")
) + theme(
  strip.placement = "outside",
  plot.title = element_text(hjust = 0.5),
  strip.text.x = element_text(size = 12),
  strip.text.y = element_text(angle = 0, size = 10, face = "italic"),
  panel.spacing = unit(0.1, "lines"),
  strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")
)


mat= fread("/Users/taoqing/Box Sync/ERPosExpression/2022-05-10_Metabric/data/brca_metabric/data_mrna_agilent_microarray_zscores_ref_all_samples.txt",data.table = F)

# import data

# Read metabric(mb) mRNA expression data from Illumina microarrays (800MB) 
#mb.expr <- read.delim("./data/data_expression.txt",h=T)
mb.expr <- fread("./data/data_expression.txt",data.table = F, colClasses=c(key='character')) 


# Read metabric(mb) clinical data file skipping first 5 rows with comments on variables
mb.clin <- read.delim("./data/data_clinical.txt",h=T, skip=5)


# import of mb.expr led to change of column headers vs SAMPLE_ID.
# Therefor we need a new variable SAMPLE_ID2 :
mb.clin$SAMPLE_ID2 <- gsub("-", ".", mb.clin$SAMPLE_ID)





sum(mb.clin$SAMPLE_ID2 %in% colnames(mb.expr))
# for all 1980 samples with clinical information, are expression data available



#########################

# A label "cohort" is defined which will be used in the script below to tag any output files

cohort <- "Metabric_"

# clin info
clin.all <- mb.clin
rownames(clin.all) <- mb.clin$SAMPLE_ID2

# mat from mb.expr:
mat.all <- as.matrix(mb.expr[,-(1:2)])
rownames(mat.all) <- mb.expr[,1]


#prepare data for genefu RS analysis by Tao
anno = cbind(Hugo_Symbol=mb.expr$Hugo_Symbol,EntrezGene.ID=mb.expr$Entrez_Gene_Id)
mat1=round(mat.all,digits=4)
save(mat1,anno,file="./data/Metabric_data_for_oncotype.RData")



################

## SELECTION OF ER POS TUMORS:

### Select 1257 ERpos tumors applying THREEGENE variable (ThreeGene model of Haibe-Kains)
clin.1 <- mb.clin[mb.clin$THREEGENE %in% c("ER+/HER2- High Prolif", "ER+/HER2- Low Prolif"), ]
mat.1 <- mat.all[ , colnames(mat.all) %in% clin.1$SAMPLE_ID2]


# Select subset of 1196  with ER-pos-IHC from clin.1
clin.2 <- na.omit(clin.1[clin.1$ER_IHC=="pos", ])
mat.2 <- mat.all[ , colnames(mat.all) %in% clin.2$SAMPLE_ID2]





###############

### First use the subset of 1196 (clin.2) for the subsequent analysis here:

clin <- clin.2
rownames(clin) <- clin$SAMPLE_ID2
expr <- mb.expr[,colnames(mb.expr)%in%clin.2$SAMPLE_ID2] # select expr data for selected samples
expr <- cbind(mb.expr[,c(1:2)], expr)  # re-add gene names and EntrezID, dataframe
mat <- expr
rownames(mat) <- mat[,1]
mat <- as.matrix(mat[,-(1:2)]) # expression data as matrix with genenames as rownames

################


##  Exclusion of high RS-score needed here  !!!!!!!!!!!!!!!!!!!
# library(genefu)
# # Syntax of function oncotypedx: # oncotypedx(data, annot, do.mapping = FALSE, mapping, verbose = FALSE) 
# mb.annot <- mb.expr[,c(1,2)]
# colnames(mb.annot) <- c("Hugo_Symbol", "EntrezGene.ID")
# rs.mb <- oncotypedx(data= expr, annot= mb.annot, do.mapping=FALSE, verbose=TRUE)
###
### DOES NOT YET WORK !!!

#####################################



### TK-EDIT:   SELECT SAMPLE SUBGROUP   ###

### Generate "type" as factor age>=55:old or age<=50:young
clin$type <- NA
clin$type[clin$AGE_AT_DIAGNOSIS>=55] <- "old"
clin$type[clin$AGE_AT_DIAGNOSIS<=50] <- "young"
clin$type <- as.factor(clin$type)

### TK-EDIT:   EXCLUDE SAMPLES with age between 51 and 54 years  (clin$type=NA)
clin <- clin[!is.na(clin$type),]
mat <- mat[,colnames(mat) %in% rownames(clin)]


#immune signatures genes
# Ori-script:
# signature=read.xlsx("./data/collected_signatures.xlsx",sheetIndex = "R")
library(readxl)
signature <- read_excel("data/tk_collected_signatures.xlsx", sheet = "R")

### !!! REMARK: It is critical that there are no white spaces between genes and commas of the genelists !
###   Otherwise after 'strsplit' (see below) the will be spaces in genename strings
###          and these genes will not be found in expression data!



sig_lists=list()
for(i in 1:dim(signature)[1]){
  sig_lists[as.character(signature$Short[i])]=list(as.character(strsplit(as.character(signature$Genes[i]),split=",")[[1]]))
}



#estimate the expression signatures by taking the average expression of all member genes in a signature for each individual patient.
# TK: na.rm=TRUE included!
mat_signature=NULL
for(sig in names(sig_lists)){
  score=apply(mat[which(rownames(mat)%in%sig_lists[[sig]]),],2,mean,na.rm=TRUE)
  mat_signature=rbind(mat_signature,score)
}


rownames(mat_signature) <- names(sig_lists)


### ADD MKI67 and ESR1 as single genes to the mat_signature
ESR1 <- mat["ESR1",]
MKI67 <- mat["MKI67",]
mat_signature <- rbind(mat_signature, ESR1, MKI67)


#zscore normalized signatures
mat_signature_zscore=as.data.frame(apply(mat_signature,1,function(x)(x-mean(x))/sd(x)),na.rm=TRUE)
mat_signature_zscore=cbind(rownames(mat_signature_zscore), mat_signature_zscore)
colnames(mat_signature_zscore)[1] <- "SAMPLE_ID2"

finalMat=merge(clin,mat_signature_zscore,by="SAMPLE_ID2")



#pairs scatter plot: pairwise comparison
#this part will generate the figure in the slide #9
library(psych)
# Ori-script:
#  pdf("./out/signature_correlation.pdf",w=8,h=8)
# TK-REPLACEMENT:  include cohort-tag in filename
fn=paste0("./out/",cohort,"signature_correlation.pdf")
pdf(fn,w=8,h=8)
# Ori-script:
## pairs.panels(finalMat%>%select("MKS","TIS","TCell","BCell","DendriticCell","MastCell","ERS","ERS_liminal","ERS_Pos_Symmans","ERS_Neg_Symmans")%>%mutate_all(~as.numeric(as.matrix(.))), 
##   TK-Edit: Changed "ERS_liminal" to "ERS_luminal" in file "tk_collected_signatures.xlsx" !!
##              Added MMKI67 and ESR1 as single genes
pairs.panels(finalMat%>%select("MKS","TIS","TCell","BCell","DendriticCell","MastCell","ERS","ERS_luminal","ERS_Pos_Symmans","ERS_Neg_Symmans","MKI67","ESR1")%>%mutate_all(~as.numeric(as.matrix(.))), 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             cex=1
)
dev.off()


############################


#compare signatures in young vs old using Mann -Whitney test, donot adjust by any covariates.
##this part will generate the table in the slide #7
ttMW=NULL
# Ori-script:
#  for(s in as.character(signature$Short)){
# TK-REPLACEMENT:  as.character(signature$Short))  REPLACED BY  rownames(mat_signature)
# because of added ESR1 single gene to mat_signature
for(s in rownames(mat_signature)){
  submat=finalMat[,c("type",s)]
  submat[,s]=as.numeric(as.matrix(submat[,s]))
  
  pvalue=wilcox.test(submat[submat$type=="young",s],submat[submat$type=="old",s])$p.value
  tmpmean=tapply(submat[,s],submat$type,mean)
  ttMW=rbind(ttMW,cbind(Signature=s,'log2 Fold Change'=as.numeric(tmpmean["young"])-as.numeric(tmpmean["old"]),'Mean in Young'=tmpmean["young"],'Mean in Old'=tmpmean["old"],Pvalue=pvalue))
}

ttMW=as.data.frame(ttMW)
ttMW$Pvalue=as.numeric(as.matrix(ttMW$Pvalue))
ttMW$FDR=p.adjust(ttMW$Pvalue,method="fdr")

ttMW$`log2 Fold Change`=as.numeric(as.matrix(ttMW$`log2 Fold Change`))
ttMW$`Mean in Young`=as.numeric(as.matrix(ttMW$`Mean in Young`))
ttMW$`Mean in Old`=as.numeric(as.matrix(ttMW$`Mean in Old`))


######## TK-EDIT: library(xlsx) fail (ERROR Java-dependency) !
#  write.xlsx(ttMW,"./out/TCGA_BRCA_Signature_YoungVsOld.xlsx")
######## TK-REPLACEMENT: using write.table() instead of write.xlsx for output
#       and include cohort-tag in filename
fn=paste0("./out/",cohort,"TCGA_BRCA_Signature_YoungVsOld.txt")
write.table(ttMW, file=fn, row.names=TRUE, col.names = NA, quote=FALSE, sep="\t")






#violion plot compare signatures in young vs old: p value was generated by Mann -Whitney test
##this part will generate the figures in the slide #8
# Ori-script:
#  for(s in as.character(signature$Short)){
# TK-REPLACEMENT:  as.character(signature$Short))  REPLACED BY  rownames(mat_signature)
# because of added ESR1 single gene to mat_signature
for(s in rownames(mat_signature)){
  submat=finalMat[,c("type",s)]
  colnames(submat)[2]="Signature"
  med <- submat %>%
    group_by(type) %>% 
    filter(!is.na(type)) %>%
    summarize(med = median(as.numeric(as.matrix(Signature))))
  
  if(ttMW$Pvalue[ttMW$Signature==s]<0.0001){
    pvalue<-formatC(ttMW$Pvalue[ttMW$Signature==s], format = "e", digits = 2)
  }else{
    pvalue=round(ttMW$Pvalue[ttMW$Signature==s],digits=4)
  }
  
  if(ttMW$FDR[ttMW$Signature==s]<0.0001){
    FDR<-formatC(ttMW$FDR[ttMW$Signature==s], format = "e", digits = 2)
  }else{
    FDR=round(ttMW$FDR[ttMW$Signature==s],digits=4)
  }
  
  p=ggplot(submat,
           aes(x=type,y=as.numeric(as.matrix(Signature)))) + 
    geom_violin() + 
    geom_jitter(alpha=.1,width = .1,size=.2)+th+ggtitle(paste0("p value=",pvalue,"\nFDR=",FDR))+xlab("")+ylab(s)+
    geom_segment(data = med, 
                 aes(y=med,yend=med,x=as.numeric(type)-.2,
                     xend=as.numeric(type)+.2),color="red")
  # Ori-script:
  # fn=paste0("./out/",s,"_wilcox.pdf")
  # TK-REPLACEMENT:  include cohort-tag in filename
  fn=paste0("./out/",cohort,s,"_wilcox.pdf")
  ggsave(fn,w = 3, h = 3, useDingbat=F,limitsize = FALSE)
}





#include clinical covariates in our analysis
#we didn't show the following analysis in our slides.
#estimate the association between age group and expression signatures, adjusted by covariates. model: gene-expression signature ~ age group (young vs old) + race + lymphnodes status + pathologic stage
tt=NULL
# Ori-script:
#  for(s in as.character(signature$Short)){
# TK-REPLACEMENT:  as.character(signature$Short))  REPLACED BY  rownames(mat_signature)
# because of added ESR1 single gene to mat_signature
for(s in rownames(mat_signature)){
  # Ori-script:
  #  model=paste0("as.numeric(",s,")~type+Race+lymphnodes_status+pathologic_stage")#+race+number_of_lymphnodes_positive_by_he+pathologic_stage
  # TK-Edit:  REPLACE WITH CLINICAL COVARIATES AVAILABLE FOR AFFY DATASET !!!!
  model=paste0("as.numeric(",s,")~type+RACE_code.2015_05_26+nodal+t1vsrest+g12vs3")  # + race + lymphnodes_pos/neg + T1vs.T2T3T4 + HistolGrade3vs.HG1HG2
  fit=glm(formula=model,data=finalMat,family=gaussian(link = "identity"))
  
  #  stats = data.frame(summary(fit)$coefficients)["typeold",]
  ######## TK-EDIT: This may be an error in the script, "typeold" not in "fit", use "typeyoung" instead:
  stats = data.frame(summary(fit)$coefficients)["typeyoung",]
  
  tt=rbind(tt,stats)
}

rownames(tt)=NULL
# Ori-script:
#  tt=cbind(Signature=as.character(signature$Short),tt)
# TK-REPLACEMENT:  as.character(signature$Short))  REPLACED BY  rownames(mat_signature)
# because of added ESR1 single gene to mat_signature
tt=cbind(Signature=rownames(mat_signature),tt)
tt$Signature=as.character(tt$Signature)
tt$FDR=p.adjust(tt$Pr...t..,method="fdr")
# Ori-script:
# write.table(tt,"./out/TCGA_ERPos_Age_Signature_Association.txt",sep="\t")
# TK-REPLACEMENT:  include cohort-tag in filename
fn=paste0("./out/",cohort,"TCGA_ERPos_Age_Signature_Association.txt")
write.table(tt,file=fn,row.names=TRUE, col.names = NA, quote=FALSE, sep="\t")



#violion plot compare signatures in young vs old: p value was generated by the regression analysis adjusted by covariates 
# Ori-script:
#  for(s in as.character(signature$Short)){
# TK-REPLACEMENT:  as.character(signature$Short))  REPLACED BY  rownames(mat_signature)
# because of added ESR1 single gene to mat_signature
for(s in rownames(mat_signature)){
  submat=finalMat[,c("type",s)]
  colnames(submat)[2]="Signature"
  pvalue=round(tt[tt$Signature==s,"Pr...t.."],digits=4)
  FDR=round(tt[tt$Signature==s,"FDR"],digits=4)
  
  med <- submat %>%
    group_by(type) %>% 
    filter(!is.na(type)) %>%
    summarize(med = median(as.numeric(as.matrix(Signature))))
  
  p=ggplot(submat,
           aes(x=type,y=as.numeric(as.matrix(Signature)))) + 
    geom_violin() + 
    geom_jitter(alpha=.1,width = .1,size=.2)+th+ggtitle(paste0("p value=",pvalue,"\nFDR=",FDR))+xlab("")+ylab(s)+
    geom_segment(data = med, 
                 aes(y=med,yend=med,x=as.numeric(type)-.2,
                     xend=as.numeric(type)+.2),color="red")
  # Ori-script:
  #  fn=paste0("./out/",cohort,s,".pdf")
  # TK-REPLACEMENT:  include cohort-tag in filename
  fn=paste0("./out/",cohort,s,".pdf")
  ggsave(fn,w = 3, h = 3, useDingbat=F,limitsize = FALSE)
}



########################################################################################################




# Cross tables needed for slide 3:

# set outfile:
fn=paste0("./out/",cohort,"CrossTables_for_slide-3.txt")
sink(fn)


# numbers old/young
table(clin$type)

# no info found in clin for ethnicity
# # Ethnicity
# # 1=Caucasian/White, 2=AfricanAmerican/Black, 3=Non_Caucas_Non_AfricAmer, 4=Asian, 5=Hispanic,   9=n.a./mixed
# table(clin$type, clin$RACE_code.2015_05_26)
# chisq.test(table(clin$type, clin$RACE_code.2015_05_26))
# # number of missings:
# table(clin$type, is.na(clin$RACE_code.2015_05_26))




# T stage
# 1=T1, 2=T2/T3/T4
table(clin$type, clin$TUMOR_STAGE)
chisq.test(table(clin$type,clin$TUMOR_STAGE))
# number of missings:
table(clin$type, is.na(clin$TUMOR_STAGE))



# lymph node
# 0=LNN 1=N+
table(clin$type, clin$LYMPH_NODES_POSITIVE)
chisq.test(table(clin$type, clin$LYMPH_NODES_POSITIVE))
# number of missings:
table(clin$type, is.na(clin$LYMPH_NODES_POSITIVE))



# age:
by(clin$AGE_AT_DIAGNOSIS, clin$type, median)
by(clin$AGE_AT_DIAGNOSIS, clin$type, range)



sink()   # stop output into file












