
##############################################################
###
###  2021-03-12
###
###  RE-DO Complete Analyses with n2007 starting cohort
###
##############################################################



# import data

# New Version using txt-file import:
#  and changes to new R versions based on DataInBrief paper and suppl.:

cd.compl <- read.delim("2021-03-12_AgeERpos_MAS5_combined.txt",h=T)
# much faster read (<1min) than read_excel (6min) for 400MB file

t.header <- read.delim("2021-03-12_AgeERpos_MAS5_SampleInfo.txt",h=T)

rownames(cd.compl) <- cd.compl[,1]     # use first column with Affy IDs as rownames
cd.compl <- cd.compl[,-1]

rownames(t.header) <- paste("X",t.header$combi_id, sep="")

n.col.cd=ncol(cd.compl)
n.probes=nrow(cd.compl)

tdsn=as.numeric(t.header$datas_new09)   # transposed dataset vector
tcdm=t(cd.compl)   # transposed chip data matrix


##### Figure TK-1:  Selection of the finding cohort from multiple datasets based on dataset comparibility

#  Select comparable samples from n=2007 BC with MAS5 Affymetrix data

#  Variables used:
#  cd.compl ==> chip-data complete (data.frame, 2007 samples in columns, 22283 probesets in rows)
#  t.header ==> transposed sample info including dataset allocation (data.frame, 2007 rows)
#  n.probes ==> number of ProbeSets (rows in cd.compl)
#  n.col.cd ==> number of columns/samples in cd.compl
#  tdsn  ==> transposed dataset allocation (numeric vector)
#  tcdm  ==> transposed chip data matrix
#


ds.mean=by(tcdm, tdsn, colMeans)  # means of each probeset within x individual datasets, 
#  List is sorted by numeric dataset tdsn
tcdm.mean=apply(tcdm,2,mean)
#   generates named list of 22283 global means, which can be indexed by tcdm.mean[probes]
tcdm.stdev=apply(tcdm,2, sd)   # the same for the StdDev


# Calculate comparability metrics for the datasets:
# (calculate sum of squared differences of dataset-mean from total-mean for all probesets)
# Define variables
n.datas=length(ds.mean)   # number of datasets
diff.to.mean= matrix(0,nrow=n.probes,ncol=n.datas) # matrix of differences from mean
nrm.diff.to.mean= matrix(0,nrow=n.probes,ncol=n.datas) # matrix of NORMALIZED diff from mean

for (probes in 1:n.probes)   # loop for all probesets
{
  
  for (i in 1:n.datas)
    # calculate for each dataset diff from global-mean of all datasets
    # and save in matrix "diff.to.mean"
  {
    diff.to.mean[probes,i]=ds.mean[[i]][probes]- tcdm.mean[probes]
    # again the same but normalize by dividing through StdDev
    nrm.diff.to.mean[probes,i]=(ds.mean[[i]][probes]- tcdm.mean[probes]) / tcdm.stdev[probes]
  }
}
# calculate squares of differences
squ.diff=diff.to.mean^2
squ.nrm.diff= nrm.diff.to.mean^2
# sum of squared differences by column
sum.squ.diff=apply(na.omit(squ.diff),2,sum)
sum.squ.nrm.diff=apply(na.omit(squ.nrm.diff),2,sum)
#  important: these vectors are still sorted by numeric dataset tdsn !

# summarize results
comparab=data.frame(sort(unique(tdsn)),sum.squ.diff,sum.squ.nrm.diff)
names(comparab)=c("dataset","sum.squ.diff","sum.squ.nrm.diff")
sort.comparab=comparab[order(comparab$sum.squ.nrm.diff),]

# integrate normalized comparab data in sample info in t.header
for (i in 1:n.col.cd) 
{t.header$comparab_nrm[i]= comparab$sum.squ.nrm.diff[comparab$dataset==tdsn[i]]}
# remove temporary variables:
rm(diff.to.mean, nrm.diff.to.mean, ds.mean,tcdm.mean, tcdm.stdev, squ.diff, squ.nrm.diff, sum.squ.diff, sum.squ.nrm.diff)


# plot comparability-figure

fn=paste0("./out/","Fig-TK-1_Array-Comparab.pdf")
pdf(fn,w=8,h=8)

  plot(sort(t.header$comparab_nrm),type="l")
  abline(5000,0,col="red")
  plot(sort(t.header$comparab_nrm),type="l",ylim =c(0,20000))
  abline(5000,0,col="red")

dev.off()


########################################################

# SELECT FINDING COHORT OF BC with comparability metric <5000

#  Select a subset of datasets with low comparability metric
# Select a subset of comparab by defining criteria:
compar.subset= subset(comparab, subset= sum.squ.nrm.diff < 5000)
# vector of corresponding datasets:
datas.subset=compar.subset$dataset
# generate logical vector FALSE/TRUE for the complete dataset of 2007 BC:
subset.index.vector=(tdsn %in% datas.subset)
# query selected samples from transposed chipdata matrix:
tcdm.find=tcdm[subset.index.vector , ]
# query  corresponding transposed dataset vector:
tdsn.find=tdsn[subset.index.vector]
# query selected samples from NOT-transposed dataset (and corresponding header and t.header):
cd.compl.find= cd.compl[,subset.index.vector]
t.header.find= t.header[subset.index.vector,]

##### Figure TK-2: Analysis of a potential dataset bias among probesets in FINDING COHORT:

# Calculate Kruskal-Wallis statistics for all probesets according to their association
#   with the dataset vector among the finding cohort samples:

# Variables used:    tcdm.find, tdsn.find, n.probes, cd.compl.find
# Newly generated variables:    kruskal.param

kruskal.result=kruskal.test(tcdm.find[,1],tdsn.find)  # Define variable as "list" by performing test once
kruskal.stat= matrix(0,nrow=n.probes)  # matrix for chi-statistics
kruskal.p=matrix(0,nrow=n.probes)  # matrix for p-values

# loop perfoming kruskal-test for each probeset (chip-data vs. dataset-vector)
# the results of each test are variables of type LIST, which are combined by indexing with [[i]]
for (i in 1:n.probes)
{
  kruskal.result[[i]]=kruskal.test(tcdm.find[,i], tdsn.find)   # results as LIST
  kruskal.stat[i]=kruskal.result[[i]]$statistic  # read statistics from LIST in matrix
  kruskal.p[i]=kruskal.result[[i]]$p.value   # read p-values from LIST in matrix
}
# summarize statistics and p-values and probeset names in one dataframe:
kruskal.param.find= data.frame(row.names=rownames(cd.compl.find),kruskal.stat,kruskal.p)  # dataframe with results

# remove temporary variables:
rm(kruskal.result, kruskal.stat, kruskal.p)

hist(kruskal.param.find$kruskal.stat,breaks=30,col="dark blue") # histogram of kruskal wallis stat for all probesets

# histogram of kruskal wallis stat for probes from hemoglobin metagene:


### optional plot figure in pdf file:
## fn=paste0("./out/","Fig-TK-2_Example-dataset-bias-hemoglobin-metagene.pdf")
## pdf(fn,w=8,h=8)

  Hemoglobin.meta.probes=c("204419_x_at", "204848_x_at", "209116_x_at", "204018_x_at", "209458_x_at", "211745_x_at", "214414_x_at", "211696_x_at", "217232_x_at", "217414_x_at", "211699_x_at", "213515_x_at")
  hist(kruskal.param.find[rownames(kruskal.param.find) %in% Hemoglobin.meta.probes,]$kruskal.stat,breaks=30,col="dark blue")

## dev.off()



##### Figure 3: Analysis of a potential dataset bias among probesets in FULL COHORT:


##### ZUSÄTZLICH KRUSKAL-WALLIS AUCH NOCH IN FULL COHORT DURCHFÜHREN !!!!!


# Calculate Kruskal-Wallis statistics for all probesets according to their association
#   with the dataset vector among the finding cohort samples:

# Variables used:    tcdm, tdsn, n.probes, cd.compl

# Newly generated variables:    kruskal.param.full2007


kruskal.result=kruskal.test(tcdm[,1],tdsn)  # Define variable as "list" by performing test once
kruskal.stat= matrix(0,nrow=n.probes)  # matrix for chi-statistics
kruskal.p=matrix(0,nrow=n.probes)  # matrix for p-values

# loop perfoming kruskal-test for each probeset (chip-data vs. dataset-vector)
# the results of each test are variables of type LIST, which are combined by indexing with [[i]]
for (i in 1:n.probes)
{
  kruskal.result[[i]]=kruskal.test(tcdm[,i], tdsn)   # results as LIST
  kruskal.stat[i]=kruskal.result[[i]]$statistic  # read statistics from LIST in matrix
  kruskal.p[i]=kruskal.result[[i]]$p.value   # read p-values from LIST in matrix
}
# summarize statistics and p-values and probeset names in one dataframe:
kruskal.param.full2007= data.frame(row.names=rownames(cd.compl),kruskal.stat,kruskal.p)  # dataframe with results

# remove temporary variables:
rm(kruskal.result, kruskal.stat, kruskal.p)

hist(kruskal.param.full2007$kruskal.stat,breaks=30,col="dark blue") # histogram of kruskal wallis stat for all probesets

# histogram of kruskal wallis stat for probes from hemoglobin metagene:

Hemoglobin.meta.probes=c("204419_x_at", "204848_x_at", "209116_x_at", "204018_x_at", "209458_x_at", "211745_x_at", "214414_x_at", "211696_x_at", "217232_x_at", "217414_x_at", "211699_x_at", "213515_x_at")
hist(kruskal.param.full2007[rownames(kruskal.param.full2007) %in% Hemoglobin.meta.probes,]$kruskal.stat,breaks=30,col="dark blue")

########################################################



#################################  RecurrenceScore von Finding cohort (n1170)

library(genefu)

# Syntax der function oncotypedx: # oncotypedx(data, annot, do.mapping = FALSE, mapping, verbose = FALSE) 
# Note that for Affymetrix HGU datasets, the mapping is not necessary.
rs.find <- oncotypedx(data= tcdm.find, annot= U133A_Affy_Entrez, do.mapping=FALSE, verbose=TRUE)

# oncotypedx-score in t.header.find integrieren
t.header.find <- cbind(t.header.find,rs.find$score)

# distribution of recurrence-scores from genefu-oncotypedx-function:
hist(rs.find$score,breaks = 50)
quantile(rs.find$score)
quantile(rs.find$score, probs = seq(0, 1, 0.2))




# clinical parameters of finding cohort
table(t.header.find$dataset)
table(t.header.find$nodal)
table(t.header.find$age50)
table(t.header.find$t1vsrest)
table(t.header.find$grade_1234)
table(t.header.find$g1vs23)
table(t.header.find$g12vs3)
table(t.header.find$PGR_CO_MDAN_min0.0078)
table(t.header.find$MSubTypeHugh)
table(t.header.find$treat3class)
table(t.header.find$Neoadj_Treat)
table(t.header.find$pCR)

# Pragmatic selection of Low/Interm-Risk:  <80%-Quantile of rs.find$score


# SELECT LowInt-subset of finding cohort

# Select a subset of datasets with rs.find$score < 80%-Quantile 
# Select a subset of t.header.find by defining criteria:
t.header.find.lowint <- subset(t.header.find, subset= rs.find$score < as.numeric(quantile(rs.find$score, probs=0.8)))
# Select the subset of corresponding samples(rows) 
#   from transposed dataset of the finding cohort (tcdm.find)
tcdm.find.lowint <- tcdm.find[rownames(tcdm.find) %in% rownames(t.header.find.lowint),]
dim(tcdm.find.lowint) # 936 samples

############################


# Adaption of script "compare_expression_signature.R" from Tao Qing for AffyDataset  


# The following data from the script above will be needed:
#   t.header.find.lowint   : sample info for 948 comparable ERposHer2neg
#   tcdm.find.lowint  : matrix of MAS5 Affy data with the 948 samples as rows(!) and 22283 Probesets in cols(!)

# We will also retain the following data for a potential second validation cohort:
#   tcdm
#   t.header
#   tcdm.find
#   t.header.find

# We delete the following data to free memory:
#   cd.compl
#   cd.compl.find
rm(cd.compl, cd.compl.find)



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


#########################

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


#immune signatures genes
# Ori-script:
# signature=read.xlsx("./data/collected_signatures.xlsx",sheetIndex = "R")
# TK-REPLACEMENT:  ALTERED file and additional sheet "R_affy"  (+ using read_excel function):
#    !NOTE: "MKI67" as "single gene" was added here as "signature" since represented by 4 probesets
#           which do not cause error when using "apply" to calculate mean below!
library(readxl)
signature <- read_excel("data/tk_collected_signatures.xlsx", sheet = "R_affy")


# Ori-script:
# sig_lists=list()
# for(i in 1:dim(signature)[1]){
#   sig_lists[as.character(signature$Short[i])]=list(as.character(strsplit(as.character(signature$Genes[i]),split=",")[[1]]))
# }

# TK-REPLACEMENT:   Replace "Genes" with "Probesets"
sig_lists=list()
for(i in 1:dim(signature)[1]){
  sig_lists[as.character(signature$Short[i])]=list(as.character(strsplit(as.character(signature$Probesets[i]),split=",")[[1]]))
}

#estimate the expression signatures by taking the average expression of all member genes in a signature for each individual patient.
mat_signature=NULL
for(sig in names(sig_lists)){
  score=apply(mat[which(rownames(mat)%in%sig_lists[[sig]]),],2,mean)
  mat_signature=rbind(mat_signature,score)
}
rownames(mat_signature)=names(sig_lists)

### TK-EDIT:  ADD ESR1 as single gene to the mat_signature
ESR1 <- mat["205225_at",]
mat_signature <- rbind(mat_signature, ESR1)


# Ori-script:
#zscore normalized signatures
# mat_signature_zscore=as.data.frame(cbind(bcr_patient_barcode=colnames(mat_signature),apply(mat_signature,1,function(x)(x-mean(x))/sd(x))))
# TK-REPLACEMENT:
mat_signature_zscore=as.data.frame(cbind(Xcombi_id=colnames(mat_signature),apply(mat_signature,1,function(x)(x-mean(x))/sd(x))))


# Ori-script:
# finalMat=merge(clin,mat_signature_zscore,by="bcr_patient_barcode")
# TK-REPLACEMENT:   merge by Xcombi_id
finalMat=merge(clin, mat_signature_zscore, by="Xcombi_id")




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
########################################################################################################




# Cross tables needed for slide 3:

# set outfile:
fn=paste0("./out/",cohort,"CrossTables_for_slide-3.txt")
sink(fn)


# numbers old/young
table(clin$type)


# Ethnicity
# 1=Caucasian/White, 2=AfricanAmerican/Black, 3=Non_Caucas_Non_AfricAmer, 4=Asian, 5=Hispanic,   9=n.a./mixed
table(clin$type, clin$RACE_code.2015_05_26)
chisq.test(table(clin$type, clin$RACE_code.2015_05_26))
# number of missings:
table(clin$type, is.na(clin$RACE_code.2015_05_26))




# T stage
# 1=T1, 2=T2/T3/T4
table(clin$type, clin$t1vsrest)
chisq.test(table(clin$type,clin$t1vsrest))
# number of missings:
table(clin$type, is.na(clin$t1vsrest))



# lymph node
# 0=LNN 1=N+
table(clin$type, clin$nodal)
chisq.test(table(clin$type, clin$nodal))
# number of missings:
table(clin$type, is.na(clin$nodal))



# age:
by(clin$age, clin$type, median)
by(clin$age, clin$type, range)



sink()   # stop output into file



########################################################################################################

#### Now we may also use those samples that have been excluded from the finding cohort
#     as a validation cohort (validation-cohort-2 for Tao)


# SELECT VALIDATION COHORT OF BC with comparability metric >=5000

#  Select a subset of datasets with low comparability metric
# Select a subset of comparab by defining criteria:
compar.subset= subset(comparab, subset= sum.squ.nrm.diff >= 5000)
# vector of corresponding datasets:
datas.subset=compar.subset$dataset
# generate logical vector FALSE/TRUE for the complete dataset of 2007 BC:
subset.index.vector=(tdsn %in% datas.subset)
# query selected samples from transposed chipdata matrix:
tcdm.valid <- tcdm[subset.index.vector , ]
# query  corresponding transposed dataset vector:
tdsn.valid <- tdsn[subset.index.vector]
# query selected samples info from t.header:
t.header.valid <-  t.header[subset.index.vector,]

# remove the temporary subsetting variables:
rm(compar.subset, datas.subset, subset.index.vector)


########  Now again we select the sample with low/int-RS from the validation cohort:

library(genefu)

# Syntax der function oncotypedx: # oncotypedx(data, annot, do.mapping = FALSE, mapping, verbose = FALSE) 
# Note that for Affymetrix HGU datasets, the mapping is not necessary.
rs.valid <- oncotypedx(data= tcdm.valid, annot= U133A_Affy_Entrez, do.mapping=FALSE, verbose=TRUE)

# oncotypedx-score in t.header.find integrieren
t.header.valid <- cbind(t.header.valid,rs.valid$score)

# distribution of recurrence-scores from genefu-oncotypedx-function:
hist(rs.valid$score,breaks = 50)
quantile(rs.valid$score)
quantile(rs.valid$score, probs = seq(0, 1, 0.2))




# clinical parameters of VALIDATION cohort
table(t.header.valid$dataset)
table(t.header.valid$nodal)
table(t.header.valid$age50)
table(t.header.valid$t1vsrest)
table(t.header.valid$grade_1234)
table(t.header.valid$g1vs23)
table(t.header.valid$g12vs3)
table(t.header.valid$PGR_CO_MDAN_min0.0078)
table(t.header.valid$MSubTypeHugh)
table(t.header.valid$treat3class)
table(t.header.valid$Neoadj_Treat)
table(t.header.valid$pCR)

# Pragmatic selection of Low/Interm-Risk:  <80%-Quantile of rs.find$score


# SELECT LowInt-subset of validation cohort

# Select a subset of datasets with rs.find$score < 80%-Quantile 
# Select a subset of t.header.find by defining criteria:
t.header.valid.lowint <- subset(t.header.valid, subset= rs.valid$score < as.numeric(quantile(rs.valid$score, probs=0.8)))
# Select the subset of corresponding samples(rows) 
#   from transposed dataset of the valid cohort (tcdm.valid)
tcdm.valid.lowint <- tcdm.valid[rownames(tcdm.valid) %in% rownames(t.header.valid.lowint),]
dim(tcdm.valid.lowint) #  669 samples

############################



# A label "cohort" is defined which will be used in the script below to tag any output files

# AffyDataset, 669 samples VALIDATION cohort with low/int RS:
cohort <- "s669val_"
# clin info
clin <- t.header.valid.lowint
clin$Xcombi_id <- rownames(clin)
# mat from AffyDataset, 669 samples validation cohort with low/int RS, cols=Xcombi_id  rows=AffyIDs:
mat <- log(t(tcdm.valid.lowint))

###########################################################################################################
###
###
#######  NOW WE CAN AGAIN RUN THE WHOLE SCRIPT FROM TAO:
###
###



### TK-EDIT:   SELECT SAMPLE SUBGROUP   ###

### TK-EDIT:   Generate "type" as factor age>=55:old or age<=50:young
clin$type <- NA
clin$type[clin$age>=55] <- "old"
clin$type[clin$age<=50] <- "young"
clin$type <- as.factor(clin$type)

### TK-EDIT:   EXCLUDE SAMPLES with age between 51 and 54 years  (clin$type=NA)
clin <- clin[!is.na(clin$type),]
mat <- mat[,colnames(mat) %in% rownames(clin)]


#immune signatures genes
# Ori-script:
# signature=read.xlsx("./data/collected_signatures.xlsx",sheetIndex = "R")
# TK-REPLACEMENT:  ALTERED file and additional sheet "R_affy"  (+ using read_excel function):
#    !NOTE: "MKI67" as "single gene" was added here as "signature" since represented by 4 probesets
#           which do not cause error when using "apply" to calculate mean below!
library(readxl)
signature <- read_excel("data/tk_collected_signatures.xlsx", sheet = "R_affy")


# Ori-script:
# sig_lists=list()
# for(i in 1:dim(signature)[1]){
#   sig_lists[as.character(signature$Short[i])]=list(as.character(strsplit(as.character(signature$Genes[i]),split=",")[[1]]))
# }

# TK-REPLACEMENT:   Replace "Genes" with "Probesets"
sig_lists=list()
for(i in 1:dim(signature)[1]){
  sig_lists[as.character(signature$Short[i])]=list(as.character(strsplit(as.character(signature$Probesets[i]),split=",")[[1]]))
}

#estimate the expression signatures by taking the average expression of all member genes in a signature for each individual patient.
mat_signature=NULL
for(sig in names(sig_lists)){
  score=apply(mat[which(rownames(mat)%in%sig_lists[[sig]]),],2,mean)
  mat_signature=rbind(mat_signature,score)
}
rownames(mat_signature)=names(sig_lists)

### TK-EDIT:  ADD ESR1 as single gene to the mat_signature
ESR1 <- mat["205225_at",]
mat_signature <- rbind(mat_signature, ESR1)


# Ori-script:
#zscore normalized signatures
# mat_signature_zscore=as.data.frame(cbind(bcr_patient_barcode=colnames(mat_signature),apply(mat_signature,1,function(x)(x-mean(x))/sd(x))))
# TK-REPLACEMENT:
mat_signature_zscore=as.data.frame(cbind(Xcombi_id=colnames(mat_signature),apply(mat_signature,1,function(x)(x-mean(x))/sd(x))))


# Ori-script:
# finalMat=merge(clin,mat_signature_zscore,by="bcr_patient_barcode")
# TK-REPLACEMENT:   merge by Xcombi_id
finalMat=merge(clin, mat_signature_zscore, by="Xcombi_id")




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
########################################################################################################



# Cross tables needed for slide 3:

# set outfile:
fn=paste0("./out/",cohort,"CrossTables_for_slide-3.txt")
sink(fn,)


# numbers old/young
table(clin$type)


# Ethnicity
# 1=Caucasian/White, 2=AfricanAmerican/Black, 3=Non_Caucas_Non_AfricAmer, 4=Asian, 5=Hispanic,   9=n.a./mixed
table(clin$type, clin$RACE_code.2015_05_26)
chisq.test(table(clin$type, clin$RACE_code.2015_05_26))
# number of missings:
table(clin$type, is.na(clin$RACE_code.2015_05_26))




# T stage
# 1=T1, 2=T2/T3/T4
table(clin$type, clin$t1vsrest)
chisq.test(table(clin$type,clin$t1vsrest))
# number of missings:
table(clin$type, is.na(clin$t1vsrest))



# lymph node
# 0=LNN 1=N+
table(clin$type, clin$nodal)
chisq.test(table(clin$type, clin$nodal))
# number of missings:
table(clin$type, is.na(clin$nodal))



# age:
by(clin$age, clin$type, median)
by(clin$age, clin$type, range)



sink()   # stop output into file



############################################


