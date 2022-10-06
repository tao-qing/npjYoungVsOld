setwd("/Users/qingtao/Box Sync/ERPosExpression/")
library(data.table)
library(dplyr)
#load data
#deg
#degmat=read.csv("./out/TCGA_BRCA_ERPos_YoungVsOld_DEGs.csv",stringsAsFactors = F)
#deg=degmat$Genes[degmat$DEGs!="None"]

clin=fread("./data/TCGA_BRCA_Clinical.txt",data.table=F,h=T)%>%mutate(type=ifelse(Age<=50,"young",ifelse(Age>=55,"old","other")))%>%filter(type%in%c("young","old")  & ( as.numeric(as.matrix(HER2))<15.17) & ERstatus=="Positive")%>%mutate(type=factor(type,levels=c("old","young")))

mat<-fread("./data/TCGA_ER_breast_cancer.txt",stringsAsFactors = F,data.table = F)
rownames(mat)=mat$V1
mat=mat[,-1]
mat=round(mat,digits = 2)


#DEG list 
#gmt_file = paste0("./GSEA/DEG.gmt")
#sink(gmt_file)
##for(i in names(deg)){
#  cat(i)
#  cat("\tna\t")
#  cat(paste(deg[[i]],collapse = "\t"))
#  cat("\n")
#}
#sink()

#Phenotype data
pheno_file = paste0("./GSEA/phenotype.cls")
sink(pheno_file)
cat("704 2 1\n")
cat("# Y O\n")
cat(c(rep("Y",table(clin$type)["young"]),rep("O",table(clin$type)["old"])))
cat("\n")
sink()

#expression data
finMat=mat[,c(clin$bcr_patient_barcode[clin$type=="young"],clin$bcr_patient_barcode[clin$type=="old"])]
finMat=cbind(NAME=rownames(finMat),DESCRIPTION="na",finMat)

write.table(finMat,"./data/TCGA_ER_breast_cancer_GSEA.txt")


library(data.table)
finMat=fread("TCGA_ER_breast_cancer_GSEA.txt",stringsAsFactors = F)
gct_file = paste0("TCGA_ER_breast_cancer_GSEA.gct")
sink(gct_file)
cat("#1.2\n")
cat(paste0(nrow(finMat), "\t", ncol(finMat), "\n"))
sink()

gct_out <- cbind(symbol = rownames(finMat), description = "na", finMat)
write.table(gct_out, gct_file, append = T, quote = F, row.names = F, sep = "\t")


library(fgsea)
load("./GSEA/MSigDB_gene_sets.RData")
degmat=read.csv("./out/TCGA_BRCA_ERPos_YoungVsOld_DEGs.csv",stringsAsFactors = F)
#degmat=degmat[which( (degmat$FC>=1.5 | degmat$FC <= 0.67) & degmat$FDR<0.05),]
#tmp=degmat[which( abs(log2(degmat$FC))>1 & degmat$FDR<0.05),]


pathways=gene_lists$MSigDBHallmark
ranks=log2(degmat$FC)
names(ranks)=degmat$Genes
ranks=sort(ranks, decreasing = T)


fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500, nperm=1000)

#The warning produced indicates that there are few genes that have the same fold change and so are ranked equally. fgsea with arbitrarily order determine which comes first in the ranked list. As long as this number is small it shouldn’t significantly effect the results. If the number is large something is suspicious about the fold change results.


library(ggplot2)
th= theme_bw()+ theme(legend.key.size =unit(.2, "cm"),legend.title = element_text(size=8) ,panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(colour="black", size=10,hjust = 0.95), axis.text.y = element_text(colour="black", size=10,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=14,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))

fgseaRes=as.data.frame(fgseaRes)

write.csv(fgseaRes[,1:7],"./GSEA/fgseaRes_out.csv",quote=F,row.names = F)

p = ggplot(fgseaRes, aes(y=-log10(fgseaRes$padj), x=fgseaRes$NES)) 
p = p + geom_point(aes(colour=ifelse(fgseaRes$padj<0.05,"significant","none")),shape=16)
p = p + geom_text_repel(aes(label=ifelse(-log10(fgseaRes$padj)>2,fgseaRes$pathway,"")),size=3,segment.alpha =0.5,segment.colour ="grey")
p = p + ylab("-log10(FDR)") + xlab("NES") + theme_bw() 
p = p + guides(color=guide_legend(override.aes=list(fill=NA)),linetype=guide_legend(override.aes=list(fill=NA)))+ theme(plot.title = element_text(size=14, face="bold"),axis.title = element_text(size=14), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14)) + labs(size = "-log10(P)") 
p = p + scale_color_manual("FDR<0.05", values = c("significant" = "red","none" = "grey"))
p=p+th
p
fn = paste0("./out/TCGA_ERPos_Old_vs_Young_DEGs.pdf")
ggsave(fn,w = 8, h = 5, useDingbat=F,limitsize = FALSE)
