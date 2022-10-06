#visulize DEG in three cohort
setwd("/Users/taoqing/Box Sync/ERPosExpression/")
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

require("biomaRt")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
affy_ensembl= c("affy_hg_u133_plus_2", "external_gene_name")
idconvert=getBM(attributes= affy_ensembl, mart= ensembl, values = "*", uniqueRows=T)

th = theme_bw()+ 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))+ 
    theme(legend.position='none',
          axis.text.x = element_text(colour="black", size=12,vjust = 0.95), 
          axis.text.y = element_text(colour="black", size=12,hjust = 0.95),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5,face="bold"),
          axis.title=element_text(size=12,face="bold"))+ 
    theme(strip.placement = "outside",
          plot.title = element_text(hjust = 0.5),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(angle = 0,size=10,face="italic"),
          panel.spacing = unit(0.1, "lines"),
          strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))

#deg 
degmat_tcga=fread("/Users/taoqing/Box Sync/ERPosExpression/out/TCGA_BRCA_ERPos_YoungVsOld_DEGs_20220530.csv",data.table = F,stringsAsFactors = F)
degmat_tcga$gene=sapply(degmat_tcga$V1,function(x)strsplit(x,split="\\|")[[1]][1])

p = ggplot(degmat_tcga, aes(y=-log10(FDR), x=log2(FC))) +
    geom_point(aes(colour=ifelse(degmat_tcga$FDR<0.05 & (degmat_tcga$FC >=1.5 | degmat_tcga$FC<=0.667),"significant","none")),shape=16) +
    geom_text_repel(aes(label=ifelse(-log10(FDR)>5 & abs(log2(FC))>0.585,gene,"")),size=3,segment.alpha =0.5,segment.colour ="grey") +
    ylab("-log10(FDR)") + 
    xlab("log2(FC)") + 
    theme_bw() +
    guides(color=guide_legend(override.aes=list(fill=NA)),linetype=guide_legend(override.aes=list(fill=NA))) + 
    theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14),
          axis.text.x = element_text(colour="black", size=14), 
          axis.text.y = element_text(colour="black", size=14)) + 
    labs(size = "-log10(P)") +
    scale_color_manual("FDR", values = c("significant" = "red","none" = "grey")) +
    th
p
fn = paste0("./out/TCGA_ERPos_Old_vs_Young_DEGs.pdf")
#ggsave(fn,w = 8, h = 5, useDingbat=F,limitsize = FALSE)
ggsave(fn,w = 5, h = 4, useDingbat=F,limitsize = FALSE)


#deg cohort A
degmat_A=fread("./ThomasDEG/DEG_Cohort_A.csv",data.table = F,stringsAsFactors = F)
degmat_A$Gene=idconvert$external_gene_name[pmatch(degmat_A$Probeset,idconvert$affy_hg_u133_plus_2)]
degmat_A=degmat_A[-which(is.na(degmat_A$Gene)),]


#deg cohort B
degmat_B=fread("./ThomasDEG/DEG_Cohort_B.csv",data.table = F,stringsAsFactors = F)
degmat_B$Gene=idconvert$external_gene_name[pmatch(degmat_B$Probeset,idconvert$affy_hg_u133_plus_2)]
degmat_B=degmat_B[-which(is.na(degmat_B$Gene)),]

#venn diagram
degmat_A$FC=as.numeric(as.matrix(degmat_A$FC))
degmat_B$FC=as.numeric(as.matrix(degmat_B$FC))

# Load library
library(VennDiagram)

# upregulated genes
tcga_up <- degmat_tcga$gene[(degmat_tcga$FC>=1.5) & degmat_tcga$FDR<0.05]
cohortA_up <- degmat_A$Gene[(degmat_A$FC>=1.5) & degmat_A$adj.P.Val<0.05]
cohortB_up <- degmat_B$Gene[(degmat_B$FC>=1.5) & degmat_B$adj.P.Val<0.05]

upintersec<-Reduce(intersect,list(tcga_up,cohortA_up,cohortB_up))

write.table(upintersec,"./out/Figure2_three_cohort_upregulated_genes_overlap.txt",quote=F,col.names = F,row.names = F)


# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(tcga_up, cohortA_up, cohortB_up),
  category.names = c("TCGA" , "Cohort-A" , "Cohort-B"),
  filename = 'up_regulated_gene_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


#down-regulated
tcga_down <- unique(degmat_tcga$gene[degmat_tcga$FC<=0.67 & degmat_tcga$FDR<0.05])
cohortA_down <- unique(degmat_A$Gene[degmat_A$FC<=0.67 & degmat_A$adj.P.Val<0.05])
cohortB_down <- unique(degmat_B$Gene[degmat_B$FC<=0.67 & degmat_B$adj.P.Val<0.05])

downintersec<-Reduce(intersect,list(tcga_down,cohortA_down,cohortB_down))

write.table(downintersec,"./out/Figure2_three_cohort_downregulated_genes_overlap.txt",quote=F,col.names = F,row.names = F)

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(tcga_down, cohortA_down, cohortB_down),
  category.names = c("TCGA" , "Cohort-A" , "Cohort-B"),
  filename = 'down_regulated_gene_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
