setwd("/Users/qingtao/Box Sync/ERPosExpression/")
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(maftools)

#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')

laml = read.maf(maf = laml.maf,
                clinicalData = laml.clin,
                verbose = FALSE)
oncoplot(maf = laml, draw_titv = TRUE)


#GISTIC results LAML
all.lesions = system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes = system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes = system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis = system.file("extdata", "scores.gistic", package = "maftools")

#Read GISTIC results along with MAF
laml.plus.gistic = read.maf(
  maf = laml.maf,
  gisticAllLesionsFile = all.lesions,
  gisticAmpGenesFile = amp.genes,
  gisticDelGenesFile = del.genes,
  gisticScoresFile = scores.gis,
  isTCGA = TRUE,
  verbose = FALSE, 
  clinicalData = laml.clin
)

library(TCGAWorkflowData)
data("elmerExample")
data("GBMnocnvhg19")
data("GBMIllumina_HiSeq")
data("LGGIllumina_HiSeq")
data("mafMutect2LGGGBM")
data("markersMatrix")
data("histoneMarks")
data("biogrid")
data("genes_GR")
library(GenomicRanges)
library(TCGAbiolinks)
##############################
## Recurrent CNV annotation ## 
##############################
# Get gene information from GENCODE using biomart
genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg19") 
genes <- genes[genes$chromosome_name %in% c(1:22,"X","Y"),]
genes[genes$chromosome_name == "X", "chromosome_name"] <- 23
genes[genes$chromosome_name == "Y", "chromosome_name"] <- 24
genes$chromosome_name <- sapply(genes$chromosome_name,as.integer)
genes <- genes[order(genes$start_position),]
genes <- genes[order(genes$chromosome_name),]
genes <- genes[,c("external_gene_name", "chromosome_name", "start_position","end_position")]
colnames(genes) <- c("GeneSymbol","Chr","Start","End")
genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)

library(RTCGAToolbox)
# Download GISTIC results
lastAnalyseDate <- getFirehoseAnalyzeDates(1)
gistic <- getFirehoseData("GBM",gistic2_Date = lastAnalyseDate, GISTIC = TRUE)

# get GISTIC results
gistic.allbygene <- getData(gistic,type = "GISTIC", platform = "AllByGene")
gistic.thresholedbygene <- getData(gistic,type = "GISTIC", platform = "ThresholdedByGene")

library(TCGAbiolinks)
query.gbm.nocnv <- GDCquery(project = "TCGA-GBM",
                            data.category = "Copy number variation",
                            legacy = TRUE,
                            file.type = "nocnv_hg19.seg",
                            sample.type = c("Primary solid Tumor"))
# to reduce time we will select only 20 samples
query.gbm.nocnv$results[[1]] <- query.gbm.nocnv$results[[1]][1:20,]


GDCdownload(query.gbm.nocnv, chunks.per.download = 100)

gbm.nocnv <- GDCprepare(query.gbm.nocnv, save = TRUE, save.filename = "GBMnocnvhg19.rda")



library(TCGAbiolinks)
# Samples: primary solid tumor w/ DNA methylation and gene expression
matched_met_exp <- function(project, n = NULL){
  # get primary solid tumor samples: DNA methylation
  message("Download DNA methylation information")
  met450k <- GDCquery(project = project,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450",
                      legacy = TRUE, 
                      sample.type = c("Primary solid Tumor"))
  met450k.tp <-  met450k$results[[1]]$cases
  
  # get primary solid tumor samples: RNAseq
  message("Download gene expression information")
  exp <- GDCquery(project = project,
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "normalized_results", 
                  sample.type = c("Primary solid Tumor"),
                  legacy = TRUE)
  exp.tp <- exp$results[[1]]$cases
  # Get patients with samples in both platforms
  patients <- unique(substr(exp.tp,1,15)[substr(exp.tp,1,12) %in% substr(met450k.tp,1,12)] )
  if(!is.null(n)) patients <- patients[1:n] # get only n samples
  return(patients)
}
lgg.samples <- matched_met_exp("TCGA-LGG", n = 10)
gbm.samples <- matched_met_exp("TCGA-GBM", n = 10)
samples <- c(lgg.samples,gbm.samples)

#-----------------------------------
# 1 - Methylation
# ----------------------------------
query.met <- GDCquery(project = c("TCGA-LGG","TCGA-GBM"),
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450",
                      legacy = TRUE, 
                      barcode = samples)
GDCdownload(query.met)
met <- GDCprepare(query.met, save = FALSE)
met <- subset(met,subset = as.character(GenomicRanges::seqnames(met)) %in% c("chr9"))

#-----------------------------------
# 2 - Expression
# ----------------------------------
query.exp <- GDCquery(project = c("TCGA-LGG","TCGA-GBM"),
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type  = "normalized_results", 
                      legacy = TRUE, 
                      barcode =  samples)
GDCdownload(query.exp)
exp <- GDCprepare(query.exp, save = FALSE)
save(exp, met, gbm.samples, lgg.samples, file = "elmerExample.rda", compress = "xz")


library(TCGAbiolinks)
LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")
GBMmut <- GDCquery_Maf(tumor = "GBM", pipelines = "mutect2")
mut <- plyr::rbind.fill(LGGmut, GBMmut)
save(mut, LGGmut, GBMmut,file = "mafMutect2LGGGBM.rda")

file <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/CNV.hg19.bypos.111213.txt"
if(!file.exists(basename(file))) downloader::download(file, basename(file))
commonCNV <- readr::read_tsv(basename(file), progress = FALSE)
commonCNV <- as.data.frame(commonCNV)




#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)





query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)


query.maf.hg19 <- GDCquery(project = "TCGA-BRCA", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           legacy = TRUE)

#somatic mutation

#copy number variation

#methylation


