setwd("/Users/taoqing/Box Sync/ERPosExpression/")

library(ggplot2)
library(ggsignif)
library(DescTools)
library(dplyr)
library(data.table)

#clinical data
covariates=fread("/Users/taoqing/Box Sync/GermlineSomatic (tao.qing@yale.edu)/Huang_lab_data/TCGA_PanCanAtlas_2018/covariates/TCGA_ancestry_PC.txt",data.table = F)

clin = fread("/Users/taoqing/Box Sync/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/clinical_PANCAN_patient_with_followup.tsv",data.table = F) %>%
       filter(acronym=="BRCA") %>%
       select("bcr_patient_barcode",
              "breast_carcinoma_estrogen_receptor_status",
              "age_at_initial_pathologic_diagnosis",
              "number_of_lymphnodes_positive_by_he",
              "pathologic_stage",
              "menopause_status",
              "her2_immunohistochemistry_level_result",
              "lab_proc_her2_neu_immunohistochemistry_receptor_status",
              "er_level_cell_percentage_category") %>%
       mutate(washu_assigned_ethnicity=covariates$washu_assigned_ethnicity[pmatch(bcr_patient_barcode,covariates$bcr_patient_barcode)]) %>%
       mutate(lymphnodes_status = case_when( number_of_lymphnodes_positive_by_he == "[Not Available]" ~ "Unknown",
                                             number_of_lymphnodes_positive_by_he == "0" ~ "negative",
                                             number_of_lymphnodes_positive_by_he > 0 ~ "positive")
       ) %>%
      mutate(menopause_status_fin = case_when( menopause_status %in% c("[Not Available]","[Not Evaluated]","[Unknown]") ~ "Unknown",
                                               menopause_status == "Indeterminate (neither Pre or Postmenopausal)" ~ "Indeterminate",
                                               menopause_status %in% c("Peri (6-12 months since last menstrual period)",
                                                                       "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)") ~ "Pre",
                                               menopause_status == "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)" ~ "Post"
                                               )
      ) %>%
       mutate(washu_assigned_ethnicity = case_when(washu_assigned_ethnicity%in%c("African","Asian","European") ~ washu_assigned_ethnicity,
                                                   !washu_assigned_ethnicity%in%c("African","Asian","European") ~ "Unknown")
       ) %>%
       mutate(pathologic_stage=ifelse(pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"),"Stage I & II",
                                             ifelse(pathologic_stage%in%c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC","Stage IV","Stage X"),"Stage III IV X","Unknown"))) %>%
       rename(ERstatus=breast_carcinoma_estrogen_receptor_status,
              HER2status=lab_proc_her2_neu_immunohistochemistry_receptor_status,
              Age=age_at_initial_pathologic_diagnosis,
              Race=washu_assigned_ethnicity)


tmp = clin[duplicated(clin$bcr_patient_barcode),]

#expression 
#mat=readRDS("/Users/qingtao/Analysis/R/YaleBRCA/TCGAData/Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.UQnormalized/PanCancer_RNA_seq_UQ_norm.rds")
mat = readRDS("./data/PanCancer_RNA_seq_UQ_norm.rds") %>%
      select(which(substr(names(.),14,16)%in%c("01A","11A"))) %>%
      select(which(!duplicated(substr(names(.),1,12)))) %>%
      `colnames<-`(gsub("\\.","-",substr(colnames(.),1,12))) %>% 
      filter(!grepl("\\?",rownames(.))) %>%
      select(which(names(.)%in%clin$bcr_patient_barcode))


anno=sapply(rownames(mat),function(x)as.matrix(strsplit(x,split="\\|")[[1]][c(1,2)])) %>%
     t() %>%
     as.data.frame() %>%
     rename_with(.,~c("symbol","EntrezGene.ID"))  %>%
     `rownames<-` (NULL)

dupindex=which(duplicated(anno$symbol))
mat=mat[-dupindex,]
anno=anno[-dupindex,]
rownames(mat)=anno$symbol

#fwrite(mat,"./out/TCGA_breast_cancer.txt",row.names = T)

#TCGA subtype
library(TCGAbiolinks)
subtypes <- PanCancerAtlas_subtypes() %>%
            as.data.frame() %>%
            filter(cancer.type == "BRCA") %>%
            mutate(bcr_patient_barcode=substr(pan.samplesID,1,12)) %>%
            filter(substr(pan.samplesID,14,16)%in%c("01A","11A")) %>%
            filter(!duplicated(bcr_patient_barcode))
            

#oncotype score
load("./out/TCGA_breast_cancer_OncotypeDX_score.RData")

overlaped = intersect(colnames(mat),clin$bcr_patient_barcode)

#final clinical matrix
clin = clin %>%
  filter(bcr_patient_barcode%in%overlaped) %>%
  mutate(ESR1=as.numeric(mat["ESR1",bcr_patient_barcode]),
         HER2=as.numeric(mat["ERBB2",bcr_patient_barcode]),
         OncotypeDx_score=rs_mb$score[bcr_patient_barcode]
  ) %>% 
  left_join(subtypes[,c("bcr_patient_barcode","cancer.type","Subtype_mRNA")],by="bcr_patient_barcode") %>%
  mutate(type=ifelse(Age<=50,"Young",ifelse(Age>=55,"Old","other"))) %>%
  filter(type%in%c("Young","Old") & HER2<15.17 & ERstatus=="Positive") %>%
  slice_min(OncotypeDx_score, prop = 0.8)

mat = mat %>% 
      select(clin$bcr_patient_barcode)


save(mat,clin,anno,file="./data/TCGA_BRCA.RData")


library(ggplot2)
th= theme_bw()+ 
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

#ERS1 expression versus ER cell composition
clin_ESR1=clin%>%
          mutate(type=ifelse(Age<=50 & ERstatus=="Positive","young_ER_Positive",
                             ifelse(Age>=55 & ERstatus=="Positive","old_ER_Positive",
                                    ifelse(ERstatus=="Negative","ER_negative","other")))) %>%
          filter(type%in%c("young_ER_Positive","old_ER_Positive","ER_negative")) %>%
          mutate(type=factor(type,levels=c("young_ER_Positive","old_ER_Positive","ER_negative"))) %>%
          filter(as.numeric(as.matrix(HER2))<15.17 & type!="ER_negative")


med <- clin_ESR1 %>%
  group_by(type) %>% 
  filter(!is.na(type)) %>%
  summarize(med = median(as.numeric(as.matrix(ESR1))))


p=ggplot(clin_ESR1,
         aes(x=type,y=as.numeric(as.matrix(ESR1)))) + 
  geom_violin() + 
  geom_jitter(width = .1,size=.2)+th+xlab("")+ylab("ESR1 expression")+#+ggtitle(paste0("p value=",pvalue,"\nFDR=",FDR))
  geom_segment(data = med, 
               aes(y=med,yend=med,x=as.numeric(type)-.2,
                   xend=as.numeric(type)+.2),color="red")


fn=paste0("./out/ESR1_expression_in_each_group.pdf")
ggsave(fn,w = 5, h = 3, useDingbat=F,limitsize = FALSE)


#ERS1 by ER percentile
p=ggplot(clin_ESR1 %>% 
           filter(grepl("%",er_level_cell_percentage_category) & !is.na(ESR1))%>%mutate(er_level_cell_percentage_category=gsub("<10%","1-9%",er_level_cell_percentage_category)),
         aes(x=er_level_cell_percentage_category,y=as.numeric(as.matrix(ESR1)))) + 
  geom_boxplot(fill="blue")+th+xlab("ER Percentage Categories")+ylab("ESR1 expression")

fn=paste0("./out/ESR1_expression_in_ER_Percent_Cate.pdf")
ggsave(fn,w = 8, h = 3, useDingbat=F,limitsize = FALSE)


library(DescTools)
tmp=clin_ESR1 %>% filter(grepl("%",er_level_cell_percentage_category) & !is.na(ESR1))

tmp$er_level_cell_percentage_category=factor(tmp$er_level_cell_percentage_category,levels=names(table(tmp$er_level_cell_percentage_category)),ordered=T)

JTres <- JonckheereTerpstraTest(tmp$ESR1,tmp$er_level_cell_percentage_category,alternative="increasing", nperm=10000)
ccJTres <- cor.test(tmp$ESR1,as.numeric(tmp$er_level_cell_percentage_category),method="k", alternative = "greater")

#tau  = 0.2745422  ,p-value = 1.508e-11


#age by ER percentage
p = ggplot(clin_ESR1 %>% 
    filter(grepl("%",er_level_cell_percentage_category) & !is.na(ESR1)) %>%
    mutate(er_level_cell_percentage_category=gsub("<10%","1-9%",er_level_cell_percentage_category)),
          aes(x=er_level_cell_percentage_category,
              y=as.numeric(as.matrix(Age)))) + 
  geom_boxplot(fill="blue") +
  th +
  xlab("ER Percentage Categories") +
  ylab("Age")

fn=paste0("./out/Age_distribution_in_ER_Percent_Cate.pdf")
ggsave(fn,w = 8, h = 3, useDingbat=F,limitsize = FALSE)

tmp=clin_ESR1 %>% filter(grepl("%",er_level_cell_percentage_category) & !is.na(Age))
tmp$er_level_cell_percentage_category=factor(tmp$er_level_cell_percentage_category,levels=names(table(tmp$er_level_cell_percentage_category)),ordered=T)
tmp$Age=as.numeric(as.matrix(tmp$Age))
JTres <- JonckheereTerpstraTest(tmp$Age,tmp$er_level_cell_percentage_category,alternative="increasing", nperm=10000)
ccJTres <- cor.test(tmp$Age,as.numeric(tmp$er_level_cell_percentage_category),method="k", alternative = "greater")
#tau 0.03626571, p-value = 0.1927


p = ggplot(clin_ESR1) + aes(x = clin_ESR1$type, fill = factor(er_level_cell_percentage_category)) + 
    geom_bar(position = "fill") +
    theme_bw() +
    ylab("Proportions") + 
    xlab("") +
    th
    #theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5,hjust = 0.95), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0,size=12,face="bold"),axis.title=element_text(size=12,face="bold"),panel.border = element_blank(),axis.line= element_line(color='white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p = p + scale_fill_manual("Grade",values = c("1" = "green",  "2" = "gray", "3"="red"))
p
ggsave("Grade_distribution_in_PDL1Cate.pdf",plot=p,w = 3, h = 3, useDingbat=F,limitsize = FALSE)



clin_ESR1$Age=as.numeric(as.matrix(clin$Age))
clin_ESR1$ESR1=as.numeric(as.matrix(clin$ESR1))

cor(clin_ESR1$ESR1,clin_ESR1$Age,method = "spearman")
cor.test(clin_ESR1$ESR1,clin_ESR1$Age,method = "spearman")

p = ggplot(clin_ESR1, aes(y=Age, x=ESR1)) +
    geom_point(shape=16) +
    ylab("Age of Diagnosis") + xlab("ESR1 expression") + 
    theme_bw() +
    guides(color=guide_legend(override.aes=list(fill=NA)),
           linetype=guide_legend(override.aes=list(fill=NA))) +
    geom_smooth(method='lm')+
    th
p
fn = paste0("./out/TCGA_ESR1_Age_correlation.pdf")
ggsave(fn,w = 5, h = 5, useDingbat=F,limitsize = FALSE)


