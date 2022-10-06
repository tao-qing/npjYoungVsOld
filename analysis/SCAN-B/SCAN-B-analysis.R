setwd("/Users/taoqing/Box Sync/ERPosExpression/")

library(ggplot2)
library(ggsignif)
library(DescTools)
library(kSamples)
library(dplyr)
library(data.table)
library(tibble)
library(ggplot2)
library(xlsx)

cohort = "SCANB"

#####load data 
load("./data/SCANB_breastcancer_data.RData")
clin = clin %>%
      mutate(Age=as.numeric(as.matrix(Age))) %>%
      mutate(type = case_when(Age>=55~"old",
                          Age<=50~"young"))


table(clin$type)

#age 
tapply(clin$Age,clin$type,summary)


table(clin$Chemo_Treated,clin$type)
chisq.test(clin$Race,clin$type)


#compare signatures in young vs old using Mann -Whitney test, donot adjust by any covariates.
##this part will generate the table in the slide #7
#immune signatures genes
signature=read.xlsx("./data/collected_signatures.xlsx",sheetIndex = "R") 
signature=signature %>% dplyr::select(c("Short","Categories","Genes"))

signature = rbind(signature,
                  c("ESR1","ESR1","ESR1"),
                  c("MKI67","Preliferation","MKI67")
                  )


sig_lists=list()
for(i in 1:dim(signature)[1]){
  sig_lists[as.character(signature$Short[i])]=list(as.character(strsplit(as.character(signature$Genes[i]),split=",")[[1]]))
}

#estimate the expression signatures by taking the average expression of al member genes in a signature for each individual patient.
mat_signature=NULL
for(sig in names(sig_lists)){
  score=apply(mat[which(rownames(mat)%in%sig_lists[[sig]]),],2,mean)
  mat_signature=rbind(mat_signature,score)
}
rownames(mat_signature)=names(sig_lists)

#zscore normalized signatures
mat_signature_zscore=as.data.frame(cbind(Samples=colnames(mat_signature),apply(mat_signature,1,function(x)(x-mean(x))/sd(x))))

finalMat=clin %>% left_join(mat_signature_zscore, by="Samples")


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



ttMW=NULL
for(s in as.character(signature$Short)){
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

write.xlsx(ttMW,"./out/SCAN-B_Signature_YoungVsOld.xlsx")



#####3D Plot#######
library(plotly)
library(MASS)

tmp=finalMat%>%filter(type=="young")
mat3d=kde2d(as.numeric(tmp$ERS_Pos_Symmans), as.numeric(tmp$BCell), n = 1000)

fig <- plot_ly(x =mat3d$x, y = mat3d$y, z = mat3d$z) %>% add_surface()%>%layout(
  title = "",
  scene = list(
    xaxis = list(title = "ERS Pos Symmans"),
    yaxis = list(title = "B-cell signature"),
    zaxis = list(title = "Density of cases")
  ))
fig


youngMat=finalMat %>% 
         filter(type=="young") %>%
         mutate(ERS_Pos_Symmans=as.numeric(as.matrix(ERS_Pos_Symmans)),
                BCell=as.numeric(as.matrix(BCell)),
                group= case_when(ERS_Pos_Symmans >= median(ERS_Pos_Symmans) & BCell >= median(BCell) ~ "Immune-High ER-High",
                                 ERS_Pos_Symmans >= median(ERS_Pos_Symmans) & BCell < median(BCell) ~ "Immune-Low ER-High",
                                 ERS_Pos_Symmans < median(ERS_Pos_Symmans) & BCell >= median(BCell) ~ "Immune-High ER-Low",
                                 ERS_Pos_Symmans < median(ERS_Pos_Symmans) & BCell < median(BCell) ~ "Immune-Low ER-Low") 
                )
         

oldMat=finalMat %>%
       filter(type=="old") %>%
       mutate(ERS_Pos_Symmans=as.numeric(as.matrix(ERS_Pos_Symmans)),
              BCell=as.numeric(as.matrix(BCell)),
              group= case_when(ERS_Pos_Symmans >= median(ERS_Pos_Symmans) & BCell >= median(BCell) ~ "Immune-High ER-High",
                               ERS_Pos_Symmans >= median(ERS_Pos_Symmans) & BCell < median(BCell) ~ "Immune-Low ER-High",
                               ERS_Pos_Symmans < median(ERS_Pos_Symmans) & BCell >= median(BCell) ~ "Immune-High ER-Low",
                               ERS_Pos_Symmans < median(ERS_Pos_Symmans) & BCell < median(BCell) ~ "Immune-Low ER-Low") 
  )




library(ggplot2)
th = theme_bw() + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black")
) + theme(
  legend.position = 'top',
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


youngMatPlot = ggplot(data=youngMat, aes(x=ERS_Pos_Symmans, y=BCell,col=group)) + geom_point() + 
  geom_hline(yintercept=median(youngMat$BCell), linetype="dashed",col="grey") + 
  geom_vline(xintercept=median(youngMat$ERS_Pos_Symmans), linetype="dashed",col="grey") +
  guides(color=guide_legend(ncol=2))+
  th
ggsave("./out/SCAN-B-young-ERS_Pos_Symmans_vs_BCell.pdf",height=4,width=5)



oldMatPlot = ggplot(data=oldMat, aes(x=ERS_Pos_Symmans, y=BCell,col=group)) + geom_point() + 
geom_hline(yintercept=median(oldMat$BCell), linetype="dashed",col="grey") + 
geom_vline(xintercept=median(oldMat$ERS_Pos_Symmans), linetype="dashed",col="grey") +
guides(color=guide_legend(ncol=2))+
th
ggsave("./out/SCAN-B-old-ERS_Pos_Symmans_vs_BCell.pdf",height=4,width=5)


############Differentially expressed gene analysis###############
library(limma)
clin = clin[pmatch(colnames(mat),clin$Samples),]


# Construct model matrix including type
X <- model.matrix(~clin$type)
# Use limma to fit gene expression in mat to model matrix
fit <- lmFit(mat,X)
# Add Gene names
fit$genes <- row.names(fit)

# Empirical Bayes Statistics for Differential Expression
fit <- eBayes(fit)
#results <- decideTests(fit)
#summary(results)

deg <- topTable(fit, coef=2, adjust="fdr",n=30865)

expValue=apply(mat,1,function(x)tapply(x,clin$type,mean))

degmat<-as.data.frame(deg) %>%
  rename("Genes"="ID",
         "P value" = "P.Value",
         "FDR" = "adj.P.Val"
  ) %>%
  mutate('Mean Expression in Young'=expValue["young",Genes],
         'Mean Expression in Old'=expValue["old",Genes],
         GeneCard=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",Genes),
         DEGs=ifelse(logFC>=0.584 & FDR <0.05,"Significantly up-regulated in young",
                     ifelse(logFC<= (-0.584) & FDR <0.05,"Significantly down-regulated in young","None")),
         FC = 2^logFC) %>%
  select(c("Genes","Mean Expression in Young","Mean Expression in Old","FC","P value","FDR","DEGs","GeneCard"))

  write.csv(degmat,file="./out/SCANB_ERPos_YoungVsOld_DEGs_20220630.csv")


###heatmap for DEGs####
require("ComplexHeatmap")
clin = clin[order(clin$type,decreasing=T),]

plotMat=mat %>% 
        filter(rownames(.) %in% degmat$Genes[which(degmat$DEGs!="None")])

plotMat = plotMat[,pmatch(clin$Samples,colnames(plotMat))]
        

library(circlize)
#col_fun = colorRamp2(c(0, , 15), c("green", "white", "red"))
allheatmap=Heatmap(as.matrix(plotMat),name = "Exp",
                   top_annotation = HeatmapAnnotation(Group=clin$type,col=list(Group=c("old"="lightgreen","young"="red"))),
                   show_row_names=FALSE,
                   cluster_columns=FALSE,
                   show_column_names = FALSE#,
                   #col=col_fun
)

pdf("./out/SCANB_allDEG_heatmap.pdf",height=4,width=5)
allheatmap
dev.off()


#survival analysis
require('ggfortify')
library(survminer)
library(survival)

source("./TCGASurvival/run_survival_function.R")
source("./TCGASurvival/samstein_functions.R")
source("./TCGASurvival/createOncoMatrix.R")
source("./TCGASurvival/ggkmTable.R")


plotReadableKM = function(data=NULL, predictor=NULL,surv_time,status,legend_lbls=NULL, ptitle=NULL){
  model <- formula(paste0("Surv(data$",surv_time,", data$",status,") ~", predictor))
  fit <- surv_fit(model, data)
  p <- ggsurvplot(fit, data, pval=TRUE, pval.method = TRUE, risk.table = TRUE,
                  tables.theme = theme_cleantable(),
                  title = ptitle,
                  xlab = "Day",
                  legend.title = "",
                  legend.labs = legend_lbls,
                  risk.table.fontsize = 8,
                  pval.size = 7,
                  font.x = c(22),
                  font.y = c(22),
                  font.tickslab = c(20)
  )
  
  p$plot <- p$plot + theme(legend.title = element_text(size=10), 
                           legend.text = element_text(size=10),
                           plot.title = element_text(size=15)
  ) + 
    #guides(color = guide_legend(override.aes=list(size=3))) +
    guides(color=guide_legend(ncol=2)) 
 # p$table <- p$table + theme(plot.title=element_text(size=22),
  #                            axis.text.y = element_text(size=22))
  
  return(p)
}
  
          
p1=plotReadableKM(data=oldMat %>% 
                       mutate(OS=as.numeric(OS),
                              OS_Event=as.numeric(OS_Event)) %>%
                       filter(group %in% c("Immune-High ER-Low","Immune-Low ER-High")) %>%
                       filter(Chemo_Treated == 1), 
                 predictor="group", 
                 surv_time = "OS",
                 status="OS_Event",
                 legend_lbls=NULL, 
                 ptitle="Old patients with Chemotherapy") 

ggsave("./out/Survival_Old_patients_with_Chemotherapy.png",plot=p1$plot,w = 7, h = 7,dpi=100,limitsize = FALSE)#, useDingbat=F,limitsize = FALSE


p2=plotReadableKM(data=oldMat %>% 
                 mutate(OS=as.numeric(OS),
                        OS_Event=as.numeric(OS_Event)) %>%
                 filter(group %in% c("Immune-High ER-Low","Immune-Low ER-High")) %>%
               filter(Chemo_Treated == 0)
               , 
               predictor="group", 
               surv_time = "OS",
               status="OS_Event", 
               legend_lbls=NULL, 
               ptitle="Old patients without Chemotherapy")
ggsave("./out/Survival_Old_patients_without_Chemotherapy.pdf",plot=p2$plot,w = 7, h = 7, useDingbat=F,limitsize = FALSE)


p3=plotReadableKM(data=oldMat %>% 
                 mutate(OS=as.numeric(OS),
                        OS_Event=as.numeric(OS_Event)) %>%
                 filter(group %in% c("Immune-High ER-Low","Immune-Low ER-High")) 
               , 
               predictor="group", 
               surv_time = "OS",
               status="OS_Event", 
               legend_lbls=NULL, 
               ptitle="All old patients")
ggsave("./out/Survival_all_Old_patients_Chemotherapy.pdf",plot=p3$plot,w = 7, h = 7, useDingbat=F,limitsize = FALSE)


p4=plotReadableKM(data=youngMat %>% 
                    mutate(OS=as.numeric(OS),
                           OS_Event=as.numeric(OS_Event)) %>%
                    filter(group %in% c("Immune-High ER-Low","Immune-Low ER-High")) %>%
                    filter(Chemo_Treated == 1)
                  , 
                  predictor="group", 
                  surv_time = "OS",
                  status="OS_Event", 
                  legend_lbls=NULL, 
                  ptitle="Young patients with Chemotherapy")
ggsave("./out/Survival_Young_patients_with_Chemotherapy.pdf",plot=p4$plot,w = 7, h = 7, useDingbat=F,limitsize = FALSE)


p5=plotReadableKM(data=youngMat %>% 
                    mutate(OS=as.numeric(OS),
                           OS_Event=as.numeric(OS_Event)) %>%
                    filter(group %in% c("Immune-High ER-Low","Immune-Low ER-High")) %>%
                    filter(Chemo_Treated == 0)
                  , 
                  predictor="group", 
                  surv_time = "OS",
                  status="OS_Event", 
                  legend_lbls=NULL, 
                  ptitle="Young patients without Chemotherapy")
ggsave("./out/Survival_Young_patients_without_Chemotherapy.pdf",plot=p5$plot,w = 7, h = 7, useDingbat=F,limitsize = FALSE)


p6=plotReadableKM(data=youngMat %>% 
                    mutate(OS=as.numeric(OS),
                           OS_Event=as.numeric(OS_Event)) %>%
                    filter(group %in% c("Immune-High ER-Low","Immune-Low ER-High")) 
                  , 
                  predictor="group", 
                  surv_time = "OS",
                  status="OS_Event", 
                  legend_lbls=NULL, 
                  ptitle="All Young patients")
ggsave("./out/Survival_all_Young_patients_Chemotherapy.pdf",plot=p6$plot,w = 7, h = 7, useDingbat=F,limitsize = FALSE)

