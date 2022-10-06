#install.packages("devtools")
#devtools::install_github("cansysbio/ConsensusTME")

setwd("/Users/taoqing/Box Sync/ERPosExpression/")
library(data.table)
library(dplyr)
library(ConsensusTME)

bulkExpMatrix <- as.matrix(read.csv("./data/TCGA_ER_breast_cancer.txt",h=T,r=1))

TMEscore=ConsensusTME::consensusTMEAnalysis(as.matrix(bulkExpMatrix), cancer = "BRCA", statMethod = "ssgsea")

TMEscore=t(TMEscore)

#clinical information
clin=fread("./data/TCGA_BRCA_Clinical.txt",data.table=F,h=T)

clin=clin %>%
  mutate(type=ifelse(Age<=50,"Young",ifelse(Age>=55,"old","other"))) %>%
  filter(type%in%c("Young","old") & HER2<15.17 & ERstatus=="Positive")


TMEscore = TMEscore %>%
            as.data.frame()  %>%
            mutate(bcr_patient_barcode=gsub("\\.","-",rownames(.))) %>%
            right_join(clin[,c("bcr_patient_barcode","ERstatus","Age","type")])


#save(TMEscore,file="./out/TCGA_BRCA_ConsensusTME_score.RData")


#Figure
library(ggplot2)
library(reshape2)

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

load("./out/TCGA_BRCA_ConsensusTME_score.RData")

TMEscore1 = melt(TMEscore,id=c("bcr_patient_barcode", "ERstatus", "Age",  "type")) %>%
            mutate(variable = as.character(variable)) %>%
            mutate( type = case_when(
                          type == "old" ~ "older",
                          type == "Young" ~ "younger"
                          ),
                    type = factor(type,levels=c("older","younger"))) %>%
           rename("cell_type" = "variable" )


stats_mat_freq <- NULL  
for(cell in unique(TMEscore1$cell_type)){

    these_data <- TMEscore1 %>%
      filter(cell_type == cell) %>%
      droplevels()
    
    new_stat <- cbind(tibble(cell_type = cell,
                             l_compare = levels(these_data$type)[1],
                             r_compare = levels(these_data$type)[2]),
                      broom::tidy(wilcox.test(value ~ type,
                                              data = these_data)))
    new_stat <- new_stat %>%
      mutate(significant = ifelse(p.value < 0.05,T,F)) %>%
      mutate(lab = if(p.value < 0.0001){
        "***"}else{
          if(p.value < 0.001){
            "**"
          }else{
            if(p.value < 0.01){"*"}else{"n.s."}  
          }
        },
        P = if(p.value < 0.001){
          "P < 0.0001"
          }else{
            #if(p.value >= 0.001){
             paste0("P = ",round(p.value,digits = 3))
            #}
        }
      ) %>%
      mutate(y_height = max(these_data$value))
    stats_mat_freq <- rbind(stats_mat_freq,new_stat)
}

#stats_mat_freq$cell_type = factor(stats_mat_freq$cell_type)


main_med <- TMEscore1 %>%
  group_by(cell_type, type) %>% 
  #filter(!is.na(neighborhood_category)) %>%
  summarize(med = median(value))


p <- 
  ggplot(TMEscore1,
         aes(x=type,y=value)) + 
  geom_violin() + 
  geom_segment(data = main_med, 
               aes(y=med,yend=med,x=as.numeric(type)-.2,
                   xend=as.numeric(type)+.2),
                   color="red") + 
  geom_jitter(alpha=.1,width = .1,size=.2) + 
  facet_wrap(~cell_type, scales = "free_y",ncol = 5) + 
  ylim(c(-0.4, 0.8)) +
  ggsignif::geom_signif(data=stats_mat_freq,
                        aes(xmin=l_compare,
                            xmax=r_compare,
                            annotations = P,
                            y_position=y_height
                        ),
                        manual=TRUE,vjust = -0.3) +
  th +
  #theme_classic() +
  #theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + 
  labs(y="Immune scores", x="")

ggsave("./out/ImmuneCell_Composition_ConsensusTME_oldvsyoung.pdf",w = 10, h = 10, useDingbat=F,limitsize = FALSE)




