setwd("/Users/zhanglei/Box Sync/ERPosExpression/")
library(maftools)
library(data.table)
library(dplyr)

maf = fread("TCGA_BRCA_maf_maftools.txt",stringsAsFactors = F,data.table = F)
GATA3maf = maf %>% 
           filter(Hugo_Symbol=="GATA3")

clin = fread("TCGA_BRCA_clin_maftools.txt",stringsAsFactors = F,data.table = F) %>%
       mutate(GATA3mut = ifelse(Tumor_Sample_Barcode%in%GATA3maf$Tumor_Sample_Barcode,1,0),
              Subtype = ifelse(Subtype_mRNA == "LumB","LumB","Other"),
              Subtype = factor(Subtype,c("Other","LumB")),
              AgeGroup = type
       )


logit <- glm(GATA3mut ~ AgeGroup + Subtype + AgeGroup * Subtype , data = clin, family = "binomial")
summary(logit)



Call:
  glm(formula = mut ~ type + subtype + type * subtype, family = "binomial", 
      data = clin)

Deviance Residuals: 
  Min       1Q   Median       3Q      Max  
-0.7422  -0.6303  -0.4379  -0.4379   2.1873  

Coefficients:
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)            -2.2963     0.1855 -12.382  < 2e-16 ***
  typeYoung               1.1477     0.2611   4.396  1.1e-05 ***
  subtypeLumB             0.7812     0.3088   2.529   0.0114 *  
  typeYoung:subtypeLumB  -0.8693     0.5220  -1.665   0.0959 .  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 562.42  on 662  degrees of freedom
Residual deviance: 540.27  on 659  degrees of freedom
AIC: 548.27

Number of Fisher Scoring iterations: 5








