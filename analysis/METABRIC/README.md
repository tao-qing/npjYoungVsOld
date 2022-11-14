# ERposAge

## This directory [*METABRIC*](https://github.com/tkarn/ERposAge/blob/master/METABRIC/) contains data related to METABRIC analyses in the paper:


1. [*METABRIC_Analyses.Rmd*](https://github.com/tkarn/ERposAge/blob/master/METABRIC/METABRIC_Analyses.Rmd):  An R-Markdown file of the METABRIC analyses in the paper

2. [*METABRIC_Analyses.html*](https://github.com/tkarn/ERposAge/blob/master/METABRIC/METABRIC_Analyses.html):  The corresponding code together with all output from METABRIC_Analyses.Rmd

3. sudirectory [*data*](https://github.com/tkarn/ERposAge/blob/master/METABRIC/data/) with the following input files for the R-Markdown script:
 - *data_clinical.txt*
 - *signatures.txt*
 - *Common39Genes-probesets.txt*

4. sudirectory [*functions*](https://github.com/tkarn/ERposAge/blob/master/METABRIC/functions/) with the R-script:
 - *oncotypedx2_withoutBAG1.R*
 (an adapted oncotypedx function from the genefu package leaving out the gene BAG1, which is not available in METABRIC expression data)
 
5. [*METABRIC_Analyses.zip*](https://github.com/tkarn/ERposAge/blob/master/METABRIC/METABRIC_Analyses.zip): A zip package for downloading all files above (1.-4.): 

**NOTE:**

One additional file necessary for the METABRIC analysis is **NOT** available in this resource due to size:

- The log expression data and clinical data for 1980 breast cancer samples can be downloaded from [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric) or as a zip file from:
https://user.uni-frankfurt.de/~tkarn/data/METABRIC.zip

 
