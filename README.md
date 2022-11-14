# ERposAge

## Supplementary Resource to:  

Qing T, Karn T et al. Molecular differences between younger versus older ER-positive and HER2-negative breast cancers.
npj Breast Cancer. 2022 Nov 7;8(1):119. doi: 10.1038/s41523-022-00492-0

*PMID 36344517*

[https://www.ncbi.nlm.nih.gov/pubmed/36344517/](https://pubmed.ncbi.nlm.nih.gov/36344517/)

************************************************************

## This resource contains the following data regarding the analyses described in the paper:

## Directory [*Microarray*](https://github.com/tkarn/ERposAge/blob/master/Microarray/) contains data related to Affymetrix microarray analyses in the paper:


1. [*Affymetrix_Analyses.Rmd*](https://github.com/tkarn/ERposAge/blob/master/Microarray/Affymetrix_Analyses.Rmd):  An R-Markdown file of the Affymetrix analyses in the paper

2. [*Affymetrix_Analyses.html*](https://github.com/tkarn/ERposAge/blob/master/Microarray/Affymetrix_Analyses.html):  The corresponding code together with all output from Affymetrix_Analyses.Rmd

3. Sudirectory [*data*](https://github.com/tkarn/ERposAge/blob/master/Microarray/data/) with the following input files for the R-Markdown script:
 - *n2007_SampleInfo.txt*
 - *signatures.txt*
 - *Common39Genes-probesets.txt*

4. [*Affymetrix_Analyses.zip*](https://github.com/tkarn/ERposAge/blob/master/Microarray/Affymetrix_Analyses.zip): A zip package for downloading all files above (1.-3.): 

**NOTE:**

One additional file necessary for the microarray analysis is **NOT** available in this resource due to size:

- The file  *n2007_MAS5.txt*  contains all Affymetrix microarray data and can be downloaded as a 187 MB zip file from:
https://user.uni-frankfurt.de/~tkarn/data/n2007_MAS5.zip

---

## Directory [*METABRIC*](https://github.com/tkarn/ERposAge/blob/master/METABRIC/) contains data related to METABRIC analyses in the paper:


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

 
