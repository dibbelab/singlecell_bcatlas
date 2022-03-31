# Single-cell Breast Cancer Atlas
[![DOI](https://zenodo.org/badge/362182268.svg)](https://zenodo.org/badge/latestdoi/362182268)

On-line single cell <b>BC Atlas available</b> at [http://bcatlas.tigem.it](http://bcatlas.tigem.it)

# A single-cell analysis of breast cancer cell lines to study tumour heterogeneity and drug response
#### G. Gambardella<sup>1,2+</sup>, G. Viscido<sup>1,2,+</sup>, B. Tumaini<sup>1</sup>, A. Isacchi<sup>3</sup>, R. Bosotti<sup>3</sup>, D. di Bernardo<sup>1,2,++</sup>

<sup>1</sup>Telethon Institute of Genetics and Medicine, Armenise/Harvard Laboratory of Integrative Genomics, Pozzuoli, Italy.  
<sup>2</sup>University of Naples Federico II, Department of Chemical, Materials and Industrial Engineering, Naples, Italy.  
<sup>3</sup>NMSsrl, Nerviano Medical Sciences, 20014, Nerviano, Milan, Italy.

## Abstract
Cancer cells within a tumour have heterogeneous phenotypes and exhibit dynamic plasticity. How to evaluate such heterogeneity and its impact on outcome and drug response is still unclear. Here, we transcriptionally profile 35,276 individual cells from 32 breast cancer cell lines to yield a single cell atlas. We find high degree of heterogeneity in the expression of biomarkers. We then train a deconvolution algorithm on the atlas to determine cell line composition from bulk gene expression profiles of tumour biopsies, thus enabling cell line-based patient stratification. Finally, we link results from large-scale in vitro drug screening in cell lines to the single cell data to computationally predict drug responses starting from single-cell profiles. We find that transcriptional heterogeneity enables cells with differential drug sensitivity to co-exist in the same population. Our work provides a framework to determine tumour heterogeneity in terms of cell line composition and drug response.

The full article [(Gambardella, Viscido et al. 2022)](https://www.nature.com/articles/s41467-022-29358-6)

## Figures

![alt text](https://github.com/dibbelab/singlecell_bcatlas/blob/main/figures/Figure_01.png?raw=true)

<b>Figure 1 – The Breast Cancer Single Cell Atlas.</b> (A) Representation of single-cell expression profiles of 35,276 cells from 32 cell lines color-coded according to cancer subtype (LA=Luminal A, LB=Luminal B, H=Her2-enirched, TNA = Triple Negative type A, TNB = Triple Negative type B). (B) Expression levels of the indicated genes in the atlas, with red indicating expression, together with their (C) distribution within the cell lines, shown as a violin plot. (D) Dotplot of literature-based biomarker genes along the columns for each of the 32 sequenced cell lines along the rows. Biomarker genes are grouped by type (Basal Epith. = Basal Epithelial, Luminal Epith. = Luminal Epithelial, L.P. = Luminal Progenitor, EMT = Epithelial to Mesenchymal Transition). (E) Graphical representation of 35,276 cells color-coded according to their cluster of origin. Clusters are numbered from 1 to 22. (F) For the indicated cluster, the corresponding pie-chart represents the cluster composition in terms of cell lines. Cell lines in the same pie-chart are distinguished by colour. Only the top 10 most heterogenous clusters are shown. Gray colour represents the cell-lines in the cluster contributing with less then 5% of total cells in the cluster. While the other colours represent a specific cell line. For example, Cluster 2 is the most heterogeneous cluster mainly composed by 3 cell-line while cluster 19 is the most homogeneous since in its mainly composed by the cells coming from one cell-line. (G) Expression levels in the atlas of the five luminal biomarkers identified as the most differentially expressed in each of the five luminal clusters (1, 2, 6, 8 and 18). (H) Expression of 22 atlas-derived biomarkers in the biopsies of 937 breast cancer patient from TCGA. (I) Accuracy in classifying tumour subtype for 937 patients from TCGA by using either the PAM50 gene signature or the 22 cluster-derived  biomarker genes (scCCL) alone or augmented with HER2 gene (scCCL + HER2).
<hr/>

![alt text](https://github.com/dibbelab/singlecell_bcatlas/blob/main/figures/Figure_02.png?raw=true)
<b>Figure 2 –Automatic classification of patients’ tumour cells</b> (A) Cancer cells from triple negative breast cancer (TNBC) biopsies of five patients were embedded in the BC atlas by means of the mapping algorithm in order to predict their tumour subtype. (B) For each patient, the pie chart shows  cell line composition obtained by mapping patient’s cells onto the atlas. (C) Tissue-slide of an oestrogen positive breast tumour biopsy sequenced by means of the 10x Genomics Visium spatial transcriptomics (top-left) and the position of the mapped tissue tiles onto the atlas (top-left). Tiles are colour-coded according to the cell line (bottom-left) and to tumour subtype (bottom-right) as predicted by the mapping algorithm. (D) Cell line composition of 937 BC patients from the TGCA database as estimated by the Bisque algorithm from their bulk RNA-seq data. For ease of interpretation, in the heatmap patients are clustered according to their cell line composition. The bottom row reports the annotated cancer subtype in TGCA. (E) Predicted cell-line composition by the Bisque algorithm for four representative patients. (F) The distribution of the 937 BC patients across the 32 cell lines. For each cell line, the stacked  bars report the percentage of patients of a given cancer subtype assigned by Bisque to that cell line. Since each patient is usually predicted by Bisque to be composed by multiple cell lines, the patient is associated to the cell line making up  the largest fraction  of the patient’s cell line composition. (G) Performance of Bisque in classifying the tumour subtype of the 937 BC patients in TGCA from bulk gene expression profiles. Since each patient is usually predicted by Bisque to be composed by multiple cell lines, the patient is associated to the tumour type of the cell line making up  the largest fraction  of the cell line composition. (PPV: Positive Predictive Value; AUC: Area Under the Curve). Source data are provided in a Source data file.
<hr/>

![alt text](https://github.com/dibbelab/singlecell_bcatlas/blob/main/figures/Figure_03.png?raw=true)
<b>Figure 3 – Transcriptional heterogeneity in breast cancer cell lines and its impact on drug response.</b> (A) Percentage of cells expressing the indicated genes in each of the 32 cell lines. (B) Fluorescence cytometry of HCC38, MDA-MB-361 and AU565 cell lines stained with a fluorescent antibody against HER2. (C) Expression of HER2 protein in MDA-MB-361  cells is dynamic and re-established in less than 3 weeks. (D) Cell cycle phase for the HER2+ and HER2- subpopulations of MDA-MB-361  cells. p-value refers to the Fisher’s exact test. (E) Enriched pathways (GSEA, FDR<10%) across differentially expressed genes between the HER2+ (orange) and HER2- (blue) MDA-MB-361 cells. (F) Gene expression versus drug potency for four anti-HER2 drugs. Each dot corresponds to a cell line with percentage of cells expressing ERBB2 or EGFR in the cell line [y-axis] versus the experimental drug potency2 as Area Under the Curve (AUC) [x-axis]. PCC (Pearson correlation coefficient) and its p-value are also shown. (G) PCC values  computed as in F for 66 drugs for which the cognate drug targets is known. The  PCC distribution when choosing a random gene is also shown. (H) Bioinformatics pipeline for the identification of drug sensitivity biomarkers for 450 drugs. (I) The top 250 most expressed genes in  a single cell are used as input for a GSEA against the ranked list of genes correlated with drug potency for each one of the 450 drugs to predict its drug sensitivity. (J) Performance of DREEP in predicting drug sensitivity of 32 cell lines in the atlas to 450 drugs in terms of PPV (Positive Predicted Value) versus Recall. (K) Dose-response curve in terms of cell viability following treatment with either afatinib or etoposide at the indicated concentrations on sorted MDA-MB-361  cells (triplicate experiment). (L) Percentage of HER2+ cells in MDA-MB-361  after 72h treatment with either afatinib or etoposide, and (M) cell viability. (N) Percentage of HER2+ cells in MDA-MB-361 cell-line at the indicated time-points either following 48h of afatinib pre-treatment (red bars) or without any afatinib pre-treatment (black bars) and (O) the relative number of cells rescaled for the number of cells at the beginning of the experiment. Source data are provided in a Source data file.
<hr/>

## How to download RAW data of the Atlas
Raw UMI counsts before to be normalized with [GFICF package](https://github.com/dibbelab/gficf) are arevailable on figshare at following DOI [10.6084/m9.figshare.15022698](https://figshare.com/articles/dataset/Single_Cell_Breast_Cancer_cell-line_Atlas/15022698). Specifically raw UMI count matrix of the 35,276 cells is available in two main formats:
1. R dataset as rds file [RAW.UMI.counts.BC.cell.lines.rds](https://figshare.com/ndownloader/files/28893384)
2. MatrixMarket format (like cell ranger output)
    1. Matrix file [matrix.mtx.gz](https://figshare.com/ndownloader/files/30469062)
    2. Cell name file file [barcodes.tsv.gz](https://figshare.com/ndownloader/files/30469065)
    3. Gene name file [features.tsv.gz](https://figshare.com/ndownloader/files/30469065)
  
## How to download Processed data of the Atlas
Processed data with [GFICF package](https://github.com/dibbelab/gficf) are arevailable on figshare at following DOI [10.6084/m9.figshare.15022698](https://figshare.com/articles/dataset/Single_Cell_Breast_Cancer_cell-line_Atlas/15022698). Specifically you have to download the file [GFICF.processed.counts.gficf](https://figshare.com/ndownloader/files/33943715)


## DREEP predictions
First dowload the following two files
1. Pre-computed top 250 expressed genes in each sequenced cell [topGFICFgenes.rds](https://www.dropbox.com/s/0oucq3dpc7zig30/topGFICFgenes.rds?dl=0)
2. Gene to Drug IC50 correlation matrix [PCC.Drug.to.Genes.rds](https://www.dropbox.com/s/b4ukgwzb36mzx2b/PCC.Drug.to.Genes.rds?dl=0)

The script below contain the code useful to covert each cell of the Atlas in a vector of Drugs and the associated p-values.

```R
require(fgsea)
library(Matrix)

pcc.all = readRDS(file = "/path/to/PCC.Drug.to.Genes.rds")
topG = readRDS("/path/to/topGFICFgenes.rds")

M = Matrix::Matrix(data = 0,nrow = ncol(data$gficf),ncol = ncol(pcc.all))
rownames(M) = colnames(data$gficf)
colnames(M) = colnames(pcc.all)
M.pv = M
for (i in 1:ncol(pcc.all))
{
  print(i)
  s = pcc.all[,i]
  res = fgsea::fgsea(pathways = topG,stats = s,nperm = 10000,gseaParam = 0,nproc = 6)
  M[,i] = res$ES
  M.pv[,i] = res$pval
}

# print the drug-by-cell matrix containing the ES scores of each Drug
head(M)

# print the drug-by-cell matrix containing the pvalues associated to each ES score in M
head(M.pv)

```
