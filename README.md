# singlecell_bcatlas
Soon available..

# DREEP prediction for reviwer

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
