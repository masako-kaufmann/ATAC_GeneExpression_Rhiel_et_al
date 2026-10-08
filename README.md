# Chromatin accessibility and gene expression at CRISPR–Cas9 off-target sites

This repository contains the analysis code used to annotate on- and off-target sites for three CRISPR–Cas9 guide RNAs (CCR5, EMX1, FANCF) with (i) chromatin accessibility from ATAC-seq and (ii) expression of the genes overlapping those sites, using public ENCODE data.


Two independent scripts, each runnable on its own:

| Script | Question it answers |
| --- | --- |
| `ATAC.R` | Is each on-/off-target site located in open chromatin? |
| `GeneExpression.R` | Are the genes at off-target sites expressed? |

