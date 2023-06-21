# VisioCell

Tutorial for integrating single cell and Visium data to predict cell type proportions using Seurat

## About

Spatial Transcriptomics data can provide a great overview of the molecular signature overlayed on histology section.

However, infering cell types can be challenging within a specific spot. Each spot can contains multiple cell types.

Therefor, in this pipeline, we are using a deconvolution method to infer the probability of having a specific cell at a certain spot. From these probabilities, we can derive the cell type proportion.  

![alt text](https://github.com/AhmedNaglah/VisioCell/blob/main/dPOD.jpg?raw=true)

Predicted dPOD proportion in each SPOT on a Visium Slide

For more information about the topic, you can review the below article

https://www.10xgenomics.com/resources/analysis-guides/integrating-single-cell-and-visium-spatial-gene-expression-data

For more information about spatial transcriptomics

https://www.10xgenomics.com/products/spatial-gene-expression

For more information about single cell RNA-Seq

https://www.10xgenomics.com/products/single-cell-gene-expression

## Usage

### Environment:

RStudio 2023.06.0+421 "Mountain Hydrangea" Release (583b465ecc45e60ee9de085148cd2f9741cc5214, 2023-06-06) for Ubuntu Jammy

R version 4.2.3 (2023-03-15)

### Dependencies:

Seurat

SeuratData

ggplot2

patchwork

dplyr

### Input:

```R
# Single Cell
KBR <- LoadH5Seurat("../Single_Cell_Reference.h5Seurat")

# Spatial Transcriptomics
sample = "../Visium_Sample"
spatial <- Load10X_Spatial(paste0(sample,'/outs/'))

```
