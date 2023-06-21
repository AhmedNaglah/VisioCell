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

System: x86_64, linux-gnu
<br>R version 4.2.3 (2023-03-15)
<br>R version 4.2.3 (2023-03-15)

### Dependencies:

Seurat
<br>SeuratData
<br>ggplot2
<br>patchwork
<br>dplyr

### Input:

Visium and Chromium Data Files

```R
# Single Cell
KBR <- LoadH5Seurat("../Single_Cell_Reference.h5Seurat")

# Spatial Transcriptomics
sample = "../Visium_Sample"
spatial <- Load10X_Spatial(paste0(sample,'/outs/'))

```

### Output:

Seurat object in RDS file

```R
saveRDS(spatial,paste0(sample,'_cell_type_seurat.RDS'))

```


