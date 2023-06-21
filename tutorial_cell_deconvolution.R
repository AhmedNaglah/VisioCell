library(Seurat)
library(SeuratDisk)

## Inputs

# Single Cell
KBR <- LoadH5Seurat("/orange/pinaki.sarder/ahmed.naglah/data/Deconvolution/Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_07302021.h5Seurat")

# Spatial Transcriptomics
sample = "/orange/pinaki.sarder/ahmed.naglah/data/Normalkidney57919085202301"
spatial <- Load10X_Spatial(paste0(sample,'/outs/'))

## Preprocess Spatial Transcriptomics

spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE)
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE)
spatial <- FindNeighbors(spatial, dims = 1:30)
spatial <- FindClusters(spatial, verbose = FALSE)
spatial <- RunUMAP(spatial, dims = 1:30)

## Proprocess Single Cell Data

Idents(KBR) <- KBR@meta.data$subclass.l2

KBR <- SCTransform(KBR, verbose = FALSE)


## Run Deconvolution

label = "subclass.l2"

KBR <- subset(KBR, idents = 'NA', invert=T)
anchors <- FindTransferAnchors(reference = KBR, query = spatial, normalization.method = "SCT", query.assay='SCT')
predictions.assay <- TransferData(anchorset = anchors, refdata = KBR@meta.data[[label]], prediction.assay = TRUE, dims = 1:30, weight.reduction = spatial[["pca"]])
spatial[["predictions"]] <- predictions.assay

DefaultAssay(spatial) <- "predictions"
df_pred <- spatial@assays[["predictions"]]@data
max_pred <- apply(df_pred,2,function(x) max.col(t(x),'first'))
max_pred_val <- apply(df_pred,2,function(x) max(t(x)))
max_pred <- as.data.frame(max_pred)
max_pred$Seurat_subset <- rownames(df_pred)[max_pred$max_pred]
max_pred$score <- max_pred_val
max_pred$Barcode <- rownames(max_pred)
write.csv(max_pred[,c('Barcode','Seurat_subset')],paste0(sample,'_type.transfer','.csv'),quote=F,row.names = F)
spatial@meta.data$transfer_subset <- max_pred$Seurat_subset
spatial@meta.data$transfer_subset_score <- max_pred$score
saveRDS(spatial,paste0(sample,'_type_seurat_only.RDS'))

###################################################################
###################################################################
###################################################################

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


sample = "/orange/pinaki.sarder/ahmed.naglah/data/Kidney_Visium"

kidney <- readRDS(paste0(sample,'_seurat_only.RDS'))
DefaultAssay(kidney) <- "predictions"
ftr <-      "dPOD"

plot <- SpatialFeaturePlot(kidney, features = ftr) + theme(legend.position = "right")

SpatialFeaturePlot(kidney, features = ftr, interactive = TRUE)
