library(SeuratDisk)

KBR <- LoadH5Seurat("./Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_07302021_SCT.h5Seurat")

library(Seurat)

# Spatial Transcriptomics
sample = "./V42D20-364_XY01_2235505.RDS"  # Input 
sample_out = "./V42D20-364_XY01_2235505_out.RDS"  # Output

spatial <- readRDS(file = sample)

## Preprocess Spatial Transcriptomics

spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE)
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE)
spatial <- FindNeighbors(spatial, dims = 1:30)
spatial <- FindClusters(spatial, verbose = FALSE)
spatial <- RunUMAP(spatial, dims = 1:30)

Idents(KBR) <- KBR@meta.data$subclass.l2

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
spatial@meta.data$transfer_subset <- max_pred$Seurat_subset
spatial@meta.data$transfer_subset_score <- max_pred$score

saveRDS(spatial,sample_out)
