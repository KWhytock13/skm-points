
####Integration comparisons
## Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(STACAS)
library(SeuratWrappers)
library(harmony)
library(lisi)
library(pdfCluster)
library(grid)

###Set working directory to where you want your files saving. 
##Recommend to link this to onedrive if you are logged in 
setwd("~/OneDrive - AdventHealth/Points/Data/Integration")


sample_name = 'decont_Points'
filePrefix <- "Points"
Participant <- "orig.ident"
Group <- 'Group'


###Load individual RDS seurat objects 
##Read RDS
M3 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M3decont.RDS")
M22 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M22decont.RDS")
M28 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M28decont.RDS")
M34 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M34decont.RDS")
M42 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M42decont.RDS")
M45 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M45decont.RDS")
M57 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M57decont.RDS")
M58 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M58decont.RDS")
M62 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M62decont.RDS")
M78 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M78decont.RDS")
M89 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M89decont.RDS")
M91 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M91decont.RDS")
M94 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M94decont.RDS")
M102 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M102decont.RDS")
M119 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M119decont.RDS")
M120 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M120decont.RDS")
M130 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M130decont.RDS")
M135 <- readRDS("~/OneDrive - AdventHealth/Points/Data/IndividualSeurat/seurat_M135decont.RDS")


#Merge data sets
Merge <- merge(M3, y = c(M22, M28, M34, M42, M45, M57, M58, M62, M78, M89, M91, M94, M102, M119, M120, M130, M135), project = "Points") 
DefaultAssay(object = Merge) <- "Pcoding"


###save Merged seurat. 
saveRDS(Merge, file = paste0('seurat_', sample_name, 'decont', '.RDS'))

rm(list = ls())

Merge <- readRDS("seurat_decont_Pointsdecont.RDS")
Participant <- "orig.ident"
Group <- 'Group'

##Run PCA with SCT
options(future.globals.maxSize = 8000 * 1024^2)
Run_PCA <- function(seurat, assay, number_of_features){
  DefaultAssay(object = seurat) <- assay
  seurat <- SCTransform(seurat, vst.flavor = "v2", assay = assay, verbose = FALSE, variable.features.n = number_of_features)
  seurat <- RunPCA(seurat, verbose = FALSE)
  return(seurat)
}

Merge <- Run_PCA(Merge, assay = 'Pcoding', number_of_features = 5000)
ElbowPlot(Merge)


##Cluster 
cluster <- function(seurat, PCA_dim, res, sample_name){
  DefaultAssay(object = seurat) <- 'SCT'
  seurat <- RunUMAP(seurat, dims = PCA_dim)
  seurat <- RunTSNE(seurat, dims = PCA_dim)
  seurat <- FindNeighbors(seurat, dims = PCA_dim)
  seurat <- FindClusters(seurat, resolution = res)
  return(seurat)
}

seurat <- cluster(seurat=Merge, PCA_dim = 1:15, res = c(0.2, 0.4, 0.6, 0.8, 1, 1.2), sample_name = sample_name)
DimPlot(seurat, split.by = 'orig.ident', group.by = 'orig.ident', ncol = 6)



##Generic muscle tissue markers 
Markers <- c('PAX7', #MuSC
             'PDGFRA', #FAP
             'TNXB', #FAP
             'ABCA8', #FAP
             'CFD', #FAP
             'COL6A1', #FAP
             'HMCN2', #FAP
             'LUM', #FAP
             'FBN1', #FAP
             'MYH7', #myonuclei
             "ATP2A2", #Type I
             "GAPDH", #Type II
             'MYH2', #myonuclei
             'MYH1', #myonuclei
             'MYH7B', #myonuclei
             'PECAM1', #endothelial
             'CD34', #endothelial
             'BTNL9', #endothelial
             'EGFL7', #endothelial
             'VWF', #endothelial
             'ACTA2', #SMC
             'MYH11', 
             'DES', ##Pericyte
             'PTPRC', #Immune
             'FCGR1A', #Macrophage
             'MRC1', #Macrophage
             'ITGAX', #Macrophage
             'TREM2', #Macrophage
             'CD86', #Macrophage,
             'GZMA', #NK
             'NKG7', #NK
             'CD3E', #Tcell
             'CD8A', #Tcell
             'KIT'
)

DefaultAssay(object = seurat) <- "RNA"
seurat <- NormalizeData(seurat)
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

###decide which resolution to use based on expression of makers
P1 <- DotPlot(seurat, features = Markers, 
              scale.min = 0,
              group.by = "SCT_snn_res.0.2"
) +
  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))

P2 <- DimPlot(seurat, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.2" )

P1 + P2


# plot data
P1 <- DimPlot(seurat, group.by = c(Participant)) 
P2 <- DimPlot(seurat, group.by = c(Group))
P3 <- DimPlot(seurat, group.by = c("SCT_snn_res.0.2")) + ggtitle('Cell Clusters')
ggp_all <- (P1 + P2 + P3)  + 
  plot_layout(ncol = 3) + # Create grid of plots with title
  plot_annotation(title = "NoBatchCorrection.SCT") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  plot_layout(widths = c(1, 1, 1), heights = unit(c(8, 1), c('cm', 'null')))
ggp_all   
ggsave("uncorrected.SCT.png")

# Add umap coordinates, extract meta and save for later
umap <- data.frame(seurat@reductions$umap@cell.embeddings)
colnames(umap) <- paste0("Uncorrected_SCT",  colnames(umap))
seurat <- AddMetaData(seurat, metadata = umap)
meta <- seurat@meta.data
write.csv(meta, 'meta.data.uncorrected.SCT.csv')

rm(list = ls())

################Run RPCA ######################
rm(list = ls())
Merge <- readRDS("seurat_decont_Pointsdecont.RDS")

start <- Sys.time()
###Split object by participant 
Merge.list <- SplitObject(Merge, split.by = "orig.ident") 


#normalize and identify variable features for each dataset independently
for (i in 1:length(Merge.list)) {
  suppressWarnings(Merge.list[[i]] <- SCTransform(Merge.list[[i]], vst.flavor = "v2", assay = 'Pcoding', method = "glmGamPoi",
                                                  variable.features.n = 5000, ##Can change this number if needed, default to 5000
                                                  verbose = FALSE, return.only.var.genes=FALSE))
}

features <- SelectIntegrationFeatures(object.list = Merge.list, nfeatures = 5000)
Merge.list <- PrepSCTIntegration(object.list = Merge.list, anchor.features = features)
Merge.list <- lapply(X = Merge.list, 
                     FUN = RunPCA,
                     features = features)

##find intergration anchors
anchors <- FindIntegrationAnchors(object.list = Merge.list, normalization.method = "SCT",
                                  anchor.features = features,
                                  dims = 1:20, 
                                  reduction = "rpca", 
                                  k.anchor = 10, 
) 
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                              #k.weight = 100,
                              dims = 1:20)

DefaultAssay(combined.sct) <- "integrated" ##should probably run it on integrated data but doesn't really affect the results
combined.sct <- ScaleData(combined.sct, 
                          #assay = "SCT", 
                          verbose = FALSE)
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
ElbowPlot(combined.sct) 

####Keep the same numner of PCAs 1:10
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:15)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:15)
combined.sct <- FindClusters(combined.sct, resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2))

end <- Sys.time()
end - start # Time difference of 

Participant <- "orig.ident"
Group <- 'Group'

DimPlot(combined.sct , split.by = 'orig.ident', group.by = 'orig.ident', ncol = 6)

DefaultAssay(object = combined.sct ) <- "RNA"
combined.sct  <- NormalizeData(combined.sct )
all.genes <- rownames(combined.sct )
combined.sct  <- ScaleData(combined.sct , features = all.genes)

Markers <- c('PAX7', #MuSC
             'PDGFRA', #FAP
             'TNXB', #FAP
             'ABCA8', #FAP
             'CFD', #FAP
             'COL6A1', #FAP
             'HMCN2', #FAP
             'LUM', #FAP
             'FBN1', #FAP
             'MYH7', #myonuclei
             "ATP2A2", #Type I
             "GAPDH", #Type II
             'MYH2', #myonuclei
             'MYH1', #myonuclei
             'MYH7B', #myonuclei
             'PECAM1', #endothelial
             'CD34', #endothelial
             'BTNL9', #endothelial
             'EGFL7', #endothelial
             'VWF', #endothelial
             'ACTA2', #SMC
             'MYH11', 
             'DES', ##Pericyte
             'PTPRC', #Immune
             'FCGR1A', #Macrophage
             'MRC1', #Macrophage
             'ITGAX', #Macrophage
             'TREM2', #Macrophage
             'CD86', #Macrophage,
             'GZMA', #NK
             'NKG7', #NK
             'CD3E', #Tcell
             'CD8A', #Tcell
             'KIT'
)


###decide which resolution to use based on expression of makers
P1 <- DotPlot(combined.sct , features = Markers, 
              scale.min = 0,
              group.by = "integrated_snn_res.0.2"
) +
  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))

P2 <- DimPlot(combined.sct , reduction = "umap", label = TRUE, group.by = "integrated_snn_res.0.2" )

P1 + P2


# plot data
P1 <- DimPlot(combined.sct, group.by = c(Participant)) 
P2 <- DimPlot(combined.sct, group.by = c(Group))
P3 <- DimPlot(combined.sct, group.by = c("integrated_snn_res.0.2")) + ggtitle('Cell Clusters')
ggp_all <- (P1 + P2 + P3)  + 
  plot_layout(ncol = 3) + # Create grid of plots with title
  plot_annotation(title = "RPCA_k.anchor.10") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  plot_layout(widths = c(1, 1, 1, 1), heights = unit(c(8, 1), c('cm', 'null')))
ggp_all   
ggsave("RPCA_k.anchor.10.png")

DimPlot(combined.sct, split.by = c(Participant))
DimPlot(combined.sct, group.by = c("integrated_snn_res.0.2"), split.by = Participant)
DimPlot(combined.sct, group.by = c("integrated_snn_res.0.2"), split.by = Group)

# Add umap coordinates, extract meta and save for later
umap <- data.frame(combined.sct@reductions$umap@cell.embeddings)
colnames(umap) <- paste0("RPCA_k.anchor.10",  colnames(umap))
combined.sct <- AddMetaData(combined.sct, metadata = umap)
meta <- combined.sct@meta.data
write.csv(meta, 'meta.RPCA_k.anchor.SCT.10.csv')

######Run STACAS #######
rm(list = ls())
Merge <- readRDS("seurat_decont_Pointsdecont.RDS")

start <- Sys.time()
SCT <- replicate(18, "SCT") ##The numner of replicates needs to be how many samples you have. 
ref.list <- SplitObject(Merge, split.by = "orig.ident")
ref.list <- lapply(ref.list, FUN = SCTransform, variable.features.n = 5000)
features.sct <- SelectIntegrationFeatures(ref.list, nfeatures = 5000)


object_integrated <- Merge %>% SplitObject(split.by = "orig.ident") %>%
  Run.STACAS(dims = 1:20, anchor.features = features.sct, assay = SCT) %>%
  RunUMAP(dims = 1:10) 
ElbowPlot(object_integrated)

object_integrated  <- FindNeighbors(object_integrated , reduction = "pca", dims = 1:10)
object_integrated  <- FindClusters(object_integrated , resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2))
object_integrated@meta.data

end <- Sys.time()
end - start # Time difference of 

DimPlot(object_integrated, split.by = 'orig.ident', group.by = 'orig.ident', ncol = 6)

DefaultAssay(object = object_integrated) <- "RNA"
object_integrated <- NormalizeData(object_integrated)
all.genes <- rownames(object_integrated)
object_integrated <- ScaleData(object_integrated, features = all.genes)

Markers <- c('PAX7', #MuSC
             'PDGFRA', #FAP
             'TNXB', #FAP
             'ABCA8', #FAP
             'CFD', #FAP
             'COL6A1', #FAP
             'HMCN2', #FAP
             'LUM', #FAP
             'FBN1', #FAP
             'MYH7', #myonuclei
             "ATP2A2", #Type I
             "GAPDH", #Type II
             'MYH2', #myonuclei
             'MYH1', #myonuclei
             'MYH7B', #myonuclei
             'PECAM1', #endothelial
             'CD34', #endothelial
             'BTNL9', #endothelial
             'EGFL7', #endothelial
             'VWF', #endothelial
             'ACTA2', #SMC
             'MYH11', 
             'DES', ##Pericyte
             'PTPRC', #Immune
             'FCGR1A', #Macrophage
             'MRC1', #Macrophage
             'ITGAX', #Macrophage
             'TREM2', #Macrophage
             'CD86', #Macrophage,
             'GZMA', #NK
             'NKG7', #NK
             'CD3E', #Tcell
             'CD8A', #Tcell
             'KIT'
)


###decide which resolution to use based on expression of makers
P1 <- DotPlot(object_integrated, features = Markers, 
              scale.min = 0,
              group.by = "integrated_snn_res.0.2"
) +
  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))

P2 <- DimPlot(object_integrated, reduction = "umap", label = TRUE, group.by = "integrated_snn_res.0.2" )

P1 + P2


Participant <- "orig.ident"
Group <- 'Group'

# plot data
P1 <- DimPlot(object_integrated, group.by = c(Participant)) 
P2 <- DimPlot(object_integrated, group.by = c(Group))
P3 <- DimPlot(object_integrated, group.by = c("integrated_snn_res.0.2")) + ggtitle('Cell Clusters')
ggp_all <- (P1 + P2 + P3)  + 
  plot_layout(ncol = 3) + # Create grid of plots with title
  plot_annotation(title = "STACAS.SCT.Anchor.integrate.dims.20") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  plot_layout(widths = c(1, 1, 1), heights = unit(c(8, 1), c('cm', 'null')))
ggp_all   
ggsave("STACAS.SCT.Anchor.integrate.dims.20.png")

DimPlot(object_integrated, split.by = c(Group), group.by = c("integrated_snn_res.0.2")) 


# Add umap coordinates, extract meta and save for later
umap <- data.frame(object_integrated@reductions$umap@cell.embeddings)
colnames(umap) <- paste0("STACAS",  colnames(umap))
object_integrated <- AddMetaData(object_integrated, metadata = umap)
meta <- object_integrated@meta.data
write.csv(meta, 'meta.STACA.SCT.Anchor.integrate.dims.20.csv')

######Run fastMNN#########
rm(list = ls())
Merge <- readRDS("seurat_decont_Pointsdecont.RDS")

start <- Sys.time()
SCT <- replicate(20, "SCT")
ref.list <- SplitObject(Merge, split.by = "orig.ident")
ref.list <- lapply(ref.list, FUN = SCTransform, variable.features.n = 5000)
features.sct <- SelectIntegrationFeatures(ref.list, nfeatures = 5000)


merge.list <- SplitObject(Merge, split.by="orig.ident")
merge.list <- lapply(X = merge.list, 
                     FUN = SCTransform, 
                     method = "glmGamPoi", 
                     return.only.var.genes = FALSE)
var.features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 5000)


Merge <- readRDS("seurat_decont_Pointsdecont.RDS")
ref.list <- SplitObject(Merge, split.by = "orig.ident")
ref.list <- lapply(ref.list, FUN = SCTransform, variable.features.n = 5000)



Merge <- RunFastMNN(object.list = ref.list, features = features.sct, assay ='SCT')
Merge <- RunUMAP(Merge, reduction = "mnn", dims = 1:10)
Merge <- FindNeighbors(Merge, reduction = "mnn", dims = 1:10)
Merge <- FindClusters(Merge)

end <- Sys.time()
end - start # Time difference of 

DimPlot(Merge, split.by = 'orig.ident', group.by = 'orig.ident', ncol = 6)

##Generic muscle tissue markers 
Markers <- c('PAX7', #MuSC
             'PDGFRA', #FAP
             'TNXB', #FAP
             'ABCA8', #FAP
             'CFD', #FAP
             'COL6A1', #FAP
             'HMCN2', #FAP
             'LUM', #FAP
             'FBN1', #FAP
             'MYH7', #myonuclei
             "ATP2A2", #Type I
             "GAPDH", #Type II
             'MYH2', #myonuclei
             'MYH1', #myonuclei
             'MYH7B', #myonuclei
             'PECAM1', #endothelial
             'CD34', #endothelial
             'BTNL9', #endothelial
             'EGFL7', #endothelial
             'VWF', #endothelial
             'ACTA2', #SMC
             'MYH11', 
             'DES', ##Pericyte
             'PTPRC', #Immune
             'FCGR1A', #Macrophage
             'MRC1', #Macrophage
             'ITGAX', #Macrophage
             'TREM2', #Macrophage
             'CD86', #Macrophage,
             'GZMA', #NK
             'NKG7', #NK
             'CD3E', #Tcell
             'CD8A', #Tcell
             'KIT'
)

DefaultAssay(object = Merge) <- "RNA"
Merge <- NormalizeData(Merge)
all.genes <- rownames(object_integrated)
Merge <- ScaleData(Merge, features = all.genes)

###decide which resolution to use based on expression of makers
P1 <- DotPlot(Merge, features = Markers, 
              scale.min = 0
) +
  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))

P2 <- DimPlot(Merge, reduction = "umap", label = TRUE )

P1 + P2

Participant <- "orig.ident"
Group <- 'Group'

# plot data
P1 <- DimPlot(Merge, group.by = c(Participant)) 
P2 <- DimPlot(Merge, group.by = c(Group))
P3 <- DimPlot(Merge) + ggtitle('Cell Clusters')
ggp_all <- (P1 + P2 + P3)  + 
  plot_layout(ncol = 3) + # Create grid of plots with title
  plot_annotation(title = "FastMNN.SCT") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  plot_layout(widths = c(1, 1, 1), heights = unit(c(8, 1), c('cm', 'null')))
ggp_all   
ggsave("FastMNN.SCT.png")


# Add umap coordinates, extract meta and save for later
umap <- data.frame(Merge@reductions$umap@cell.embeddings)
colnames(umap) <- paste0("FastMNN.SCT",  colnames(umap))
Merge <- AddMetaData(Merge, metadata = umap)
meta <- Merge@meta.data
write.csv(meta, 'meta.FastMNN.SCT.csv')

#####run harmony ######
rm(list = ls())
Merge <- readRDS("seurat_decont_Pointsdecont.RDS")

start <- Sys.time()

merge.list <- SplitObject(Merge, split.by="orig.ident")
merge.list <- lapply(X = merge.list, 
                     FUN = SCTransform, 
                     method = "glmGamPoi", 
                     return.only.var.genes = FALSE)
var.features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 5000)

merge.sct <- merge(x = merge.list[[1]], y = merge.list[2:length(merge.list)], merge.data=TRUE)
VariableFeatures(merge.sct) <- var.features
merge.sct <- RunPCA(merge.sct, verbose = FALSE)
merge.sct <- RunHarmony(merge.sct, assay.use="SCT", group.by.vars = "orig.ident")
ElbowPlot(merge.sct)
merge.sct <- RunUMAP(merge.sct, reduction = "harmony", dims = 1:15)
merge.sct <- FindNeighbors(merge.sct, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4)) 
end <- Sys.time()
end - start # Time difference of 

Participant <- "orig.ident"
Group <- 'Group'


DimPlot(merge.sct, split.by = 'orig.ident', group.by = 'orig.ident', ncol = 6)

##Generic muscle tissue markers 
Markers <- c('PAX7', #MuSC
             'PDGFRA', #FAP
             'TNXB', #FAP
             'ABCA8', #FAP
             'CFD', #FAP
             'COL6A1', #FAP
             'HMCN2', #FAP
             'LUM', #FAP
             'FBN1', #FAP
             'MYH7', #myonuclei
             "ATP2A2", #Type I
             "GAPDH", #Type II
             'MYH2', #myonuclei
             'MYH1', #myonuclei
             'MYH7B', #myonuclei
             'PECAM1', #endothelial
             'CD34', #endothelial
             'BTNL9', #endothelial
             'EGFL7', #endothelial
             'VWF', #endothelial
             'ACTA2', #SMC
             'MYH11', 
             'DES', ##Pericyte
             'PTPRC', #Immune
             'FCGR1A', #Macrophage
             'MRC1', #Macrophage
             'ITGAX', #Macrophage
             'TREM2', #Macrophage
             'CD86', #Macrophage,
             'GZMA', #NK
             'NKG7', #NK
             'CD3E', #Tcell
             'CD8A', #Tcell
             'KIT'
)


DefaultAssay(object = merge.sct) <- "RNA"
merge.sct <- NormalizeData(merge.sct)
all.genes <- rownames(object_integrated)
merge.sct <- ScaleData(merge.sct, features = all.genes)

###decide which resolution to use based on expression of makers
P1 <- DotPlot(merge.sct, features = Markers, 
              scale.min = 0, group.by = 'SCT_snn_res.0.2'
) +
  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))

P2 <- DimPlot(merge.sct, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.2" )

P1 + P2

# plot data
P1 <- DimPlot(merge.sct, group.by = c(Participant)) 
P2 <- DimPlot(merge.sct, group.by = c(Group))
P3 <- DimPlot(merge.sct, group.by ="SCT_snn_res.0.2" ) + ggtitle('Cell Clusters')
ggp_all <- (P1 + P2 + P3)  + 
  plot_layout(ncol = 3) + # Create grid of plots with title
  plot_annotation(title = "harmony.SCT") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  plot_layout(widths = c(1, 1, 1), heights = unit(c(8, 1), c('cm', 'null')))
ggp_all   
ggsave("harmony.SCT.png")


# Add umap coordinates, extract meta and save for later
umap <- data.frame(merge.sct@reductions$umap@cell.embeddings)
colnames(umap) <- paste0("harmony.SCT",  colnames(umap))
merge.sct <- AddMetaData(merge.sct, metadata = umap)
meta <- merge.sct@meta.data
write.csv(meta, 'meta.harmony.SCT.csv')

DimPlot(merge.sct, split.by = c(Group), group.by = c("SCT_snn_res.0.2")) 
saveRDS(merge.sct, 'Points.harmony.sct.RDS')

##Compute batch correction metrics using saved meta.data csv files. Examples below
######unadjusted
df <- read.csv('meta.data.uncorrected.SCT.csv')
head(df)

##compute Lisi
cols <- c('Uncorrected_SCTUMAP_1', 'Uncorrected_SCTUMAP_2')
df2 <- df[, colnames(df) %in% cols]

cols <- c('SCT_snn_res.0.2', 'orig.ident')
df3 <- df[, colnames(df) %in% cols]
table(df3$SCT_snn_res.0.3)
table(df3$orig.ident)

res <- compute_lisi(df2, df3, c('SCT_snn_res.0.2', 'orig.ident'))
head(res)
mean(res$orig.ident)

#compute ARI
adj.rand.index(df$SCT_snn_res.0.2, df$orig.ident)

#####RPCA######
df <- read.csv('meta.RPCA_k.anchor.SCT.10.csv')
head(df)

##compute Lisi
cols <- c('RPCA_k.anchor.10UMAP_1', 'RPCA_k.anchor.10UMAP_2')
df2 <- df[, colnames(df) %in% cols]

cols <- c('integrated_snn_res.0.2', 'orig.ident')
df3 <- df[, colnames(df) %in% cols]
table(df3$integrated_snn_res.0.2)
table(df3$orig.ident)

res <- compute_lisi(df2, df3, c('integrated_snn_res.0.2', 'orig.ident'))
head(res)
mean(res$orig.ident)

#compute ARI
adj.rand.index(df$integrated_snn_res.0.2, df$orig.ident)

#####harmony######
df <- read.csv('meta.harmony.SCT.csv')
head(df)

##compute Lisi
cols <- c('harmony.SCTUMAP_1', 'harmony.SCTUMAP_2')
df2 <- df[, colnames(df) %in% cols]

cols <- c('SCT_snn_res.0.2', 'orig.ident')
df3 <- df[, colnames(df) %in% cols]
table(df3$SCT_snn_res.0.2)
table(df3$orig.ident)

res <- compute_lisi(df2, df3, c('SCT_snn_res.0.2', 'orig.ident'))
head(res)
mean(res$orig.ident)

#compute ARI
adj.rand.index(df$SCT_snn_res.0.2, df$orig.ident)

#####MNN######
df <- read.csv('meta.FastMNN.SCT.csv')
head(df)

##compute Lisi
cols <- c('FastMNN.SCTUMAP_1', 'FastMNN.SCTUMAP_2')
df2 <- df[, colnames(df) %in% cols]

cols <- c('seurat_clusters', 'orig.ident')
df3 <- df[, colnames(df) %in% cols]
table(df3$seurat_clusters)
table(df3$orig.ident)

res <- compute_lisi(df2, df3, c('seurat_clusters', 'orig.ident'))
head(res)
mean(res$orig.ident)

#compute ARI
adj.rand.index(df$seurat_clusters, df$orig.ident)


####STACAS######
df <- read.csv('meta.STACA.SCT.Anchor.integrate.dims.20.csv')
head(df)

##compute Lisi
cols <- c('STACASUMAP_1', 'STACASUMAP_2')
df2 <- df[, colnames(df) %in% cols]

cols <- c('integrated_snn_res.0.2', 'orig.ident')
df3 <- df[, colnames(df) %in% cols]
table(df3$integrated_snn_res.0.2)
table(df3$orig.ident)

res <- compute_lisi(df2, df3, c('integrated_snn_res.0.2', 'orig.ident'))
head(res)
mean(res$orig.ident)

#compute ARI
adj.rand.index(df$integrated_snn_res.0.2, df$orig.ident)






