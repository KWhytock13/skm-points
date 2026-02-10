###Performing QC, clustering and decontx on each sample 
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(celda)
library(SingleCellExperiment) 
library(simpleSingleCell) 
library(scran) 
library(scuttle)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)
library(decontX)

setwd("~/Library/CloudStorage/OneDrive-AdventHealth/Clinical Studies Projects/Points/Data")

my_theme <- theme_classic() + theme(text = element_text(size=10, colour = "black"),
                                    axis.text.x = element_text(colour = "black", size = 10),
                                    axis.text.y = element_text(colour = "black", size = 10),
                                    axis.title.x = element_text(colour = "black", size = 10),
                                    axis.title.y = element_text(colour = "black", size = 10),
)



remove_duplicate_genes<-function(eset, column_of_symbol, method = "mean"){
  eset<-as.data.frame(eset)
  rownames(eset)<-NULL
  dups <- dim(eset)[1] - length(unique(eset[,column_of_symbol]))
  if(dups==0){
    eset<-tibble:: column_to_rownames(eset,var = column_of_symbol)
    return(eset)
  }else{
    if(method=="mean"){
      order_index=apply(eset[,setdiff(colnames(eset),column_of_symbol)],1,function(x) mean(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>%dplyr:: distinct(!!sym(column_of_symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = column_of_symbol)
      return(eset)
    }else if(method == "sd"){
      order_index = apply(eset[,setdiff(colnames(eset),column_of_symbol)],1,function(x) sd(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>% distinct(!!sym(column_of_symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = column_of_symbol)
      return(eset)
    }
  }
}


preprocess2 <- function(raw_counts, stats, well, gene_info, ambient_control, samples, sample_name){
  print(dim(raw_counts)) 
  raw_counts$GeneID <- gsub("_.*","",raw_counts$GeneID) ###Remove gene symbol to just show ENSG 
  sample <- subset(stats, Sample %in% c(samples)) 
  sample_barcodes <- sample$Barcode  
  sample_barcodes <- c(sample_barcodes, "GeneID") ##include GeneID as a column header for filtering 
  counts <- raw_counts[(names(raw_counts) %in%  sample_barcodes)] ##Filter the gene df so we use only sample barcodes  
  print(dim(counts)) 
  subset <- gene_info %>% dplyr::select(2,3) #subset gene info to just ensemble_Id and gene_name 
  colnames(subset) <- c("GeneID", 'symbol') ##change column names 
  counts <- left_join(counts, subset, by= c('GeneID')) ##add gene information to counts 
  countsb <- remove_duplicate_genes(counts, column_of_symbol = 'symbol', method = 'mean') ##remove duplicate genes that have the same symbol based on which one has the greatest mean 
  print(summary(duplicated(rownames(countsb))))  
  counts <- countsb[,-1]  
  print(dim(counts)) 
  well$Barcode <- gsub("\\+", "", well$Barcode) ##remove the + from the barcode in the wellist so it will match the genedf 
  ambientbarcodes <- well[well$Comment == ambient_control,]  ##subsets wells that are ambient controls 
  amb_bar <- c(ambientbarcodes$Barcode, 'GeneID') ##gets list of barcodes that belong to ambient wells, include GeneID as a column header for filtering 
  ambient_matrix <- counts[(names(counts) %in% amb_bar)] ##create ambient df to save and maybe use later 
  write.csv(ambient_matrix, file = paste0('ambient_matrix_', sample_name, '.csv')) ##and save it 
  amb_bar <- c(ambientbarcodes$Barcode) ##create list of barcodes without geneID to subset count df 
  counts <- counts%>% dplyr::select(-one_of(amb_bar)) ##remove ambient wells from gene df  
  print(dim(counts))
  keep <- colnames(counts)   #update stats file so same barcodes match gene matrix file 
  stats <- filter(stats, Barcode %in% keep)   #update stats file so same barcodes match gene matrix file  
  row <- rownames(counts)   ##Order gene info to match same order as count matrix 
  gene_info <- gene_info[match(row, gene_info$Gene_Name), ]     ###Convert gene matrices to dgcMatrix format ready to make single cell experiment 
  mat <- as.matrix(counts) 
  mat <- as(mat, "dgCMatrix")## Paste in the experiment name into the column names to make barcodes/cell IDs unique.
  colnames(mat) <- paste(sample_name, colnames(mat), sep="_")   ###Paste in the experiment name into the column names to make barcodes/cell IDs unique for stats file  
  stats2 <- stats[,-1] 
  rownames(stats2) <- stats[,1]  
  rownames(stats2) <- paste(sample_name, rownames(stats2), sep = "_")  
  preprocessed <- list(mat = mat, gene_info = gene_info, stats2 = stats2)  
  return(preprocessed)
}


create_SCE <- function(data){
  sce <- SingleCellExperiment(list(counts=data$mat), colData = data$stats2, rowData = data$gene_info) ##create single cell experiment
  sce <- addPerCellQCMetrics(sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(sce), value = FALSE)))
  pcodingmito <- c("MT-ND1", "MT-ND2",  "MT-ND3",  "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",   "MT-CO1", "MT-CO2",  "MT-CO3",  "MT-CYB",  "MT-ATP8", "MT-ATP6")
  sce <- addPerCellQCMetrics(sce, subsets=list(pcodingmito=pcodingmito, x = rownames(sce), value = FALSE))
  colData(sce)$percent_ribo <- (sce@colData$Ribosomal_Reads/sce@colData$sum) *100
  colData(sce)$Complexity <- log10(sce@colData$detected) / log10(sce@colData$sum) ##Add complexity measurement
  return(sce)
}

filter_cells <- function(sce,  sum_no, detected_low, detected_high, complexity_no, mito_percent_no, pc_mito_percent_no, ribo_percent_no, sample_name){
  sce <- sce[,(sce@colData$sum >= sum_no)] ##Remove nanowells with number of reads below this amount, we have sometimes not used this filer, or used 20,000
  sce <- sce[,(sce@colData$detected >= detected_low & sce@colData$detected <= detected_high)]  ##Remove nanowells if they have a low number of a high number of genes or high number of genes
  sce <- sce[,(sce@colData$Complexity >= complexity_no)] ##Remove nanowells with a low complexity, this means they have a lot of reads but not that many genes detected, for full-length snRNA-seq 0.65 is typical
  sce <- sce[,!(sce@colData$subsets_Mito_percent == "NaN")] 
  sce <- sce[,!(sce@colData$subsets_Mito_percent > mito_percent_no)] ####Remove nanowells with a mitochondrial % above this number, 20 is typically default
  sce <- sce[,!(sce@colData$subsets_pcodingmito_percent > pc_mito_percent_no)]
  sce <- sce[,! (sce@colData$percent_ribo > ribo_percent_no)]
  df <- as.data.frame(colData(sce))
  G1 <- ggplot(df, aes(x=sum)) + geom_histogram(bins = 50, colour = 'black', fill = 'white') + labs(title = 'Overall Counts', x = 'Counts') + my_theme
  G2 <- ggplot(df, aes(x=detected)) + geom_histogram(bins = 50, colour = 'black', fill = 'white') + labs(title = 'Genes Detected', x = 'Genes Detected') + my_theme
  G3 <- ggplot(df, aes(x=Complexity)) + geom_histogram(bins = 50, colour = 'black', fill = 'white') + labs(title = 'Complexity', x = 'Complexity') + my_theme
  G4 <- ggplot(df, aes(x=subsets_Mito_percent)) + geom_histogram(bins = 50, colour = 'black', fill = 'white') + labs(title = 'Mito %', x = 'Mito %') + my_theme
  G5 <- ggplot(df, aes(x=subsets_pcodingmito_percent)) + geom_histogram(bins = 50, colour = 'black', fill = 'white') + labs(title = 'PC Mito %', x = 'Protein Coding mito %') + my_theme
  G6 <- ggplot(df, aes(x=percent_ribo)) + geom_histogram(bins = 50, colour = 'black', fill = 'white') + labs(title = 'Ribo %', x = 'Ribosomal %') + my_theme
  G7 <-  ggplot(df, aes(x=sum, y=detected, color = subsets_pcodingmito_percent)) + geom_point() + scale_colour_gradient(low = "gray90", high = "black") + stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    labs(x = 'Counts', y = 'Genes Detected', colour = '% of PC Mito Genes') + my_theme
  G8 <-  ggplot(df, aes(x=sum, y=detected, color = percent_ribo)) + geom_point() + scale_colour_gradient(low = "gray90", high = "black") + stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    labs(x = 'Counts', y = 'Genes Detected', colour = '% of Ribo Genes') + my_theme
  G9 <- ggplot(df, aes(x = subsets_Mito_percent, subsets_pcodingmito_percent)) + geom_point() +  stat_smooth(method=lm) +
    labs(x = '% all Mito Genes', y = '% PC Mito Genes') + my_theme
  g <- arrangeGrob(G1, G2, G3, G4, G5, G6, G7, G8, G9,  nrow = 3, top = sample_name) #generates g
  ggsave(g, filename = paste0('AfterQCmetrics/AfterQC_', sample_name, '.png'))
  return(sce)
}

generate_seurat <- function(sce, counts_cut_off){
  ave.counts <- rowMeans(counts(sce))
  keep <- ave.counts >= 0.1  ##filter genes based on default >= 0.1. Log10(0.1) = -1
  sce1 <- sce
  sce1 <- sce1[keep,]
  seurat <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))  ####Convert to seurat object
  filtered <- assay(sce1) 
  filtered <- CreateAssayObject(counts = filtered) ##add unfiltered assay for data exploration, this is not used for clustering etc
  seurat[["filtered"]] <- filtered
  ##Add feature level meta data to both assay slots
  seurat[["filtered"]] <- AddMetaData(seurat[["filtered"]],
                                      as.data.frame(SingleCellExperiment::rowData(sce1))) 
  seurat[["RNA"]] <- AddMetaData(seurat[["RNA"]],
                                 as.data.frame(SingleCellExperiment::rowData(sce)))
  
  ##Get  protein coding genes only for clustering
  genes <- seurat@assays$filtered@meta.features 
  PC <- genes %>% dplyr::filter(Gene_Biotype == 'protein_coding') ##Protein coding genes DF
  PCoding <- PC$Gene_Name #Protein coding genes list
  ##Make a new assay with protein coding genes only 
  filtered <- as.data.frame(seurat@assays$filtered@counts) 
  pcode <- filtered[rownames(filtered) %in% PCoding,]
  pcoding <- CreateAssayObject(counts = pcode)
  seurat[["Pcoding"]] <- pcoding
  print(Assays(seurat))
  return(seurat)
}

generate_seurat2 <- function(sce, counts_cut_off, am_genes){
  ave.counts <- rowMeans(counts(sce))
  keep <- ave.counts >= 0.1  ##filter genes based on default >= 0.1. Log10(0.1) = -1
  sce1 <- sce
  sce1 <- sce1[keep,]
  seurat <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))  ####Convert to seurat object
  filtered <- assay(sce1) 
  filtered <- CreateAssayObject(counts = filtered) ##add unfiltered assay for data exploration, this is not used for clustering etc
  seurat[["filtered"]] <- filtered
  ##Add feature level meta data to both assay slots
  seurat[["filtered"]] <- AddMetaData(seurat[["filtered"]],
                                      as.data.frame(SingleCellExperiment::rowData(sce1))) 
  seurat[["RNA"]] <- AddMetaData(seurat[["RNA"]],
                                 as.data.frame(SingleCellExperiment::rowData(sce)))
  
  ##Get  protein coding genes only for clustering
  genes <- seurat@assays$filtered@meta.features 
  PC <- genes %>% dplyr::filter(Gene_Biotype == 'protein_coding') ##Protein coding genes DF
  PCoding <- PC$Gene_Name #Protein coding genes list
  ##Make a new assay with protein coding genes only 
  filtered <- as.data.frame(seurat@assays$filtered@counts) 
  pcode <- filtered[rownames(filtered) %in% PCoding,]
  pcode <- pcode[! rownames(pcode) %in% am_genes,] ###Remove ambient genes from pcoding assay so data isnt clustered by it 
  pcoding <- CreateAssayObject(counts = pcode)
  seurat[["Pcoding"]] <- pcoding
  print(Assays(seurat))
  return(seurat)
}


mito_pseudo <- readLines("mito_pseudogenes.csv")


remove_genes <- function(seurat){
  seurat <- seurat[!grepl("^MT-", rownames(seurat)), ] #filter mitochondrial genes
  seurat <- seurat[!grepl("^HB[^(P)]", rownames(seurat)), ] 
  #seurat <- seurat[!grepl("MALAT1", rownames(seurat)), ] 
  #seurat <- seurat[!grepl("NEAT1", rownames(seurat)), ] 
  seurat <- subset(seurat,features=setdiff(rownames(seurat),mito_pseudo))
  return(seurat)
}

Run_PCA <- function(seurat, assay, number_of_features, sample_name){
  DefaultAssay(object = seurat) <- assay
  seurat <- SCTransform(seurat, vst.flavor = "v2", assay = assay, verbose = FALSE, variable.features.n = number_of_features)
  seurat <- RunPCA(seurat, verbose = FALSE)
  G1 <- DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
  ggsave(G1, filename = paste0('PCA_plots/PCA_Heatmap_', sample_name, '.png'))
  G2 <- ElbowPlot(seurat)
  ggsave(G2, filename = paste0('PCA_plots/PCA_Elbow_', sample_name, '.png'))
  return(seurat)
}
cluster <- function(seurat, PCA_dim, res, sample_name){
  DefaultAssay(object = seurat) <- 'SCT'
  seurat <- RunUMAP(seurat, dims = PCA_dim)
  seurat <- RunTSNE(seurat, dims = PCA_dim)
  seurat <- FindNeighbors(seurat, dims = PCA_dim)
  seurat <- FindClusters(seurat, resolution = res)
  return(seurat)
}


run_decontx_with_background <- function(seurat, data, raw){
  sce1 <- SingleCellExperiment(list(counts = raw))
  sce <- SingleCellExperiment(list(counts = seurat@assays$RNA@counts), 
                              colData = seurat@meta.data, rowData = data$gene_info)
  sce <- decontX(sce, z = colData(sce)$CellType, background = sce1)
  return(sce)
}

create_seurat_with_decont <- function(sce, counts_cut_off){
  seurat1 <- CreateSeuratObject(counts = decontXcounts(sce), meta.data = as.data.frame(colData(sce)))
  seurat1[["RNA"]] <- AddMetaData(seurat1[["RNA"]],
                                  as.data.frame(SingleCellExperiment::rowData(sce)))
  cols <- c('sum.1', 'detected.1', 'total.1', 'nCount_filtered', 'nFeature_filtered', 'nCount_Pcoding', 'nCount_SCT', 'nFeature_SCT', 
            'SCT_snn_res.0.2', 'SCT_snn_res.0.3', 'SCT_snn_res.0.4', 'SCT_snn_res.0.5', 'SCT_snn_res.0.6', 'SCT_snn_res.0.7', 
            'SCT_snn_res.0.8', 'SCT_snn_res.0.9', 'seurat_clusters', 'CellType', 'ident', 'sum.2', 'detected.2', 'subsets_Mito_sum.1', 
            'subsets_Mito_detected.1', 'subsets_Mito_percent.1', 'altexps_filtered_sum', 'altexps_filtered_detected', 'altexps_filtered_percent',
            'altexps_Pcoding_sum', 'altexps_Pcoding_detected','altexps_Pcoding_percent',  "altexps_SCT_sum", "altexps_SCT_detected", "altexps_SCT_percent",
            "total.2",  "decontX_clusters")
  seurat1@meta.data <- seurat1@meta.data[, ! colnames(seurat1@meta.data) %in% cols]
  ave.counts <- rowMeans(decontXcounts(sce))
  keep <- ave.counts >= counts_cut_off  ##filter genes based on default >= 0.1. Log10(0.1) = -1
  sce1 <- sce
  sce1 <- sce1[keep,]
  counts = decontXcounts(sce1)
  filtered <- CreateAssayObject(counts = counts) ##add unfiltered assay for data exploration, this is not used for clustering etc
  seurat1[["filtered"]] <- filtered
  seurat1@assays$filtered@counts[1:6, 1:6]
  ##Add feature level meta data to both assay slots
  seurat1[["filtered"]] <- AddMetaData(seurat1[["filtered"]],
                                       as.data.frame(SingleCellExperiment::rowData(sce1))) 
  seurat1@assays$filtered@counts[rownames(seurat1@assays$filtered@counts) %in% 'CD36', ][1:10] 
  ##Get  protein coding genes only for clustering
  genes <- seurat1@assays$filtered@meta.features 
  PC <- genes %>% dplyr::filter(Gene_Biotype == 'protein_coding') ##Protein coding genes DF
  PCoding <- PC$Gene_Name #Protein coding genes list
  ##Make a new assay with protein coding genes only 
  filtered <- as.data.frame(seurat1@assays$filtered@counts) 
  pcode <- filtered[rownames(filtered) %in% PCoding,]
  pcoding <- CreateAssayObject(counts = pcode)
  seurat1[["Pcoding"]] <- pcoding
  return(seurat1)
}
Add_Well_To_Meta <- function(well, seurat, sample_name){
  well$Barcode <- gsub("\\+", "", well$Barcode)
  well$Barcode  <- paste(sample_name, well$Barcode, sep = "_")
  well2 <- well %>% dplyr::select(Barcode, Sample, Size1, Circularity1, Integ.Signal1, SampleWell, Dispense.tip)
  barcodes <- rownames(seurat@meta.data)
  data <- filter(well2, Barcode %in% barcodes) 
  rownames(data) <- data$Barcode
  seurat <- AddMetaData(seurat, data)
  return(seurat)
}


##Load filtered data from cogentAP for preliminary analysis,  ambient matrix pre-generated by running preprocess2
raw_counts <- read.table("~/OneDrive - AdventHealth/Clinical Studies Projects/Points/Data/CogentAP/M62/M62_analysis_incl_introns_genematrix.csv", header = T, sep = ',')
stats <- read.csv("~/OneDrive - AdventHealth/Clinical Studies Projects/Points/Data/CogentAP/M62/M62_analysis_incl_introns_stats.csv")
well <- read.table("~/OneDrive - AdventHealth/Clinical Studies Projects/Points/Data/WellList files/141021_MIRAGE0062_V2_WellList.TXT",  header = T,  sep="\t")
gene_info <- read.csv("~/OneDrive - AdventHealth/Clinical Studies Projects/Points/Data/CogentAP/M62/gene_info_incl_introns.csv")
raw <- read.csv("~/OneDrive - AdventHealth/Clinical Studies Projects/Points/Data/Ambient/ambient_matrix_M62.csv", row.names = 1)
am_genes <- read.csv("~/OneDrive - AdventHealth/Clinical Studies Projects/Points/Data/Ambient/ambient_genes_M62.csv")
raw <- as.matrix(raw) 
raw <- as(raw, "dgCMatrix")


ambient_process <- function(df){
df <- subset(df, Row_means > 5) ##remove genes that have average expression lower than 5
df <- df[-ncol(df)]
df <- df %>% remove_rownames %>% column_to_rownames(var="X")
df$per_NA <- rowSums(is.na(df)) / ncol(df)
df <- subset(df, per_NA < 0.2) ##subset to be genes that are present in atleast 80% nanowells
am_genes <- rownames(df)
return(am_genes)
}

am_genes <- ambient_process(am_genes)

sample_name <- "M62"


unique(well$Comment)

data <- preprocess2(raw_counts = raw_counts, 
                    stats = stats,
                    well = well, 
                    gene_info = gene_info,
                    ambient_control = 'Ambient Control', 
                    samples = c('sample'),
                    sample_name = sample_name)


###Create a single cell experiment and plots QC metrics
sce <- create_SCE(data = data)

##Filter cells based on metric
sce <- filter_cells(sce = sce, 
                    sum_no = 10000, ##Remove nanowells with number of reads below this amount, we have sometimes not used this filter, or used 20,000
                    detected_low = 500, ##remove nanowells with that detect no of genes lower than this number, may indicate poor quality nuclei,
                    detected_high = 15000, ##remove nanowells with that detect no of genes higher than this number, could be an indicator of doublet or a lot of ambient RNA,
                    mito_percent_no = 30, ####Remove nanowells with a mitochondrial % above this number, 20 is typically default
                    pc_mito_percent_no = 15, ####Remove nanowells with a protein coding mitchondrial % above this number.
                    ribo_percent_no = 30, ###Remove nanowells with rRNA read % above this nmber
                    complexity_no = 0.6, ##remove nanowells with low complexity, i.e. a lot of reads but not many genes detected, for full-length snRNA-seq 0.65 is typical
                    sample_name = sample_name
)



###Take a look at how frequently the genes are expressed
ave.counts <- rowMeans(counts(sce))
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
     xlab=expression(Log[10]~"average count"))

seurat <- generate_seurat(sce = sce,
                          counts_cut_off = 0.01 ##filter genes based on default >= 0.1. Log10(0.1) = -1

)

##Run and plot PCA
seurat <- Run_PCA(seurat = seurat,
                  assay = "Pcoding", ##usually use used Protein coding assay but can also use RNA or unfiltered
                  number_of_features =  3000, ##number of features used for clustering
                  sample_name = sample_name)

##Cluster your seurat object 
seurat <- cluster(seurat = seurat,
                  PCA_dim = 1:7, ##number of PCA dimensions used
                  res = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), ##Different resolutions to determine the different clusters
                  sample_name = sample_name
)

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
              group.by = "SCT_snn_res.0.7"
) +
  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))

P2 <- DimPlot(seurat, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.7" )

P1 + P2

#Top 5 genes for each cluster

Idents(seurat) <- seurat@meta.data$SCT_snn_res.0.7


markersd <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Top5 <- markersd %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print(Top5, n = 60)

top5 <- Top5$gene
str(top5)
top5 <- unique(top5)

avg <- AverageExpression(seurat, group.by = "SCT_snn_res.0.7")
avg <- as.data.frame(avg$RNA)

top5_df <- subset(avg, rownames(avg) %in% top5)
top5_df <- as.data.frame(t(top5_df))
str(top5_df)
#top5_df_S <- standardize(top5_df) ##Standardize     
top5_df_S <- scale(top5_df) ##Standardize    
top5_df_S <- as.data.frame(top5_df_S)
str(top5_df_S)

top5_df_S <- tibble::rownames_to_column(top5_df_S, "SCT_snn_res.0.7")

df.long <- pivot_longer(data = top5_df_S, cols = -SCT_snn_res.0.7, names_to = "Gene", values_to = "Abundance")
head(df.long)

df.long$Gene <- factor(df.long$Gene, levels= top5)

ggplot(data = df.long, mapping = aes(x = SCT_snn_res.0.7, y = Gene, fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option="plasma") +
  coord_flip() +
  theme_classic() +
  theme(axis.title=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.ticks.x =element_blank(),
        axis.text.x = element_text(colour = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 10)) 



###Provisionally, identify cell-types, this will change from sample to sample.Example given below 

seurat@meta.data <- seurat@meta.data %>% mutate(CellType = case_when(
  SCT_snn_res.0.7 == "0"  ~ "CKM+ Myonuclei", 
  SCT_snn_res.0.7 == "1"  ~ "EGF+ Myonuclei",
  SCT_snn_res.0.7 == "2"  ~ "MYH7B+ Myonuclei",
  SCT_snn_res.0.7 == "3"  ~ "FAPs",  
  SCT_snn_res.0.7 == "4"  ~ "EC",
  SCT_snn_res.0.7 == "5"  ~ "SMC/Pericytes",
  SCT_snn_res.0.7 == "6"  ~ "Immune",
  SCT_snn_res.0.7 == "7"  ~ "MuSC"
)) 

seurat <- Add_Well_To_Meta(well = well, seurat = seurat, sample_name = sample_name)



# Create a SingleCellExperiment object and run decontX
sce <- run_decontx_with_background(seurat, data, raw)
df <- as.data.frame(colData(sce))
ggplot(df, aes(x=decontX_contamination)) + geom_histogram(bins = 50, colour = 'black', fill = 'white')


##Create a seurat object with the decontx counts
seurat1 <- create_seurat_with_decont(sce, counts_cut_off = 0.01)
seurat1 <- remove_genes(seurat1)

seurat1 <- Run_PCA(seurat = seurat1,
                   assay = "Pcoding", ##usually use used Protein coding assay but can also use RNA or unfiltered
                   number_of_features = 3000, ##number of features used for clustering
                   sample_name = sample_name)

##Cluster your seurat object based on the number of dimensions used, we can explore different resolutions
seurat1 <- cluster(seurat = seurat1,
                   PCA_dim = 1:10, ##number of PCA dimensions used
                   res = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35), ##Different resolutions to determine the different clusters
                   sample_name = sample_name
)



DefaultAssay(object = seurat1) <- "RNA"

seurat1 <- NormalizeData(seurat1)
all.genes <- rownames(seurat1)
seurat1 <- ScaleData(seurat1, features = all.genes)

P3 <- DimPlot(seurat1, reduction = "umap", group.by = 'SCT_snn_res.0.3')


P4 <- DotPlot(seurat1, features = Markers, scale.min = 0, group.by = 'SCT_snn_res.0.3')  +  coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG')))
P3 + P4


Idents(seurat1) <- seurat1@meta.data$SCT_snn_res.0.3


markersd <- FindAllMarkers(seurat1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Top5 <- markersd %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print(Top5, n = 60)

top5 <- Top5$gene
str(top5)
top5 <- unique(top5)

avg <- AverageExpression(seurat1, group.by = "SCT_snn_res.0.3")
avg <- as.data.frame(avg$RNA)

top5_df <- subset(avg, rownames(avg) %in% top5)
top5_df <- as.data.frame(t(top5_df))
str(top5_df)
#top5_df_S <- standardize(top5_df) ##Standardize     
top5_df_S <- scale(top5_df) ##Standardize    
top5_df_S <- as.data.frame(top5_df_S)
str(top5_df_S)

top5_df_S <- tibble::rownames_to_column(top5_df_S, "SCT_snn_res.0.3")

df.long <- pivot_longer(data = top5_df_S, cols = -SCT_snn_res.0.3, names_to = "Gene", values_to = "Abundance")
head(df.long)

df.long$Gene <- factor(df.long$Gene, levels= top5)

ggplot(data = df.long, mapping = aes(x = SCT_snn_res.0.3, y = Gene, fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option="plasma") +
  coord_flip() +
  theme_classic() +
  theme(axis.title=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.ticks.x =element_blank(),
        axis.text.x = element_text(colour = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 10)) 



##Identify cell types, this will change from sample to sample, examples given below

##Provisionally, identify cell-types, this will change from sample to sample.Example given below 

seurat@meta.data <- seurat@meta.data %>% mutate(CellType = case_when(
  SCT_snn_res.0.7 == "0"  ~ "CKM+ Myonuclei", 
  SCT_snn_res.0.7 == "1"  ~ "EGF+ Myonuclei",
  SCT_snn_res.0.7 == "2"  ~ "MYH7B+ Myonuclei",
  SCT_snn_res.0.7 == "3"  ~ "FAPs",  
  SCT_snn_res.0.7 == "4"  ~ "EC",
  SCT_snn_res.0.7 == "5"  ~ "SMC/Pericytes",
  SCT_snn_res.0.7 == "6"  ~ "Immune",
  SCT_snn_res.0.7 == "7"  ~ "MuSC"
)) 


##Save Seurat RDS file
saveRDS(seurat1, file = paste0('seurat_', sample_name, 'decont', '.RDS'))

sessionInfo()

