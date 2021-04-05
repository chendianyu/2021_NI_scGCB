library(Seurat)
library(Matrix)
library(tidyverse)
library(Matrix.utils)
library(pheatmap)
library(multcomp)
library(CountClust)
library(RColorBrewer)
library(rstatix)
library(patchwork)
library(ggbeeswarm)

read_10x_matrix <- function(matrix_dir, project_id){
    barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
    features.path <- paste0(matrix_dir, "features.tsv.gz")
    matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
    mat <- readMM(file = matrix.path)
    feature.names = read.delim(features.path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    rownames(mat) = feature.names$V2
    ## some ensembl gene id correspond to same gene symbol
    ## thus we need to combine these rows 
    mat <- aggregate.Matrix(mat, rownames(mat), FUN="sum")
    seurat_object <- CreateSeuratObject(counts = mat, 
                                        project = project_id, 
                                        min.cells = 3, 
                                        min.features = 200)
    return(seurat_object)
}

gcb_np_batch1 <- read_10x_matrix("data/filtered_feature_bc_matrix_np_batch1/", "gcb_np_batch1")
gcb_np_batch2 <- read_10x_matrix("data/filtered_feature_bc_matrix_np_batch2/", "gcb_np_batch2")

## Integrate two batches
gcb_np_list <- list(gcb_np_batch1, gcb_np_batch2)
cell_prefix_gcb_np <- c("gcb_np_batch1", "gcb_np_batch2")
for (i in 1:length(gcb_np_list)) {
    gcb_np_list[[i]] <- RenameCells(gcb_np_list[[i]], 
                                     add.cell.id = cell_prefix_gcb_np[i])
    gcb_np_list[[i]][["percent.mt"]] <- PercentageFeatureSet(gcb_np_list[[i]], pattern = "^mt-")
    gcb_np_list[[i]] <- subset(gcb_np_list[[i]], 
                                subset = nFeature_RNA > 500 & nCount_RNA < 20000 & percent.mt < 10)
    gcb_np_list[[i]] <- NormalizeData(gcb_np_list[[i]], verbose = FALSE)
    gcb_np_list[[i]] <- FindVariableFeatures(gcb_np_list[[i]], selection.method = "vst", 
                                              nfeatures = 2000, verbose = FALSE)
}

gcb_np_anchors <- FindIntegrationAnchors(object.list = gcb_np_list, dims = 1:30)
gcb_np_integrated <- IntegrateData(anchorset = gcb_np_anchors, 
                                    dims = 1:30)

DefaultAssay(gcb_np_integrated) <- "integrated"
gcb_np_integrated <- ScaleData(gcb_np_integrated)

gcb_np_integrated <- RunPCA(gcb_np_integrated, features = VariableFeatures(object = gcb_np_integrated))
DimPlot(gcb_np_integrated, reduction = "pca")
ElbowPlot(gcb_np_integrated)

gcb_np_integrated <- FindNeighbors(gcb_np_integrated, dims = 1:15)
gcb_np_integrated <- FindClusters(gcb_np_integrated, resolution = 0.2)
gcb_np_integrated <- RunUMAP(gcb_np_integrated, dims = 1:15)
DimPlot(gcb_np_integrated,
        pt.size = 0.01,
        label = T) +
    coord_equal(expand = FALSE)+
    lims(x=c(-7,15), y=c(-7,15)) +
    NoLegend()

DimPlot(gcb_np_integrated,
        group.by = "orig.ident",
        pt.size = 0.01) +
    coord_equal(expand = FALSE)+
    lims(x=c(-7,15), y=c(-7,15)) +
    NoLegend()

## change assay to orignal data
DefaultAssay(gcb_np_integrated) <- "RNA"
gcb_np_integrated <- ScaleData(gcb_np_integrated)
gcb_np_integrated$orig.ident <- factor(gcb_np_integrated$orig.ident)
gcb_np_integrated.markers <- FindAllMarkers(gcb_np_integrated,
                                             test.use = "MAST",
                                             only.pos = TRUE, 
                                             min.pct = 0.25, 
                                             logfc.threshold = 0.25, 
                                             latent.vars = "orig.ident")
top10_gcb_np_integrated <- gcb_np_integrated.markers %>% 
    group_by(cluster) %>% 
    arrange(p_val_adj) %>% 
    slice(1:10)
DoHeatmap(gcb_np_integrated, 
          features = top10_gcb_np_integrated$gene, 
          #group.colors = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", 
          #                 "#3C5488FF", "#7E6148FF", "#9632B8FF"),
          angle = 0,
          draw.lines = FALSE) +
    scale_fill_gradientn(colors = c("#1565C0", "white", "#b92b27")) +
    scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", 
                                  "#3C5488FF", "#7E6148FF", "#9632B8FF"))

## DEG between each pair of clusters
deg_pairwise <- apply(combn(levels(gcb_np_integrated$seurat_clusters), 2), 
                      2, 
                      function(x) FindMarkers(gcb_np_integrated,
                                              ident.1 = x[1],
                                              ident.2 = x[2],
                                              test.use = "MAST",
                                              min.pct = 0.25,
                                              logfc.threshold = 0.25,
                                              latent.vars = "orig.ident"))
names(deg_pairwise) <- apply(combn(levels(gcb_np_integrated$seurat_clusters), 2), 
                             2, 
                             function(x) paste0("cluster_", x[1],"_vs_",x[2]))
lapply(names(deg_pairwise), function(df) write.csv(deg_pairwise[[df]], file=paste0("deg_", df, ".csv")))

########
## vlnplot of marker gene
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
    p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
        geom_quasirandom(size=0.01, alpha = 0.5) +
        xlab("") + ylab(feature) + ggtitle("") + 
        theme(legend.position = "none", 
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.title.y = element_text(size = rel(1), angle = 0), 
              axis.text.y = element_text(size = rel(1)), 
              plot.margin = plot.margin ) 
    return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
    ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
}
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
    
    plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
    
    # Add back x-axis title to bottom plot. patchwork is going to support this?
    plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
        theme(axis.text.x=element_text(), axis.ticks.x = element_line())
    
    # change the y-axis tick to only max value 
    #ymaxs<- purrr::map_dbl(plot_list, extract_max)
    #plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
    #                            scale_y_continuous(breaks = c(y)) + 
    #                            expand_limits(y = y))
    
    p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
    return(p)
}

features <- c("Bcl6", "Aicda", "Irf8", "Bach2", "Irf4", "Cd38", "Pax5", "Prdm1", "Mki67", "Sdc1",
              "Myc", "Mif", "Pa2g4", "Nolc1", "Fbl", "Ldha")
StackedVlnPlot(obj = gcb_np_integrated, 
               features = features,
               idents = 0:5)
