library(monocle)

cds <- as.CellDataSet(subset(gcb_np_integrated, seurat_clusters != 6))
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

clustering_DEG_genes <- differentialGeneTest(cds,
                                             fullModelFormulaStr = '~seurat_clusters',
                                             cores = 1)
cds_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
cds  <- setOrderingFilter(cds, ordering_genes = cds_ordering_genes)
cds  <- reduceDimension(cds, method = 'DDRTree')
cds  <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "seurat_clusters")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")

## pseudotime by cluster
pData(cds) %>% 
    mutate(seurat_clusters = factor(seurat_clusters, levels = c(5, 0, 1, 4, 2, 3))) %>% 
    ggplot(aes(seurat_clusters, Pseudotime, color = seurat_clusters)) +
    geom_violin(draw_quantiles=c(0.25,0.75), scale = "width") +
    geom_quasirandom(dodge.width = 1, size=0.01, alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=2, color = "black") +
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.ticks.x = element_line(size = 1))

fit <- aov(Pseudotime~seurat_clusters, pData(cds))
TukeyHSD(fit)

## gene expression by trajectory
plot_cell_trajectory(cds,
                     markers = c("Bcl6", "Cxcr4", "Cd83", "Cd86",
                                 "Cox6b1", "Ndufa4", "Cox8a", "Uqcrb",
                                 "Atp5g2", "Ndufb5", "Cox6c", "Atp5e",
                                 "Atp5h"),
                     use_color_gradient = TRUE) +
    scale_color_gradient(low = "white", 
                         high = "blue")

## gene expression by pseudotime
my_genes <- row.names(subset(fData(cds),
                             gene_short_name %in% c("Mki67")))
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, 
                         color_by = "seurat_clusters")

## DEG by pseudotime
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
pseudotime_heatmap <- plot_pseudotime_heatmap(cds[sig_gene_names,],
                                              num_clusters = 5,
                                              cores = 1,
                                              show_rownames = F,
                                              return_heatmap = T)
gene_by_cluster <- as.data.frame(cutree(pseudotime_heatmap$tree_row, k=5))
colnames(gene_by_cluster) <- "Cluster"
gene_by_cluster$Gene <- rownames(gene_by_cluster)

gene_by_cluster$Cluster <- factor(gene_by_cluster$Cluster, levels = c(5, 4, 1, 2, 3))
gene_by_cluster <- gene_by_cluster[order(gene_by_cluster$Cluster), ]
plot_pseudotime_heatmap(cds[gene_by_cluster$Gene,],
                        cluster_rows = F,
                        cores = 1,
                        show_rownames = F)
