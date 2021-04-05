library(rstatix)

## topic weight by affinity
gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters != 5,
           affinity %in% c("Low_affinity", "High_affinity")) %>% 
    gather(key = "topic", value = "weight", paste0("Topic_", 1:16)) %>% 
    mutate(topic=factor(topic, levels = paste0("Topic_", 1:16)),
           affinity=factor(affinity, levels = c("Low_affinity", "High_affinity"))) %>% 
    ggplot(aes(affinity, weight, color=affinity)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~topic, scales = "free") +
    scale_color_manual(values = c("Low_affinity" = "#3C5488FF", 
                                  "High_affinity" = "#E64B35FF")) +
    theme(panel.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(colour = "black"),
          axis.ticks.x = element_line(size = 1))

topic_affinity_df <- gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters != 5,
           affinity %in% c("Low_affinity", "High_affinity"))
pvalue_topic_affinity <- list()
for (i in 1:16) {
    pvalue_topic_affinity[[i]] <- t.test(as.formula(paste(paste0("Topic_", i), "~ affinity")), topic_affinity_df)
}

########
## Considering other factors that may influence the signature score
## topic weight by cluster and affinity
gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters != 6,
           affinity %in% c("Low_affinity", "High_affinity")) %>% 
    gather(key = "topic", value = "weight", paste0("Topic_", 1:16)) %>% 
    mutate(topic=factor(topic, levels = paste0("Topic_", 1:16)),
           affinity=factor(affinity, levels = c("Low_affinity", "High_affinity"))) %>% 
    ggplot(aes(seurat_clusters, weight, color=affinity, group = interaction(seurat_clusters, affinity))) +
    geom_violin(draw_quantiles=c(0.25,0.75), trim = FALSE, width = 1) +
    geom_quasirandom(dodge.width = 1, size=0.01, width = 0.1, alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=2, position = position_dodge(width = 1), color = "black") +
    facet_wrap(~topic, scales = "free", ncol = 2) +
    scale_color_manual(values = c("Low_affinity" = "#3C5488FF", 
                                  "High_affinity" = "#E64B35FF")) +
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.ticks.x = element_line(size = 1))

topic_cluster_affinity_df <- gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters != 6,
           affinity %in% c("Low_affinity", "High_affinity")) %>% 
    mutate(affinity=factor(affinity, levels = c("Low_affinity", "High_affinity")))
pvalue_topic_cluster_affinity <- list()
for (i in 1:16) {
    topic_fit <- aov(as.formula(paste(paste0("Topic_", i), "~ seurat_clusters*affinity")), topic_cluster_affinity_df)
    pvalue_topic_cluster_affinity[[i]] <- summary(topic_fit)
}
topic_cluster_affinity_df %>% 
    group_by(seurat_clusters) %>% 
    pairwise_wilcox_test(Topic_1~affinity, p.adjust.method = "BH") %>% 
    adjust_pvalue(method = "BH")


## topic weight by mutation number and affinity
metadata_all %>% 
    filter(seurat_clusters != 6,
           affinity %in% c("Low_affinity", "High_affinity")) %>% 
    gather(key = "topic", value = "weight", paste0("Topic_", 1:16)) %>% 
    mutate(topic=factor(topic, levels = paste0("Topic_", 1:16)),
           affinity=factor(affinity, levels = c("Low_affinity", "High_affinity")),
           Total_diff=factor(Total_diff)) %>% 
    ggplot(aes(Total_diff, weight, color=affinity, group = interaction(Total_diff, affinity))) +
    geom_violin(draw_quantiles=c(0.25,0.75), trim = FALSE, width = 1) +
    geom_quasirandom(dodge.width = 1, size=0.01, width = 0.1, alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=2, position = position_dodge(width = 1), color = "black") +
    facet_wrap(~topic, scales = "free", ncol = 1) +
    scale_color_manual(values = c("Low_affinity" = "#3C5488FF", 
                                  "High_affinity" = "#E64B35FF")) +
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.ticks.x = element_line(size = 1))

topic_mutation_affinity_df <- metadata_all %>% 
    filter(seurat_clusters != 6,
           affinity %in% c("Low_affinity", "High_affinity")) %>% 
    mutate(affinity=factor(affinity, levels = c("Low_affinity", "High_affinity")),
           Total_diff=factor(Total_diff))
pvalue_topic_mutation_affinity <- list()
for (i in 1:16) {
    topic_fit <- aov(as.formula(paste(paste0("Topic_", i), "~ Total_diff*affinity")), topic_mutation_affinity_df)
    pvalue_topic_mutation_affinity[[i]] <- summary(topic_fit)
}

## topic weight by cluster and isotype
metadata_all %>% 
    filter(seurat_clusters != 5,
           affinity %in% c("Low_affinity", "High_affinity")) %>% 
    gather(key = "topic", value = "weight", paste0("Topic_", 1:16)) %>% 
    mutate(topic=factor(topic, levels = paste0("Topic_", 1:16)),
           affinity=factor(affinity, levels = c("Low_affinity", "High_affinity")),
           c_gene=factor(c_gene)) %>% 
    ggplot(aes(c_gene, weight, color=affinity, group = interaction(c_gene, affinity))) +
    geom_violin(draw_quantiles=c(0.25,0.75), trim = FALSE, width = 1) +
    geom_quasirandom(dodge.width = 1, size=0.01, width = 0.1, alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=2, position = position_dodge(width = 1), color = "black") +
    facet_wrap(~topic, scales = "free", ncol = 2) +
    scale_color_manual(values = c("Low_affinity" = "#3C5488FF", 
                                  "High_affinity" = "#E64B35FF")) +
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.ticks.x = element_line(size = 1))

topic_isotype_affinity_df <- metadata_all %>% 
    filter(seurat_clusters != 5,
           affinity %in% c("Low_affinity", "High_affinity")) %>% 
    mutate(affinity=factor(affinity, levels = c("Low_affinity", "High_affinity")),
           c_gene=factor(c_gene))
pvalue_topic_isotype_affinity <- list()
for (i in 1:16) {
    topic_fit <- aov(as.formula(paste(paste0("Topic_", i), "~ c_gene*affinity")), topic_isotype_affinity_df)
    pvalue_topic_isotype_affinity[[i]] <- summary(topic_fit)
}
