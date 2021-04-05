# topic model
## B cells 
b.cells = 0:5
b.cells.names = rownames(gcb_np_integrated@meta.data)[which(gcb_np_integrated@meta.data$seurat_clusters %in% b.cells)]

type <- "b.cells.names"
## Exclude ribosomal proteins
ribo.loc <- grepl("Rpl|Rps", rownames(gcb_np_integrated[["RNA"]]@counts))

## topic model parameters
n.topics <- 16
tolerance <- 0.5

## create folder to store result
save.path <- 'topic_model/integrate/'
sub.dir <- paste0(n.topics, "topics_tol", tolerance)
dir.create(paste0(save.path, sub.dir), recursive = TRUE)

## fit topic model, saves as .rda file
FitGoM(t(as.matrix(gcb_np_integrated[["RNA"]]@counts)[!ribo.loc, eval(as.symbol(type))]),
       K=n.topics, tol=tolerance,
       path_rda=paste0(save.path, sub.dir, '/FitGoM_k', n.topics,
                       '_tol', tolerance, '.rda'))

# load the topic modeling results 
load(paste0(save.path, sub.dir, '/FitGoM_k', n.topics,
            '_tol', tolerance, '.rda'))

## omega contains topic modeling scores for each cell 
## add these scores to seurat object
omega <- as.data.frame(Topic_clus$omega)
colnames(omega)<-paste0("Topic_", colnames(omega))
gcb_np_integrated <- AddMetaData(gcb_np_integrated, omega)

## list of topic names
topic_id <- colnames(gcb_np_integrated@meta.data[paste0("Topic_", seq(1:n.topics))])
## find the outliers 
outlier <- lapply(topic_id, 
                  function(x) which(gcb_np_integrated@meta.data[,x]>
                                        (mean(gcb_np_integrated@meta.data[,x], na.rm=TRUE)+
                                             3*sd(gcb_np_integrated@meta.data[,x], na.rm=TRUE))))                
names(outlier) <- topic_id
## find the max value of the data when outliers are ignored 
max.nooutlier <- lapply(topic_id, function(x) if(length(unlist(outlier[x])))
{ceiling(10*max(gcb_np_integrated@meta.data[-unlist(outlier[x]),x], na.rm=TRUE))/10}
else {ceiling(10*max(gcb_np_integrated@meta.data[,x], na.rm=TRUE))/10})
names(max.nooutlier) <- topic_id

## Make a Seruat object where outlier values are reset to the max value for non-outlier cells - used only for plotting 
gcb_np_integrated.override <- gcb_np_integrated
for (xx in topic_id) {
    gcb_np_integrated.override@meta.data[unlist(outlier[xx]),xx] = rep(unlist(max.nooutlier[xx]))
}

## generate list of UMAP plots 
gcb_np_integrated.type <- subset(gcb_np_integrated, cells = eval(as.symbol(type)))
gcb_np_integrated.type.override <- subset(gcb_np_integrated.override, cells = eval(as.symbol(type)))
topic_featureplot <- lapply(topic_id, function(x) FeaturePlot(gcb_np_integrated.type.override, 
                                                              x, 
                                                              min.cutoff = 0, 
                                                              pt.size = 2,
                                                              max.cutoff = unlist(max.nooutlier[x])) +
                                scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10)))
topic_featureplot <- lapply(X = topic_featureplot, FUN = AugmentPlot)
CombinePlots(plots = topic_featureplot)

# collect gene scores for each topic - the amount that gene contributes to that topic 
theta <- as.data.frame(Topic_clus$theta)
colnames(theta) <- paste0("Topic_", colnames(theta))

# Run CountClust function ExtractTopFeatures to find the top 200 genes driving each topic
n.features <- 50
features <- ExtractTopFeatures(theta, top_features=n.features)
features.genes <- as.data.frame(sapply(1:n.topics, function(x) rownames(theta)[features$indices[x,]]))
colnames(features.genes) <- paste0("topic_", seq(1:n.topics))
features.scores <- as.data.frame(features$scores)
rownames(features.scores) <- paste0("topic_", seq(1:n.topics))
write.csv(features.genes, 
          paste0(save.path, sub.dir, '/top', n.features, 'features_k', n.topics, '_tol', tolerance, '.csv'), 
          row.names = FALSE)

# Visualize the top 20 genes identified by ExtractTopFeatures on a bar plot 
# These plots do not use log scale, as was used for the main text figures for improved visualization 
plot.features <- data.frame()
for (j in 1:nrow(features.scores)){
    
    current.topic <- as.data.frame(t(features.scores[j,]))
    colnames(current.topic) <- "score"
    rownames(current.topic) <- features.genes[,j]
    current.topic$genes <- features.genes[,j]
    current.topic.top <- current.topic[1:50,]
    current.topic.top$topic <- rep(colnames(features.genes)[j], nrow(current.topic.top))
    current.topic.top$genes <- factor(current.topic.top$genes, levels=rev(current.topic.top$genes))
    current.topic.top$order <- paste0(current.topic.top$genes, "_", current.topic.top$topic)
    current.topic.top$order <- factor(current.topic.top$order, levels=rev(current.topic.top$order))
    
    plot.features <- rbind(plot.features,current.topic.top)
}

# turn levels entry into a factor so that plots stay in numerical rather than alaphabetical order 
plot.features$topic <- factor(plot.features$topic, levels=unique(plot.features$topic))

# generate bar plots  
ggplot(data=plot.features, aes(x=order, y=score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    facet_wrap(~topic, scales="free", ncol=4) + 
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)


########
## Highlight specific genes
## oxidative phosphorylation
plot.features %>% 
    filter(topic == "topic_1",
           genes %in% c("Atp5e", "Atp5g2", "Cox6b1", 
                        "Cox6c", "Cox8a", "Uqcrb",
                        "Ndufa4", "Ndufb5", "Atp5h")) %>% 
    ggplot(aes(order, score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)

## antigen presentation
plot.features %>% 
    filter(topic == "topic_3",
           genes %in% c("Ms4a1", "H2-DMa", "H2-DMb2",
                        "H2-Oa", "H2-Ob", "Cd74",
                        "Cd81", "H2-Eb1", "H2-Ab1", "H2-Eb2")) %>% 
    ggplot(aes(order, score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)

## cell cycle
plot.features %>% 
    filter(topic == "topic_6",
           genes %in% c("Tyms", "Rnaseh2b", "Rbl1",
                        "Dhfr", "Psmc3", "Mcm7",
                        "Pole3", "Mcm10", "Anapc5", "Hsp90aa1")) %>% 
    ggplot(aes(order, score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)

## cell cycle 15
plot.features %>% 
    filter(topic == "topic_6",
           genes %in% c("Tyms", "Rnaseh2b", "Rbl1",
                        "Dhfr", "Psmc3", "Mcm7",
                        "Cenps", "Pole3", "Mcm10",
                        "Kif18b", "Rfc5", "Cenpk",
                        "Anapc5", "Fignl1", "Hsp90aa1")) %>% 
    ggplot(aes(order, score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)

## nuclear division in topic 9
plot.features %>% 
    filter(topic == "topic_9",
           genes %in% c("Fbxo5", "Mad2l1", "Ndc80", "Spc24", "Aurkb",
                        "Cdca5", "Cenpw", "Ccnf", "Cdk1", "Incenp",
                        "Kif11", "Ncaph", "Nusap1", "Ska1", "Smc4")) %>% 
    ggplot(aes(order, score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)

## mitotic sister chromatid segregation in topic 9
plot.features %>% 
    filter(topic == "topic_9",
           genes %in% c("Nusap1", "Fbxo5", "Incenp", "Ndc80",
                        "Smc4", "Ncaph", "Ttk",
                        "Cdca5", "Mad2l1", "Aurkb")) %>% 
    ggplot(aes(order, score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)

## adaptive immune in topic 5
plot.features %>% 
    filter(topic == "topic_5",
           genes %in% c("Fcer2a", "Swap70", "Icam1",
                        "Nfkb2", "Cd40", "Nfkbid",
                        "Bcl3", "Slamf1")) %>% 
    ggplot(aes(order, score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)

## adaptive immune in topic 5
plot.features %>% 
    filter(topic == "topic_5",
           genes %in% c("Mybbp1a", "Ppan", "Gar1",
                        "Nop16", "Fbl", "Rcl1",
                        "Gnl2")) %>% 
    ggplot(aes(order, score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)

## adaptive immune in topic 7
plot.features %>% 
    filter(topic == "topic_7",
           genes %in% c("Calm3", "Tpx2", "Dynll1",
                        "Hmmr", "Skp1a", "Ccnb1")) %>% 
    ggplot(aes(order, score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)

## adaptive immune in topic 6
plot.features %>% 
    filter(topic == "topic_6",
           genes %in% c("Tyms", "Rbl1", "Dhfr",
                        "Mcm7", "Pole3", "Mcm10")) %>% 
    ggplot(aes(order, score)) +
    geom_bar(stat="identity") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=8), text=element_text(size=10), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(breaks=plot.features$order, labels=plot.features$genes)


########
## Connect clustering and topic model
cluster_topic_mtx <- gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters != 6) %>% 
    gather(key = "topic", value = "weight", paste0("Topic_", 1:16), factor_key = T) %>%
    group_by(seurat_clusters, topic) %>% 
    summarise(weight = mean(weight)) %>%
    spread(key = topic, value = weight) 
cluster_topic_mtx <- column_to_rownames(cluster_topic_mtx, var = "seurat_clusters") 
topic_order <- order(max.col(t(apply(cluster_topic_mtx, 2, scale))))
pheatmap(t(cluster_topic_mtx)[c(2, 11, 5, 12, 1, 9, 6, 8, 7, 14, 16, 15, 10, 4, 3, 13), ],
         cluster_rows = F,
         cluster_cols = F,
         scale = "row")

########
## topic weight across cluster
plot_topic_cluster <- function(topic){
    gcb_np_integrated@meta.data %>% 
        filter(seurat_clusters!=6) %>% 
        ggplot(aes_string("seurat_clusters", paste0("Topic_", topic), color = "seurat_clusters")) +
        geom_violin(draw_quantiles=c(0.25,0.75)) +
        geom_quasirandom(size = 0.01, width = 0.3, alpha = 0.5) +
        stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              panel.border = element_rect(fill = NA))
}
## Topic 1
gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters!=6) %>% 
    ggplot(aes(seurat_clusters, Topic_1, color = seurat_clusters)) +
    geom_violin(draw_quantiles=c(0.25,0.75)) +
    geom_quasirandom(size = 0.01, width = 0.3, alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA))
fit <- aov(Topic_1 ~ seurat_clusters, gcb_np_integrated@meta.data %>% filter(seurat_clusters!=6))

## Topic 3
gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters!=6) %>% 
    ggplot(aes(seurat_clusters, Topic_3, color = seurat_clusters)) +
    geom_violin(draw_quantiles=c(0.25,0.75)) +
    geom_quasirandom(size = 0.01, width = 0.15, alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA))
fit <- aov(Topic_3 ~ seurat_clusters, gcb_np_integrated@meta.data %>% filter(seurat_clusters!=6))

## Topic 6
gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters!=6) %>% 
    ggplot(aes(seurat_clusters, Topic_6, color = seurat_clusters)) +
    geom_violin(draw_quantiles=c(0.25,0.75)) +
    geom_quasirandom(size = 0.01, width = 0.15, alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA))
fit <- aov(Topic_6 ~ seurat_clusters, gcb_np_integrated@meta.data %>% filter(seurat_clusters!=6))

## Topic 9
gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters!=6) %>% 
    ggplot(aes(seurat_clusters, Topic_9, color = seurat_clusters)) +
    geom_violin(draw_quantiles=c(0.25,0.75)) +
    geom_quasirandom(size = 0.01, width = 0.15, alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA))
fit <- aov(Topic_9 ~ seurat_clusters, gcb_np_integrated@meta.data %>% filter(seurat_clusters!=6))

## Topic 5
gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters!=6) %>% 
    ggplot(aes(seurat_clusters, Topic_5, color = seurat_clusters)) +
    geom_violin(draw_quantiles=c(0.25,0.75)) +
    geom_quasirandom(size = 0.01, width = 0.15, alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA))
fit <- aov(Topic_5 ~ seurat_clusters, gcb_np_integrated@meta.data %>% filter(seurat_clusters!=6))

## Topic 7
gcb_np_integrated@meta.data %>% 
    filter(seurat_clusters!=6) %>% 
    ggplot(aes(seurat_clusters, Topic_7, color = seurat_clusters)) +
    geom_violin(draw_quantiles=c(0.25,0.75)) +
    geom_quasirandom(size = 0.01, width = 0.15, alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA))
fit <- aov(Topic_7 ~ seurat_clusters, gcb_np_integrated@meta.data %>% filter(seurat_clusters!=6))


