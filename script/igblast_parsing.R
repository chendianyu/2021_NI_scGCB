## parsing igblast_vh186.2 result
igblast_batch1 <- read.delim("data/gcb_np_batch1.tsv", stringsAsFactors = FALSE)
igblast_batch1$sequence_id <- paste0("gcb_np_batch1_", igblast_batch1$sequence_id)
igblast_batch2 <- read.delim("data/gcb_np_batch2.tsv", stringsAsFactors = FALSE)
igblast_batch2$sequence_id <- paste0("gcb_np_batch2_", igblast_batch2$sequence_id)
igblast <- bind_rows(igblast_batch1, igblast_batch2) 

## vgene mutation position offset
igblast_vh186.2 <- igblast[igblast$v_call=="IGHV1-72*01", ]
vgene_mutation_offset <- (igblast_vh186.2$v_germline_start + 1) %/% 3

## compare each amino acid in sequence
vgene_aa_match <- sapply(1:nrow(igblast_vh186.2), function(x) {
    apply(do.call(rbind, strsplit(c(igblast_vh186.2[x,]$v_sequence_alignment_aa, igblast_vh186.2[x,]$v_germline_alignment_aa), "")), 
          2, 
          function(i){
              i[1] == i[2]
          }
    )
}
)
vgene_mutation_pos <- lapply(vgene_aa_match, function(x) which(!x))
## add offset for position
vgene_mutation_pos <- lapply(1:nrow(igblast_vh186.2), function(i) vgene_mutation_pos[[i]] + vgene_mutation_offset[i])

mutation_pos <- data.frame(pos=unlist(vgene_mutation_pos))
             
ggplot(mutation_pos, aes(pos)) +
    geom_bar(aes(y = (..count..)/length(vgene_mutation_pos)), fill = "#3C548899") + 
    scale_x_continuous(limits = c(0,99),
                       breaks = seq(0,100,1),
                       labels = c(0, rep("",9), 10, rep("",9), 20, rep("",9), 30,
                                  rep("",9), 40, rep("",9), 50, rep("",9), 60,
                                  rep("",9), 70, rep("",9), 80, rep("",9), 
                                  90, rep("",9), 100)) +
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.ticks.x = element_line(size = 1), 
          axis.ticks.length=unit(.25, "cm")) +
    coord_cartesian(expand = FALSE)

as.data.frame(table(str_split(igblast_vh186.2$v_sequence_alignment_aa, "", simplify = T)[, 33])) %>%
    arrange(Freq) %>% 
    mutate(Var1 = factor(Var1, levels = Var1)) %>% 
    mutate(fraction = Freq/sum(Freq),
           ymax = cumsum(fraction),
           ymin = c(0, head(ymax, n=-1))) %>%
    ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
    geom_rect(color = "black") +
    coord_polar(theta="y") + 
    xlim(c(2, 4)) +
    theme_void() +
    scale_fill_manual(values = rev(brewer.pal(12, "Paired")[1:11]))





