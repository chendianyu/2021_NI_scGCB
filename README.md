# 2021_NI_scGCB
This repository contains some scripts used for our Nature Immunology paper entitled ["Coupled analysis of transcriptome and BCR mutations reveals role of OXPHOS in affinity maturation"](https://www.nature.com/articles/s41590-021-00936-y), including:  
* **script**:  
    * `GCB_NP.Rmd`: Classical analysis of scRNA-seq data of NP-specific GC B cells  
    * `Signature_score.Rmd`: Signature score calculation based on gene list  
    * `Topic_model.Rmd`: Topic modeling analysis of scRNA-seq data  
    * `Igblast_parsing.Rmd`: Detailed code for parsing output of CellRanger VDJ and IgBlast  
    * `BCR.Rmd`: Combine BCR (affinity) info with scRNA-seq to identify difference between high- and low-affinity groups  
    * `monocle2.R`: Pseudotime analysis using Monocle 2

* **html**: Corresponding html output of `.Rmd` file in directory `script` 
