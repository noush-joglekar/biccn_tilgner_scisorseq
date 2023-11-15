#!/bin/R
# Code for reproducing panels from figure 5
## Relating to the similarity in exon and short-read
## gene expression between glial cell types

## Set up -----
library(dplyr)
library(tidyr)
library(tibble)
library(superheat)
library(viridis)
library(plotly)
library(ggsignif)
library(MetBrewer)
library(ggforce)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Seurat)
library(slingshot)
library(harmony)

## Fig 5a preprocessing -------
subtype_alt <- read.table('../data/altExonUsage_devel_subtype.gz',
                          sep = "\t", header = TRUE)

cP <- read.csv('../data/colorPalette.csv',header = FALSE)
colnames(cP) <- c("color","Subtype")

my_pal <- cP$color
names(my_pal) <- cP$Subtype

subtype_glia <- subtype_alt %>% select(c(Exon,Gene,contains(c("Astro","_OPC","COP","MFOL","MOL")))) %>% 
    unite(GE,c("Exon","Gene"),sep = "::") %>% column_to_rownames("GE")

threshold <- 0.95
subtype_glia2 <- subtype_glia %>%
    select(which(colMeans(is.na(.)) < 0.95)) %>%
    filter(rowMeans(is.na(.)) <= threshold)


C_st_G = cor(subtype_glia2,use = "complete",method = "spearman")

colnames(C_st_G) <- rownames(C_st_G)

colnames(C_st_G) <- gsub("VisCortex","VIS",gsub("Hippocampus","HIPP",colnames(C_st_G)))
rownames(C_st_G) <- colnames(C_st_G)

cnn <- data.frame("CT" = rownames(C_st_G)) %>% 
    separate(CT,into = c("Age","Region","Type"))

regions <- c("VIS","HIPP")
reg_cols <- c("#ddbbc6","#fdc681")

cts <- names(my_pal)
ct_cols <- my_pal

ages <- c("P14","P21","P28","P56")
age_cols <- c("#ae3b24","#c8812a","#56632c","#f68d62")

ha = rowAnnotation(
  Celltype = cnn$Type,
  Age = cnn$Age,
  Region = cnn$Region,
  
  col = list(Celltype = structure(ct_cols,names = cts),
             Age = structure(age_cols, names = ages),
             Region = structure(reg_cols,names = regions)),
  
  annotation_legend_param = list(
    Celltype = list(
      title = "cellType",
      at = cts,
      labels = cts
    ),
    Age = list(
      title = "Age",
      at = ages,
      labels = ages
    ),
    Region = list(
      title = "Region",
      at = regions,
      labels = regions
    )
  )
)


col_func <- colorRamp2(seq(0.8,1,length.out = 7),brewer.pal(7,"BrBG"))

options(repr.plot.width=15, repr.plot.height=11)

## Fig 5a ------
Heatmap(C_st_G, col = col_func, show_column_names = FALSE,
        show_column_dend = FALSE,
        left_annotation = ha)

## Fig 5c preprocessing -----
lrBC_glia <- readRDS('R_objects/lrBC_devel_glia_baseSeurat.rds') ## file too large to deposit. Please get in touch if you want it
lrBC_glia

set.seed(10)
subsetObj <- subset(lrBC_glia, subset = Subtype != "DivOPCs")
subsetObj <- subset(subsetObj, subset = seurat_clusters != 10)
dim(subsetObj)
subsetCells <- colnames(subsetObj)[sample(45114,10000,replace = FALSE)]
subsetObj <- subset(lrBC_glia,cells = subsetCells)

subsetObj <- RunUMAP(subsetObj,dims = 1:5)

sds <- slingshot(Embeddings(subsetObj, "umap"),start.clus = "OPCs",
                 clusterLabels = subsetObj$Subtype)

options(repr.plot.width=6, repr.plot.height=6)
pal_use <- my_pal[subsetObj$Subtype]

## Fig 5c ------
plot(reducedDim(sds), col = pal_use, pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')

sessionInfo()
# R version 4.0.1 (2020-06-06)
# Platform: x86_64-conda_cos6-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /pbtech_mounts/homes059/anj2026/miniconda3/envs/RV4_Seurat323/lib/libopenblasp-r0.3.9.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#   [1] harmony_0.1.0         Rcpp_1.0.8.3          slingshot_1.6.1      
# [4] princurve_2.1.6       Seurat_3.2.3          RColorBrewer_1.1-3   
# [7] circlize_0.4.15       ComplexHeatmap_2.13.1 patchwork_1.1.1      
# [10] ggforce_0.3.3         MetBrewer_0.2.0       ggsignif_0.6.3       
# [13] plotly_4.10.0         ggplot2_3.3.5         viridis_0.6.2        
# [16] viridisLite_0.4.0     superheat_1.0.0       tibble_3.1.8         
# [19] tidyr_1.2.0           dplyr_1.0.7          
# 
# loaded via a namespace (and not attached):
#   [1] uuid_0.1-4                  plyr_1.8.6                 
# [3] igraph_1.5.1                repr_1.1.3                 
# [5] lazyeval_0.2.2              splines_4.0.1              
# [7] listenv_0.8.0               scattermore_0.7            
# [9] GenomeInfoDb_1.24.2         digest_0.6.29              
# [11] foreach_1.5.2               htmltools_0.5.3            
# [13] fansi_0.5.0                 magrittr_2.0.3             
# [15] tensor_1.5                  cluster_2.1.0              
# [17] doParallel_1.0.17           ROCR_1.0-11                
# [19] globals_0.14.0              matrixStats_0.62.0         
# [21] colorspace_2.0-3            ggrepel_0.9.1              
# [23] crayon_1.5.1                RCurl_1.98-1.2             
# [25] jsonlite_1.8.7              spatstat_1.64-1            
# [27] spatstat.data_1.7-0         ape_5.4-1                  
# [29] survival_3.1-12             zoo_1.8-8                  
# [31] iterators_1.0.14            glue_1.6.2.9000            
# [33] polyclip_1.10-0             gtable_0.3.0               
# [35] zlibbioc_1.36.0             XVector_0.28.0             
# [37] leiden_0.3.7                DelayedArray_0.14.1        
# [39] GetoptLong_1.0.5            future.apply_1.7.0         
# [41] shape_1.4.6                 SingleCellExperiment_1.10.1
# [43] BiocGenerics_0.36.1         abind_1.4-5                
# [45] scales_1.2.1                DBI_1.1.3                  
# [47] miniUI_0.1.1.1              xtable_1.8-4               
# [49] clue_0.3-61                 reticulate_1.24            
# [51] rsvd_1.0.3                  stats4_4.0.1               
# [53] htmlwidgets_1.5.3           httr_1.4.2                 
# [55] ellipsis_0.3.2              ica_1.0-2                  
# [57] pkgconfig_2.0.3             farver_2.1.0               
# [59] uwot_0.1.10                 deldir_0.2-9               
# [61] utf8_1.2.2                  tidyselect_1.1.1           
# [63] rlang_1.0.2                 reshape2_1.4.4             
# [65] later_1.1.0.1               munsell_0.5.0              
# [67] tools_4.0.1                 cli_3.3.0                  
# [69] generics_0.1.2              ggridges_0.5.3             
# [71] evaluate_0.22               stringr_1.4.0              
# [73] fastmap_1.1.0               goftest_1.2-2              
# [75] fitdistrplus_1.1-3          purrr_0.3.4                
# [77] RANN_2.6.1                  pbapply_1.5-0              
# [79] future_1.21.0               nlme_3.1-148               
# [81] mime_0.9                    compiler_4.0.1             
# [83] png_0.1-7                   spatstat.utils_2.0-0       
# [85] tweenr_1.0.2                stringi_1.7.8              
# [87] lattice_0.20-41             IRdisplay_1.0              
# [89] Matrix_1.5-3                vctrs_0.4.1                
# [91] pillar_1.8.1                lifecycle_1.0.1            
# [93] lmtest_0.9-38               GlobalOptions_0.1.2        
# [95] RcppAnnoy_0.0.18            data.table_1.13.6          
# [97] cowplot_1.1.1               bitops_1.0-6               
# [99] irlba_2.3.3                 httpuv_1.5.5               
# [101] GenomicRanges_1.40.0        R6_2.5.1                   
# [103] promises_1.1.1              KernSmooth_2.23-17         
# [105] gridExtra_2.3               IRanges_2.24.1             
# [107] parallelly_1.23.0           codetools_0.2-16           
# [109] MASS_7.3-51.6               assertthat_0.2.1           
# [111] SummarizedExperiment_1.18.2 rjson_0.2.21               
# [113] withr_2.5.0                 sctransform_0.3.2          
# [115] S4Vectors_0.28.1            GenomeInfoDbData_1.2.4     
# [117] mgcv_1.8-31                 parallel_4.0.1             
# [119] rpart_4.1-15                IRkernel_1.2               
# [121] Cairo_1.5-12.2              Rtsne_0.15                 
# [123] pbdZMQ_0.3-4                Biobase_2.50.0             
# [125] shiny_1.5.0                 base64enc_0.1-3   R version 4.0.1 (2020-06-06)
# Platform: x86_64-conda_cos6-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /pbtech_mounts/homes059/anj2026/miniconda3/envs/RV4_Seurat323/lib/libopenblasp-r0.3.9.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#   [1] harmony_0.1.0         Rcpp_1.0.8.3          slingshot_1.6.1      
# [4] princurve_2.1.6       Seurat_3.2.3          RColorBrewer_1.1-3   
# [7] circlize_0.4.15       ComplexHeatmap_2.13.1 patchwork_1.1.1      
# [10] ggforce_0.3.3         MetBrewer_0.2.0       ggsignif_0.6.3       
# [13] plotly_4.10.0         ggplot2_3.3.5         viridis_0.6.2        
# [16] viridisLite_0.4.0     superheat_1.0.0       tibble_3.1.8         
# [19] tidyr_1.2.0           dplyr_1.0.7          
# 
# loaded via a namespace (and not attached):
#   [1] uuid_0.1-4                  plyr_1.8.6                 
# [3] igraph_1.5.1                repr_1.1.3                 
# [5] lazyeval_0.2.2              splines_4.0.1              
# [7] listenv_0.8.0               scattermore_0.7            
# [9] GenomeInfoDb_1.24.2         digest_0.6.29              
# [11] foreach_1.5.2               htmltools_0.5.3            
# [13] fansi_0.5.0                 magrittr_2.0.3             
# [15] tensor_1.5                  cluster_2.1.0              
# [17] doParallel_1.0.17           ROCR_1.0-11                
# [19] globals_0.14.0              matrixStats_0.62.0         
# [21] colorspace_2.0-3            ggrepel_0.9.1              
# [23] crayon_1.5.1                RCurl_1.98-1.2             
# [25] jsonlite_1.8.7              spatstat_1.64-1            
# [27] spatstat.data_1.7-0         ape_5.4-1                  
# [29] survival_3.1-12             zoo_1.8-8                  
# [31] iterators_1.0.14            glue_1.6.2.9000            
# [33] polyclip_1.10-0             gtable_0.3.0               
# [35] zlibbioc_1.36.0             XVector_0.28.0             
# [37] leiden_0.3.7                DelayedArray_0.14.1        
# [39] GetoptLong_1.0.5            future.apply_1.7.0         
# [41] shape_1.4.6                 SingleCellExperiment_1.10.1
# [43] BiocGenerics_0.36.1         abind_1.4-5                
# [45] scales_1.2.1                DBI_1.1.3                  
# [47] miniUI_0.1.1.1              xtable_1.8-4               
# [49] clue_0.3-61                 reticulate_1.24            
# [51] rsvd_1.0.3                  stats4_4.0.1               
# [53] htmlwidgets_1.5.3           httr_1.4.2                 
# [55] ellipsis_0.3.2              ica_1.0-2                  
# [57] pkgconfig_2.0.3             farver_2.1.0               
# [59] uwot_0.1.10                 deldir_0.2-9               
# [61] utf8_1.2.2                  tidyselect_1.1.1           
# [63] rlang_1.0.2                 reshape2_1.4.4             
# [65] later_1.1.0.1               munsell_0.5.0              
# [67] tools_4.0.1                 cli_3.3.0                  
# [69] generics_0.1.2              ggridges_0.5.3             
# [71] evaluate_0.22               stringr_1.4.0              
# [73] fastmap_1.1.0               goftest_1.2-2              
# [75] fitdistrplus_1.1-3          purrr_0.3.4                
# [77] RANN_2.6.1                  pbapply_1.5-0              
# [79] future_1.21.0               nlme_3.1-148               
# [81] mime_0.9                    compiler_4.0.1             
# [83] png_0.1-7                   spatstat.utils_2.0-0       
# [85] tweenr_1.0.2                stringi_1.7.8              
# [87] lattice_0.20-41             IRdisplay_1.0              
# [89] Matrix_1.5-3                vctrs_0.4.1                
# [91] pillar_1.8.1                lifecycle_1.0.1            
# [93] lmtest_0.9-38               GlobalOptions_0.1.2        
# [95] RcppAnnoy_0.0.18            data.table_1.13.6          
# [97] cowplot_1.1.1               bitops_1.0-6               
# [99] irlba_2.3.3                 httpuv_1.5.5               
# [101] GenomicRanges_1.40.0        R6_2.5.1                   
# [103] promises_1.1.1              KernSmooth_2.23-17         
# [105] gridExtra_2.3               IRanges_2.24.1             
# [107] parallelly_1.23.0           codetools_0.2-16           
# [109] MASS_7.3-51.6               assertthat_0.2.1           
# [111] SummarizedExperiment_1.18.2 rjson_0.2.21               
# [113] withr_2.5.0                 sctransform_0.3.2          
# [115] S4Vectors_0.28.1            GenomeInfoDbData_1.2.4     
# [117] mgcv_1.8-31                 parallel_4.0.1             
# [119] rpart_4.1-15                IRkernel_1.2               
# [121] Cairo_1.5-12.2              Rtsne_0.15                 
# [123] pbdZMQ_0.3-4                Biobase_2.50.0             
# [125] shiny_1.5.0                 base64enc_0.1-3   

## Short read Fig 5b
load('seuratObjects/harmony_hippDev.Robj')
load('seuratObjects/harmony_visDev.Robj')

aE_Hipp <- AverageExpression(hipp)
aE_Vis <- AverageExpression(vis)

glia_hipp <- colnames(aE_Hipp$RNA)
aE_Hipp_all_glia <- aE_Hipp$RNA[,glia_hipp[grepl("Astro|COP|OPC|MOL|MFOL",glia_hipp)]]
colnames(aE_Hipp_all_glia) <- paste("Hipp",colnames(aE_Hipp_all_glia),sep = "_")

glia_vis <- colnames(aE_Vis$RNA)
aE_Vis_all_glia <- aE_Vis$RNA[,glia_vis[grepl("Astro|COP|OPC|MOL|MFOL",glia_vis)]]
colnames(aE_Vis_all_glia) <- paste("Vis",colnames(aE_Vis_all_glia),sep = "_")

aE_all_glia <- merge(aE_Hipp_all_glia,aE_Vis_all_glia,by = 0) %>% column_to_rownames('Row.names')
gliaNames <- colnames(aE_all_glia)

## to remove:
dfCountsH <- as.data.frame(table(hipp$Subtype,hipp$Age))
colnames(dfCountsH)[1:2] <- c("Subtype","Age")
lowDFCountsHG <- dfCountsH %>% filter(grepl("Astro|COP|OPC|MOL|MFOL",Subtype )) %>%
  filter(Freq <= 50)

dfCountsV <- as.data.frame(table(vis$Subtype,vis$Age))
colnames(dfCountsV)[1:2] <- c("Subtype","Age")
lowDFCountsVG <- dfCountsV %>% filter(grepl("Astro|COP|OPC|MOL|MFOL",Subtype )) %>%
  filter(Freq <= 50)

gliaNames <- gliaNames[!grepl("Div|Vis_MOLs_P14|Vis_MFOLs_P56",gliaNames)]

aE_all_glia <- aE_all_glia[,gliaNames]

cg <- data.frame(CT = colnames(aE_all_glia))
cg <- cg %>% separate(CT, into = c("Region","Subtype","Age"))

colPal <- read.csv('../colorPalette.csv', header = FALSE)
colnames(colPal) <- c("Color","Subtype")

subtypes <- colPal %>% filter(grepl("Astro|COP|OPC|MOL|MFOL",Subtype ) & Subtype != "DivOPCs") %>% .$Subtype
subtype_cols <- colPal %>% filter(grepl("Astro|COP|OPC|MOL|MFOL",Subtype ) & Subtype != "DivOPCs") %>% .$Color
  
regions <- c("Vis","Hipp")
reg_cols <- c("#ddbbc6","#fdc681")

ages <- c("P14","P21","P28","P56")
age_cols <- c("#ae3b24","#c8812a","#56632c","#f68d62")

ha = rowAnnotation(
  Subtype = cg$Subtype,
  Age = cg$Age,
  Region = cg$Region,
  
  col = list(Subtype = structure(subtype_cols,names = subtypes),
             Age = structure(age_cols, names = ages),
             Region = structure(reg_cols,names = regions)),
  
  annotation_legend_param = list(
    Subtype = list(
      title = "Subtype",
      at = subtypes,
      labels = subtypes
    ),
    Age = list(
      title = "Age",
      at = ages,
      labels = ages
    ),
    Region = list(
      title = "Region",
      at = regions,
      labels = regions
    )
  )
)




c_g_all <- cor(aE_all_glia,use = "complete",method = "spearman")
superheat(c_g_all,
          bottom.label.text.angle = 90,
          bottom.label.text.size = 3,
          left.label.text.size = 3,
          col.dendrogram = T,
          pretty.order.rows = T,
          heat.pal = viridis::magma(100))


col_func <- colorRamp2(seq(0.8,1,length.out = 7),
                       brewer.pal(7,"BrBG"))

## Fig 5b --------
Heatmap(c_g_all,col = col_func,left_annotation = ha)
