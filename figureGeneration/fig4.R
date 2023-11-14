#!/bin/R
# Code for reproducing panels from figure 4.
# Querying changes in neuronal exon inclusion
# within HIPP and VIS during development

## Setup-----
library(dplyr)
library(tidyr)
library(tibble)
library(viridis)
library(MetBrewer)
library(ggforce)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


## Read in cell-type PSI values ----
type_alt <- read.table('../data/altExonUsage_devel_type.gz',sep = "\t", header = TRUE)

## Pre-processing ----
type_neurons <- type_alt %>% select(c(Exon,Gene,contains("Neuron"))) #%>% 
dim(type_neurons)

type_neurons <- type_neurons %>% unite(GE,c("Exon","Gene"),sep = "::") %>% column_to_rownames("GE")

## Get correlation
CN = cor(type_neurons,use = "complete",method = "spearman")

a <- gsub("Neuron","",rownames(CN))
b <- gsub("Hippocampus","HIPP",a)
rownames(CN) <- gsub("VisCortex","VIS",b)
colnames(CN) <- rownames(CN)

cnn <- data.frame("CT" = rownames(CN)) %>% 
    separate(CT,into = c("Age","Region","Type"))

regions <- c("VIS","HIPP")
reg_cols <- c("#ddbbc6","#fdc681")

cts <- c("Excite","Inhib")
ct_cols <- c("#ffcc5c","#65c3ba")

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

#col_func <- colorRamp2(seq(0.8,1,length.out = 7),brewer.pal(7,"BrBG"))
col_func <- colorRamp2(seq(0.85,1,length.out = 7),viridis::magma(7))


options(repr.plot.width=8, repr.plot.height=6)

## Fig 4a ----
Heatmap(CN, col = col_func, show_column_names = FALSE,
        show_column_dend = FALSE,
        left_annotation = ha)

## separate into Excite and Inhib
tmp_mat <- CN %>% as.data.frame() %>% select(contains("Excite"))
tmp_mat2 <- tmp_mat[grep("Excite",rownames(tmp_mat)),]
CE <- tmp_mat2[upper.tri(tmp_mat2)]

tmp_mat <- CN %>% as.data.frame() %>% select(contains("Inhib"))
tmp_mat2 <- tmp_mat[grep("Inhib",rownames(tmp_mat)),]
CI <- tmp_mat2[upper.tri(tmp_mat2)]


df <- data.frame("CorValue" = c(CE,CI),"Type" = rep(c("Excite","Inhib"),each = length(CE)) )

## Plot boxplot
options(repr.plot.width=6, repr.plot.height=6)
g0 = ggplot(df, aes(x = Type, y = CorValue, fill = Type)) +
    geom_boxplot() + 
    scale_fill_manual(values = c("#ffcc5c","#65c3ba")) + 
    theme_classic(base_size = 20) +
    theme(legend.position = "bottom")

g0


## Read in subtype PSI values -----
subtype_alt <- read.table('../data/altExonUsage_devel_subtype,gz',
                          sep = "\t", header = TRUE)

subtype_neurons <- subtype_alt %>% select(c(Exon,Gene,contains(c("Inh","Excite","Granule","NIPC")))) %>% 
    unite(GE,c("Exon","Gene"),sep = "::") %>% column_to_rownames("GE")

threshold <- 0.95
subtype_neurons2 <- subtype_neurons %>%
    select(which(colMeans(is.na(.)) < 0.95)) %>%
    filter(rowMeans(is.na(.)) <= threshold)

## get correlation ----
C_st_N = cor(subtype_neurons2,use = "complete",method = "spearman")

## separate out excite, inhib, and cajal retzius
tmp_mat_s <- C_st_N %>% as.data.frame() %>% select(contains(c("Excite","Granule")))
tmp_mat_s2 <- tmp_mat_s[grep("Excite|Granule",rownames(tmp_mat_s)),]
CE <- tmp_mat_s2[upper.tri(tmp_mat_s2)]


tmp_mat <- C_st_N %>% as.data.frame() %>% select(contains("Inhib"))
tmp_mat2 <- tmp_mat[grep("Inhib",rownames(tmp_mat)),]
CI <- tmp_mat2[upper.tri(tmp_mat2)]

tmp_mat <- C_st_N %>% as.data.frame() %>% select(contains("Inh"))
tmp_mat2 <- tmp_mat[grep("Inh",rownames(tmp_mat)),]
CR <- tmp_mat2[upper.tri(tmp_mat2)]

df <- data.frame("CorValue" = c(CE,CI,CR),
                 "Type" = c(rep("Excite",length(CE)), 
                            rep("Inhib", length(CI)),
                           rep("Inh_wCR",length(CR))))

df$Type <- factor(df$Type, levels <- c("Excite","Inh_wCR","Inhib"))

## Fig 4b ----
options(repr.plot.width=6, repr.plot.height=6)
g1 = ggplot(df, aes(x = Type, y = CorValue, fill = Type)) +
    geom_boxplot() + 
    scale_fill_manual(values = c("#ffcc5c","#65c3ba","#82d188")) + 
    theme_classic(base_size = 20) +
    theme(legend.position = "bottom")


g1

## Preprocessing for 4c-d -----

timepoints <- c("P14","P21","P28","P56")
regions <- c("Hippocampus","VisCortex")
combs <- matrix(c(timepoints[1:3],timepoints[2:4]),2,3,byrow = T)

cor_neurons_ci <- NULL
for (region in regions){
    for (c in 1:ncol(combs)){
        tps = combs[,c]
        for (ct in c("Excite","Inhib")){
            v1 <- type_neurons %>% as.data.frame() %>% select(contains(paste(tps[1],region,ct, sep = "_")))
            v2 <- type_neurons %>% as.data.frame() %>% select(contains(paste(tps[2],region,ct, sep = "_")))
            cor = cor.test(v1[,1],v2[,1])
            cor_neurons_ci[[paste(region,ct,c,sep = "_")]] <- c(unname(cor$estimate),
                                                                cor$conf.int,region,ct,paste(tps,collapse = "_"))
        }
    }
}


cor_neu_reg <- NULL
for (tp in timepoints){
    for (ct in c("Excite","Inhib")){
        v1 <- type_neurons %>% as.data.frame() %>% select(contains(paste(tp,"Hippocampus",ct, sep = "_")))
        v2 <- type_neurons %>% as.data.frame() %>% select(contains(paste(tp,"VisCortex",ct, sep = "_")))
        cor = cor.test(v1[,1],v2[,1])
        cor_neu_reg[[paste(tp,ct,sep = "_")]] <- c(unname(cor$estimate),
                                                   cor$conf.int,tp,ct)
    }
}


options(repr.plot.width=16, repr.plot.height=8)

cor_neurons_df <- as.data.frame(do.call('rbind',cor_neurons_ci)) %>% remove_rownames()
colnames(cor_neurons_df) <- c("Estimate","LCI","UCI","Region","Celltype","Transition")
cor_neurons_df <- cor_neurons_df %>% mutate_at(c("Estimate","LCI","UCI"),as.double)

## Fig 4c ----
g2 = ggplot(cor_neurons_df %>% filter(Region == "VisCortex"), 
       aes(x = Transition, y = Estimate, fill = Celltype)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, position = position_dodge(.9)) +
    theme_classic(base_size = 20) + 
    scale_fill_manual(values = c("#ffcc5c","#65c3ba")) + 
    facet_zoom(ylim = c(0.85, 0.95))


g2

cor_neu_reg_df <- as.data.frame(do.call('rbind',cor_neu_reg)) %>% remove_rownames()
colnames(cor_neu_reg_df) <- c("Estimate","LCI","UCI","Timepoint","Celltype")
cor_neu_reg_df <- cor_neu_reg_df %>% mutate_at(c("Estimate","LCI","UCI"),as.double)

## Fig 4d ----
g3 = ggplot(cor_neu_reg_df, 
       aes(x = Timepoint, y = Estimate, fill = Celltype)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1, position = position_dodge(.9)) +
    theme_classic(base_size = 20) + 
    scale_fill_manual(values = c("#ffcc5c","#65c3ba")) + 
    facet_zoom(ylim = c(0.85, 0.95))


g3

## sessionInfo()
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
# 
# other attached packages:
#   [1] RColorBrewer_1.1-3    circlize_0.4.15       ComplexHeatmap_2.13.1
# [4] patchwork_1.1.1       ggforce_0.3.3         ggplot2_3.3.5        
# [7] MetBrewer_0.2.0       viridis_0.6.2         viridisLite_0.4.0    
# [10] tibble_3.1.8          tidyr_1.2.0           dplyr_1.0.7          
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.8.3        png_0.1-7           assertthat_0.2.1   
# [4] digest_0.6.29       foreach_1.5.2       utf8_1.2.2         
# [7] IRdisplay_1.0       R6_2.5.1            repr_1.1.3         
# [10] stats4_4.0.1        evaluate_0.22       pillar_1.8.1       
# [13] GlobalOptions_0.1.2 rlang_1.0.2         uuid_0.1-4         
# [16] S4Vectors_0.28.1    GetoptLong_1.0.5    labeling_0.4.2     
# [19] polyclip_1.10-0     munsell_0.5.0       compiler_4.0.1     
# [22] pkgconfig_2.0.3     BiocGenerics_0.36.1 base64enc_0.1-3    
# [25] shape_1.4.6         htmltools_0.5.3     tidyselect_1.1.1   
# [28] gridExtra_2.3       IRanges_2.24.1      codetools_0.2-16   
# [31] matrixStats_0.62.0  fansi_0.5.0         crayon_1.5.1       
# [34] withr_2.5.0         MASS_7.3-51.6       jsonlite_1.8.7     
# [37] gtable_0.3.0        lifecycle_1.0.1     DBI_1.1.3          
# [40] magrittr_2.0.3      scales_1.2.1        cli_3.3.0          
# [43] farver_2.1.0        doParallel_1.0.17   ellipsis_0.3.2     
# [46] generics_0.1.2      vctrs_0.4.1         IRkernel_1.2       
# [49] rjson_0.2.21        Cairo_1.5-12.2      iterators_1.0.14   
# [52] tools_4.0.1         glue_1.6.2.9000     tweenr_1.0.2       
# [55] purrr_0.3.4         parallel_4.0.1      fastmap_1.1.0      
# [58] clue_0.3-61         colorspace_2.0-3    cluster_2.1.0      
# [61] pbdZMQ_0.3-4 
# 
# 
# 
