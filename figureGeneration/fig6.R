#!/bin/R
# Code for reproducing panels from figure 6
## Relating to the variability in exon inclusion
## levels over developmental time. Related
## panel for GO analysis in a separate file bc of
## sessioninfo considerations

## Setup ------
library(dplyr)
library(tidyr)
library(tibble)
library(superheat)
library(viridis)
library(MetBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

## Read in cell-type PSI values ----
type_alt <- read.table('../data/altExonUsage_devel_type.gz',sep = "\t", header = TRUE)


## Preprocessing ------
samples <- as.vector(sapply(c("Hippocampus","VisCortex"), 
       function(i) paste(c("P14","P21","P28","P56"),i,sep = "_")))

ctVar <- list()

for (sample in samples){
    ctVar[[sample]] <- type_alt %>% unite(GE,c("Exon","Gene"),sep = "::") %>%
    select(GE,c(contains(sample))) %>%
     filter(rowMeans(is.na(.[2:5])) <= 0.5) %>% 
     rowwise() %>%
      mutate(eVar = max(c_across(contains(sample)),na.rm = T) - 
                    min(c_across(contains(sample)),na.rm = T)) %>%
     select(GE,eVar)
}


meltedCTVar <- purrr::map_df(ctVar, ~as.data.frame(.x), .id="id")
head(meltedCTVar)
                            
meltedCTVar2 <- meltedCTVar %>% separate(id, into = c("Age","Region")) %>%
    mutate_at(c("Age"), as.factor)
                            
meltedCTVar3 <- meltedCTVar2 %>% group_by(Region,Age) %>% 
     mutate(m = median(eVar),d = abs(eVar - m)) %>% 
    select(-m)


## Find exons that change status. -------
## Mostly we identify the "constant exons"
changeStatusDF <- meltedCTVar2 %>% 
     pivot_wider(names_from = Age, values_from = eVar) %>% 
     drop_na() %>% mutate(s1 = P21-P14, s2 = P28-P21, s3 = P56-P28) %>% 
     mutate(s1s = case_when(abs(s1) < 0.1 ~ 0,s1 <= -0.1 ~ -1, s1 >= 0.1 ~ 1),
            s2s = case_when(abs(s2) < 0.1 ~ 0,s2 <= -0.1 ~ -1, s2 >= 0.1 ~ 1),
            s3s = case_when(abs(s3) < 0.1 ~ 0,s3 <= -0.1 ~ -1, s3 >= 0.1 ~ 1))


changeStatusDF2 <- changeStatusDF %>% select(Region,s1s,s2s,s3s) %>% 
    group_by(Region,s1s,s2s,s3s) %>% add_count() %>% distinct() %>% 
    group_by(Region) %>% mutate(x1 = 1, x2 = 2, x3 = 3, x4 = 4, 
                                y1 = x1, y2 = y1+s1s, y3 = y2+s2s, y4 = y3+s3s) %>%
    mutate(pchange = n*100/sum(n)) %>% 
    unite(ID,c(s1s,s2s,s3s),sep = ":") %>%
    unite(xy1,c('x1','y1'),sep ="_") %>% unite(xy2,c('x2','y2'),sep ="_") %>% 
    unite(xy3,c('x3','y3'),sep ="_") %>% unite(xy4,c('x4','y4'),sep ="_") 

tmpDF <- changeStatusDF2 %>% filter(n >= 100) %>% select(ID)

changeStatusDF3 <- changeStatusDF2 %>% filter(ID %in% tmpDF$ID) %>%
    pivot_longer(cols = c(xy1,xy2,xy3,xy4),values_to = "val") %>%
    separate(val, into = c("x","y")) %>% mutate_at(c('x','y'),as.double) %>%
    select(-name) %>% mutate(Age = case_when(x == 1 ~"P14", x == 2 ~ "P21",
                                      x == 3 ~"P28", x == 4 ~ "P56"))


ordering_changeStatusDF3 <- changeStatusDF3 %>% select(Region,ID, pchange) %>% 
    distinct() %>% group_by(Region) %>% arrange(-pchange,.by_group = T)

changeStatusDF3$ID <- factor(changeStatusDF3$ID, levels <- 
                             ordering_changeStatusDF3$ID[1:(nrow(ordering_changeStatusDF3)/2)])

options(repr.plot.width=10, repr.plot.height=5)
ggplot(changeStatusDF3, aes(x = x, y = y, color = Region, size = pchange/100)) + 
    geom_line() + facet_grid(Region~ID) +
    theme_bw(base_size = 15)

idMap <- changeStatusDF %>% unite(ID,c(s1s,s2s,s3s),sep = ":") %>%
    filter(ID %in% changeStatusDF3$ID) %>% select(GE,Region,ID)

constt_exons <- idMap %>% filter(ID == "0:0:0")

A = changeStatusDF %>% select(GE) %>% distinct() %>%
     mutate(category = case_when(GE %in% constt_exons$GE ~ "Constant", 
                                 TRUE ~ "Changing")) %>%
    group_by(category) %>% select(category) %>% add_count() %>% distinct()


A$category <- factor(A$category, levels = c("Constant","Changing"))

options(repr.plot.width=8, repr.plot.height=6)


## Obtain (and plot) constant versus changing exons -----
ggplot(A, aes(category,n)) +
    geom_bar(stat = "identity", width = 0.7) +
    theme_classic(base_size = 20)

## Row-normalize and cluster -----
set.seed(10)
options(repr.plot.width=4, repr.plot.height=12)
exonVar_mat <- changeStatusDF %>% filter(!(GE %in% constt_exons$GE)) %>%
    unite(GER,c(GE,Region),sep="%%") %>%
    select(GER,contains("P")) %>% column_to_rownames('GER')

exonVar_mat2 <- t(scale(t(exonVar_mat))) %>% as.data.frame() %>% drop_na() %>% as.matrix()

set.seed(1)

col_func <- colorRamp2(c(-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),viridis::viridis(9))

options(repr.plot.width=4, repr.plot.height=12)
H <- Heatmap(exonVar_mat2, show_row_names = FALSE, 
             cluster_columns = FALSE,
             km = 9, show_row_dend = FALSE,
             col = col_func, show_column_dend = FALSE,
            name = "eVar")

## Fig 6a ----
H <- draw(H)

## Save this cluster order for future analysis. This forms G1-G9 -----
ro <- row_order(H)
df <- NULL

df = lapply(1:length(ro), function(i) 
    data.frame(empID = i, 
               GER = rownames(exonVar_mat2)[as.vector(unname(unlist(ro[[i]])))]))
allDFs <- do.call('rbind',df)
              
exonVar_mat0 <- changeStatusDF %>% filter(GE %in% constt_exons$GE) %>%
    unite(GER,c(GE,Region),sep="%%") %>%
    select(GER,contains("P")) %>%
    mutate(empID = 0) %>% select(empID,GER)

allDFs <- rbind(allDFs,exonVar_mat0)

empiricalIDs <- allDFs %>% 
    separate(GER, into = c("GE","Region"),sep = "%%") %>% 
    group_by(empID,Region) %>% add_count() %>%
    ungroup() %>%
    mutate(pchange = n*100/sum(n))

eVar_wIDs <- inner_join(meltedCTVar2 %>% pivot_wider(names_from = Age, values_from = eVar),
                             empiricalIDs,by = c("Region","GE"))
eVar_wIDs$empID <- factor(eVar_wIDs$empID,levels = 0:9)


## Initialize preprocessing for figs 6b-c and associated supplemental figure ----
heatPlots <- NULL
heatPlots2 <- NULL
heatDFs <- NULL

for (region in c("Hippocampus","VisCortex")){
    heatPlots[[region]] <- NULL

    region_changeStatus <- empiricalIDs %>% filter(Region == region)
    oDF <- region_changeStatus %>% select(empID,pchange) %>% 
        distinct() %>% arrange(-pchange,.by_group = T)

    region_changeStatus$empID <- factor(region_changeStatus$empID, 
            levels = oDF$empID[1:length(unique(empiricalIDs$empID))])

N = nrow(eVar_wIDs %>% filter(Region == region))
    
    for (ix in unique(eVar_wIDs$empID)){
        h <- eVar_wIDs %>% 
            filter(empID == ix & Region == region) %>%
            column_to_rownames('GE') %>%
            select(contains("P",ignore.case = F))

        hp <- superheat(h,
            col.dendrogram = FALSE,
            pretty.order.rows = TRUE,
            #n.clusters.rows = 2,
            heat.pal.values = seq(0,1,by = 0.2),
            print.plot = F)

        ## heatmap of exon variability
        ho <- h[hp$order.rows,hp$order.cols]
        hom <- reshape2::melt(ho %>% rownames_to_column('GE'),"GE")
        hom$GE <- factor(hom$GE, levels <- rownames(ho))
        heatDFs[[region]][[ix]] <- hom

        heatPlots[[region]][[ix]] <- ggplot(hom,aes(x = variable,y = GE, fill = value)) +
            geom_tile() +
            #scale_fill_viridis(limits = c(0, 1)) +
            scale_fill_distiller(palette = "RdYlGn",limits = c(0,1)) +
            theme_classic(base_size = 16) + 
            theme(axis.text.y = element_blank(),axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                  legend.position = "none")
        
        ## lineplots2
        hom3 <- hom %>% group_by(variable) %>% mutate(mVal = mean(value)) ## Originally had hom2
        
        if(ix != 0){
            linePlots2[[region]][[ix]] <- ggplot(hom3,aes(x = variable,y = value,group = GE)) +
                geom_point(col = "lightgrey",alpha = 0.8) +
                geom_line(col = "lightgrey",alpha = 0.5) +
                scale_fill_viridis() +
                theme_classic(base_size = 16) + 
                ylim(0,1) +
                ggtitle(label = paste(nrow(ho2),"exons")) + 
                theme(axis.text.y = element_blank(),axis.title.x = element_blank(),
                axis.title.y = element_blank()) +
                geom_line(aes(x = variable,y = mVal),col = "forestgreen",alpha = 0.6,size = 10*nrow(h)/N) +
                geom_point(aes(x = variable,y = mVal),col = "forestgreen",size = 10*nrow(h)/N)
        } else {
            linePlots2[[region]][[ix]] <- ggplot(hom3,aes(x = variable,y = value,group = GE)) +
                geom_point(col = "lightgrey",alpha = 0.8) +
                geom_line(col = "lightgrey",alpha = 0.5) +
                scale_fill_viridis() +
                theme_classic(base_size = 16) + 
                ylim(0,1) +
                ggtitle(label = paste(nrow(ho2),"exons")) + 
                theme(axis.text.y = element_blank(),axis.title.x = element_blank(),
                axis.title.y = element_blank()) +
                geom_line(aes(x = variable,y = mVal),col = "forestgreen",alpha = 0.6,size = 2) +
                geom_point(aes(x = variable,y = mVal),col = "forestgreen",size = 2)
        } 
    }
}

## Get plots for Hippocampus ------
empID_CH <- factor(eVar_wIDs$empID,levels <- c(0:9))

options(repr.plot.width=16, repr.plot.height=3)
patchwork::wrap_plots(linePlots2$Hippocampus[levels(empID_CH)],nrow = 1)

options(repr.plot.width=16, repr.plot.height=4)
patchwork::wrap_plots(heatPlots$Hippocampus[levels(empID_CH)],nrow = 1)


## Get plots for Visual cortex ------
empID_CH <- factor(c(0:9),levels <- c(0:9))

options(repr.plot.width=16, repr.plot.height=3)
patchwork::wrap_plots(linePlots2$VisCortex[levels(empID_CH)],nrow = 1)

options(repr.plot.width=16, repr.plot.height=4)
patchwork::wrap_plots(heatPlots$VisCortex[levels(empID_CH)],nrow = 1)


#sessionInfo()
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
#   [1] RColorBrewer_1.1-3    circlize_0.4.15       ComplexHeatmap_2.13.1
# [4] ggplot2_3.3.5         MetBrewer_0.2.0       viridis_0.6.2        
# [7] viridisLite_0.4.0     superheat_1.0.0       tibble_3.1.8         
# [10] tidyr_1.2.0           dplyr_1.0.7          
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.8.3        png_0.1-7           assertthat_0.2.1   
# [4] digest_0.6.29       foreach_1.5.2       utf8_1.2.2         
# [7] IRdisplay_1.0       plyr_1.8.6          R6_2.5.1           
# [10] repr_1.1.3          stats4_4.0.1        evaluate_0.22      
# [13] pillar_1.8.1        GlobalOptions_0.1.2 rlang_1.0.2        
# [16] uuid_0.1-4          S4Vectors_0.28.1    GetoptLong_1.0.5   
# [19] labeling_0.4.2      stringr_1.4.0       munsell_0.5.0      
# [22] compiler_4.0.1      pkgconfig_2.0.3     BiocGenerics_0.36.1
# [25] base64enc_0.1-3     shape_1.4.6         htmltools_0.5.3    
# [28] tidyselect_1.1.1    gridExtra_2.3       IRanges_2.24.1     
# [31] codetools_0.2-16    matrixStats_0.62.0  fansi_0.5.0        
# [34] crayon_1.5.1        withr_2.5.0         jsonlite_1.8.7     
# [37] gtable_0.3.0        lifecycle_1.0.1     DBI_1.1.3          
# [40] magrittr_2.0.3      scales_1.2.1        cli_3.3.0          
# [43] stringi_1.7.8       reshape2_1.4.4      farver_2.1.0       
# [46] doParallel_1.0.17   ellipsis_0.3.2      generics_0.1.2     
# [49] vctrs_0.4.1         IRkernel_1.2        rjson_0.2.21       
# [52] iterators_1.0.14    tools_4.0.1         Cairo_1.5-12.2     
# [55] glue_1.6.2.9000     purrr_0.3.4         parallel_4.0.1     
# [58] fastmap_1.1.0       clue_0.3-61         colorspace_2.0-3   
# [61] cluster_2.1.0       pbdZMQ_0.3-4        patchwork_1.1.1    
# 
