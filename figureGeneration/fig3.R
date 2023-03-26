library(ComplexHeatmap)
library(circlize)
library(simplifyEnrichment)
library(dplyr)
library(tidyr)
library(tibble)
library(fastcluster)
library(MetBrewer)
library(clusterProfiler)
library("org.Mm.eg.db")
library(biomaRt)
library(igraph)
library(ggplot2)
library("stringr")


## hVEX preprocessing -------

fullMatrix <- data.table::fread('../data/HVEx_pairwise_matrix',header = T,sep = "\t")

subs <- abs(fullMatrix[,-c(1:2)])
subs_m <- as.matrix(subs)

DR <- dist(subs_m)
DR <- as.data.frame(as.matrix(DR))

DR2 <- DR %>% mutate(numNA = ncol(DR) - rowMeans(is.na(.))*ncol(DR))
tooManyNAs_R <- unname(which(DR2$numNA <= ncol(DR)-200 ))

subs_m <- as.matrix(subs)
subs_m <- subs_m[-c(tooManyNAs_R),]

## Plot hVEX --------
col_func <- colorRamp2(seq(1,0,length.out = 7),
                       c("#EA7580","#F4999E","#F7BE9F","#89C1A5","#14A7B3","#0A7AAF","#172869"))

regions <- c("VisCortex","Hippocampus","Striatum","Thalamus","Cerebellum")
reg_cols <- c("#ddbbc6","#fdc681","#9c9648","#7e9a66","#50939a")

cts <- c("Astro","Oligo","ExciteNeuron","InhibNeuron")
ct_cols <- c("#d86472","#841f1f","#939ed0","#5782c2")

ages <- c("P14","P21","P28","P56")
age_cols <- c("#ae3b24","#c8812a","#56632c","#f68d62")

detectChange <- function(cName){
  components <- unlist(strsplit(cName,"-|_"))
  if(components[1] == components[4]){
    Age = components[1] }
  else{Age = "Different"}
  if(components[2] == components[5]){
    Region = components[2] }
  else{Region = "Different"}
  if(components[3] == components[6]){
    Celltype = components[3] }
  else{Celltype = "Different"}
  return(c(Age,Region,Celltype))
}

consensusDF <- as.data.frame(t(unname(sapply(colnames(subs_m), 
                                             function(c) detectChange(c)))))

ha0 = HeatmapAnnotation(
  Celltype = consensusDF$V3,
  Region = consensusDF$V2,
  Age = consensusDF$V1,

  col = list(Age = c(structure(age_cols, names = ages),"Different" = "#A7A7A7"),
             Region = c(structure(reg_cols,names = regions),"Different" = "#E6E6E6"),
             Celltype = c(structure(ct_cols,names = cts),"Different" = "#A7A7A7")),
  
  annotation_legend_param = list(
    Age = list(
      title = "Age",
      at = c(ages,"Different"),
      labels = c(ages,"Different")
    ),
    Region = list(
      title = "Region",
      at = c(regions,"Different"),
      labels = c(regions,"Different")
    ),
    Celltype = list(
      title = "cellType",
      at = c(cts,"Different"),
      labels = c(cts,"Different")
    )
  )
)


fh = function(x) fastcluster::hclust(dist(x),method = "ward.D")

pa = cluster::pam(subs_m, k = 4)
pc = cluster::pam(t(subs_m), k = 2)

ht = Heatmap(subs_m,name = "dPSI", na_col = "#f4f4f5",
        cluster_rows = fh, cluster_columns = fh,
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE,show_column_dend = FALSE,
        col = col_func, top_annotation = ha0, 
        row_split = paste0("A", pa$clustering),
        column_split = pc$clustering,
        use_raster = TRUE)

ht <- draw(ht)

pdf('../Plots/hVEX_CompHeatamp.pdf',width = 14,height = 10)
ht
dev.off()

### Graphical representation on top -----

col_clust <- as.data.frame(pc$clustering) %>% 
  tibble::rownames_to_column("Comp")
colnames(col_clust)[2] <- "Cluster"
col_clust_ct <- col_clust %>% separate(Comp, into = c("from","to"), sep = "-") %>% 
  separate(from, into = c("age1","reg1","ct1"),sep = "_") %>%
  separate(to, into = c("age2","reg2","ct2"),sep = "_") %>% 
  dplyr::select(ct1,ct2,Cluster) %>% 
  rowwise() %>%
  mutate(s = paste(sort(c(ct1,ct2)), collapse = "_")) %>% 
  dplyr::select(-c(ct1,ct2)) %>% 
  group_by_all() %>%
  summarize(importance = n()) %>% separate(s, into = c("source","target")) %>%
  ungroup() %>%
  group_by(source) %>% mutate(imp = importance / sum(importance)) %>% ungroup()

links_1 <- col_clust_ct %>% filter(Cluster == 1) %>% ungroup() %>% 
  dplyr::select(-c(Cluster,importance)) %>% as.data.frame()

links_2 <- col_clust_ct %>% filter(Cluster == 2) %>% ungroup() %>% 
  dplyr::select(-c(Cluster,importance)) %>% as.data.frame()

nodes <- data.frame(
  name=c("Astro","Oligo","ExciteNeuron","InhibNeuron"))

network_1 <- graph_from_data_frame(d=links_1, vertices=nodes, directed=F)
network_2 <- graph_from_data_frame(d=links_2, vertices=nodes, directed=F)

coul <- c("#d86472","#841f1f","#939ed0","#5782c2")
my_color <- coul[as.numeric(factor(nodes$name, levels = c("Astro","Oligo","ExciteNeuron","InhibNeuron")))]

A = plot(network_1, edge.width=E(network_1)$imp*10, vertex.color = my_color)

B = plot(network_2, edge.width=E(network_2)$imp*10, vertex.color = my_color)

legend("bottomleft", legend=levels(as.factor(factor(nodes$name, 
                                                    levels = c("Astro","Oligo","ExciteNeuron","InhibNeuron")))) , 
       col = coul , bty = "n", 
       pch=20 , pt.cex = 2, cex = 1.5, text.col=coul , horiz = FALSE, inset = c(0, 0))


## Correlation of PSI : supp figure ------
fullDF <- full_join(develDF,adultDF) %>% unite(GE,c("Exon","Gene"),sep = "::") %>% 
    column_to_rownames("GE")
C = cor(fullDF,use = "complete",method = "spearman")
Cnames <- data.frame(CT = colnames(C)) %>% separate(CT, into = c("Age","Region","Celltype"))

regions <- c("VisCortex","Hippocampus","Striatum","Thalamus","Cerebellum")
reg_cols <- c("#ddbbc6","#fdc681","#9c9648","#7e9a66","#50939a")

cts <- c("Astro","Oligo","ExciteNeuron","InhibNeuron")
ct_cols <- c("#d86472","#841f1f","#939ed0","#5782c2")

ages <- c("P14","P21","P28","P56")
age_cols <- c("#ae3b24","#c8812a","#56632c","#f68d62")

ha0 = HeatmapAnnotation(
  Celltype = Cnames$Celltype,
  Region = Cnames$Region,
  Age = Cnames$Age,
  col = list(Age = c(structure(age_cols, names = ages)),
             Region = c(structure(reg_cols,names = regions)),
             Celltype = c(structure(ct_cols,names = cts))),
  
  annotation_legend_param = list(
    Age = list(
      title = "Age",
      at = ages,
      labels = ages
    ),
    Region = list(
      title = "Region",
      at = regions,
      labels = regions
    ),
    Celltype = list(
      title = "cellType",
      at = cts,
      labels = cts
    )
  ),
  height = unit(8, "mm")
)

col_func <- colorRamp2(seq(0,1,length.out = 7),c("#EA7580","#F4999E","#F7BE9F","#89C1A5","#14A7B3","#0A7AAF","#172869"))

options(repr.plot.width=12, repr.plot.height=10)
ht = Heatmap(C,name = "Spearman",
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = TRUE,show_column_dend = FALSE,
        col = col_func, 
        top_annotation = ha0,row_km = 2,column_km = 2,
        use_raster = FALSE)
ht <- draw(ht)

pdf('../Plots/adultAndDevel_spearmanCorr.pdf',width = 12,height = 10,useDingbats = F)
ht
dev.off()

## Circos plot -------
### Adult only - visualize CELL TYPE variability ------

all_eVex <- read.table('../data/threeAxes_eVex_df_forCircos_adult',sep = "\t", header = T)

mat_adult <- all_eVex %>% filter(Axis %in% c('Adult_BRspec','Adult_CTspec')) %>% 
    select(-Axis) %>%
    pivot_wider(names_from = To, values_from = n) %>%
    column_to_rownames('From') %>% as.matrix()

ct_list_adult <- unique(unlist(dimnames(mat_adult)))

ct_split_adult <- matrix(unlist(strsplit(ct_list_adult,split = ":")),
                         nrow = length(ct_list_adult),byrow = T) %>% as.data.frame()
ct_split_adult <- ct_split_adult %>% group_by(V1 = factor(V1,levels = c("VIS","HIPP","STRI","THAL","CEREB"))) %>% 
        mutate(reg_id = cur_group_id()) %>% 
    ungroup() %>% group_by(V3 = factor(V3,levels = c("Astro","Oligo","ExciteNeuron","InhibNeuron"))) %>% 
        mutate(ct_id = cur_group_id()) %>%
    rename(Region = V1, Age = V2, CT = V3) %>% select(-Age)

regions <- ct_split_adult$Region
cts <- ct_split_adult$CT

group_ct_adult = structure(cts, names = ct_list_adult)
group_reg_adult = structure(regions, names = ct_list_adult)

ct_split_adult$ct_id <- plyr::mapvalues(ct_split_adult$ct_id, 
                                        1:4, met.brewer("Egypt", 4)[1:4])
col_ct_adult = structure(ct_split_adult$ct_id,names = ct_list_adult)

ct_split_adult$reg_id <- plyr::mapvalues(ct_split_adult$reg_id, 
                                         1:5, c("#debbc7","#fdc681","#9c9648","#7e9a66","#50939a"))
col_reg_adult = structure(ct_split_adult$reg_id,names = ct_list_adult)

options(repr.plot.width=12, repr.plot.height=12)

rn_list_adult <- rownames(mat_adult)
rn_split_adult <- matrix(unlist(strsplit(rn_list_adult,split = ":")),
                         nrow = length(rn_list_adult),byrow = T) %>% as.data.frame()
rn_list_adult <- rn_split_adult %>% arrange(V3,V1,V2) %>% unite(rct,c(V1,V2,V3),sep = ":")

cn_list_adult <- colnames(mat_adult)
cn_split_adult <- matrix(unlist(strsplit(cn_list_adult,split = ":")),
                         nrow = length(cn_list_adult),byrow = T) %>% as.data.frame()
cn_list_adult <- cn_split_adult %>% arrange(V3,V1,V2) %>% unite(cct,c(V1,V2,V3),sep = ":")

### Plot circos ------
circos.clear()
pdf(file = "../Plots/adultOnly_CTvar.pdf", width = 12, height = 12, useDingbats = F)
chordDiagram(mat_adult,directional = 1,diffHeight = mm_h(5),
             group = group_reg_adult, 
             grid.col = col_reg_adult,
             order = unique(c(rn_list_adult$rct,cn_list_adult$cct)), big.gap = 15,
             annotationTrack = NULL,
            preAllocateTracks = list(list(track.height = mm_h(5)),
                             list(track.height = mm_h(5))))

highlight.sector(ct_list_adult[which(ct_split_adult$CT == "Astro")], track.index = 2, col = "#d86472", 
    text = "Astro", cex = 2, text.col = "white", niceFacing = TRUE)
highlight.sector(ct_list_adult[which(ct_split_adult$CT == "Oligo")], track.index = 2, col = "#841f1f", 
    text = "Oligo", cex = 2, text.col = "white", niceFacing = TRUE)
highlight.sector(ct_list_adult[which(ct_split_adult$CT == "ExciteNeuron")], track.index = 2, col = "#939ed0", 
    text = "Excite", cex = 2, text.col = "white", niceFacing = TRUE)
highlight.sector(ct_list_adult[which(ct_split_adult$CT == "InhibNeuron")], track.index = 2, col = "#5782c2", 
    text = "Inhib", cex = 2, text.col = "black", niceFacing = TRUE)

highlight.sector(ct_list_adult[which(ct_split_adult$Region == "HIPP")], 
                 track.index = 1, col = "#fdc681", 
    text = "HIPP", cex = 2, text.col = "black", niceFacing = TRUE)
highlight.sector(ct_list_adult[which(ct_split_adult$Region == "VIS")], 
                 track.index = 1, col = "#debbc7", 
    text = "VIS", cex = 2, text.col = "white", niceFacing = TRUE)
highlight.sector(ct_list_adult[which(ct_split_adult$Region == "CEREB")], 
                 track.index = 1, col = "#50939a", 
    text = "CEREB", cex = 2, text.col = "white", niceFacing = TRUE)
highlight.sector(ct_list_adult[which(ct_split_adult$Region == "THAL")], 
                 track.index = 1, col = "#7e9a66", 
    text = "THAL", cex = 2, text.col = "white", niceFacing = TRUE)
highlight.sector(ct_list_adult[which(ct_split_adult$Region == "STRI")], 
                 track.index = 1, col = "#9c9648", 
    text = "STRI", cex = 2, text.col = "black", niceFacing = TRUE)
dev.off()

## Quantify within vs. across brain region changes ------
within_vs_acrossBR <- hv_df_forCircos_adult %>% separate(From,into = c("FR","FAge","FCT"),sep = ":") %>% 
     separate(To,into = c("TR","TAge","TCT"),sep = ":") %>% ungroup() %>%
     select(FR,TR,n) %>% 
     mutate(Status = case_when(FR == TR ~ "Within", TRUE ~ "Across")) %>%
    group_by(Status) %>%
    summarize(across(n, sum),.groups="keep") %>% ungroup() %>%
    mutate(Perc = n*100/sum(n))

options(repr.plot.width=8, repr.plot.height=6)
pdf('../Plots/within_vs_acrossBR.pdf',8,6, useDingbats = FALSE)
ggplot(within_vs_acrossBR, aes(x = Status, y = Perc, fill = Status)) +
    geom_bar(stat = "identity", width = 0.6) +
    theme_classic(base_size = 20) + 
    scale_fill_brewer(palette = "Greys")
dev.off()

## Compare developmental age vs. brain region specific changes for the same cell type ------
### Developmental -----
cts <- c("Astro","Oligo","Excite","Inhib")

exonVar <- list()

for (sample in cts){
    exonVar[[sample]] <- develDF %>% #select(-contains("VisCortex")) %>% ## keeping brain region constant
    unite(GE,c("Exon","Gene"),sep = "::") %>%
    select(GE,c(contains(sample))) %>%
     filter(rowMeans(is.na(.[2:5])) <= 0.5) %>% 
     rowwise() %>%
      mutate(eVar = max(c_across(contains(sample)),na.rm = T) - 
                    min(c_across(contains(sample)),na.rm = T)) %>%
     select(GE,eVar)
}


meltedExonVar_devel <- purrr::map_df(exonVar, ~as.data.frame(.x), .id="id")
meltedExonVar_devel$Type <- "Devel"

### Adult -----
exonVar <- list()

for (sample in cts){
    exonVar[[sample]] <- adultDF %>% unite(GE,c("Exon","Gene"),sep = "::") %>%
    select(GE,contains(sample)) %>%
     filter(rowMeans(is.na(.[2:6])) <= 0.6) %>% 
     rowwise() %>%
      mutate(eVar = max(c_across(contains(sample)),na.rm = T) - 
                    min(c_across(contains(sample)),na.rm = T)) %>%
     select(GE,eVar)
}


meltedExonVar_adult <- purrr::map_df(exonVar, ~as.data.frame(.x), .id="id")
meltedExonVar_adult$Type <- "Adult"

combDF_5BR_4TP <- rbind(meltedExonVar_adult,meltedExonVar_devel)



### Plot comparisons -----
options(repr.plot.width=10, repr.plot.height=8)
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
}

a <- ggplot(combDF_5BR_4TP, aes(x = id, y = eVar, fill = Type)) +
    geom_boxplot(outlier.shape = NA) + theme_classic(base_size = 15) +
     stat_summary(fun = median, fun.data = give.n, geom = "text",
                position = position_dodge(width = 0.75)) +
    coord_cartesian(ylim = c(0,0.6)) +
    scale_fill_brewer(palette = "Greys") +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")


pdf('../Plots/exonVariability_devel_vs_adult.pdf',10,8,useDingbats = F)
a 
dev.off()

set.seed(1)

## Plot EVEx categories -------
A <- read.table('../data/EVEx_4axes', sep = "\t", header = TRUE)
options(repr.plot.width=5, repr.plot.height=9)

col_func <- colorRamp2(seq(1,0,length.out = 7),
                       c("#EA7580","#F4999E","#F7BE9F","#89C1A5","#14A7B3","#0A7AAF","#172869"))

H <- Heatmap(as.matrix(trueHigh), show_row_names = FALSE, km = 5, show_row_dend = FALSE,
             col = col_func, show_column_dend = FALSE,
            name = "dPSI")
H <- draw(H)

pdf('../Plots/EVEx_5categories_complexHeatmap.pdf', width= 5, height = 9, useDingbats = F)
H
dev.off()

ro <- row_order(H)       


## Get length distributions of EVEx categories ------
A2 <- A %>% separate("Exon", into = c("chr","start","end","strand"),sep = "_") %>%
    mutate_at(c("start","end"), as.integer) %>% mutate(length = end - start) %>% 
    select(Gene,Cluster,length)

A2$Cluster <- factor(A2$Cluster,levels <- c("develAll","devTime","allCT","adultAll","adultCT"))

options(repr.plot.width=8, repr.plot.height=8)

pdf('../Plots/EVEx_categoryLengths_boxplot.pdf',8,7,useDingbats = FALSE)
ggplot(A2, aes(x = Cluster, y = length, color = Cluster)) +
    geom_boxplot(alpha = 0.2, size=2) + theme_classic(base_size = 20) +
    scale_color_manual(values = met.brewer("Renoir", 5)[c(4,5,3,1,2)]) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) + 
    geom_signif(comparisons = list(c("adultCT","adultAll"),
                                  c("devTime","allCT")),
               y_position = c(200,210,190),tip_length = 0.01) +
    ylim(0,300)
dev.off()

## Get protein-coding capacity of EVEx categories ------
protCod <- as.data.frame(table(clusterInfo_PC %>% select(Cluster,type))) %>% 
        group_by(Cluster) %>% mutate(Prop = Freq *100/ sum(Freq))

protCod$Cluster <- factor(protCod$Cluster,levels <- c("develAll","devTime","allCT","adultAll","adultCT"))

pdf('../Plots/EVEx_categoryNonProtCoding_barplot.pdf',8,6,useDingbats = FALSE)
ggplot(protCod %>% filter(type == "nonPC"), aes(x = type, y = Prop, fill = Cluster)) +
    geom_bar(stat = "identity",position = "dodge") +
    theme_classic(base_size = 20) +
    scale_fill_manual(values = met.brewer("Renoir", 5)[c(4,5,3,1,2)]) +
    theme(legend.position = "bottom")
dev.off()

## Superfamily annotations ------
supFamDF_wCats_filtered <- read.table('../data/supFamDF_wEVExCats_filtered',
                       sep = "\t", header = T)

supFamDF_wCats_filtered_mat <- supFamDF_wCats_filtered %>% select(-n) %>%
  pivot_wider(names_from = Cluster, values_from = perc) %>% 
  ungroup() %>% filter(rowSums(is.na(.[2:7])) <5 & !is.na(lowVar_NULL) | 
                         is.na(lowVar_NULL)) %>%
  replace(is.na(.),0) %>%
  column_to_rownames('exon_associated_domains_superfamilies_meaning') %>%
  as.matrix()


### Take only interesting superfamilies ---------
supFamDF_wCats_filtered_mat_reduced <- supFamDF_wCats_filtered_mat %>% 
  as.data.frame() %>% 
  filter_at(vars(-lowVar_NULL),any_vars(. > 4)) %>% 
  as.matrix()

ord <- hclust( dist(supFamDF_wCats_filtered_mat_reduced, method = "euclidean"), method = "ward.D" )$order
ord2 <- hclust( dist(t(supFamDF_wCats_filtered_mat_reduced), method = "euclidean"), method = "ward.D" )$order

supFamDF_wCats_filtered_reduced <- supFamDF_wCats_filtered %>% 
  filter(exon_associated_domains_superfamilies_meaning %in% rownames(supFamDF_wCats_filtered_mat_reduced))
supFamDF_wCats_filtered_reduced$Cluster <- factor(supFamDF_wCats_filtered_reduced$Cluster,
                                                  levels = c("lowVar_NULL","develAll","devTime",
                                                             "allCT","adultAll","adultCT"))

colnames(supFamDF_wCats_filtered_reduced)[2] <- "Superfamily"
supFamDF_wCats_filtered_reduced$Superfamily <- factor(supFamDF_wCats_filtered_reduced$Superfamily,
                                                      levels = rownames(supFamDF_wCats_filtered_mat_reduced)[ord])

g <- ggplot(supFamDF_wCats_filtered_reduced, 
       aes(x = Cluster, y = Superfamily, fill = perc)) +
  geom_tile() +
  geom_text(aes(label= round(perc,2))) +
  scale_fill_distiller(palette = 3,limits = c(0,10),oob = squish) + theme_classic() +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 12))


pdf('../Plots/EVex_superFamilies_percent.pdf',8,8)
g
dev.off()


## GO analysis adult brain-region specificity ------
evexGenes_uniq <- evex %>%
  separate(Gene, into = c("Stable","version"), sep = "[.]") %>%
  dplyr::select(-Exon) %>% distinct() %>%
  group_by(Stable,Cluster) %>% add_count() %>% 
  filter(n <= 1)

egoList <- lapply(unique(evexGenes_uniq$Cluster),
                  function(i) enrichGO(gene = as.vector(evexGenes_uniq %>% 
                                                          filter(Cluster == i))$Stable, 
                                       universe = genesMouse$Gene.stable.ID,
                                       OrgDb = org.Mm.eg.db, keyType = 'ENSEMBL',
                                       ont = "BP",readable = TRUE))

names(egoList) <- unique(evexGenes_uniq$Cluster)

e = enrichplot::dotplot(egoList$adultAll) + scale_color_viridis_c(option = "plasma")
pdf('../Plots/EVEx_adultAll_dotPlot_GO_BP.pdf',
    width = 8, height = 8,useDingbats = F)
print(e)
dev.off()

