#! /bin/R 
# Ternary / simplex plots for representing isoform variability

# Set up ----------

library(dplyr)
library(tidyr)
library(data.table)
library(cowplot)
library(ggplot2)
library(ggtern)
library(MetBrewer)
library(ggrastr)
library(ggsignif)
library(igraph)
library(ComplexHeatmap)
library(circlize)

# Cell type -------
## Read in cell type data --------
gW <- fread('../data/genomeWide_withContribs.gz',
                        sep = "\t",header = T) 

gW <- gW %>% separate("isoID",c("Gene","ID"),sep = "-",extra = "merge",remove = FALSE) 
gW <- gW %>% mutate(explainability = case_when(abs(missingVar) <= 0.1 ~ "Explainable",
                                                             abs(missingVar) > 0.1 ~ "Missing"))
head(gW)
plotDir <- '../Plots/'
if(!dir.exists){dir.create(plotDir)}

## Cell-type specific ternary plot --------

gW$Celltype <- factor(gW$Celltype, 
                             levels = c("Progenitor","InhibNeuron","ExciteNeuron",
                                        "Astro","Oligo","Immune"))


pdf(file.path(plotDir,"ternary_coloredByCT.pdf"),15,15,useDingbats = FALSE)
ggtern(data=gW %>% filter(Age >= 0.1 | Region >= 0.1 | Subtype >= 0.1),
       aes(x=Age, y=Region, z=Subtype, col = Celltype),
       aes(x,y,z)) + 
  facet_wrap(~Celltype,nrow = 3, ncol = 3) +
  geom_point_rast(size = 0.01, alpha = 0.3) +
  theme_arrowdefault() +
  theme(legend.position = "none") + 
  labs( title= "Isoform variability contribution") +
  geom_Tline(Tintercept=.5,colour='black') + 
  geom_Lline(Lintercept=.5, colour='black') +
  geom_Rline(Rintercept=.5, colour='black') +
  scale_color_met_d("Renoir")
dev.off()

### Excitatory neuron sub plots ----- 

ct <- "ExciteNeuron"

densInset = ggplot(getMean %>% filter(Celltype == ct), 
                 aes(x = m, color = Celltype)) + 
    geom_density(size = 2) +
    theme_classic() +
    scale_color_manual(values = colPalette) +
    theme(legend.position = "none") +
    xlab("Mean variability")

boxInset <- ggplot(getMeanMelted %>% filter(Celltype == ct), 
                 aes(x = Axis, y = Variability, fill = Axis)) + 
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    scale_fill_manual(values = met.brewer("Peru1",4)[-1]) +
    theme(legend.position = "none")


gEN = ggtern(data=normalizedData %>% filter(Celltype == ct),
           aes(x=Age, y=Region, z=Subtype, col = Status, label = Gene),
           aes(x,y,z)) + 
    geom_point_rast(size = 0.01, alpha = 0.3) +
    theme_arrowdefault() +
    theme(legend.position = "none") + 
    ggtitle(ct) +
    geom_Tline(Tintercept=.5,colour='black') + 
    geom_Lline(Lintercept=.5, colour='black') +
    geom_Rline(Rintercept=.5, colour='black') +
    scale_color_met_d("Peru1")

pdf(file.path(plotDir,paste0(ct,"_tern.pdf")),12,8,useDingbats = FALSE)
print(gEN)
dev.off()

pdf(file.path(plotDir,paste0(ct,"_dens.pdf")),12,8,useDingbats = FALSE)
print(densInset)
dev.off()

pdf(file.path(plotDir,paste0(ct,"_box.pdf")),12,8,useDingbats = FALSE)
print(boxInset)
dev.off()


## Normalized data for classification --------

normalizedData <- gW %>% filter(Age >= 0.1 | Region >= 0.1 | Subtype >= 0.1) %>% 
  mutate(sum = (Age+Region+Subtype)) %>% 
  mutate_at(c("Age","Region","Subtype"),~ ./sum) %>% select(-sum) %>%
  mutate(Status = case_when(Age > 0.5 & Region < 0.5 & Subtype < 0.5 ~ "AgeSpec",
                            Age < 0.5 & Region > 0.5 & Subtype < 0.5 ~ "RegSpec",
                            Age < 0.5 & Region < 0.5 & Subtype > 0.5 ~ "SubtypeSpec",
                            TRUE ~ "Complex"))

normalizedData$Celltype <- factor(normalizedData$Celltype, 
                                  levels = c("Progenitor","InhibNeuron","ExciteNeuron",
                                             "Astro","Oligo","Immune"))

normalizedData$Status <- factor(normalizedData$Status, 
                                levels = c("Complex","AgeSpec","RegSpec","SubtypeSpec"))

### Ternary all with classification --------

pdf(file.path(plotDir,"AllCT_4Quadrants.pdf"),15,15,useDingbats = FALSE)
ggtern(data=normalizedData ,
       aes(x=Age, y=Region, z=Subtype, col = Status, label = Gene),
       aes(x,y,z)) + 
  facet_wrap(~Celltype,nrow = 3,ncol = 3) +
  geom_point_rast(size = 0.01, alpha = 0.3) +
  theme_arrowdefault() +
  theme(legend.position = "none") + 
  labs( title= "Ternary Plot and Filled Contour") +
  geom_Tline(Tintercept=.5,colour='black') + 
  geom_Lline(Lintercept=.5, colour='black') +
  geom_Rline(Rintercept=.5, colour='black') +
  scale_color_met_d("Peru1")
dev.off()

## Cross quadrants preprocessing --------

countQuadrants <- normalizedData %>% filter(Celltype != "Vasc") %>%
  group_by(Celltype, Status, Gene) %>% add_count(name = "multTranPerQuadrant") %>% 
  select(Celltype, Status,Gene,multTranPerQuadrant) %>% ungroup() %>% distinct() %>% 
  group_by(Celltype,Gene) %>%
  add_count(name = "AcrossQuadrants") %>% distinct()

summarizeQuadrantInfo <- countQuadrants %>% ungroup() %>%
  select(Celltype,Gene,AcrossQuadrants) %>% distinct()

df <- as.data.frame(table(summarizeQuadrantInfo[,c(1,3)])) %>%
  group_by(Celltype) %>%
  mutate(Perc = Freq*100/sum(Freq)) %>% ungroup()

### Cross quadrants plot  --------

pdf(file.path(plotDir,"percGenesAcrossQuadrants.pdf"),12,8,useDingbats = FALSE)
ggplot(df, aes(x = AcrossQuadrants, y = Perc, fill = Celltype)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_met_d("Renoir") +
  theme_classic() +
  xlab("Number of quadrants") +
  ylab("Percentage of genes")
dev.off()

## Cumulative plot -----

getMean = gW %>% filter(isoID %in% (normalizedData %>% filter(Status == "Complex") %>% .$isoID)) %>%
  rowwise() %>% mutate(m = mean(c(Age,Region,Subtype)))

getMeanMelted <- getMean %>% select(isoID,Age,Region,Subtype,Celltype) %>% 
  pivot_longer(c(Age,Region,Subtype),names_to = "Axis",values_to = "Variability")

getMean <- getMean %>% mutate(Broad = case_when(Celltype %in% c("ExciteNeuron","InhibNeuron") ~ "Neuron",
                                                Celltype == "Progenitor" ~ "Progenitor",
                                                Celltype %in% c("Astro","Oligo","Immune") ~ "Glia",
                                                Celltype == "Vasc" ~ "Vasc"))

getMeanCumu <- NULL
for(ct in unique(getMean$Broad)){
    ctDF <- getMean %>% filter(Broad == ct)
    for (i in seq(0,1,length.out = 51)){
    perc <- nrow(ctDF %>% filter(m <= i))/nrow(ctDF)
    getMeanCumu[[ct]] <- c(getMeanCumu[[ct]],perc)
  }
}

cumuDF <- as.data.frame(do.call('cbind',getMeanCumu))
cumuDF$meanVar <- seq(0,1,length.out = 51)
cumuDF2 <- reshape2::melt(cumuDF,id.vars = "meanVar",variable.name = "Broad",
                          value.name = "percIsos")

pdf(file.path(plotDir,"CDFplot_broadCategs_meanVar_ComplexQuadrant.pdf"),8,6,useDingbats = FALSE)
ggplot(cumuDF2 %>% filter(Broad != "Vasc"), 
       aes(x = meanVar, y = percIsos, color = Broad)) +
  geom_line(size = 1.2) +
  theme_classic() +
  scale_color_met_d("Hokusai2")
dev.off()

## Get inter-quadrant genes and represent ------
numQuadrantsPerGeneAndCT <- normalizedData %>% 
  group_by(Gene,Celltype) %>% 
  summarize(numQuadrants = n_distinct(Status),.groups = "keep") %>% 
  ungroup() %>% select(-Gene)

numQuadrantsPerGene <- normalizedData %>% filter(Celltype == ct) %>% 
  group_by(Gene) %>% summarize(numQuadrants = n_distinct(Status))

getQuadrantConnections <- function(normalizedData,ct) {
  connections <- combn(levels(normalizedData$Status),2)
  conName <- combn(levels(normalizedData$Status),2, paste, collapse = "-")
  ct_df <- normalizedData %>% filter(Celltype == ct)
  lengthVec <- NULL
  for(i in 1:ncol(connections)){
    tmpDF <- ct_df %>% filter(Status %in% connections[,i])
    genesIn2 <- unique(tmpDF %>% filter(Status %in% connections[,i]) %>% 
                         select(Gene,Status) %>% distinct() %>% 
                         group_by(Gene) %>% add_count() %>% 
                         filter(n == 2) %>% .$Gene)
    lengthVec[[i]] <- length(genesIn2)
  }
  names(lengthVec) <- conName
  lengthVec$nGenes <- length(unique(ct_df$Gene))
  lengthVec$CT <- ct
  return(lengthVec)
}
  
ctQuadConn <- lapply(levels(normalizedData$Celltype), 
                     function(ct) getQuadrantConnections(normalizedData,ct))

ctQuadConnDF <- do.call('rbind',ctQuadConn) %>% as.data.frame()

ctQuadConnFlat <- ctQuadConnDF %>% pivot_longer(cols = contains("-"),
                                                names_to = "Connection") %>%
  mutate_at(c("nGenes","value"),as.integer) %>%
  mutate(Importance = value*100/nGenes) %>% 
  separate(Connection,into = c("from","to"),sep = "-")

## Network diagrams -----

getNetworkDiag <- function(ct){
  links_en <- ctQuadConnFlat %>% filter(CT == ct) %>%
    select(from,to,Importance) %>%
    as.data.frame()
  nodes <- data.frame(name = as.vector(unique(normalizedData$Status)))
  
  
  network_1 <- graph_from_data_frame(d=links_en, vertices=nodes, directed=F)
  
  coul <- met.brewer("Peru1",4)[1:4]
  my_color <- coul[as.numeric(factor(nodes$name, 
                                     levels = levels(normalizedData$Status)))]
  pdf(paste0('Plot_Simplex/networkPlots_simplex/network_',ct,'.pdf'),4,4,useDingbats = FALSE)
  plot(network_1, edge.width=E(network_1)$Importance/10, vertex.color = my_color)
  dev.off()
  return()
}

cts <- as.vector(unique(normalizedData$Celltype))
A = lapply(cts, function(ct) getNetworkDiag(ct))

## Barplot - supp fig --------

barPlot_quadrants <- normalizedData %>% group_by(Status,Celltype) %>% 
  count() %>% ungroup() %>% group_by(Celltype) %>% mutate(Perc = n*100/sum(n))

barPlot_quadrants$Status <- factor(barPlot_quadrants$Status, 
                                   levels = c("Complex","AgeSpec","RegSpec","SubtypeSpec"))

pdf(file.path(plotDir,"AllCT_4Quadrants_barPlot.pdf"),15,15,useDingbats = FALSE)
ggplot(barPlot_quadrants,aes(x = Status, y = Perc, fill = Status)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~Celltype) +
  theme_classic() +
  scale_fill_met_d("Peru1")
dev.off()

## Supp fig example --------

pdf(file.path(plotDir,"AllCT_Rufy3.pdf"),15,4,useDingbats = FALSE)
ggtern(data=gW %>% filter(Gene == "Rufy3"),
       aes(x=Age,y=Region,z=Subtype, label = ID, col = ID),
       aes(x,y,z)) + 
  facet_wrap(~Celltype,nrow = 1) +
  geom_Tline(Tintercept=.5,colour='black') + 
  geom_Lline(Lintercept=.5, colour='black') + 
  geom_Rline(Rintercept=.5, colour='black') +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(  title= "Rufy3")
dev.off()

# Hyper variable genes preprocessing ------

normalizedData <- gW %>% filter(Age >= 0.25 | Region >= 0.25 | Subtype >= 0.25) %>% ## make note
  mutate(sum = (Age+Region+Subtype)) %>% 
  mutate_at(c("Age","Region","Subtype"),~ ./sum) %>% select(-sum) %>%
  mutate(Status = case_when(Age > 0.5 & Region < 0.5 & Subtype < 0.5 ~ "AgeSpec",
                            Age < 0.5 & Region > 0.5 & Subtype < 0.5 ~ "RegSpec",
                            Age < 0.5 & Region < 0.5 & Subtype > 0.5 ~ "SubtypeSpec",
                            TRUE ~ "Complex"))

gW_mod <- right_join(gW,normalizedData %>% select(isoID, Status),by = "isoID")
gW_mod <- gW_mod %>% filter(Status != "Complex" | 
                              (Status == "Complex" & Age >= 0.25 & Region >= 0.25 & Subtype >= 0.25) )
highlyVariableIsos <- normalizedData %>% filter(isoID %in% gW_mod$isoID) %>%
  separate(isoID,into = c("Gene","id"),sep = "-",remove = FALSE,extra = "merge")
contributingGenes <- highlyVariableIsos %>% select(Gene,Status,Celltype)

diversityOfStatus <- contributingGenes %>% distinct() %>% group_by(Celltype,Gene) %>% add_count() %>% 
  filter(Celltype != "Vasc") %>%
  select(-Status) %>% distinct() %>% pivot_wider(names_from = Celltype, values_from = n,values_fill = 0) %>%
  column_to_rownames("Gene") %>% 
  as.matrix()

rS <- rowSums(diversityOfStatus)

## Hyper variable genes plot - supp -----

set.seed(1)

colors = structure(rev(met.brewer("Hokusai3",5)))

mat <- diversityOfStatus
row_ha = rowAnnotation(foo = anno_lines(rS),width = unit(3, "in"))
ht <- Heatmap(mat, right_annotation = row_ha, row_km = 7, 
              show_row_dend = FALSE, show_column_dend = FALSE,
              show_row_names = FALSE, col = colors,use_raster = TRUE)
ht <- draw(ht)

options(repr.plot.width=12, repr.plot.height=12)
pdf(file.path(plotDir,'hyperVariableGenes.pdf'),12,12)
ht
dev.off()

# Bulk ------
## Preprocessing bulk ----
gW_bulk <- fread('../data/genomeWide_bulk.gz', header = T)

gW_bulk <- gW_bulk %>% separate("isoID",c("Gene","ID"),sep = "-",extra = "merge",remove = FALSE) 

normalizedData_bulk <- gW_bulk %>% filter(Age >= 0.1 | Region >= 0.1 | Celltype >= 0.1) %>% 
  mutate(sum = (Age+Region+Celltype)) %>% 
  mutate_at(c("Age","Region","Celltype"),~ ./sum) %>% select(-sum) %>%
  mutate(Status = case_when(Age > 0.5 & Region < 0.5 & Celltype < 0.5 ~ "AgeSpec",
                            Age < 0.5 & Region > 0.5 & Celltype < 0.5 ~ "RegSpec",
                            Age < 0.5 & Region < 0.5 & Celltype > 0.5 ~ "CelltypeSpec",
                            TRUE ~ "Complex"))
gW_bulk <- right_join(gW_bulk, normalizedData_bulk %>% select(isoID,Status))
gW_bulk <- gW_bulk %>% rowwise() %>%
  mutate(avg = case_when(Status == "Complex" ~ mean(c(Age,Region,Celltype)),
                         Status == "CelltypeSpec" ~ Celltype,
                         Status == "RegSpec" ~ Region,
                         Status == "AgeSpec" ~ Age))

## Ternary plot bulk -----
pdf(file.path(plotDir,"bulkData_coloredByAvgVar.pdf"),8,8,useDingbats = FALSE)
ggtern(data=gW_bulk %>% filter(Age >= 0.1 | Region >= 0.1 | Celltype >= 0.1),
       aes(x=Age, y=Region, z=Celltype, col = avg),
       aes(x,y,z)) + 
  geom_point_rast(size = 0.3, alpha = 0.6) +
  theme_arrowdefault() + 
  labs( title= "Isoform variability contribution") +
  geom_Tline(Tintercept=.5,colour='black') + 
  geom_Lline(Lintercept=.5, colour='black') +
  geom_Rline(Rintercept=.5, colour='black') +
  scale_color_distiller(palette = 4,direction="horizontal")
dev.off()

allBound <- fread('../data/flIso_oneRegionVsAll_perCellType.gz')

numSig <- allBound %>%
    mutate(sigStatus = case_when(FDR <= 0.05 & abs(dPI) >= 0.1 ~ "Sig",
                                                              TRUE ~ "NonSig")) %>%
    select(Region,Type,sigStatus) %>% 
    group_by_all() %>% add_count() %>% distinct() %>% 
    ungroup() %>% group_by(Region,Type) %>%
    mutate(perc = n*100/sum(n))

numSig$Region <- factor(numSig$Region, 
                        levels <- c("VisCortex","Hippocampus","Striatum","Thalamus","Cerebellum"))
numSig$Type <- factor(numSig$Type, levels <- c("Astro","Oligo","Immune","InhibNeuron","ExciteNeuron"))

## Plot number of significant genes per brain region -----
options(repr.plot.width=10, repr.plot.height=8)

g1 = ggplot(numSig %>% filter(sigStatus == "Sig"),
      aes(x = Type, y = perc, fill = Region)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    theme_classic(base_size = 20) + 
    scale_fill_manual(values = c("#debbc6","#fdc681","#9c9648","#7e9a66","#50939a")) +
    theme(legend.position = "bottom")

pdf('../Plots/PercSignificantPerBR.pdf',10,8,useDingbats = FALSE)
g1
dev.off()

## Preprocessing for getting the number of overlaps between regions -----
uniqPerRegion <- allBound %>%
    mutate(sigStatus = case_when(FDR <= 0.05 & abs(dPI) >= 0.1 ~ "Sig",
                                                              TRUE ~ "NonSig")) %>%
      filter(sigStatus == "Sig") %>% ungroup() %>%
      select(Gene,Type) %>%
      group_by_all() %>% add_count() %>% distinct()

numOverlaps <- table(uniqPerRegion[,c(2:3)]) %>% as.data.frame() %>% 
    group_by(Type) %>% mutate(Perc = Freq/sum(Freq))

numOverlaps$Type <- factor(numOverlaps$Type, 
                           levels <- c("Astro","Oligo","Immune","ExciteNeuron","InhibNeuron"))
numOverlaps$n <- factor(numOverlaps$n,5:1)

### Plot number of overlaps -----
g3 <- ggplot(numOverlaps, aes(x = Type, y = Perc, fill = n)) +
    geom_bar(stat = "identity", width = 0.7) +
    theme_classic(base_size = 20) +
    scale_fill_brewer(name = "# of regions", palette = "Greys",direction = -1) +
    ylab("% of sig genes")
    
pdf('../Plots/numRegOverlaps_perCT.pdf',10,8,useDingbats = FALSE)
g3
dev.off()

## Cumulative matrix pre-processing ------
cumuMat <- NULL
for (i in seq(0,1,length.out = 11)){
    cumuMat[[as.character(i)]] <- allBound %>% select(Gene,FDR,dPI,Region,Type) %>% 
     mutate(sigStatus = case_when(FDR <= 0.05 & abs(dPI) >= i ~ "Sig",
                                  TRUE ~ "NonSig")) %>%
     select(Region,Type,sigStatus) %>% group_by_all() %>%
     add_count() %>% distinct() %>% ungroup() %>%
     group_by(Region,Type) %>% mutate(perc = n/sum(n), cutoff = i)
}

cumuDF <- do.call('rbind',cumuMat)

cumuDF$Region <- factor(cumuDF$Region, 
                        levels <- c("VisCortex","Hippocampus","Striatum","Thalamus","Cerebellum"))
cumuDF$Type <- factor(cumuDF$Type, levels <- c("Astro","Oligo","Immune","InhibNeuron","ExciteNeuron"))

### Plot cumulative stats ------
options(repr.plot.width=16, repr.plot.height=5)

cp <- ggplot(cumuDF %>% filter(sigStatus == "Sig" & cutoff >0 & cutoff <1),
       aes(x = cutoff,y = perc, col = Region)) +
    geom_point(size = 3) + geom_line(size =1) +
    facet_wrap(~Type,nrow = 1) +
    scale_x_reverse(breaks = seq(1,0,length.out = 11)) + 
    theme_classic(base_size = 20) +
    scale_color_manual(values = c("#debbc7","#fdc681","#9c9648","#7e9a66","#50939a")) +
    ylab("Signif genes (%)") + theme(axis.text.x = element_text(angle = 90),
                                     legend.position = "bottom")
                                     
pdf('../Plots/cumulativePlots_perCT.pdf',16,5,useDingbats = FALSE)
cp
dev.off()

## Rug plot of significant genes x cell type per region -----
drawRugPlot <- function(ct){
    wideDF <- allBound %>% mutate(sigStatus = case_when(FDR <= 0.05 & abs(dPI) >= 0.1 ~ 1,
                                                                  TRUE ~ 0.5)) %>%
         select(Gene,Region,Type,sigStatus) %>% filter(Type == ct) %>%
         pivot_wider(names_from = Region, values_from = sigStatus) %>% 
        replace(is.na(.),0)

    wideDF <- wideDF %>% filter_all(any_vars(. == 1))

    wideDF2 <- scale(wideDF[,-c(1:2)])
    ord <- hclust( dist(wideDF2, method = "euclidean"), method = "ward.D" )$order
    ord2 <- hclust( dist(t(wideDF2), method = "euclidean"), method = "ward.D" )$order

    ctDF <- reshape2::melt(wideDF,c("Gene","Type"))
    colnames(ctDF)[3:4] <- c("Region","SigStatus")

    ctDF$SigStatus <- plyr::mapvalues(ctDF$SigStatus,from = c(0.5,0,1),to = c("NonSig","Untested","Sig"))

    ctDF$Gene <- factor(ctDF$Gene, levels = wideDF$Gene[ord])
    ctDF$Region <- factor(ctDF$Region, levels = colnames(wideDF2)[ord2])

    options(repr.plot.width=8, repr.plot.height=6)

    ro <- ggplot(ctDF,aes(x = Region,y = Gene, fill = SigStatus)) +
        geom_tile() +
        theme_classic(base_size = 20) +
        theme(axis.text.y = element_blank()) +
        scale_fill_manual(values = c("#1D1E24","#921A1D","#FCEEDF"))

    pdf(paste0('../Plots/rugPlot_',ct,'.pdf'),8,6,useDingbats = FALSE)
    print(ro)
    dev.off()
    
    return()

}

### Rug plot construction -----
cts <- c("Astro","Oligo","Immune","InhibNeuron","ExciteNeuron")
rugPlots <- lapply(cts, function(ct) drawRugPlot(ct))


# TSS and PolyA site representation -------

## Read in tss-polya site usage info -----
allBound <- fread('../data/endSites_oneRegionVsAll_perCellType.gz')

cumuMat <- NULL
for (i in seq(0,1,length.out = 11)){
    cumuMat[[as.character(i)]] <- allBound %>% select(Gene,FDR,dPI,Region,Type,endSite) %>% 
     mutate(sigStatus = case_when(FDR <= 0.05 & abs(dPI) >= i ~ "Sig",
                                  TRUE ~ "NonSig")) %>%
     select(Region,Type,endSite,sigStatus) %>% group_by_all() %>%
     add_count() %>% distinct() %>% ungroup() %>%
     group_by(Region,Type,endSite) %>% mutate(perc = n/sum(n), cutoff = i)
}

cumuDF <- do.call('rbind',cumuMat)

cumuDF$Region <- factor(cumuDF$Region, 
                        levels <- c("VisCortex","Hippocampus","Striatum","Thalamus","Cerebellum"))
cumuDF$Type <- factor(cumuDF$Type, levels <- c("Astro","Oligo","Immune","InhibNeuron","ExciteNeuron"))
cumuDF$endSite <- factor(cumuDF$endSite, levels <- c("TSS","PolyA"))

## Plot TSS and PolyA site differences per region -----
options(repr.plot.width=12, repr.plot.height=6)

pdf('../Plots/barPlot_numSig_tssPolyA.pdf',12,6,useDingbats = FALSE)
ggplot(cumuDF %>% filter(sigStatus == "Sig" & cutoff == 0.1),
      aes(x = Region, y = perc*100, fill = interaction(endSite,Region) )) +
    facet_wrap(~Type,nrow = 1) +
    geom_bar(stat = "identity", position = "dodge") + 
    theme_classic(base_size = 18) +
    theme(axis.text.x = element_text(angle = 90)) +
    ylab("Significant genes (%)") + 
    scale_fill_manual(values = c("#debbc7","#AF8494","#fdc681","#DEA05E",
                                  "#f68d61","#685E27","#7e9a66","#4E7134","#50939a","#226C70"))
dev.off()
           


# SessionInfo ------

# R version 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] igraph_1.2.6    ggsignif_0.6.4  ggrastr_1.0.1   MetBrewer_0.2.0 ggtern_3.3.5    ggplot2_3.3.6  
# [7] tidyr_1.2.0     dplyr_1.0.9  circlize_0.4.15       ComplexHeatmap_2.13.1    cowplot_1.1.1 data.table_1.13.6
