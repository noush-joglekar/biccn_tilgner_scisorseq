#! /bin/R
# UMAPs and QC from short read scRNA data

# Set up ----------

library(Seurat)
library(patchwork)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(harmony)

## Load in integrated data -----

load('../data/fullObject_harmonyIntegrated.Robj') ## Seurat object too big. Individual objects for all samples deposited on NeMO

## Read in color palette -----
cP <- read.csv('data/colorPalette.csv',header = FALSE)
colnames(cP) <- c("color","Subtype")

my_pal <- cP$color
names(my_pal) <- cP$Subtype

## Fig 1a ------

fig1a <- DimPlot(harmony_merged, group.by = "Subtype", label = T, cols = my_pal) +
    NoLegend()

pdf('../Plots/Subtype_UMAP.pdf',12,12,useDingbats = FALSE)
fig1a
dev.off()

## Pre-process for Fig 1b ------

h <- subset(harmony_merged,subset = Region == "Hippocampus")
v <- subset(harmony_merged,subset = Region == "VisCortex")
p56 <- subset(harmony_merged,subset = Age == "P56")
p56$Region <- factor(p56$Region, levels = rev(c("VisCortex","Hippocampus","Striatum","Thalamus","Cerebellum")))

hplot <- DimPlot(object = subset(h, Age != "P56"), reduction = "umap", split.by = "Age", label = F, 
        group.by = "Age",cols = met.brewer("Morgenstern",8)[8:5], pt.size = 2,
        ncol = 1) + theme_void() +
  theme(legend.position = "none")

vplot <- DimPlot(object = subset(v, Age != "P56"), reduction = "umap", split.by = "Age", label = F,
                 group.by = "Age",cols = met.brewer("Morgenstern",8)[1:4], pt.size = 2,
                 ncol = 1) + theme_void() + 
  theme(legend.position = "none") 

adultPlot <- DimPlot(object = p56, reduction = "umap", split.by = "Region", label = F, group.by = "Region",
        cols = rev(c("#ddbbc6","#fdc681","#9c9648","#7e9a66","#50939a")), pt.size = 2) +theme_void() + 
    theme(legend.position = "none")

## Fig 1b ------

options(repr.plot.width=15, repr.plot.height=15)

fig1b <- (plot_spacer() | plot_spacer() | plot_spacer() | hplot | vplot ) / (adultPlot) +
    plot_layout(heights = unit(c(8.5, 3), c('inches', 'inches')))

pdf('../Plots/spaceAndTime.pdf',15,15,useDingbats = FALSE)
fig1b
dev.off()

## Preprocess for Fig 1c -----

percCells <- harmony_merged@meta.data %>% dplyr::filter(Type != "Doublets") %>% 
        select(Age,Region,Replicate,Subtype) %>% group_by_all() %>% summarize(n = n()) %>%
        mutate(perc = n*100/sum(n))

percCells <- percCells %>% unite(group, c("Age","Region","Replicate"), sep = ":", remove = FALSE) %>%
    unite(CG, c("Age","Region"), sep = "_")
percCells$ColorMap <- numCells$CG
percCells$ColorMap <- factor(numCells$ColorMap,levels <- unique(numCells$ColorMap)[c(1:6,8,11,9,10,7)] )
percCells$ColorMap <- plyr::mapvalues(numCells$ColorMap, from = unique(numCells$ColorMap)[c(1:6,8,11,9,10,7)], 
  to = c(met.brewer("Morgenstern",8)[c(8,1,7,2,6,3,5,4)],"#9c9648","#7e9a66","#50939a"))

uniqCols <- percCells %>% ungroup() %>% select(group,ColorMap) %>% distinct()
uniqCols$ColorMap <- as.character(uniqCols$ColorMap)
colorPal <- uniqCols$ColorMap
names(colorPal) <- uniqCols$group

percCells$group <- factor(percCells$group, levels <- c(grep("Vis",unique(percCells$group), value = T),
                                                     grep("Hipp",unique(percCells$group), value = T),
                                                     grep("Stri",unique(percCells$group), value = T),
                                                     grep("Thal",unique(percCells$group), value = T),
                                                     grep("Cereb",unique(percCells$group), value = T)))

percCells$Subtype <- factor(percCells$Subtype, levels = cP$Subtype)

newLabels <- gsub("Cerebellum","CEREB",
    gsub("Thalamus","THAL",
         gsub("Striatum","STRI",
              gsub("Hippocampus","HIPP",
                   gsub(pattern = "VisCortex","VIS",levels(numCells$group))))))


## Fig 1c --------

options(repr.plot.width=15, repr.plot.height=12)
fig1c = ggplot(percCells, aes(x = group, y = Subtype, col = group)) +
    geom_point(aes(size = perc)) + 
    scale_size_continuous(range = c(1, 8)) +
    theme_classic(base_size = 20) + 
    theme(axis.text.x = element_text(angle = 90)) +
    scale_color_manual(values = colorPal, guide = "none") +
    scale_x_discrete(labels= newLabels)

pdf('../Plots/percCells_subType.pdf',15,12,useDingbats = FALSE)
fig1c
dev.off()


## Preprocessing for fig 1d -------
cpS <- harmony_merged@meta.data %>% 
        select(Age,Region,Replicate) %>% group_by_all() %>% summarize(n = n())

cpS <- cpS %>% unite(group, c("Age","Region","Replicate"), sep = ":", remove = FALSE) %>%
    unite(CG, c("Age","Region"), sep = "_")
cpS$ColorMap <- cpS$CG
cpS$ColorMap <- factor(cpS$ColorMap,levels <- unique(cpS$ColorMap)[c(1:6,8,11,9,10,7)] )
cpS$ColorMap <- plyr::mapvalues(cpS$ColorMap, from = unique(cpS$ColorMap)[c(1:6,8,11,9,10,7)], 
  to = c(met.brewer("Morgenstern",8)[c(8,1,7,2,6,3,5,4)],"#9c9648","#7e9a66","#50939a"))

cpS$group <- factor(cpS$group, levels <- c(grep("Vis",unique(cpS$group), value = T),
                                                     grep("Hipp",unique(cpS$group), value = T),
                                                     grep("Stri",unique(cpS$group), value = T),
                                                     grep("Thal",unique(cpS$group), value = T),
                                                     grep("Cereb",unique(cpS$group), value = T)))

uniqCols <- numCells %>% ungroup() %>% select(group,ColorMap) %>% distinct()
uniqCols$ColorMap <- as.character(uniqCols$ColorMap)
colorPal <- uniqCols$ColorMap
names(colorPal) <- uniqCols$group

## Fig 1d --------
options(repr.plot.width=10, repr.plot.height=12)

fig1d <- ggplot(cpS, aes(x = group, y = n, fill = group)) +
    geom_bar(stat = "identity") + 
    theme_classic(base_size = 20) + 
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values = colorPal, guide = "none") +
    xlab("Sample") + ylab("Number of cells") +
    scale_x_discrete(labels= newLabels)

pdf('../Plots/numberOfCells.pdf',10,12,useDingbats = FALSE)
fig1d
dev.off()
