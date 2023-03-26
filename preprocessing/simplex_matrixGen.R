# Setup ------
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(parallel)

## Read in and preprocess ------
fullDataset <- '../data/LR_annotTranscripts_byCT/'
countData <- Read10X(fullDataset)
dim(countData)

CN <- data.frame(CT = colnames(countData)) %>% 
    tidyr::separate(CT, into = c("Age","Sample","Region","Broad","Type","Subtype"),remove = FALSE)

regions <- unique(CN$Region)
ages <- unique(CN$Age)
cts <- unique(CN$Type)

isoInfo <- data.frame(Isoform=rownames(countData))
isoInfo <- isoInfo %>% separate(Isoform, into = c("Gene","ID"), sep = "-", extra = "merge",remove = FALSE)
isoCountsPerGene <- isoInfo %>% group_by(Gene) %>% add_count(name = "numPerGene")
multPerGene <- isoCountsPerGene %>% filter(numPerGene > 1) %>% ungroup()

# Bulk -------
## Cell type functions -------
getRegSpecVector <- function(reg,df,ct){
    filt <- CN %>% filter(Region == reg & Type == ct) %>% select(CT)
    mat <- as.matrix(df[,filt$CT])
    pi <- rowSums(mat)/sum(mat)
    return(pi)
}

getAgeSpecVector <- function(age,df,ct){
    filt <- CN %>% filter(Age == age & Type == ct) %>% select(CT)
    mat <- as.matrix(df[,filt$CT])
    pi <- rowSums(mat)/sum(mat)
    return(pi)
}


getSubtypeSpecVector <- function(df,ct){
    filt <- CN %>% filter(Type == ct)
    st <- unique(filt$Subtype)
    A = lapply(st, function(i) df[,filt$CT] %>% as.data.frame() %>% 
        mutate(a = rowSums(across(filt %>% filter(Subtype == i) %>% .$CT))) %>% .$a)
    B = as.data.frame(do.call('cbind',A))
    rownames(B) <- rownames(df)
    colnames(B) <- st
    #B <- B[,which(colSums(B) >= 10)] ## no need for threshold
    return(B)
}

getVarEstimatesPerGene <- function(gene,ct){
    regions <- unique(CN$Region)
    ages <- unique(CN$Age)
    cts <- unique(CN$Type)
    
    dummy <- multPerGene %>% filter(Gene == gene) %>% select(Isoform)
    df <- countData[dummy$Isoform,]
    filt <- CN %>% filter(Type == ct) %>% select(CT)
    df <- df[,filt$CT]
    
    ## bulk
    
    mat <- as.matrix(df)
    suff <- names(which(colSums(mat) > 10))
    if(length(suff) >= 3){
        bulkDF <- sweep(mat,2,colSums(mat),`/`) %>% as.data.frame() %>%
                rowwise() %>%
                mutate(bulkVar = max(c_across(all_of(suff)),na.rm = T) - 
                        min(c_across(all_of(suff)),na.rm = T)) %>%
                select(bulkVar)
    } else {bulkDF = NULL}
    
    ## region sp
    tmpReg <- lapply(regions, function(reg) getRegSpecVector(reg,df,ct))
    names(tmpReg) <- regions
    regDF <- do.call('cbind',tmpReg) %>% as.data.frame() %>% select(where(function(x) any(!is.na(x))))
    if(ncol(regDF) >=2){
        cn <- colnames(regDF)
        regVec <- regDF %>% mutate(M = pmax(!!!rlang::syms(cn)),
                                   m = pmin(!!!rlang::syms(cn)),
                                  regVar = M -m, 
                                  med = median(c_across(cn),na.rm = T),
                                  MID = which.max(c_across(cols = cn)),
                                  mID = which.min(c_across(cols = cn))) %>% 
        rowwise() %>%
        mutate(name = case_when((M - med) > (med - m) ~ colnames(mat)[MID],
                               TRUE ~ colnames(mat)[mID])) %>%
        separate(name, into = c("a","s","regName","b","t","st"), sep = "::") %>%
        select(-c(M,m,med,MID,mID,a,s,b,t,st)) %>%
        as.data.frame() %>% select(regVar,regName)
    } else {regVec = NULL}
                                 
    ## age sp
    tmpAge <- lapply(ages, function(age) getAgeSpecVector(age,df,ct))
    names(tmpAge) <- ages
    ageDF <- do.call('cbind',tmpAge) %>% as.data.frame() %>% select(where(function(x) any(!is.na(x))))     
    if(ncol(ageDF) >=2){
        cn <- colnames(ageDF)
        ageVec <- ageDF %>% mutate(M = pmax(!!!rlang::syms(cn)),
                                   m = pmin(!!!rlang::syms(cn)),
                                  ageVar = M -m, 
                                  med = median(c_across(cn),na.rm = T),
                                  MID = which.max(c_across(cols = cn)),
                                  mID = which.min(c_across(cols = cn))) %>% 
        rowwise() %>%
        mutate(name = case_when((M - med) > (med - m) ~ colnames(mat)[MID],
                               TRUE ~ colnames(mat)[mID])) %>%
        separate(name, into = c("ageName","s","r","b","t","st"), sep = "::") %>%
        select(-c(M,m,med,MID,mID,s,r,b,t,st)) %>%
        as.data.frame() %>%
        select(ageVar,ageName)
        rownames(ageVec) <- rownames(df)
    } else {ageVec = NULL}
                                                                          
    ## cell subtype sp
    mat <- getSubtypeSpecVector(df,ct)
    
    if(is.data.frame(mat)){
    if(ncol(mat) >=2 & sum(mat) != 0){
        stVec <- sweep(mat,2,colSums(mat),`/`) %>% as.data.frame() %>%
                rowwise() %>%
                mutate(M = max(c_across(cols = colnames(mat)),na.rm = T),
                    m = min(c_across(cols = colnames(mat)),na.rm = T),
                    stVar = M - m, 
                       med = median(c_across(cols = colnames(mat)),na.rm = T),
                      MID = which.max(c_across(cols = colnames(mat))),
                      mID = which.min(c_across(cols = colnames(mat)))) %>%
        rowwise() %>%
        mutate(stName = case_when((M - med) > (med - m) ~ colnames(mat)[MID],
                               TRUE ~ colnames(mat)[mID])) %>%
        select(-c(M,m,med,MID,mID)) %>%
        as.data.frame() %>% select(stVar,stName)
        rownames(stVec) <- rownames(df)
    } else {stVec = NULL}
    } else {stVec = NULL}

    ## combine
    if(is.null(regVec) | is.null(ageVec) | is.null(stVec) | is.null(bulkDF) ){
        return()
    } else {
        fullDF <- do.call('cbind',list(bulkDF$bulkVar,ageVec,regVec,stVec))
        colnames(fullDF) <- c("bulk","Age","AgeName","Region","RegName","Subtype","SubtypeName")
        fullDF <- fullDF[rowSums(fullDF[,c(1,2,4,6)])>0,]
        fullDF <- fullDF %>% rownames_to_column("isoID") %>%
                rowwise() %>% mutate(missingVar = abs(bulk - (Age + Region + Subtype)))
        fullDF$Celltype <- ct
        return(fullDF)
    }
}

## Cell type execution -------
genomeWide_perGene <- lapply(unique(multPerGene$Gene),function(gene) do.call('rbind',lapply(cts, 
                                                           function(ct) getVarEstimatesPerGene(gene,ct))))

gW <- as.data.frame(do.call('rbind',genomeWide_perGene) %>% drop_na())

write.table(gW,'../data/genomeWide_withContribs',sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


## Bulk functions -------

getRegSpecVector <- function(reg,df){
    filt <- CN %>% filter(Region == reg) %>% select(CT)
    mat <- as.matrix(df[,filt$CT])
    pi <- rowSums(mat)/sum(mat)
    return(pi)
}

getAgeSpecVector <- function(age,df){
    filt <- CN %>% filter(Age == age) %>% select(CT)
    mat <- as.matrix(df[,filt$CT])
    pi <- rowSums(mat)/sum(mat)
    return(pi)
}


getCelltypeSpecVector <- function(ct,df){
    filt <- CN %>% filter(Type == ct) %>% select(CT)
    mat <- as.matrix(df[,filt$CT])
    pi <- rowSums(mat)/sum(mat)
    return(pi)
}

getVarEstimatesPerGene <- function(gene){
    regions <- unique(CN$Region)
    ages <- unique(CN$Age)
    
    dummy <- multPerGene %>% filter(Gene == gene) %>% select(Isoform)
    df <- countData[dummy$Isoform,]
    filt <- CN %>% select(CT)
    df <- df[,filt$CT]
    
    ## region sp
    tmpReg <- lapply(regions, function(reg) getRegSpecVector(reg,df))
    names(tmpReg) <- regions
    regDF <- do.call('cbind',tmpReg) %>% as.data.frame() %>% select(where(function(x) any(!is.na(x))))
    if(ncol(regDF) >=2){
        cn <- colnames(regDF)
        regVec <- regDF %>% mutate(M = pmax(!!!rlang::syms(cn)),
                                   m = pmin(!!!rlang::syms(cn)),
                                  regVar = M -m) %>%
                rowwise() %>% mutate(med = median(c_across(all_of(cn)),na.rm = T),
                                  MID = which.max(c_across(cols = cn)),
                                  mID = which.min(c_across(cols = cn))) %>% 
        mutate(regName = case_when((M - med) > (med - m) ~ regions[MID],
                               TRUE ~ cn[mID])) %>%
          select(-c(M,m,med,MID,mID)) %>%
          as.data.frame() %>% select(regVar,regName)
    } else {regVec = NULL}
                                                                          
                                 
    ## age sp
    tmpAge <- lapply(ages, function(age) getAgeSpecVector(age,df))
    names(tmpAge) <- ages
    ageDF <- do.call('cbind',tmpAge) %>% as.data.frame() %>% select(where(function(x) any(!is.na(x))))     
    if(ncol(ageDF) >=2){
        cn <- colnames(ageDF)
        ageVec <- ageDF %>% mutate(M = pmax(!!!rlang::syms(cn)),
                                   m = pmin(!!!rlang::syms(cn)),
                                  ageVar = M -m) %>%
                rowwise() %>% mutate(med = median(c_across(all_of(cn)),na.rm = T),
                                  MID = which.max(c_across(cols = cn)),
                                  mID = which.min(c_across(cols = cn))) %>% 
        mutate(x = (M - med), y = (med - m)) %>%
        mutate(ageName = case_when((M - med) > (med - m) ~ cn[MID],
                               TRUE ~ ages[mID])) %>%
        select(-c(M,m,med,MID,mID)) %>%
        as.data.frame() %>%
        select(ageVar,ageName)
        rownames(ageVec) <- rownames(df)
    } else {ageVec = NULL}

                                                                          
    ## cell type sp
    tmpCT <- lapply(cts, function(ct) getCelltypeSpecVector(ct,df))
    names(tmpCT) <- cts
    ctDF <- do.call('cbind',tmpCT) %>% as.data.frame() %>% select(where(function(x) any(!is.na(x))))
    if(ncol(ctDF) >=2){
        cn <- colnames(ctDF)
        ctVec <- ctDF %>% mutate(M = pmax(!!!rlang::syms(cn)),
                                   m = pmin(!!!rlang::syms(cn)),
                                  ctVar = M -m) %>%
                rowwise() %>% mutate(med = median(c_across(all_of(cn)),na.rm = T),
                                  MID = which.max(c_across(cols = cn)),
                                  mID = which.min(c_across(cols = cn))) %>% 
        mutate(ctName = case_when((M - med) > (med - m) ~ cn[MID],
                               TRUE ~ cts[mID])) %>%
          select(-c(M,m,med,MID,mID)) %>%
          as.data.frame() %>% select(ctVar,ctName)
    } else {ctVec = NULL}

                                                                    
    ## combine
    if(is.null(regVec) | is.null(ageVec) | is.null(ctVec) ){
        return()
    } else {
        fullDF <- do.call('cbind',list(ageVec,regVec,ctVec))
        colnames(fullDF) <- c("Age","AgeName","Region","RegName","Celltype","CelltypeName")
        fullDF <- fullDF[rowSums(fullDF[,c(1,3,5)])>0,]
        fullDF <- fullDF %>% rownames_to_column("isoID")
        return(fullDF)
    }
}

# Bulk -------
## Bulk execution -------

genomeWide_perGeneBulk <- lapply(unique(multPerGene$Gene),function(gene) getVarEstimatesPerGene(gene))
gWB <- as.data.frame(do.call('rbind',genomeWide_perGeneBulk) %>% drop_na())
write.table(gWB,'genomeWide_bulk',sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)



# SessionInfo ------
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /pbtech_mounts/homes059/anj2026/miniconda3/envs/hdWGCNA/lib/libopenblasp-r0.3.20.so
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
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#   [1] Seurat_3.2.3 tibble_3.1.8 tidyr_1.2.1  dplyr_1.0.10
