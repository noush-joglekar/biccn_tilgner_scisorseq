## Setup
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(parallel)

workingDir <- '../'
if(!dir.exists(workingDir)){dir.create(workingDir)}
setwd(workingDir)

print("Reading in data")

fullDataset <- '../data/LR_annotTranscripts_byCT/'
countData <- Read10X(fullDataset)
dim(countData)

print("Processing all genes")

args <- commandArgs(trailing = TRUE)
threshold <- as.integer(args[1])
sampleIx <- args[2]
print(paste0("Sample:", sampleIx))

CN <- data.frame(CT = colnames(countData)) %>% 
    tidyr::separate(CT, into = c("Age","Sample","Region","Broad","Type","Subtype"),remove = FALSE)


regions <- unique(CN$Region)
ages <- unique(CN$Age)
cts <- unique(CN$Type)

isoInfo <- data.frame(Isoform=rownames(countData))
isoInfo <- isoInfo %>% separate(Isoform, into = c("Gene","ID"), sep = "-", extra = "merge",remove = FALSE)

isoCountsPerGene <- isoInfo %>% group_by(Gene) %>% add_count(name = "numPerGene")

multPerGene <- isoCountsPerGene %>% filter(numPerGene > 1) %>% ungroup()

dsVec <- function(vec,size){
    vec = sapply(1:length(vec), function(i) round(vec[i],0))
    if(sum(vec) > size){
        prob = size/sum(vec)
        splitUnif = split(runif(sum(vec),0,1),rep(1:length(vec),vec) )
        dS = sapply(1:length(splitUnif), function(i) sum(splitUnif[[i]] < prob))
        newVec = vec
        newVec[newVec >= 1] = dS
    } else {newVec = vec}
    return(newVec)
}

downsample <- function(mat,size){
    tmpMat <- lapply(1:ncol(mat), function(c) dsVec(mat[,c],size))
    newMat <- do.call('cbind',tmpMat)
    colnames(newMat) <- colnames(mat)
    return(newMat)
}

getRegSpecVector <- function(reg,df,ct,threshold){
    filt <- CN %>% filter(Region == reg & Type == ct) %>% select(CT)
    mat <- as.matrix(df[,filt$CT])
    if(sum(mat) >= threshold){
        dsMat = downsample(mat,threshold)
        pi <- rowSums(dsMat)/sum(dsMat)
        return(pi)
    } else {return(NULL)}
}

getAgeSpecVector <- function(age,df,ct,threshold){
    filt <- CN %>% filter(Age == age & Type == ct) %>% select(CT)
    mat <- as.matrix(df[,filt$CT])
    if(sum(mat) >= threshold){
        dsMat = downsample(mat,threshold)
        pi <- rowSums(dsMat)/sum(dsMat)
        return(pi)
    } else {return(NULL)}
}

getSubtypeSpecVector <- function(df,ct,threshold){
    filt <- CN %>% filter(Type == ct)
    st <- unique(filt$Subtype)
    A = lapply(st, function(i) df[,filt$CT] %>% as.data.frame() %>% 
        mutate(a = rowSums(across(filt %>% filter(Subtype == i) %>% .$CT))) %>% .$a)
    B = as.data.frame(do.call('cbind',A))
    rownames(B) <- rownames(df)
    colnames(B) <- st
    B <- B[,which(colSums(B) >= threshold)]
    if(is.data.frame(B) && ncol(B) > 1){
        dsMat = as.data.frame(downsample(B,threshold))
        return(dsMat)
    } else{return(B)}
}

getVarEstimatesPerGene <- function(gene,ct,threshold){
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
    tmpReg <- lapply(regions, function(reg) getRegSpecVector(reg,df,ct,threshold))
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
    tmpAge <- lapply(ages, function(age) getAgeSpecVector(age,df,ct,threshold))
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
    mat <- getSubtypeSpecVector(df,ct,threshold)
    
    if(is.data.frame(mat)){
    if(ncol(mat) >=2 && sum(mat) != 0){
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


threshold = 10

genomeWide_perGene <- mclapply(unique(multPerGene$Gene),function(gene) do.call("rbind",
	lapply(cts,function(ct) {
#		print(gene)
		getVarEstimatesPerGene(gene,ct,threshold)})),mc.cores = 6)

rbd = do.call('rbind',genomeWide_perGene)
if (!is.null(rbd)){
	gW <- as.data.frame(rbd %>% drop_na())
	write.table(gW,paste0('genomeWide_withContribs_threshold_10_subsample_10_ix_',sampleIx),
		sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
} else{print("No output")}
