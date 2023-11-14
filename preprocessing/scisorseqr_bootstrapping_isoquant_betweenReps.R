## Setup
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(ggplot2)

workingDir <- '../bootstrapping/'
setwd(workingDir)

fullDataset <- '../data//LR_annotTranscripts_byCT/'

print("Reading input -----")
countData <- Read10X(fullDataset)
dim(countData)

args <- commandArgs(trailing = TRUE)
seed = as.integer(args[1])
set.seed(seed)

CN <- data.frame(CT = colnames(countData)) %>% 
    tidyr::separate(CT, into = c("Age","Sample","Region","Broad","Type","Subtype"),remove = FALSE)

CN2 = CN %>% mutate(Rep = case_when(Sample == "M1" ~ "Rep1",
                                   Sample == "M5" & Age == "P21" ~ "Rep2",
                                   Sample == "M5" & Age != "P21" ~ "Rep1", 
                                   Sample == "M8" ~ "Rep1",
                                   TRUE ~ "Rep2"))

CN = CN2

## functions:

getCorr <- function(mat){
    evalMat <- sweep(mat[,4:5],2,colSums(mat[,4:5]),`/`)
    if(nrow(evalMat) > 2){
        cv = cor(evalMat[,1],evalMat[,2])
    } else {
        cv = as.vector(lm(evalMat[,1] ~ evalMat[,2])$coef[2])
    }
    return(c(cv,nrow(evalMat)))
}

deltaPI <- function(mat){
  pos <- sum(c(rev(mat$delta[mat$delta > 0]),0,0)[1:2])
  neg <- abs(sum(c(mat$delta[mat$delta < 0],0,0)[1:2]))
  change=max(pos,neg)
  if (pos >= neg){
    x = c(rev(mat$delta[mat$delta > 0]),0)[1:2]
    y = mat %>% dplyr::filter(delta %in% x) %>% dplyr::select(ID) %>% as.matrix()
    return(c(change,y[1],y[2]))}
  else {
    x = c(mat$delta[mat$delta < 0],0,0)[1:2]
    y = mat %>% dplyr::filter(delta %in% x) %>% dplyr::select(ID) %>% as.matrix()
    return(c(-change,y[1],y[2]))}
}

Get_Pval_DeltaPI <- function(mat){
  cc = getCorr(mat)
  cv = cc[1]
  numIso = cc[2]
  dPI_calc <- deltaPI(mat)
  max_ix1 <- dPI_calc[2]
  max_ix2 <- dPI_calc[3]
  pval <- chisq.test(mat[,4:5])$p.value
  d_pi <- dPI_calc[1]
  return(list(pval,d_pi,max_ix1,max_ix2,cv,numIso))
}


dsVec <- function(vec,size){
    vec = sapply(1:length(vec), function(i) round(vec[i],0))
    prob = size/sum(vec)
    splitUnif = split(runif(sum(vec),0,1),rep(1:length(vec),vec) )
    dS = sapply(1:length(splitUnif), function(i) sum(splitUnif[[i]] < prob))
    newVec = vec
    newVec[newVec >= 1] = dS
    return(newVec)
}

downsample <- function(mat,size){
    newMat <- mat
    newMat[,4] <- dsVec(mat[,4],size)
    newMat[,5] <- dsVec(mat[,5],size)
    return(newMat)
}
                    
downsampleAndCalcSig <- function(mat,size){
    newMat <- downsample(mat,size)
    newMat <- newMat %>% rowwise() %>% 
        filter(Rep1 > 0 & Rep2 > 0) %>% as.data.frame()
    if(nrow(newMat) >= 2){
        newMat <- newMat %>% mutate(pi1 = Rep1/sum(Rep1),
            pi2 = Rep2/sum(Rep2), delta = pi1-pi2) %>%
            arrange(delta, .by_group = TRUE) %>% as.data.frame()
        outMat <- Get_Pval_DeltaPI(newMat)
        return(outMat)
    } else return()
}

oneRegion_betweenReps <- function(REG,ctoi,size,numCores){
    
    cols <- CN %>% filter(Age == "P56" & Type == ctoi & Region == REG) %>% select(CT)
        
    subsetDF <- left_join(as.data.frame(countData[,cols$CT]) %>% mutate(iso = rownames(countData)) %>%
        separate(iso, into = c("Gene","ID"), sep = "-", extra = "merge") %>%
        pivot_longer(cols = contains("::"),names_to = "CT",values_to = "counts"), 
                    CN, by = "CT")
    
    processedDF <- subsetDF %>% group_by(Gene, !!sym("Type"), ID, Rep) %>%
            summarise(Sum = sum(counts)) %>%
            filter(sum(Sum) >= threshold ) %>%
            ungroup() %>% group_by(Gene) %>% 
            filter(length(unique(Rep)) == 2) %>%
            pivot_wider(names_from = "Rep", values_from = "Sum") %>% 
            replace(is.na(.), 0) %>% 
            filter(Rep1 > 0 & Rep2 > 0) %>% 
            filter(n() >= 2) %>% as.data.frame()

    if(nrow(processedDF) == 0){return (NULL)}
    
    PerGene = split(processedDF,processedDF$Gene)
    
    uncorrected_output <- parallel::mclapply(names(PerGene),
                                         function(geneName) 
                                             downsampleAndCalcSig(PerGene[[geneName]],size),
                                         mc.cores=numCores)
    
    if(nrow(processedDF) == 0){return (NULL)}
    names(uncorrected_output) <- names(PerGene)
    output_DF <- as.data.frame(do.call(rbind, uncorrected_output))
    colnames(output_DF) <- c("pvals","dPI","maxDeltaPI_ix1","maxDeltaPI_ix2","corr","numIso")
    output_DF <- output_DF %>% rownames_to_column("Gene")
    output_DF$FDR <- p.adjust(output_DF$pvals, method = "BH")
    output_DF$dPI <- as.numeric(output_DF$dPI)
    output_DF$Region <- REG
    output_DF$Type <- ctoi

    cat("Finished processing",REG,"-",ctoi,"\n")
    cat("Total genes tested:",nrow(output_DF),"\n")
    cat("Number of significant genes by fdr:",length(which(output_DF$FDR <= 0.05)),"\n")
    cat("Number of fdr genes with deltaPI >= 0.2:",length(which(output_DF$FDR <= 0.05 & abs(output_DF$dPI) >= 0.2)),"\n")
    sigs <- output_DF[which(output_DF$FDR <= 0.05 & abs(output_DF$dPI) >= 0.1),]
    cat(dim(sigs))
    
    return(output_DF) 
}


getSummary <- function(df,ct,reg){
    ctReg <- df %>% filter(Type == ct & Region == reg)
    
    summSt = NULL

    for (thresh in c(0.05,0.1,0.2,0.3,0.4,0.5)){
            sig = ctReg %>% filter(abs(dPI)> thresh & FDR <= 0.05) %>% .$Gene
            tot = nrow(ctReg)
            summSt[[as.character(thresh)]] <- c(percSig = length(sig)*100/tot)
        }

    df = as.data.frame(do.call('rbind',summSt)) %>% 
        mutate(Region = reg, CT = ct) %>% tibble::rownames_to_column("thresh")
    return(df)
}


regions <- c('Hippocampus','VisCortex','Cerebellum','Striatum','Thalamus')
cts <- c("Oligo","Astro","Immune","ExciteNeuron","InhibNeuron")

## User defined parameters
threshold <- 100 ## threshold for number of reads per sample
size <- 50 ## threshold for number of reads to subsample
numCores <- 8 ## number of cores

print("Running rep 1 vs rep 2  -------")
allCT_DF <- do.call('rbind',lapply(cts, function(ct) do.call('rbind',
                                             lapply(regions, function(reg) oneRegion_betweenReps(reg,ct,size,numCores)))))


DF = do.call('rbind',lapply(regions, function(reg)
    do.call('rbind',lapply(cts, function(ct) getSummary(allCT_DF,ct,reg)))))


print("Writing output 1 ------")
outdir_betn = "bootstrapping_betweenReps_10050/"
if(!dir.exists(outdir_betn)){dir.create(outdir_betn)}
write.table(DF, paste0(outdir_betn,"summaryPerThresh_betweenReps_",seed,".tab"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
