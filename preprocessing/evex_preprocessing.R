#! /bin/R
# Setup --------
library(dplyr)
library(tidyr)
library(tibble)

develDF <- read.table('../data/altExonUsage_devel_psiPerSample_avgOverReps_type',
                      sep = "\t",header = T)
adultDF <- read.table('../data/altExonUsage_P56_psiPerSample_avgOverReps_type',
                      sep = "\t",header = T)


## Identify highly variable exons in development keeping cell type constant ------
regions <- c("Hippocampus","VisCortex")
celltypes <- c("Astro","Oligo","ExciteNeuron","InhibNeuron")
timepoints <- c("P14","P21","P28","P56")

perReg <- list()
perReg_unf <- list()
for (reg in regions){
    ctDF <- list()
    ctDF_unf <- list()
    for (ct in celltypes){
        rct <- paste(reg,ct,sep="_")
        devel_rct <- develDF %>% select(Exon,Gene,contains(rct))
        colnames(devel_rct)[3:6] <- gsub(colnames(devel_rct)[3:6],pattern = paste0("_",rct),replacement = "")
        cols <- combn(timepoints, 2, paste, collapse = "-")

        devel_rct[cols] <- t(apply(devel_rct[3:6], 1, function(x) {
                     out <- combn(x, 2, function(x) x[2] - x[1])
                     c(out, -out)
                    }))  
        
        devel_rct <- devel_rct %>% select(-timepoints) %>% unite("GE",c("Exon","Gene"),sep = "::")
        devel_rct$Region <- reg
        devel_rct$cellType <- ct
        ctDF_unf[[ct]] <- devel_rct
        devel_rct <- devel_rct %>% pivot_longer(colnames(devel_rct)[2:7],
                                                names_to = "Transition",values_to = "dPSI") %>%
                drop_na() %>% filter(abs(dPSI) >= 0.25)
        ctDF[[ct]] <- devel_rct
    }
    perReg[[reg]] <- do.call('rbind',ctDF)
    perReg_unf[[reg]] <- do.call('rbind',ctDF_unf)
}

highVar_devel <- do.call('rbind',perReg) %>% remove_rownames()
highVar_devel$Axis <- "Devel"

highVar_devel_unf <- do.call('rbind',perReg_unf) %>% remove_rownames()
highVar_devel_unf$Axis <- "Devel"

wide_highVar_devel <- highVar_devel %>% pivot_wider(names_from = cellType, values_from = dPSI) 
                                   
                                   
## Timepoint constant, pairwise CT variability ----------
regions <- c("Hippocampus","VisCortex")
celltypes <- c("Astro","Oligo","ExciteNeuron","InhibNeuron")
timepoints <- c("P14","P21","P28")

perReg_tp <- list()
perReg_tp_unf <- list()
for (reg in regions){
    tpDF <- list()
    tpDF_unf <- list()
    for (tp in timepoints){
        tr <- paste(tp,reg,sep= "_")
        devel_tr <- develDF %>% select(Exon,Gene,contains(tr))
        colnames(devel_tr)[3:6] <- gsub(colnames(devel_tr)[3:6],pattern = paste0(tr,"_"),replacement = "")
        cols <- combn(celltypes, 2, paste, collapse = "-")

        devel_tr[cols] <- t(apply(devel_tr[3:6], 1, function(x) {
                     out <- combn(x, 2, function(x) x[1] - x[2])
                     c(out, -out)
                    }))  
        
        devel_tr <- devel_tr %>% select(-all_of(celltypes)) %>% unite("GE",c("Exon","Gene"),sep = "::")
        devel_tr$Region <- reg
        devel_tr$Age <- tp
        tpDF_unf[[tp]] <- devel_tr
        devel_tr <- devel_tr %>% pivot_longer(contains("-"), names_to = "Comparison",values_to = "dPSI") %>% 
            drop_na() %>% filter(abs(dPSI) >= 0.25) 
        tpDF[[tp]] <- devel_tr
    }
    perReg_tp[[reg]] <- do.call('rbind',tpDF)
    perReg_tp_unf[[reg]] <- do.call('rbind',tpDF_unf)
}

highVar_devel_tp <- do.call('rbind',perReg_tp) %>% remove_rownames()
highVar_devel_tp$Axis <- "Devel_CTspec"

highVar_devel_tp_unf <- do.call('rbind',perReg_tp_unf) %>% remove_rownames()
highVar_devel_tp_unf$Axis <- "Devel_CTspec"

           
## Identify highly variable exons in P56 brain regions keeping cell type constant ------                                  
celltypes <- c("Astro","Oligo","ExciteNeuron","InhibNeuron")

perCT <- list()
perCT_unf <- list()

for(ct in celltypes){
    adult_ct <- adultDF %>% select(Exon,Gene,contains(ct))
    regions <- c("CEREB","HIPP","STRI","THAL","VIS")
    cols <- combn(regions, 2, paste, collapse = "-")

    adult_ct[cols] <- t(apply(adult_ct[3:7], 1, function(x) {
                     out <- combn(x, 2, function(x) x[1] - x[2])
                     c(out, -out)
                    }))
    
     adult_ct <- adult_ct %>% select(-contains("_")) %>% unite("GE",c("Exon","Gene"),sep = "::")
     adult_ct$cellType <- ct
     perCT_unf[[ct]] <- adult_ct                             
     adult_ct <- adult_ct %>% 
            pivot_longer(contains("-"),names_to = "Comparison",values_to = "dPSI") %>% 
            drop_na() %>% filter(abs(dPSI) >= 0.25) 
    perCT[[ct]] <- adult_ct
}
                                  
highVar_BRspec_adult <- do.call("rbind",perCT)
highVar_BRspec_adult$Axis <- "Adult_BRspec"                
                                  
                                  
highVar_BRspec_adult_unf <- do.call("rbind",perCT_unf) %>% remove_rownames()
highVar_BRspec_adult_unf$Axis <- "Adult_BRspec"

                              
## Identify highly variable exons in P56 keeping brain-region constant ------   
regions <- c("VIS","HIPP","STRI","THAL","CEREB")

perReg <- list()
perReg_unf <- list()


for(reg in regions){
    adult_reg <- adultDF %>% select(Exon,Gene,contains(reg))
    celltypes <- c("InhibNeuron","ExciteNeuron","Astro","Oligo")
    cols <- combn(celltypes, 2, paste, collapse = "-")
   
    orig <- colnames(adult_reg)[-c(1:2)]
    

    adult_reg[cols] <- t(apply(adult_reg[3:6], 1, function(x) {
             out <- combn(x, 2, function(x) x[1] - x[2])
             c(out, -out)
            }))   
                                     
    adult_reg <- adult_reg %>% select(-contains("_")) %>% unite("GE",c("Exon","Gene"),sep = "::")
    adult_reg$Region <- reg
    perReg_unf[[reg]] <- adult_reg
    adult_reg <- adult_reg %>% pivot_longer(contains("-"), 
                                            names_to = "Comparison",values_to = "dPSI") %>% 
            drop_na() %>% filter(abs(dPSI) >= 0.25) 
    perReg[[reg]] <- adult_reg
}

highVar_CTspec_adult <- do.call("rbind",perReg)
highVar_CTspec_adult$Axis <- "Adult_CTspec"
dim(highVar_CTspec_adult)
head(highVar_CTspec_adult,2)

highVar_CTspec_adult_unf <- do.call("rbind",perReg_unf) %>% remove_rownames()
highVar_CTspec_adult_unf$Axis <- "Adult_CTspec"

## Putting it all together ------                               
highlyVarExons <- unique(c(highVar_CTspec_adult$GE,highVar_BRspec_adult$GE,highVar_devel$GE,highVar_devel_tp))                               
### Select exon with max abs dPSI ------

BR_all <- highVar_BRspec_adult_unf %>% filter(GE %in% highlyVarExons) %>%
    pivot_longer(contains("-"),names_to = "Comparison",values_to = "dPSI") %>%
    group_by(GE) %>% slice_max(n = 1, order_by = abs(dPSI), with_ties = F)     

CT_all <- highVar_CTspec_adult_unf %>% filter(GE %in% highlyVarExons) %>%
    pivot_longer(contains("-"),names_to = "Comparison",values_to = "dPSI") %>%
    group_by(GE) %>% slice_max(n = 1, order_by = abs(dPSI), with_ties = F)

devel_all <- highVar_devel_unf %>% filter(GE %in% highlyVarExons) %>%
     pivot_longer(contains("-"),names_to = "Comparison",values_to = "dPSI") %>%
    group_by(GE) %>% slice_max(n = 1, order_by = abs(dPSI), with_ties = F)

devel_tp_all <- highVar_devel_tp_unf %>% filter(GE %in% highlyVarExons) %>%
     pivot_longer(contains("-"),names_to = "Comparison",values_to = "dPSI") %>%
    group_by(GE) %>% slice_max(n = 1, order_by = abs(dPSI), with_ties = F)                               
                               
a <- BR_all %>% select(GE, dPSI, Axis)
b <- CT_all %>% select(GE, dPSI, Axis)
c <- devel_all %>% select(GE, dPSI, Axis)
d <- devel_tp_all %>% select(GE, dPSI, Axis)

allAxes <- do.call(rbind,list(a,b,c,d)) %>% group_by(GE,Axis) %>%
    slice_max(n=1,order_by = abs(dPSI),with_ties = F) %>% 
    mutate(dPSI = abs(dPSI)) %>% ungroup() %>%
    pivot_wider(names_from = Axis,values_from = dPSI) %>% 
    replace(is.na(.),-0.05) %>%
    column_to_rownames("GE") 


trueHigh <- allAxes %>% filter_all(any_vars(. >= 0.75) )

complexEvents <- allAxes %>% filter_all(all_vars(. < 0.75) ) %>%
    filter(Adult_BRspec + Devel >= 0.75 | Adult_CTspec + Devel >= 0.75 | Adult_BRspec + Adult_CTspec >= 0.75 |
          Devel + Devel_CTspec >= 0.75 | Adult_BRspec + Devel_CTspec >= 0.75 | Adult_CTspec + Devel_CTspec >= 0.75)

## Determine hVEX ------
HVEx <- allAxes %>% rownames_to_column("GE") %>% separate(GE,into = c("Exon","Gene"), sep = "::")                               
allDataCombined <- left_join(adultDF2,develDF2, by = c("Gene","Exon"))

allDataCombined <- allDataCombined %>% filter(Exon %in% HVEx$Exon)
                               
                               
cols <- combn(colnames(allDataCombined)[3:46], 2, paste, collapse = "-")
allDataCombined[cols] <- t(apply(allDataCombined[3:46], 1, function(x) {
     out <- combn(x, 2, function(x) x[1] - x[2])
     c(out, -out)
}))
allDataCombined <- allDataCombined %>% select(c(Gene,Exon,contains("-")))
write.table(allDataCombined,"../data/HVEx_pairwise_matrix",row.names = FALSE, col.names = T, quote = F, sep = "\t")    
                               
                               
## Preprocessing for circos plot -------
ctspec_hv_df_forCircos <- highVar_CTspec_adult %>% 
    filter(GE %in% rownames(complexEvents) | GE %in% rownames(trueHigh) ) %>% 
    select(Region,Comparison,Axis) %>% 
    group_by_all() %>% add_count() %>% distinct() %>%
    separate(Comparison, into = c("F","T"), sep = "-") %>% 
    mutate(Age = "P56") %>% 
    unite("From", c("Region","Age","F"),sep = ":",remove = F) %>%
    unite("To", c("Region","Age","T"),sep = ":",remove = T) %>% select(-F)
                               
brspec_hv_df_forCircos <- highVar_BRspec_adult %>% 
    filter(GE %in% rownames(complexEvents) | GE %in% rownames(trueHigh) ) %>% 
    select(cellType,Comparison,Axis) %>% 
    group_by_all() %>% add_count() %>% distinct() %>% 
    separate(Comparison, into = c("F","T"), sep = "-") %>%
    mutate(Age = "P56") %>%
     unite("From", c("F","Age","cellType"),sep = ":",remove = F) %>%
     unite("To", c("T","Age","cellType"),sep = ":",remove = T) %>% select(-F)

hv_df_forCircos_adult <- do.call('rbind',list(brspec_hv_df_forCircos,ctspec_hv_df_forCircos))     
                               
write.table(hv_df_forCircos_adult,file = '../data/threeAxes_eVex_df_forCircos_adult', sep = "\t",
           row.names = F, col.names = T, quote = F)
                                                                                      

## Write out EVEx -------
A = data.frame(GE = c(rownames(trueHigh)[ro$`1`],
                    rownames(trueHigh)[ro$`2`],
                    rownames(trueHigh)[ro$`3`],
                    rownames(trueHigh)[ro$`4`],
                    rownames(trueHigh)[ro$`5`]), 
               Cluster = c(rep("devTime",length(ro$`1`)),
                          rep("develAll",length(ro$`2`)),
                          rep("adultCT",length(ro$`3`)),
                          rep("adultAll",length(ro$`4`)),
                          rep("allCT",length(ro$`5`)))) %>%
    separate("GE",into = c("Exon","Gene"),sep = "::") 

write.table(A, '../data/EVEx_4axes',sep = "\t",quote = FALSE,row.names = FALSE, col.names = TRUE)
                               
                               
