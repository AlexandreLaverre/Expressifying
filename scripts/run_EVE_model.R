#!/usr/bin/env Rscript
library(evemodel)
library(ape)

args = commandArgs(trailingOnly=TRUE)
tissue = args[1]

path <- "/Users/alaverre/Documents/Expressifying/"

################################################################################
##### DATA ##### 
# Gene orthology
ortho <- read.csv(paste0(path, "data/gene_orthologies/one2one_orthogroups.csv"))
human.genes <- ortho$Homo_sapiens
names(human.genes) <- ortho$Orthogroups

# Gene expression
gene.express <- read.csv(paste0(path, "results/gene_expression/mammals_", tissue, "_gene_expression_log2TPM_orthogroups.csv"))

gene.express.test <- t(gene.express[,3:ncol(gene.express)])
colnames(gene.express.test) <- paste0(gene.express$species, "_", gene.express$sample)
exprMat <- as.matrix(gene.express.test)
colSpecies <- sub("_SRX.*$","",colnames(exprMat))
colSpecies <- sub("_ERX.*$","",colSpecies)
colSpecies <- sub("_DRX.*$","",colSpecies)
message("Original data: ", nrow(exprMat), " genes, ", ncol(exprMat), " samples and ", length(unique(colSpecies)), " species.")

par(mai=c(2.2, 1, 0.5, 1))
plot(table(colSpecies), xlab="", ylab="Nb samples", las=2, main="")

# Filters
# Reduce missing values
exprMat <- exprMat[,colSums(is.na(exprMat)) < 2000] # remove samples with more than 2000 NA genes 

# Subsample to reduce the over-represented species
max_samples <- 25
for (sp in colSpecies){
  sp_samples <- grep(sp, colnames(exprMat), value = TRUE)
  if (length(sp_samples) > max_samples){
    
    # Remove samples with the higher number of NA genes
    nb_to_remove <- length(sp_samples)-max_samples
    NA_count <- apply(exprMat[,sp_samples], 2, function(x) sum(is.na(x)))
    samples_to_remove <- names(NA_count[order(NA_count, decreasing=T)[1:nb_to_remove]])
    samples_to_keep <- !colnames(exprMat) %in% samples_to_remove
    
    exprMat <- exprMat[,samples_to_keep]
  }
}

exprMat <- exprMat[rowSums(is.na(exprMat)) < 40, ] # remove gene with more than 40 NA samples 

colSpecies <- sub("_SRX.*$","",colnames(exprMat))
colSpecies <- sub("_ERX.*$","",colSpecies)
colSpecies <- sub("_DRX.*$","",colSpecies)
message("Filtered data: ", nrow(exprMat), " genes, ", ncol(exprMat), " samples and ", length(unique(colSpecies)), " species.")

plot(table(colSpecies), xlab="", ylab="Nb samples", las=2, main="Distribution of liver samples filtered")
table(colSpecies)

# Reformat Gene Expression to our input for comparison
filtered.gene.express <- as.data.frame(t(exprMat))

orthogroups <- colnames(filtered.gene.express)
species <- sub("_SRX.*$","",row.names(filtered.gene.express))
species <- sub("_ERX.*$","",species)
filtered.gene.express$species <- species
filtered.gene.express$sample <- gsub("^.*_", "", row.names(filtered.gene.express))
filtered.gene.express <- filtered.gene.express[,c("sample", "species", orthogroups)]

write.csv(filtered.gene.express, file=paste0(path, "results/gene_expression/mammals_", tissue, "_gene_expression_log2TPM_orthogroups_EVE_filtered.csv"), row.names = F)

################################################################################
# Species Tree from Zoonomia
speciesTree <- read.tree("/Users/alaverre/Documents/Detecting_positive_selection/data/species_trees/241-mammals.nk")

# keep only species with gene expression data
common.species <- intersect(unique(colSpecies), speciesTree$tip.label)
speciesTree <- keep.tip(speciesTree, common.species)

plot(speciesTree)

################################################################################
# Running beta shared test

res <- betaSharedTest(tree = speciesTree, gene.data = exprMat, colSpecies = colSpecies)
saveRDS(res, file=paste0(path, "results/EVE_model/", tissue, "_log2TPM_betaSharedTest.RDS"))

################################################################################
# Format output
EVE_results <- as.data.frame(res[["indivBetaRes"]][["par"]])
EVE_results$SharedBeta <- res$sharedBeta
EVE_results$LL <-  res[["indivBetaRes"]][["ll"]]
EVE_results$LRT <- res$LRT
EVE_results$iterations <-  res[["indivBetaRes"]][["iterations"]]
EVE_results$pval <- pchisq(res$LRT,df = 1,lower.tail = F)
EVE_results$p.adjusted <- p.adjust(EVE_results$pval)
rownames(EVE_results) <- rownames(exprMat)
EVE_results$Human_ID <- human.genes[rownames(EVE_results)]

write.table(EVE_results, file=paste0(path, "results/EVE_model/", tissue, "_log2TPM_betaSharedTest.txt"))

################################################################################