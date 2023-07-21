library(BgeeDB)
library(XML)

################################################################################

path <- "/Users/alaverre/Documents/Expressifying/"
source(paste0(path, "/scripts/modified_bgee_functions.R"))
species <- read.csv(paste0(path, "/data/Bgee_libraries/species_list.csv"), header=T, fill=T)
mammals <- species[which(species$Large_Taxa == "Mammal" | species$Large_Taxa == "Marsupial" ), "complete_name"]
#mammals[mammals != "Monodelphis_domestica"]

# Get numbers of experiments per IDs
IDs <- list()
for (sp in mammals){
  bgee <- Bgee$new(species=sp, dataType = "rna_seq", pathToData=paste0(path, "/data/Bgee_libraries/"))
  
  annotation_bgee <- getAnnotation(bgee)$sample.annotation
  annotation_bgee$ID <- paste(annotation_bgee$Anatomical.entity.name, annotation_bgee$Anatomical.entity.ID)
  
  IDs[[sp]] <- as.data.frame(table(annotation_bgee$ID))
  colnames(IDs[[sp]]) <- c("ID", sp)
}

# Merging tables
IDs_count <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), IDs)
rownames(IDs_count) <- IDs_count$ID
IDs_count <- IDs_count[, -1]

# Counting the number of species with at least 2 experiments
IDs_count$Nb_sp_var <- apply(IDs_count, 1, function(x) sum(x > 1, na.rm = T))
IDs_count <- IDs_count[order(IDs_count$Nb_sp_var, decreasing=T),]
head(IDs_count)

# Get species list for tissues with the highest variance
liver_variance <- head(colnames(IDs_count[which(IDs_count[1,] > 1)]), -1) # UBERON:0002107 
adult_kidney_variance <- head(colnames(IDs_count[which(IDs_count[2,] > 1)]), -1) # UBERON:0000082  
cerebellum_variance <- head(colnames(IDs_count[which(IDs_count[3,] > 1)]), -1) #UBERON:0002037
testis_variance  <- head(colnames(IDs_count[which(IDs_count[4,] > 1)]), -1) # UBERON:0000473

species_list <- list(liver_variance, adult_kidney_variance, cerebellum_variance, testis_variance)
names(species_list) <- c("liver", "kidney", "cerebellum", "testis")
UBERON <- c("UBERON:0002107", "UBERON:0000082", "UBERON:0002037", "UBERON:0000473")
names(UBERON) <- names(species_list) 

################################################################################
# Retrieve gene expression data for each tissue
for (tissu in names(species_list)){
  print(tissu)
  output.expression <- paste0(path, "/data/gene_expression/mammals_", tissu, "_gene_expression_logTPM.Rds")
  output.score <- paste0(path, "/data/gene_expression/mammals_", tissu, "_gene_expression_score.Rds")
  
  if (file.exists(output.expression)){
    print("Already done!")
  }else{
  
    gene.expression <- list()
    gene.score <- list()
    for (sp in species_list[[tissu]]){
      print(sp)
      bgee <- Bgee$new(species=sp, dataType = "rna_seq", pathToData=paste0(path, "/data/Bgee_libraries/"))
      data <- getData(bgee, anatEntityId = UBERON[[tissu]])
      data$Rank <- as.numeric(data$Rank)
      simplified.data <- formatData(bgee, data, callType = "present", stats = "tpm")
      all_gene_expression <- log(simplified.data@assayData[["exprs"]])
      
      simplified.data <- formatData(bgee, data, callType = "present", stats = "rank")
      all_gene_rank <- simplified.data@assayData[["exprs"]]
      all_gene_rank_normalised <- apply(all_gene_rank, 2, function(x) 100-(x*100)/max(x, na.rm=T))
      
      # Remove genes with less than 2 measures
      gene.expression[[sp]] <- all_gene_expression[rowSums(!is.na(all_gene_expression)) >= 2,]
      gene.score[[sp]] <- all_gene_rank_normalised[rowSums(!is.na(all_gene_rank_normalised)) >= 2,]
      
    }
    saveRDS(gene.expression, file=output.expression)
    saveRDS(gene.score, file=output.score)
  }
}

sp="Bos_taurus"
plot(gene.expression[[sp]][,1]~gene.score[[sp]][,1], cex=0.1,
     xlab="Expression score", ylab="log(TPM)", main=paste(sp, tissu, "in", colnames(gene.score[[sp]])[1]))

par(mfrow=c(2,1))
par(mai=c(0.4,0.8,0.3,0.3), mgp=c(2.2,0.8,0))
boxplot(t(gene.expression[[sp]][1:20,]), outline=F, 
        ylab="log(TPM)", las=1, main="20 genes in liver Bos_taurus")
boxplot(t(gene.score[[sp]][1:20,]), outline=F,
        ylab="Score", las=1, main="")

Nsample = length(colnames(gene.expression[[sp]]))
median_expression <- apply(gene.expression[[sp]][,1:Nsample], 1, function(x) median(log(x), na.rm=T))
median_score <- apply(gene.score[[sp]][,1:Nsample], 1, function(x) median(x, na.rm=T))

plot(median_expression~median_score, cex=0.1, 
     xlab="Median score", ylab="median log(TPM)", main=paste("Median across all samples in", tissu, sp))

# SD
sd_expression <- apply(gene.expression[[sp]][,1:Nsample], 1, function(x) sd(log(x), na.rm=T))
sd_score <- apply(gene.score[[sp]][,1:Nsample], 1, function(x) sd(x, na.rm=T))

plot(sd_expression~sd_score, cex=0.1, 
     xlab="Standard Deviation Score", ylab="Standard Deviation log(TPM)", main=paste("SD across all samples in", tissu, sp))



################################################################################
