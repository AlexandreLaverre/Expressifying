#!/usr/bin/env Rscript
library(BgeeDB)

################################################################################
# Input arguments
args = commandArgs(trailingOnly=TRUE)
path <- ifelse(length(args)>0, getwd(), dirname(rstudioapi::getSourceEditorContext()$path))
path <- paste0(path, "/../../")

tissue = ifelse(length(args)>0, args[1], "liver")   # default = "liver"
stage = ifelse(length(args)>1, args[2], "post-juvenile") # default = "post-juvenile"

#example.organs <- c("kidney", "liver", "heart", "lung", "testis", "cerebellum",
#                  "spinal cord", "colon", "pituitary gland", "uterus", "muscle tissue",
#                  "blood",  "zone of skin", "esophagus", "stomach", "intestine",
#                  "bone element", "spleen", "brain", "ovary") 

#example.stage <- c("post-juvenile", "embryonic", "fully formed")

################################################################################
# BgeeDB functions to properly retrieve gene expression rank
source(paste0(path, "/scripts/get_dataset/modified_bgee_functions.R"))

# Bgee species list manually annotated for "Taxa"
species <- read.csv(paste0(path, "/data/Bgee_libraries/species_list.csv"), header=T, fill=T)
mammals <- species[which(species$Large_Taxa == "Mammal" | species$Large_Taxa == "Marsupial" ), "complete_name"]

# Retrieve UBERON ID for anatomical entity
Anatomical.UBERON <- read.csv(paste0(path, "/data/AnatomicalEntity.to.ID.csv"), h = T, sep="\t")
Anat.Name2ID <- Anatomical.UBERON$Anatomical.entity.ID
names(Anat.Name2ID) <-  Anatomical.UBERON$Anatomical.entity.name
tissue.ID = Anat.Name2ID[[tissue]]

# Retrieve UBERON ID for stage name
UBERON.stage <- c("UBERON:0000113", "UBERON:0000068", "UBERON:0000066")
names(UBERON.stage) <- c("post-juvenile", "embryonic", "fully formed") 
stage.ID = UBERON.stage[[stage]]

################################################################################
# Retrieve gene expression (TPM and rank score) in selected tissue 
dir.create(paste0(path, "/data/gene_expression/log2TPM/", stage, "_", tissue, "/"), recursive=T, showWarnings=F)
dir.create(paste0(path, "/data/gene_expression/rank_score/", stage, "_", tissue, "/"), recursive=T, showWarnings=F)

for (sp in mammals){
  print(paste0("#####################", sp, "#####################"))
  output.expression <- paste0(path, "/data/gene_expression/log2TPM/", stage, "_", tissue, "/", sp, "gene_expression.csv")
  output.score <- paste0(path, "/data/gene_expression/rank_score/", stage, "_", tissue, "/",  sp, "_gene_expression.csv")
  
  if (!file.exists(output.expression)){
    bgee <- Bgee$new(species=sp, dataType = "rna_seq", pathToData=paste0(path, "/data/Bgee_libraries/"))
    data <- getData(bgee, anatEntityId = tissue.ID, withDescendantAnatEntities = TRUE, 
                    stageId = stage.ID, withDescendantStages = TRUE)
    
    # At least 2 libraries per condition for variance
    if (length(unique(data$Library.ID)) >= 2){
      
      data$Rank <- as.numeric(data$Rank)
      simplified.data <- formatData(bgee, data, callType = "present", stats = "tpm")
      all_gene_expression <- log2(0.01+simplified.data@assayData[["exprs"]])
      
      simplified.data <- formatData(bgee, data, callType = "present", stats = "rank")
      all_gene_rank <- simplified.data@assayData[["exprs"]]
      
      # New score based on rank normalized per library
      all_gene_rank_normalised <- apply(all_gene_rank, 2, function(x) 100-(x*100)/max(x, na.rm=T))
      
      # Remove genes with less than 2 measures
      gene.expression <- all_gene_expression[rowSums(!is.na(all_gene_expression)) >= 2,]
      gene.score <- all_gene_rank_normalised[rowSums(!is.na(all_gene_rank_normalised)) >= 2,]
      
      # Write output
      write.table(gene.expression, file=output.expression, quote=F, sep="\t")
      write.table(gene.score, file=output.score, quote=F, sep="\t")
    }
  }else{print("Already done!")}
}

################################################################################
