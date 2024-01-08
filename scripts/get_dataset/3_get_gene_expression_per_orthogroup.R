#!/usr/bin/env Rscript
library(progress)

################################################################################
# Input arguments
args = commandArgs(trailingOnly=TRUE)
path <-  if (length(args)>0) getwd() else dirname(rstudioapi::getSourceEditorContext()$path)
path <- paste0(path, "/../../")

min_species = if (length(args)>0) args[1] else 10  # minimum number of species to keep an orthogroup, can be "no-missing" for complete case (default=10)
gene_list = if (length(args)>1) args[2] else FALSE # write gene list for GO Enrichment (default=FALSE)

################################################################################
all_orthogroups = read.csv(paste0(path, "data/gene_orthologies/one2one_orthogroups_OMA22_mammals.csv"), row.names = 1)
conditions <- list.dirs(paste0(path, "/data/gene_expression/log2TPM/"), full.names = FALSE, recursive = FALSE)

dir.create(paste0(path, "results/gene_expression_orthogroups/log2TPM/"), recursive=T, showWarnings=F)
dir.create(paste0(path, "results/gene_expression_orthogroups/rank_score/"), recursive=T, showWarnings=F)

for(condition in conditions){
  print(condition)
  output.file.expression = paste0(path, "results/gene_expression_orthogroups/log2TPM/", condition, ".csv")
  output.file.rank = paste0(path, "results/gene_expression_orthogroups/rank_score/", condition, ".csv")
  
  # Check if not already done
  if (!file.exists(output.file.expression)){
    ############################################################################
    print("Reading data for each species...")
    species <- list.files(paste0(path, "data/gene_expression/log2TPM/", condition, "/"), full.names = FALSE)
    species <- sub("\\.csv$", "", species)
    
    gene.expression <- list()
    gene.score <- list()
    for (sp in species){
      gene.expression[[sp]] = read.table(paste0(path, "data/gene_expression/log2TPM/", condition, "/", sp, ".csv"))
      gene.score[[sp]] = read.table(paste0(path, "data/gene_expression/rank_score/", condition, "/", sp, ".csv"))
    }

    # Remove the naked-mole rat because OMA annotations (male) are not the same as Bgee (female) 
    gene.expression[["Heterocephalus_glaber"]] <- NULL
    gene.score[["Heterocephalus_glaber"]] <- NULL
    
    orthogroups = all_orthogroups[,species]
    
    # Select orthogroups with genes for at least 10 species or all species
    if (min_species == "no_missing"){
      orthogroups_all_sp <- orthogroups[complete.cases(orthogroups),]
    }else{
      orthogroups_all_sp = orthogroups[rowSums(!is.na(orthogroups)) >= min_species,]
    }
    
    if (gene_list){
      # Background genes
      orthogroups_human <- orthogroups_all_sp$Homo_sapiens
      orthogroups_human <- orthogroups_human[!is.na(orthogroups_human)]
      write.table(orthogroups_human, file=paste0(path, "results/gene_list/human_genes_", condition, "_1-1_orthogroups.txt"), row.names = F, quote=F, col.names=F)
    }

    ############################################################################
    print("Combining gene expression in each orthogroup...")
    
    # Create empty data.frame for gene expression in each orthogroup
    samples <- unlist(sapply(gene.expression, function(x) colnames(x)))
    orthogroup_IDs <- row.names(orthogroups_all_sp)
    expression_ortho <- data.frame(matrix(ncol=length(orthogroup_IDs)+1, nrow=length(samples)))
    
    colnames(expression_ortho) <- c("species", orthogroup_IDs)
    row.names(expression_ortho) <- samples
    score_ortho <- expression_ortho
    
    ### Complete the data.frame
    # Create a progress bar
    pb <- progress_bar$new(total = length(orthogroup_IDs), format = "[:bar] :percent :elapsed")
    
    for (orthogroup in orthogroup_IDs){
      for (sp in colnames(orthogroups_all_sp)){
        gene = as.character(orthogroups_all_sp[orthogroup,sp])
        samples = colnames(gene.expression[[sp]])
        
        # Get the expression of this gene if present in data
        if (gene %in% rownames(gene.expression[[sp]])){
          expression = t(gene.expression[[sp]][gene,])
          score = t(gene.score[[sp]][gene,])
        }else{
          expression = NA
          score = NA
        }
      
        expression_ortho[samples, orthogroup] = expression
        expression_ortho[samples, "species"] = sp
        score_ortho[samples, orthogroup] = score
        score_ortho[samples, "species"] = sp
      }
      pb$tick()
    }
    
    ################################################################################
    print("Filtering and writting output...")
    
    ## Filters to get gene expression for each species in each orthogroup
    # Count the number of non-NA values per species
    species_without_na <- aggregate(!is.na(expression_ortho[, -1]), by = list(expression_ortho$species), FUN = sum)
    
    # Keep orthogroup with data for at least min_species
    nb_sp_ortho <- apply(species_without_na[, -1], 2, function(x) sum(x >= 2, na.rm = T))
    ortho_all_sp <- names(nb_sp_ortho[which(nb_sp_ortho >= min_species)])
    
    express_ortho_all_sp <- expression_ortho[,c("species", ortho_all_sp)]
    score_ortho_all_sp <- score_ortho[,c("species", ortho_all_sp)] 
    
    # Print summary and save final data.frame per tissue
    nb_ortho = ncol(express_ortho_all_sp)
    print(paste(length(species), "species;", nrow(express_ortho_all_sp), "samples;", nb_ortho-1, "genes."))
    
    if (gene_list){
      # Background expressed genes
      orthogroups_human <- orthogroups_all_sp[colnames(express_ortho_all_sp)[2:nb_ortho],"Homo_sapiens"]
      orthogroups_human <- orthogroups_human[!is.na(orthogroups_human)]
      write.table(orthogroups_human, file=paste0(path, "results/gene_list/mammals_", tissue, "_expressed_filtered.txt"), row.names = F, quote=F, col.names=F)
    }
    
    # add sample as column
    express_ortho_all_sp <- cbind(sample=row.names(express_ortho_all_sp), express_ortho_all_sp)
    score_ortho_all_sp <- cbind(sample=row.names(score_ortho_all_sp), score_ortho_all_sp)
    
    write.csv(express_ortho_all_sp, file=output.file.expression, row.names = F)
    write.csv(score_ortho_all_sp, file=output.file.rank, row.names = F)
  
  }else{print("Already done!")}
}

################################################################################