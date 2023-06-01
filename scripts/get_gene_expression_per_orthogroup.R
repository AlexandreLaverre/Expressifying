path = "/Users/alaverre/Documents/Expressifying/"

tissues <- c("liver", "cerebellum", "kidney", "testis")

for(tissue in tissues){
  print(tissue)
  gene.expression = readRDS(paste0(path, "data/gene_expression/mammals_", tissue, "_gene_expression.Rds"))
  
  # Correct dog species name 
  dog_index <- which(names(gene.expression) == "Canis_lupus familiaris")
  names(gene.expression)[dog_index] <- "Canis_lupus_familiaris"
  
  # Remove the naked-mole rat because OMA annotations (male) are not the same as Bgee (female) 
  gene.expression[["Heterocephalus_glaber"]] <- NULL
  
  species = names(gene.expression)
  
  orthogroups = read.csv(paste0(path, "data/gene_orthologies/one2one_orthogroups.csv"), row.names = 1)
  orthogroups = orthogroups[,species]
  
  # Select orthogroups with genes for all species
  #orthogroups_15_sp = orthogroups[rowSums(!is.na(orthogroups)) > 15,]
  orthogroups_all_sp =  orthogroups[complete.cases(orthogroups),]
  
  # Create empty data.frame for gene expression in each orthogroup
  samples <- unlist(sapply(gene.expression, function(x) colnames(x)))
  orthogroup_IDs <- row.names(orthogroups_all_sp)
  expression_ortho <- data.frame(matrix(ncol=length(orthogroup_IDs)+1, nrow=length(samples)))
  colnames(expression_ortho) <- c("species", orthogroup_IDs)
  row.names(expression_ortho) <- samples
  
  # Complete the data.frame
  for (orthogroup in orthogroup_IDs){
    for (sp in colnames(orthogroups_all_sp)){
      
      gene = as.character(orthogroups_all_sp[orthogroup,sp])
      samples = colnames(gene.expression[[sp]])
      
      # Get the expression of this gene if present in data
      if (gene %in% rownames(gene.expression[[sp]])){
        expression = gene.expression[[sp]][gene,]
      }else{expression = NA}
      
      expression_ortho[samples, orthogroup] = expression
      expression_ortho[samples, "species"] = sp
    }
  }
  
  ## Filters to get gene expression for each species in each orthogroup
  # Count the number of non-NA values per species
  species_without_na <- aggregate(!is.na(expression_ortho[, -1]), by = list(expression_ortho$species), FUN = sum)
  
  # Get orthogroups with at least 2 samples for each species
  nb_sp_ortho <- apply(species_without_na[, -1], 2, function(x) sum(x >= 2, na.rm = T))
  ortho_all_sp <- names(nb_sp_ortho[which(nb_sp_ortho == length(species))])
  
  express_ortho_all_sp <- expression_ortho[,c("species", ortho_all_sp)] 
  
  # Print summary and save final data.frame per tissue
  print(paste(length(species), "species;", 
              nrow(express_ortho_all_sp), "samples;", 
              ncol(express_ortho_all_sp)-1, "genes."))

  write.csv(express_ortho_all_sp, file=paste0(path, "results/mammals_", tissue, "_gene_expression_orthogroups.csv"))
}

# liver = 24 species, 585 samples, 481 genes.
# cerebellum = 21 species, 112 samples, 577 genes.
# kidney = 21 species, 153 samples, 600 genes. 
# testis = 19 species, 160 samples, 750 genes.
