library(BgeeDB)
library(XML)

################################################################################
formatData <- function (myBgeeObject, data, stats = NULL, callType = "all") 
{
  if (length(myBgeeObject$quantitativeData) == 0 | length(myBgeeObject$dataType) == 
      0) {
    stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object is valid.")
  }
  if (myBgeeObject$quantitativeData != TRUE) {
    stop("ERROR: formatting quantitative data is not possible for the species and data type of the input Bgee class object.")
  }
  if (length(stats) == 0) {
    if (myBgeeObject$dataType == "affymetrix") {
      cat("\nWARNING: stats parameter set to \"intensities\" for the next steps.\n")
      stats <- "intensities"
    }
    else if (myBgeeObject$dataType == "rna_seq") {
      stop("Please specify the stats parameters. Should be set to \"rpkm\", \"counts\" or \"tpm\"")
    }
  }
  else if (myBgeeObject$dataType == "affymetrix" & stats != 
           "intensities") {
    cat("\nWARNING: For Affymetrix microarray data, stats parameter can only be set to \"intensities\". This will be used for the next steps.\n")
    stats <- "intensities"
  }
  else if (myBgeeObject$dataType == "rna_seq" & compareVersion(gsub("_", 
                                                                    ".", myBgeeObject$release), "13.2") <= 0 & !(stats %in% 
                                                                                                                 c("rpkm", "counts"))) {
    stop("Choose whether data formatting should create a matrix of RPKMs or read counts, with stats option set as \"rpkm\" or \"counts\"")
  }
  else if (myBgeeObject$dataType == "rna_seq" & compareVersion(gsub("_", 
                                                                    ".", myBgeeObject$release), "13.2") > 0 & !(stats %in% 
                                                                                                                c("fpkm", "counts", "tpm", "rank"))) {
    stop("Choose whether data formatting should create a matrix of FPKMs, TPMs or read counts, with stats option set as \"fpkm\", \"tpm\" or \"counts\"")
  }
  if (!(callType %in% c("present", "present high quality", 
                        "all"))) {
    stop("Choose between displaying intensities for present genes, present high quality genes or all genes, e.g., 'present', 'present high quality', 'all' ")
  }
  if (length(data) == 1) 
    data[[1]]
  else data
  cat("\nExtracting expression data matrix...\n")
  if (stats == "rpkm") {
    columns <- c("Library.ID", "Gene.ID", "RPKM")
    expr <- .extract.data(data, columns, callType)
  }
  else if (stats == "fpkm") {
    columns <- c("Library.ID", "Gene.ID", "FPKM")
    expr <- .extract.data(data, columns, callType)
  }
  else if (stats == "tpm") {
    columns <- c("Library.ID", "Gene.ID", "TPM")
    expr <- .extract.data(data, columns, callType)
  }
  else if (stats == "counts") {
    columns <- c("Library.ID", "Gene.ID", "Read.count")
    expr <- .extract.data(data, columns, callType)
  }
  else if (stats == "rank") {
    columns <- c("Library.ID", "Gene.ID", "Rank")
    expr <- .extract.data(data, columns, callType)
  }
  else {
    columns <- c("Chip.ID", "Probeset.ID", "Log.of.normalized.signal.intensity", 
                 "Gene.ID")
    expr <- .extract.data(data, columns, callType)
  }
  if (is.data.frame(expr$assayData)) {
    expr$assayData <- expr$assayData[, order(names(expr$assayData))]
    expr$pheno <- expr$pheno[order(sampleNames(expr$pheno))]
    eset <- new("ExpressionSet", exprs = as.matrix(expr$assayData), 
                phenoData = expr$pheno, featureData = expr$features)
  }
  else if (is.list(expr$assayData)) {
    eset <- mapply(function(x, y, z) {
      x <- x[, order(names(x))]
      y <- y[order(sampleNames(y))]
      new("ExpressionSet", exprs = as.matrix(x), phenoData = y, 
          featureData = z)
    }, expr$assayData, expr$pheno, expr$features)
  }
  return(eset)
}

.extract.data <- function (data, columns, callType) 
{
  if (class(data) == "list") {
    calls <- lapply(data, function(x) .calling(x, callType, 
                                               columns[3]))
    expr <- lapply(calls, function(x) {
      xt <- x[, columns[1:3]]
      xtt <- xt %>% spread_(columns[1], columns[3])
      rownames(xtt) <- xtt[, columns[2]]
      xtt[, -1, drop = FALSE]
    })
    cat("\nExtracting features information...\n")
    features <- mapply(.extract.data.feature, calls, expr, 
                       rep(list(columns), times = length(calls)))
    cat("\nExtracting samples information...\n")
    phenos <- mapply(.extract.data.pheno, calls, rep(list(columns[1]), 
                                                     times = length(calls)))
  }
  else {
    calls <- .calling(data, callType, columns[3])
    xt <- calls[, columns[1:3]]
    xtt <- xt %>% spread_(columns[1], columns[3])
    rownames(xtt) <- xtt[, columns[2]]
    expr <- xtt[, -1, drop = FALSE]
    cat("\nExtracting features information...\n")
    features <- .extract.data.feature(calls, expr, columns)
    cat("\nExtracting samples information...\n")
    phenos <- .extract.data.pheno(calls, columns[1])
  }
  return(list(assayData = expr, pheno = phenos, features = features, 
              calls = calls))
}

.calling <- function (x, callType, column) 
{
  if (callType == "present") {
    cat("  Keeping only present genes.\n")
    x[(x$Detection.flag == "absent"), column] <- NA
  }
  else if (callType == "present high quality") {
    cat("  Keeping only present high quality genes.\n")
    x[which(x$Detection.flag == "absent" | x$Detection.quality == 
              "poor quality"), column] <- NA
  }
  return(x)
}

.extract.data.feature <- function (calls, expr, columns) 
{
  if (length(columns) == 3) {
    fdata <- calls[match(rownames(expr), calls[, columns[2]]), 
                   columns[2], drop = FALSE]
  }
  else if (length(columns) == 4) {
    fdata <- calls[match(rownames(expr), calls[, columns[2]]), 
                   columns[c(2, 4)], drop = FALSE]
  }
  rownames(fdata) <- fdata[, columns[2]]
  fdata <- as(fdata, "AnnotatedDataFrame")
  return(fdata)
}

.extract.data.pheno <- function (calls, column) 
{
  phdata <- calls[, c(column, "Anatomical.entity.ID", "Anatomical.entity.name", 
                      "Stage.ID", "Stage.name")]
  phdata <- phdata[!duplicated(phdata[, column]), ]
  rownames(phdata) <- phdata[, column]
  phdata <- as.data.frame(phdata)
  metadata <- data.frame(labelDescription = colnames(phdata), 
                         row.names = colnames(phdata))
  phenodata <- new("AnnotatedDataFrame", data = phdata, varMetadata = metadata)
  return(phenodata)
}
################################################################################

path <- "/Users/alaverre/Documents/Expressifying/"
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
