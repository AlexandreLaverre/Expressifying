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