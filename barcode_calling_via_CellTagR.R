# CellTagR installed in local home directory
# This must be run in R version 4.2
library("CellTagR")
setwd("/fh/fast/greenberg_p/user/jwlee/S211/S211_rd_4/")

# Modifying CellTagPatternCalling to make a more specific pattern for S211 V1
# For the custom sequencing primer, there is no V1 prefix. Assume that all of
# the UTR-aligned sequences contain it though, since the sequencing primer binds
# to the V1 sequence.
CellTagPatternCalling <- function (celltag.version) {
  celltag.df <- data.frame(version = c("v1", "v2", 
                                       "v3"), nt.before.tag = c("ACCGGT", "GTGATG", 
                                                                "TGTACG"), stringsAsFactors = F)
  rownames(celltag.df) <- celltag.df$version
  short.nt.before.tag <- celltag.df[celltag.version, "nt.before.tag"]
  short.nt.after.tag <- "GAATTC"
  pattern <- paste0(short.nt.before.tag, "[ATCG]{8}", 
                    short.nt.after.tag)
  return(c(pattern, short.nt.before.tag, short.nt.after.tag))
}
# Copying CellTagExtraction to here so that it can use CellTagPatternCalling
# as defined in this working environment rather than in the CellTagR namespace
CellTagExtraction <- function (celltag.obj, celltag.version, technique = "10x") {
  celltag.obj@curr.version <- celltag.version
  if (file_test("-f", celltag.obj@fastq.bam.dir)) {
    fastq.bam.input <- celltag.obj@fastq.bam.dir
  }
  else {
    fastq.bam.input <- list.files(celltag.obj@fastq.bam.dir, 
                                  full.names = T)
  }
  file.extension.unique <- unique(file_ext(fastq.bam.input))
  if (length(celltag.obj@celltag.version) > 0) {
    if (celltag.obj@curr.version %in% celltag.obj@celltag.version) {
      print("This CellTag has already been processed!")
    }
    else {
      celltag.obj@celltag.version <- c(celltag.obj@celltag.version, 
                                       celltag.obj@curr.version)
    }
  }
  else {
    celltag.obj@celltag.version <- celltag.obj@curr.version
  }
  p.calling <- CellTagPatternCalling(celltag.version)
  if (endsWith(file.extension.unique, "fastq") || endsWith(file.extension.unique, 
                                                           "fq")) {
    if (length(fastq.bam.input) > 1) {
      stop("Please process the whitelist files one at a time!")
    }
    rslt <- fastq.process(fastq.file = fastq.bam.input, 
                          pattern = p.calling[1], p.calling[2], p.calling[3])
    celltag.obj@fastq.full.celltag[[celltag.version]] <- rslt[[1]]
    celltag.obj@fastq.only.celltag[[celltag.version]] <- rslt[[2]]
  }
  if (endsWith(file.extension.unique, "bam")) {
    rslt <- NULL
    for (i in 1:length(fastq.bam.input)) {
      curr.rslt <- bam.process(bam.file = fastq.bam.input[i], 
                               pattern = p.calling[1], p.calling[2], p.calling[3], 
                               technique)
      if (length(fastq.bam.input) > 1) 
        curr.rslt$Cell.BC <- paste0("Sample-", i, "_", 
                                    curr.rslt$Cell.BC)
      if (is.null(rslt)) {
        rslt <- curr.rslt
      }
      else {
        rslt <- rbind(rslt, curr.rslt, fill = TRUE)
      }
    }
    celltag.obj@bam.parse.rslt[[celltag.version]] <- rslt
  }
  return(celltag.obj)
}
# Modyfing SingleCellDataWhitelist to account for potentially faulty whitelist
# No usage of whitelist
SingleCellDataWhitelist <- function (celltag.obj, whitels.cell.tag.file, replace.option = FALSE) {
  CellTags <- as.matrix(GetCellTagCurrentVersionWorkingMatrix(celltag.obj, 
                                                              "binary.mtx"))
  cell.names <- rownames(CellTags)
  CellTags <- t(CellTags)
  celltag.rownames <- row.names(CellTags)
  if (endsWith(whitels.cell.tag.file, ".csv")) {
    separator <- ","
  }
  else {
    if (endsWith(whitels.cell.tag.file, ".txt") | endsWith(whitels.cell.tag.file, 
                                                           ".tsv")) {
      separator <- "\t"
    }
    else {
      separator <- " "
    }
  }
  whitelist <- read.delim(whitels.cell.tag.file, sep = separator, 
                          header = T, stringsAsFactors = F)
  whitelist.names <- whitelist[, 1]
  whitelist <- Reduce(intersect, list(whitelist.names, celltag.rownames))
  celltags.whitelisted <- CellTags #[whitelist, ]
  colnames(celltags.whitelisted) <- cell.names
  new.obj <- SetCellTagCurrentVersionWorkingMatrix(celltag.obj, 
                                                   "whitelisted.count", as(t(as.matrix(celltags.whitelisted)), 
                                                                           "dgCMatrix"))
  return(new.obj)
}
GetCellTagCurrentVersionWorkingMatrix <- function(celltag.obj, slot.to.select) {
  curr.mtx <- slot(celltag.obj, slot.to.select)
  if (nrow(curr.mtx) <= 0) {
    return(curr.mtx)
  } else {
    curr.version <- celltag.obj@curr.version
    curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
    colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
    full.mtx.sub <- curr.mtx.sub[Matrix::rowSums(is.na(curr.mtx.sub)) != ncol(curr.mtx.sub),]
    
    return(full.mtx.sub)
  }
}
SetCellTagCurrentVersionWorkingMatrix <- function(celltag.obj, slot.to.set, final.to.set) {
  cop.final <- final.to.set
  colnames(cop.final) <- paste0(celltag.obj@curr.version, ".", colnames(cop.final))
  curr.version.existing.mtx <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, slot.to.set)
  
  if (sum(dim(slot(celltag.obj, slot.to.set))) <= 0) {
    slot(celltag.obj, slot.to.set) <- cop.final
  } else  {
    curr.existing.mtx <- slot(celltag.obj, slot.to.set)
    if (ncol(curr.version.existing.mtx) > 0) {
      curr.ver.exist.colnames <- paste0(celltag.obj@curr.version, ".", colnames(curr.version.existing.mtx))
      indx <- which(colnames(curr.existing.mtx) %in% curr.ver.exist.colnames)
      curr.existing.mtx <- curr.existing.mtx[, -indx]
    }
    new.rownames <- unique(c(rownames(curr.existing.mtx), rownames(cop.final)))
    
    diff.rnms <- setdiff(new.rownames, rownames(cop.final))
    cop.comp.mtx <- matrix(NA, nrow = length(diff.rnms), ncol = ncol(cop.final))
    rownames(cop.comp.mtx) <- diff.rnms
    colnames(cop.comp.mtx) <- colnames(cop.final)
    
    diff.rnms.2 <- setdiff(new.rownames, rownames(curr.existing.mtx))
    cem.comp.mtx <- matrix(NA, nrow = length(diff.rnms.2), ncol = ncol(curr.existing.mtx))
    rownames(cem.comp.mtx) <- diff.rnms.2
    colnames(cem.comp.mtx) <- colnames(curr.existing.mtx)
    
    to.merge.mtx.cop <- rbind(cop.final, cop.comp.mtx)
    to.merge.mtx.cem <- rbind(curr.existing.mtx, cem.comp.mtx)
    
    if (ncol(to.merge.mtx.cem) <= 0) {
      new.mtx <- to.merge.mtx.cop[,colnames(cop.final)]
    } else {
      new.mtx <- cbind(to.merge.mtx.cop[,colnames(cop.final)], to.merge.mtx.cem[,colnames(cop.final)])
    }
    
    slot(celltag.obj, slot.to.set) <- new.mtx
  }
  
  return(celltag.obj)
} 
# Alternate function for getting linked list and nodes
# This is accounting for using only V1 and V2 without V3.
convertCellTagMatrix2LinkList <- function (celltag.obj) {
  celltag.dt <- celltag.obj@clone.composition
  v1.df <- as.data.frame(celltag.dt$v1)
  v2.df <- as.data.frame(celltag.dt$v2)
  rownames(v1.df) <- v1.df$cell.barcode
  rownames(v2.df) <- v2.df$cell.barcode
  all.cells <- unique(c(celltag.dt$v1$cell.barcode, celltag.dt$v2$cell.barcode))
  celltag_data <- data.frame(row.names = all.cells)
  celltag_data[rownames(v1.df), "CellTagV1"] <- v1.df[rownames(v1.df), 
                                                      "clone.id"]
  celltag_data[rownames(v2.df), "CellTagV2"] <- v2.df[rownames(v2.df), 
                                                      "clone.id"]
  celltag.obj@celltag.aggr.final <- celltag_data
  message("Preprocessing data..")
  Cells_with_tag <- rownames(celltag_data)[!(is.na(celltag_data$CellTagV1) & 
                                               is.na(celltag_data$CellTagV2))] # & is.na(celltag_data$CellTagV3))]
  message(paste0(" Cells that have CellTagV1: ", sum(!is.na(celltag_data$CellTagV1))))
  message(paste0(" Cells that have CellTagV2: ", sum(!is.na(celltag_data$CellTagV2))))
  celltag_data <- celltag_data[Cells_with_tag, ]
  tags <- c("CellTagV1", "CellTagV2")#, "CellTagV3")
  for (i in tags) {
    celltag_data[is.na(celltag_data[, i]), i] <- "e"
  }
  message("Constructing link list..")
  findRoot <- function(cell_id, tag) {
    tagid <- celltag_data[cell_id, tag]
    tmp <- as.data.frame(t(c(paste0(tag, "_", tagid), 
                             cell_id, tag)), stringsAsFactors = F)
    rownames(tmp) <- NULL
    colnames(tmp) <- c("source", "target", "tag")
    return(tmp)
  }
  all_cell_id <- rownames(celltag_data)
  remaining_cell_id <- all_cell_id
  tags <- c("CellTagV2", "CellTagV1") #"CellTagV3", "CellTagV2", "CellTagV1")
  linkList <- data.frame()
  for (tag in tags) {
    remaining_cells <- celltag_data[remaining_cell_id, ]
    subcells <- remaining_cells[remaining_cells[, tag] != 
                                  "e", ]
    tmp <- foreach(i = rownames(subcells), .combine = rbind, 
                   .packages = "foreach") %do% {
                     findRoot(i, tag)
                   }
    linkList <- rbind(linkList, tmp)
    done_id <- rownames(subcells)
  }
  hiddenlink_D3 <- foreach(i = (unique(celltag_data$CellTagV1)[-1]), 
                           .combine = rbind, .packages = "foreach") %do% {
                             sub_cells <- celltag_data[celltag_data$CellTagV1 == i, 
                             ]
                             prev_tag <- sub_cells$CellTagV1
                             prev_tag <- prev_tag[prev_tag != "e"]
                             prev_tag <- names(which.max(table(prev_tag)))
                             if (class(prev_tag) != "NULL") {
                               tmp <- as.data.frame(t(c(paste0("CellTagV1", 
                                                               "_", prev_tag), paste0("CellTagV1", 
                                                                                      "_", i), "CellTagV1")), stringsAsFactors = F)
                               rownames(tmp) <- NULL
                               colnames(tmp) <- c("source", "target", 
                                                  "tag")
                               return(tmp)
                             }
                           }
  rm(remaining_cells, remaining_cell_id, sub_cells, subcells, 
     all_cell_id, done_id, prev_tag, tag, tags)
  modifyCellName <- function(linkList) {
    linkList$target_unmodified <- linkList$target
    node_cell <- grep("-", linkList$target)
    linkList[node_cell, "target"] <- paste0(linkList[node_cell, 
                                                     "target"], "_", stringr::str_split_fixed(linkList[node_cell, 
                                                                                                       "tag"], "g", 2)[, 2])
    return(linkList)
  }
  linkList <- rbind(linkList, hiddenlink_D3)
  linkList <- modifyCellName(linkList)
  message("finished")
  celltag.obj@network.link.list <- linkList
  return(celltag.obj)
}

# Get list of files to process. Note the first entry is this directory
dirs <- c(list.dirs('.', recursive=FALSE))
# Put the celltag processing to loop through the dirs
for (dir in dirs[1:length(dirs)]){
  # Set up CellTag object
  bam.test.obj <- CellTagObject(object.name = "bam.cell.tag.obj",
                                fastq.bam.directory = paste(dir, "/alignment.bam", sep=''))
  # Extract the CellTag information, V1 first
  bam.test.obj <- CellTagExtraction(bam.test.obj, celltag.version = "v1")
  # Check the bam file result
  print(head(bam.test.obj@bam.parse.rslt[["v1"]]))
  # Generate the sparse count matrix
  bam.test.obj <- CellTagMatrixCount(celltag.obj = bam.test.obj, barcodes.file = paste(dir, '/barcodes.tsv.gz', sep=''))  
  # Check the dimension of the raw count matrix
  print(dim(bam.test.obj@raw.count))
  # Save the object
  saveRDS(bam.test.obj, file = paste(dir, "/bam.test.obj.lower_threshold.rds", sep=''))
  # Celltag error correction
  # Generating the collapsing file
  bam.test.obj <- CellTagDataForCollapsing(celltag.obj = bam.test.obj, output.file = 
                                             paste(dir, "/collapsing_lower_threshold.txt", sep=''))
  # Ensure that starcode/1.4-GCC-11.2.0 is loaded in the session first before running the next line
  system(paste('../starcode/starcode -s --print-clusters -t 36 ', dir, '/collapsing_lower_threshold.txt > ',
                dir, '/collapsing_result_lower_threshold.txt', 
               sep=''))
  # Recount and generate collapsed matrix
  bam.test.obj <- CellTagDataPostCollapsing(celltag.obj = bam.test.obj, 
                                            collapsed.rslt.file = paste(dir, "/collapsing_result_lower_threshold.txt", sep=''))
  # Calling binarization
  bam.test.obj <- SingleCellDataBinarization(bam.test.obj, 1)
  # Read the RDS file and get the object
  bam.test.obj <- SingleCellDataWhitelist(bam.test.obj, "../diversity_check/S211-BFP_V1.BFP_v1.whitelist.csv")
  # Lowering threshold to 1 to see if it improves traceability
  bam.test.obj <- MetricBasedFiltering(bam.test.obj, 1, comparison = "greater")
  # Clone calling with Jaccard analysis first
  bam.test.obj <- JaccardAnalysis(bam.test.obj, fast=TRUE)
  # Call clones
  bam.test.obj <- CloneCalling(celltag.obj = bam.test.obj, correlation.cutoff=0.7)
  # Write the raw celltags out per cell too
  write.table(bam.test.obj@bam.parse.rslt, file=paste(dir, "/output/bam.parse.rslt.v1.txt", sep=''), quote=FALSE, 
              row.names=FALSE)
  write.table(bam.test.obj@clone.size.info[["v1"]], file=paste(dir, "/output/clone.size.info.v1.txt", sep=''),
              quote=FALSE, row.names=FALSE)
}
# Only V1 was utilized due to its better coverage and trace-ability further back in time
