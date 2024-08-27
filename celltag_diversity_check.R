library("CellTagR")
# extract the dataset of interest, will need to set to appropriate working directory
# the FASTQ file name changes from file to file as needed
extract.dir <- "."
full.fpath <- paste0(extract.dir, "/S211.V1-BFP/", "S211-CellTag-V1-diversity_S1_L001_R1_001.fastq")
# set up the CellTag Object
test.obj <- CellTagObject(object.name = "v1.whitelist", fastq.bam.directory = full.fpath)
# extract the CellTags from the FASTQ
test.obj <- CellTagExtraction(celltag.obj = test.obj, celltag.version = "v1")
# count and sort the CellTags in descending order of occurrence
test.obj <- AddCellTagFreqSort(test.obj)
# visually check the statistics the stats
test.obj@celltag.freq.stats
# generate the whitelist
test.obj <- CellTagWhitelistFiltering(celltag.obj = test.obj, percentile = 0.95, output.dir = NULL)