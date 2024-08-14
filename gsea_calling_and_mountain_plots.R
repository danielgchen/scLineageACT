### R
# load gene sets
# > import package
library(msigdbr)
print(paste0('MSigDb Version: ', packageVersion("msigdbr")))
# > gather gene sets
categories <- c("H", "C2", "C5", "C6", "C7")
gene_sets <- c()  # instantiate list of all gene sets
msigdbr_t2gs <- c()
for (category in categories) {
  msigdbr_df <- msigdbr(species = "Homo sapiens", category = category)
  msigdbr_list <- split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
  gene_sets <- c(gene_sets, msigdbr_list)
  msigdbr_t2g <- msigdbr_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
  msigdbr_t2gs <- rbind(msigdbr_t2gs, msigdbr_t2g)
}
rm(msigdbr_df, msigdbr_list)  # remove to free space

# retrieve the gsea object
library(GSEABase) 
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(ggplot2)
library(stringr)

# > read in the rankings
rankedGenes <- read.csv('/fh/fast/greenberg_p/user/dchen2/LINEAGE_TRACING/gsea/d16_vs_d0.csv', row.names = 1)
geneNames <- rownames(rankedGenes)
geneList <- rankedGenes[, 1]
names(geneList) <- geneNames
geneList <- sort(geneList, decreasing=TRUE)
geneList <- geneList[geneList != 0]
egmt <- GSEA(geneList, TERM2GENE=msigdbr_t2gs, verbose=F, pvalueCutoff=1, minGSSize=0)
egmtd <- data.frame(egmt)
write.csv(egmtd, file="/fh/fast/greenberg_p/user/dchen2/LINEAGE_TRACING/gsea/d16_vs_d0.gsea.csv")

# > read in the rankings
rankedGenes <- read.csv('/fh/fast/greenberg_p/user/dchen2/LINEAGE_TRACING/gsea/persists_vs_vanishes.csv', row.names = 1)
geneNames <- rownames(rankedGenes)
geneList <- rankedGenes[, 1]
names(geneList) <- geneNames
geneList <- sort(geneList, decreasing=TRUE)
geneList <- geneList[geneList != 0]
egmt <- GSEA(geneList, TERM2GENE=msigdbr_t2gs, verbose=F, pvalueCutoff=1, minGSSize=0)
egmtd <- data.frame(egmt)
write.csv(egmtd, file="/fh/fast/greenberg_p/user/dchen2/LINEAGE_TRACING/gsea/persists_vs_vanishes.gsea.csv")

# > read in the rankings
rankedGenes <- read.csv('/fh/fast/greenberg_p/user/dchen2/LINEAGE_TRACING/gsea/futureprolif_vs_nonprolif.csv', row.names = 1)
geneNames <- rownames(rankedGenes)
geneList <- rankedGenes[, 1]
names(geneList) <- geneNames
geneList <- sort(geneList, decreasing=TRUE)
geneList <- geneList[geneList != 0]
egmt <- GSEA(geneList, TERM2GENE=msigdbr_t2gs, verbose=F, pvalueCutoff=1, minGSSize=0)
egmtd <- data.frame(egmt)
write.csv(egmtd, file="/fh/fast/greenberg_p/user/dchen2/LINEAGE_TRACING/gsea/futureprolif_vs_nonprolif.gsea.csv")

# > filter the pathways for visualization
pathways <- c("GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP","HALLMARK_G2M_CHECKPOINT")
tmp <- msigdbr_t2gs[msigdbr_t2gs$gs_name %in% pathways, ]

library(enrichplot)
# > read in the rankings
rankedGenes <- read.csv('/fh/fast/greenberg_p/user/dchen2/LINEAGE_TRACING/gsea/persists_vs_vanishes.csv', row.names = 1)
geneNames <- rownames(rankedGenes)
geneList <- rankedGenes[, 1]
names(geneList) <- geneNames
geneList <- sort(geneList, decreasing=TRUE)
geneList <- geneList[geneList != 0]
egmt <- GSEA(geneList, TERM2GENE=tmp, verbose=F, pvalueCutoff=1, minGSSize=0, eps=0)
dev.off()
for (idx in 1:2) {
  fn <- paste0("/fh/fast/greenberg_p/user/dchen2/LINEAGE_TRACING/gsea/persists_vs_vanishes.",
               pathways[idx], ".svg")
  p <- gseaplot2(egmt, geneSetID = idx, title = pathways[idx],
                 base_size=10)
  print(p)
  ggsave(fn, width=20, units = "cm", scale = 1)
}

# > read in the rankings
rankedGenes <- read.csv('/fh/fast/greenberg_p/user/dchen2/LINEAGE_TRACING/gsea/futureprolif_vs_nonprolif.csv', row.names = 1)
geneNames <- rownames(rankedGenes)
geneList <- rankedGenes[, 1]
names(geneList) <- geneNames
geneList <- sort(geneList, decreasing=TRUE)
geneList <- geneList[geneList != 0]
egmt <- GSEA(geneList, TERM2GENE=tmp, verbose=F, pvalueCutoff=1, minGSSize=0, eps=0)
dev.off()
for (idx in 1:2) {
  fn <- paste0("/fh/fast/greenberg_p/user/dchen2/LINEAGE_TRACING/gsea/futureprolif_vs_nonprolif.",
               pathways[idx], ".svg")
  p <- gseaplot2(egmt, geneSetID = idx, title = pathways[idx],
                 base_size=10)
  print(p)
  ggsave(fn, width=20, units = "cm", scale = 1)
}