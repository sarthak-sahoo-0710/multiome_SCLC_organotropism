# import required_packages
library(dplyr)
library(Signac)
library(Seurat)
library(Matrix)
library(Rmagic)
library(AUCell)
library(presto)
library(Biobase)
library(chromVAR)
library(GSEABase)
library(GEOquery)
library(ggplot2)
library(patchwork)
library(TFBSTools)
library(JASPAR2020)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# setup multiprocessing (if required)
plan("multicore", workers = 30)
options(future.globals.maxSize = 50000 * 1024^2, future.seed = NULL)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# for comparison of test case: Gen1LMD with reference: Gen1BMD
WT_obj <- Read10X_h5("./03_FilteredMatricesH5/SCAF3364_RA_22_6_filtered_feature_bc_matrix.h5")
WT_fragpath <- "./05_ATACFragments/SCAF3364_RA_22_6_atac_fragments.tsv.gz"
Bra1_obj <- Read10X_h5("./03_FilteredMatricesH5/SCAF3366_406_Brain_filtered_feature_bc_matrix.h5")
Bra1_fragpath <- "./05_ATACFragments/SCAF3366_406_Brain_atac_fragments.tsv.gz"
Liv1_obj <- Read10X_h5("./03_FilteredMatricesH5/SCAF3365_404_Liver_filtered_feature_bc_matrix.h5")
Liv1_fragpath <- "./05_ATACFragments/SCAF3365_404_Liver_atac_fragments.tsv.gz"
WT <- CreateSeuratObject(counts = WT_obj$`Gene Expression`,assay = "RNA",project = "RA22-6")
WT[["ATAC"]] <- CreateChromatinAssay(counts = WT_obj$Peaks,sep = c(":", "-"),fragments = WT_fragpath,annotation = annotation)
Bra1 <- CreateSeuratObject(counts = Bra1_obj$`Gene Expression`,assay = "RNA",project = "Gen1BMD")
Bra1[["ATAC"]] <- CreateChromatinAssay(counts = Bra1_obj$Peaks,sep = c(":", "-"),fragments = Bra1_fragpath,annotation = annotation)
Liv1 <- CreateSeuratObject(counts = Liv1_obj$`Gene Expression`,assay = "RNA",project = "Gen1LMD")
Liv1[["ATAC"]] <- CreateChromatinAssay(counts = Liv1_obj$Peaks,sep = c(":", "-"),fragments = Liv1_fragpath,annotation = annotation)

reference_data <- WT
test_data <- Bra1
test_data2 <- Liv1

# combine the data sets together
combined <- merge(x = reference_data, y=list(test_data, test_data2), add.cell.ids = c("RA22-6", "Gen1BMD", "Gen1LMD"))

# RNA based computes of cell features
combined <- PercentageFeatureSet(combined, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")

# ATAC based computes of cell features
peaks <- reduce(unlist(as(c(reference_data@assays$ATAC@ranges, test_data@assays$ATAC@ranges), "GRangesList")))
peakwidths <- width(peaks)
peaks <- peaks[peakwidths < 10000 & peakwidths > 20]
counts_atac_merged <- FeatureMatrix(combined@assays$ATAC@fragments,features = peaks,cells = colnames(combined), process_n=10000)
combined[['ATAC']] <- CreateChromatinAssay(counts_atac_merged,fragments = combined@assays$ATAC@fragments,annotation = combined@assays$ATAC@annotation,sep = c(":","-"))
combined <- NucleosomeSignal(combined, assay = "ATAC", process_n=10000)
combined <- TSSEnrichment(combined, assay = "ATAC")

# filtering based on both RNA and ATAC counts
combined <- subset(combined, subset = nFeature_RNA > 2500 & nFeature_RNA < 7500 & percent.mt < 30 & 
                     nFeature_ATAC > 1000 & nFeature_ATAC < 30000 & TSS.enrichment > 1 & nucleosome_signal < 2)

# perform RNA specific analysis
DefaultAssay(combined) <- "RNA"
combined <- SCTransform(combined)

genes_to_check <- c("ARHGEF9","ASCL1","NEUROD1","HNF4A","HNF1A","FOXA1","FOXA2","LMNA","SUN1",'SUN2','SPAG4','SYNE2','SYNE3','SYNE4','LMNB1','LMNB2','ACTB','TUBB','PLEC','TOR1AIP1','LBR','LEMD3','EMD','TMPO')
up_genes_downstream <- c('ARHGEF9', 'RUNX1T1', 'TRPM3', 'EDNRB-AS1', 'KANK4', 'TCERG1L', 'EML6', 'TNS3', 'CHRNA3', 'ANTXR2', 'LINC00698', 'TGFBR3', 'PPARGC1A', 'DHRS2', 'MIR100HG', 'CADM2', 'SOX6', 'SYNPO2', 'SEMA3C', 'CTNNA3', 'PPFIBP2', 'ZNF90', 'CYYR1', 'ATP2B3', 'CCDC152', 'ST6GALNAC5', 'SLC18A1', 'COLCA2', 'PDZRN4', 'GRM8', 'LINC01307', 'PRPH2', 'STMN2', 'ERICH3', 'CHL1', 'KCNK10', 'CERKL', 'ERVMER61-1', 'SLAIN1', 'UGT2A1', 'PAPPA', 'EPHA4', 'CADPS', 'MAP2', 'PARD3B', 'C1orf21', 'PPP4R4', 'OLFM2', 'PRKCE', 'IGFBPL1', 'COL4A3', 'BMP7', 'LINC00461', 'NTM', 'ADARB2', 'MFAP4', 'PLXNA2', 'PPP1R17', 'NEUROD1', 'ID3', 'WNT4', 'MAPK10', 'RYR3', 'ARRB1', 'ZEB2', 'SEMA3A', 'NHS', 'NR2F1', 'CDH1', 'SAMD11', 'OPRD1', 'LGR5', 'NTS', 'BRINP3', 'KALRN', 'CLSTN2', 'TSPAN5', 'CACNA1E', 'MET', 'ACVR1C', 'PTPRQ', 'RASSF9', 'GAS2', 'EPHA5', 'ELN', 'PTH2R', 'HEPACAM2', 'LINC00867', 'EYS', 'LGR6', 'COL4A4', 'SOX3', 'CTDSPL', 'GRIA3', 'FNDC5', 'FRMD6', 'CNTNAP4', 'ARHGEF3', 'LINC00616', 'YBX3', 'SLC7A11-AS1', 'ABCG1', 'DOCK5', 'SERPINI1', 'TLL1', 'L3MBTL4', 'PKD1L3', 'POU3F3', 'MIR137HG', 'PTPRM', 'ARHGAP22', 'SLC1A5', 'COL4A6', 'PCSK1', 'MYO5B', 'CR1L', 'DISP3', 'ZNF486', 'TMEM255A', 'LSAMP', 'SSTR2', 'DSP', 'EYA2', 'MAP6', 'HPCA', 'MAN1A1', 'LAMA2', 'LRFN5', 'TGFBR2', 'CDH20', 'CNTN1', 'KCND3', 'NHLH1', 'GRIN2B', 'DMKN', 'SYTL4', 'COL1A2', 'LRRTM3', 'PHLDA1', 'CNTN2', 'SNTB1', 'PKHD1', 'LURAP1L', 'PRSS12', 'GAP43', 'OSBP2', 'AFF3')
subset_df_SCT <- as.data.frame(combined@assays$SCT@data) %>% dplyr::filter(row.names(combined@assays$SCT@data) %in% genes_to_check)
subset_df_SCT_up <- as.data.frame(combined@assays$SCT@data) %>% dplyr::filter(row.names(combined@assays$SCT@data) %in% genes_to_check)
#write.table(subset_df_SCT, "./experimental_github/Gen1LMD_BMD/upregulated_genes_imputed_subset_genes_SCT.tsv", sep='\t', col.names=TRUE, row.names = TRUE,quote = FALSE)


combined <- RunPCA(combined)
combined <- RunUMAP(object = combined, reduction.name = 'UMAP_RNA', reduction.key = 'UMAPRNA_', dims = 1:40)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.2, cluster.name = 'RNA_snn_res_0.2')

# compute specific imputed genes
DefaultAssay(combined) <- "SCT"
data_matrix <- as.matrix(combined@assays$SCT@data)
data_MAGIC <- magic(data_matrix,solver='approximate')
df_result <- as.data.frame(data_MAGIC$result)
dim(data_MAGIC$result)
subset_df<-df_result %>% dplyr::filter(row.names(df_result) %in% genes_to_check)
subset_df_up<-df_result %>% dplyr::filter(row.names(df_result) %in% up_genes_downstream)
#write.table(subset_df_up, "./experimental_github/Gen1LMD_BMD/upregulated_genes_imputed_subset_genes.tsv", sep='\t', col.names=TRUE, row.names = TRUE,quote = FALSE)


DefaultAssay(combined) <- "SCT"
# compute pathway scores on SCT normalised values
geneSets <- getGmt("./data_27_10_24/IMPOWER_pathways.gmt")
geneSets <- subsetGeneSets(geneSets, rownames(combined)) 
cbind(nGenes(geneSets))
counts_matrix <- LayerData(combined, assay = "SCT")
cells_rankings <- AUCell_buildRankings(counts_matrix, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
AUCmatrix <-  as.data.frame(cells_AUC@assays@data$AUC)

# run preprocessing for ATAC - elbow plot
DefaultAssay(combined) <- "ATAC"
combined <- FindTopFeatures(combined, min.cutoff = 5)
combined <- RunTFIDF(combined)
combined <- RunSVD(combined)

# UMAP and cluster analysis for ATAC
combined <- RunUMAP(object = combined, reduction = 'lsi', reduction.name = 'UMAP_ATAC', reduction.key = 'UMAPATAC_', dims = 2:40)
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:20)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3, resolution = 0.2, cluster.name = 'ATAC_snn_res_0.2')

rm(Liv1,Liv1_obj,Bra1,Bra1_obj, WT, WT_obj, reference_data, test_data)

combined <- FindMultiModalNeighbors(combined, reduction.list = list("pca", "lsi"),
dims.list = list(1:ncol(Embeddings(combined,"pca")),1:ncol(Embeddings(combined,"lsi"))),
modality.weight.name = c("RNA.weight","ATAC.weight"),verbose = TRUE)

combined <- RunUMAP(combined, nn.name = "weighted.nn", assay = "RNA")
combined <- FindClusters(combined, graph.name = "wsnn", resolution = 0.2)

DefaultAssay(combined) <- "SCT"

combined <- PrepSCTFindMarkers(combined, assay = 'SCT')

Idents(combined) <- combined$wsnn_res.0.2

# Perform differential expression analysis using the SCT assay
DE_ct <- FindAllMarkers(combined, 
                     test.use = "wilcox",   # Use Wilcoxon rank-sum test
                     assay = "SCT",         # Use SCT assay for normalization
                     min.pct = 0.1)         # Keep markers expressed in at least 10% of cells

# Filter based on AUC, p-value, and other conditions
top_markers_ct <- DE_ct %>% dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.01)

DefaultAssay(combined) <- "ATAC"
DA_ct <- FindAllMarkers(combined, 
                        test.use = "wilcox",   # Use Wilcoxon rank-sum test
                        assay = "ATAC",         # Use SCT assay for normalization
                        min.pct = 0.1)         # Keep markers expressed in at least 10% of cells

# Filter based on AUC, p-value, and other conditions
top_markers_ct_ATAC <- DA_ct %>% dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.01)

combined <- RegionStats(combined, genome = BSgenome.Hsapiens.UCSC.hg38)
combined <- LinkPeaks(combined, peak.assay = "ATAC", expression.assay = "SCT", genes.use = top_markers_ct$gene)


pfm <- getMatrixSet(x = JASPAR2020,opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
df_pfm <- data.frame(t(sapply(pfm, function(x) c(id=x@ID, name=x@name, symbol=ifelse(!is.null(x@tags$symbol),x@tags$symbol,NA)))))
combined <- AddMotifs(combined, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

open_peaks <- AccessiblePeaks(combined)
peaks_matched <- MatchRegionStats(meta.feature = combined[['ATAC']]@meta.features[open_peaks, ],
                                  query.feature = combined[['ATAC']]@meta.features[top_markers_ct_ATAC$gene, ],
                                  n = 50000)

closest_genes <- ClosestFeature(combined, regions = row.names(top_markers_ct_ATAC))

unique_clusters <- unique(top_markers_ct_ATAC$cluster)

# Initialize an empty list to store results for each cluster
enriched_motif_list <- list()

# Loop over each cluster
for (cluster_id in unique_clusters) {
  # Perform motif enrichment for the current cluster
  motif_enrichment <- FindMotifs(combined,
                                 features = top_markers_ct_ATAC$gene[top_markers_ct_ATAC$cluster == cluster_id],
                                 background = peaks_matched) %>%
    mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
    mutate(padj = p.adjust(pvalue, method="BH"))
  
  # Store the result in the list with the cluster_id as the name
  enriched_motif_list[[as.character(cluster_id)]] <- motif_enrichment
  
  # Define the output file name
  output_file <- paste0("./experimental_github/Gen1LMD_BMD/motif_enrichment_cluster_", cluster_id, ".txt")
  
  # Write the results to a file
  #write.table(motif_enrichment, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Optionally, print a message to track progress
  message("Written results for cluster ", cluster_id, " to ", output_file)
  
}

#write.table(combined@assays$ATAC@links, file = "./experimental_github/Gen1LMD_BMD/links.txt", quote = FALSE, sep = "\t", row.names = TRUE)
#write.table(AUCmatrix, file = "./experimental_github/Gen1LMD_BMD/IMPOWER.txt", quote = FALSE, sep = "\t", row.names = TRUE)
#write.table(subset_df, "./experimental_github/Gen1LMD_BMD/LINC_imputed_subset_genes_subset_genes.tsv", sep='\t', col.names=TRUE, row.names = TRUE,quote = FALSE)
#write.table(subset_df_SCT, "./experimental_github/Gen1LMD_BMD/LINC_imputed_subset_genes_SCT.tsv", sep='\t', col.names=TRUE, row.names = TRUE,quote = FALSE)
#write.csv(combined@meta.data, file = "./experimental_github/Gen1LMD_BMD/metadata.txt", quote = FALSE)
#write.csv(combined@reductions$umap@cell.embeddings, file = "./experimental_github/Gen1LMD_BMD/reductions_umap.txt", quote = FALSE)
#write.csv(combined@reductions$UMAP_RNA@cell.embeddings, file = "./experimental_github/Gen1LMD_BMD/reductions_umap_rna.txt", quote = FALSE)
#write.csv(combined@reductions$UMAP_ATAC@cell.embeddings, file = "./experimental_github/Gen1LMD_BMD/reductions_umap_atac.txt", quote = FALSE)
#write.table(DE_ct, file = "./experimental_github/Gen1LMD_BMD/RNA_cluster_markers.txt", quote = FALSE, sep = "\t", row.names = TRUE)
#write.table(DA_ct, file = "./experimental_github/Gen1LMD_BMD/ATAC_cluster_markers.txt", quote = FALSE, sep = "\t", row.names = TRUE)
#write.table(closest_genes, file = "./experimental_github/Gen1LMD_BMD/closest_genes.txt", quote = FALSE, sep = "\t", row.names = TRUE)

gene_of_interest <- 'MET'
p1 <- DimPlot(combined,
        group.by = "orig.ident",
        reduction = "umap") & NoAxes()
p2 <- DimPlot(combined,
              group.by = "wsnn_res.0.2",
              reduction = "umap") & NoAxes()
p3 <- FeaturePlot(combined,
              c(gene_of_interest),
              reduction = "umap") & NoAxes()

DefaultAssay(combined) <- "ATAC"
gene_of_interest <- 'ARHGEF9'
p4<-CoveragePlot(combined,
             region = gene_of_interest,
             features = gene_of_interest,
             group.by = "wsnn_res.0.2",
             extend.upstream = 300000,
             ymax=200,
             extend.downstream = 50, expression.assay = 'SCT')
p5 <- p4 & scale_fill_manual(values = c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", 
                                        "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5"))
p5 & theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
           axis.text.x = element_text(size=10), axis.text.y = element_text(size=8), title = element_text(size=12))
p5

(p1|p2|p3)/p4

MotifPlot(combined, motifs = c('Sox6','HAND2','HES1','NEUROG2','KLF15','PDX1','HOXA1','HOXA2'), ncol=4)

