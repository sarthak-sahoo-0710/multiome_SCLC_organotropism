**Multi-Omics Analysis: RNA and ATAC-seq Integration of SCLC cell lines**

**Overview**
This analysis compares reference cell line (RA22-6) to test cases (Gen1LMD and Gen1BMD) to:
  - Normalize and preprocess RNA and ATAC sequencing data.
  - Identify differentially expressed genes (DEGs) and accessible chromatin regions.
  - Perform motif enrichment and pathway analysis.
  - Visualize UMAP embeddings, feature plots, and coverage plots.

**Citation**

**Title:** "Metastatic organotropism in small cell lung cancer"
**DOI:** 10.1101/2024.10.07.617066v2

**Features**

**Data Integration:** Merges RNA and ATAC-seq data into a single multi-modal Seurat object.

**Preprocessing:**
  - RNA: Normalized with SCTransform.
  - ATAC: Features filtered by nucleosome signal and TSS enrichment.

**Dimensionality Reduction:**
  - PCA for RNA.
  - LSI for ATAC.

**Cluster Analysis:** Identifies clusters with RNA and ATAC metrics.

**Differential Expression:**
  - RNA: Wilcoxon rank-sum test.
  - ATAC: Identifies regions linked to DEGs.

**Motif Enrichment:** Analyzes enriched motifs in accessible peaks using JASPAR2020.

**Pathway Analysis:** Computes AUC scores for custom gene sets.

**Visualization:**
  - UMAP projections for RNA and ATAC.
  - Feature plots for genes of interest.
  - Coverage plots for ATAC peaks.

**File Outputs**

**Expression and Accessibility Analysis:**
  - RNA_cluster_markers.txt: Differentially expressed RNA markers.
  - ATAC_cluster_markers.txt: Differentially accessible ATAC markers.

**Motif Enrichment:**
  - motif_enrichment_cluster_<cluster_id>.txt: Motif enrichment results for each cluster.

**Meta and Reduction Data:**
  - metadata.txt: Meta information for all cells.
  - reductions_umap.txt: UMAP coordinates (multi-modal).
  - reductions_umap_rna.txt: UMAP coordinates (RNA).
  - reductions_umap_atac.txt: UMAP coordinates (ATAC).

**Gene Links:**
  - links.txt: Peaks linked to DEGs.
  - closest_genes.txt: Closest genes to differential peaks.

**Visualizations**

**UMAP Plots:**
  - Group cells by their original identity (orig.ident).
  - Display RNA and ATAC clusters (wsnn_res.0.2).

**Feature Plots:**
  - Highlight gene expression for a given gene_of_interest.

**Coverage Plots:**
  - Show chromatin accessibility for specific regions and genes.
