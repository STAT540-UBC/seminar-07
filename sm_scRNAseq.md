STAT540 - Seminar 7: Single-Cell RNA-seq Analysis with Seurat
================

# Introduction

Single-cell RNA sequencing (scRNA-seq) allows us to profile gene
expression at the resolution of individual cells, revealing cellular
heterogeneity that bulk RNA-seq cannot capture. In this seminar, we’ll
analyze a dataset of Peripheral Blood Mononuclear Cells (PBMCs) using
the Seurat package in R.

We will cover:

- Quality control and filtering of cells
- Normalization and identification of highly variable features
- Dimensionality reduction (PCA, UMAP)
- Clustering and cell type identification
- Finding differentially expressed marker genes

## Learning Objectives

By the end of this seminar, you will be able to:

1.  Perform quality control on scRNA-seq data
2.  Normalize and scale gene expression data
3.  Identify cell populations through clustering
4.  Visualize single-cell data using dimensionality reduction
5.  Identify marker genes that define cell populations

## Dataset

We’ll use a downsampled version of the PBMC3k dataset, which contains
approximately 2,700 single cells isolated from a healthy donor. This
dataset has been downsampled to make it feasible to run on a typical
laptop with 8 GB RAM.

# Setup

First, let’s load the required libraries:

``` r
# Install Seurat if needed
#install.packages('Seurat')
#install.packages('dplyr')
#install.packages('ggplot2')
#install.packages('remotes')
library(Seurat)
library(dplyr)
library(ggplot2)
library(remotes)
# Install SeuratData from GitHub
remotes::install_github('satijalab/seurat-data')
```

# Loading the Data

Seurat works with expression matrices where rows represent genes and
columns represent cells. We’ll load the PBMC dataset directly from the
SeuratData package.

``` r
library(SeuratData)
# Install the dataset if not already installed
InstallData("pbmc3k")
# Load the dataset
pbmc3k <- LoadData("pbmc3k")
# Downsample the data to ~1000 cells for faster processing
set.seed(42)
pbmc3k <- subset(pbmc3k, cells = sample(colnames(pbmc3k), 1000))
```

Alternatively, if you have your own data in the 10X Genomics format:

``` r
# Read in 10X data
# pbmc.data <- Read10X(data.dir = "path/to/filtered_gene_bc_matrices/hg19/")
# pbmc3k <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```

For this seminar, we’ll work with the pre-loaded `pbmc3k` object.

# Quality Control

Not all cells in our dataset are high quality. Some may be dying cells,
doublets (two cells captured together), or have low RNA content. We need
to filter these out.

## QC Metrics

Common quality control metrics include:

- *nFeature_RNA*: Number of genes detected per cell
- *nCount_RNA*: Total number of molecules (UMIs) detected per cell
- *Percent mitochondrial genes*: Proportion of reads mapping to
  mitochondrial genes (high values indicate dying cells)

``` r
# Calculate percentage of mitochondrial genes
# Mitochondrial genes start with "MT-" in human data
pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc3k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

![](sm_scRNAseq_files/figure-gfm/qc_metrics-1.png)<!-- -->

We can also visualize relationships between QC metrics:

``` r
# Feature-count relationship
plot1 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

![](sm_scRNAseq_files/figure-gfm/qc_scatter-1.png)<!-- -->

``` r
# Calculate correlation between nCount and nFeature
cor_value <- cor(pbmc3k$nCount_RNA, pbmc3k$nFeature_RNA)
cat("Correlation between nCount_RNA and nFeature_RNA:", round(cor_value, 3), "\n")
```

    ## Correlation between nCount_RNA and nFeature_RNA: 0.953

**Key observation**: There is typically a strong positive correlation
(\>0.95) between `nCount_RNA` (total UMIs) and `nFeature_RNA` (number of
genes detected). This makes sense biologically - cells with more total
UMIs generally detect more unique genes.

Some exceptions include:

- **Abnormally high UMIs and high genes**: Likely doublets/multiplets
  (two or more cells captured together)
- **Abnormally low UMIs and low genes**: Empty droplets, dying cells, or
  ambient RNA contamination

Additionally, cells that deviate from the linear trend (e.g., unusually
high UMIs for a given number of genes) may indicate technical artifacts
or biological outliers worth investigating.

## Doublet Detection

Before filtering based on QC thresholds, we can use computational
methods to identify doublets. DoubletFinder works best when run on the
full dataset before aggressive filtering. We’ll run it first, then apply
QC filtering.

``` r
# Install DoubletFinder if needed
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

# Pre-process for doublet detection
pbmc3k <- NormalizeData(pbmc3k)
pbmc3k <- FindVariableFeatures(pbmc3k, selection.method = "vst", nfeatures = 2000)
pbmc3k <- ScaleData(pbmc3k)
pbmc3k <- RunPCA(pbmc3k)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10)
pbmc3k <- FindClusters(pbmc3k)

# Find optimal pK value
sweep.res <- paramSweep(pbmc3k, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Use pK with highest BCmetric
pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Estimate doublet rate (assume 7.5% for ~1000 cells)
nExp_poi <- round(0.075 * nrow(pbmc3k@meta.data))

# Run DoubletFinder
pbmc3k <- doubletFinder(pbmc3k, PCs = 1:10, pN = 0.25, pK = pK_value,
                        nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Get the doublet classification column name (it's auto-generated)
doublet_col <- grep("^DF.classifications", colnames(pbmc3k@meta.data), value = TRUE)

# Visualize doublets on UMAP
DimPlot(pbmc3k, reduction = "umap", group.by = doublet_col)

# Count doublets
table(pbmc3k@meta.data[[doublet_col]])

# Keep doublet classification in metadata for filtering later
pbmc3k$doublet_class <- pbmc3k@meta.data[[doublet_col]]
```

**Note**: DoubletFinder should be run **before** aggressive QC
filtering. It identifies doublets by creating artificial doublets and
comparing them to real cells. Running it on the full dataset provides
better context for doublet identification. We’ll filter out doublets
along with low-quality cells in the next step.

**This step is set to `eval=FALSE` by default** since it requires
additional package installation and significantly increases runtime. For
this tutorial, we’ll proceed with the feature count-based doublet
filtering.

## Filtering Cells

Based on the QC metrics (and optionally doublet detection), we filter
out low-quality cells:

``` r
# Filter cells based on QC metrics
# If you ran doublet detection above, also add: & doublet_class == "Singlet"
pbmc3k <- subset(pbmc3k, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

These thresholds mean we keep cells with: \* More than 200 genes
detected (removes empty droplets) \* Fewer than 2500 genes (removes
potential doublets) \* Less than 5% mitochondrial content (removes dying
cells)

# Normalization

Raw count data must be normalized to account for differences in
sequencing depth between cells. Seurat uses a log-normalization method
by default.

``` r
# Normalize the data
pbmc3k <- NormalizeData(pbmc3k, normalization.method = "LogNormalize", scale.factor = 10000)
```

This approach normalizes gene expression for each cell by total
expression, multiplies by a scale factor (10,000), and log-transforms
the result.

# Feature Selection

Not all genes are informative for distinguishing cell types. We identify
highly variable features (genes with high cell-to-cell variation) to
focus downstream analysis on the most informative genes.

``` r
# Identify highly variable features
pbmc3k <- FindVariableFeatures(pbmc3k, selection.method = "vst", nfeatures = 2000)

# Identify the top 10 most variable genes
top10 <- head(VariableFeatures(pbmc3k), 10)
top10
```

    ##  [1] "S100A9"  "LYZ"     "GNLY"    "FTL"     "FTH1"    "PPBP"    "S100A8" 
    ##  [8] "HLA-DRA" "IGLL5"   "CD74"

``` r
# Plot variable features
plot1 <- VariableFeaturePlot(pbmc3k)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
```

![](sm_scRNAseq_files/figure-gfm/find_variable_features-1.png)<!-- -->

# Scaling

Before dimensionality reduction, we scale the data so that highly
expressed genes don’t dominate.

``` r
# Scale the data
all.genes <- rownames(pbmc3k)
pbmc3k <- ScaleData(pbmc3k, features = all.genes)
```

Scaling transforms the data so each gene has mean expression of 0 and
variance of 1 across cells.

# Dimensionality Reduction: PCA

scRNA-seq data is high-dimensional (thousands of genes). Principal
Component Analysis (PCA) reduces this complexity while retaining
important variation.

``` r
# Run PCA
pbmc3k <- RunPCA(pbmc3k, features = VariableFeatures(object = pbmc3k))

# Examine the PCA results
print(pbmc3k[["pca"]], dims = 1:5, nfeatures = 5)
```

    ## PC_ 1 
    ## Positive:  CST3, TYROBP, LST1, AIF1, FTH1 
    ## Negative:  MALAT1, IL32, CD3D, LTB, AES 
    ## PC_ 2 
    ## Positive:  CD79A, HLA-DQA1, MS4A1, HLA-DQB1, HLA-DRA 
    ## Negative:  NKG7, GZMB, GZMA, PRF1, CST7 
    ## PC_ 3 
    ## Positive:  SDPR, GP9, PF4, PPBP, TMEM40 
    ## Negative:  S100A4, MALAT1, CYBA, S100A10, TSPO 
    ## PC_ 4 
    ## Positive:  CD79A, CD79B, HLA-DQA1, HLA-DQB1, CD74 
    ## Negative:  CD3D, NOSIP, VIM, S100A6, FYB 
    ## PC_ 5 
    ## Positive:  S100A8, GZMB, MS4A6A, LGALS2, NKG7 
    ## Negative:  LTB, HES4, MS4A4A, CDKN1C, CKB

``` r
# Visualize PCA
VizDimLoadings(pbmc3k, dims = 1:2, reduction = "pca")
```

![](sm_scRNAseq_files/figure-gfm/pca-1.png)<!-- -->

``` r
DimPlot(pbmc3k, reduction = "pca")
```

![](sm_scRNAseq_files/figure-gfm/pca-2.png)<!-- -->

## Determining Dimensionality

How many principal components should we use? We can use an elbow plot to
identify where additional PCs contribute little variation:

``` r
ElbowPlot(pbmc3k)
```

![](sm_scRNAseq_files/figure-gfm/elbow-1.png)<!-- -->

The “elbow” in the plot suggests where to make the cutoff. For this
dataset, PCs 1-10 capture most of the variation.

# Clustering

Now we cluster cells based on their gene expression profiles in PCA
space. Similar cells should cluster together, potentially representing
the same cell type.

``` r
# Find neighbors in PCA space
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10)

# Find clusters
pbmc3k <- FindClusters(pbmc3k, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 977
    ## Number of edges: 34678
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8524
    ## Number of communities: 6
    ## Elapsed time: 0 seconds

``` r
# Look at cluster assignments
head(Idents(pbmc3k), 5)
```

    ## AAACATACAACCAC AAACGCTGTAGCCA AAAGAGACGAGATA AAAGAGACGCGAGA AAAGAGACGGACTT 
    ##              0              3              0              1              3 
    ## Levels: 0 1 2 3 4 5

The `resolution` parameter controls cluster granularity. Higher values
produce more clusters.

# Non-linear Dimensionality Reduction: UMAP

While PCA is useful, visualizing data in 2D requires further
dimensionality reduction. UMAP (Uniform Manifold Approximation and
Projection) preserves both local and global structure.

``` r
# Run UMAP
pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)

# Visualize clusters
DimPlot(pbmc3k, reduction = "umap", label = TRUE)
```

![](sm_scRNAseq_files/figure-gfm/umap-1.png)<!-- -->

Each point represents a cell, colored by cluster assignment. Cells of
the same type should group together.

# Finding Marker Genes

To understand what cell types our clusters represent, we identify genes
that are differentially expressed in each cluster (marker genes).

## Markers for All Clusters

``` r
# Find markers for every cluster compared to all remaining cells
pbmc3k.markers <- FindAllMarkers(pbmc3k, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top markers for each cluster
pbmc3k.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
```

    ## # A tibble: 12 x 7
    ## # Groups:   cluster [6]
    ##       p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene     
    ##       <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>    
    ##  1 2.97e-33       4.84 0.275 0.017  4.07e-29 0       MAL      
    ##  2 1.05e-31       3.67 0.315 0.039  1.43e-27 0       AQP3     
    ##  3 8.43e-51       8.33 0.28  0.001  1.16e-46 1       FOLR3    
    ##  4 1.29e-48       7.61 0.275 0.003  1.77e-44 1       S100A12  
    ##  5 1.26e-36       7.79 0.25  0.011  1.72e-32 2       IGLL5    
    ##  6 5.16e-92       7.57 0.525 0.011  7.08e-88 2       LINC00926
    ##  7 2.67e-52       4.34 0.539 0.059  3.67e-48 3       GZMK     
    ##  8 3.55e-34       3.97 0.441 0.065  4.87e-30 3       GZMH     
    ##  9 1.54e-40       6.48 0.258 0.005  2.12e-36 4       SH2D1B   
    ## 10 7.35e-70       6.27 0.919 0.131  1.01e-65 4       GNLY     
    ## 11 1.34e-80       5.65 0.475 0.007  1.83e-76 5       CDKN1C   
    ## 12 3.23e-48       5.56 0.311 0.007  4.44e-44 5       CKB

Parameters: \* `only.pos = TRUE`: Only report positive markers
(upregulated genes) \* `min.pct = 0.25`: Gene must be detected in at
least 25% of cells in either group \* `logfc.threshold = 0.25`: Minimum
log fold-change threshold

## Visualizing Marker Expression

``` r
# Violin plot for specific markers
VlnPlot(pbmc3k, features = c("MS4A1", "CD79A"))
```

![](sm_scRNAseq_files/figure-gfm/visualize_markers-1.png)<!-- -->

``` r
# Feature plot showing expression on UMAP
FeaturePlot(pbmc3k, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
```

![](sm_scRNAseq_files/figure-gfm/visualize_markers-2.png)<!-- -->

## Heatmap of Top Markers

``` r
# Get top 10 markers per cluster
top10_markers <- pbmc3k.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

DoHeatmap(pbmc3k, features = top10_markers$gene) + NoLegend()
```

![](sm_scRNAseq_files/figure-gfm/marker_heatmap-1.png)<!-- -->

# Cell Type Annotation

Based on known marker genes, we can assign cell type identities to
clusters:

- Cluster 0: Naive CD4+ T cells (IL7R, CCR7)
- Cluster 1: CD14+ Monocytes (CD14, LYZ)
- Cluster 2: Memory CD4+ T cells (IL7R, S100A4)
- Cluster 3: B cells (MS4A1/CD20)
- Cluster 4: CD8+ T cells (CD8A)
- Cluster 5: FCGR3A+ Monocytes (FCGR3A, MS4A7)
- Cluster 6: NK cells (GNLY, NKG7)
- Cluster 7: Dendritic Cells (FCER1A, CST3)
- Cluster 8: Platelets (PPBP)

``` r
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc3k)
pbmc3k <- RenameIdents(pbmc3k, new.cluster.ids)
DimPlot(pbmc3k, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

![](sm_scRNAseq_files/figure-gfm/annotate-1.png)<!-- -->

# Finding Markers Between Specific Clusters

Sometimes we want to compare two specific cell populations:

``` r
# Find markers distinguishing cluster 1 from cluster 5 (two monocyte populations)
cluster1.markers <- FindMarkers(pbmc3k, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", min.pct = 0.25)
head(cluster1.markers, n = 5)
```

    ##               p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## FCGR3A 7.407410e-39  -5.169201 0.125 0.934 1.015852e-34
    ## LYZ    5.851566e-30   2.792958 1.000 0.967 8.024837e-26
    ## CDKN1C 1.742245e-23  -4.494179 0.005 0.475 2.389315e-19
    ## S100A8 2.470944e-23   4.289250 0.935 0.492 3.388653e-19
    ## S100A9 1.370778e-22   3.551349 0.975 0.820 1.879886e-18

# Deliverables

Complete the following two deliverables and submit your knitted R
Markdown file (both .Rmd and .md files).

## Deliverable 1: Quality Control and Parameter Exploration (1 point)

Perform the following analyses and provide written interpretations:

1.  **Dimensionality Selection**: Based on the elbow plot, at
    approximately which PC does the curve start to flatten? Would you
    choose 10, 15, or 20 PCs for downstream analysis? Explain your
    reasoning (2-3 sentences).

**Your written response here:**

## Deliverable 2: Marker Gene Analysis and Biological Interpretation (1 point)

1.  **Marker Gene Validation**: Examine the top marker genes for cluster
    0 from the `pbmc3k.markers` object. Use `FeaturePlot()` to visualize
    expression of the top 3 marker genes from this cluster across the
    UMAP. Based on their expression patterns, do these genes seem to be
    good markers? Explain what makes a gene a “good” marker (3-4
    sentences).

``` r
# Your code here for marker visualization and analysis
```

**Your written response here:**

# Saving Your Analysis

``` r
# Save the Seurat object for later use
saveRDS(pbmc3k, file = "pbmc3k_final.rds")

# To load it later:
# pbmc3k <- readRDS("pbmc3k_final.rds")
```

# Summary

In this seminar, we covered a standard scRNA-seq analysis workflow:

1.  **Quality Control**: Filtering low-quality cells based on gene
    counts and mitochondrial content
2.  **Normalization**: Accounting for sequencing depth differences
    between cells
3.  **Feature Selection**: Identifying highly variable genes
4.  **Dimensionality Reduction**: Using PCA and UMAP to visualize and
    analyze high-dimensional data
5.  **Clustering**: Grouping similar cells together
6.  **Marker Identification**: Finding genes that define each cluster
7.  **Annotation**: Assigning biological meaning to clusters

These steps form the foundation for more advanced analyses like
trajectory inference, cell-cell communication, and integration of
multiple datasets.

# References

1.  Hao, Y., et al. (2021). Integrated analysis of multimodal
    single-cell data. *Cell*, 184(13), 3573-3587.
2.  Stuart, T., et al. (2019). Comprehensive integration of single-cell
    data. *Cell*, 177(7), 1888-1902.
3.  Satija Lab. Seurat PBMC3k Tutorial.
    <https://satijalab.org/seurat/articles/pbmc3k_tutorial>
4.  Satija Lab. Integration and Label Transfer.
    <https://satijalab.org/seurat/articles/integration_mapping>

# Acknowledgments

This tutorial was developed with assistance from [Claude
Code](https://claude.com/claude-code), an AI-powered coding assistant by
Anthropic.
