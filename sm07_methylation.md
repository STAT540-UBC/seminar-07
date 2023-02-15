STAT540 - Seminar 7: DNA methylation analysis with Illumina Infinium DNA
microarrays
================

## Attributions

This seminar was developed by Keegan Korthauer, and is modeled after the
Bioconductor Workflow [“**A cross-package Bioconductor workflow for
analysing methylation array
data**”](https://bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html)
by Jovana Maksimovic, Belinda Phipson and Alicia Oshlack. We skip over
some details that are included in the workflow, and add in some details
not included in the workflow. *As we’ve seen with other analysis types
in this course (e.g. RNA-seq), there are several sensible options for
analyzing DNA microarray data*.

## Learning objectives

By the end of this tutorial, you should be able to  
- Identify the basic steps in a methylation array analysis pipeline  
- Appreciate the role of normalization prior to analysis in methylation
datasets  
- Explain the difference between beta and M-values, and which are
appropriate for various steps of the anlaysis  
- Undertake exploratory visual analysis of methylation data (plotting
Beta value densities, PCA plots, sample correlation plots)  
- Undertake differential methylation analysis of **probes and regions**
and annotate probes with genomic information (chromosome, position, CpG
island status, nearest gene, etc) for interpretation  
- Plot results from differential analysis (differentially methylated
probes and regions) to highlight areas under regulatory influence

## Setup - install and load packages

First we load necessary packages into our session. If any of these
packages are not already installed on your system, you’ll first need to
install them with `BiocManager::install()`. You likely don’t already
have the Bioconductor Workflow package (`methylationArrayAnalysis`),
`DMRcate`, or `DMRcatedata` pre-installed, so we’ve included those
commands in the following chunk. The first package includes the dataset
we’ll work with, so it involves a sizeable download (on the order of
hundreds of MB) and may take several minutes depending on your
connection. We’ll also make heavy use of the `minfi`, `limma`, and
`DMRcate` Bioconductor packages, which are automatically loaded by the
workflow package.

First install `methylationArrayAnalysis` and `DMRcatedata` by running
this command in your console (only need to run it once, not every time
you knit, so chunk is set to `eval = FALSE`).

``` r
BiocManager::install("methylationArrayAnalysis")
BiocManager::install("DMRcate")
BiocManager::install("DMRcatedata")
```

Install any other packages below before loading the libraries if you
don’t already have them. Note that the loading of these packages may
take a couple of minutes, so be patient!

``` r
library(methylationArrayAnalysis)
library(DMRcate)
library(DMRcatedata)
library(tidyverse)
library(ggplot2)
theme_set(theme_bw())
library(ComplexHeatmap)
library(circlize)
set.seed(540)
```

## Introduction

The Illumina Infiniumn Array is a microarray-based high throughput
platform for methylation profiling on more than 450,000 pre-selected
probes for CpGs across the human genome. This platform is an aggregation
of the earlier Illumina HumanMethylation27 (“27K”) and
HumanMethylation450 (“450K”) arrays. There also exists an “850K” array
called the Illumina MethylationEPIC that houses more than 850,000
probes.

On these arrays, each CpG is measured by two values: a methylated
intensity (M) and an unmethylated intensity (U). Methylation levels are
commonly reported as either beta values (beta=M/(M+U)) or M-values
(Mvalue=log2(M/U)=log2(beta/(1-beta))), which are the logit transformed
beta values. A small offset is often added to the denominator of the
beta value to avoid dividing by small values. Beta values are readily
interpretable as a methylation percentage, and are generally preferable
for visualization and effect size interpretation. M-values, on the other
hand, are more appropriate for statistical testing due to their
distributional properties (see [Du et
al. 2010](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587)
for an excellent primer on why).

## Read in 450K methylation array data

A great number of these types of datasets have been made publicly
available through Gene Expression Omnibus (GEO). However, they are often
available most conveniently as Bioconductor objects in the form of
processed datasets, e.g. beta values already computed. We would like to
demonstrate the analysis steps from the raw signal intensity data (the M
and U values). To avoid the hassle of downloading the raw signal
intensity data and formatting it into a Bioconductor object, we’ll use
the data provided in the `methylationArrayAnalysis` package.

In this seminar, we are going to perform differential methylation
analysis on sorted T-cell types from [Zhang et
al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3798997/). This
data is publicly available at GEO accession
[GSE49667](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49667),
but as mentioned above, we’ll save time by loading the raw data in
already formatted objects. The dataset contains 10 samples: naive,
rTreg, act_naive, act_rTreg, collected from 3 different individuals
(M28, M29, M30).

Let’s take a look at the data files we have.

``` r
# set up a path to the data directory
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
# list the files
list.files(dataDirectory, recursive = TRUE)
```

    ##  [1] "48639-non-specific-probes-Illumina450k.csv"
    ##  [2] "5975827018/5975827018_R06C02_Grn.idat"     
    ##  [3] "5975827018/5975827018_R06C02_Red.idat"     
    ##  [4] "6264509100/6264509100_R01C01_Grn.idat"     
    ##  [5] "6264509100/6264509100_R01C01_Red.idat"     
    ##  [6] "6264509100/6264509100_R01C02_Grn.idat"     
    ##  [7] "6264509100/6264509100_R01C02_Red.idat"     
    ##  [8] "6264509100/6264509100_R02C01_Grn.idat"     
    ##  [9] "6264509100/6264509100_R02C01_Red.idat"     
    ## [10] "6264509100/6264509100_R02C02_Grn.idat"     
    ## [11] "6264509100/6264509100_R02C02_Red.idat"     
    ## [12] "6264509100/6264509100_R03C01_Grn.idat"     
    ## [13] "6264509100/6264509100_R03C01_Red.idat"     
    ## [14] "6264509100/6264509100_R03C02_Grn.idat"     
    ## [15] "6264509100/6264509100_R03C02_Red.idat"     
    ## [16] "6264509100/6264509100_R04C01_Grn.idat"     
    ## [17] "6264509100/6264509100_R04C01_Red.idat"     
    ## [18] "6264509100/6264509100_R04C02_Grn.idat"     
    ## [19] "6264509100/6264509100_R04C02_Red.idat"     
    ## [20] "6264509100/6264509100_R05C01_Grn.idat"     
    ## [21] "6264509100/6264509100_R05C01_Red.idat"     
    ## [22] "6264509100/6264509100_R05C02_Grn.idat"     
    ## [23] "6264509100/6264509100_R05C02_Red.idat"     
    ## [24] "6264509100/6264509100_R06C01_Grn.idat"     
    ## [25] "6264509100/6264509100_R06C01_Red.idat"     
    ## [26] "6264509100/6264509100_R06C02_Grn.idat"     
    ## [27] "6264509100/6264509100_R06C02_Red.idat"     
    ## [28] "ageData.RData"                             
    ## [29] "human_c2_v5.rdata"                         
    ## [30] "model-based-cpg-islands-hg19-chr17.txt"    
    ## [31] "SampleSheet.csv"                           
    ## [32] "wgEncodeRegDnaseClusteredV3chr17.bed"

We can see several pairs of ‘Red’ and ‘Grn’ IDAT (Intensity Data) files.
These are the raw signal intensities for the read and green channels
(which can be parsed to obtain the relative methylated and unmethylated
signals). There are more than 10 pairs, since there are three extra
samples included that we will not consider here. There is also a file
called `SampleSheet.csv`, which contains sample metadata linking the
IDAT files to samples. This file is formatted with one row per sample,
and one column per metadata/phenotypic data variable. We will use the
`read.metharray.sheet` sheet from the `minfi` package to read this info
in. We’ll then remove an extra sample from another dataset that is
included but we’ll not consider in this seminar.

``` r
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
```

    ## [read.metharray.sheet] Found the following CSV files:

    ## [1] "/Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/SampleSheet.csv"

``` r
targets
```

    ##    Sample_Name Sample_Well   Sample_Source Sample_Group Sample_Label Pool_ID
    ## 1            1          A1             M28        naive        naive    <NA>
    ## 2            2          B1             M28        rTreg        rTreg    <NA>
    ## 3            3          C1             M28    act_naive    act_naive    <NA>
    ## 4            4          D1             M29        naive        naive    <NA>
    ## 5            5          E1             M29    act_naive    act_naive    <NA>
    ## 6            6          F1             M29    act_rTreg    act_rTreg    <NA>
    ## 7            7          G1             M30        naive        naive    <NA>
    ## 8            8          H1             M30        rTreg        rTreg    <NA>
    ## 9            9          A2             M30    act_naive    act_naive    <NA>
    ## 10          10          B2             M30    act_rTreg    act_rTreg    <NA>
    ## 11          11         H06 VICS-72098-18-B        birth        birth    <NA>
    ##     Array      Slide
    ## 1  R01C01 6264509100
    ## 2  R02C01 6264509100
    ## 3  R03C01 6264509100
    ## 4  R04C01 6264509100
    ## 5  R05C01 6264509100
    ## 6  R06C01 6264509100
    ## 7  R01C02 6264509100
    ## 8  R02C02 6264509100
    ## 9  R03C02 6264509100
    ## 10 R04C02 6264509100
    ## 11 R06C02 5975827018
    ##                                                                                                                        Basename
    ## 1  /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/6264509100/6264509100_R01C01
    ## 2  /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/6264509100/6264509100_R02C01
    ## 3  /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/6264509100/6264509100_R03C01
    ## 4  /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/6264509100/6264509100_R04C01
    ## 5  /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/6264509100/6264509100_R05C01
    ## 6  /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/6264509100/6264509100_R06C01
    ## 7  /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/6264509100/6264509100_R01C02
    ## 8  /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/6264509100/6264509100_R02C02
    ## 9  /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/6264509100/6264509100_R03C02
    ## 10 /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/6264509100/6264509100_R04C02
    ## 11 /Library/Frameworks/R.framework/Versions/4.2/Resources/library/methylationArrayAnalysis/extdata/5975827018/5975827018_R06C02

``` r
# only keep files from Slide ID 6264509100 (remove one extra sample from another study)
targets <- filter(targets, Slide == "6264509100")
```

Next, we read the IDAT files into R using the `read.metharray.exp`
function, also from the `minfi` package. This creates an `RGChannelSet`
object that contains all the raw intensity data, from both the red and
green colour channels, for each sample. Note that it takes in our
targets data frame as an argument - the function will look in the
`Basename` column for the filenames of the idat files to grab.

``` r
# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets = targets)
rgSet
```

    ## class: RGChannelSet 
    ## dim: 622399 10 
    ## metadata(0):
    ## assays(2): Green Red
    ## rownames(622399): 10600313 10600322 ... 74810490 74810492
    ## rowData names(0):
    ## colnames(10): 6264509100_R01C01 6264509100_R02C01 ... 6264509100_R03C02
    ##   6264509100_R04C02
    ## colData names(10): Sample_Name Sample_Well ... Basename filenames
    ## Annotation
    ##   array: IlluminaHumanMethylation450k
    ##   annotation: ilmn12.hg19

``` r
assays(rgSet)$Green[1:5, 1:5]
```

    ##          6264509100_R01C01 6264509100_R02C01 6264509100_R03C01
    ## 10600313               108               165               238
    ## 10600322              3855              5324              6653
    ## 10600328              1239              1263              1755
    ## 10600336               836              1094              1214
    ## 10600345               629               651               922
    ##          6264509100_R04C01 6264509100_R05C01
    ## 10600313               276               289
    ## 10600322              9404             10689
    ## 10600328              2821              3363
    ## 10600336              2242              2618
    ## 10600345              2769              2681

As you can see, this object looks a lot like an `ExpressionSet` object,
with two assays - one for Green and one for Red channel. Indeed, it
behaves in much the same way, but has several useful
methylation-specific methods we can perform on it.

At this stage, it can be useful to rename the samples with more
descriptive names.

``` r
# give the samples descriptive names
targets$ID <- paste(targets$Sample_Group, targets$Sample_Name, sep = ".")
sampleNames(rgSet) <- targets$ID
sampleNames(rgSet)
```

    ##  [1] "naive.1"      "rTreg.2"      "act_naive.3"  "naive.4"      "act_naive.5" 
    ##  [6] "act_rTreg.6"  "naive.7"      "rTreg.8"      "act_naive.9"  "act_rTreg.10"

## Normalization

To minimise the unwanted variation within and between samples, we need
to **normalize** our data. Many different types of normalization have
been developed for methylation arrays and it is beyond the scope of this
seminar to compare and contrast all of them. Several methods are built
into the `minfi` package. Although there is no single normalization
method that is universally considered best, the `preprocessFunnorm`
([Fortin et
al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0503-2))
function is most appropriate for datasets with **global methylation
differences** such as cancer/normal or vastly different tissue types,
whilst the `preprocessQuantile` function ([Touleimat and Tost
2012](https://www.futuremedicine.com/doi/full/10.2217/epi.12.21)) is
more suitable for datasets where you **don’t expect global differences**
(e.g. for example a single tissue). For further discussion on
appropriate choice of normalization, including data-driven tests for the
assumptions of quantile normalization, see [Hicks and Irizarry
(2015)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0679-0).

In this dataset, we are comparing different blood cell types, which we
expect to be globally relatively similar. Therefore, we will apply the
`preprocessQuantile` method. This function implements quantile
normalization on the methylated and unmethylated signal intensities
separately, and takes into account the different probe types
(e.g. background probes). We specify `sex = "M"` since all samples in
the dataset are male.

``` r
# normalize the data; this results in a GenomicRatioSet object
mSetNorm <- preprocessQuantile(rgSet, sex = "M")
```

    ## [preprocessQuantile] Mapping to genome.

    ## [preprocessQuantile] Fixing outliers.

    ## [preprocessQuantile] Quantile normalizing.

``` r
mSetNorm
```

    ## class: GenomicRatioSet 
    ## dim: 485512 10 
    ## metadata(0):
    ## assays(2): M CN
    ## rownames(485512): cg13869341 cg14008030 ... cg08265308 cg14273923
    ## rowData names(0):
    ## colnames(10): naive.1 rTreg.2 ... act_naive.9 act_rTreg.10
    ## colData names(10): Sample_Name Sample_Well ... Basename filenames
    ## Annotation
    ##   array: IlluminaHumanMethylation450k
    ##   annotation: ilmn12.hg19
    ## Preprocessing
    ##   Method: Raw (no normalization or bg correction)
    ##   minfi version: 1.44.0
    ##   Manifest version: 0.4.0

Note that after normalization, the data is housed in a `GenomicRatioSet`
object. This is a much more compact representation of the data as the
two colour channel information has been discarded and the M and U
intensity information has been converted to M-values, together with
associated genomic coordinates. We can pull out beta values with
`getBeta(mSetNorm)`.

Let’s plot the per-sample distributions of the raw versus normalized
beta values.

``` r
# vizualise what the data looks like before and after normalization
par(mfrow = c(1,2))
densityPlot(rgSet, sampGroups = targets$Sample_Group, main = "Raw", legend = FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetNorm), sampGroups = targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
```

![](sm07_methylation_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

As we can see, the distributions of beta values appear much more similar
after normalization.

Let’s also check out some PCA plots based on the top 1000 most variable
genes using `plotMDS` in limma.

``` r
# MDS plots to look at largest sources of variation
par(mfrow=c(1,2))
pal <- brewer.pal(8,"Dark2")

plotMDS(getM(mSetNorm), top = 1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetNorm), top = 1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)
```

![](sm07_methylation_files/figure-gfm/PCA-1.png)<!-- -->

Notice that PCs 1 and 2 largely separate samples by donor/individual.
Let’s look at higher dimensions.

``` r
# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,2))
plotMDS(getM(mSetNorm), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal, 
       cex=0.7, bg="white")

plotMDS(getM(mSetNorm), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
```

![](sm07_methylation_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

PC 3 seems to capture the cell type differences.

Let’s also look at a sample-sample correlation plot based on the top
1000 most variable probes.

``` r
top1000 <- which(rank(-rowVars(getM(mSetNorm))) <= 1000)
cc <- as.matrix(cor(getM(mSetNorm[top1000,]), 
                     use = "pairwise.complete.obs"), 
                     row.names = colnames(mSetNorm))

annot <- HeatmapAnnotation(df = data.frame(Individual = colData(mSetNorm)$Sample_Source,
                                           CellType = colData(mSetNorm)$Sample_Group),
                           col = list(Group = c("M28" =  pal[8], "M29" = pal[7], "M30" = pal[6]),
                                      CellType = c("naive" = pal[1], "rTreg" = pal[2], 
                                                   "act_naive" = pal[3], "act_rTreg" = pal[4])))

bcols <- colorRamp2(c(0, 1), c("red", "white"))

Heatmap(cc, name = "Corr", col = bcols,
        column_title = "Correlation of M values",
        cluster_rows = TRUE, cluster_columns = TRUE, 
        top_annotation = annot, 
        row_names_gp = gpar(fontsize = 8), 
        column_names_gp = gpar(fontsize = 8))
```

![](sm07_methylation_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Samples cluster strongly by donor/individual.

## Quality Control (Filtering)

First, let’s remove poor quality probes. We define poor quality probes
as those which have a total signal (M+U) that is similar to the
background signal. We evaluate this by calculating ‘detection p-values’
for each probe that test whether the total signal is greater than the
background signal. If there are any samples with low mean detection
p-values, they should be removed. We will also remove probes with
p-values larger than 0.01 in at least one sample (since small p-values
indicate reliable signal). We can calculate these detection p-values
with the `detectionP` function in `minfi`.

``` r
# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)
```

    ##                  naive.1       rTreg.2  act_naive.3     naive.4   act_naive.5
    ## cg00050873  0.000000e+00  0.000000e+00 0.000000e+00 0.00000e+00  0.000000e+00
    ## cg00212031  0.000000e+00  0.000000e+00 0.000000e+00 0.00000e+00  0.000000e+00
    ## cg00213748  2.139652e-88  4.213813e-31 1.181802e-12 1.29802e-47  8.255482e-15
    ## cg00214611  0.000000e+00  0.000000e+00 0.000000e+00 0.00000e+00  0.000000e+00
    ## cg00455876 1.400696e-234 9.349236e-111 4.272105e-90 0.00000e+00 3.347145e-268
    ## cg01707559  0.000000e+00  0.000000e+00 0.000000e+00 0.00000e+00  0.000000e+00
    ##              act_rTreg.6      naive.7       rTreg.8  act_naive.9  act_rTreg.10
    ## cg00050873  0.000000e+00  0.00000e+00  0.000000e+00 0.000000e+00  0.000000e+00
    ## cg00212031  0.000000e+00  0.00000e+00  0.000000e+00 0.000000e+00  0.000000e+00
    ## cg00213748  2.592206e-23  1.16160e-28  1.469801e-05 1.543654e-21  1.365951e-08
    ## cg00214611  0.000000e+00  0.00000e+00  0.000000e+00 0.000000e+00  0.000000e+00
    ## cg00455876 4.690740e-308 1.08647e-219 5.362780e-178 0.000000e+00 7.950724e-295
    ## cg01707559  0.000000e+00  0.00000e+00  0.000000e+00 0.000000e+00  0.000000e+00

``` r
# flag probes with pvals > 0.01 in any sample
keep <- rowSums(detP < 0.01) == ncol(mSetNorm) 
table(keep)
```

    ## keep
    ##  FALSE   TRUE 
    ##    977 484535

``` r
# remove poor quality probes from GenomicRatioSet
mSetNormFilt <- mSetNorm[keep,]
mSetNormFilt
```

    ## class: GenomicRatioSet 
    ## dim: 484535 10 
    ## metadata(0):
    ## assays(2): M CN
    ## rownames(484535): cg13869341 cg14008030 ... cg08265308 cg14273923
    ## rowData names(0):
    ## colnames(10): naive.1 rTreg.2 ... act_naive.9 act_rTreg.10
    ## colData names(10): Sample_Name Sample_Well ... Basename filenames
    ## Annotation
    ##   array: IlluminaHumanMethylation450k
    ##   annotation: ilmn12.hg19
    ## Preprocessing
    ##   Method: Raw (no normalization or bg correction)
    ##   minfi version: 1.44.0
    ##   Manifest version: 0.4.0

It is also advisable to consider removing **probes on sex chromosomes**
(if they are not of interest in your study), **probes that may interact
with known SNPs** (they may inflate inter-individual variation), and/or
**probes that map to multiple places in the genome**. Here, we’ll remove
probes that may interact with known SNPs with the convenient `minfi`
function `dropLociWithSnps` - for more info on other types of filtering,
see the [Bioconductor workflow on which this seminar is
based](https://bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html).

``` r
mSetNormFilt <- dropLociWithSnps(mSetNormFilt)
mSetNormFilt
```

    ## class: GenomicRatioSet 
    ## dim: 467023 10 
    ## metadata(0):
    ## assays(2): M CN
    ## rownames(467023): cg13869341 cg14008030 ... cg08265308 cg14273923
    ## rowData names(0):
    ## colnames(10): naive.1 rTreg.2 ... act_naive.9 act_rTreg.10
    ## colData names(10): Sample_Name Sample_Well ... Basename filenames
    ## Annotation
    ##   array: IlluminaHumanMethylation450k
    ##   annotation: ilmn12.hg19
    ## Preprocessing
    ##   Method: Raw (no normalization or bg correction)
    ##   minfi version: 1.44.0
    ##   Manifest version: 0.4.0

Let’s recheck the PCA plots.

``` r
# MDS plots to look at largest sources of variation
par(mfrow=c(1,2))
pal <- brewer.pal(8,"Dark2")

plotMDS(getM(mSetNormFilt), top = 1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetNormFilt), top = 1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("topright", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)
```

![](sm07_methylation_files/figure-gfm/PCA2-1.png)<!-- -->

We will also filter out probes that have shown to be cross-reactive
(i.e. probes that have been demonstrated to map to multiple places in
the genome). This list was originally published by [Chen et
al. (2013)](https://doi.org/10.4161/epi.23470).

``` r
# exclude cross reactive probes 
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetNormFilt) %in% xReactiveProbes$TargetID)
table(keep)
```

    ## keep
    ##  FALSE   TRUE 
    ##  27428 439595

``` r
mSetNormFilt <- mSetNormFilt[keep,] 
mSetNormFilt
```

    ## class: GenomicRatioSet 
    ## dim: 439595 10 
    ## metadata(0):
    ## assays(2): M CN
    ## rownames(439595): cg13869341 cg24669183 ... cg08265308 cg14273923
    ## rowData names(0):
    ## colnames(10): naive.1 rTreg.2 ... act_naive.9 act_rTreg.10
    ## colData names(10): Sample_Name Sample_Well ... Basename filenames
    ## Annotation
    ##   array: IlluminaHumanMethylation450k
    ##   annotation: ilmn12.hg19
    ## Preprocessing
    ##   Method: Raw (no normalization or bg correction)
    ##   minfi version: 1.44.0
    ##   Manifest version: 0.4.0

Great! After filtering poor quality probes and probes near SNPs, we see
that PC1 and PC2 are not as strongly associated with individual
variation (especially PC2).

Let’s recheck the sample-sample correlation plot.

``` r
top1000 <- which(rank(-rowVars(getM(mSetNormFilt))) <= 1000)
cc <- as.matrix(cor(getM(mSetNormFilt[top1000,]), 
                     use = "pairwise.complete.obs"), 
                     row.names = colnames(mSetNorm))

Heatmap(cc, name = "Corr", col = bcols,
        column_title = "Correlation of M values",
        cluster_rows = TRUE, cluster_columns = TRUE, 
        top_annotation = annot, 
        row_names_gp = gpar(fontsize = 8), 
        column_names_gp = gpar(fontsize = 8))
```

![](sm07_methylation_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Samples still cluster by individual/donor, but the correlation patterns
are not nearly as striking as before filtering.

## Probe-level differential methylation analysis

We are interested in pairwise comparisons between the four cell types,
taking into account individual to individual variation. We perform this
analysis on the matrix of M-values in `limma`, obtaining moderated
t-statistics and associated p-values for each CpG site (this part will
be very familiar to you!). We’ll use a design matrix with no intercept
term (recall the **cell means parameterization**).

A convenient way to set up the model when the user has many comparisons
of interest that they would like to test is to use a **contrast matrix**
in conjunction with the design matrix. A contrast matrix will take
linear combinations of the columns of the design matrix corresponding to
the comparisons of interest. Here we’ll construct a contrast matrix of
the 4 pairwise comparisons we are interested in.

``` r
# use the above to create a design matrix (cell means parameterization)
design <- model.matrix(~ 0 + Sample_Group + Sample_Source, 
                       data = targets)
# simplify column names of design matrix
colnames(design) <- gsub("Sample_Group|Sample_Source", "", colnames(design))
design
```

    ##    act_naive act_rTreg naive rTreg M29 M30
    ## 1          0         0     1     0   0   0
    ## 2          0         0     0     1   0   0
    ## 3          1         0     0     0   0   0
    ## 4          0         0     1     0   1   0
    ## 5          1         0     0     0   1   0
    ## 6          0         1     0     0   1   0
    ## 7          0         0     1     0   0   1
    ## 8          0         0     0     1   0   1
    ## 9          1         0     0     0   0   1
    ## 10         0         1     0     0   0   1
    ## attr(,"assign")
    ## [1] 1 1 1 1 2 2
    ## attr(,"contrasts")
    ## attr(,"contrasts")$Sample_Group
    ## [1] "contr.treatment"
    ## 
    ## attr(,"contrasts")$Sample_Source
    ## [1] "contr.treatment"

``` r
# fit the linear model 
fit <- lmFit(getM(mSetNormFilt), design)

# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(naive - rTreg,
                            naive - act_naive,
                            rTreg - act_rTreg,
                            act_naive - act_rTreg,
                            levels = design)
contMatrix
```

    ##            Contrasts
    ## Levels      naive - rTreg naive - act_naive rTreg - act_rTreg
    ##   act_naive             0                -1                 0
    ##   act_rTreg             0                 0                -1
    ##   naive                 1                 1                 0
    ##   rTreg                -1                 0                 1
    ##   M29                   0                 0                 0
    ##   M30                   0                 0                 0
    ##            Contrasts
    ## Levels      act_naive - act_rTreg
    ##   act_naive                     1
    ##   act_rTreg                    -1
    ##   naive                         0
    ##   rTreg                         0
    ##   M29                           0
    ##   M30                           0

``` r
# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))
```

    ##        naive - rTreg naive - act_naive rTreg - act_rTreg act_naive - act_rTreg
    ## Down            1607               397                 0                   557
    ## NotSig        436597            438977            439595                438128
    ## Up              1391               221                 0                   910

Now, we’ll show how to annotate the significant probes for genomic
context. We use the `IlluminaHumanMethylation450kmanifest` which
contains all of the annotation information for each of the CpG probes on
the 450k array. This information can be combined with the `topTable`
output. Here, we’ll demonstrate this for the first contrast (naive vs
rTreg).

``` r
# read in annotation and keep only particularly relevant columns (we've preselected some)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)[,c(1:4,18:19,24:30)]
head(ann450k)
```

    ## DataFrame with 6 rows and 13 columns
    ##                    chr       pos      strand        Name           Islands_Name
    ##            <character> <integer> <character> <character>            <character>
    ## cg00050873        chrY   9363356           -  cg00050873   chrY:9363680-9363943
    ## cg00212031        chrY  21239348           -  cg00212031 chrY:21238448-21240005
    ## cg00213748        chrY   8148233           -  cg00213748   chrY:8147877-8148210
    ## cg00214611        chrY  15815688           -  cg00214611 chrY:15815488-15815779
    ## cg00455876        chrY   9385539           -  cg00455876   chrY:9385471-9385777
    ## cg01707559        chrY   6778695           +  cg01707559   chrY:6778574-6780028
    ##            Relation_to_Island UCSC_RefGene_Name UCSC_RefGene_Accession
    ##                   <character>       <character>            <character>
    ## cg00050873            N_Shore    TSPY4;FAM197Y2 NM_001164471;NR_001553
    ## cg00212031             Island            TTTY14              NR_001543
    ## cg00213748            S_Shore                                         
    ## cg00214611             Island     TMSB4Y;TMSB4Y    NM_004202;NM_004202
    ## cg00455876             Island                                         
    ## cg01707559             Island TBL1Y;TBL1Y;TBL1Y NM_134259;NM_033284;..
    ##              UCSC_RefGene_Group     Phantom         DMR    Enhancer
    ##                     <character> <character> <character> <character>
    ## cg00050873         Body;TSS1500                                    
    ## cg00212031               TSS200                                    
    ## cg00213748                                                         
    ## cg00214611        1stExon;5'UTR                                    
    ## cg00455876                                                         
    ## cg01707559 TSS200;TSS200;TSS200                                    
    ##                     HMM_Island
    ##                    <character>
    ## cg00050873   Y:9973136-9976273
    ## cg00212031 Y:19697854-19699393
    ## cg00213748   Y:8207555-8208234
    ## cg00214611 Y:14324883-14325218
    ## cg00455876   Y:9993394-9995882
    ## cg01707559   Y:6838022-6839951

``` r
# overall Island annotation across all significant 
table(ann450k$Relation_to_Island) / nrow(ann450k)
```

    ## 
    ##     Island    N_Shelf    N_Shore    OpenSea    S_Shelf    S_Shore 
    ## 0.30947536 0.05117072 0.12949216 0.36260072 0.04593089 0.10133014

``` r
# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mSetNormFilt), ann450k$Name), ]
DMPs <- topTable(fit2, n = Inf, p.value = 0.05, coef = "naive - rTreg", genelist = ann450kSub)
head(DMPs)
```

    ##              chr       pos strand       Name            Islands_Name
    ## cg07499259  chr1  12188502      + cg07499259                        
    ## cg26992245  chr8  29848579      - cg26992245                        
    ## cg09747445 chr15  70387268      - cg09747445 chr15:70387929-70393206
    ## cg18808929  chr8  61825469      - cg18808929  chr8:61822358-61823028
    ## cg25015733  chr2  99342986      - cg25015733  chr2:99346882-99348177
    ## cg21179654  chr3 114057297      + cg21179654                        
    ##            Relation_to_Island                                UCSC_RefGene_Name
    ## cg07499259            OpenSea                                  TNFRSF8;TNFRSF8
    ## cg26992245            OpenSea                                                 
    ## cg09747445            N_Shore                                   TLE3;TLE3;TLE3
    ## cg18808929            S_Shelf                                                 
    ## cg25015733            N_Shelf                                           MGAT4A
    ## cg21179654            OpenSea ZBTB20;ZBTB20;ZBTB20;ZBTB20;ZBTB20;ZBTB20;ZBTB20
    ##                                                                             UCSC_RefGene_Accession
    ## cg07499259                                                                     NM_152942;NM_001243
    ## cg26992245                                                                                        
    ## cg09747445                                                        NM_001105192;NM_020908;NM_005078
    ## cg18808929                                                                                        
    ## cg25015733                                                                               NM_012214
    ## cg21179654 NM_001164343;NM_001164346;NM_001164345;NM_001164342;NM_001164344;NM_001164347;NM_015642
    ##                                   UCSC_RefGene_Group Phantom DMR Enhancer
    ## cg07499259                                5'UTR;Body                     
    ## cg26992245                                                           TRUE
    ## cg09747445                            Body;Body;Body                     
    ## cg18808929                                                           TRUE
    ## cg25015733                                     5'UTR                     
    ## cg21179654 3'UTR;3'UTR;3'UTR;3'UTR;3'UTR;3'UTR;3'UTR                     
    ##                     HMM_Island     logFC     AveExpr         t      P.Value
    ## cg07499259 1:12111023-12111225  3.654104  2.46652171  18.73274 7.455335e-08
    ## cg26992245                      4.450696 -0.09180715  18.32038 8.864345e-08
    ## cg09747445                     -3.337299 -0.25201484 -18.24876 9.138351e-08
    ## cg18808929                     -2.990263  0.77522878 -17.91096 1.056603e-07
    ## cg25015733                     -3.054336  0.83280190 -17.33212 1.363424e-07
    ## cg21179654                      2.859016  1.32460816  17.28740 1.391009e-07
    ##              adj.P.Val        B
    ## cg07499259 0.005180072 7.430392
    ## cg26992245 0.005180072 7.334745
    ## cg09747445 0.005180072 7.317706
    ## cg18808929 0.005180072 7.235565
    ## cg25015733 0.005180072 7.087729
    ## cg21179654 0.005180072 7.075921

``` r
# Look at Island annotation across all significant 
table(DMPs$Relation_to_Island) / nrow(DMPs)
```

    ## 
    ##     Island    N_Shelf    N_Shore    OpenSea    S_Shelf    S_Shore 
    ## 0.05537025 0.06937959 0.14142762 0.57204803 0.06070714 0.10106738

Now we’ll plot the top 5 differentially methylated probes.

``` r
top5 <- getBeta(mSetNormFilt)[rownames(DMPs[1:5,]),] %>%
  data.frame() %>%
  rownames_to_column(var = "probe") %>%
  pivot_longer(cols = -probe, names_to = c("Sample"), values_to = "Beta") %>%
  left_join(colData(mSetNormFilt) %>% data.frame() %>% rownames_to_column(var = "Sample"), 
            by = "Sample") 

top5 %>%
   ggplot(aes(x = Sample_Group, y = Beta, color = Sample_Group)) + 
   geom_point(position = position_jitter(width = 0.05), na.rm = T) + 
   facet_grid(. ~ probe) + 
   ggtitle("Probe beta values within top 5 DM (naive-rTreg) probes") + 
   xlab("Sample Type") + 
   ylab("Beta") +
   theme(axis.text.x = element_text(angle = 90))
```

![](sm07_methylation_files/figure-gfm/stripplot-1.png)<!-- -->

All of these show a clear difference between naive and rTreg.

Finally, plot location of differentially methylated probes along each
chromosome.

``` r
# fetch chromosome length info
chrlen = getChromInfoFromUCSC("hg19") %>%
  rename(chr = chrom) %>%
  filter(chr %in% c("chrX", "chrY", paste("chr",1:22,sep=""))) %>%
  mutate(chr = factor(chr, levels = c(paste("chr",1:22,sep=""), "chrX", "chrY")))

DMPs %>%
  ggplot() +
    geom_linerange(aes(x = chr, ymin = 0, ymax = size), data = chrlen, alpha = 0.5) + 
    geom_point(aes(x = chr, y = pos), alpha = 0.5,
              position = position_jitter(width = 0.03), na.rm = TRUE) + 
    ggtitle("DMP (naive-rTreg) positions along the genome") + 
    ylab("Position of DMPs") +
    xlab("chr") +
    coord_flip() 
```

![](sm07_methylation_files/figure-gfm/coord-1.png)<!-- -->

## Region level differential methylation analysis

Although performing a probe-wise analysis is useful and informative,
sometimes we are interested in knowing whether several neighboring CpGs
are concordantly differentially methylated. In other words, we want to
identify **differentially methylated regions** (DMRs). There are several
Bioconductor packages that have functions for identifying differentially
methylated regions from 450k data, including `bumphunter` and `DMRcate`.
Here, we will demonstrate an analysis using the `DMRcate`. As it is
based on limma, we can directly use the design and contrast matrices we
previously defined.

First, our matrix of M-values is annotated with the relevant information
about the probes such as their genomic position, gene annotation, etc.
By default, this is done using the `ilmn12.hg19` annotation, but this
can be modified. The limma pipeline is then used for differential
methylation analysis to calculate moderated t-statistics.

``` r
myAnnotation <- cpg.annotate(object = getM(mSetNormFilt), datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "naive - rTreg", arraytype = "450K")
```

    ## Your contrast returned 2998 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

``` r
myAnnotation
```

    ## CpGannotated object describing 439595 CpG sites, with independent
    ## CpG threshold indexed at fdr=0.05 and 2998 significant CpG sites.

Now that we have the relevant statistics for the individual CpGs, we can
then use the `dmrcate` function to combine them to identify DMRs. The
parameter `lambda` here represents a smoothing bandwidth (in basepairs).
Here we use the default value of 1000. `C` is a related parameter that
is recommended to be set to 2.

We then need to run `extractRanges` to pull out a `GRanges` object with
our results.

``` r
# run DMRcate
DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2)
```

    ## Fitting chr1...

    ## Fitting chr2...

    ## Fitting chr3...

    ## Fitting chr4...

    ## Fitting chr5...

    ## Fitting chr6...

    ## Fitting chr7...

    ## Fitting chr8...

    ## Fitting chr9...

    ## Fitting chr10...

    ## Fitting chr11...

    ## Fitting chr12...

    ## Fitting chr13...

    ## Fitting chr14...

    ## Fitting chr15...

    ## Fitting chr16...

    ## Fitting chr17...

    ## Fitting chr18...

    ## Fitting chr19...

    ## Fitting chr20...

    ## Fitting chr21...

    ## Fitting chr22...

    ## Fitting chrX...

    ## Fitting chrY...

    ## Demarcating regions...

    ## Done!

``` r
DMRs
```

    ## DMResults object with 538 DMRs.
    ## Use extractRanges() to produce a GRanges object of these.

``` r
# pull out results
results.ranges <- extractRanges(DMRs, genome = "hg19")
```

    ## snapshotDate(): 2022-10-31

    ## see ?DMRcatedata and browseVignettes('DMRcatedata') for documentation

    ## loading from cache

``` r
results.ranges
```

    ## GRanges object with 538 ranges and 8 metadata columns:
    ##         seqnames              ranges strand |   no.cpgs min_smoothed_fdr
    ##            <Rle>           <IRanges>  <Rle> | <integer>        <numeric>
    ##     [1]    chr17   57915665-57918682      * |        12      5.90700e-91
    ##     [2]     chr3 114012316-114012912      * |         5     1.97313e-180
    ##     [3]    chr18   21452730-21453131      * |         7     5.67432e-115
    ##     [4]    chr17   74639731-74640078      * |         6      9.20242e-90
    ##     [5]     chrX   49121205-49122718      * |         6      7.56676e-84
    ##     ...      ...                 ...    ... .       ...              ...
    ##   [534]     chr2   43454761-43455103      * |        14      1.25073e-25
    ##   [535]     chr6   31832650-31833452      * |        18      2.50956e-28
    ##   [536]     chr6 144385771-144387124      * |        22      2.78436e-60
    ##   [537]     chrX   43741310-43742501      * |        10      2.74791e-60
    ##   [538]     chr2   25141532-25142229      * |         8      3.94046e-25
    ##            Stouffer      HMFDR      Fisher   maxdiff   meandiff
    ##           <numeric>  <numeric>   <numeric> <numeric>  <numeric>
    ##     [1] 7.38597e-10 0.02385182 7.27170e-08  0.398286   0.313161
    ##     [2] 1.63569e-07 0.00722324 1.50179e-06  0.543428   0.425162
    ##     [3] 8.12889e-07 0.01257976 1.97231e-06 -0.386747  -0.254609
    ##     [4] 1.61913e-07 0.01419742 2.74827e-06 -0.252864  -0.195190
    ##     [5] 3.15833e-07 0.01186578 3.82468e-06  0.452909   0.300624
    ##     ...         ...        ...         ...       ...        ...
    ##   [534]    0.967997  0.1544075    0.623920 -0.218836 -0.0427390
    ##   [535]    0.887299  0.2322842    0.650632  0.153367  0.0490080
    ##   [536]    0.996676  0.0806609    0.694295  0.325422  0.0449451
    ##   [537]    0.924207  0.0486641    0.724933  0.413832  0.0234547
    ##   [538]    0.992465  0.0576007    0.750823  0.282058  0.0314244
    ##         overlapping.genes
    ##               <character>
    ##     [1]       VMP1, MIR21
    ##     [2]             TIGIT
    ##     [3]             LAMA3
    ##     [4]        ST6GALNAC1
    ##     [5]             FOXP3
    ##     ...               ...
    ##   [534]             THADA
    ##   [535]           SLC44A4
    ##   [536]              <NA>
    ##   [537]              MAOB
    ##   [538]             ADCY3
    ##   -------
    ##   seqinfo: 23 sequences from an unspecified genome; no seqlengths

Let’s visualize the results for DMR number 8 (which I selected because
it has a relatively high number of CpGs (14), and a relatively high mean
beta difference). The `DMR.plot` function from `DMRcate` plots
sample-level heatmaps and smoothed group means (only useful for larger
numbers of CpGs) of methylation values in the context of genomic
location.

``` r
# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 8, CpGs = getBeta(mSetNormFilt), phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")
```

    ## snapshotDate(): 2022-10-31

    ## see ?DMRcatedata and browseVignettes('DMRcatedata') for documentation

    ## loading from cache

    ## Warning in updateObjectFromSlots(object, ..., verbose = verbose): dropping
    ## slot(s) 'columns' from object = 'GeneRegionTrack'

![](sm07_methylation_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## Exercise

There are no graded deliverables for this seminar, but as an exercise,
repeat the analysis in the previous sections (plot top 5 differentially
methylated probes, plot all differentially methylated probes across the
genome, and obtain DMRs) for one of the other contrasts. In addition,
add the name(s) of nearest genes to the plot of the top 5 differentially
methylated probes.

## Interpretation and functional enrichment analysis

The next step in the analysis pipeline would be to interpret our result
by associating these differentially methylated regions (DMRs) with
biological features. DNA methylation in different genomic regions,
i.e. promoter regions, enhancers, gene body etc., has different impact
on gene transcription. A common practice in DMR analysis is to separate
DMRs associated with different types of genomic regions and study them
separately. Our key interest is to see how these DMRs affect gene
expression level and consequently biological function. Gene set
enrichment analyses will be explored later in the course (Seminar 10).

## Methylation sequencing

We’ve only scratched the surface of DNA methylation analysis. This
seminar focused entirely on analysis of DNA methylation microarrays (in
particular the 450K array). As sequencing technologies are becoming
popular, techniques like Whole Genome Bisulfite Sequencing (WGBS) and
Methylated DNA Immunoprecipitation followed by sequencing (MEDIP-seq)
are becoming more widely used. These bring with them new analysis
challenges. Refer to [lecture
11](https://stat540-ubc.github.io/lectures/lect11-epigenetics/lect11-epigenetics.html)
for an overview.

## Session info

``` r
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] circlize_0.4.15                                    
    ##  [2] ComplexHeatmap_2.14.0                              
    ##  [3] forcats_0.5.2                                      
    ##  [4] dplyr_1.0.10                                       
    ##  [5] purrr_1.0.0                                        
    ##  [6] readr_2.1.3                                        
    ##  [7] tidyr_1.2.1                                        
    ##  [8] tibble_3.1.8                                       
    ##  [9] ggplot2_3.4.0                                      
    ## [10] tidyverse_1.3.2                                    
    ## [11] DMRcatedata_2.16.0                                 
    ## [12] ExperimentHub_2.6.0                                
    ## [13] AnnotationHub_3.6.0                                
    ## [14] BiocFileCache_2.6.0                                
    ## [15] dbplyr_2.2.1                                       
    ## [16] methylationArrayAnalysis_1.22.0                    
    ## [17] FlowSorted.Blood.450k_1.36.0                       
    ## [18] stringr_1.5.0                                      
    ## [19] DMRcate_2.12.0                                     
    ## [20] Gviz_1.42.0                                        
    ## [21] minfiData_0.44.0                                   
    ## [22] missMethyl_1.32.0                                  
    ## [23] IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
    ## [24] RColorBrewer_1.1-3                                 
    ## [25] IlluminaHumanMethylation450kmanifest_0.4.0         
    ## [26] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1 
    ## [27] minfi_1.44.0                                       
    ## [28] bumphunter_1.40.0                                  
    ## [29] locfit_1.5-9.7                                     
    ## [30] iterators_1.0.14                                   
    ## [31] foreach_1.5.2                                      
    ## [32] Biostrings_2.66.0                                  
    ## [33] XVector_0.38.0                                     
    ## [34] SummarizedExperiment_1.28.0                        
    ## [35] Biobase_2.58.0                                     
    ## [36] MatrixGenerics_1.10.0                              
    ## [37] matrixStats_0.63.0                                 
    ## [38] GenomicRanges_1.50.2                               
    ## [39] GenomeInfoDb_1.34.7                                
    ## [40] IRanges_2.32.0                                     
    ## [41] S4Vectors_0.36.1                                   
    ## [42] BiocGenerics_0.44.0                                
    ## [43] limma_3.54.0                                       
    ## [44] BiocStyle_2.26.0                                   
    ## [45] rmarkdown_2.19                                     
    ## [46] knitr_1.41                                         
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rappdirs_0.3.3                rtracklayer_1.58.0           
    ##   [3] R.methodsS3_1.8.2             bit64_4.0.5                  
    ##   [5] DelayedArray_0.24.0           R.utils_2.12.2               
    ##   [7] data.table_1.14.6             rpart_4.1.19                 
    ##   [9] doParallel_1.0.17             KEGGREST_1.38.0              
    ##  [11] RCurl_1.98-1.9                GEOquery_2.66.0              
    ##  [13] AnnotationFilter_1.22.0       generics_0.1.3               
    ##  [15] GenomicFeatures_1.50.4        preprocessCore_1.60.2        
    ##  [17] RSQLite_2.2.20                bit_4.0.5                    
    ##  [19] tzdb_0.3.0                    xml2_1.3.3                   
    ##  [21] lubridate_1.9.0               httpuv_1.6.7                 
    ##  [23] assertthat_0.2.1              gargle_1.2.1                 
    ##  [25] xfun_0.36                     hms_1.1.2                    
    ##  [27] evaluate_0.19                 promises_1.2.0.1             
    ##  [29] fansi_1.0.3                   restfulr_0.0.15              
    ##  [31] scrime_1.3.5                  progress_1.2.2               
    ##  [33] readxl_1.4.1                  DBI_1.1.3                    
    ##  [35] htmlwidgets_1.6.0             reshape_0.8.9                
    ##  [37] googledrive_2.0.0             ellipsis_0.3.2               
    ##  [39] backports_1.4.1               permute_0.9-7                
    ##  [41] annotate_1.76.0               biomaRt_2.54.0               
    ##  [43] deldir_1.0-6                  sparseMatrixStats_1.10.0     
    ##  [45] vctrs_0.5.1                   ensembldb_2.22.0             
    ##  [47] cachem_1.0.6                  withr_2.5.0                  
    ##  [49] BSgenome_1.66.2               checkmate_2.1.0              
    ##  [51] GenomicAlignments_1.34.0      prettyunits_1.1.1            
    ##  [53] mclust_6.0.0                  cluster_2.1.4                
    ##  [55] lazyeval_0.2.2                crayon_1.5.2                 
    ##  [57] genefilter_1.80.3             labeling_0.4.2               
    ##  [59] edgeR_3.40.2                  pkgconfig_2.0.3              
    ##  [61] nlme_3.1-160                  ProtGenerics_1.30.0          
    ##  [63] nnet_7.3-18                   rlang_1.0.6                  
    ##  [65] lifecycle_1.0.3               filelock_1.0.2               
    ##  [67] modelr_0.1.10                 dichromat_2.0-0.1            
    ##  [69] cellranger_1.1.0              rngtools_1.5.2               
    ##  [71] base64_2.0.1                  Matrix_1.5-1                 
    ##  [73] Rhdf5lib_1.20.0               reprex_2.0.2                 
    ##  [75] base64enc_0.1-3               GlobalOptions_0.1.2          
    ##  [77] googlesheets4_1.0.1           png_0.1-8                    
    ##  [79] rjson_0.2.21                  bitops_1.0-7                 
    ##  [81] R.oo_1.25.0                   rhdf5filters_1.10.0          
    ##  [83] blob_1.2.3                    DelayedMatrixStats_1.20.0    
    ##  [85] shape_1.4.6                   doRNG_1.8.6                  
    ##  [87] nor1mix_1.3-0                 jpeg_0.1-10                  
    ##  [89] scales_1.2.1                  memoise_2.0.1                
    ##  [91] magrittr_2.0.3                plyr_1.8.8                   
    ##  [93] zlibbioc_1.44.0               compiler_4.2.2               
    ##  [95] BiocIO_1.8.0                  illuminaio_0.40.0            
    ##  [97] clue_0.3-64                   Rsamtools_2.14.0             
    ##  [99] cli_3.5.0                     DSS_2.46.0                   
    ## [101] htmlTable_2.4.1               Formula_1.2-4                
    ## [103] MASS_7.3-58.1                 tidyselect_1.2.0             
    ## [105] stringi_1.7.8                 highr_0.10                   
    ## [107] yaml_2.3.6                    askpass_1.1                  
    ## [109] latticeExtra_0.6-30           VariantAnnotation_1.44.0     
    ## [111] tools_4.2.2                   timechange_0.1.1             
    ## [113] rstudioapi_0.14               foreign_0.8-83               
    ## [115] bsseq_1.34.0                  gridExtra_2.3                
    ## [117] farver_2.1.1                  digest_0.6.31                
    ## [119] BiocManager_1.30.19           shiny_1.7.4                  
    ## [121] quadprog_1.5-8                Rcpp_1.0.9                   
    ## [123] siggenes_1.72.0               broom_1.0.2                  
    ## [125] BiocVersion_3.16.0            later_1.3.0                  
    ## [127] org.Hs.eg.db_3.16.0           httr_1.4.4                   
    ## [129] AnnotationDbi_1.60.0          biovizBase_1.46.0            
    ## [131] colorspace_2.0-3              rvest_1.0.3                  
    ## [133] XML_3.99-0.13                 fs_1.5.2                     
    ## [135] splines_4.2.2                 statmod_1.5.0                
    ## [137] multtest_2.54.0               xtable_1.8-4                 
    ## [139] jsonlite_1.8.4                R6_2.5.1                     
    ## [141] Hmisc_4.8-0                   pillar_1.8.1                 
    ## [143] htmltools_0.5.4               mime_0.12                    
    ## [145] glue_1.6.2                    fastmap_1.1.0                
    ## [147] BiocParallel_1.32.5           interactiveDisplayBase_1.36.0
    ## [149] beanplot_1.3.1                codetools_0.2-18             
    ## [151] utf8_1.2.2                    lattice_0.20-45              
    ## [153] curl_4.3.3                    gtools_3.9.4                 
    ## [155] openssl_2.0.5                 interp_1.1-3                 
    ## [157] survival_3.4-0                munsell_0.5.0                
    ## [159] GetoptLong_1.0.5              rhdf5_2.42.0                 
    ## [161] GenomeInfoDbData_1.2.9        HDF5Array_1.26.0             
    ## [163] haven_2.5.1                   gtable_0.3.1
