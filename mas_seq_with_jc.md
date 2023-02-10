scRNA-seq Analyis Using MAS-Seq with CRISPRClean
================
2022-12-02

-   [1 Plotting median UMIs/cell and
    Genes/cell](#1-plotting-median-umiscell-and-genescell)
-   [2 Setup the seurat object](#2-setup-the-seurat-object)
-   [3 Data pre-processing and QC
    workflow](#3-data-pre-processing-and-qc-workflow)
    -   [3.1 Filtering low-quality
        cells](#31-filtering-low-quality-cells)
        -   [3.1.1 Scatter plot of genes vs
            counts](#311-scatter-plot-of-genes-vs-counts)
        -   [3.1.2 Remove cells with too few genes detected and cells
            that were sequenced too
            high](#312-remove-cells-with-too-few-genes-detected-and-cells-that-were-sequenced-too-high)
        -   [3.1.3 Filtering low-complexity
            cells](#313-filtering-low-complexity-cells)
    -   [3.2 Doublet filtering](#32-doublet-filtering)
        -   [3.2.1 clustering using SCTransform
            workflow](#321-clustering-using-sctransform-workflow)
        -   [3.2.2 scDblFinder doublet
            calling](#322-scdblfinder-doublet-calling)
    -   [3.3 Cell cycle evaluation](#33-cell-cycle-evaluation)
-   [4 Cell clustering using SCTransform
    workflow](#4-cell-clustering-using-sctransform-workflow)
-   [5 Harmony data projection](#5-harmony-data-projection)

<style type="text/css">
  body{
  font-size: 14pt;
}
</style>

# 1 Plotting median UMIs/cell and Genes/cell

We first want to assess some basic metrics such as UMIs/cell and
Genes/cell to assess the gain in transcriptomic read information of
untargeted genes.

We created a series of functions to assist with the workflow in
assessing benefits of depletion.

You will use a directory containing the barcodes.tsv, genes.tsv, and
matrix.mtx file as output from PacBio’s Pigeon toolkit.

You will also need the gene target list provided by Jumpcode Genomics.
Because we are targeting \~350 genes for removal, we experimentally
remove UMIs associated with the targeted transcripts. We want to
evaluate how reads are being redistributed for genes that have not been
targeted with CRISPR. Thus, we should see a boost in both UMIs/cell and
Genes/cell when considering all non-targeted transcripts.

``` r
#We are first going to create a function to evaluate both umis/cell and genes/cell for each sample
umis_genes_per_cell <- function(matrix, gene_list, sample) {
    x<-c("Seurat", "dplyr", "patchwork")
    lapply(x, require, character.only = TRUE)
    if(dir.exists(paths = matrix)) {
        mtx <- Read10X(matrix)
    }
    if(file.exists(gene_list)) {
      targets <- read.delim("~/targets.txt", header = F)
      targets <- targets$V1
    }
    #remove the targeted genes from the control matrix file
    matrix.removed.genes <- mtx[!rownames(mtx) %in% targets, ]
    
    #here we only consider genes that are expressed in at least 5 cells
    keep <- rowSums(matrix.removed.genes>0)>=5
    matrix.removed.genes <- matrix.removed.genes[keep,]
    
    #create a dataframe to plot umis/cell
    df <- as.data.frame(colSums(matrix.removed.genes))
    df <- df %>% dplyr::rename(c("umi_counts"=1)) %>% mutate(gene_counts = colSums(matrix.removed.genes>0), condition = sample)
}
```

``` r
#matrix = PATH to directory containing barcodes.tsv, genes.tsv, and matrix.mtx from PacBio pipeline generating Seurat inputs
#gene list = gene list provided by Jumpcode containing gene targets
#sample = whichever name you would like to name the conditions to compare (i.e., MAS-Seq vs Jumpcode Depletion)

#control
control_dirs <- c("~/R/mas_seq/control3_genes_seurat/mas_untreated_genes_seurat/genes_seurat/")

control <- umis_genes_per_cell(matrix = control_dirs, gene_list = "~/targets.txt", sample = "MAS-Seq")
```

    ## Loading required package: Seurat

    ## Attaching SeuratObject

    ## Attaching sp

    ## Loading required package: dplyr

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Loading required package: patchwork

``` r
#depleted
depleted_dirs <- c("~/R/mas_seq/depleted3_genes_seurat/mas_treated_genes_seurat/genes_seurat/")

depleted <- umis_genes_per_cell(matrix = depleted_dirs, gene_list = "~/targets.txt", sample = "MAS-Seq + JC")
```

``` r
#plot histogram of UMIs and genes/cell
depletion_benefit <- function(control,depleted) {
    require("ggplot2")
    df.m <- rbind(control, depleted)
    median <- df.m %>% dplyr::group_by(condition) %>% dplyr::summarize(median_umi=median(umi_counts), median_genes=median(gene_counts))
    
    p1 <- ggplot(df.m, aes(x =umi_counts , fill=condition, ..scaled..)) + geom_density(alpha=0.2) + 
        xlab("nUMIs/cell") + ylab("Scaled Density") + xlim(0,median(df.m$umi_counts) + 8*mad(df.m$umi_counts)) +
        geom_vline(data = median, aes(xintercept = median_umi, color=condition), linetype='dashed') + geom_text(data = median, aes(x = round(median(df.m$umi_counts) + 1*mad(df.m$umi_counts)), y = c(1,0.9), label=paste0(condition, ":", median_umi), color=condition), size=4, hjust=0, fontface=2, show.legend = F)
        
      #this will generate plot comparing genes/cell
    p2 <- ggplot(df.m, aes(x =gene_counts , fill=condition, ..scaled..)) + geom_density(alpha=0.2) + 
        xlab("nGenes/cell") + ylab(NULL) +
        geom_vline(data = median, aes(xintercept = median_genes, color=condition), linetype='dashed') + xlim(0,median(df.m$gene_counts) + 8*mad(df.m$gene_counts)) + geom_text(data = median, aes(x = round(median(df.m$gene_counts) + 1*mad(df.m$gene_counts)), y = c(1,0.9), label=paste0(condition, ":", median_genes), color=condition), size=4, hjust=0, fontface=2, show.legend = F)
    p1 + p2 + plot_layout(guides = "collect", widths = c(2,2)) & theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(), axis.title.y = element_text(size = 14), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white'))
}
```

``` r
#generate the plots
depletion_benefit(control = control, depleted = depleted)
```

    ## Loading required package: ggplot2

    ## Warning: Removed 195 rows containing non-finite values (stat_density).

    ## Warning: Removed 35 rows containing non-finite values (stat_density).

![](mas_seq_with_jc_files/figure-gfm/plot%20the%20function-1.png)<!-- -->

Here we can see that we’re boosting the reads of the untargeted
transcripts of the 10x libraries.

# 2 Setup the seurat object

We are going to use the Seurat toolkit to perform all subsequent
downstream analysis for this tutorial.

We start by reading the data using the Read10x() function which reads
the directory containing the Seurat inputs generated by Pigeon: PacBio’s
Transcript Toolkit and generates a count matrix. We next use the count
matrix to create a seurat object.

``` r
library(Seurat)
#control
pbmc_control.mtx <- Read10X("~/R/mas_seq/control3_genes_seurat/mas_untreated_genes_seurat/genes_seurat/")
#depleted
pbmc_depleted.mtx <- Read10X("~/R/mas_seq/depleted3_genes_seurat/mas_treated_genes_seurat/genes_seurat/")
```

We are only going to include genes that have at least 3 total isoform
counts. Later, we will use data-driven techniques to determine
appropriate thresholds for filtering.

``` r
#control
pbmc_control.so <- CreateSeuratObject(pbmc_control.mtx, min.cells = 3, project = "MAS-Seq")
#depleted
pbmc_depleted.so <- CreateSeuratObject(pbmc_depleted.mtx, min.cells = 3, project = "JC")

#create a list containing the control and depleted conditions
list.so <- list(pbmc_control.so,pbmc_depleted.so)
#rename the indices of the list
names(list.so) <- c("mas","jc")
#remove the separate files to conserve memory
rm(pbmc_control.so,pbmc_depleted.so)
```

# 3 Data pre-processing and QC workflow

The steps below demonstrate an example of data pre-processing and QC for
scRNA-seq data generated with MAS-Seq sequencing combined with Jumpcode
depletion.

## 3.1 Filtering low-quality cells

We are going to use a data-driven technique to evaluate the number of
unique genes and UMIs detected in each cell.

1.  Empty droplets or low-quality cells typically express few genes and
    have low library complexity.
2.  Doublets exhibit higher UMI and gene counts which affects downstream
    interpretation.
3.  The total number of genes detected corresponds linearly with unique
    genes. Over-sequenced cells will have a disruption in this
    linearity.

### 3.1.1 Scatter plot of genes vs counts

``` r
library(ggplot2)
library(ggExtra)
```

    ## Warning: package 'ggExtra' was built under R version 4.2.2

``` r
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:patchwork':
    ## 
    ##     align_plots

``` r
#control
p1 <- ggplot(list.so$mas@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + geom_smooth(method="lm") + ggtitle("MAS-Seq")
p1 <- ggMarginal(p1, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

``` r
p2 <- ggplot(list.so$mas@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
p2 <- ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

``` r
plot_grid(plotlist = list(p1,p2), ncol=2, align='h', rel_widths = c(1, 1))
```

    ## Warning: Graphs cannot be horizontally aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

![](mas_seq_with_jc_files/figure-gfm/plotting%20relationship%20between%20detected%20genes%20and%20counts%20control-1.png)<!-- -->

``` r
#depleted 
p1 <- ggplot(list.so$jc@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + geom_smooth(method="lm") + ggtitle("MAS-Seq + JC")
p1 <- ggMarginal(p1, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

``` r
p2 <- ggplot(list.so$jc@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
p2 <- ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

``` r
plot_grid(plotlist = list(p1,p2), ncol=2, align='h', rel_widths = c(1, 1))
```

    ## Warning: Graphs cannot be horizontally aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

![](mas_seq_with_jc_files/figure-gfm/plotting%20relationship%20between%20detected%20genes%20and%20counts%20depleted-1.png)<!-- -->

### 3.1.2 Remove cells with too few genes detected and cells that were sequenced too high

We are going to work with log10 transformed data because it better
recapitulates the linear relationship as an effect of sequencing.

Using a data-driven technique, we are going use the median absolute
deviation (mad) as cutoff points. Typically, this can range from 3-5 mad
from the median as a cutoff. Make sure to test different options with
your data to determine an appropriate cutoff based on distributions.

For this dataset, we are going to apply several thresholds: 1. We are
going to set the maximum UMI threshold to 4 mads because we do not have
many cells over-sequenced (or approaching sequence saturation) as
determined by the upper limit of our distribution and would like to keep
as many cells for downstream analysis as possible 2. We are going to set
a minimum gene threshold to 3 mads of the gene counts/cell as evident by
the long left-tail of the distribution to remove these low-quality cells

``` r
# Gene/UMI scatter plot for control

p1 <- ggplot(list.so$mas@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm") + geom_hline(aes(yintercept = median(log10(list.so$mas$nFeature_RNA)) - 3*mad(log10(list.so$mas$nFeature_RNA))), colour = "green", linetype = 2) +
  geom_vline(aes(xintercept = median(log10(list.so$mas$nCount_RNA)) + 4*mad(log10(list.so$mas$nCount_RNA))), colour = "red", linetype = 2)

ggMarginal(p1, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](mas_seq_with_jc_files/figure-gfm/low%20quality%20cell%20removal%20control-1.png)<!-- -->

For replication purposes, we are going to apply the same thresholds to
the depleted samples. Keep in mind that with depletion, by increasing
counts of untargeted transcripts, we can achieve a higher sequence
saturation at the same sequence depth. As a result, these thresholds can
be changed depending on the experiment.

``` r
# Gene/UMI scatter plot for depleted

p1 <- ggplot(list.so$jc@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm") + geom_hline(aes(yintercept = median(log10(list.so$jc$nFeature_RNA)) - 3*mad(log10(list.so$jc$nFeature_RNA))), colour = "green", linetype = 2) +
  geom_vline(aes(xintercept = median(log10(list.so$jc$nCount_RNA)) + 4*mad(log10(list.so$jc$nCount_RNA))), colour = "red", linetype = 2)

ggMarginal(p1, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](mas_seq_with_jc_files/figure-gfm/low%20quality%20cell%20removal%20depleted-1.png)<!-- -->

``` r
#filter the cells from the control

min.gene.thresh <- median(log10(list.so$mas$nFeature_RNA)) - 3*mad(log10(list.so$mas$nFeature_RNA))
max.umi.thresh <- median(log10(list.so$mas$nCount_RNA)) + 3*mad(log10(list.so$mas$nCount_RNA))

cells.keep <- rownames(list.so$mas@meta.data %>% filter(log10(nFeature_RNA) > min.gene.thresh) %>% filter(log10(nCount_RNA) < max.umi.thresh))

list.so$mas <- subset(list.so$mas, cells = cells.keep)
```

``` r
#filter the cells from the depleted

min.gene.thresh.jc <- median(log10(list.so$jc$nFeature_RNA)) - 3*mad(log10(list.so$jc$nFeature_RNA))
max.umi.thresh.jc <- median(log10(list.so$jc$nCount_RNA)) + 3*mad(log10(list.so$jc$nCount_RNA))

cells.keep.jc <- rownames(list.so$jc@meta.data %>% filter(log10(nFeature_RNA) > min.gene.thresh.jc) %>% filter(log10(nCount_RNA) < max.umi.thresh.jc))

list.so$jc <- subset(list.so$jc, cells = cells.keep.jc)
```

### 3.1.3 Filtering low-complexity cells

For this workflow, we are going to perform a QC filtering based on
complexity or the ratio of Genes/UMI. Typically, dead or decaying cells
will express very few genes contributing to their total UMI count. As
previously mentioned, the number of genes/cell should scale with an
increase in sequencing depth. As a result, cells with a library
complexity outside the expected ratio are deemed lower-quality and
should be removed for QC.

Due to the linear relationship between the log10(UMI) and log10(gene
counts), we can use a linear model to calculate the residuals in
relation to the regression line. For this example, we are going to
exclude cells that have residuals with \>20% variance below the linear
regression to exclude low complexity cells. Keep in mind that this
threshold can change depending on the cell-type and experiment.

``` r
#control
lm.model = lm(data = list.so$mas@meta.data, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
list.so$mas@meta.data$residuals <- residuals(lm.model)
list.so$mas@meta.data <- list.so$mas@meta.data %>% mutate(complexity = ifelse(test = abs(list.so$mas@meta.data$residuals) <= 0.20, yes = "high" , no = "low"))
p2 <- ggplot(list.so$mas@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point(aes(color = complexity)) + geom_abline(intercept = lm.model$coefficients[1] - 0.20 , slope = lm.model$coefficients[2], color="orange", linetype=2) + geom_smooth(method="lm")
ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](mas_seq_with_jc_files/figure-gfm/plotting%20low%20complexity%20cells%20control-1.png)<!-- -->

``` r
#depleted
lm.model = lm(data = list.so$jc@meta.data, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
list.so$jc@meta.data$residuals <- residuals(lm.model)
list.so$jc@meta.data <- list.so$jc@meta.data %>% mutate(complexity = ifelse(test = abs(list.so$jc@meta.data$residuals) <= 0.20, yes = "high" , no = "low"))
p2 <- ggplot(list.so$jc@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point(aes(color = complexity)) + geom_abline(intercept = lm.model$coefficients[1] - 0.20 , slope = lm.model$coefficients[2], color="orange", linetype=2) + geom_smooth(method="lm")
ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](mas_seq_with_jc_files/figure-gfm/plotting%20low%20complexity%20cells%20depleted-1.png)<!-- -->

We can see that we already did a good job of filtering low-complexity
cells with our prior filtering, but there are a few outlier cells we are
going to exclude from our further analysis

``` r
#filter the cells from the control

list.so$mas <- list.so$mas[, list.so$mas@meta.data[, "complexity"] == 'high']

p2 <- ggplot(list.so$mas@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](mas_seq_with_jc_files/figure-gfm/filter%20low%20complexity%20cells%20control-1.png)<!-- -->

``` r
#filter the cells from the depleted

list.so$jc <- list.so$jc[, list.so$jc@meta.data[, "complexity"] == 'high']

p2 <- ggplot(list.so$jc@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](mas_seq_with_jc_files/figure-gfm/filter%20low%20complexity%20cells%20depleted-1.png)<!-- -->

Next,we are going to use the doublet removal toolkit scDblFinder
(Germain et. al., 2022) to remove doublets from the data. This tool uses
a machine-learning algorithm to simulate artificial doublets from the
data based on cell clustering. From there, real cells get assigned a
doublet score probability to which we will filter cells called as
doublets.

For clustering, we are going to use the Seurat SCTransform workflow.

Here, we want to perform clustering using residual default cutoff of 1.3
(default) rather than selecting a fixed number of highly variable genes.

We will also be scoring cell cycle genes to observe their impact on
clustering as well.

## 3.2 Doublet filtering

### 3.2.1 clustering using SCTransform workflow

``` r
#SCTransform
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- SCTransform(x, variable.features.n = NULL, variable.features.rv.th = 1.3, verbose = F)
})

#cell cyclce scoring
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- CellCycleScoring(x, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
})
```

    ## Warning: The following features are not present in the object: TYMS, DTL, RRM2,
    ## EXO1, DSCC1, E2F8, not searching for symbol synonyms

    ## Warning: The following features are not present in the object: PIMREG, CKAP2L,
    ## HJURP, CDCA3, TTK, CDC25C, KIF2C, CDCA2, NEK2, not searching for symbol synonyms

    ## Warning: The following features are not present in the object: TYMS, DTL, UHRF1,
    ## CDC45, EXO1, E2F8, not searching for symbol synonyms

    ## Warning: The following features are not present in the object: CDK1, PIMREG,
    ## CKAP2L, HJURP, CDCA3, TTK, CDC25C, CDCA2, NEK2, CENPA, not searching for symbol
    ## synonyms

``` r
#PCA
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- RunPCA(x, assay = "SCT")
})
```

    ## PC_ 1 
    ## Positive:  S100A9, S100A8, FTL, LYZ, FTH1, FCN1, S100A12, CST3, CTSS, MNDA 
    ##     IFI30+PIK3R2, TYROBP, S100A6, PSAP, ODF3B+TYMP, CD14, AIF1, LST1, S100A4, CD74 
    ##     HLA-DRB1, VCAN, MS4A6A, SERPINA1, CYBB, FCER1G, HLA-DRA, COTL1, TKT, LGALS1 
    ## Negative:  GNLY, NKG7, IFITM1, IL32, LTB, CCL5, IL7R, MALAT1+NEAT1, CD3E, CD3D 
    ##     GZMA, HLA-A, LDHB, CORO1B+PTPRCAP, PCED1B-AS1, NPM1, EEF1B2, CST7, BTG1, ENSG00000196656 
    ##     GZMB, CD3G, GIMAP7, CD247, KLRB1, MYL12A, CD2, GZMK, ENSG00000237550, ENSG00000268093+RPL18 
    ## PC_ 2 
    ## Positive:  GNLY, NKG7, S100A9, S100A8, CCL5, GZMA, GZMB, CST7, FGFBP2, KLRD1 
    ##     HOPX, GZMH, TRDC, PRF1, CTSW, KLRF1, SPON2, KLRB1, CCL3-AS1+CCL4, KLRC1 
    ##     CLIC3, LAIR2, S100A12, XCL2, HCST, CD160, KLRG1, GZMM, CD247, KLRC3 
    ## Negative:  HLA-DRA, CD74, HLA-DRB1, HLA-DPB1, HLA-DQB1, HLA-DPA1, HLA-DQA1, CD79A, LTB, MS4A1 
    ##     PIK3C2A+RPS13, CD79B, BANK1, CD37, JCHAIN, ENSG00000237550, FTH1, ENSG00000268093+RPL18, EEF1B2, CST3 
    ##     CCR7, LINC00926, CD52, VPREB3, AFF3.1, FTL, LDHB, CD1C, NPM1, CXCR4 
    ## PC_ 3 
    ## Positive:  S100A9, S100A8, IL7R, LDHB, S100A12, CD3E, PIK3C2A+RPS13, LTB, CCR7, ENSG00000237550 
    ##     CD3D, ENSG00000196656, NPM1, LINC02446, GIMAP7, CD8B, BTG1, EEF1B2, ENSG00000268093+RPL18, IL32 
    ##     NOSIP, CD3G, PCED1B-AS1, GAS5, SARAF, MAL, TCF7, LEF1, CD27, PIK3IP1 
    ## Negative:  HLA-DRA, GNLY, HLA-DRB1, CD74, HLA-DPB1, NKG7, HLA-DPA1, HLA-DQB1, FCGR3A, FTH1 
    ##     GZMB, FTL, HLA-DQA1, CST3, FCER1G, SAT1, LST1, GZMA, IFITM3, FGFBP2 
    ##     CD79A, TRDC, IFI30+PIK3R2, IFITM2, KLRD1, AIF1, KLRF1, HOPX, CST7, CCL5 
    ## PC_ 4 
    ## Positive:  S100A9, HLA-DRA, S100A8, CD74, CD79A, HLA-DPB1, HLA-DQB1, HLA-DRB1, MS4A1, JCHAIN 
    ##     LYZ, BANK1, HLA-DQA1, S100A12, CD37, LINC00926, CD79B, VPREB3, IGKC, AFF3.1 
    ##     HLA-DPA1, GZMB, TSPAN13, GNLY, CCDC50, ITM2C, SPIB, CD79A+RPS19, RALGPS2, CD24 
    ## Negative:  FTL, FCGR3A, FTH1, LST1, SAT1, AIF1, IFITM3, LYPD2, MS4A7, PELATON 
    ##     FCER1G, IFI30+PIK3R2, COTL1, HMOX1, SERPINA1, PSAP, CTSS, RHOC, IFITM2, CFD 
    ##     S100A4, NEAT1, TYROBP, SPI1, C1QA, ODF3B+TYMP, CST3, MT2A, S100A6, ACTB 
    ## PC_ 5 
    ## Positive:  PPBP, ENSG00000288796, ENSG00000284874, TUBB1, NRGN, GNG11, CAVIN2, CLU, GP9, CMTM5 
    ##     TREML1, H2AC6, PF4, ITGA2B, CD9, ACRBP, MYL9, SPARC, MAP3K7CL, TSC22D1 
    ##     PTCRA, ENSG00000288869, TMEM40, C2orf88, GP1BB, MYL4, ENSG00000254614, RGS18, LY6G6F, RUFY1 
    ## Negative:  FTL, GNLY, FCGR3A, PIK3C2A+RPS13, IFITM1, LINC02446, EEF1B2, ENSG00000268093+RPL18, IFI30+PIK3R2, LST1 
    ##     IFITM2, ENSG00000196656, HLA-DRB1, ENSG00000237550, CD52, BTG1, CD8B, VIM, GIMAP7, CD3D 
    ##     S100A6, IL7R, S100A4, NPM1, GZMB, IL32, CST3, LYZ, MALAT1+NEAT1, LGALS1

    ## PC_ 1 
    ## Positive:  S100A9, HLA-DRA, S100A12, IFI30, VCAN, MNDA, CD74, AIF1, HLA-DRB1, FTH1 
    ##     S100A8, LYZ, CD14, COTL1, ODF3B+TYMP, APLP2, MS4A6A, LGALS2, HLA-DPB1, CYBB 
    ##     GRN, CTSS, MS4A7, TKT, VIM, SPI1, S100A6, AP1S2, NPC2, NCF2 
    ## Negative:  NKG7, GNLY, GZMA, GZMB, CST7, KLRB1, CD3D, GZMH, HOPX, FGFBP2 
    ##     LTB, PRF1, LAIR2, PTPRCAP, IL7R, GZMM, KLRD1, CLIC3, TRDC, GZMK 
    ##     HCST, CD3G, XCL2, GIMAP7, CTSW, CCL3-AS1+CCL4, CD3E, MALAT1, CD247, S100B 
    ## PC_ 2 
    ## Positive:  GNLY, NKG7, GZMB, S100A9, GZMA, CST7, GZMH, FGFBP2, S100A12, HOPX 
    ##     PRF1, CLIC3, TRDC, KLRD1, VCAN, CCL3-AS1+CCL4, XCL2, LAIR2, KLRF1, CTSW 
    ##     MNDA, IFI30, KLRB1, EFHD2, S100A8, FCGR3A, KLRC1, LYZ, ITGB2, CD160 
    ## Negative:  LTB, IL7R, CD3D, LINC02446, C12orf57, LEF1, CD52, CCR7, MAL, CD3E 
    ##     PRKCQ-AS1, GIMAP7, NOSIP, CD3G, CD8B, OXNAD1, MALAT1, CD27, RGCC, RGS10 
    ##     SNHG6, SNHG8, SNHG32, ISG20, FCMR, RCAN3, SNRPD2, CAMK4, GYPC, MIF 
    ## PC_ 3 
    ## Positive:  S100A9, S100A12, VCAN, S100A8, MNDA, CD14, CD3D, IL7R, GIMAP7, LYZ 
    ##     LINC02446, RBP7, APLP2, CTSD, RGS2, CD3E, CD3G, CTSS, FYB1, GIMAP4 
    ##     PLBD1, VIM, TKT, CEBPD, CD8B, ANXA1, PRKCQ-AS1, LEF1, S100A6, NCF2 
    ## Negative:  HLA-DRA, CD74, HLA-DPB1, HLA-DRB1, CD79A, HLA-DQB1, HLA-DPA1, MS4A1, HLA-DQA1, BANK1 
    ##     JCHAIN, VPREB3, CD79B, AFF3, HLA-DMA, FCER1A, CD37, IGKC, GZMB, RALGPS2 
    ##     SPIB, PPP1R14A, CD1C, CCDC50, TSPAN13, TNFRSF13C, FCER2, EAF2, CD22, IRF8 
    ## PC_ 4 
    ## Positive:  FCGR3A, AIF1, LYPD2, MS4A7, PELATON, FTH1, IFI30, CDKN1C, RHOC, HMOX1 
    ##     COTL1, TIMP1, MTSS1, HES4, BCL2A1, WARS1, C1QA, IFITM2, CFD, PPM1N 
    ##     VMO1, SPI1, CALHM6, TUBA1B, DUSP6, NPC2, STXBP2, CX3CR1, S100A11, ISG15 
    ## Negative:  S100A9, S100A12, VCAN, S100A8, CD79A, JCHAIN, LYZ, MS4A1, MNDA, CD74 
    ##     BANK1, CD14, VPREB3, MS4A6A, IGKC, GZMB, ITM2C, LTB, MZB1, AFF3 
    ##     HLA-DRA, HLA-DQB1, CCDC50, APLP2, NCF1, TSPAN13, SPIB, PLBD1, CD36, RALGPS2 
    ## PC_ 5 
    ## Positive:  PF4, ENSG00000284874, GNG11, GP9, NRGN, CMTM5, ACRBP, PPBP, CLEC1B, CLU 
    ##     PRKAR2B, H2AC6, TMEM40, PTCRA, C2orf88, ITGA2B, HRAT92, TSC22D1, ENSG00000288869, TUBB1 
    ##     LY6G6F, TREML1, CD9, SPARC, ENSG00000236304, GP1BB, ITGB3, MTURN, RGS18, MYL9 
    ## Negative:  S100A9, S100A12, LTB, HLA-DRA, VCAN, CD74, HLA-DRB1, CD3D, S100A8, HLA-DPB1 
    ##     CD52, MNDA, IL7R, S100B, LINC02446, LYZ, CD79A, VIM, HCST, CD14 
    ##     CD8B, HLA-DQB1, PTPRCAP, GIMAP7, ISG20, MIF, CD3G, GZMK, KLRB1, MALAT1

``` r
#generate UMAP coordinates
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- RunUMAP(x, dims = 1:30)
})
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 15:40:03 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 15:40:03 Read 9846 rows and found 30 numeric columns

    ## 15:40:03 Using Annoy for neighbor search, n_neighbors = 30

    ## 15:40:03 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 15:40:04 Writing NN index file to temp file C:\Users\dante\AppData\Local\Temp\RtmpAvdq4I\file6cbc4a55b2d
    ## 15:40:05 Searching Annoy index using 1 thread, search_k = 3000
    ## 15:40:08 Annoy recall = 100%
    ## 15:40:09 Commencing smooth kNN distance calibration using 1 thread
    ## 15:40:10 Initializing from normalized Laplacian + noise
    ## 15:40:11 Commencing optimization for 500 epochs, with 426142 positive edges
    ## 15:40:42 Optimization finished
    ## 15:40:42 UMAP embedding parameters a = 0.9922 b = 1.112
    ## 15:40:42 Read 9580 rows and found 30 numeric columns
    ## 15:40:42 Using Annoy for neighbor search, n_neighbors = 30
    ## 15:40:42 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 15:40:44 Writing NN index file to temp file C:\Users\dante\AppData\Local\Temp\RtmpAvdq4I\file6cbc38502824
    ## 15:40:44 Searching Annoy index using 1 thread, search_k = 3000
    ## 15:40:47 Annoy recall = 100%
    ## 15:40:48 Commencing smooth kNN distance calibration using 1 thread
    ## 15:40:49 Initializing from normalized Laplacian + noise
    ## 15:40:49 Commencing optimization for 500 epochs, with 409506 positive edges
    ## 15:41:20 Optimization finished

``` r
#find k-nearest neighbors
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- FindNeighbors(x, dims = 1:30)
})
```

    ## Computing nearest neighbor graph
    ## Computing SNN
    ## Computing nearest neighbor graph
    ## Computing SNN

``` r
#Find clusters using the louvain algorithm with multilevel refinement. It is recommended to overcluster the data first when using scDblFinder
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- FindClusters(x, resolution = 1.2, algorithm = 2)
})
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 9846
    ## Number of edges: 385869
    ## 
    ## Running Louvain algorithm with multilevel refinement...
    ## Maximum modularity in 10 random starts: 0.8401
    ## Number of communities: 19
    ## Elapsed time: 2 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 9580
    ## Number of edges: 365481
    ## 
    ## Running Louvain algorithm with multilevel refinement...
    ## Maximum modularity in 10 random starts: 0.8562
    ## Number of communities: 22
    ## Elapsed time: 2 seconds

### 3.2.2 scDblFinder doublet calling

We use the natural log normalized features to simulate artificial
doublets. The expected doublet rate is assumed to be 1% per thousand
cells captured which is appropriate for 10x datasets.

``` r
library(scDblFinder)
library(SingleCellExperiment)
```

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:sp':
    ## 
    ##     %over%

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomeInfoDb

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'SummarizedExperiment'

    ## The following object is masked from 'package:SeuratObject':
    ## 
    ##     Assays

    ## The following object is masked from 'package:Seurat':
    ## 
    ##     Assays

``` r
#natural log normalize the raw counts data. SCTransform counts data uses pearson residuals which can only be used for clustering/visualization
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- NormalizeData(x, assay = 'RNA')
}) 

#Run scDblFinder. We have to first convert to a single cell experiment object for this tool
mas <- as.SingleCellExperiment(list.so$mas, assay = 'RNA')

jc <- as.SingleCellExperiment(list.so$jc, assay = 'RNA')

#convert to a list
list.dbr <- list(mas,jc)

rm(mas,jc)

d <- lapply(list.dbr, FUN = function(x) {
  x <- scDblFinder(x, clusters = 'seurat_clusters', dbr = NULL, dims = 30, includePCs = 30, returnType = "table", processing = "normFeatures", nfeatures = 2000)
})
```

    ## Warning in .checkSCE(sce): Some cells in `sce` have an extremely low read
    ## counts; note that these could trigger errors and might best be filtered out

    ## 19 clusters

    ## Creating ~7877 artificial doublets...

    ## Dimensional reduction

    ## Evaluating kNN...

    ## Training model...

    ## iter=0, 704 cells excluded from training.

    ## iter=1, 430 cells excluded from training.

    ## iter=2, 355 cells excluded from training.

    ## Threshold found:0.411

    ## 374 (3.8%) doublets called

    ## Warning in .checkSCE(sce): Some cells in `sce` have an extremely low read
    ## counts; note that these could trigger errors and might best be filtered out

    ## 22 clusters

    ## Creating ~7664 artificial doublets...

    ## Dimensional reduction

    ## Evaluating kNN...

    ## Training model...

    ## iter=0, 644 cells excluded from training.

    ## iter=1, 403 cells excluded from training.

    ## iter=2, 376 cells excluded from training.

    ## Threshold found:0.446

    ## 389 (4.1%) doublets called

For cells called doublets, we should see twice the amount of UMI counts
and more genes/cell. We have view this with Violin plots to check the
doublet calls.

``` r
#plotting control doublet calls
list.so$mas$class <- as.data.frame(d[[1]]) %>% filter(type == "real") %>% select(class)
list.so$mas$class <- as.factor(list.so$mas$class)
list.so$mas$class <- factor(list.so$mas$class, levels = c("singlet","doublet"))

DimPlot(list.so$mas, group.by = c("class"), order = T)
```

![](mas_seq_with_jc_files/figure-gfm/doublet%20calls%20control-1.png)<!-- -->

``` r
VlnPlot(list.so$mas, features = c("nCount_RNA", "nFeature_RNA"), group.by = "class")
```

![](mas_seq_with_jc_files/figure-gfm/doublet%20calls%20control-2.png)<!-- -->

``` r
#plotting depleted doublet calls
list.so$jc$class <- as.data.frame(d[[2]]) %>% filter(type == "real") %>% select(class)
list.so$jc$class <- as.factor(list.so$jc$class)
list.so$jc$class <- factor(list.so$jc$class, levels = c("singlet","doublet"))

DimPlot(list.so$jc, group.by = c("class"), order = T)
```

![](mas_seq_with_jc_files/figure-gfm/doublet%20calls%20depleted-1.png)<!-- -->

``` r
VlnPlot(list.so$jc, features = c("nCount_RNA", "nFeature_RNA"), group.by = "class")
```

![](mas_seq_with_jc_files/figure-gfm/doublet%20calls%20depleted-2.png)<!-- -->

The doublet call information looks like what we would expect. Nowe we
will filter these cells from each sample.

``` r
#remove doublets control
list.so$mas <- list.so$mas[, list.so$mas@meta.data[, "class"] == "singlet"]

#remove doublets depleted
list.so$jc <- list.so$jc[, list.so$jc@meta.data[, "class"] == "singlet"]
```

## 3.3 Cell cycle evaluation

Depending on the experiment, cell cycle related influence may contribute
to uninteresting variation in clustering. As a result, we can choose to
regress the influence of cell cycle related genes in clustering.

We will look at the S Phase and G2 Phase scores separately.
Additionally, we can separate all cycling cells (S + G2 phase) cells
from all non-cycling cells by subtracting the G2M and S scores.

``` r
#we quantified cells is S and G2 phase earlier, so we will look to see the proportion of cycling cells
#control
list.so$mas$cc.difference <- list.so$mas$S.Score - list.so$mas$G2M.Score
VlnPlot(list.so$mas, features = c('S.Score','G2M.Score','cc.difference'), group.by = 'orig.ident')
```

![](mas_seq_with_jc_files/figure-gfm/cell%20cycling%20control-1.png)<!-- -->

``` r
#we quantified cells is S and G2 phase earlier, so we will look to see the proportion of cycling cells
#depleted
list.so$jc$cc.difference <- list.so$jc$S.Score - list.so$jc$G2M.Score
VlnPlot(list.so$jc, features = c('S.Score','G2M.Score','cc.difference'), group.by = 'orig.ident')
```

![](mas_seq_with_jc_files/figure-gfm/cell%20cycling%20depleted-1.png)<!-- -->

We can see there is little influence of these PBMC cells. We can always
choose to regress this influence out to our choosing when selecting
vars.to.regress during SCTransform. For now, we won’t regress cell cycle
genes from the data.

# 4 Cell clustering using SCTransform workflow

Now that we have done our cell QC filtering, we will repeat the
SCTransform workflow to produce clustring results without the influence
of doublets

``` r
#SCTransform
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- SCTransform(x, variable.features.n = NULL, variable.features.rv.th = 1.3, verbose = F)
})

#PCA
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- RunPCA(x, assay = "SCT")
})
```

    ## PC_ 1 
    ## Positive:  S100A9, S100A8, FTL, LYZ, FTH1, FCN1, CST3, S100A12, CTSS, IFI30+PIK3R2 
    ##     MNDA, TYROBP, S100A6, PSAP, ODF3B+TYMP, CD14, AIF1, LST1, S100A4, CD74 
    ##     HLA-DRB1, VCAN, SERPINA1, MS4A6A, FCER1G, HLA-DRA, CYBB, COTL1, TKT, LGALS1 
    ## Negative:  GNLY, NKG7, IFITM1, MALAT1, IL32, LTB, IL7R, MALAT1+NEAT1, CCL5, CD3E 
    ##     CD3D, HLA-A, LDHB, GZMA, NPM1, CORO1B+PTPRCAP, PCED1B-AS1, EEF1B2, BTG1, ENSG00000196656 
    ##     GIMAP7, CST7, CD3G, GZMB, MYL12A, CD247, KLRB1, ENSG00000237550, CD2, GZMK 
    ## PC_ 2 
    ## Positive:  GNLY, NKG7, S100A9, S100A8, GZMA, CCL5, GZMB, CST7, FGFBP2, KLRD1 
    ##     HOPX, GZMH, TRDC, PRF1, CTSW, KLRF1, SPON2, KLRB1, CCL3-AS1+CCL4, KLRC1 
    ##     CLIC3, LAIR2, XCL2, CD160, HCST, KLRG1, GZMM, IFITM2, S100A12, KLRC3 
    ## Negative:  HLA-DRA, CD74, HLA-DRB1, HLA-DPB1, LTB, HLA-DQB1, EEF1A1, CD79A, HLA-DPA1, PIK3C2A+RPS13 
    ##     HLA-DQA1, ENSG00000237550, MS4A1, ENSG00000268093+RPL18, EEF1B2, CD79B, BANK1, EEF1G, CCR7, CD37 
    ##     LDHB, NPM1, JCHAIN, ENSG00000196656, GAS5, IL7R, CD52, CXCR4, MALAT1, LINC00926 
    ## PC_ 3 
    ## Positive:  S100A9, S100A8, IL7R, S100A12, LDHB, CD3E, PIK3C2A+RPS13, CCR7, CD3D, LTB 
    ##     ENSG00000237550, GIMAP7, LINC02446, IL32, ENSG00000196656, CD8B, NPM1, BTG1, EEF1A1, LYZ 
    ##     CD3G, PCED1B-AS1, EEF1B2, NOSIP, ENSG00000268093+RPL18, MAL, VCAN, SARAF, TCF7, GAS5 
    ## Negative:  HLA-DRA, HLA-DRB1, CD74, HLA-DPB1, GNLY, HLA-DPA1, HLA-DQB1, FCGR3A, HLA-DQA1, NKG7 
    ##     FTH1, FTL, CST3, GZMB, CD79A, SAT1, LST1, FCER1G, MS4A1, CD79B 
    ##     IFITM3, BANK1, AIF1, IFI30+PIK3R2, JCHAIN, RHOC, IFITM2, GZMA, MS4A7, COTL1 
    ## PC_ 4 
    ## Positive:  S100A9, S100A8, HLA-DRA, CD74, CD79A, HLA-DPB1, HLA-DQB1, HLA-DRB1, MS4A1, LYZ 
    ##     JCHAIN, BANK1, HLA-DQA1, S100A12, CD79B, CD37, LINC00926, VPREB3, AFF3.1, IGKC 
    ##     GNLY, TSPAN13, CD79A+RPS19, GZMB, RALGPS2, SPIB, HLA-DPA1, CCDC50, CD24, FCER2 
    ## Negative:  FTL, FCGR3A, FTH1, LST1, SAT1, AIF1, IFITM3, LYPD2, MS4A7, PELATON 
    ##     FCER1G, IFI30+PIK3R2, COTL1, HMOX1, SERPINA1, RHOC, PSAP, CTSS, IFITM2, S100A4 
    ##     CFD, NEAT1, SPI1, TYROBP, C1QA, ODF3B+TYMP, CST3, ACTB, MT2A, S100A6 
    ## PC_ 5 
    ## Positive:  EEF1A1, GNLY, FTL, FCGR3A, PIK3C2A+RPS13, IFITM1, LINC02446, EEF1B2, IFI30+PIK3R2, IFITM2 
    ##     MALAT1, HLA-DRB1, ENSG00000268093+RPL18, LST1, GZMB, ENSG00000196656, ENSG00000237550, BTG1, HLA-DRA, CD8B 
    ##     S100A6, CD52, EEF1G, VIM, IL7R, S100A4, NPM1, GIMAP7, CST3, CD74 
    ## Negative:  PPBP, ENSG00000288796, ENSG00000284874, TUBB1, NRGN, GNG11, CAVIN2, GP9, CLU, CMTM5 
    ##     TREML1, H2AC6, PF4, ITGA2B, CD9, ACRBP, SPARC, MAP3K7CL, TSC22D1, MYL9 
    ##     TMEM40, ENSG00000288869, PTCRA, C2orf88, RGS18, ENSG00000254614, GP1BB, LY6G6F, MYL4, RUFY1

    ## PC_ 1 
    ## Positive:  S100A9, HLA-DRA, S100A12, VCAN, IFI30, MNDA, CD74, AIF1, FTH1, HLA-DRB1 
    ##     S100A8, LYZ, CD14, COTL1, ODF3B+TYMP, MS4A6A, APLP2, LGALS2, CYBB, HLA-DPB1 
    ##     GRN, CTSS, MS4A7, TKT, VIM, SPI1, S100A6, NCF2, AP1S2, CPVL 
    ## Negative:  NKG7, GNLY, GZMA, GZMB, CST7, KLRB1, CD3D, GZMH, HOPX, FGFBP2 
    ##     PRF1, LAIR2, LTB, PTPRCAP, IL7R, GZMM, KLRD1, CLIC3, TRDC, GZMK 
    ##     CD3G, HCST, XCL2, GIMAP7, CD3E, CTSW, CCL3-AS1+CCL4, CD247, MALAT1, KLRF1 
    ## PC_ 2 
    ## Positive:  GNLY, NKG7, GZMB, S100A9, GZMA, CST7, GZMH, FGFBP2, HOPX, PRF1 
    ##     CLIC3, TRDC, S100A12, KLRD1, VCAN, LAIR2, CCL3-AS1+CCL4, XCL2, KLRF1, CTSW 
    ##     MNDA, KLRB1, IFI30, EFHD2, S100A8, FCGR3A, KLRC1, LYZ, ITGB2, CD160 
    ## Negative:  LTB, IL7R, CD3D, LINC02446, C12orf57, LEF1, CD52, CCR7, CD3E, MAL 
    ##     PRKCQ-AS1, NOSIP, GIMAP7, CD3G, CD8B, OXNAD1, MALAT1, CD27, RGCC, RGS10 
    ##     SNHG6, SNHG8, SNHG32, FCMR, RCAN3, CAMK4, ISG20, SNRPD2, GYPC, OCIAD2 
    ## PC_ 3 
    ## Positive:  S100A9, S100A12, VCAN, S100A8, MNDA, CD14, CD3D, LYZ, IL7R, GIMAP7 
    ##     RBP7, LINC02446, CTSD, APLP2, RGS2, CD3E, CTSS, CD3G, VIM, CEBPD 
    ##     TKT, FYB1, PLBD1, GIMAP4, ANXA1, CD8B, S100A6, NCF2, PRKCQ-AS1, LEF1 
    ## Negative:  HLA-DRA, CD74, HLA-DPB1, HLA-DRB1, CD79A, HLA-DQB1, HLA-DPA1, MS4A1, HLA-DQA1, BANK1 
    ##     JCHAIN, VPREB3, CD79B, AFF3, HLA-DMA, CD37, IGKC, FCER1A, GZMB, RALGPS2 
    ##     SPIB, PPP1R14A, TSPAN13, CD1C, TNFRSF13C, FCER2, CCDC50, CD22, EAF2, IGHV5-78 
    ## PC_ 4 
    ## Positive:  FCGR3A, AIF1, LYPD2, MS4A7, PELATON, FTH1, IFI30, CDKN1C, RHOC, HMOX1 
    ##     COTL1, TIMP1, MTSS1, BCL2A1, WARS1, C1QA, IFITM2, CFD, PPM1N, CALHM6 
    ##     SPI1, TUBA1B, VMO1, DUSP6, NPC2, STXBP2, CX3CR1, BRI3, ISG15, C1QB 
    ## Negative:  S100A9, S100A12, VCAN, S100A8, CD79A, JCHAIN, LYZ, MNDA, MS4A1, CD74 
    ##     BANK1, CD14, VPREB3, MS4A6A, IGKC, LTB, ITM2C, MZB1, GZMB, AFF3 
    ##     HLA-DQB1, APLP2, NCF1, HLA-DRA, PLBD1, CD36, CCDC50, TSPAN13, SPIB, NCF1C 
    ## PC_ 5 
    ## Positive:  S100A9, S100A12, LTB, HLA-DRA, VCAN, CD74, CD3D, HLA-DRB1, S100A8, CD52 
    ##     MNDA, HLA-DPB1, IL7R, CD79A, S100B, CD14, VIM, LINC02446, LYZ, HCST 
    ##     PTPRCAP, CD8B, ISG20, GIMAP7, HLA-DQB1, MIF, KLRB1, CD3G, IFI30, GZMK 
    ## Negative:  PF4, ENSG00000284874, GNG11, GP9, NRGN, CMTM5, ACRBP, PPBP, CLU, CLEC1B 
    ##     PRKAR2B, TMEM40, H2AC6, PTCRA, C2orf88, HRAT92, TSC22D1, ITGA2B, LY6G6F, ENSG00000288869 
    ##     TUBB1, CD9, TREML1, SPARC, GP1BB, ENSG00000236304, ITGB3, MTURN, RGS18, F13A1

``` r
#generate UMAP coordinates
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- RunUMAP(x, dims = 1:30)
})
```

    ## 15:52:41 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 15:52:41 Read 9472 rows and found 30 numeric columns

    ## 15:52:41 Using Annoy for neighbor search, n_neighbors = 30

    ## 15:52:41 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 15:52:43 Writing NN index file to temp file C:\Users\dante\AppData\Local\Temp\RtmpAvdq4I\file6cbc5f256dc5
    ## 15:52:43 Searching Annoy index using 1 thread, search_k = 3000
    ## 15:52:46 Annoy recall = 100%
    ## 15:52:47 Commencing smooth kNN distance calibration using 1 thread
    ## 15:52:49 Found 2 connected components, falling back to 'spca' initialization with init_sdev = 1
    ## 15:52:49 Initializing from PCA
    ## 15:52:49 Using 'irlba' for PCA
    ## 15:52:49 PCA: 2 components explained 43.94% variance
    ## 15:52:49 Commencing optimization for 500 epochs, with 405848 positive edges
    ## 15:53:20 Optimization finished
    ## 15:53:20 UMAP embedding parameters a = 0.9922 b = 1.112
    ## 15:53:20 Read 9191 rows and found 30 numeric columns
    ## 15:53:20 Using Annoy for neighbor search, n_neighbors = 30
    ## 15:53:20 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 15:53:21 Writing NN index file to temp file C:\Users\dante\AppData\Local\Temp\RtmpAvdq4I\file6cbc1fe46f4e
    ## 15:53:21 Searching Annoy index using 1 thread, search_k = 3000
    ## 15:53:24 Annoy recall = 100%
    ## 15:53:25 Commencing smooth kNN distance calibration using 1 thread
    ## 15:53:27 Found 2 connected components, falling back to 'spca' initialization with init_sdev = 1
    ## 15:53:27 Initializing from PCA
    ## 15:53:27 Using 'irlba' for PCA
    ## 15:53:27 PCA: 2 components explained 40.5% variance
    ## 15:53:27 Commencing optimization for 500 epochs, with 389032 positive edges
    ## 15:53:57 Optimization finished

``` r
#find k-nearest neighbors. For consistency, are going to use the same number of neighbors used in RunUMAP()
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- FindNeighbors(x, dims = 1:30, k.param = 30)
})
```

    ## Computing nearest neighbor graph
    ## Computing SNN
    ## Computing nearest neighbor graph
    ## Computing SNN

``` r
#Find clusters using the louvain algorithm with multilevel refinement. Each clustering resolution will be stored in the metadata
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- FindClusters(x, resolution = seq(0.1,1.2,0.1), algorithm = 2, verbose=F)
})
```

We are going to initially visualize the UMAP plots of each sample
separately. We are going to choose a relatively high resolution to
showcase that with Jumpcode depletion we can extract more heterogeneity
with samples that haven’t approached sequence saturation.

``` r
Idents(list.so$mas) <- "SCT_snn_res.1"
Idents(list.so$jc) <- "SCT_snn_res.1"

#visualize UMAP of the control
DimPlot(list.so$mas, label = T, pt.size = 0.75) + theme_classic() + ggtitle('MAS')
```

![](mas_seq_with_jc_files/figure-gfm/UMAP%20plots%20of%20clustering%20results-1.png)<!-- -->

``` r
#visualize UMAP of the depleted
DimPlot(list.so$jc, label = T, pt.size = 0.75) + theme_classic() + ggtitle('JC')
```

![](mas_seq_with_jc_files/figure-gfm/UMAP%20plots%20of%20clustering%20results-2.png)<!-- -->

# 5 Harmony data projection

From the UMAP plots, we see an additional 4 clusters!

Depending on the seed and CPU processing, UMAPs can visually be
generated differently between users. Additionally, UMAPs between control
and depleted can sometimes be difficult to visually compare.

As a result, we also want to ensure that the control and depleted
conditions can be projected together so that experiments used with
depletion can be compared to other previously sequenced samples.

To do so, we are going to use the Harmony algorithm (Korsunsky et. al.,
2019) to project the control and depleted samples onto one another.

DO NOT attempt to find nearest neighbors or perform clustering on the
harmony dimensions. This will effectively “batch correct” the observable
clustering differences between the control and depleted conditions which
we ultimately want to examine.

Thus, use the clustering results generated prior running harmony.

``` r
library(harmony)
```

    ## Warning: package 'harmony' was built under R version 4.2.2

    ## Loading required package: Rcpp

``` r
#add condition variables
list.so$mas$condition <- "MAS-Seq"
list.so$jc$condition <- "MAS-Seq + JC"

#set the default assay to the SCT slot
DefaultAssay(list.so$mas) <- "SCT"
DefaultAssay(list.so$jc) <- "SCT"

#select variable features that are in common between each sample to aid in projection
features <- SelectIntegrationFeatures(object.list = list.so, nfeatures = 2000)

#merge seurat
temp <- merge(x = list.so[[1]], y = list.so[[2]], merge.data=TRUE)
```

    ## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
    ## duplicated across objects provided. Renaming to enforce unique cell names.

``` r
#again make sure the default assay is the SCT slot after merging data
DefaultAssay(temp) <- "SCT"

#manually set variable features previously calculated
VariableFeatures(temp) <- features

#Run PCA on the merged data
temp <- RunPCA(temp)
```

    ## PC_ 1 
    ## Positive:  S100A9, S100A8, LYZ, FTL, S100A12, FTH1, MNDA, HLA-DRA, VCAN, CD74 
    ##     AIF1, HLA-DRB1, FCN1, CTSS, ODF3B+TYMP, CD14, IFI30, CST3, S100A6, TYROBP 
    ##     COTL1, MS4A6A, CYBB, APLP2, TKT, SPI1, LGALS2, VIM, NCF2, PSAP 
    ## Negative:  GNLY, NKG7, GZMA, CD3D, LTB, GZMB, CST7, IL7R, KLRB1, MALAT1 
    ##     CD3E, GZMH, CCL5, CD3G, GIMAP7, FGFBP2, HOPX, GZMK, IL32, GZMM 
    ##     PRF1, KLRD1, CD247, LAIR2, CD2, C12orf57, TRDC, CD52, CTSW, IFITM1 
    ## PC_ 2 
    ## Positive:  GNLY, NKG7, S100A9, GZMB, GZMA, S100A8, CST7, GZMH, FGFBP2, HOPX 
    ##     PRF1, KLRD1, TRDC, CLIC3, CTSW, KLRF1, KLRB1, CCL3-AS1+CCL4, LAIR2, XCL2 
    ##     CCL5, KLRC1, S100A12, SPON2, EFHD2, LYZ, CD160, HCST, TYROBP, CMC1 
    ## Negative:  LTB, IL7R, CCR7, LINC02446, CD52, LEF1, CD79A, MAL, CD3E, MALAT1 
    ##     C12orf57, CD3D, NOSIP, EEF1A1, HLA-DRA, CD8B, OXNAD1, CD27, SARAF, PRKCQ-AS1 
    ##     ENSG00000196656, RGS10, MS4A1, ENSG00000237550, LDHB, RGCC, CXCR4, SNHG8, CD3G, GIMAP7 
    ## PC_ 3 
    ## Positive:  HLA-DRA, CD74, HLA-DRB1, HLA-DPB1, HLA-DQB1, HLA-DPA1, CD79A, HLA-DQA1, MS4A1, BANK1 
    ##     CD79B, JCHAIN, FCGR3A, GNLY, VPREB3, GZMB, HLA-DMA, CD37, CD1C, IGKC 
    ##     NKG7, CCDC50, SPIB, IRF8, TSPAN13, LINC00926, RALGPS2, CLEC10A, FCER2, CALHM6 
    ## Negative:  S100A9, S100A8, S100A12, IL7R, VCAN, CD3D, GIMAP7, CD3E, LINC02446, LYZ 
    ##     CD14, CD3G, MNDA, CD8B, LEF1, NOSIP, PRKCQ-AS1, MAL, CCR7, C12orf57 
    ##     OXNAD1, GIMAP4, FCN1, IL32, RBP7, GIMAP5, CD27, RGCC, CD2, CAMK4 
    ## PC_ 4 
    ## Positive:  S100A9, S100A8, CD79A, S100A12, CD74, LYZ, JCHAIN, MS4A1, HLA-DRA, BANK1 
    ##     VCAN, HLA-DQB1, HLA-DPB1, VPREB3, LTB, IGKC, ITM2C, CD37, TSPAN13, RALGPS2 
    ##     GZMB, SPIB, CD79B, CCDC50, LINC00926, MZB1, CD14, MNDA, HLA-DRB1, HLA-DQA1 
    ## Negative:  FCGR3A, AIF1, FTH1, LYPD2, MS4A7, FTL, PELATON, LST1, COTL1, IFITM3 
    ##     CDKN1C, RHOC, HMOX1, SAT1, IFI30, FCER1G, IFITM2, TIMP1, CFD, WARS1 
    ##     BCL2A1, C1QA, MTSS1, SPI1, SERPINA1, DUSP6, PPM1N, TUBA1B, CTSS, VMO1 
    ## PC_ 5 
    ## Positive:  ENSG00000284874, PF4, GNG11, PPBP, NRGN, GP9, CMTM5, TUBB1, CLU, ACRBP 
    ##     H2AC6, ITGA2B, ENSG00000288796, TREML1, CAVIN2, PTCRA, TMEM40, CD9, TSC22D1, PRKAR2B 
    ##     C2orf88, CLEC1B, ENSG00000288869, SPARC, LY6G6F, GP1BB, RGS18, MYL9, HRAT92, ENSG00000236304 
    ## Negative:  HLA-DRB1, HLA-DRA, LTB, CD52, CD3D, LINC02446, IL7R, FCGR3A, CD74, HLA-DPB1 
    ##     AIF1, CD8B, S100A9, S100B, MALAT1, GIMAP7, EEF1A1, IFI30, FTL, VIM 
    ##     S100A6, HLA-DPA1, HLA-DQB1, HCST, MIF, CD3E, LST1, COTL1, IFITM2, S100A8

``` r
#run harmony to correct the PCA
harmonized_so <- RunHarmony(temp, group.by.vars = "condition", reduction.save="harmony", reduction="pca")
```

    ## Harmony 1/10

    ## Harmony 2/10

    ## Harmony 3/10

    ## Harmony 4/10

    ## Harmony 5/10

    ## Harmony 6/10

    ## Harmony 7/10

    ## Harmony 8/10

    ## Harmony 9/10

    ## Harmony 10/10

    ## Harmony converged after 10 iterations

    ## Warning: Invalid name supplied, making object name syntactically valid. New
    ## object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details
    ## on syntax validity

``` r
#remove the temp merged object
rm(temp)

#run umap on the corrected PCA
harmonized_so <- RunUMAP(harmonized_so, dims=1:30, reduction = "harmony", n.neighbors=30)
```

    ## 15:56:53 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 15:56:53 Read 18663 rows and found 30 numeric columns

    ## 15:56:53 Using Annoy for neighbor search, n_neighbors = 30

    ## 15:56:53 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 15:56:56 Writing NN index file to temp file C:\Users\dante\AppData\Local\Temp\RtmpAvdq4I\file6cbc690f304b
    ## 15:56:56 Searching Annoy index using 1 thread, search_k = 3000
    ## 15:57:03 Annoy recall = 100%
    ## 15:57:04 Commencing smooth kNN distance calibration using 1 thread
    ## 15:57:07 Initializing from normalized Laplacian + noise
    ## 15:57:09 Commencing optimization for 200 epochs, with 811482 positive edges
    ## 15:57:33 Optimization finished

``` r
#visualize the projected data split by condition
DimPlot(harmonized_so, label = T, reduction = "umap", split.by = "condition", repel = T, label.size = 3, pt.size = 0.75) & NoLegend()
```

![](mas_seq_with_jc_files/figure-gfm/haromy%20projection-1.png)<!-- -->

Now we can visually compare the clustering results between control and
depletion!
