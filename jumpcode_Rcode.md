CRISPRclean Data Analysis Using R
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
        -   [3.1.2 Quantify % of reads contributed by mitochondrial,
            ribosomoal, and all jumpcode gene
            targets](#312-quantify--of-reads-contributed-by-mitochondrial-ribosomoal-and-all-jumpcode-gene-targets)
        -   [3.1.3 Dead cell removal](#313-dead-cell-removal)
        -   [3.1.4 Remove cells with too few genes
            detected](#314-remove-cells-with-too-few-genes-detected)
        -   [3.1.5 Filtering low-complexity
            cells](#315-filtering-low-complexity-cells)
    -   [3.2 Doublet filtering](#32-doublet-filtering)
        -   [3.2.1 clustering using SCTransform
            workflow](#321-clustering-using-sctransform-workflow)
        -   [3.2.2 scDblFinder doublet
            calling](#322-scdblfinder-doublet-calling)
    -   [3.3 Cell cycle evaluation](#33-cell-cycle-evaluation)
-   [4 Cell clustering using SCTransform
    workflow](#4-cell-clustering-using-sctransform-workflow)

# 1 Plotting median UMIs/cell and Genes/cell

We first want to assess some basic metrics such as UMIs/cell and
Genes/cell to assess the gain in transcriptomic read information of
untargeted genes.

We created a series of functions to assist with the workflow in
assessing benefits of depletion.

You can use a directory containing the barcodes.tsv, genes.tsv, and
matrix.mtx file as output from CellRanger or
filtered_feature_barcode_matrix.h5 file from CellRanger

You will also need the gene target list provided by Jumpcode Genomics.
Because we are targeting \~350 genes for removal, we experimentally
remove UMIs associated with the targeted transcripts. We want to
evaluate how reads are being redistributed for genes that have not been
targeted with CRISPR. Thus, we should see a boost in both UMIs/cell and
Genes/cell when considering all non-targeted transcripts.

``` r
#We are first going to create a function to evaluate both umis/cell and genes/cell for each sample
umis_genes_per_cell <- function(matrix, gene_list, sample = c('10x-v3','CRISPRClean')){
    x<-c("Seurat", "dplyr", "patchwork")
    lapply(x, require, character.only = TRUE)
    if(dir.exists(paths = matrix)) {
        mtx <- Read10X(matrix)
    } else {
        mtx <- Read10X_h5(matrix)
    }
    if(file.exists(gene_list)) {
        targets <- read.delim(gene_list, header = F)
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
#matrix = PATH to directory containing .h5 file
#gene list = gene list provided by Jumpcode containing gene targets
#sample = whichever name you would like to name the conditions to compare (i.e., 10x-v3 vs Jumpcode Depletion)

#control
control <- umis_genes_per_cell(matrix = "~/R/control_filtered_no_mask.h5", gene_list = "~/targets.txt", sample = "10x-v3")
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
depleted <- umis_genes_per_cell(matrix = "~/R/depleted_filtered_no_mask_rep3.h5", gene_list = "~/targets.txt", sample = "CRISPRclean")
```

``` r
#plot histogram of UMIs and genes/cell
depletion_benefit <- function(control,depleted) {
    require("ggplot2")
    df.m <- rbind(control, depleted)
    median <- df.m %>% dplyr::group_by(condition) %>% dplyr::summarize(median_umi=round(median(umi_counts)), median_genes=round(median(gene_counts)))
    
    p1 <- ggplot(df.m, aes(x =umi_counts , fill=condition, ..scaled..)) + geom_density(alpha=0.2) + 
        xlab("nUMIs/cell") + ylab("Scaled Density") + xlim(0,median(df.m$umi_counts) + 8*mad(df.m$umi_counts)) +
        geom_vline(data = median, aes(xintercept = median_umi, color=condition), linetype='dashed') + geom_text(data = median, aes(x = round(median(df.m$umi_counts) + 1*mad(df.m$umi_counts)), y = c(1,0.9), label=paste0(condition, ":", median_umi), color=condition), size=4, hjust=0, fontface=2, show.legend = F)
        
      #this will generate plot comparing genes/cell
    p2 <- ggplot(df.m, aes(x =gene_counts , fill=condition, ..scaled..)) + geom_density(alpha=0.2) + 
        xlab("nGenes/cell") + ylab(NULL) +
        geom_vline(data = median, aes(xintercept = median_genes, color=condition), linetype='dashed') + xlim(0,median(df.m$gene_counts) + 8*mad(df.m$gene_counts)) + geom_text(data = median, aes(x = round(median(df.m$gene_counts) + 1*mad(df.m$gene_counts)), y = c(1,0.9), label=paste0(condition, ":", median_genes), color=condition), size=4, hjust=0, fontface=2, show.legend = F)
    p1 + p2 + plot_layout(guides = "collect", widths = c(2,2)) & theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(), axis.title.y = element_text(size = 14), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white'), legend.title = element_blank())
}
```

``` r
#generate the plots
depletion_benefit(control = control, depleted = depleted)
```

    ## Loading required package: ggplot2

    ## Warning: Removed 278 rows containing non-finite values (stat_density).

    ## Warning: Removed 15 rows containing non-finite values (stat_density).

![](jumpcode_Rcode_files/figure-gfm/plot%20the%20function-1.png)<!-- -->

Here we can see a clear benefit in depletion via a boost in both
UMIs/cell and Genes/cell.

# 2 Setup the seurat object

We are going to use the Seurat toolkit to perform all subsequent
downstream analysis for this tutorial.

We start by reading the data using the Read10x_h5() function which reads
the .h5 file containing the Seurat inputs generated by CellRanger and
generates a count matrix. We next use the count matrix to create a
seurat object.

``` r
library(Seurat)
#control
pbmc_control.mtx <- Read10X_h5("~/R/control_filtered_no_mask.h5")
#depleted
pbmc_depleted.mtx <- Read10X_h5("~/R/depleted_filtered_no_mask_rep3.h5")
```

We are only going to include genes that have at least 3 total in all
cells. And we will start by removing cells that are expressing a minimum
of 300 features. Later, we will use data-driven techniques to determine
appropriate thresholds for filtering.

``` r
#control
pbmc_control.so <- CreateSeuratObject(pbmc_control.mtx, min.cells = 3, project = "10X-V3")
#depleted
pbmc_depleted.so <- CreateSeuratObject(pbmc_depleted.mtx, min.cells = 3, project = "CRISPRclean")

#create a list containing the control and depleted conditions
list.so <- list(pbmc_control.so,pbmc_depleted.so)
#rename the indices of the list
names(list.so) <- c("10x","jc")
#remove the separate files to conserve memory
rm(pbmc_control.so,pbmc_depleted.so)
```

# 3 Data pre-processing and QC workflow

The steps below demonstrate an example of data pre-processing and QC for
scRNA-seq data generated with 10x-v3 3’ gene expression combined with
Jumpcode depletion.

## 3.1 Filtering low-quality cells

We are going to use a data-driven technique to evaluate the number of
unique genes and UMIs detected in each cell.

1.  Dead or decaying cells typically have low complexity with a high %
    of UMIs contributed by mitochondrial genes
2.  Empty droplets or low-quality cells typically express few genes and
    have low library complexity.
3.  Doublets exhibit higher UMI and gene counts which affects downstream
    interpretation.
4.  The total number of genes detected corresponds linearly with unique
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
p1 <- ggplot(list.so[[1]]@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + geom_smooth(method="lm") + ggtitle("10x-v3")
p1 <- ggMarginal(p1, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

``` r
p2 <- ggplot(list.so[[1]]@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
p2 <- ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

``` r
plot_grid(plotlist = list(p1,p2), ncol=2, align='h', rel_widths = c(1, 1))
```

    ## Warning: Graphs cannot be horizontally aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

![](jumpcode_Rcode_files/figure-gfm/plotting%20relationship%20between%20detected%20genes%20and%20counts%20control-1.png)<!-- -->

``` r
#depleted 
p1 <- ggplot(list.so[[2]]@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + geom_smooth(method="lm") + ggtitle("CRISPRclean")
p1 <- ggMarginal(p1, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

``` r
p2 <- ggplot(list.so[[2]]@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
p2 <- ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

``` r
plot_grid(plotlist = list(p1,p2), ncol=2, align='h', rel_widths = c(1, 1))
```

    ## Warning: Graphs cannot be horizontally aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

![](jumpcode_Rcode_files/figure-gfm/plotting%20relationship%20between%20detected%20genes%20and%20counts%20depleted-1.png)<!-- -->

### 3.1.2 Quantify % of reads contributed by mitochondrial, ribosomoal, and all jumpcode gene targets

``` r
#percent mito
list.so <- lapply(list.so, FUN = function(x) {
  x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = 'percent.mt')
})

#percent ribo
list.so <- lapply(list.so, FUN = function(x) {
  x <- PercentageFeatureSet(x, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = 'percent.rb')
})

#create regex for target list
targets <- read.delim("~/targets.txt", header = F)
targets <- targets$V1
targets <- paste0("^", targets, "$", collapse = "|")

#percent all targets
list.so <- lapply(list.so, FUN = function(x) {
  x <- PercentageFeatureSet(x, pattern = targets, col.name = 'percent.all')
})
```

``` r
#here we can visualize the depletion for mitochondrial and ribosomal genes
temp=rbind(
  list.so[[1]]@meta.data %>% select(percent.mt, percent.rb, percent.all, orig.ident),
  list.so[[2]]@meta.data %>% select(percent.mt, percent.rb, percent.all, orig.ident)
  )
ggplot(temp, aes(x = orig.ident, y=percent.mt)) + geom_boxplot() + ggplot(temp, aes(x = orig.ident, y=percent.rb)) + geom_boxplot() + ggplot(temp, aes(x = orig.ident, y=percent.all)) + geom_boxplot() & theme_classic()
```

![](jumpcode_Rcode_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

we can see that we removed a large fraction of targets reads with
CRISPRclean

``` r
temp=rbind(
  list.so[[1]]@meta.data %>% select(nCount_RNA, nFeature_RNA, orig.ident),
  list.so[[2]]@meta.data %>% select(nCount_RNA, nFeature_RNA, orig.ident)
  )
p1 <- ggplot(temp, aes(x = orig.ident, y=nCount_RNA)) + geom_boxplot() + theme_classic()
p2 <- ggplot(temp, aes(x = orig.ident, y=nFeature_RNA)) + geom_boxplot() + theme_classic()
p3 <- ggplot(temp, aes(x = nFeature_RNA/nCount_RNA, fill=orig.ident)) + geom_density() + theme_classic()
(p1 | p2) / 
  p3
```

![](jumpcode_Rcode_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

we can see that we have an increase in library complexity as well

### 3.1.3 Dead cell removal

Cells with high fraction of mitochondrial content are typically of lower
quality and should be filtered out. Keep in mind the cell type as
fraction of mitochondrial reads can vary between cell types. Here We
will keep cells within the 99th quantile and remove outlier cells from
the distribution.

``` r
#control sample
p1 <- ggplot(list.so[[1]]@meta.data, aes(x=nFeature_RNA, y=percent.mt)) +
      geom_point(aes(color="red")) +
      geom_hline(aes(yintercept = quantile(list.so[[1]]$percent.mt, probs = 0.99)), colour = "blue", linetype = 2)
ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100)
```

![](jumpcode_Rcode_files/figure-gfm/plot%20dead%20cell%20threshold%20control-1.png)<!-- -->

If you feel the threshold is too conservative, you can always change
using the probs argument to change the cutoff of the distribution

``` r
list.so[[1]]@meta.data <- list.so[[1]]@meta.data %>% mutate(miQC.keep = ifelse(test = list.so[[1]]@meta.data$percent.mt <= quantile(list.so[[1]]$percent.mt, probs = 0.99), yes = 'keep', no = 'discard'))

FeatureScatter(list.so[[1]], feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = "miQC.keep", jitter = T)
```

![](jumpcode_Rcode_files/figure-gfm/dead%20cell%20removal%20control-1.png)<!-- -->

``` r
#depleted sample
p3 <- ggplot(list.so[[2]]@meta.data, aes(x=nFeature_RNA, y=percent.mt)) +
      geom_point(aes(color="red")) +
      geom_hline(aes(yintercept = quantile(list.so[[2]]$percent.mt, probs = 0.99)), colour = "blue", linetype = 2)
ggMarginal(p3, type = "histogram", fill="lightgrey", bins=100)
```

![](jumpcode_Rcode_files/figure-gfm/dead%20cell%20removal-1.png)<!-- -->

``` r
list.so[[2]]@meta.data <- list.so[[2]]@meta.data %>% mutate(miQC.keep = ifelse(test = list.so[[2]]@meta.data$percent.mt <= quantile(list.so[[2]]$percent.mt, probs=0.99), yes = 'keep', no = 'discard'))

FeatureScatter(list.so[[2]], feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = "miQC.keep", jitter = T)
```

![](jumpcode_Rcode_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#filter the cells
list.so[[1]] <- list.so[[1]][, list.so[[1]]@meta.data[, "miQC.keep"] == 'keep']
list.so[[2]] <- list.so[[2]][, list.so[[2]]@meta.data[, "miQC.keep"] == 'keep']
```

### 3.1.4 Remove cells with too few genes detected

We are going to work with log10 transformed data because it better
preserves the linear relationship as an effect of sequencing.

Using a data-driven technique, we are going use the median absolute
deviation (mad) as cutoff points. Typically, this can range from 3-5 mad
from the median as a cutoff. Make sure to test different options with
your data to determine an appropriate cutoff based on distributions.

For this dataset, we are going to removed cells with significantly low
amount of features (or genes) detected. These could represent empty
droplets and low-quality cells

The goal here is to get a representative normal distribution of the data

``` r
# Gene/UMI scatter plot before filtering for control
p1 <- ggplot(list.so[[1]]@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
  geom_point(aes(color="red")) +
  geom_smooth(method="lm") + 
  geom_hline(aes(yintercept = median(log10(list.so[[1]]$nFeature_RNA)) - 5*mad(log10(list.so[[1]]$nFeature_RNA))), colour = "green", linetype = 2)

ggMarginal(p1, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](jumpcode_Rcode_files/figure-gfm/low%20quality%20cell%20removal-1.png)<!-- -->

``` r
# Gene/UMI scatter plot before filtering for control
p3 <- ggplot(list.so[[2]]@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
  geom_point(aes(color="red")) +
  geom_smooth(method="lm") +
  geom_hline(aes(yintercept = median(log10(list.so[[2]]$nFeature_RNA)) - 5*mad(log10(list.so[[2]]$nFeature_RNA))), colour = "green", linetype = 2)

ggMarginal(p3, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](jumpcode_Rcode_files/figure-gfm/low%20quality%20cell%20removal%20depleted-1.png)<!-- -->

``` r
#filter the cells from the control

min.gene.thresh <- median(log10(list.so[[1]]$nFeature_RNA)) - 5*mad(log10(list.so[[1]]$nFeature_RNA))

cells.keep <- rownames(list.so[[1]]@meta.data %>% filter(log10(nFeature_RNA) > min.gene.thresh))
```

``` r
#control
list.so[[1]] <- subset(list.so[[1]], cells = cells.keep)
```

``` r
#filter the cells from the depleted

min.gene.thresh.d <- median(log10(list.so[[2]]$nFeature_RNA)) - 5*mad(log10(list.so[[2]]$nFeature_RNA))

cells.keep.d <- rownames(list.so[[2]]@meta.data %>% filter(log10(nFeature_RNA) > min.gene.thresh.d))
```

``` r
#depleted
list.so[[2]] <- subset(list.so[[2]], cells = cells.keep.d)
```

### 3.1.5 Filtering low-complexity cells

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
exclude cells that have residuals with \>40% variance below the linear
regression to exclude low complexity cells. Keep in mind that this
threshold can change depending on the cell-type and experiment.

``` r
#control
lm.model = lm(data = list.so[[1]]@meta.data, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
list.so[[1]]@meta.data$residuals <- residuals(lm.model)
list.so[[1]]@meta.data <- list.so[[1]]@meta.data %>% mutate(complexity = ifelse(test = list.so[[1]]@meta.data$residuals >= -0.4, yes = "high" , no = "low"))

p2 <- ggplot(list.so[[1]]@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point(aes(color = complexity)) + geom_abline(intercept = lm.model$coefficients[1] - 0.4 , slope = lm.model$coefficients[2], color="orange", linetype=2) + geom_smooth(method="lm")
ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](jumpcode_Rcode_files/figure-gfm/plotting%20low%20complexity%20cells%20control-1.png)<!-- -->

``` r
#control
lm.model = lm(data = list.so[[2]]@meta.data, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
list.so[[2]]@meta.data$residuals <- residuals(lm.model)
list.so[[2]]@meta.data <- list.so[[2]]@meta.data %>% mutate(complexity = ifelse(test = list.so[[2]]@meta.data$residuals >= -0.4, yes = "high" , no = "low"))

p2 <- ggplot(list.so[[2]]@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point(aes(color = complexity)) + geom_abline(intercept = lm.model$coefficients[1] - 0.4 , slope = lm.model$coefficients[2], color="orange", linetype=2) + geom_smooth(method="lm")
ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](jumpcode_Rcode_files/figure-gfm/plotting%20low%20complexity%20cells%20depleted-1.png)<!-- -->

``` r
#filter the cells from the control
list.so[[1]] <- list.so[[1]][, list.so[[1]]@meta.data[, "complexity"] == 'high']

p2 <- ggplot(list.so[[1]]@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
ggMarginal(p2, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](jumpcode_Rcode_files/figure-gfm/filter%20cells%20from%20control%20and%20view-1.png)<!-- -->

``` r
#filter the cells from the depleted
list.so[[2]] <- list.so[[2]][, list.so$jc@meta.data[, "complexity"] == 'high']

p4 <- ggplot(list.so[[2]]@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
ggMarginal(p4, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](jumpcode_Rcode_files/figure-gfm/filter%20cells%20from%20depleted%20and%20view-1.png)<!-- -->

Next, we are going to use the doublet removal toolkit scDblFinder
(Germain et. al., 2022) to remove doublets from the data. This tool uses
a machine-learning algorithm to simulate artificial doublets from the
data based on cell clustering. From there, real cells get assigned a
doublet score probability to which we will filter cells called as
doublets.

For clustering, we are going to use the Seurat SCTransform workflow.

Here, we want to perform clustering using residual default cutoff of 1.3
(default) rather than selecting a fixed number of highly variable genes.

We will also be demonstrating how to score cell cycle related genes

## 3.2 Doublet filtering

### 3.2.1 clustering using SCTransform workflow

``` r
#SCTransform
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- SCTransform(x, verbose = T, vars.to.regress = c("percent.mt","percent.rb"), variable.features.n = NULL, variable.features.rv.th = 1.3)
})
```

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 19948 by 11193

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 5000 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## Found 74 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 19948 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   2%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  12%  |                                                                              |==========                                                            |  15%  |                                                                              |============                                                          |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  32%  |                                                                              |========================                                              |  35%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  52%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  85%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  98%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 19948 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   2%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  12%  |                                                                              |==========                                                            |  15%  |                                                                              |============                                                          |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  32%  |                                                                              |========================                                              |  35%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  52%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  85%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  98%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 2.406948 mins

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Regressing out percent.mt, percent.rb

    ## Centering data matrix

    ## Set default assay to SCT

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 19968 by 11325

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 5000 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## Found 101 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 19968 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   2%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  12%  |                                                                              |==========                                                            |  15%  |                                                                              |============                                                          |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  32%  |                                                                              |========================                                              |  35%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  52%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  85%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  98%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 19968 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   2%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  12%  |                                                                              |==========                                                            |  15%  |                                                                              |============                                                          |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  32%  |                                                                              |========================                                              |  35%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  52%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  85%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  98%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 3.309266 mins

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Regressing out percent.mt, percent.rb

    ## Centering data matrix

    ## Set default assay to SCT

``` r
#cell cyclce scoring
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- CellCycleScoring(x, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
})

#PCA
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- RunPCA(x, assay = "SCT")
})
```

    ## PC_ 1 
    ## Positive:  LYZ, S100A9, S100A8, CTSS, VCAN, FCN1, CST3, HLA-DRA, MNDA, NEAT1 
    ##     FTL, AIF1, IFI30, S100A12, FOS, PSAP, LST1, CD74, CD14, TYROBP 
    ##     S100A6, CYBB, FGL2, COTL1, S100A11, MS4A6A, HLA-DPA1, FTH1, S100A4, SAT1 
    ## Negative:  GNLY, NKG7, CCL5, GZMB, GZMA, FGFBP2, KLRD1, CST7, PRF1, HOPX 
    ##     CLIC3, PPBP, KLRB1, PF4, GNG11, GP1BB, CAVIN2, TUBB1, HIST1H2AC, SPON2 
    ##     KLRF1, CCL4, CTSW, TSC22D1, GP9, PTCRA, TRDC, ACRBP, IL32, CMTM5 
    ## PC_ 2 
    ## Positive:  GNLY, NKG7, GZMB, GZMA, FGFBP2, KLRD1, PRF1, CST7, HOPX, CLIC3 
    ##     SPON2, KLRF1, KLRB1, CCL4, CTSW, TRDC, GZMH, CCL5, IL2RB, KLRK1 
    ##     CD160, CMC1, GZMM, CD7, FCGR3A, XCL2, IFITM2, AKR1C3, MATK, CD247 
    ## Negative:  PPBP, NRGN, PF4, GP1BB, TUBB1, CAVIN2, GNG11, HIST1H2AC, RGS18, GP9 
    ##     ACRBP, PTCRA, TSC22D1, HLA-DRA, PRKAR2B, CMTM5, TMEM40, ODC1, C2orf88, CLDN5 
    ##     MMD, RUFY1, MAP3K7CL, F13A1, LGALSL, PDLIM1, LIMS1, AC127502.2, MTURN, CLU 
    ## PC_ 3 
    ## Positive:  HLA-DRA, CD74, IGKC, IGHM, MS4A1, CD79A, HLA-DPB1, HLA-DPA1, HLA-DQA1, TCL1A 
    ##     IGHD, HLA-DQB1, BANK1, CD79B, HLA-DRB1, RALGPS2, IGLC2, JCHAIN, TCF4, LINC02397 
    ##     TNFRSF13C, BCL11A, LINC00926, FCRL1, SPIB, FCER2, IGLC3, HLA-DRB5, VPREB3, CD22 
    ## Negative:  S100A8, S100A9, VCAN, LYZ, S100A12, IL7R, FCN1, FOS, MNDA, CTSS 
    ##     NRGN, NEAT1, PPBP, PF4, CCL5, GP1BB, CAVIN2, GNG11, CD14, TUBB1 
    ##     RGS18, S100A4, S100A6, GP9, ACRBP, RGS10, F13A1, THBS1, TSC22D1, CMTM5 
    ## PC_ 4 
    ## Positive:  IL32, GZMK, IL7R, LTB, TRAC, CD3D, CD3G, TRBC2, BCL11B, CD52 
    ##     CD3E, MALAT1, TRBC1, ITGB1, LIME1, RTKN2, PTTG1, CD27, AQP3, CD2 
    ##     CRIP1, MAF, PCLAF, SYNE2, GPR183, CD28, MKI67, CCDC167, LDHB, CD84 
    ## Negative:  GNLY, NKG7, GZMB, FGFBP2, KLRD1, LYZ, CCL5, HLA-DRA, PRF1, CLIC3 
    ##     S100A9, HOPX, TYROBP, GZMA, KLRF1, SPON2, FCER1G, S100A8, CST7, CST3 
    ##     TRDC, CD74, CTSW, CCL4, FCN1, FCGR3A, CD160, KLRB1, CD7, XCL2 
    ## PC_ 5 
    ## Positive:  FCGR3A, CDKN1C, LST1, AIF1, CST3, MS4A7, HLA-DPA1, IFITM3, SAT1, IFI30 
    ##     FCER1G, FTL, COTL1, TCF7L2, SMIM25, C1QA, FTH1, HLA-DPB1, LYPD2, CSF1R 
    ##     RHOC, CALHM6, HLA-DRB1, HLA-DRA, TMSB4X, CEBPB, HMOX1, HES4, SERPINA1, MTSS1 
    ## Negative:  S100A8, VCAN, S100A9, S100A12, FOS, IGHM, LYZ, IGKC, MS4A1, CD79A 
    ##     MNDA, IGHD, TCL1A, CD14, CYP1B1, THBS1, BANK1, RALGPS2, NCF1, CSF3R 
    ##     IGLC2, CD36, PLBD1, MEGF9, LINC02397, SLC2A3, FCRL1, LINC00926, TNFRSF13C, RGS2

    ## PC_ 1 
    ## Positive:  LYZ, S100A9, VCAN, HLA-DRA, CD74, MNDA, AIF1, S100A12, IFI30, CD14 
    ##     CTSS, CYBB, HLA-DPA1, S100A8, HLA-DRB1, FCN1, COTL1, FGL2, MS4A6A, HLA-DPB1 
    ##     FTH1, GRN, TKT, FOS, FTL, MPEG1, NCF2, TYMP, CST3, VIM 
    ## Negative:  GNLY, NKG7, GZMA, GZMB, CCL5, FGFBP2, CST7, PRF1, KLRB1, HOPX 
    ##     KLRD1, CLIC3, KLRF1, CTSW, SPON2, TRDC, CCL4, KLRK1, GZMH, GZMM 
    ##     CMC1, GZMK, CD7, IL2RB, CD160, ARL4C, MATK, CD247, AKR1C3, IL32 
    ## PC_ 2 
    ## Positive:  NRGN, GP1BB, PF4, TUBB1, CAVIN2, GNG11, HIST1H2AC, GP9, RGS18, ACRBP 
    ##     PTCRA, TSC22D1, CMTM5, TMEM40, PRKAR2B, PPBP, C2orf88, MMD, CLDN5, MAP3K7CL 
    ##     CLU, ODC1, LGALSL, SPARC, F13A1, RUFY1, MTURN, TREML1, CLEC1B, AC147651.1 
    ## Negative:  GNLY, NKG7, GZMB, GZMA, FGFBP2, PRF1, CST7, HOPX, KLRD1, CLIC3 
    ##     KLRB1, KLRF1, SPON2, CTSW, TRDC, CCL4, FCGR3A, GZMH, S100A9, KLRK1 
    ##     LYZ, CMC1, GZMM, CD160, IL2RB, IFITM2, MATK, CD7, HCST, AKR1C3 
    ## PC_ 3 
    ## Positive:  HLA-DRA, CD74, MS4A1, CD79A, IGHM, IGKC, HLA-DPB1, HLA-DPA1, TCL1A, HLA-DQA1 
    ##     IGHD, HLA-DQB1, BANK1, CD79B, FCRL1, RALGPS2, LINC02397, HLA-DRB1, IGLC2, TNFRSF13C 
    ##     FCER2, LINC00926, VPREB3, IGLC3, BCL11A, CD22, AFF3, JCHAIN, TCF4, CD37 
    ## Negative:  S100A9, LYZ, VCAN, S100A12, MNDA, S100A8, IL7R, CD14, AIF1, FOS 
    ##     VIM, S100A6, FCN1, NRGN, THBS1, FYB1, CTSS, CEBPD, RGS2, RGS18 
    ##     FTH1, ANXA1, TKT, S100A11, CYP1B1, GIMAP7, PLBD1, PF4, MS4A6A, FGL2 
    ## PC_ 4 
    ## Positive:  GNLY, NKG7, GZMB, LYZ, FGFBP2, S100A9, GZMA, PRF1, CLIC3, CCL5 
    ##     HOPX, KLRD1, HLA-DRA, SPON2, KLRF1, CST7, VCAN, FCGR3A, CTSW, CD74 
    ##     TRDC, CCL4, IFI30, AIF1, NRGN, MNDA, CD160, GZMH, AKR1C3, EFHD2 
    ## Negative:  IL7R, LTB, TRAC, CD52, CD3D, BCL11B, CD3E, CD3G, IL32, TRBC1 
    ##     GZMK, CRIP1, MALAT1, ITGB1, LIME1, INPP4B, CD27, RCAN3, MAL, TRAT1 
    ##     AQP3, CD2, BCL2, GPR183, RTKN2, LEF1, ETS1, BIRC3, PBXIP1, SPOCK2 
    ## PC_ 5 
    ## Positive:  MS4A1, CD79A, IGHM, IGHD, VCAN, S100A9, S100A12, FCRL1, BANK1, CD79B 
    ##     TCL1A, IGLC2, FCER2, S100A8, RALGPS2, LINC02397, LINC00926, TNFRSF13C, VPREB3, IGLC3 
    ##     CD22, MNDA, PAX5, LYZ, GZMK, HLA-DOB, CD37, CD24, CD72, NKG7 
    ## Negative:  PLD4, SERPINF1, LILRA4, GZMB, IL3RA, ITM2C, CCDC50, UGCG, JCHAIN, IRF8 
    ##     TPM2, IRF7, PPP1R14B, LRRC26, DNASE1L3, SCT, CLEC4C, TCF4, MZB1, DERL3 
    ##     LINC00996, SMPD3, MAP1A, SCN9A, EPHB1, BCL11A, SCAMP5, FCGR3A, AC023590.1, PACSIN1

``` r
#generate UMAP coordinates
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- RunUMAP(x, dims = 1:30, n.neighbors=20)
})
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 12:30:14 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 12:30:14 Read 11193 rows and found 30 numeric columns

    ## 12:30:14 Using Annoy for neighbor search, n_neighbors = 20

    ## 12:30:14 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 12:30:15 Writing NN index file to temp file C:\Users\dante\AppData\Local\Temp\Rtmp82HhsV\file3a9c20c61ce
    ## 12:30:16 Searching Annoy index using 1 thread, search_k = 2000
    ## 12:30:18 Annoy recall = 100%
    ## 12:30:19 Commencing smooth kNN distance calibration using 1 thread
    ## 12:30:20 Initializing from normalized Laplacian + noise
    ## 12:30:22 Commencing optimization for 200 epochs, with 311714 positive edges
    ## 12:30:36 Optimization finished
    ## 12:30:36 UMAP embedding parameters a = 0.9922 b = 1.112
    ## 12:30:36 Read 11325 rows and found 30 numeric columns
    ## 12:30:36 Using Annoy for neighbor search, n_neighbors = 20
    ## 12:30:36 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 12:30:38 Writing NN index file to temp file C:\Users\dante\AppData\Local\Temp\Rtmp82HhsV\file3a9c7f445fe7
    ## 12:30:38 Searching Annoy index using 1 thread, search_k = 2000
    ## 12:30:41 Annoy recall = 100%
    ## 12:30:41 Commencing smooth kNN distance calibration using 1 thread
    ## 12:30:42 Initializing from normalized Laplacian + noise
    ## 12:30:43 Commencing optimization for 200 epochs, with 318216 positive edges
    ## 12:30:55 Optimization finished

``` r
#find k-nearest neighbors
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- FindNeighbors(x, dims = 1:30, k.param=20)
})
```

    ## Computing nearest neighbor graph
    ## Computing SNN
    ## Computing nearest neighbor graph
    ## Computing SNN

``` r
#Find clusters using the louvain algorithm with multilevel refinement. It is recommended to overcluster the data first when using scDblFinder
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- FindClusters(x, resolution = 2, algorithm = 2)
})
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 11193
    ## Number of edges: 387484
    ## 
    ## Running Louvain algorithm with multilevel refinement...
    ## Maximum modularity in 10 random starts: 0.8121
    ## Number of communities: 37
    ## Elapsed time: 2 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 11325
    ## Number of edges: 395736
    ## 
    ## Running Louvain algorithm with multilevel refinement...
    ## Maximum modularity in 10 random starts: 0.8203
    ## Number of communities: 33
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
control <- as.SingleCellExperiment(list.so[[1]], assay = 'RNA')

jc <- as.SingleCellExperiment(list.so[[2]], assay = 'RNA')

#convert to a list
list.dbr <- list(control,jc)

rm(control,jc)

d <- lapply(list.dbr, FUN = function(x) {
  x <- scDblFinder(x, clusters = 'seurat_clusters', dbr = NULL, dims = 30, includePCs = 30, returnType = "table", k = 10, processing = "normFeatures")
})
```

    ## 37 clusters

    ## Creating ~13690 artificial doublets...

    ## Dimensional reduction

    ## Evaluating kNN...

    ## Training model...

    ## iter=0, 630 cells excluded from training.

    ## iter=1, 625 cells excluded from training.

    ## iter=2, 627 cells excluded from training.

    ## Threshold found:0.601

    ## 687 (6.1%) doublets called

    ## 33 clusters

    ## Creating ~10890 artificial doublets...

    ## Dimensional reduction

    ## Evaluating kNN...

    ## Training model...

    ## iter=0, 549 cells excluded from training.

    ## iter=1, 684 cells excluded from training.

    ## iter=2, 755 cells excluded from training.

    ## Threshold found:0.541

    ## 814 (7.2%) doublets called

For cells called doublets, we should see twice the amount of UMI counts
and more genes/cell. We have view this with Violin plots to check the
doublet calls.

``` r
#plotting control doublet calls
list.so[[1]]$class <- as.data.frame(d[[1]]) %>% filter(type == "real") %>% select(class)
list.so[[1]]$class <- as.factor(list.so[[1]]$class)
list.so[[1]]$class <- factor(list.so[[1]]$class, levels = c("singlet","doublet"))

DimPlot(list.so[[1]], group.by = c("class"), order = T)
```

![](jumpcode_Rcode_files/figure-gfm/doublet%20calls%20control-1.png)<!-- -->

``` r
VlnPlot(list.so[[1]], features = c("nCount_RNA", "nFeature_RNA"), group.by = "class")
```

![](jumpcode_Rcode_files/figure-gfm/doublet%20calls%20control-2.png)<!-- -->

``` r
FeatureScatter(list.so[[1]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "class")
```

![](jumpcode_Rcode_files/figure-gfm/doublet%20calls%20control-3.png)<!-- -->

``` r
#plotting depleted doublet calls
list.so[[2]]$class <- as.data.frame(d[[2]]) %>% filter(type == "real") %>% select(class)
list.so[[2]]$class <- as.factor(list.so[[2]]$class)
list.so[[2]]$class <- factor(list.so[[2]]$class, levels = c("singlet","doublet"))

DimPlot(list.so[[2]], group.by = c("class"), order = T)
```

![](jumpcode_Rcode_files/figure-gfm/doublet%20calls%20depleted-1.png)<!-- -->

``` r
VlnPlot(list.so[[2]], features = c("nCount_RNA", "nFeature_RNA"), group.by = "class")
```

![](jumpcode_Rcode_files/figure-gfm/doublet%20calls%20depleted-2.png)<!-- -->

``` r
FeatureScatter(list.so[[2]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "class")
```

![](jumpcode_Rcode_files/figure-gfm/doublet%20calls%20depleted-3.png)<!-- -->

The doublet call information looks like what we would expect. Now we
will filter these cells from each sample.

``` r
#remove doublets control
list.so[[1]] <- list.so[[1]][, list.so[[1]]@meta.data[, "class"] == "singlet"]

#remove doublets depleted
list.so[[2]] <- list.so[[2]][, list.so[[2]]@meta.data[, "class"] == "singlet"]
```

Even after doublet filtering, we can see that cells with high UMI counts
are retained. These could be sources of homotypic doublets which doublet
filtering algorithms are less efficient at removing. We will remove
these to potentially remove any sources of technical variation.

``` r
# Gene/UMI scatter plot before filtering for control
p1 <- ggplot(list.so[[1]]@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
  geom_point() +
  geom_smooth(method="lm") + 
  geom_vline(aes(xintercept = quantile(log10(list.so[[1]]$nCount_RNA), probs=0.99)), colour = "red", linetype = 2)

ggMarginal(p1, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](jumpcode_Rcode_files/figure-gfm/low%20quality%20cell%20removal%20from%20control-1.png)<!-- -->

``` r
# Gene/UMI scatter plot before filtering for control
p1 <- ggplot(list.so[[2]]@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
  geom_point() +
  geom_smooth(method="lm") + 
  geom_vline(aes(xintercept = quantile(log10(list.so[[2]]$nCount_RNA), probs=0.99)), colour = "red", linetype = 2)

ggMarginal(p1, type = "histogram", fill="lightgrey")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](jumpcode_Rcode_files/figure-gfm/low%20quality%20cell%20removal%20from%20depleted-1.png)<!-- -->

``` r
#filter the cells from the control

max.umi.thresh <- quantile(log10(list.so[[1]]$nCount_RNA), probs=0.99)[[1]]

cells.keep <- rownames(list.so[[1]]@meta.data %>% filter(log10(nCount_RNA) < max.umi.thresh))
```

``` r
#control
list.so[[1]] <- subset(list.so[[1]], cells = cells.keep)
```

``` r
#filter the cells from the depleted

max.umi.thresh.d <- quantile(log10(list.so[[2]]$nCount_RNA), probs=0.99)[[1]]

cells.keep.d <- rownames(list.so[[2]]@meta.data %>% filter(log10(nCount_RNA) < max.umi.thresh.d))
```

``` r
#depleted
list.so[[2]] <- subset(list.so[[2]], cells = cells.keep.d)
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
list.so[[1]]$cc.difference <- list.so[[1]]$S.Score - list.so[[1]]$G2M.Score
VlnPlot(list.so[[1]], features = c('S.Score','G2M.Score','cc.difference'), group.by = 'orig.ident')
```

![](jumpcode_Rcode_files/figure-gfm/cell%20cycling%20control-1.png)<!-- -->

``` r
#we quantified cells is S and G2 phase earlier, so we will look to see the proportion of cycling cells
#depleted
list.so[[2]]$cc.difference <- list.so[[2]]$S.Score - list.so[[2]]$G2M.Score
VlnPlot(list.so[[2]], features = c('S.Score','G2M.Score','cc.difference'), group.by = 'orig.ident')
```

![](jumpcode_Rcode_files/figure-gfm/cell%20cycling%20depleted-1.png)<!-- -->

We can see there is little influence of these PBMC cells. We can always
choose to regress this influence out to our choosing when selecting
vars.to.regress during SCTransform. For now, we won’t regress cell cycle
genes from the data because there could be some interesting variation
(i.e., proliferating cell populations).

# 4 Cell clustering using SCTransform workflow

Now that we have done our cell QC filtering, we will repeat the
SCTransform workflow to produce clustering results without the influence
of doublets

``` r
#SCTransform
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- SCTransform(x, verbose = T, vars.to.regress = c("percent.mt","percent.rb"), variable.features.n = NULL, variable.features.rv.th = 1.3)
})
```

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 19439 by 10400

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 5000 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## Found 75 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 19439 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  21%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  38%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  49%  |                                                                              |====================================                                  |  51%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  62%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  69%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  79%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  85%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 19439 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  21%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  38%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  49%  |                                                                              |====================================                                  |  51%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  62%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  69%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  79%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  85%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 3.128791 mins

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Regressing out percent.mt, percent.rb

    ## Centering data matrix

    ## Set default assay to SCT

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 19452 by 10405

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 5000 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## Found 91 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 19452 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  21%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  38%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  49%  |                                                                              |====================================                                  |  51%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  62%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  69%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  79%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  85%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 19452 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  21%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  38%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  49%  |                                                                              |====================================                                  |  51%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  62%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  69%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  79%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  85%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  92%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 3.170854 mins

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Regressing out percent.mt, percent.rb

    ## Centering data matrix

    ## Set default assay to SCT

``` r
#PCA
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- RunPCA(x, assay = "SCT")
})
```

    ## PC_ 1 
    ## Positive:  GNLY, NKG7, CCL5, GZMB, GZMA, KLRD1, FGFBP2, CST7, PRF1, HOPX 
    ##     KLRB1, CLIC3, PPBP, KLRF1, PF4, SPON2, GNG11, GP1BB, CAVIN2, CTSW 
    ##     TUBB1, HIST1H2AC, CCL4, TRDC, IL32, TSC22D1, GP9, PTCRA, ACRBP, GZMM 
    ## Negative:  LYZ, S100A9, S100A8, CTSS, VCAN, FCN1, CST3, HLA-DRA, NEAT1, MNDA 
    ##     FTL, AIF1, S100A12, IFI30, FOS, CD74, PSAP, LST1, S100A6, CD14 
    ##     TYROBP, CYBB, S100A11, FGL2, COTL1, MS4A6A, FTH1, S100A4, SAT1, HLA-DPA1 
    ## PC_ 2 
    ## Positive:  GNLY, NKG7, GZMB, GZMA, FGFBP2, KLRD1, PRF1, CST7, HOPX, CLIC3 
    ##     KLRB1, KLRF1, SPON2, CTSW, CCL4, TRDC, GZMH, IL2RB, CMC1, KLRK1 
    ##     CCL5, GZMM, CD160, FCGR3A, CD7, IFITM2, XCL2, AKR1C3, CD247, ID2 
    ## Negative:  PPBP, NRGN, PF4, GP1BB, TUBB1, CAVIN2, GNG11, HIST1H2AC, RGS18, ACRBP 
    ##     GP9, PTCRA, HLA-DRA, TSC22D1, PRKAR2B, CMTM5, TMEM40, ODC1, C2orf88, CLDN5 
    ##     MMD, RUFY1, MAP3K7CL, LIMS1, F13A1, PDLIM1, LGALSL, RGS10, CD74, AC127502.2 
    ## PC_ 3 
    ## Positive:  HLA-DRA, CD74, IGHM, IGKC, MS4A1, CD79A, HLA-DPB1, HLA-DPA1, HLA-DQA1, TCL1A 
    ##     IGHD, HLA-DQB1, BANK1, CD79B, HLA-DRB1, RALGPS2, IGLC2, TNFRSF13C, LINC02397, BCL11A 
    ##     TCF4, LINC00926, FCRL1, SPIB, JCHAIN, FCER2, VPREB3, HLA-DRB5, CD22, IGLC3 
    ## Negative:  S100A8, S100A9, VCAN, LYZ, S100A12, IL7R, FCN1, FOS, MNDA, CTSS 
    ##     NRGN, NEAT1, S100A4, PPBP, S100A6, PF4, GP1BB, GNG11, CAVIN2, CD14 
    ##     RGS18, TUBB1, RGS10, CCL5, LIMS1, GZMK, ACRBP, GP9, F13A1, IL32 
    ## PC_ 4 
    ## Positive:  GNLY, NKG7, GZMB, FGFBP2, LYZ, KLRD1, S100A9, CCL5, S100A8, PRF1 
    ##     CLIC3, HOPX, TYROBP, GZMA, KLRF1, SPON2, CST7, VCAN, HLA-DRA, FCER1G 
    ##     TRDC, CTSW, FCN1, CTSS, CCL4, CST3, CD160, KLRB1, CD7, CD74 
    ## Negative:  IL32, GZMK, IL7R, LTB, TRAC, MALAT1, CD3D, CD52, TRBC2, CD3G 
    ##     BCL11B, CD3E, ITGB1, TRBC1, RTKN2, CD2, LIME1, SYNE2, MAF, GPR183 
    ##     CD27, CRIP1, AQP3, CD84, TIGIT, INPP4B, CDC14A, BIRC3, ETS1, BCL2 
    ## PC_ 5 
    ## Positive:  FCGR3A, CDKN1C, LST1, AIF1, CST3, MS4A7, IFITM3, SAT1, HLA-DPA1, FCER1G 
    ##     FTL, IFI30, COTL1, TCF7L2, FTH1, SMIM25, C1QA, HLA-DPB1, CSF1R, LYPD2 
    ##     HLA-DRA, CALHM6, RHOC, IFITM2, HES4, SERPINA1, TMSB4X, HMOX1, CEBPB, HLA-DRB1 
    ## Negative:  S100A8, VCAN, S100A12, S100A9, FOS, IGHM, IGKC, LYZ, MS4A1, CD79A 
    ##     MNDA, IGHD, TCL1A, CD14, CYP1B1, THBS1, NCF1, BANK1, CSF3R, CD36 
    ##     PLBD1, MEGF9, RALGPS2, SLC2A3, IGLC2, RGS2, LINC02397, VNN2, TNFRSF13C, FCRL1

    ## PC_ 1 
    ## Positive:  NKG7, GNLY, GZMA, GZMB, CCL5, FGFBP2, CST7, PRF1, KLRB1, HOPX 
    ##     KLRD1, CLIC3, KLRF1, CTSW, SPON2, TRDC, CCL4, KLRK1, GZMM, GZMH 
    ##     CMC1, GZMK, CD7, IL2RB, ARL4C, CD160, MATK, CD247, IL32, AKR1C3 
    ## Negative:  LYZ, S100A9, VCAN, HLA-DRA, MNDA, AIF1, IFI30, CD74, S100A12, CD14 
    ##     CTSS, CYBB, FCN1, S100A8, HLA-DPA1, HLA-DRB1, COTL1, FGL2, FTH1, MS4A6A 
    ##     GRN, TKT, FOS, MPEG1, TYMP, HLA-DPB1, FTL, NCF2, VIM, CST3 
    ## PC_ 2 
    ## Positive:  NKG7, GNLY, GZMB, GZMA, FGFBP2, PRF1, CST7, HOPX, KLRD1, CLIC3 
    ##     KLRB1, KLRF1, SPON2, CTSW, TRDC, S100A9, LYZ, CCL4, FCGR3A, GZMH 
    ##     KLRK1, CMC1, GZMM, VCAN, CD160, IFITM2, IL2RB, MATK, HCST, CD7 
    ## Negative:  NRGN, GP1BB, PF4, CAVIN2, TUBB1, GNG11, HIST1H2AC, GP9, ACRBP, RGS18 
    ##     PTCRA, TSC22D1, CMTM5, TMEM40, PRKAR2B, PPBP, C2orf88, MMD, ODC1, CLDN5 
    ##     CLU, MAP3K7CL, LGALSL, SPARC, F13A1, RUFY1, MTURN, TREML1, PDLIM1, AC147651.1 
    ## PC_ 3 
    ## Positive:  S100A9, LYZ, VCAN, S100A12, MNDA, IL7R, S100A8, CD14, AIF1, VIM 
    ##     NRGN, FOS, FTH1, RGS18, FCN1, S100A6, FYB1, PF4, CEBPD, CTSS 
    ##     ANXA1, RGS2, TUBB1, THBS1, TKT, GP1BB, CAVIN2, GNG11, GIMAP7, S100A11 
    ## Negative:  HLA-DRA, CD74, MS4A1, CD79A, IGHM, IGKC, HLA-DPB1, HLA-DPA1, HLA-DQA1, HLA-DQB1 
    ##     IGHD, TCL1A, BANK1, CD79B, FCRL1, HLA-DRB1, LINC02397, RALGPS2, IGLC2, TNFRSF13C 
    ##     FCER2, LINC00926, VPREB3, BCL11A, CD22, AFF3, IGLC3, CD37, TCF4, PAX5 
    ## PC_ 4 
    ## Positive:  GNLY, NKG7, GZMB, FGFBP2, GZMA, PRF1, CCL5, CLIC3, LYZ, HOPX 
    ##     KLRD1, S100A9, CST7, KLRF1, SPON2, FCGR3A, HLA-DRA, CTSW, TRDC, VCAN 
    ##     CCL4, CD74, IFI30, NRGN, CD160, AIF1, KLRB1, EFHD2, AKR1C3, GZMH 
    ## Negative:  IL7R, LTB, CD52, TRAC, CD3D, BCL11B, CD3E, CD3G, IL32, TRBC1 
    ##     MALAT1, ITGB1, CRIP1, GZMK, INPP4B, RCAN3, LIME1, CD27, MAL, BCL2 
    ##     TRAT1, CD2, GPR183, AQP3, RTKN2, LEF1, ETS1, BIRC3, PBXIP1, FYB1 
    ## PC_ 5 
    ## Positive:  FCGR3A, PLD4, CDKN1C, HLA-DPA1, IL3RA, AIF1, SERPINF1, LILRA4, PPP1R14B, UGCG 
    ##     GZMB, ITM2C, CCDC50, IRF7, TPM2, LRRC26, SCT, CLEC4C, DNASE1L3, MS4A7 
    ##     IRF8, NPC2, JCHAIN, TCF4, SMPD3, IFI30, MZB1, DERL3, CST3, LINC00996 
    ## Negative:  VCAN, S100A12, MS4A1, CD79A, S100A9, IGHM, S100A8, IGHD, FCRL1, LYZ 
    ##     MNDA, BANK1, TCL1A, FCER2, CD79B, IGLC2, LINC02397, RALGPS2, CD14, FOS 
    ##     LINC00926, VPREB3, TNFRSF13C, THBS1, CYP1B1, CD22, IGLC3, PLBD1, PAX5, GNLY

``` r
#generate UMAP coordinates
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- RunUMAP(x, dims = 1:30, n.neighbors=20)
})
```

    ## 12:48:54 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 12:48:54 Read 10400 rows and found 30 numeric columns

    ## 12:48:54 Using Annoy for neighbor search, n_neighbors = 20

    ## 12:48:54 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 12:48:56 Writing NN index file to temp file C:\Users\dante\AppData\Local\Temp\Rtmp82HhsV\file3a9c2ab47f9b
    ## 12:48:56 Searching Annoy index using 1 thread, search_k = 2000
    ## 12:49:00 Annoy recall = 100%
    ## 12:49:02 Commencing smooth kNN distance calibration using 1 thread
    ## 12:49:04 Initializing from normalized Laplacian + noise
    ## 12:49:04 Commencing optimization for 200 epochs, with 287768 positive edges
    ## 12:49:17 Optimization finished
    ## 12:49:17 UMAP embedding parameters a = 0.9922 b = 1.112
    ## 12:49:17 Read 10405 rows and found 30 numeric columns
    ## 12:49:17 Using Annoy for neighbor search, n_neighbors = 20
    ## 12:49:17 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 12:49:19 Writing NN index file to temp file C:\Users\dante\AppData\Local\Temp\Rtmp82HhsV\file3a9c24cd859
    ## 12:49:19 Searching Annoy index using 1 thread, search_k = 2000
    ## 12:49:21 Annoy recall = 100%
    ## 12:49:22 Commencing smooth kNN distance calibration using 1 thread
    ## 12:49:24 Initializing from normalized Laplacian + noise
    ## 12:49:25 Commencing optimization for 200 epochs, with 289580 positive edges
    ## 12:49:37 Optimization finished

``` r
#find k-nearest neighbors. For consistency, are going to use the same number of neighbors used in RunUMAP()
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- FindNeighbors(x, dims = 1:30, k.param=20)
})
```

    ## Computing nearest neighbor graph
    ## Computing SNN
    ## Computing nearest neighbor graph
    ## Computing SNN

``` r
#Find clusters using the louvain algorithm with multilevel refinement. Each clustering resolution will be stored in the metadata
list.so <- lapply(X = list.so, FUN = function(x) {
  x <- FindClusters(x, resolution = seq(0.1,1,0.1), algorithm = 2, verbose=F)
})
```

``` r
#set the Ident of the object to resolution
Idents(list.so[[1]]) <- "SCT_snn_res.0.2"
Idents(list.so[[2]]) <- "SCT_snn_res.0.2"

#visualize UMAP of the control
DimPlot(list.so[[1]], label = T, pt.size = 0.75) + theme_classic() + ggtitle('10x')
```

![](jumpcode_Rcode_files/figure-gfm/UMAP%20plots%20of%20clustering%20results-1.png)<!-- -->

``` r
#visualize UMAP of the depleted
DimPlot(list.so[[2]], label = T, pt.size = 0.75) + theme_classic() + ggtitle('JC')
```

![](jumpcode_Rcode_files/figure-gfm/UMAP%20plots%20of%20clustering%20results-2.png)<!-- -->
