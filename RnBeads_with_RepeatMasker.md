Mouse Infinium methylation data processing with RnBeads
================
Seung June Kim

#### Introduction

The goal is to display some key data frame structures and visuals that
you can expect to get from running the code in
RnBeads_with_RepeatMasker.R

Using RnBeads dev version (Müller et al. 2019).

#### Setup[^1][^2]

What should RepeatMasker table and \$SAMPLE.csv file look like?
`sample.table` was created solely for demonstration purposes; it’s not
needed to run the code.

``` r
head(tab.rmsk)
```

    ##   mean.mean.g1 mean.mean.g2 mean.mean.diff mean.mean.quot.log2 comb.p.val comb.p.adj.fdr
    ## 1    0.7539077    0.7713189   -0.017411198        -0.032513175  0.1091943      0.4557239
    ## 2    0.8585011    0.8533401    0.005160963         0.008598615  0.5222993      0.8408950
    ## 3    0.8657933    0.8744245   -0.008631208        -0.014148592  0.2042860      0.6075644
    ## 4    0.7520166    0.7590930   -0.007076353        -0.013335547  0.6216694      0.8860324
    ## 5    0.8108180    0.8189176   -0.008099642        -0.014166401  0.4723224      0.8144906
    ## 6    0.8954512    0.8973859   -0.001934682        -0.003079324  0.7580065      0.9374626
    ##   num.sites mean.num.na.g1 mean.num.na.g2 mean.mean.covg.g1 mean.mean.covg.g2
    ## 1         1              0              0          22.33333          25.00000
    ## 2         1              0              0          27.66667          20.00000
    ## 3         1              0              0          24.66667          24.33333
    ## 4         1              0              0          23.66667          23.33333
    ## 5         1              0              0          19.00000          19.33333
    ## 6         1              0              0          25.66667          29.33333
    ##   mean.nsamples.covg.thresh.g1 mean.nsamples.covg.thresh.g2 combinedRank
    ## 1                            3                            3         6985
    ## 2                            3                            3        18576
    ## 3                            3                            3        13776
    ## 4                            3                            3        20454
    ## 5                            3                            3        16904
    ## 6                            3                            3        25046

``` r
head(sample.table,14)
```

    ##            X.Header.                              X          X.1          X.2     X.3
    ## 1  Investigator Name                                                                 
    ## 2       Project Name                                                                 
    ## 3    Experiment Name                                                                 
    ## 4               Date                     02-02-2023                                  
    ## 5        [Manifests]                                                                 
    ## 6                  A MouseMethylation-12v1-0_A2.bpm                                  
    ## 7             [Data]                                                                 
    ## 8        Sample_Name                    Sample_Well Sample_Plate Sample_Group Pool_ID
    ## 9            B1 DMSO                            A01     PAS23023                    0
    ## 10          NB1 DMSO                            C01     PAS23023                    0
    ## 11           B2 DMSO                            E01     PAS23023                    0
    ## 12          NB2 DMSO                            G01     PAS23023                    0
    ## 13           B3 DMSO                            A02     PAS23023                    0
    ## 14          NB3 DMSO                            C02     PAS23023                    0
    ##             X.4              X.5         X.6              X.7              X.8
    ## 1                                                                             
    ## 2                                                                             
    ## 3                                                                             
    ## 4                                                                             
    ## 5                                                                             
    ## 6                                                                             
    ## 7                                                                             
    ## 8    Sentrix_ID Sentrix_Position Sample_Type Sample_Treatment Sample_Replicate
    ## 9  205832350006           R01C01           B                D                1
    ## 10 205832350006           R03C01          NB                D                1
    ## 11 205832350006           R05C01           B                D                2
    ## 12 205832350006           R01C02          NB                D                2
    ## 13 205832350006           R03C02           B                D                3
    ## 14 205832350006           R05C02          NB                D                3
    ##                X.9
    ## 1                 
    ## 2                 
    ## 3                 
    ## 4                 
    ## 5                 
    ## 6                 
    ## 7                 
    ## 8  Sample_Identity
    ## 9              B_D
    ## 10            NB_D
    ## 11             B_D
    ## 12            NB_D
    ## 13             B_D
    ## 14            NB_D

#### Interesting but confusing functions

`get.table()` returns a useful, unsorted data frame of statistics for
each \$ANNOTATION probe(rows) from `diffmeth.\$ANNOTATION` S4 object.
`annotation()` returns row-matching annotation for each probe in
tab.\$ANNOTATION.

``` r
diffmeth.rmsk <- rnb.execute.computeDiffMeth(rnb.set, cmp.cols, region.types="rmsk")
tab.rmsk <- get.table(diffmeth.rmsk,comparison, region.type="rmsk",return.data.frame=T)
aa.rmsk <- annotation(rnb.set,type="rmsk")
```

``` r
head(tab.rmsk)
```

    ##   mean.mean.g1 mean.mean.g2 mean.mean.diff mean.mean.quot.log2 comb.p.val comb.p.adj.fdr
    ## 1    0.7539077    0.7713189   -0.017411198        -0.032513175  0.1091943      0.4557239
    ## 2    0.8585011    0.8533401    0.005160963         0.008598615  0.5222993      0.8408950
    ## 3    0.8657933    0.8744245   -0.008631208        -0.014148592  0.2042860      0.6075644
    ## 4    0.7520166    0.7590930   -0.007076353        -0.013335547  0.6216694      0.8860324
    ## 5    0.8108180    0.8189176   -0.008099642        -0.014166401  0.4723224      0.8144906
    ## 6    0.8954512    0.8973859   -0.001934682        -0.003079324  0.7580065      0.9374626
    ##   num.sites mean.num.na.g1 mean.num.na.g2 mean.mean.covg.g1 mean.mean.covg.g2
    ## 1         1              0              0          22.33333          25.00000
    ## 2         1              0              0          27.66667          20.00000
    ## 3         1              0              0          24.66667          24.33333
    ## 4         1              0              0          23.66667          23.33333
    ## 5         1              0              0          19.00000          19.33333
    ## 6         1              0              0          25.66667          29.33333
    ##   mean.nsamples.covg.thresh.g1 mean.nsamples.covg.thresh.g2 combinedRank
    ## 1                            3                            3         6985
    ## 2                            3                            3        18576
    ## 3                            3                            3        13776
    ## 4                            3                            3        20454
    ## 5                            3                            3        16904
    ## 6                            3                            3        25046

``` r
head(aa.rmsk)
```

    ##      Chromosome   Start     End Strand name
    ## 699        chr1 3035552 3035969      *  LTR
    ## 843        chr1 3121629 3121965      *  LTR
    ## 1368       chr1 3469565 3469891      *  LTR
    ## 2176       chr1 3984251 3985340      *  LTR
    ## 2489       chr1 4208633 4209515      * LINE
    ## 3132       chr1 4584673 4585628      * LINE

These two data frames are bound together.

``` r
annotated.tab.rmsk <- data.frame(tab.rmsk,aa.rmsk, row.names=NULL)
head(annotated.tab.rmsk)
```

    ##   mean.mean.g1 mean.mean.g2 mean.mean.diff mean.mean.quot.log2 comb.p.val comb.p.adj.fdr
    ## 1    0.7539077    0.7713189   -0.017411198        -0.032513175  0.1091943      0.4557239
    ## 2    0.8585011    0.8533401    0.005160963         0.008598615  0.5222993      0.8408950
    ## 3    0.8657933    0.8744245   -0.008631208        -0.014148592  0.2042860      0.6075644
    ## 4    0.7520166    0.7590930   -0.007076353        -0.013335547  0.6216694      0.8860324
    ## 5    0.8108180    0.8189176   -0.008099642        -0.014166401  0.4723224      0.8144906
    ## 6    0.8954512    0.8973859   -0.001934682        -0.003079324  0.7580065      0.9374626
    ##   num.sites mean.num.na.g1 mean.num.na.g2 mean.mean.covg.g1 mean.mean.covg.g2
    ## 1         1              0              0          22.33333          25.00000
    ## 2         1              0              0          27.66667          20.00000
    ## 3         1              0              0          24.66667          24.33333
    ## 4         1              0              0          23.66667          23.33333
    ## 5         1              0              0          19.00000          19.33333
    ## 6         1              0              0          25.66667          29.33333
    ##   mean.nsamples.covg.thresh.g1 mean.nsamples.covg.thresh.g2 combinedRank Chromosome
    ## 1                            3                            3         6985       chr1
    ## 2                            3                            3        18576       chr1
    ## 3                            3                            3        13776       chr1
    ## 4                            3                            3        20454       chr1
    ## 5                            3                            3        16904       chr1
    ## 6                            3                            3        25046       chr1
    ##     Start     End Strand name
    ## 1 3035552 3035969      *  LTR
    ## 2 3121629 3121965      *  LTR
    ## 3 3469565 3469891      *  LTR
    ## 4 3984251 3985340      *  LTR
    ## 5 4208633 4209515      * LINE
    ## 6 4584673 4585628      * LINE

Order the rows by “combinedRank”. The rownames of the
`annotated.tab.\$ANNOTATION.order` data frame refers to the original row
numbers of the annotated.tab.sites data frame, which is not ordered.

``` r
annotated.tab.rmsk.order <- annotated.tab.rmsk[order(annotated.tab.rmsk[,"combinedRank"]),]
head(annotated.tab.rmsk.order,3)
```

    ##       mean.mean.g1 mean.mean.g2 mean.mean.diff mean.mean.quot.log2   comb.p.val
    ## 12842    0.1096282    0.5779636     -0.4683353           -2.297169 2.902269e-13
    ## 2839     0.1322201    0.4972435     -0.3650233           -1.842680 3.274019e-14
    ## 22365    0.1141823    0.4741471     -0.3599648           -1.962986 6.936504e-13
    ##       comb.p.adj.fdr num.sites mean.num.na.g1 mean.num.na.g2 mean.mean.covg.g1
    ## 12842   1.031741e-09         1              0              0          24.00000
    ## 2839    9.113735e-10         2              0              0          21.16667
    ## 22365   1.103012e-09         1              0              0          21.00000
    ##       mean.mean.covg.g2 mean.nsamples.covg.thresh.g1 mean.nsamples.covg.thresh.g2
    ## 12842          28.33333                            3                            3
    ## 2839           28.00000                            3                            3
    ## 22365          20.00000                            3                            3
    ##       combinedRank Chromosome     Start       End Strand name
    ## 12842            4       chr7 110920614 110920788      * SINE
    ## 2839            14       chr2  90676949  90677308      * LINE
    ## 22365           17      chr13 102691299 102691567      * LINE

``` r
rownames(annotated.tab.rmsk.order)[1:10]
```

    ##  [1] "12842" "2839"  "22365" "24294" "6831"  "20269" "17342" "17343" "7604"  "8467"

`annotated.tab.\$ANNOTATION.order` data frame does not contain any
information about beta values. These are extracted from rnb.set using
`meth()`, annotation name, and row index from
`annotated.tab.\$ANNOTATION.order`.

``` r
topDiffBeta.rmsk.full <- meth(rnb.set.full, type="rmsk")[as.integer(rownames(annotated.tab.rmsk.order)),]
```

#### Heatmap(s)

Differentially methylated probes annotated with RepeatMasker

``` r
colSide <- annotated.tab.rmsk.order$name

colSide <-gsub("LTR","slateblue3",colSide)
colSide <- gsub("LINE","green3",colSide)
colSide <- gsub("SINE","green3",colSide)
colSide <- gsub("DNA","black",colSide)
colSide <- gsub("Simple_repeat","yellow3",colSide)
colSide <- gsub("Other","blue",colSide)
colSide <- gsub("Unknown","blue",colSide)

par(font=2, font.axis=2, font.lab=2,cex.lab=1.5, cex.axis=1.5, lty=1)
heatmap.2(topDiffBeta.rmsk.full[1:200,c(1,5,9,2,6,10,3,7,11,4,8,12)],
          col=colfunc(10), scale="none", Rowv=F, Colv=T,
          cexCol=1, labCol=NA, labRow=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "column",
          symkey=FALSE,symbreaks=F,revC = FALSE, RowSideColors = colSide[1:200],
          breaks=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
          lmat=rbind(c(5,0,4),c(3,1,2)),
          lhei=c(1,4),
          lwid=c(1.75,0.1, 3.5),margins=c(5,12),cexRow=1.2)
```

    ## Error in plot.new(): figure margins too large

![](RnBeads_with_RepeatMasker_files/figure-gfm/heatmap-1.png)<!-- -->

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-müller2019" class="csl-entry">

Müller, Fabian, Michael Scherer, Yassen Assenov, Pavlo Lutsik, Jörn
Walter, Thomas Lengauer, and Christoph Bock. 2019. “RnBeads 2.0:
Comprehensive Analysis of DNA Methylation Data.” *Genome Biology* 20
(1): 55. <https://doi.org/10.1186/s13059-019-1664-9>.

</div>

</div>

[^1]: Run rmarkdown::render(“RnBeads_with_RepeatMasker.Rmd”) in the same
    workspace as .R to avoid worrying about global environment.

[^2]: rnb.set object is not stored well in the global environment. If
    there are errors, re-run the script from top to bottom (remove
    output folders “result” and “resultFull”). There is a way to write
    theses objects to disk instead of RAM, but I didn’t want to lose any
    disk space.
