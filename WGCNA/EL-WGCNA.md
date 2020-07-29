EL-WGCNA
================
Maggie Schedl
7/21/2020

``` r
library(WGCNA)
```

    ## Loading required package: dynamicTreeCut

    ## Loading required package: fastcluster

    ##
    ## Attaching package: 'fastcluster'

    ## The following object is masked from 'package:stats':
    ##
    ##     hclust

    ##

    ##
    ## Attaching package: 'WGCNA'

    ## The following object is masked from 'package:stats':
    ##
    ##     cor

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.3     ✓ dplyr   1.0.0
    ## ✓ tidyr   1.1.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ##
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ##
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

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
    ##     union, unique, unsplit, which, which.max, which.min

    ##
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ##
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ##
    ##     expand

    ## The following object is masked from 'package:base':
    ##
    ##     expand.grid

    ## Loading required package: IRanges

    ##
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ##
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ##
    ##     reduce

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ##
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ##
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ##
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ##
    ##     count

    ##
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ##
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:purrr':
    ##
    ##     simplify

    ## The following objects are masked from 'package:base':
    ##
    ##     aperm, apply, rowsum

# Data input and cleaning <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-DataInput.pdf>

Already have a counts file that is varience stabilizing transformed from
the stage design from the DESeq2 analysis That will be the input

``` r
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

vst_s_counts <- read.table("vst_stage_counts.txt", sep = "\t")
# change names
colnames(vst_s_counts) <- c("29_4cell_rep_1" ,"29_4cell_rep_2" ,"29_4cell_rep_3", "29_blast_rep_1", "29_blast_rep_2" ,"29_blast_rep_3", "29_gast_rep_1",  "29_gast_rep_2",  "29_gast_rep_3", "29_larv_rep_1",  "29_larv_rep_2",  "29_larv_rep_3", "33_4cell_rep_1", "33_4cell_rep_2", "33_4cell_rep_3", "33_blast_rep_1", "33_blast_rep_2", "33_blast_rep_3", "33_gast_rep_1", "33_gast_rep_2",  "33_gast_rep_3",  "33_larv_rep_1",  "33_larv_rep_2",  "33_larv_rep_3", "eggs_rep_1" , "eggs_rep_2",  "eggs_rep_3")
head(vst_s_counts)
```

    ##                        29_4cell_rep_1 29_4cell_rep_2 29_4cell_rep_3
    ## TRINITY_DN100003_c0_g1       6.841787       6.728011       6.728011
    ## TRINITY_DN10000_c0_g1        8.752377       9.299948       9.625846
    ## TRINITY_DN100013_c0_g1       6.981198       6.847842       6.975697
    ## TRINITY_DN100032_c0_g1       6.728011       6.728011       6.728011
    ## TRINITY_DN100034_c0_g1       7.424897       6.969098       9.404068
    ## TRINITY_DN100035_c0_g1       6.842057       6.728011       6.728011
    ##                        29_blast_rep_1 29_blast_rep_2 29_blast_rep_3
    ## TRINITY_DN100003_c0_g1       6.728011       6.941719       6.728011
    ## TRINITY_DN10000_c0_g1       10.246692      10.018470      10.111578
    ## TRINITY_DN100013_c0_g1       6.728011       6.728011       6.728011
    ## TRINITY_DN100032_c0_g1       6.951739       6.728011       6.728011
    ## TRINITY_DN100034_c0_g1       7.115040       6.728011       7.973263
    ## TRINITY_DN100035_c0_g1       6.728011       6.728011       6.728011
    ##                        29_gast_rep_1 29_gast_rep_2 29_gast_rep_3 29_larv_rep_1
    ## TRINITY_DN100003_c0_g1      6.937392      6.728011      7.309218      6.728011
    ## TRINITY_DN10000_c0_g1       9.908478      9.546428      9.454904      9.806620
    ## TRINITY_DN100013_c0_g1      6.875572      6.728011      6.938358      6.728011
    ## TRINITY_DN100032_c0_g1      7.090861      6.728011      7.244085      7.148334
    ## TRINITY_DN100034_c0_g1      6.728011      6.728011      7.132576      6.887480
    ## TRINITY_DN100035_c0_g1      7.119796      7.124233      7.073104      7.117383
    ##                        29_larv_rep_2 29_larv_rep_3 33_4cell_rep_1
    ## TRINITY_DN100003_c0_g1      7.206892      7.064670       6.728011
    ## TRINITY_DN10000_c0_g1       8.935433      9.152730       8.609209
    ## TRINITY_DN100013_c0_g1      6.728011      6.865192       7.003605
    ## TRINITY_DN100032_c0_g1      7.157699      7.276963       6.728011
    ## TRINITY_DN100034_c0_g1      6.728011      6.728011       7.523401
    ## TRINITY_DN100035_c0_g1      7.420164      7.309879       6.852180
    ##                        33_4cell_rep_2 33_4cell_rep_3 33_blast_rep_1
    ## TRINITY_DN100003_c0_g1       6.728011       6.728011       6.929967
    ## TRINITY_DN10000_c0_g1        9.207888       9.636638      10.014947
    ## TRINITY_DN100013_c0_g1       6.728011       6.919949       6.728011
    ## TRINITY_DN100032_c0_g1       6.728011       6.728011       6.871203
    ## TRINITY_DN100034_c0_g1       6.907339       9.352878       7.078306
    ## TRINITY_DN100035_c0_g1       6.728011       6.728011       6.930449
    ##                        33_blast_rep_2 33_blast_rep_3 33_gast_rep_1
    ## TRINITY_DN100003_c0_g1       6.886602       7.038232      7.027424
    ## TRINITY_DN10000_c0_g1        9.688650       9.636674      9.840236
    ## TRINITY_DN100013_c0_g1       6.885991       7.004515      6.728011
    ## TRINITY_DN100032_c0_g1       6.728011       7.095649      6.940402
    ## TRINITY_DN100034_c0_g1       6.952878       8.721305      6.878374
    ## TRINITY_DN100035_c0_g1       6.728011       6.867299      7.245988
    ##                        33_gast_rep_2 33_gast_rep_3 33_larv_rep_1 33_larv_rep_2
    ## TRINITY_DN100003_c0_g1      7.148118      7.220518      6.881662      7.253414
    ## TRINITY_DN10000_c0_g1       9.226925      9.402862      9.645214      8.779761
    ## TRINITY_DN100013_c0_g1      6.728011      6.864707      6.881084      6.920057
    ## TRINITY_DN100032_c0_g1      7.061331      7.161419      7.546698      7.199519
    ## TRINITY_DN100034_c0_g1      6.728011      6.865632      6.882127      7.033382
    ## TRINITY_DN100035_c0_g1      7.221143      6.922436      6.994519      7.288240
    ##                        33_larv_rep_3 eggs_rep_1 eggs_rep_2 eggs_rep_3
    ## TRINITY_DN100003_c0_g1      6.912144   6.728011   6.728011   6.728011
    ## TRINITY_DN10000_c0_g1       9.162756   8.534400   8.775307   9.732666
    ## TRINITY_DN100013_c0_g1      6.728011   6.728011   6.728011   6.906318
    ## TRINITY_DN100032_c0_g1      7.522915   6.728011   6.728011   6.728011
    ## TRINITY_DN100034_c0_g1      6.858652   7.489706   6.921235   9.801128
    ## TRINITY_DN100035_c0_g1      7.349524   6.728011   6.728011   6.728011

“each row corresponds to a gene and column to a sample or auxiliary
information”

``` r
# check for genes with too many missing values
gsg = goodSamplesGenes(vst_s_counts, verbose = 3)
```

    ##  Flagging genes and samples with too many missing values...
    ##   ..step 1

``` r
gsg$allOK
```

    ## [1] TRUE

``` r
# good!
```

Do not need to remove any samples to do the analysis

Need to subset these genes to a more managable number, 5000 is good (out
of the 13296 in the full dataset), these should be the genes with the
most variance (most signal\!)

``` r
matrix <- as.matrix(vst_s_counts) # make into a matrix so the vars fucntion will take it as an input and
vars <- rowVars(matrix) # apply the varience function to it to get the varience in each row
vars <- as.data.frame(vars) # change varience into data frame
vst_counts <- as.data.frame(matrix) # turn the counts data back into a dataframe
vst_s_counts$Genes <- rownames(vst_s_counts) # create a new column in the original dataframe that makes the row names an actual column
rownames(vars) <- vst_s_counts$Names # use that column as the template for the rownames of the vars dataframe
vst_counts_var <- cbind(vst_counts, vars) #add the row varience to the counts dataframe , binds by having the same rownames
vst_counts_var <- vst_counts_var[order(-vars),] # order this dataframe by the varience column and in decreasing order, so highest varience should be first
orderedvst_counts_var <- head(vst_counts_var, 5000) # then only take the top 5000 amount for a scaled down dataset
orderedvst_counts_var$vars <- NULL #remove the last column because it's not needed anymore
torderedvst_counts_var <- as.data.frame(t(as.matrix(orderedvst_counts_var))) #transpose to right format where samples are rows
# head(torderedvst_counts_var)
```

# Exploratory sample clustering

``` r
sampleTree = hclust(dist(torderedvst_counts_var), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
```
![1](/images/1.png)

There are no ovbsious outliers but they completely cluster by stage
except for 4 cells and eggs. There is also the large break between early
and later development

I might need to split the dataset between those two big groups

# Step by step network construction <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-man.pdf>

Choosing the soft-thresholding power: analysis of network topology

“Constructing a weighted gene network entails the choice of the soft
thresholding power β to which co-expression similarity is raised to
calculate adjacency”

``` r
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(torderedvst_counts_var, powerVector = powers, verbose = 5, networkType = "signed") #make sure to specify network type
```

    ## pickSoftThreshold: will use block size 5000.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 5000 of 5000

    ## Warning: executing %dopar% sequentially: no parallel backend registered

    ##    Power SFT.R.sq    slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1 0.039000  2.01000        0.97700    2580      2760   2940
    ## 2      2 0.247000  1.64000        0.04060    1920      1940   2490
    ## 3      3 0.682000  1.58000        0.81900    1580      1500   2230
    ## 4      4 0.744000  1.26000        0.81400    1360      1340   2030
    ## 5      5 0.737000  0.96900        0.76800    1190      1210   1870
    ## 6      6 0.657000  0.73700        0.64600    1060      1090   1720
    ## 7      7 0.588000  0.57800        0.56200     949       994   1600
    ## 8      8 0.502000  0.45100        0.46000     859       907   1490
    ## 9      9 0.382000  0.34700        0.30700     783       830   1400
    ## 10    10 0.242000  0.25500        0.14700     717       763   1310
    ## 11    12 0.049200  0.10600       -0.06770     610       648   1170
    ## 12    14 0.000019  0.00204       -0.09120     526       554   1050
    ## 13    16 0.030300 -0.08320       -0.00489     459       479    957
    ## 14    18 0.093400 -0.15900        0.09610     405       417    875
    ## 15    20 0.167000 -0.23100        0.20200     359       365    807
    ## 16    22 0.242000 -0.30200        0.30700     321       323    747
    ## 17    24 0.302000 -0.35800        0.38400     288       286    695
    ## 18    26 0.358000 -0.39500        0.47200     261       256    648
    ## 19    28 0.407000 -0.43400        0.54400     236       230    607
    ## 20    30 0.449000 -0.47800        0.58900     215       207    570

``` r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```
![2](/image/2.png)

This is not how the scale independence graph should look, I need a
singed R^2 above 0.8 and below 30 for the power I want to do a signed
network. “If the scale-free topology fit index fails to reach values
above 0.8 for reasonable powers (less than 15 for unsigned or signed
hybrid networks, and less than 30 for signed networks) and the mean
connectivity remains relatively high (in the hundreds or above), chances
are that the data exhibit a strong driver that makes a subset of the
samples globally different from the rest. The difference causes high
correlation among large groups of genes which invalidates the assumption
of the scale-free topology
approximation.”<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html>

I could change the way I normalized the data, try a TMM normalization,
or I could subset the data into the two major groups.

I think I will try subsetting first?

Select out rows in dataframe with only 4cell and egg samples

``` r
# need the first 3, a middle 3 and the last 3

fourcell1 <- torderedvst_counts_var[1:3,]
fourcell2 <- torderedvst_counts_var[13:15,]
eggs <- torderedvst_counts_var[25:27,]
earlydevcounts <- rbind(eggs, fourcell1, fourcell2)
```

I am worried now though that this may be too few samples..

Select out rows in dataframe with blastula, gastrula, and larvae

``` r
latedevcounts1 <- torderedvst_counts_var[4:12,]
latedevcounts12 <- torderedvst_counts_var[16:24,]
latedevcounts <- rbind(latedevcounts1, latedevcounts12)
```

Try picking soft thresholding power again

early development

``` r
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(earlydevcounts, powerVector = powers, verbose = 5, networkType = "signed") #make sure to specify network type
```

    ## pickSoftThreshold: will use block size 5000.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 5000 of 5000

    ## Warning in eval(xpr, envir = envir): Some correlations are NA in block 1 :
    ## 5000 .

    ##    Power SFT.R.sq    slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1 0.043000  2.09000          0.849  2530.0    2540.0   2800
    ## 2      2 0.110000  3.50000          0.820  1590.0    1610.0   1860
    ## 3      3 0.245000  4.71000          0.599  1110.0    1130.0   1380
    ## 4      4 0.299000  4.59000          0.365   832.0     854.0   1100
    ## 5      5 0.329000  3.67000          0.198   652.0     674.0    906
    ## 6      6 0.844000  1.68000          0.800   529.0     550.0    773
    ## 7      7 0.870000  0.98000          0.835   439.0     459.0    672
    ## 8      8 0.847000  1.08000          0.861   372.0     390.0    592
    ## 9      9 0.738000  0.92000          0.827   320.0     335.0    528
    ## 10    10 0.620000  0.69300          0.729   279.0     292.0    476
    ## 11    12 0.341000  0.38900          0.515   218.0     227.0    397
    ## 12    14 0.085300  0.17300          0.365   176.0     182.0    340
    ## 13    16 0.000279  0.00967          0.352   145.0     149.0    297
    ## 14    18 0.032900 -0.10500          0.435   122.0     123.0    263
    ## 15    20 0.114000 -0.20700          0.517   105.0     104.0    235
    ## 16    22 0.217000 -0.32500          0.595    90.7      88.3    214
    ## 17    24 0.298000 -0.42200          0.666    79.4      75.7    197
    ## 18    26 0.361000 -0.50400          0.717    70.1      65.6    182
    ## 19    28 0.418000 -0.58300          0.751    62.4      57.2    169
    ## 20    30 0.456000 -0.65800          0.775    55.9      50.5    158

``` r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```
![3](/images/3.png)

This still does not look good. I don’t know why there is the v shape in
the scale independence. It partly could have to do with the smaller
number of samples. To accurately run the program it has to look like the
one below

later development

``` r
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(latedevcounts, powerVector = powers, verbose = 5, networkType = "signed") #make sure to specify network type
```

    ## pickSoftThreshold: will use block size 5000.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 5000 of 5000
    ##    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1 0.038500  1.5600         0.9060  2520.0    2530.0   2710
    ## 2      2 0.009140 -0.3070        -0.2520  1640.0    1590.0   2060
    ## 3      3 0.000387  0.0223         0.0944  1200.0    1180.0   1750
    ## 4      4 0.001630  0.0305         0.2820   937.0     929.0   1540
    ## 5      5 0.017300 -0.0723         0.3840   763.0     751.0   1390
    ## 6      6 0.182000 -0.1980         0.5400   639.0     621.0   1270
    ## 7      7 0.463000 -0.3090         0.7050   546.0     520.0   1170
    ## 8      8 0.680000 -0.4050         0.8170   474.0     441.0   1090
    ## 9      9 0.796000 -0.4800         0.8850   417.0     379.0   1020
    ## 10    10 0.851000 -0.5390         0.9150   370.0     328.0    959
    ## 11    12 0.884000 -0.6330         0.9100   299.0     249.0    856
    ## 12    14 0.890000 -0.7000         0.9080   248.0     193.0    771
    ## 13    16 0.899000 -0.7560         0.9130   209.0     152.0    701
    ## 14    18 0.898000 -0.7980         0.9110   179.0     121.0    641
    ## 15    20 0.887000 -0.8350         0.8990   155.0      97.7    590
    ## 16    22 0.883000 -0.8580         0.8960   136.0      79.4    545
    ## 17    24 0.885000 -0.8760         0.8980   120.0      65.6    505
    ## 18    26 0.900000 -0.8850         0.9150   107.0      54.2    470
    ## 19    28 0.909000 -0.9000         0.9220    95.9      45.2    439
    ## 20    30 0.921000 -0.9090         0.9330    86.3      37.8    411

``` r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

[4](/images/4.png)

This one looks really great\! power 10 has an R^2 of 0.9150

I will go forward with only the late development samples I suppose.

creating the network adjacency matrix, using the power picked and
specifying the sign. The adjacency matrix is for all samples and is
between each gene.

``` r
vsoftPower = 10 # seting the soft power to 10
vadjacency = adjacency(latedevcounts, power = vsoftPower, type= "signed")
```

Create topological overlap matrix from adjacency matrix

``` r
# Turn adjacency into topological overlap
vTOM = TOMsimilarity(vadjacency)
```

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
vdissTOM = 1-vTOM # this is the dissimilarity, which is 1 minus the similarity
```

hierarchical clustering tree (dendrogram) of genes

``` r
# Call the hierarchical clustering function
vgeneTree = hclust(as.dist(vdissTOM), method = "average")

# Average linkage: average the dissimilarities between all objects
#single (minimum dissimilarity) and complete (maximum dissimilarity) are other options, but seem too harsh

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(vgeneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04)
```

![5](/images/5.png)

Each one of the lines on this are genes. “Branches of the dendrogram
group together densely interconnected, highly co-expressed genes.” Now I
need to cut the branches into coexpression modules

These are “dynamically” cut because no specific height works for each
individual branch “DeepSplit controls how finely the branches should be
split. Higher values give more smaller modules, low values (0) give
fewer larger modules”

``` r
# chose a number of genes to be the minium in a module, the tutorial said 30 as an example
minModuleSize = 30; #must have 30 genes in it to be a module
# Module identification using dynamic tree cut:
vdynamicMods = cutreeDynamic(dendro = vgeneTree, distM = vdissTOM,
deepSplit = 2, pamRespectsDendro = FALSE, # deep split set at 2, might want to change...
minClusterSize = minModuleSize)
```

    ##  ..cutHeight not given, setting it to 0.985  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

``` r
table(vdynamicMods)
```

    ## vdynamicMods
    ##    1    2    3    4    5    6    7    8    9   10   11   12   13   14
    ## 1283  880  717  516  311  255  251  214  148  128  121   84   58   34

This gives 14 co-expression modules

Visualize modules on the tree

``` r
# Convert numeric lables into colors
vdynamicColors = labels2colors(vdynamicMods)
table(vdynamicColors)
```

    ## vdynamicColors
    ##       black        blue       brown        cyan       green greenyellow
    ##         251         880         717          34         311         121
    ##     magenta        pink      purple         red      salmon         tan
    ##         148         214         128         255          58          84
    ##   turquoise      yellow
    ##        1283         516

``` r
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(vgeneTree, vdynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
```

![6](/images/6.png)

Some modules my need to be collapsed into each other if they are also
co-expressed. “The Dynamic Tree Cut may identify modules whose
expression profiles are very similar. It may be prudent to merge such
modules since their genes are highly co-expressed. To quantify
co-expression similarity of entire modules, we calculate their
eigengenes and cluster them on their correlation:”

``` r
# Calculate eigengenes
vMEList = moduleEigengenes(latedevcounts, colors = vdynamicColors, trapErrors = FALSE, excludeGrey = TRUE)
vMEs = vMEList$eigengenes
# Calculate dissimilarity of module eigengenes
vMEDiss = 1-cor(vMEs);
# Cluster module eigengenes
vMETree = hclust(as.dist(vMEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(vMETree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
```


``` r
vMEList = moduleEigengenes(latedevcounts, colors = vdynamicColors, trapErrors = FALSE, excludeGrey = TRUE)
vMEs = vMEList$eigengenes
# Calculate dissimilarity of module eigengenes
vMEDiss = 1-cor(vMEs);
# Cluster module eigengenes
vMETree = hclust(as.dist(vMEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(vMETree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
vMEDissThres = 0.15 # this means to merge modlues if they 85% correlated to eachother
# Plot the cut line into the dendrogram
abline(h=vMEDissThres, col = "red")
# Call an automatic merging function
vmerge = mergeCloseModules(latedevcounts, vdynamicColors, cutHeight = vMEDissThres, verbose = 3)
```

    ##  mergeCloseModules: Merging modules whose distance is less than 0.15
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 14 module eigengenes in given set.
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 13 module eigengenes in given set.
    ##    Calculating new MEs...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 13 module eigengenes in given set.

``` r
# The merged module colors
vmergedColors = vmerge$colors
# Eigengenes of the new merged modules:
vmergedMEs = vmerge$newMEs
```

![7](/images/7.png)

This would mean 2 modules are merged together

``` r
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(vgeneTree, cbind(vdynamicColors, vmergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
```

![8](/images/8.png)

# Relating modules to stage and temp treatment information

Changed the treatment/trait information into numbers

``` r
Traitdat <- as.data.frame(read.csv("el-treatment-data.csv"))
head(Traitdat)
```

    ##           Sample Stage Temp
    ## 1 29_4cell_rep_1     2    1
    ## 2 29_4cell_rep_2     2    1
    ## 3 29_4cell_rep_3     2    1
    ## 4 29_blast_rep_1     3    1
    ## 5 29_blast_rep_2     3    1
    ## 6 29_blast_rep_3     3    1

``` r
Traitdat2 <- Traitdat[,-1] #make column 1 into row names
rownames(Traitdat2) <- Traitdat[,1]
print(Traitdat2)
```

    ##                Stage Temp
    ## 29_4cell_rep_1     2    1
    ## 29_4cell_rep_2     2    1
    ## 29_4cell_rep_3     2    1
    ## 29_blast_rep_1     3    1
    ## 29_blast_rep_2     3    1
    ## 29_blast_rep_3     3    1
    ## 29_gast_rep_1      4    1
    ## 29_gast_rep_2      4    1
    ## 29_gast_rep_3      4    1
    ## 29_larv_rep_1      5    1
    ## 29_larv_rep_2      5    1
    ## 29_larv_rep_3      5    1
    ## 33_4cell_rep_1     2    2
    ## 33_4cell_rep_2     2    2
    ## 33_4cell_rep_3     2    2
    ## 33_blast_rep_1     3    2
    ## 33_blast_rep_2     3    2
    ## 33_blast_rep_3     3    2
    ## 33_gast_rep_1      4    2
    ## 33_gast_rep_2      4    2
    ## 33_gast_rep_3      4    2
    ## 33_larv_rep_1      5    2
    ## 33_larv_rep_2      5    2
    ## 33_larv_rep_3      5    2
    ## eggs_rep_1         1    1
    ## eggs_rep_2         1    1
    ## eggs_rep_3         1    1

Think I need to subset again by only the rows that are for the later
stages

``` r
Trait1 <- Traitdat2[4:12,]
Trait2 <- Traitdat2[16:24,]
Traitdat3 <- rbind(Trait1, Trait2)
print(Traitdat3)
```

    ##                Stage Temp
    ## 29_blast_rep_1     3    1
    ## 29_blast_rep_2     3    1
    ## 29_blast_rep_3     3    1
    ## 29_gast_rep_1      4    1
    ## 29_gast_rep_2      4    1
    ## 29_gast_rep_3      4    1
    ## 29_larv_rep_1      5    1
    ## 29_larv_rep_2      5    1
    ## 29_larv_rep_3      5    1
    ## 33_blast_rep_1     3    2
    ## 33_blast_rep_2     3    2
    ## 33_blast_rep_3     3    2
    ## 33_gast_rep_1      4    2
    ## 33_gast_rep_2      4    2
    ## 33_gast_rep_3      4    2
    ## 33_larv_rep_1      5    2
    ## 33_larv_rep_2      5    2
    ## 33_larv_rep_3      5    2

``` r
#set up for using the module colors in further code
# Rename to moduleColors
vmoduleColors = vmergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
vmoduleLabels = match(vmoduleColors, colorOrder)-1
vMEs = vmergedMEs
```

“Robust correlation. The default correlation method in all functions in
WGCNA is standard Pearson correlation. In general, unless there is good
reason to believe that there are no outlier measurements, we recommend
(and use ourselves) the biweight mid-correlation as a robust
alternative. This is implemented in WGCNA function bicor. Many WGCNA
functions take the argument corFnc that allows one to specify an
alternative correlation function to the standard cor and bicor is one
option. Additional arguments to the correlation function can be
specified using the argument corOptions (depending on function, this
argument may require one of two alternate forms, please see the help for
each function for details). In certain functions, notably the of the
blockwise family, correlation function cannot be specified directly as a
function; rather, one must use the argument corType to specify either
Pearson or biweight mid-correlation.”

I think I want to use the robust correlation method. However the temp
treatment is bianary but the stage treatment is not “Dealing with binary
data. When relating high-throughput data x to binary variable y such as
sample traits, one can use argument robustY = FALSE to turn off the
robust treatment for the y argment of bicor. This results in a hybrid
robust-Pearson correlation as described in Langfelder and Horvath
(2011)”

“The options robustX, robustY allow the user to revert the calculation
to standard correlation calculation. This is important, for example, if
any of the variables is binary (or, more generally, discrete) as in such
cases the robust methods produce meaningless results”

``` r
# Define numbers of genes and samples
vnGenes = ncol(latedevcounts);
vnSamples = nrow(latedevcounts);
# Recalculate MEs with color labels
vMEs0 = moduleEigengenes(latedevcounts, vmoduleColors)$eigengenes
vMEs = orderMEs(vMEs0)
vmoduleTraitCor = bicor(vMEs, Traitdat3, robustX = TRUE, robustY = FALSE, use = "all.obs"); #try pearson robust correlation
vmoduleTraitPvalue = corPvalueStudent(vmoduleTraitCor, vnSamples)
```

Visualize

``` r
sizeGrWindow(10,6)
# Will display correlations and their p-values
vtextMatrix = paste(signif(vmoduleTraitCor, 2), "\n(",
signif(vmoduleTraitPvalue, 1), ")", sep = "");
dim(vtextMatrix) = dim(vmoduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = vmoduleTraitCor,
xLabels = names(Traitdat3),
yLabels = names(vMEs),
ySymbols = names(vMEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = vtextMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
```

![9](/images/9.png)

The turquoise module is highly positively correlated with the stage. The
yellow module is highly negatively correlated with the stage. There are
a couple other singificant ones. There aren’t any serious correlations
with temperature but the green yellow module is slight (but not
significant 0.1)

I am not sure what the negative correlation means exactly.

# Running Model on All Samples

I could potentially analyze all of the stages in one
<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html>
"Lack of scale-free topology fit by itself does not invalidate the data,
but should be looked into carefully. It always helps to plot the sample
clustering tree and any technical or biological sample information below
it as in Figure 2 of Tutorial I, section 1; strong clusters in the
clustering tree indicate globally different groups of samples. It could
be the result a technical effect such as a batch effect, biological
heterogeneity (e.g., a data set consisting of samples from 2 different
tissues), or strong changes between conditions (say in a time series).
One should investigate carefully whether there is sample heterogeneity,
what drives the heterogeneity, and whether the data should be adjusted
(see previous point).

If the lack of scale-free topology fit turns out to be caused by an
interesting biological variable that one does not want to remove (i.e.,
adjust the data for), the appropriate soft-thresholding power can be
chosen based on the number of samples as in the table below. This table
has been updated in December 2017 to make the resulting networks
conservative."

I would say that stage is like a time series, and clustering is probably
the main reason there isn’t scale free topology, which is an interesting
biological variable. So I potentially could keep going with the full
dataset and just use the acceptable soft threshold power for a signed
network of 16.

Going forward below with Full Dataset

Via the above website, for a signed network with 20-30 samples the
recommended soft power value is 16 Putting an f in front of all the
terms because it’s with the full dataset

``` r
fsoftPower = 16
fadjacency = adjacency(torderedvst_counts_var, power = fsoftPower, type= "signed")
# Turn adjacency into topological overlap
fOM = TOMsimilarity(fadjacency)
```

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
fdissTOM = 1-fOM # this is the dissimilarity, which is 1 minus the similarity
```

Hierarchical clustering for module dectection

``` r
# Call the hierarchical clustering function
fgeneTree = hclust(as.dist(fdissTOM), method = "average");

# Average linkage: average the dissimilarities between all objects
#single (minimum dissimilarity) and complete (maximum dissimilarity) are other options, but seem too harsh

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(fgeneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04)
```

![10](/images/10.png)

This looks pretty good, doesn’t seem to be a problem will running will
all samples

Cutting the dendrogram into modules

Right now trying with 30 genes as the minimum module size

``` r
minModuleSize = 30;
# Module identification using dynamic tree cut:
fdynamicMods = cutreeDynamic(dendro = fgeneTree, distM = fdissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize)
```

    ##  ..cutHeight not given, setting it to 0.993  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

``` r
table(fdynamicMods)
```

    ## fdynamicMods
    ##    1    2    3    4    5    6    7    8    9   10   11
    ## 1464  810  747  599  501  259  248  214   70   51   37

``` r
# Convert numeric lables into colors
fdynamicColors = labels2colors(fdynamicMods)
table(fdynamicColors)
```

    ## fdynamicColors
    ##       black        blue       brown       green greenyellow     magenta
    ##         248         810         747         501          37          70
    ##        pink      purple         red   turquoise      yellow
    ##         214          51         259        1464         599

``` r
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(fgeneTree, fdynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
```

![11](/images/11.png)

11 modules

Correlating modules to see if any should be merged

There aren’t a lot of modules and they split pretty nicely, so not sure
about if this is necessary

``` r
# Calculate eigengenes
fMEList = moduleEigengenes(torderedvst_counts_var, colors = fdynamicColors, trapErrors = FALSE, excludeGrey = TRUE) # grey is a non-module
fMEs = fMEList$eigengenes
# Calculate dissimilarity of module eigengenes
fMEDiss = 1-cor(fMEs);
# Cluster module eigengenes
fMETree = hclust(as.dist(fMEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(fMETree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

fMEDissThres = 0.1 # merge modules if they are 90% similar
# Plot the cut line into the dendrogram
abline(h=fMEDissThres, col = "red")
# Call an automatic merging function
fmerge = mergeCloseModules(torderedvst_counts_var, fdynamicColors, cutHeight = fMEDissThres, verbose = 3)
```

    ##  mergeCloseModules: Merging modules whose distance is less than 0.1
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 11 module eigengenes in given set.
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 9 module eigengenes in given set.
    ##    Calculating new MEs...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 9 module eigengenes in given set.

``` r
# The merged module colors
fmergedColors = fmerge$colors;
# Eigengenes of the new merged modules:
fmergedMEs = fmerge$newMEs
```

![12](/images/12.png)

Visualizing merged modules

``` r
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(fgeneTree, cbind(fdynamicColors, fmergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
```

![13](/images/13.png)

Correlating modules with treatment and life stage information

Already loaded in the trait data, Traitdat2 should be good

``` r
#set up for using the module colors in further code
# Rename to moduleColors
fmoduleColors = fmergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
fmoduleLabels = match(fmoduleColors, colorOrder)-1
fMEs = fmergedMEs
```

``` r
# Define numbers of genes and samples
fnGenes = ncol(torderedvst_counts_var);
fnSamples = nrow(torderedvst_counts_var);
# Recalculate MEs with color labels
fMEs0 = moduleEigengenes(torderedvst_counts_var, fmoduleColors)$eigengenes
fMEs = orderMEs(fMEs0)
fmoduleTraitCor = bicor(fMEs, Traitdat2, robustX = TRUE, robustY = FALSE, use = "all.obs"); #try pearson robust correlation
fmoduleTraitPvalue = corPvalueStudent(fmoduleTraitCor, fnSamples)
```

``` r
sizeGrWindow(10,6)
# Will display correlations and their p-values
ftextMatrix = paste(signif(fmoduleTraitCor, 2), "\n(",
signif(fmoduleTraitPvalue, 1), ")", sep = "");
dim(ftextMatrix) = dim(fmoduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = fmoduleTraitCor,
xLabels = names(Traitdat2),
yLabels = names(fMEs),
ySymbols = names(fMEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = ftextMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
```

![14](/images/14.png)

Black module and green module are significantly positively correlated
with the life stage.

I am still not sure what the negative correlations mean but there are
two signigicantly negatively correlated ones as well

Visualizing Networks (using the full model)

``` r
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
 # dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# I did do that!
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = fdissTOM^7; # how do you decide this
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, fgeneTree, fmoduleColors, main = "Network heatmap plot, all genes")
```

![15](/images/15.png)

There looks like a hugely correlated block of modules. I will have to
think about this one and what it means. I wouldn’t be surprised if a lot
of genes are correlated with development and are an overwhelming amount
in the dataset

Investigating this further I can look at the adjacency of the module
eigengens as well as with the stage variable to see what the
relationship is

``` r
Stage = as.data.frame(Traitdat2$Stage);
names(Stage) = "stage"
MET = orderMEs(cbind(fMEs, Stage))
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)
```

![16](/images/16.png)
![17](/images/17.png)

This looks like the green, black, and blue modules are highly related
with the stage variable, which makes sense that it is the most visible
in the dataset. And the Purple, red and turquoise are negatively related
to it.

# Saving the genes in each module

``` r
# Black module
BlackGeneNames <- names(torderedvst_counts_var)[fmoduleColors=="black"]
BlackModGeneNames <- as.data.frame(BlackGeneNames) # make into dataframe
write.csv(BlackModGeneNames, "BlackModGeneNames.csv") # save as csv

#green module
greenGeneNames <- names(torderedvst_counts_var)[fmoduleColors=="green"]
greenModGeneNames <- as.data.frame(greenGeneNames) # make into dataframe
write.csv(greenModGeneNames, "greenModGeneNames.csv") # save as csv

# blue module
BlueGeneNames <- names(torderedvst_counts_var)[fmoduleColors=="blue"]
BlueModGeneNames <- as.data.frame(BlueGeneNames) # make into dataframe
write.csv(BlueModGeneNames, "BlueModGeneNames.csv") # save as csv

# pink module
pinkGeneNames <- names(torderedvst_counts_var)[fmoduleColors=="pink"]
pinkModGeneNames <- as.data.frame(pinkGeneNames) # make into dataframe
write.csv(pinkModGeneNames, "pinkModGeneNames.csv") # save as csv

# green yellow module
gyGeneNames <- names(torderedvst_counts_var)[fmoduleColors=="greenyellow"]
gyModGeneNames <- as.data.frame(gyGeneNames) # make into dataframe
write.csv(gyModGeneNames, "gyModGeneNames.csv") # save as csv

# magenta module
magGeneNames <- names(torderedvst_counts_var)[fmoduleColors=="magenta"]
magModGeneNames <- as.data.frame(magGeneNames) # make into dataframe
write.csv(magModGeneNames, "magentaModGeneNames.csv") # save as csv

# purple module
prupGeneNames <- names(torderedvst_counts_var)[fmoduleColors=="purple"]
purpModGeneNames <- as.data.frame(prupGeneNames) # make into dataframe
write.csv(purpModGeneNames, "purpleModGeneNames.csv") # save as csv

# red module
redGeneNames <- names(torderedvst_counts_var)[fmoduleColors=="red"]
redModGeneNames <- as.data.frame(redGeneNames) # make into dataframe
write.csv(redModGeneNames, "redModGeneNames.csv") # save as csv

# turquoise module
turqGeneNames <- names(torderedvst_counts_var)[fmoduleColors=="turquoise"]
turqModGeneNames <- as.data.frame(turqGeneNames) # make into dataframe
write.csv(turqModGeneNames, "turquoiseModGeneNames.csv") # save as csv
```
