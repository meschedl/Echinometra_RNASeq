DEGs
================
Maggie Schedl
6/24/2020

``` r
#libraries
library(naniar)
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.1     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.1     ✓ dplyr   1.0.0
    ## ✓ tidyr   1.1.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(dplyr)
library(stringr)
```

Read i orthogroups file and do same file manipulations as in annotation
script

``` r
orthologs <- read.delim("LORF_EL90.pep__v__SPU_peptide.tsv", sep = "\t")
head(orthologs)
```

    ##   Orthogroup
    ## 1  OG0000000
    ## 2  OG0000000
    ## 3  OG0000001
    ## 4  OG0000001
    ## 5  OG0000002
    ## 6  OG0000003
    ##                                                                           LORF_EL90.pep
    ## 1 TRINITY_DN56081_c1_g2_i2.p1, TRINITY_DN56081_c1_g3_i1.p1, TRINITY_DN56081_c1_g1_i1.p1
    ## 2                              TRINITY_DN56081_c0_g1_i2.p1, TRINITY_DN56081_c0_g1_i5.p1
    ## 3                                                          TRINITY_DN108939_c0_g1_i1.p1
    ## 4                                                           TRINITY_DN10826_c0_g1_i6.p1
    ## 5                              TRINITY_DN17786_c0_g1_i2.p1, TRINITY_DN17786_c0_g1_i5.p1
    ## 6                             TRINITY_DN49892_c0_g1_i11.p1, TRINITY_DN35058_c0_g1_i6.p1
    ##                                                                                                                                                                                                                                                                                                  SPU_peptide
    ## 1                                                                                                                                                                                                                                                                         SPU_027035, SPU_004012, SPU_027207
    ## 2                         SPU_023183, SPU_022564, SPU_002888, SPU_008707, SPU_026036, SPU_008498, SPU_000863, SPU_018384, SPU_025179, SPU_016759, SPU_003303, SPU_003553, SPU_015033, SPU_003539, SPU_003247, SPU_006610, SPU_022001, SPU_005383, SPU_003797, SPU_017505, SPU_009659, SPU_013038, SPU_016060
    ## 3                                                                                                                                                                                                                                                                                     SPU_027574, SPU_019249
    ## 4 SPU_009978, SPU_019937, SPU_018594, SPU_001189, SPU_018839, SPU_007888, SPU_007889, SPU_027473, SPU_012635, SPU_020572, SPU_020571, SPU_009630, SPU_006586, SPU_008153, SPU_017855, SPU_024683, SPU_017108, SPU_024682, SPU_002327, SPU_011979, SPU_013146, SPU_010581, SPU_009699, SPU_011131, SPU_013233
    ## 5 SPU_001711, SPU_010740, SPU_015261, SPU_000464, SPU_014175, SPU_002486, SPU_002905, SPU_002772, SPU_007666, SPU_005645, SPU_007820, SPU_009181, SPU_015741, SPU_002828, SPU_003576, SPU_024346, SPU_014085, SPU_000892, SPU_000474, SPU_016236, SPU_013389, SPU_012142, SPU_006563, SPU_013690, SPU_015560
    ## 6                                                                                                                                                                                                                                                                                                 SPU_002006

``` r
EL.orthologs <- orthologs %>% 
  select(1,2)
head(EL.orthologs)
```

    ##   Orthogroup
    ## 1  OG0000000
    ## 2  OG0000000
    ## 3  OG0000001
    ## 4  OG0000001
    ## 5  OG0000002
    ## 6  OG0000003
    ##                                                                           LORF_EL90.pep
    ## 1 TRINITY_DN56081_c1_g2_i2.p1, TRINITY_DN56081_c1_g3_i1.p1, TRINITY_DN56081_c1_g1_i1.p1
    ## 2                              TRINITY_DN56081_c0_g1_i2.p1, TRINITY_DN56081_c0_g1_i5.p1
    ## 3                                                          TRINITY_DN108939_c0_g1_i1.p1
    ## 4                                                           TRINITY_DN10826_c0_g1_i6.p1
    ## 5                              TRINITY_DN17786_c0_g1_i2.p1, TRINITY_DN17786_c0_g1_i5.p1
    ## 6                             TRINITY_DN49892_c0_g1_i11.p1, TRINITY_DN35058_c0_g1_i6.p1

``` r
# collapse df in to one column multiple rows instead of one row multiple columns
tr.orthogroups <- separate_rows(EL.orthologs, "LORF_EL90.pep", sep = ",", convert = FALSE)
head(tr.orthogroups)
```

    ## # A tibble: 6 x 2
    ##   Orthogroup LORF_EL90.pep                 
    ##   <chr>      <chr>                         
    ## 1 OG0000000  "TRINITY_DN56081_c1_g2_i2.p1" 
    ## 2 OG0000000  " TRINITY_DN56081_c1_g3_i1.p1"
    ## 3 OG0000000  " TRINITY_DN56081_c1_g1_i1.p1"
    ## 4 OG0000000  "TRINITY_DN56081_c0_g1_i2.p1" 
    ## 5 OG0000000  " TRINITY_DN56081_c0_g1_i5.p1"
    ## 6 OG0000001  "TRINITY_DN108939_c0_g1_i1.p1"

``` r
# remove extra characters of .p#s
tr.orthogroups <-mapply(gsub, pattern = ".p1", replacement = "", tr.orthogroups)
tr.orthogroups <- as.data.frame(tr.orthogroups)
# I don't know why this always removes the column names but it does
# remove .p2
tr.orthogroups <-mapply(gsub, pattern = ".p2", replacement = "", tr.orthogroups)
tr.orthogroups <- as.data.frame(tr.orthogroups)
# remove .p3
tr.orthogroups <-mapply(gsub, pattern = ".p3", replacement = "", tr.orthogroups)
tr.orthogroups <- as.data.frame(tr.orthogroups)
# remove .p4
tr.orthogroups <-mapply(gsub, pattern = ".p4", replacement = "", tr.orthogroups)
tr.orthogroups <- as.data.frame(tr.orthogroups)
# remove .p5
tr.orthogroups <-mapply(gsub, pattern = ".p5", replacement = "", tr.orthogroups)
tr.orthogroups <- as.data.frame(tr.orthogroups)
# add column names
colnames(tr.orthogroups) <- c("Orthogroup", "transcript_id")
head(tr.orthogroups)
```

    ##   Orthogroup             transcript_id
    ## 1  OG0000000  TRINITY_DN56081_c1_g2_i2
    ## 2  OG0000000  TRINITY_DN56081_c1_g3_i1
    ## 3  OG0000000  TRINITY_DN56081_c1_g1_i1
    ## 4  OG0000000  TRINITY_DN56081_c0_g1_i2
    ## 5  OG0000000  TRINITY_DN56081_c0_g1_i5
    ## 6  OG0000001 TRINITY_DN108939_c0_g1_i1

read in SPU annotation file with the orthogroups label in it

``` r
annot.ortho.SPU <- read.delim("annot.ortho.SPU.txt", sep = "\t")
head(annot.ortho.SPU)
```

    ##       spu_id Orthogroup    family_member   common_name
    ## 1 SPU_000002  OG0001426     LNB-7TM GPCR Cub/Ldla/Gpcr
    ## 2 SPU_000003  OG0002596             <NA>     Hypp_1239
    ## 3 SPU_000006  OG0012715     RNA helicase        RigIL4
    ## 4 SPU_000007  OG0002845              WNT          Wnt7
    ## 5 SPU_000008  OG0010762 PLCc superfamily       Plcxd3L
    ## 6 SPU_000009  OG0002329             <NA>         Mfsd1
    ##                                                               synonyms
    ## 1                                                                 <NA>
    ## 2                                            hypothetical protein-1239
    ## 3                                          MDA-5 like 12, LGP2 like 12
    ## 4                             wingless-related MMTV integration site 7
    ## 5 phosphatidylinositol-specific phospholipase C, X domain containing 3
    ## 6                    major facilitator superfamily domain containing 1
    ##   best_genbank_hit
    ## 1             <NA>
    ## 2   XP_002121667.1
    ## 3             <NA>
    ## 4         AAC80433
    ## 5   XP_001371514.1
    ## 6         27503032

read in lists of DEGS

``` r
sig_LH_LA <- read.delim("sig_LH_LA.txt", sep = "\t")
head(sig_LH_LA)
```

    ##                       baseMean log2FoldChange     lfcSE      stat       pvalue
    ## TRINITY_DN253_c0_g1  25635.895       1.034572 0.2158924  4.792073 1.650668e-06
    ## TRINITY_DN4175_c0_g1  2049.318      -1.849489 0.4025724 -4.594177 4.344613e-06
    ## TRINITY_DN824_c0_g1  22515.409       1.052742 0.2290497  4.596130 4.304095e-06
    ##                            padj
    ## TRINITY_DN253_c0_g1  0.01918581
    ## TRINITY_DN4175_c0_g1 0.01918581
    ## TRINITY_DN824_c0_g1  0.01918581

``` r
sig_BH_BA <- read.delim("sig_BH_BA.txt", sep = "\t")
head(sig_BH_BA)
```

    ##                        baseMean log2FoldChange     lfcSE      stat       pvalue
    ## TRINITY_DN10934_c0_g1 1121.0271     -0.6573609 0.1837711 -3.577065 3.474742e-04
    ## TRINITY_DN109_c0_g2   8979.2463      0.5815509 0.1273794  4.565501 4.983028e-06
    ## TRINITY_DN1127_c0_g1  4289.5233      0.4157550 0.1037463  4.007420 6.138565e-05
    ## TRINITY_DN11294_c0_g1  787.5860      0.5226225 0.1301763  4.014728 5.951439e-05
    ## TRINITY_DN11320_c0_g1  959.4128      0.9636948 0.2546909  3.783782 1.544634e-04
    ## TRINITY_DN11407_c0_g1  850.0212     -0.9163361 0.2709343 -3.382133 7.192519e-04
    ##                             padj
    ## TRINITY_DN10934_c0_g1 0.03192111
    ## TRINITY_DN109_c0_g2   0.00305971
    ## TRINITY_DN1127_c0_g1  0.01155689
    ## TRINITY_DN11294_c0_g1 0.01152595
    ## TRINITY_DN11320_c0_g1 0.02026460
    ## TRINITY_DN11407_c0_g1 0.04511487

``` r
sig_GH_GA <- read.delim("sig_GH_GA.txt", sep = "\t")
head(sig_GH_GA)
```

    ##                         baseMean log2FoldChange     lfcSE      stat
    ## TRINITY_DN1032_c0_g1  6995.91574      -1.539359 0.3525111 -4.366838
    ## TRINITY_DN12430_c0_g1  381.64580      -3.179922 0.7708568 -4.125178
    ## TRINITY_DN14659_c0_g1   89.58115      -3.679411 0.9416779 -3.907293
    ## TRINITY_DN14955_c0_g1  447.35252      -1.205378 0.2636935 -4.571132
    ## TRINITY_DN16066_c0_g1   65.20442      -1.981359 0.5115532 -3.873222
    ## TRINITY_DN16997_c0_g1   91.10238      -2.630511 0.6797065 -3.870068
    ##                             pvalue        padj
    ## TRINITY_DN1032_c0_g1  1.260581e-05 0.012345934
    ## TRINITY_DN12430_c0_g1 3.704482e-05 0.029478417
    ## TRINITY_DN14659_c0_g1 9.333600e-05 0.049475113
    ## TRINITY_DN14955_c0_g1 4.850967e-06 0.008823216
    ## TRINITY_DN16066_c0_g1 1.074061e-04 0.049475113
    ## TRINITY_DN16997_c0_g1 1.088048e-04 0.049475113

What are the orthogroups of the DEGs?

``` r
# give it the column to merge by
sig_LH_LA$transcript_id <- rownames(sig_LH_LA)
head(sig_LH_LA)
```

    ##                       baseMean log2FoldChange     lfcSE      stat       pvalue
    ## TRINITY_DN253_c0_g1  25635.895       1.034572 0.2158924  4.792073 1.650668e-06
    ## TRINITY_DN4175_c0_g1  2049.318      -1.849489 0.4025724 -4.594177 4.344613e-06
    ## TRINITY_DN824_c0_g1  22515.409       1.052742 0.2290497  4.596130 4.304095e-06
    ##                            padj        transcript_id
    ## TRINITY_DN253_c0_g1  0.01918581  TRINITY_DN253_c0_g1
    ## TRINITY_DN4175_c0_g1 0.01918581 TRINITY_DN4175_c0_g1
    ## TRINITY_DN824_c0_g1  0.01918581  TRINITY_DN824_c0_g1

``` r
# make dataframe of only LH_LA sig othrogroups
LH_LA_orthogroups <- merge(sig_LH_LA, tr.orthogroups, by= "transcript_id", sort = TRUE )
head(LH_LA_orthogroups)
```

    ## [1] transcript_id  baseMean       log2FoldChange lfcSE          stat          
    ## [6] pvalue         padj           Orthogroup    
    ## <0 rows> (or 0-length row.names)

``` r
# ok this is nothing because the orthogroup file is all by isoforms not genes
```

I think what I’ll need to do is just remove the \_i\# and then remove
the duplicate rows…

``` r
# I can use . as any character! 
tr.orthogroups.genes <-mapply(gsub, pattern = "_i.", replacement = "", tr.orthogroups)
tr.orthogroups.genes <- as.data.frame(tr.orthogroups.genes)
# there might be random spaces in the file??
tr.orthogroups.genes <-mapply(gsub, pattern = " ", replacement = "", tr.orthogroups.genes)
tr.orthogroups.genes <- as.data.frame(tr.orthogroups.genes)
colnames(tr.orthogroups.genes) <- c("Orthogroup", "transcript_id")
tr.orthogroups.genes.d <- unique(tr.orthogroups.genes)

tr.orthogroups.genes[1:100,]
```

    ##     Orthogroup          transcript_id
    ## 1    OG0000000  TRINITY_DN56081_c1_g2
    ## 2    OG0000000  TRINITY_DN56081_c1_g3
    ## 3    OG0000000  TRINITY_DN56081_c1_g1
    ## 4    OG0000000  TRINITY_DN56081_c0_g1
    ## 5    OG0000000  TRINITY_DN56081_c0_g1
    ## 6    OG0000001 TRINITY_DN108939_c0_g1
    ## 7    OG0000001  TRINITY_DN10826_c0_g1
    ## 8    OG0000002  TRINITY_DN17786_c0_g1
    ## 9    OG0000002  TRINITY_DN17786_c0_g1
    ## 10   OG0000003 TRINITY_DN49892_c0_g11
    ## 11   OG0000003  TRINITY_DN35058_c0_g1
    ## 12   OG0000003 TRINITY_DN114943_c0_g1
    ## 13   OG0000003  TRINITY_DN78128_c0_g1
    ## 14   OG0000004   TRINITY_DN6826_c0_g1
    ## 15   OG0000005 TRINITY_DN110841_c0_g1
    ## 16   OG0000005  TRINITY_DN46534_c0_g1
    ## 17   OG0000005  TRINITY_DN92340_c0_g1
    ## 18   OG0000006  TRINITY_DN42618_c0_g1
    ## 19   OG0000006  TRINITY_DN26758_c0_g1
    ## 20   OG0000006  TRINITY_DN72314_c0_g1
    ## 21   OG0000006  TRINITY_DN76548_c0_g1
    ## 22   OG0000006  TRINITY_DN74097_c0_g1
    ## 23   OG0000007  TRINITY_DN91752_c0_g1
    ## 24   OG0000008  TRINITY_DN35209_c0_g1
    ## 25   OG0000008  TRINITY_DN35209_c0_g2
    ## 26   OG0000009  TRINITY_DN68318_c0_g1
    ## 27   OG0000009  TRINITY_DN98384_c0_g1
    ## 28   OG0000009  TRINITY_DN79843_c0_g1
    ## 29   OG0000010  TRINITY_DN72839_c0_g1
    ## 30   OG0000011  TRINITY_DN23187_c0_g1
    ## 31   OG0000011  TRINITY_DN34527_c0_g1
    ## 32   OG0000011  TRINITY_DN42403_c0_g1
    ## 33   OG0000011  TRINITY_DN17592_c0_g1
    ## 34   OG0000011 TRINITY_DN33021_c0_g10
    ## 35   OG0000011  TRINITY_DN31699_c0_g1
    ## 36   OG0000011  TRINITY_DN31699_c0_g1
    ## 37   OG0000011  TRINITY_DN42166_c0_g1
    ## 38   OG0000011  TRINITY_DN79436_c0_g1
    ## 39   OG0000011  TRINITY_DN17106_c1_g1
    ## 40   OG0000011  TRINITY_DN60091_c0_g1
    ## 41   OG0000011  TRINITY_DN52646_c0_g1
    ## 42   OG0000011  TRINITY_DN18134_c0_g1
    ## 43   OG0000011  TRINITY_DN20561_c0_g1
    ## 44   OG0000011  TRINITY_DN61326_c0_g1
    ## 45   OG0000011  TRINITY_DN43346_c0_g1
    ## 46   OG0000011  TRINITY_DN45540_c0_g1
    ## 47   OG0000011  TRINITY_DN48062_c0_g1
    ## 48   OG0000012  TRINITY_DN82227_c0_g1
    ## 49   OG0000012  TRINITY_DN82227_c0_g1
    ## 50   OG0000012  TRINITY_DN45193_c0_g1
    ## 51   OG0000012  TRINITY_DN45193_c0_g1
    ## 52   OG0000012  TRINITY_DN51900_c0_g1
    ## 53   OG0000013  TRINITY_DN64607_c0_g1
    ## 54   OG0000013  TRINITY_DN40302_c0_g1
    ## 55   OG0000013  TRINITY_DN40302_c0_g1
    ## 56   OG0000014   TRINITY_DN3740_c0_g1
    ## 57   OG0000015 TRINITY_DN110830_c0_g1
    ## 58   OG0000017  TRINITY_DN15019_c0_g1
    ## 59   OG0000017  TRINITY_DN15019_c0_g1
    ## 60   OG0000017 TRINITY_DN109560_c0_g1
    ## 61   OG0000018  TRINITY_DN64986_c0_g1
    ## 62   OG0000018 TRINITY_DN106894_c0_g1
    ## 63   OG0000018  TRINITY_DN64986_c0_g1
    ## 64   OG0000019 TRINITY_DN102067_c0_g1
    ## 65   OG0000019 TRINITY_DN115068_c0_g1
    ## 66   OG0000020  TRINITY_DN65263_c0_g1
    ## 67   OG0000020  TRINITY_DN65263_c0_g1
    ## 68   OG0000020  TRINITY_DN65263_c0_g1
    ## 69   OG0000020  TRINITY_DN78837_c0_g1
    ## 70   OG0000020  TRINITY_DN59168_c0_g1
    ## 71   OG0000020   TRINITY_DN7000_c0_g2
    ## 72   OG0000020  TRINITY_DN51152_c0_g1
    ## 73   OG0000020  TRINITY_DN72197_c0_g1
    ## 74   OG0000020  TRINITY_DN33746_c0_g1
    ## 75   OG0000020  TRINITY_DN78842_c0_g1
    ## 76   OG0000020  TRINITY_DN34255_c0_g1
    ## 77   OG0000020  TRINITY_DN54484_c0_g1
    ## 78   OG0000020  TRINITY_DN62929_c0_g1
    ## 79   OG0000020  TRINITY_DN46633_c0_g1
    ## 80   OG0000021  TRINITY_DN35207_c0_g2
    ## 81   OG0000021  TRINITY_DN33378_c0_g2
    ## 82   OG0000021  TRINITY_DN33378_c0_g2
    ## 83   OG0000021  TRINITY_DN42287_c0_g1
    ## 84   OG0000021  TRINITY_DN42287_c0_g1
    ## 85   OG0000021 TRINITY_DN104561_c0_g1
    ## 86   OG0000021  TRINITY_DN84167_c0_g1
    ## 87   OG0000021 TRINITY_DN24349_c0_g10
    ## 88   OG0000021  TRINITY_DN24349_c0_g1
    ## 89   OG0000021  TRINITY_DN24349_c0_g1
    ## 90   OG0000021  TRINITY_DN40403_c0_g1
    ## 91   OG0000021  TRINITY_DN75469_c0_g1
    ## 92   OG0000021 TRINITY_DN107161_c0_g1
    ## 93   OG0000021  TRINITY_DN73480_c0_g1
    ## 94   OG0000022  TRINITY_DN59104_c0_g1
    ## 95   OG0000023  TRINITY_DN46228_c0_g1
    ## 96   OG0000023  TRINITY_DN63184_c0_g1
    ## 97   OG0000023  TRINITY_DN67584_c0_g1
    ## 98   OG0000023  TRINITY_DN53859_c0_g1
    ## 99   OG0000023  TRINITY_DN10093_c0_g1
    ## 100  OG0000023  TRINITY_DN10093_c0_g1

``` r
tr.orthogroups.genes.d[1:100,]
```

    ##     Orthogroup          transcript_id
    ## 1    OG0000000  TRINITY_DN56081_c1_g2
    ## 2    OG0000000  TRINITY_DN56081_c1_g3
    ## 3    OG0000000  TRINITY_DN56081_c1_g1
    ## 4    OG0000000  TRINITY_DN56081_c0_g1
    ## 6    OG0000001 TRINITY_DN108939_c0_g1
    ## 7    OG0000001  TRINITY_DN10826_c0_g1
    ## 8    OG0000002  TRINITY_DN17786_c0_g1
    ## 10   OG0000003 TRINITY_DN49892_c0_g11
    ## 11   OG0000003  TRINITY_DN35058_c0_g1
    ## 12   OG0000003 TRINITY_DN114943_c0_g1
    ## 13   OG0000003  TRINITY_DN78128_c0_g1
    ## 14   OG0000004   TRINITY_DN6826_c0_g1
    ## 15   OG0000005 TRINITY_DN110841_c0_g1
    ## 16   OG0000005  TRINITY_DN46534_c0_g1
    ## 17   OG0000005  TRINITY_DN92340_c0_g1
    ## 18   OG0000006  TRINITY_DN42618_c0_g1
    ## 19   OG0000006  TRINITY_DN26758_c0_g1
    ## 20   OG0000006  TRINITY_DN72314_c0_g1
    ## 21   OG0000006  TRINITY_DN76548_c0_g1
    ## 22   OG0000006  TRINITY_DN74097_c0_g1
    ## 23   OG0000007  TRINITY_DN91752_c0_g1
    ## 24   OG0000008  TRINITY_DN35209_c0_g1
    ## 25   OG0000008  TRINITY_DN35209_c0_g2
    ## 26   OG0000009  TRINITY_DN68318_c0_g1
    ## 27   OG0000009  TRINITY_DN98384_c0_g1
    ## 28   OG0000009  TRINITY_DN79843_c0_g1
    ## 29   OG0000010  TRINITY_DN72839_c0_g1
    ## 30   OG0000011  TRINITY_DN23187_c0_g1
    ## 31   OG0000011  TRINITY_DN34527_c0_g1
    ## 32   OG0000011  TRINITY_DN42403_c0_g1
    ## 33   OG0000011  TRINITY_DN17592_c0_g1
    ## 34   OG0000011 TRINITY_DN33021_c0_g10
    ## 35   OG0000011  TRINITY_DN31699_c0_g1
    ## 37   OG0000011  TRINITY_DN42166_c0_g1
    ## 38   OG0000011  TRINITY_DN79436_c0_g1
    ## 39   OG0000011  TRINITY_DN17106_c1_g1
    ## 40   OG0000011  TRINITY_DN60091_c0_g1
    ## 41   OG0000011  TRINITY_DN52646_c0_g1
    ## 42   OG0000011  TRINITY_DN18134_c0_g1
    ## 43   OG0000011  TRINITY_DN20561_c0_g1
    ## 44   OG0000011  TRINITY_DN61326_c0_g1
    ## 45   OG0000011  TRINITY_DN43346_c0_g1
    ## 46   OG0000011  TRINITY_DN45540_c0_g1
    ## 47   OG0000011  TRINITY_DN48062_c0_g1
    ## 48   OG0000012  TRINITY_DN82227_c0_g1
    ## 50   OG0000012  TRINITY_DN45193_c0_g1
    ## 52   OG0000012  TRINITY_DN51900_c0_g1
    ## 53   OG0000013  TRINITY_DN64607_c0_g1
    ## 54   OG0000013  TRINITY_DN40302_c0_g1
    ## 56   OG0000014   TRINITY_DN3740_c0_g1
    ## 57   OG0000015 TRINITY_DN110830_c0_g1
    ## 58   OG0000017  TRINITY_DN15019_c0_g1
    ## 60   OG0000017 TRINITY_DN109560_c0_g1
    ## 61   OG0000018  TRINITY_DN64986_c0_g1
    ## 62   OG0000018 TRINITY_DN106894_c0_g1
    ## 64   OG0000019 TRINITY_DN102067_c0_g1
    ## 65   OG0000019 TRINITY_DN115068_c0_g1
    ## 66   OG0000020  TRINITY_DN65263_c0_g1
    ## 69   OG0000020  TRINITY_DN78837_c0_g1
    ## 70   OG0000020  TRINITY_DN59168_c0_g1
    ## 71   OG0000020   TRINITY_DN7000_c0_g2
    ## 72   OG0000020  TRINITY_DN51152_c0_g1
    ## 73   OG0000020  TRINITY_DN72197_c0_g1
    ## 74   OG0000020  TRINITY_DN33746_c0_g1
    ## 75   OG0000020  TRINITY_DN78842_c0_g1
    ## 76   OG0000020  TRINITY_DN34255_c0_g1
    ## 77   OG0000020  TRINITY_DN54484_c0_g1
    ## 78   OG0000020  TRINITY_DN62929_c0_g1
    ## 79   OG0000020  TRINITY_DN46633_c0_g1
    ## 80   OG0000021  TRINITY_DN35207_c0_g2
    ## 81   OG0000021  TRINITY_DN33378_c0_g2
    ## 83   OG0000021  TRINITY_DN42287_c0_g1
    ## 85   OG0000021 TRINITY_DN104561_c0_g1
    ## 86   OG0000021  TRINITY_DN84167_c0_g1
    ## 87   OG0000021 TRINITY_DN24349_c0_g10
    ## 88   OG0000021  TRINITY_DN24349_c0_g1
    ## 90   OG0000021  TRINITY_DN40403_c0_g1
    ## 91   OG0000021  TRINITY_DN75469_c0_g1
    ## 92   OG0000021 TRINITY_DN107161_c0_g1
    ## 93   OG0000021  TRINITY_DN73480_c0_g1
    ## 94   OG0000022  TRINITY_DN59104_c0_g1
    ## 95   OG0000023  TRINITY_DN46228_c0_g1
    ## 96   OG0000023  TRINITY_DN63184_c0_g1
    ## 97   OG0000023  TRINITY_DN67584_c0_g1
    ## 98   OG0000023  TRINITY_DN53859_c0_g1
    ## 99   OG0000023  TRINITY_DN10093_c0_g1
    ## 101  OG0000023  TRINITY_DN40415_c0_g1
    ## 102  OG0000023 TRINITY_DN102989_c0_g1
    ## 103  OG0000024  TRINITY_DN45722_c0_g1
    ## 104  OG0000024 TRINITY_DN35920_c0_g17
    ## 105  OG0000024  TRINITY_DN73627_c0_g1
    ## 106  OG0000024   TRINITY_DN8604_c0_g1
    ## 107  OG0000025  TRINITY_DN36347_c0_g1
    ## 110  OG0000025  TRINITY_DN29599_c0_g1
    ## 111  OG0000026  TRINITY_DN76575_c0_g1
    ## 112  OG0000026 TRINITY_DN76575_c0_g15
    ## 116  OG0000026 TRINITY_DN76575_c0_g11
    ## 117  OG0000026 TRINITY_DN76575_c0_g14
    ## 120  OG0000026 TRINITY_DN76575_c0_g13
    ## 123  OG0000027  TRINITY_DN72186_c0_g1

``` r
sig_LH_LA$transcript_id <- rownames(sig_LH_LA)
LH_LA_orthogroups <- merge(sig_LH_LA, tr.orthogroups.genes.d, by= "transcript_id", sort = TRUE )
head(LH_LA_orthogroups)
```

    ##          transcript_id  baseMean log2FoldChange     lfcSE      stat
    ## 1  TRINITY_DN253_c0_g1 25635.895       1.034572 0.2158924  4.792073
    ## 2 TRINITY_DN4175_c0_g1  2049.318      -1.849489 0.4025724 -4.594177
    ## 3  TRINITY_DN824_c0_g1 22515.409       1.052742 0.2290497  4.596130
    ##         pvalue       padj Orthogroup
    ## 1 1.650668e-06 0.01918581  OG0001022
    ## 2 4.344613e-06 0.01918581  OG0009774
    ## 3 4.304095e-06 0.01918581  OG0001029

``` r
sig_BH_BA$transcript_id <- rownames(sig_BH_BA)
sig_BH_BA_orthogroups <- merge(sig_BH_BA, tr.orthogroups.genes.d, by= "transcript_id", sort = TRUE )
head(sig_BH_BA_orthogroups)
```

    ##           transcript_id  baseMean log2FoldChange     lfcSE      stat
    ## 1   TRINITY_DN109_c0_g2 8979.2463      0.5815509 0.1273794  4.565501
    ## 2 TRINITY_DN10934_c0_g1 1121.0271     -0.6573609 0.1837711 -3.577065
    ## 3  TRINITY_DN1127_c0_g1 4289.5233      0.4157550 0.1037463  4.007420
    ## 4 TRINITY_DN11294_c0_g1  787.5860      0.5226225 0.1301763  4.014728
    ## 5 TRINITY_DN11320_c0_g1  959.4128      0.9636948 0.2546909  3.783782
    ## 6 TRINITY_DN11407_c0_g1  850.0212     -0.9163361 0.2709343 -3.382133
    ##         pvalue       padj Orthogroup
    ## 1 4.983028e-06 0.00305971  OG0000133
    ## 2 3.474742e-04 0.03192111  OG0007619
    ## 3 6.138565e-05 0.01155689  OG0003277
    ## 4 5.951439e-05 0.01152595  OG0004847
    ## 5 1.544634e-04 0.02026460  OG0004871
    ## 6 7.192519e-04 0.04511487  OG0004876

``` r
sig_GH_GA$transcript_id <- rownames(sig_GH_GA)
sig_GH_GA_orthogroups <- merge(sig_GH_GA, tr.orthogroups.genes.d, by= "transcript_id", sort = TRUE )
head(sig_GH_GA_orthogroups)
```

    ##           transcript_id   baseMean log2FoldChange     lfcSE      stat
    ## 1  TRINITY_DN1032_c0_g1 6995.91574      -1.539359 0.3525111 -4.366838
    ## 2 TRINITY_DN12430_c0_g1  381.64580      -3.179922 0.7708568 -4.125178
    ## 3 TRINITY_DN14659_c0_g1   89.58115      -3.679411 0.9416779 -3.907293
    ## 4 TRINITY_DN14955_c0_g1  447.35252      -1.205378 0.2636935 -4.571132
    ## 5 TRINITY_DN16066_c0_g1   65.20442      -1.981359 0.5115532 -3.873222
    ## 6 TRINITY_DN16997_c0_g1   91.10238      -2.630511 0.6797065 -3.870068
    ##         pvalue        padj Orthogroup
    ## 1 1.260581e-05 0.012345934  OG0003197
    ## 2 3.704482e-05 0.029478417  OG0007708
    ## 3 9.333600e-05 0.049475113  OG0011063
    ## 4 4.850967e-06 0.008823216  OG0005154
    ## 5 1.074061e-04 0.049475113  OG0011101
    ## 6 1.088048e-04 0.049475113  OG0005235

remove the stats from the DEG lists

``` r
sig_GH_GA_orthogroup <- sig_GH_GA_orthogroups %>% select(1,8)
head(sig_GH_GA_orthogroup) #29
```

    ##           transcript_id Orthogroup
    ## 1  TRINITY_DN1032_c0_g1  OG0003197
    ## 2 TRINITY_DN12430_c0_g1  OG0007708
    ## 3 TRINITY_DN14659_c0_g1  OG0011063
    ## 4 TRINITY_DN14955_c0_g1  OG0005154
    ## 5 TRINITY_DN16066_c0_g1  OG0011101
    ## 6 TRINITY_DN16997_c0_g1  OG0005235

``` r
sig_BH_BA_orthogroup <- sig_BH_BA_orthogroups %>% select(1,8)
head(sig_BH_BA_orthogroup)
```

    ##           transcript_id Orthogroup
    ## 1   TRINITY_DN109_c0_g2  OG0000133
    ## 2 TRINITY_DN10934_c0_g1  OG0007619
    ## 3  TRINITY_DN1127_c0_g1  OG0003277
    ## 4 TRINITY_DN11294_c0_g1  OG0004847
    ## 5 TRINITY_DN11320_c0_g1  OG0004871
    ## 6 TRINITY_DN11407_c0_g1  OG0004876

``` r
dim(sig_BH_BA_orthogroup) #138 DEGS 
```

    ## [1] 138   2

``` r
LH_LA_orthogroup <- LH_LA_orthogroups %>% select(1,8)
head(LH_LA_orthogroup) #3
```

    ##          transcript_id Orthogroup
    ## 1  TRINITY_DN253_c0_g1  OG0001022
    ## 2 TRINITY_DN4175_c0_g1  OG0009774
    ## 3  TRINITY_DN824_c0_g1  OG0001029

Look for the orthogroups from the DEG lists in the annotation file
that’s merged with orthogroups.

``` r
# blastula comparison 
Borthogroups <- sig_BH_BA_orthogroup$Orthogroup

Blastula_temp_DEG_annot <- subset(annot.ortho.SPU, Orthogroup %in% Borthogroups)
dim(Blastula_temp_DEG_annot) # 140 
```

    ## [1] 140   6

``` r
# this is looking good that there are actually more, because there should be a few orthogroups with 2 genes
Blastula_temp_DEG_annot[1:10,]
```

    ##          spu_id Orthogroup                                 family_member
    ## 11   SPU_000021  OG0005252          HSCB_C superfamily, DnaJ superfamily
    ## 91   SPU_000188  OG0004572                                          <NA>
    ## 221  SPU_000478  OG0009714                                          <NA>
    ## 352  SPU_000764  OG0011406                                          <NA>
    ## 380  SPU_000832  OG0010152                                        myosin
    ## 723  SPU_001583  OG0003461 alpha-crystallin superfamily, SGS superfamily
    ## 818  SPU_001799  OG0007001              C2 superfamily, FerI superfamily
    ## 1057 SPU_002333  OG0001100                                          <NA>
    ## 1100 SPU_002431  OG0010533                                          ABCB
    ## 1118 SPU_002476  OG0004880                                          <NA>
    ##      common_name
    ## 11         HscbL
    ## 91          <NA>
    ## 221       Glod4L
    ## 352       Sarmr4
    ## 380          MyV
    ## 723       Sugt1L
    ## 818         Otof
    ## 1057      Cnpy2h
    ## 1100      Abcb1b
    ## 1118        Hype
    ##                                                              synonyms
    ## 11      similar to co-chaperone protein HscB, mitochondrial precursor
    ## 91                                                               <NA>
    ## 221                                    glyoxalase domain containing 4
    ## 352                                                    Sarm-related 4
    ## 380                                                          myosin V
    ## 723        SGT1, suppressor of G2 allele of SKP1 (S. cerevisiae)-like
    ## 818  otoferlin, brain otoferlin long isoform, otoferlin isoform CRA_b
    ## 1057                                                 canopy 2 homolog
    ## 1100                                                   P-glycoprotein
    ## 1118                                 Huntingtin interacting protein E
    ##      best_genbank_hit
    ## 11         ACI66807.1
    ## 91               <NA>
    ## 221    NP_001014249.1
    ## 352              <NA>
    ## 380          47550963
    ## 723    NP_001039668.1
    ## 818       NP_919224.1
    ## 1057      NP_064337.1
    ## 1100        XP_418707
    ## 1118       NP_009007

``` r
# gastrula comparison
Gorthogroups <- sig_GH_GA_orthogroup$Orthogroup
Gastrula_temp_DEG_annot <- subset(annot.ortho.SPU, Orthogroup %in% Gorthogroups)
Gastrula_temp_DEG_annot[1:10,]
```

    ##          spu_id Orthogroup             family_member common_name
    ## 517  SPU_001138  OG0006664        carbonic anhydrase    Cara14LA
    ## 762  SPU_001669  OG0003991 HCO3_cotransp superfamily     Slc4a11
    ## 784  SPU_001722  OG0006655     prefoldin superfamily        Uri1
    ## 956  SPU_002120  OG0011063                      <NA>        <NA>
    ## 1847 SPU_004105  OG0008208                      <NA>        <NA>
    ## 1954 SPU_004345  OG0009117                      <NA>     Fabp2L2
    ## 2280 SPU_005005  OG0003954                      <NA>   Uhrf1bp1L
    ## 2737 SPU_006003  OG0009137                      <NA>        <NA>
    ## 3770 SPU_008319  OG0007092                      <NA>        <NA>
    ## 3890 SPU_008575  OG0009250                      <NA>       Pdia4
    ##                                                                          synonyms
    ## 517                                                                          <NA>
    ## 762      solute carrier family 4, sodium bicarbonate cotransporter-like member 11
    ## 784  unconventional prefoldin RPB5 interactor 1, uri, NNX3, PPP1R19, Rmp, C19orf2
    ## 956                                                                          <NA>
    ## 1847                                                                         <NA>
    ## 1954                                            fatty acid binding protein 2-like
    ## 2280                                                 UHRF1 binding protein 1-like
    ## 2737                                                                         <NA>
    ## 3770                                                                         <NA>
    ## 3890            Protein disulfide-isomerase A4 precursor (Protein ERp-72) (ERp72)
    ##      best_genbank_hit
    ## 517          AAH46995
    ## 762        EDL28289.1
    ## 784      XP_002735303
    ## 956              <NA>
    ## 1847             <NA>
    ## 1954     XP_005495341
    ## 2280      XP_416170.1
    ## 2737             <NA>
    ## 3770             <NA>
    ## 3890           119531

``` r
dim(Gastrula_temp_DEG_annot)
```

    ## [1] 29  6

``` r
# larvae comparison
Lorthogroups <- LH_LA_orthogroup$Orthogroup
Larvae_temp_DEG_annot <- subset(annot.ortho.SPU, Orthogroup %in% Lorthogroups)
head(Larvae_temp_DEG_annot)
```

    ##          spu_id Orthogroup family_member common_name
    ## 906  SPU_002003  OG0001022          <NA>   Hnrpa2B1L
    ## 2194 SPU_004836  OG0001029          <NA>    Ewsr1l_1
    ## 9974 SPU_022236  OG0009774        HSP40B      Hsp40B
    ##                                                      synonyms best_genbank_hit
    ## 906  heterogeneous nuclear ribonucleoprotein A2/B1 isoform A2        NP_002128
    ## 2194            Ewing sarcoma breakpoint region 1 isoform EWS             <NA>
    ## 9974                                                     DNAJ         AAH84307

``` r
write.csv(Blastula_temp_DEG_annot,"Blastula_temp_DEG_annot.csv",col.names=TRUE,row.names=TRUE,sep=",")
```

    ## Warning in write.csv(Blastula_temp_DEG_annot, "Blastula_temp_DEG_annot.csv", :
    ## attempt to set 'col.names' ignored

    ## Warning in write.csv(Blastula_temp_DEG_annot, "Blastula_temp_DEG_annot.csv", :
    ## attempt to set 'sep' ignored

``` r
write.csv(Gastrula_temp_DEG_annot,"Gastrula_temp_DEG_annot.csv",col.names=TRUE,row.names=TRUE,sep=",")
```

    ## Warning in write.csv(Gastrula_temp_DEG_annot, "Gastrula_temp_DEG_annot.csv", :
    ## attempt to set 'col.names' ignored

    ## Warning in write.csv(Gastrula_temp_DEG_annot, "Gastrula_temp_DEG_annot.csv", :
    ## attempt to set 'sep' ignored

``` r
write.csv(Larvae_temp_DEG_annot,"Larvae_temp_DEG_annot.csv",col.names=TRUE,row.names=TRUE,sep=",")
```

    ## Warning in write.csv(Larvae_temp_DEG_annot, "Larvae_temp_DEG_annot.csv", :
    ## attempt to set 'col.names' ignored

    ## Warning in write.csv(Larvae_temp_DEG_annot, "Larvae_temp_DEG_annot.csv", :
    ## attempt to set 'sep' ignored

GRN some of the gene names had / in them to delimiate two gene names I
think pertaining to the same gene I separated them out into separate
lines

``` r
GRN <- read.csv("journal.pbio.1002391.s013.GRN.csv", header = TRUE)
head(GRN)
```

    ##   Gene.Name    Signal Cluster.Assignment..L.v.. Membership.Score..L.v..
    ## 1    Blimp1 Conserved                         6               0.5031831
    ## 2      Bmp2 Conserved                         6               0.9951625
    ## 3       Bra Conserved                         6               0.9805601
    ## 4    Endo16 Conserved                         6               0.9149388
    ## 5   Chordin Conserved                         6               0.9366175
    ## 6    Cycpln Conserved                         6               0.9659888
    ##   Cluster.Assignment..H.t.. Membership.Score..H.t.. Cluster.Assignment..H.e..
    ## 1                         6               0.9854622                         6
    ## 2                         6               0.9775867                         6
    ## 3                         6               0.9911317                         6
    ## 4                         6               0.9726211                         6
    ## 5                         6               0.9915541                         6
    ## 6                         6               0.6135961                         6
    ##   Membership.Score..H.e..
    ## 1               0.9376113
    ## 2               0.9685458
    ## 3               0.9819882
    ## 4               0.9803611
    ## 5               0.9718851
    ## 6               0.9876851

``` r
dim(GRN) #126
```

    ## [1] 126   8

Hopefully these are the same common names as in the annotation file…

Start with Blastula because it has the most…

``` r
GRN.genes <- GRN$Gene.Name

BlastulaGRN <- subset(Blastula_temp_DEG_annot, common_name %in% GRN.genes)
# none.

# what about the whole annotation?

GRNs <- subset(annot.ortho.SPU, common_name %in% GRN.genes)
print(GRNs) 
```

    ##           spu_id Orthogroup                   family_member common_name
    ## 205   SPU_000438  OG0005119                            <NA>       p58-b
    ## 206   SPU_000439  OG0000526                            <NA>       p58-a
    ## 454   SPU_000975  OG0012269                        Forkhead        FoxF
    ## 938   SPU_002088  OG0009280                   MSP130 family      Msp130
    ## 963   SPU_002129  OG0002162                    homeobox-NKL         Not
    ## 1185  SPU_002592  OG0004855                    homeobox-NKL         Emx
    ## 1191  SPU_002603  OG0009832                         sox-hmg        SoxC
    ## 1207  SPU_002634  OG0011416                   homeobox-HOXL        Hox7
    ## 1297  SPU_002815  OG0007600                    homeobox-NKL         Dlx
    ## 1316  SPU_002874  OG0000751         Ets, SAM/Pointed domain      Ets1/2
    ## 1443  SPU_003166  OG0004343                       Myc, bHLH         Myc
    ## 1858  SPU_004136  OG0012251                            <NA>         P19
    ## 2061  SPU_004551  OG0005903                        Forkhead        FoxB
    ## 2083  SPU_004599  OG0007858                    homeobox-PRD       Pitx2
    ## 2246  SPU_004924  OG0008773                             WNT       Wnt10
    ## 2272  SPU_004983  OG0006673                         chordin     Chordin
    ## 2730  SPU_005990  OG0005092          Spicule Matrix Protein        Sm29
    ## 2944  SPU_006462  OG0003225                             GCM         Gcm
    ## 3038  SPU_006676  OG0007507                        Forkhead        FoxA
    ## 3072  SPU_006753  OG0006475                          E2F/DP        E2f3
    ## 3404  SPU_007484  OG0004537                            <NA>      Cycpln
    ## 3587  SPU_007882  OG0009679                            <NA>    C-lectin
    ## 3680  SPU_008117  OG0010002        nuclear hormone receptor        Shr2
    ## 4429  SPU_009771  OG0005382                        Forkhead        FoxG
    ## 4490  SPU_009911  OG0010829 Transforming growth factor beta       Lefty
    ## 4678  SPU_010351  OG0004607                   homeobox-TALE        IrxA
    ## 4697  SPU_010403  OG0010228                        Forkhead        FoxY
    ## 4706  SPU_010424  OG0003558                    homeobox-PRD         Otx
    ## 4790  SPU_010635  OG0010319                GATA zinc finger       Gatae
    ## 4974  SPU_011029  OG0003133     Ras superfamily, Rho family        RhoA
    ## 4978  SPU_011038  OG0006375                     Reeler-like      Endo16
    ## 4991  SPU_011064  OG0006643 Transforming growth factor beta       Nodal
    ## 5020  SPU_011130  OG0004987                             WNT       Wnt16
    ## 5272  SPU_011756  OG0004886                             WNT        Wnt1
    ## 5489  SPU_012238  OG0011437   Hint (Hedgehog/Intein) domain          Hh
    ## 5497  SPU_012253  OG0003098                   homeobox-HOXL         Eve
    ## 5605  SPU_012491  OG0008153                    homeobox-NKL         Nk1
    ## 5621  SPU_012518  OG0001798              carbonic anhydrase     Cara7LA
    ## 5847  SPU_013015  OG0011242                           T-box         Bra
    ## 6334  SPU_014131  OG0006595                     Notch/Lin12       Notch
    ## 6997  SPU_015640  OG0005540                         zf-C2H2    ScratchX
    ## 7161  SPU_015982  OG0010470                    homeobox-PRD         Gsc
    ## 7230  SPU_016128  OG0005164                             DSL       Delta
    ## 7556  SPU_016881  OG0008519                         sox-hmg        SoxE
    ## 8105  SPU_018126  OG0007463                   homeobox-TALE        Tgif
    ## 8271  SPU_018483  OG0001475                             Ets         Erg
    ## 8416  SPU_018811  OG0009254         Spicule Matrix Proteins        Sm50
    ## 8417  SPU_018813  OG0003551         Spicule Matrix Proteins        Sm37
    ## 8462  SPU_018908  OG0010862                   homeobox-SINE        Six3
    ## 8505  SPU_019002  OG0001074                        Forkhead     FoxQ2_1
    ## 9125  SPU_020371  OG0006735                             WNT        Wnt8
    ## 9214  SPU_020565  OG0005721                    homeobox-NKL        Msxl
    ## 9697  SPU_021608  OG0009524                            bHLH        HesC
    ## 10526 SPU_023463  OG0005124                             WNT        Wnt4
    ## 10839 SPU_024139  OG0008813                        Forkhead        FoxC
    ## 10861 SPU_024189  OG0004057                         zf-C2H2         Sp5
    ## 11167 SPU_024903  OG0010291                             Ets         Ese
    ## 11245 SPU_025068  OG0001729                     tetraspanin   Ttrspn_19
    ## 11345 SPU_025302  OG0007582                    homeobox-PRD        Alx1
    ## 11468 SPU_025584  OG0003814                           T-box         Tbr
    ## 11724 SPU_026099  OG0008735                   homeobox-HOXL         Lox
    ## 11812 SPU_026277  OG0004631                             WNT        Wnt5
    ## 11957 SPU_026620  OG0007778                SIP1 superfamily        Sip1
    ## 12120 SPU_027015  OG0005470                         zf-GATA       GataC
    ## 12202 SPU_027235  OG0003129             zf-C2H2, SET domain      Blimp1
    ## 12718 SPU_028395  OG0003727                            <NA>        Pks2
    ## 12760 SPU_028479  OG0003942                             Ets         Tel
    ## 12845 SPU_028698  OG0002555                        Forkhead      FoxO_1
    ## 12983 SPU_030148  OG0007686                PDGF/VEGF domain       Vegf3
    ##                                                                                        synonyms
    ## 205                                                   P58-B biomineralization protein, Hypp_302
    ## 206                                                     P58-A biomineralization protein, FcgbpL
    ## 454                                                                              forkhead box F
    ## 938                                              mesenchyme specific protein, 130 kDa; Msp130_1
    ## 963                                                          notochord homeobox-like, Noto, flh
    ## 1185                                                           empty spiracles homeobox, Emx1/2
    ## 1191                                      SRY (sex determining region Y)-box C, Sox4/11/12-like
    ## 1207                                                                         homeobox 7, Hox7_1
    ## 1297                                                                  distal-less homeobox, Dll
    ## 1316          v-ets erythroblastosis virus E26 oncogene homolog 1/2-like, pointed (pnt) homolog
    ## 1443                                myelocytomatosis viral related oncogene, Myc/Mycn/Mycl-like
    ## 1858                                  biomineralization protein SpP19, tooth matrix protein P19
    ## 2061                                                                             forkhead box B
    ## 2083                                                                  paired-like homeodomain 2
    ## 2246                                      wingless-type MMTV integration site family, member 10
    ## 2272                                            chordin, short gastrulation (sog) homolog, Chrd
    ## 2730                                                                                       <NA>
    ## 2944                                                   glial cells missing homolog, Gcm1/2-like
    ## 3038                                                                  forkhead box A, fkh, HNF3
    ## 3072                                E2F transcription factor 3-like, E2f1/2/3-like, SpE2F3, E2F
    ## 3404                                                                                       <NA>
    ## 3587                                                                                       <NA>
    ## 3680   steroid hormone receptor 2, SpSHR2, nuclear receptor subfamily 2 group C-like, Nr2c-like
    ## 4429                                                         Forkhead box G, Bf1, BF-1,\\nFoxG1
    ## 4490                                               Left right determination factor, lft,antivin
    ## 4678                                                                        iroquois homeobox A
    ## 4697                                                          forkhead box Y, forkhead\\nC-like
    ## 4706                                                                      orthodenticle homolog
    ## 4790                                  GATA binding protein E, SpGatae, SpGata-e, Gata4/5/6-like
    ## 4974                                                            RhoA/B/C;  Rho1, RAS Homology A
    ## 4978                                                                                     Endo16
    ## 4991                                                                              nodal homolog
    ## 5020                                      wingless-type MMTV integration site family, member 16
    ## 5272                                       wingless-type MMTV integration site family, member 1
    ## 5489                                                                      hedgehog (hh) homolog
    ## 5497                                                                 even-skipped homeobox, Evx
    ## 5605          NK1 homeobox, NK1 transcription factor-related protein, Nkx1, Sax, slouch homolog
    ## 5621                                                                                       <NA>
    ## 5847                                       brachyury, T, TA, brachyenteron (byn) homolog, SpBra
    ## 6334                                                           neurogenic locus Notch homolog,N
    ## 6997                                                           scratch subfamily member X, z191
    ## 7161                                                                        goosecoid homeobox 
    ## 7230                            Delta (Dl) homolog, neurogenic locus protein delta, Dll1/4-like
    ## 7556                                            SRY (sex determining region Y)-box E, Sox8/9/10
    ## 8105                                                               TGFB-induced factor homeobox
    ## 8271               v-ets oncogene related (Erg)-like, Friend leukemia integration 1 (Fli1)-like
    ## 8416                                                                                       <NA>
    ## 8417                                                                spicule matrix protein SM37
    ## 8462                                   SIX homeobox 3, sine oculis-related homeobox 3/6, Six3/6
    ## 8505                                                                   forkhead box Q2 (copy 2)
    ## 9125                                       wingless-type MMTV integration site family, member 8
    ## 9214                                                                   Msh homeobox-like, Msxlx
    ## 9697                                                    hairy and enhancer of split-like C, Hes
    ## 10526                       wingless-type MMTV integration site family member 4, Wnt4_1, Wnt4_2
    ## 10839                                                                            forkhead box C
    ## 10861                                                       Sp5 transcription factor-like, z199
    ## 11167                                epithelium-specific ets factor-like, Ehf-like, Elf3/5-like
    ## 11245                                                                                      <NA>
    ## 11345             aristaless-like homeobox 1-like, Alx1(Cart1)/Alx3/Alx4 subfamily-like, Alx1_1
    ## 11468                                           T-box brain-like, Eomes-like, Tbx21-like, ske-T
    ## 11724                                             pancreas/duodenum homeobox, Pdx1, IPF1, Splox
    ## 11812                                       wingless-type MMTV integration site family member 5
    ## 11957                                    survival of motor neuron protein interacting protein 1
    ## 12120                                           GATA binding protein C, SpGATAc, Gata1/2/3-like
    ## 12202            B lymphocyte induced maturation protein-like, PRDM1-like, blimp1/krox, SpKrox1
    ## 12718                                                                                      <NA>
    ## 12760                  ets variant 6 (TEL oncogene) and 7-like (partial), Etv6/7-like, Yan-like
    ## 12845                                                                                      <NA>
    ## 12983 vascular endothelial growth factor-like 3, PDGF- and VEGF-related factor 3 (Pvf3) homolog
    ##       best_genbank_hit
    ## 205          XP_799905
    ## 206    \\tXP_003724499
    ## 454       NP_001158437
    ## 938       NP_001116986
    ## 963         BAD91047.1
    ## 1185   \\tNP_001179152
    ## 1191      NP_001158501
    ## 1207                  
    ## 1297      NP_001158371
    ## 1316                  
    ## 1443                  
    ## 1858        AAM70484.1
    ## 2061        ABX89143.1
    ## 2083          AAW51825
    ## 2246        AHY22362.1
    ## 2272      NP_001158390
    ## 2730         NP_999804
    ## 2944       NP_999826.1
    ## 3038        CAY90194.1
    ## 3072      XP_005106760
    ## 3404        AE003458.3
    ## 3587         NP_999805
    ## 3680      XP_002739506
    ## 4429          AAP79301
    ## 4490         AAS00535 
    ## 4678     \\tADW95342.1
    ## 4697        BAE45343.1
    ## 4706                  
    ## 4790      NP_001164701
    ## 4974          BAA75688
    ## 4978         NP_999684
    ## 4991          AAS00534
    ## 5020          BAD12590
    ## 5272          CAA38991
    ## 5489      XP_003439270
    ## 5497         NM_214651
    ## 5605        ADW95337.1
    ## 5621          AAH94913
    ## 5847        BAD74048.1
    ## 6334        AAB82088.1
    ## 6997                  
    ## 7161        AAR17089.1
    ## 7230          AAL71862
    ## 7556       NP_989612.1
    ## 8105      NP_001158452
    ## 8271      NP_001161529
    ## 8416         NP_999775
    ## 8417      \\tNP_999776
    ## 8462      NP_001158378
    ## 8505      XP_003725384
    ## 9125         NP_999832
    ## 9214      XP_002597186
    ## 9697      XP_003974264
    ## 10526     XP_002613927
    ## 10839         AAH46028
    ## 10861     NP_001158478
    ## 11167       ACZ65094.1
    ## 11245     NP_001020405
    ## 11345                 
    ## 11468                 
    ## 11724        AAN17337 
    ## 11812         CAA51916
    ## 11957       CAG32122.1
    ## 12120       ABX71821.1
    ## 12202       AEO92035.1
    ## 12718         63086968
    ## 12760                 
    ## 12845     NP_001009988
    ## 12983     XP_003965794

``` r
# only 69 of these
# not sure if the names are all matching up correctly 
GRN.genes                  
```

    ##   [1] "Blimp1"     "Bmp2"       "Bra"        "Endo16"     "Chordin"   
    ##   [6] "Cycpln"     "Delta"      "E2f3"       "Emx"        "Erg"       
    ##  [11] "Ese"        "FoxA"       "FoxF"       "FoxG"       "FoxY"      
    ##  [16] "GataC"      "Gatae"      "Gcm"        "hbn"        "Hex"       
    ##  [21] "Hox11"      "IrxA"       "Lefty"      "lim1"       "Lox"       
    ##  [26] "Msp130rel1" "Msxl"       "Nk1"        "Nk2-1"      "Nkx2.2"    
    ##  [31] "Nkx3.2"     "Not"        "p58-a"      "Ployksl"    "Six3"      
    ##  [36] "Sm27"       "Sm50"       "SoxE"       "Sp5"        "Tbx2-3"    
    ##  [41] "Unc4-1"     "Wnt1"       "Wnt16"      "Alx1"       "atbf1"     
    ##  [46] "egr"        "FoxN2"      "Hh"         "Hnf6"       "lasp1"     
    ##  [51] "Notch"      "P16rel1"    "Patched"    "Sip1"       "Smooth"    
    ##  [56] "Wnt10"      "cdx"        "ets4"       "FoxQ2"      "Nodal"     
    ##  [61] "Otx"        "Shr2"       "Six1"       "Tgif"       "Wnt4"      
    ##  [66] "Apobec"     "Bmp5"       "Brn1-2-4"   "C-lectin"   "Cara7LA"   
    ##  [71] "Dlx"        "Eve"        "Fgf9"       "Ficolin"    "FoxB"      
    ##  [76] "FoxC"       "Gsc"        "HesC"       "Hox7"       "Msp130"    
    ##  [81] "Myc"        "P133"       "P16rel2"    "P19"        "p58-b"     
    ##  [86] "Pax2-5-8"   "Pdgfr"      "Pks2"       "RhoA"       "ScratchX"  
    ##  [91] "Sm29"       "Sm37"       "Sm49"       "Soxb1"      "SoxC"      
    ##  [96] "Tel"        "Ttrspn_19"  "Vegf3"      "Wnt5"       "Wnt8"      
    ## [101] "Dr/Dr_1"    "Ets1/2"     "FoxO"       "Pitx2"      "Tbr"       
    ## [106] "Unvn"       "Bmp4"       "Bra_1"      "Endo16"     "Cyp"       
    ## [111] "Hox13b"     "pks1"       "z60"        "Fox3"       "FoxQ2_1"   
    ## [116] "Six2"       "Bmp6"       "Bmp7"       "Bmp8"       "Can1"      
    ## [121] "Fgf16"      "Fgf20"      "vegfrl"     "PMC1"       "net7"      
    ## [126] "FoxO_1"

``` r
GasturalaGRN <- subset(Gastrula_temp_DEG_annot, common_name %in% GRN.genes)        
# also none
```
