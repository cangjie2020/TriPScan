
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TriPScan

<!-- badges: start -->
<!-- badges: end -->

Translation elongation is a crucial step in the translation process, and
tRNA is an indispensable component during this process. If there are
abnormalities in tRNA expression, it may lead to changes in the ribosome
translation speed during the elongation process. Based on tRNA-seq and
Ribo-seq data, we have identified differential ribosome stalling sites
caused by tRNA expression differences during translation elongation.
Furthermore, we can calculate the changes in translation efficiency of
the genes at these differential ribosome stalling sites

## Dependencies

TriPScan requires R version \>= 3.3.0 and the following packages:

    data.table,
    DESeq2,
    dplyr,
    GenomicFeatures,
    ggplot2,
    ggrepel,
    NOISeq,
    rtracklayer,
    tidyr,
    Biostrings

All dependencies are automatically installed running the code in the
next section.

## Installation

You can install the development version of TriPScan from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cangjie2020/TriPScan")
```

## Loading

To load TriPScan run

``` r
library(TriPScan)
```

## Genome annotation matrix
Extract information from the GTF annotation file and generate a matrix. This sequence annotation file has 14 columns 
and is a crucial file in this package. It will help screen different transcripts of the same gene and assist in calculating 
the codons corresponding to the E, P, and A sites for each RPF.

### input files

gtf.path：Corresponding species gtf annotation file

``` r
geneinfo <- extractGeneInfo( gtf.path = "path/to/gtf")
#> Import genomic features from the file as a GRanges object ... OK
#> Prepare the 'metadata' data frame ... OK
#> Make the TxDb object ...
#> OK
```
The reference annotation file 

``` r
head(geneinfo)
#>               gene_id            tx_name  tx_id nexon tx_len cds_len utr5_len
#>                <char>             <char>  <int> <int>  <int>   <int>    <int>
#> 1: ENSMUSG00000000001 ENSMUST00000000001  26866     9   3262    1065      141
#> 2: ENSMUSG00000000003 ENSMUST00000000003 145808     7    902     525      140
#> 3: ENSMUSG00000000003 ENSMUST00000114041 145809     6    697     414       47
#> 4: ENSMUSG00000000028 ENSMUST00000000028 124151    20   2143    1701      313
#> 5: ENSMUSG00000000028 ENSMUST00000096990 124152    18   1747    1563       61
#> 6: ENSMUSG00000000028 ENSMUST00000115585 124153     5    832     410      422
#>    utr3_len gene_name   gene_biotype transcript_biotype         protein_id
#>       <int>    <char>         <char>             <char>             <char>
#> 1:     2056     Gnai3 protein_coding     protein_coding ENSMUSP00000000001
#> 2:      237      Pbsn protein_coding     protein_coding ENSMUSP00000000003
#> 3:      236      Pbsn protein_coding     protein_coding ENSMUSP00000109675
#> 4:      129     Cdc45 protein_coding     protein_coding ENSMUSP00000000028
#> 5:      123     Cdc45 protein_coding     protein_coding ENSMUSP00000094753
#> 6:        0     Cdc45 protein_coding     protein_coding ENSMUSP00000111248
#>     chrom strand
#>    <fctr> <fctr>
#> 1:      3      -
#> 2:      X      -
#> 3:      X      -
#> 4:     16      -
#> 5:     16      -
#> 6:     16      -
```

## Identifying differential tRNA-index

``` r
AC_count <- "../Anticodon_counts_raw.txt"
condition <- list(N1 = c(3,4),N2 = c(5,6))
tRNA_res <- diff_AC(AC_count = AC_count,condition = condition,replicates = T)
```

``` r
head(tRNA_res[[1]])
#>   Anticodon  AA Codon ..5.result_mim.mmu_more20.FB1.fq.gz.unpaired_uniq.bam
#> 1       AAA Phe   TTT                                           0.002575701
#> 2       AAC Val   GTT                                           0.024378951
#> 3       AAG Leu   CTT                                           0.005665191
#> 4       AAT Ile   ATT                                           0.013916428
#> 5       ACA Cys   TGT                                           0.002425710
#> 6       ACC Gly   GGT                                           0.025603641
#>   ..5.result_mim.mmu_more20.FB2.fq.gz.unpaired_uniq.bam
#> 1                                           0.001763311
#> 2                                           0.025105997
#> 3                                           0.006533785
#> 4                                           0.014840934
#> 5                                           0.002788376
#> 6                                           0.027870850
#>   ..5.result_mim.mmu_more20.FH1.fq.gz.unpaired_uniq.bam
#> 1                                           0.002312578
#> 2                                           0.028820882
#> 3                                           0.007336867
#> 4                                           0.015528162
#> 5                                           0.001952492
#> 6                                           0.027949406
#>   ..5.result_mim.mmu_more20.FH2.fq.gz.unpaired_uniq.bam
#> 1                                           0.001908306
#> 2                                           0.028331843
#> 3                                           0.008485709
#> 4                                           0.016511510
#> 5                                           0.001995792
#> 6                                           0.027787390
#>   ..5.result_mim.mmu_more20.FL1.fq.gz.unpaired_uniq.bam
#> 1                                           0.002380757
#> 2                                           0.035323424
#> 3                                           0.009122034
#> 4                                           0.015110099
#> 5                                           0.002308543
#> 6                                           0.023285067
#>   ..5.result_mim.mmu_more20.FL2.fq.gz.unpaired_uniq.bam
#> 1                                           0.002051268
#> 2                                           0.035521336
#> 3                                           0.009132332
#> 4                                           0.012458845
#> 5                                           0.001814826
#> 6                                           0.027289718
```

``` r
head(tRNA_res[[2]])
#>      baseMean log2FoldChange      lfcSE       stat      pvalue       padj codon
#> TTT  2161.336    -0.05688517 0.19312714 -0.2945478 0.768339377 0.78864483   TTT
#> GTT 31985.649    -0.29800067 0.09687794 -3.0760427 0.002097678 0.01678143   GTT
#> CTT  8517.411    -0.19405452 0.12812507 -1.5145711 0.129881063 0.33060634   CTT
#> ATT 14895.626     0.23128915 0.11953777  1.9348625 0.053007175 0.19861475   ATT
#> TGT  2015.648    -0.04796306 0.19618271 -0.2444816 0.806857869         NA   TGT
#> GGT 26600.685     0.15006894 0.11640318  1.2892168 0.197322714 0.39464543   GGT
```

``` r
tRNA_res[[3]]
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

## Extract codon information at E, P, and A sites for each read from ribo-seq data

``` r
load("C:/Users/96453/Desktop/ribowalts/reads_psite_list.rda")
EPA_result<-EPA_res(cdna_path = "../Mus_musculus.GRCm39.cdna.all.fa",reads_psite_list = reads_psite_list,geneinfo,tx_scale ="no")
```

``` r
head(reads_psite_list[[1]])
#>            transcript  end5 psite  end3 length cds_start cds_stop
#>                <fctr> <int> <num> <int>  <int>     <num>    <num>
#> 1: ENSMUST00000070533   196   208   224     29       151     2094
#> 2: ENSMUST00000070533   196   208   224     29       151     2094
#> 3: ENSMUST00000070533   196   208   224     29       151     2094
#> 4: ENSMUST00000070533   196   208   224     29       151     2094
#> 5: ENSMUST00000070533   697   709   722     26       151     2094
#> 6: ENSMUST00000070533   933   945   962     30       151     2094
#>    psite_from_start psite_from_stop psite_region
#>               <num>           <num>       <char>
#> 1:               57           -1886          cds
#> 2:               57           -1886          cds
#> 3:               57           -1886          cds
#> 4:               57           -1886          cds
#> 5:              558           -1385          cds
#> 6:              794           -1149          cds
```

``` r
head(EPA_result[[1]])
#>              tx_name trans.site E_Codon P_Codon A_Codon
#> 1 ENSMUST00000070533        208     TTC     ACC     CCG
#> 2 ENSMUST00000070533        208     TTC     ACC     CCG
#> 3 ENSMUST00000070533        208     TTC     ACC     CCG
#> 4 ENSMUST00000070533        208     TTC     ACC     CCG
#> 5 ENSMUST00000070533        709     GGA     GCA     GAT
#> 9 ENSMUST00000070533       1264     ATT     CAG     TTC
```

## Combine the results of EPA_res and diff_AC to identify the differential peaks caused by tRNA

``` r
EPA_path <- list(S1 = list(EPA_result[[1]],EPA_result[[2]]),
                 S2 = list(EPA_result[[3]],EPA_result[[4]]))
condition <- c("Heart","Liver")
peak_pic <- tRNA_diff_peak(EPA_path = EPA_path,tRNA = tRNA_res,condition = condition,replicates = T,
          geneinfo = geneinfo,min_diff = 8,diff_tRNA_num = 5)
```

``` r
names(peak_pic)
#> [1] "diff_peak"          "Heart_VS_Liver"     "res_Heart_VS_Liver"
#> [4] "Liver_VS_Heart"     "res_Liver_VS_Heart"
```

``` r
head(peak_pic[[1]])
#>                    name ratio_Heart ratio_Liver       diff            tx_name
#>                  <char>       <num>       <num>      <num>             <char>
#> 1: ENSMUST00000000028_6  1.43712575   1.1304348  0.3066910 ENSMUST00000000028
#> 2: ENSMUST00000000080_6  1.12096774   0.7625899  0.3583778 ENSMUST00000000080
#> 3: ENSMUST00000000090_6  0.03302522   0.1664112 -0.1333860 ENSMUST00000000090
#> 4: ENSMUST00000000109_6  0.58288770   0.0000000  0.5828877 ENSMUST00000000109
#> 5: ENSMUST00000000356_6  2.53808752   2.6436498 -0.1055623 ENSMUST00000000356
#> 6: ENSMUST00000000451_6  0.48770053   1.1278953 -0.6401947 ENSMUST00000000451
#>     site A_Codon
#>    <num>  <char>
#> 1:     6     AAG
#> 2:     6     CCG
#> 3:     6     CGC
#> 4:     6     AGC
#> 5:     6     CCA
#> 6:     6     TGG
```

``` r
head(peak_pic[[3]])
#>                     name ratio_Heart ratio_Liver     diff            tx_name
#> 1  ENSMUST00000020768_71    23.43879   0.5714286 22.86737 ENSMUST00000020768
#> 2 ENSMUST00000035269_378    22.11942   2.7096774 19.40974 ENSMUST00000035269
#> 3 ENSMUST00000037491_562    17.48366   0.8800000 16.60366 ENSMUST00000037491
#> 4  ENSMUST00000110167_99    30.24879  14.8032723 15.44552 ENSMUST00000110167
#> 5 ENSMUST00000021750_322    13.90916   0.0000000 13.90916 ENSMUST00000021750
#> 6 ENSMUST00000109313_720    13.04525   0.0000000 13.04525 ENSMUST00000109313
#>   site A_Codon
#> 1   71     GTT
#> 2  378     GTG
#> 3  562     GTA
#> 4   99     GTG
#> 5  322     GTA
#> 6  720     GTC
```

``` r
peak_pic[[2]][[14]]
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" /> \##

## To calculate the translation efficiency of differentially expressed peaks (diff-peak) genes

``` r
RNA<-read.table("../CDS/RNA.txt",sep = "\t",header = T)
CDS<-read.table("../CDS/CDS.txt",sep = "\t",header = T)
condition <- c("Heart","Liver")
TE_pic<- TE_peak(RNA_count = RNA,CDS_count = CDS,gene = geneinfo,peak_pic = peak_pic ,condition = condition ,replicates = 2 ,gene_id = F)
```

``` r
head(RNA)
#>              tx_name SRR70 SRR71 SRR74  SRR75
#> 1 ENSMUST00000000001   819   616  2199   2845
#> 2 ENSMUST00000000003     0     0     0      0
#> 3 ENSMUST00000000010     1     0     0      2
#> 4 ENSMUST00000000028   132   101    74     46
#> 5 ENSMUST00000000049     8     3 89238 121118
#> 6 ENSMUST00000000058  2855  2130   179    206
head(CDS)
#>              tx_name SRR86 SRR87 SRR90 SRR91
#> 1 ENSMUST00000000001   515   486  1589  2160
#> 2 ENSMUST00000000028   116   116    18    19
#> 3 ENSMUST00000000049   176   242 66319 84344
#> 4 ENSMUST00000000058   772   609    79   131
#> 5 ENSMUST00000000080   509   455   147   212
#> 6 ENSMUST00000000087   121   123    38   102
```

``` r
head(TE_pic[[1]])
#>                         SRR86       SRR87     SRR90     SRR91
#> ENSMUST00000000001  3.5659668   4.3060758 2.9450376 2.6428407
#> ENSMUST00000000028  2.0498474   2.5783759 0.4077711 0.5913965
#> ENSMUST00000000049 46.6973152 164.7917563 1.1336933 0.9073195
#> ENSMUST00000000058  2.7980894   2.8474770 3.2821820 4.0392432
#> ENSMUST00000000080  1.8232094   1.4480030 1.5559547 1.6698457
#> ENSMUST00000000087  0.2478779   0.3607961 0.1876004 0.3920122
```

``` r
TE_pic[[3]][[14]]
```

<img src="man/figures/README-unnamed-chunk-19-1.png" width="100%" />
