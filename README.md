# nanotail

<!-- badges: start -->
<!-- badges: end -->

The goal of **NanoTail** is to provide a set of functions to manipulate and analyze data coming from polyA lengths estimations done using Oxford Nanopore Direct RNA sequencing and Nanopolish software. The software is still in thee development phase so all suggestions are welcome. Please also expect the code to be changed frequently, so use it with caution.

## Installation

You can install the developmental version of Nanotail with

``` r
install.packages("devtools")
devtools::install_github('smaegol/nanotail')
library(nanotail)
```

## Input data

NanoTail needs an output from `nanopolish polya` to work. It can read single output file with `read_polya_single`:

``` r
path <- "/location/of/nanopolish/output"
polya_data <- read_polya_single(path)
```

It can also read multiple samples at once and associate any metadata with them. 
Let's assume we have performed an experiment, targeting one of the polyA polymerases. 2 replicates were sequenced for control samples, and 2 replicates sequences for samples with mutant PAP, therefore after all analysis we have 4 files with `nanopolish polya` output. To read all of them at once, we can use command `read_polya_multiple` and associate metadata using samples_table data.frame:

``` r
samples_table <- data.frame(polya_path = c(path1,path2,path3,path4),
                            sample_name =c("wt1","mu1","wt2","mut2"),
                            group = c("wt","mut","wt","mut"))
polya_data_multiple <- read_polya_multiple(samples_table)
```

## Shiny App

Once data are imported they can be processed in the R environment using `nanotail` functions described below or, more convinent, the interactive Shiny app can be launched, allowing for easy exploration of obtained data. To launch the app for the data imported above:

``` r
nanoTailApp(polya_table = polya_data_multiple)

```


## Citation

Please cite NanoTail as:
Krawczyk PS et al., NanoTail - R package for exploratory analysis of Nanopore Direct RNA based polyA lengths estimations
