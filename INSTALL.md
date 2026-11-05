R Notebook
================

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.23")
```

    ## 'getOption("repos")' replaces Bioconductor standard repositories, see
    ## 'help("repositories", package = "BiocManager")' for details.
    ## Replacement repositories:
    ##     CRAN: https://p3m.dev/cran/__linux__/resolute/latest

    ## Bioconductor version 3.23 (BiocManager 1.30.27), R 4.6.2 (2026-10-29)

    ## Warning: package(s) not installed when version(s) same as or greater than current; use
    ##   `force = TRUE` to re-install: 'dada2'

    ## Installation paths not writeable, unable to update packages
    ##   path: /usr/local/lib/R/library
    ##   packages:
    ##     mgcv, nlme

autre chose à installer : Pour installer le fichier : attention dans le
terminal

wget <https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip>
unzip miseqsopdata.zip

    wget -O silva_nr99_v138.2_toGenus_trainset.fa.gz "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1"
