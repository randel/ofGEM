ofGEM
=====

### A test for omnibus-filtering-based gene-environment interactions meta-analysis

This R package offers a gene-based meta-analysis test with filtering to detect gene-environment interactions (GxE) with association data. 
It first conducts a meta-filtering test to filter out unpromising SNPs by combining all samples in the consortia data. 
It then runs a test of omnibus-filtering-based GxE meta-analysis (ofGEM) that combines the strengths of the fixed- and 
random-effects meta-analysis with meta-filtering. It can also analyze data from multiple ethnic groups.

### Installation (requiring the `devtools` package):
```r
devtools::install_github('randel/ofGEM')
```
