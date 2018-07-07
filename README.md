ofGEM
=====

### A Meta-Analysis Approach with Filtering for Identifying Gene-Level Gene-Environment Interactions with Genetic Association Data 


This `R` package offers a gene-based meta-analysis test with filtering to 
  detect gene-environment interactions (GxE) with association 
  data, proposed by Wang et al. (2018) 
  <doi:10.1002/gepi.22115>. It first conducts a meta-filtering 
  test to filter out unpromising SNPs by combining all samples 
  in the consortia data. It then runs a test of 
  omnibus-filtering-based GxE meta-analysis (ofGEM) that 
  combines the strengths of the fixed- and random-effects 
  meta-analysis with meta-filtering. It can also analyze data 
  from multiple ethnic groups. 


### Installation
- For the stable version from [CRAN](https://cran.r-project.org/web/packages/ofGEM/index.html):
```r
install.packages('ofGEM')
```
- For the development version (requiring the `devtools` package):
```r
devtools::install_github('randel/ofGEM')
```

### Reference
Wang, J., Liu, Q., Pierce, B. L., Huo, D., Olopade, O. I., Ahsan, H., & Chen, L. S. (2018). A meta‐analysis approach with filtering for identifying gene‐level gene–environment interactions. *Genetic epidemiology*.
