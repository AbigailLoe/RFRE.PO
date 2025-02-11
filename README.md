# RFRE.PO

Data application and analysis of the Azithromycin for the Prevention of COPD Cohort.

Dependency `htree` was archived by CRAN, so it must be downloaded from the CRAN archives using: 
```
install.packages("https://cran.r-project.org/src/contrib/Archive/htree/htree_2.0.0.tar.gz", type = "source", repos = NULL)
```

## Install 
```
#install.packages(devtools)
install_github("AbigailLoe/RFRE.PO")
```

## Use 
```
load("dat.rda")
ia_mod = imputed_analysis(dat, arg1, arg2, ...)
print(ia_mod)
plot(ia_mod)
```