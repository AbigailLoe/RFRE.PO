# RFRE.PO
Code for RFRE.PO paper

The file **wald-selection-loop.R** is the code used for finding the forward selection model, or Model 2 in the RFRE.PO paper.

**imputed-Analysis.R** combines all of the data analysis results. It requires uploading and reading in files from the imputedCat folder; alternatively these may be found in the imputedCat.zip file.

Simulation code was performed on a high performance cluster at the University of Michigan. The necessary code is found in the **SimSetting** folder. It is currently configured for a cluster, but comments can be found where users might modify code if they wish to run it on local computers.
