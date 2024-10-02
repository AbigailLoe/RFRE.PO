# RFRE.PO
Code for RFRE.PO paper

The file **wald-selection-loop.R** is the code used for finding the forward selection model, or Model 2 in the RFRE.PO paper.

**imputed-Analysis.R** combines all of the data analysis results. It requires uploading and reading in files from the imputedCat folder; alternatively these may be found in the imputedCat.zip file.

Simulation code was performed on a high performance cluster at the University of Michigan. The necessary code is found in the **SimSetting** folder. It is currently configured for a cluster, but comments can be found where users might modify code if they wish to run it on local computers.

Files under the folder **Data Application** contain all of the necessary R code for replicating the results on a local computer found in the Data Application Section of the RFRE.PO paper. The filed ImputedAnalysis.R analyzes the data using a universal C-statistic and MSE, while imputedAnalysisTimeVarying.R is for analysis using history variables, and calculting C and MSE at the start of each follow-up window. Similarly for imputedAnalysis-noHist.R, just without the history variables.
