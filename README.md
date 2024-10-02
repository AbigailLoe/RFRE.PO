# RFRE.PO
Code for the data application and analysis of the Azithromycin for the Prevention of COPD Cohort.

Files under the folder **Data Application** contain all of the necessary R code for replicating the results on a local computer found in the Data Application Section of the RFRE.PO paper. The filed ImputedAnalysis.R analyzes the data using a universal C-statistic and MSE, while imputedAnalysisTimeVarying.R is for analysis using history variables, and calculting C and MSE at the start of each follow-up window. Similarly for imputedAnalysis-noHist.R, just without the history variables.
