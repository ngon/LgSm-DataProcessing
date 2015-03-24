READ ME FILE
/LgSm-DataProcessing/dataFiles/inspect outliers
-------------------------------------------------
author: Natalia Gonzales
date: 3-24-15
_________________________________________________

This folder contains lists of outliers identified using the script
/LgSm-DataProcessing/cookieTime/inspect.outliers.R

There is one list for each trait inspected. The first column contains the
row names in the phenotype file and the second column contains the mouse id. 

There is also a list containing all traits in one document; it is called
ALLtraits.outlierList. 

ppi.outliers contains a list of samples that were removed from all ppi traits
(ppi3, ppi6, ppi12, startle, and habituation - also, the logit forms of the 
first 3 variables). The procedure for identifying outliers is different for ppi
traits because ppi traits are not normal. See the script /inspect.outliers.R for
a detailed description of how this was done. 

Since 'wild' will ultimately be used as a categorical trait, no samples were
removed from the wild data set. 