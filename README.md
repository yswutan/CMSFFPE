# CMSFFPE
A CMS classifier for formalin-fixed paraffin-embeded colorectal cancer tissue

## installation
devtools::install_github("yswutan/CMSFFPE")

## usage
library(CMSFFPE);
library(randomForest);
library(caret)

labels <- CMSFFPEClassifier(Expr = Expr, minPosterior = 0.5, log2transfrom=FALSE)
