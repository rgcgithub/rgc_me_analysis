######################################################################################################################
## Goal: estimate the number of genes with X or more heterozygous or homozygous pLoF carriers at some sample size Y ##
## given the number of heterozygous or homozygous pLOF carriers observed per gene at current sample sizes	    ##
## Date: Aug-8-2023                                                                                                 ##
## Authors: Jonathan Marchini, Goncalo Abecasis, Josh Backman, Daren Liu, Kathie Sun                       	    ##	
######################################################################################################################

## BetaBinomial.R contains all functions
source("BetaBinomial.R")

maxSize = 8       ## maximum number of components to test best fit model; we suggest an integer between 6-10
total_n = 985000  ## sample size used to compute input pLOF counts

df = read.table("input_df.txt", sep="\t", header=T)
## input data was randomly generated based on distribution of het/hom cumulative LoF count in RGCME

ngenes = nrow(df) #19644
## total number of protein coding genes

PerGeneLofCounts <- df$n_homAA
## PerGeneLofCounts is a vector of heterozygous (or homozygous) pLOF counts from dataset of sample size total_n
### length(PerGeneLofCounts) = the total number of protein coding genes

## To find best model fit, first call LoopThroughModels:
## this finds the ideal number of components for the mixture model
model_loop = LoopThroughModels( PerGeneLofCounts, total_n, maxSize )

## assign bestFitSize to the number of components that return the lowest (most probable) AIC or BIC
bestFitSize = 2

## Once best fit is determined, you can run FindBetaBinomialModel
model_fin <- FindBetaBinomialModel( PerGeneLofCounts, total_n, bestFitSize, refineWithOptim = T )

num_carriers <-  c(1, 5, 10, 50, 100, 1000)
## num_carriers is a vector of number of carriers to return projections for, e.g. y genes were observed with >= 1, 5, or 10 homozygous LOF carriers

Xs = c(50000,250000,500000,750000,1000000,5000000)
## Xs is a vector of projected sample sizes; we went up to 5M with the caveat that the further from your actual sample size, the lower predictive validity of the estimates

all_curves = MatrixPrediction( num_carriers, Xs, model_fin )
## returns length(num_carriers) x length(Xs) matrix representing the proportion of the total number of genes observed for each carrier count threshold (row) for each projected sample size (column)

projections = apply(all_curves, 2, function(x) x*ngenes)
## returns the number of genes expected at each carrier count threshold per projected sample size
