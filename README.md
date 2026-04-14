# Multiplex Random Walk with Restart (RWR) Framework

This repository implements an optimized Random Walk with Restart (RWR) framework on a heterogeneous multiplex biological network, integrating gene, disease, and metabolite layers. The goal is to prioritize disease associations based on network topology and diffusion dynamics.

This is a fork with some improvements of the [MDAKM Lab](https://mdakm.ceid.upatras.gr/) RWR_MH algorithm in the framework of my thesis in patient classification using GCNNs for patients with hypercholesterolemia and diabetes.

Original code can be found [here](https://github.com/mariarevythi/RandomWalkLab) from Maria Revythi.


## Pipeline

### Data Preparation

Input networks are loaded as adjacency matrices (gene, metabolite, etc.)
Seed nodes (genes and diseases) are defined
Node naming is standardized across layers

### Multiplex Network Construction

Individual layers are combined into a supra-adjacency matrix
Inter-layer connections are controlled by parameter delta
Identity matrices are used to link corresponding nodes across layers

### Transition Matrix Computation

Column-normalized matrix used for random walk propagation
Incorporates inter-layer navigation via parameter lambda

### Random Walk with Restart (RWR)

Iterative diffusion process
Controlled by restart probability r
Converges to a steady-state proximity vector

### Disease Ranking

Nodes are ranked based on steady-state probabilities
Seed diseases are excluded from final ranking

## Performance Optimizations over the original


### Edge Pruning (Top-K Filtering / ElasticNet Filtering)

Introduction of additional functions in this implementation, an alternative TopK percentage pruning and an [ElasticNet](https://en.wikipedia.org/wiki/Elastic_net_regularization) inspired Lasso and Ridge pruning to impove disease separability and classification results

functions can be found [here](https://github.com/Feasao/RWRLab/blob/main/Random_walk/Functions/All_Functions.R)
get.supra.adj.multiplexWeightedTopK , get.supra.adj.multiplexWeightedEN

### Vectorized Transition Matrix

The original implementation used explicit loops:

```
for (j in 1:n) {
  Transition[,j] <- ...
}
```

This was replaced with a fully vectorized approach that leverages BLAS/LAPACK routines for up to 5000 times speed increase:

```
Transition <- SupraAdjacencyMatrix %*% Diagonal(x = scale_factors)
```

functions can be found [here](https://github.com/Feasao/RWRLab/blob/main/Random_walk/Functions/All_Functions.R)
get.transition.multiplex , get.transition.multiplexFAST


## Results

Results can be found:
For the tuning of the algorithm in [dataset_tuning_results](https://github.com/Feasao/RWRLab/tree/main/Random_walk/dataset_tuning_results)

For the final results of all patient classification with 3 methods each in [diabetes_results](https://github.com/Feasao/RWRLab/tree/main/Random_walk/diabetes_results) and [hyperchol_results](https://github.com/Feasao/RWRLab/tree/main/Random_walk/hyperchol_results)

For a held out test set of hypercholesterolemia patients in [Test_resultsHyper](https://github.com/Feasao/RWRLab/tree/main/Random_walk/Test_resultsHyper)

## Key Insights

Network pruning significantly improves both separability and prediction quality
Sparse matrix operations are critical for scalability
Vectorization provides massive performance gains without changing algorithmic behavior

## Reproduction

Input the seeds txt file with the seeds of the patient in "Input_Files/SeedFolder"
Input the gene_network metabolite_network, diseases_network and metabolite_disease relation files in "Networks" folder
set the working directory to the correct one in the setwd() command
Run the file ([RStudio](https://posit.co/downloads/) recommended).

All functions for modular reproducibility for unweighted, weighted, topK and Elasticnet
SupraAdjacencyMatrix as well as the vectorised and original multiplex transition matrix creation can be found in [All_Functions.R](https://github.com/Feasao/RWRLab/blob/main/Random_walk/Functions/All_Functions.R)

