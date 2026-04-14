
## ----------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
##only change this path (setwd) and it should work. If you are having issues also see
##other places with paths and input your full path of the files
setwd("C:/Users/your_path/RWRLab/Random_walk")
Network_List <- c("GENES","METABOL")
#PPI=GENES
#PATH=METABOL
L <- length(Network_List)
L


## ----------------------------------------------------------------------------------------------------------------------------------
#Install packages
library(igraph)
library(Matrix)
library(biomaRt)
library(hexbin)
library(supraHex)
library(dnet)
library(Rcpp)


## ----------------------------------------------------------------------------------------------------------------------------------
source("Functions/All_Functions.R")
source("Functions/CreateNetworks_TopMultiplexHeterogeneous.R") 


## ----------------------------------------------------------------------------------------------------------------------------------
# Read input files
print("Reading arguments...")
args <- commandArgs(trailingOnly = TRUE)

SEED_FOLDER <- "Input_files/SeedFolder"
files <- list.files(SEED_FOLDER, full.names = TRUE)
for (file in files){
#All_Seeds <- read.csv("Input_files/Seeds_ExampleNew.txt",header=FALSE,sep="\t",dec=".",stringsAsFactors = FALSE)
All_Seeds <- read.csv(file,header=FALSE,sep="\t",dec=".",stringsAsFactors = FALSE)
All_Seeds <- All_Seeds$V1


## ----------------------------------------------------------------------------------------------------------------------------------
Parameters_File <- read.csv("Input_files/Parameters_Example.txt", header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)
Parameters <- check.parameters(Parameters_File,L )


## ----------------------------------------------------------------------------------------------------------------------------------
# Read the multiplex networks
size <- length(Network_List)
Layers <- vector("list", size)

for (i in 1:size) {
  if (Network_List[i] == "GENES") {
    
    GENES_table <- read.table("Networks/gene_network.tsv", sep = "\t", header = TRUE)
    GENES_Network <- graph_from_data_frame(GENES_table, directed = FALSE)
    GENES_Network <- igraph::simplify(GENES_Network, remove.multiple = TRUE, remove.loops = TRUE)
    
    Layers[[1]] <- GENES_Network
    names(Layers)[1] <- "GENES_NETWORK"
    
  } else {
    if (Network_List[i] == "METABOL") {
      
      METABOLOMICS_table <- read.delim("Networks/metabolite_network.tsv", header = TRUE)
      METABOLOMICS_Network <- graph_from_data_frame(METABOLOMICS_table, directed = FALSE)
      METABOLOMICS_Network <- igraph::simplify(METABOLOMICS_Network, remove.multiple = TRUE, remove.loops = TRUE)
      
      Layers[[2]] <- METABOLOMICS_Network
      names(Layers)[2] <- "METABOLOMICS_NETWORK"
      
    }
  }
}




## ----------------------------------------------------------------------------------------------------------------------------------
List_Layers <-Layers
List_Layers


## ----------------------------------------------------------------------------------------------------------------------------------
# Pool of nodes
pool_nodes <- pool.of.nodes(List_Layers)
pool_nodes_sorted <- sort(pool_nodes)

List_Layers_Allnodes <- add.missing.nodes(List_Layers, pool_nodes)

N <- Get.Number.Nodes(List_Layers_Allnodes)
N

lambda<- 0.01 #lambda and alpha values for the elasticnet implementation
alpha=1
TopK = 0.95 # Top K percentage of edges to keep in the multiplex network
## ----------------------------------------------------------------------------------------------------------------------------------
# Adjacency matrix of the multiplex network.
# If you want to use the unweighted version, you can change the function call to get.supra.adj.multiplex
print("Generating Multiplex Adjacency Matrix...")
SupraAdjacencyMatrix <- get.supra.adj.multiplexWeightedEN(List_Layers_Allnodes,Parameters$delta,N,lambda,alpha)


## ----------------------------------------------------------------------------------------------------------------------------------
# Heterogenous (Disease) network

print("Reading the disease-similarity network...")
Disease_table <- read.table("Networks/diseases_network.csv",sep=",")
Disease_Network <- graph_from_data_frame(Disease_table,directed=FALSE)
Disease_Network <- igraph::simplify(Disease_Network, remove.multiple = TRUE, remove.loops = TRUE)
AdjMatrix_Diseases <- as_adjacency_matrix(Disease_Network,sparse = TRUE)


## ----------------------------------------------------------------------------------------------------------------------------------
M <- nrow(AdjMatrix_Diseases)
M # Number of diseases


## ----------------------------------------------------------------------------------------------------------------------------------
print("Checking Seed nodes...")
Seed_File_list <- check.seeds(All_Seeds,pool_nodes,rownames(AdjMatrix_Diseases))


## ----------------------------------------------------------------------------------------------------------------------------------
# Generation of bipartite graph

print("Reading Ensemble file with Gene-disease relations from OMIM...")
Gene_Phenotype_relation <- get.disease.gene.relations(pool_nodes)

Metabolite_Phenotype_relation <- read.delim("Input_files/metabolite_disease.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(Metabolite_Phenotype_relation) <- colnames(Gene_Phenotype_relation)
Gene_Phenotype_relation <- rbind(Gene_Phenotype_relation, Metabolite_Phenotype_relation)


## ----------------------------------------------------------------------------------------------------------------------------------
print("Generating Bipartite Graph...")

Bipartite_matrix_and_report <- get.bipartite.graph(pool_nodes_sorted, colnames(AdjMatrix_Diseases), Gene_Phenotype_relation,N,M)
Bipartite_matrix <- Bipartite_matrix_and_report[[1]]
Error_log <- Bipartite_matrix_and_report[[2]]


## ----------------------------------------------------------------------------------------------------------------------------------
print("Adapting Bipartite Graph to the Multiplex...")
SupraBipartiteMatrix <- expand.bipartite.graph(N,L,M,Bipartite_matrix)


## ----------------------------------------------------------------------------------------------------------------------------------
# Transition matrices
Parameters$lambda=0.5
# inter-subnetworks links
Transition_Protein_Disease <- get.transition.protein.disease(N,L,M,SupraBipartiteMatrix,Parameters$lambda)
Transition_Disease_Protein <- get.transition.disease.protein(N,L,M,SupraBipartiteMatrix,Parameters$lambda)


## ----------------------------------------------------------------------------------------------------------------------------------
# Intra-subnet links
Transition_Multiplex_Network <- get.transition.multiplexFAST(N,L,Parameters$lambda,SupraAdjacencyMatrix,SupraBipartiteMatrix)
Transition_Disease_Network <- get.transition.disease(M,Parameters$lambda,AdjMatrix_Diseases,SupraBipartiteMatrix)


## ----------------------------------------------------------------------------------------------------------------------------------
# Global Transition matrices
Transition_Multiplex_Heterogeneous_Matrix_1 <- cbind(Transition_Multiplex_Network, Transition_Protein_Disease)
Transition_Multiplex_Heterogeneous_Matrix_2 <- cbind(Transition_Disease_Protein, Transition_Disease_Network)
Transition_Multiplex_Heterogeneous_Matrix <- rbind(Transition_Multiplex_Heterogeneous_Matrix_1,Transition_Multiplex_Heterogeneous_Matrix_2)


## ----------------------------------------------------------------------------------------------------------------------------------
# Scores

Gene_Seeds <- Seed_File_list[[1]]
Disease_Seeds <- Seed_File_list[[2]]

tau <- Parameters$tau/L


## ----------------------------------------------------------------------------------------------------------------------------------
# We compute the restart probability of each seed, based on eta and tau.
eta=0.5
Seeds_Score <- get.seed.scores(Gene_Seeds,Disease_Seeds,eta,L,tau)
print(Seeds_Score)


## ----------------------------------------------------------------------------------------------------------------------------------
print("Performing Random Walk...")

Random_Walk_Results <- Random_Walk_Restart(Transition_Multiplex_Heterogeneous_Matrix,Parameters$r,Seeds_Score)
sourceCpp("Functions/Geometric_Mean.cpp")

## ----------------------------------------------------------------------------------------------------------------------------------
Geometric_Mean <- function(Scores, L, N) {
  FinalScores <- rep(1.0, N)  # Initialize a vector of ones
  
  for (i in 1:N) {
    for (j in 1:L) {
      FinalScores[i] <- FinalScores[i] * Scores[i + (j - 1) * N]
    }
    FinalScores[i] <- FinalScores[i]^(1.0 / L)  # Compute geometric mean
  }
  
  return(FinalScores)
}

diseases_map <- c(
  "125853" = "TYPE 2 DIABETES MELLITUS",
  "125850" = "MATURITY-ONSET DIABETES OF THE YOUNG",
  "222100" = "TYPE 1 DIABETES MELLITUS",
  "603813" = "HYPERCHOLESTEROLEMIA, FAMILIAL, 4",
  "603776" = "HYPERCHOLESTEROLEMIA, FAMILIAL, 3",
  "144010" = "HYPERCHOLESTEROLEMIA, FAMILIAL, 2",
  "143890" = "HYPERCHOLESTEROLEMIA, FAMILIAL, 1",
  
  "610370" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 1",
  "125851" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 2",
  "600496" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 3",
  "606392" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 4",
  "606394" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 6",
  "610508" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 7",
  "609812" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 8",
  "612225" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 9",
  "613370" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 10",
  "613375" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 11",
  "614521" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 13",
  "616511" = "MATURITY-ONSET DIABETES OF THE YOUNG, TYPE 14",
  
  "601410" = "TYPE 1 DIABETES MELLITUS 2",
  "612522" = "TYPE 1 DIABETES MELLITUS 6",
  "614456" = "TYPE 1 DIABETES MELLITUS 10",
  "616461" = "TYPE 1 DIABETES MELLITUS 15",
  "617219" = "TYPE 1 DIABETES MELLITUS 20",
  
  "601410" = "DIABETES MELLITUS, TRANSIENT NEONATAL, 1",
  "618856" = "DIABETES MELLITUS, TRANSIENT NEONATAL, 2",
  "618857" = "DIABETES MELLITUS, TRANSIENT NEONATAL, 3",
  
  "304800" = "DIABETES INSIPIDUS, NEPHROGENIC, 1, X-LINKED",
  "125800" = "DIABETES INSIPIDUS, NEUROHYPOPHYSEAL, X-LINKED",
  "222000" = "DIABETES INSIPIDUS, NEPHROGENIC, 2, AUTOSOMAL"
)

## ----------------------------------------------------------------------------------------------------------------------------------
final_rank_proteins <- rank_proteins(N, L,Random_Walk_Results,Gene_Seeds)
final_rank_diseases <- rank_diseases(N,L,M,Random_Walk_Results,Disease_Seeds)
final_rank_diseases$DiseaseName <-diseases_map[as.character(final_rank_diseases$DiseaseID)]

out_name <- paste0( "result_", basename(file))
out_path <- file.path("Test_results/", out_name)
## ----------------------------------------------------------------------------------------------------------------------------------
write.table(final_rank_proteins,file="final_rank_proteins.txt",sep="\t",row.names = FALSE, dec=".",quote=FALSE)
write.table(final_rank_diseases,file=out_path,sep="\t",row.names = FALSE, dec=".",quote=FALSE)


## ----------------------------------------------------------------------------------------------------------------------------------
print("Creating Network file with the top candidates...")
Top_Results_Network <- CreateNetworks_TopMultiplexHeterogeneous(Network_List,c(Gene_Seeds,Disease_Seeds),
                                                                final_rank_proteins$GeneNames[1:Parameters$k],final_rank_diseases$DiseaseID[1:Parameters$k])


## ----------------------------------------------------------------------------------------------------------------------------------
write.table(as_data_frame(Top_Results_Network, what = c("edges", "vertices", "both")), file = "Final_network.txt", 
            sep=" ", quote=FALSE,row.names=FALSE)
}
