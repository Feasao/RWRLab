### Functions called from RWR-M.R

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# GENERAL FUNCTIONS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 1.- Package installation and load of libraries.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

install.packages.if.necessary <- function(CRAN.packages=c(), bioconductor.packages=c()) {
  if (length(bioconductor.packages) > 0) {
    source("http://bioconductor.org/biocLite.R")
  }
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      biocLite(p) 
      library(p, character.only=T)
    }
  }
  
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p) 	
      library(p, character.only=T)  	
    }	
  }
  
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 2.- We check the input parameters
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
check.parameters <- function(input_file,Nr_Networks){
  
  if (ncol(input_file) != 2) { stop("Incorrect Parameters Input File")}
  if (nrow(input_file) != 4) { stop("Incorrect Parameters Input File")}
  
  r <- as.numeric(input_file[1,2] )
  print(paste("Global restart probability:", r))
  if ((r > 1 || r <= 0)){ stop("Incorrect r, it must be between 0 and 1")}
  
  delta <- as.numeric(input_file[2,2] )
  print(paste("Inter-Layer Jump Probability:", delta))
  if ((delta > 1 || delta < 0)){ stop("Incorrect delta, it must be between 0 and 1")}
  
  tau <- as.numeric(unlist(strsplit(input_file[3,2],",")))
  print(paste("Layers Restart Probability: ", input_file[3,2]))
  if (sum(tau)/L != 1) {stop("Incorrect tau, the sum of its component divided by 3 must be 1")}
  
  k <- as.numeric(input_file[4,2] )
  print(paste("Number of Top results to be displayed in the output network:", k))
  if ((k > 200 || k <= 0)){ stop("Incorrect k, it must be between 0 and 200")}
  
  parameters <- list(r,delta,tau,k)
  names(parameters) <- c("r", "delta", "tau", "k")
  return(parameters)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 3.- We check if the seed nodes are in our network.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
check.seeds <- function(Seeds, All_proteins){
  
  Genes_Seeds_OK <- Seeds[which(Seeds %in% All_proteins)]
  Genes_Seeds_KO <- Seeds[which(!Seeds %in% All_proteins)]
  
  if (length(Genes_Seeds_OK) == 0) {
    stop("Seeds not found in our network")
  } else {
    if (length(Genes_Seeds_KO) > 0){
      print("Some Seed genes no present in the network: ")
      for (i in 1:length(Genes_Seeds_KO)){
        print(Genes_Seeds_KO[i])
      }
    }
    return(Genes_Seeds_OK)
  }
  
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# MULTIPLEX RELATED FUNCTIONS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 4.- We read the different layers that integrate the multiplex network.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

read.layers <- function(vec_layers){
  
  ## "pre-allocate" an empty list of length size (Number of Layers)
  size <- length(vec_layers)
  Layers <- vector("list", size)
  
  ## We read the different Networks (Layers) of our Multiplex network. We also simplify the networks
  ## by removing possible self loops and possible multiple nodes.
  
  for (i in 1:size){
    if (vec_layers[i]=="GENES"){
      
      GENES_table <- read.table("Networks/gene_network.tsv",sep="\t")
      GENES_Network <- graph_from_data_frame(GENES_table,directed=FALSE)
      GENES_Network <- igraph::simplify(GENES_Network, remove.multiple = TRUE, remove.loops = TRUE)
      
      Layers[[i]] <- GENES_Network
      names(Layers)[i] <- "GENES_NETWORK"
    } else {
      if (vec_layers[i]=="METABOL"){
        
        METABOLway_table <- read.delim("Networks/metabolomics_network.tsv")
        METABOLway_Network <- graph_from_data_frame(METABOLway_table,directed=FALSE)
        METABOLway_Network <- igraph::simplify(METABOLway_Network, remove.multiple = TRUE, remove.loops = TRUE)
        
        Layers[[i]] <- METABOLway_Network
        names(Layers)[i] <- "METABOLWAY_NETWORK"
        
      } 
    }
  }
  ## We return a list containing an igraph object for each layer.
  return(Layers)
} 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 5.- We generate a pool of nodes. We merge the nodes present in every layer.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
pool.of.nodes <- function(Layers){
  
  ## We get the number of layers
  Nr_Layers <- length(Layers)
  
  ## We get the nodes of all the layers of the multiplex network. We save them into a vector.
  Node_Names_all <- character()
  for (i in 1:Nr_Layers) {
    Node_Names_Layer <- V(Layers[[i]])$name
    Node_Names_all <-c(Node_Names_all,Node_Names_Layer)
  }
  
  ## We remove duplicates.
  Node_Names_all <- unique(Node_Names_all)
  
  return(Node_Names_all)
} 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 6.- From the pool of nodes we add the missing proteins to each layer as 
#     isolated nodes.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
add.missing.nodes <- function (Layers,NodeNames) {
  
  ## We get the number of layers
  Nr_Layers <- length(Layers)
  
  ## We generate a new list of layers.
  Layers_New <- vector("list", Nr_Layers)
  
  ## We add to each layer the missing nodes of the total set of nodes, of the pool of nodes.
  for (i in 1:Nr_Layers){
    Node_Names_Layer <- V(Layers[[i]])$name
    Missing_Nodes <- NodeNames[which(!NodeNames %in% Node_Names_Layer)]
    Layers_New[[i]] <- add_vertices(Layers[[i]] ,length(Missing_Nodes), name=Missing_Nodes)
  }
  return(Layers_New)
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 7.- We check the total number of nodes in every layer and we return it. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
Get.Number.Nodes <- function(Layers_Allnodes) {
  
  ## We get the number of layers
  Nr_Layers <- length(Layers_Allnodes)
  vector_check <- numeric(length = Nr_Layers)  
  
  for (i in 1:Nr_Layers){
    vector_check[i] <- vcount(Layers_Allnodes[[i]])  
  }
  
  if (all(vector_check == vector_check[1])){
    print("Number of nodes in every layer updated...")
    return(vector_check[1])  
  } else {
    stop("Not correct number of nodes in each Layer...")
  }
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 8.- We build the supra adjacency matrix (not normalised yet) without edge weights. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

get.supra.adj.multiplex <- function(Layers,delta,N){
  
  ## IDEM_MATRIX.
  Idem_Matrix <- Diagonal(N, x = 1)
  L <- length(Layers)
  
  SupraAdjacencyMatrix <- Matrix(0,ncol=N*L,nrow=N*L,sparse = TRUE)
  
  Col_Node_Names <- character()
  Row_Node_Names <- character()
  
  for (i in 1:L){
    Adjacency_Layer <-  as_adjacency_matrix(Layers[[i]],sparse = TRUE)
    
    ## We order the matrix by the node name. This way all the matrix will have the same. Additionally we include a label with the layer number for each node name.
    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),order(colnames(Adjacency_Layer))]
    Layer_Col_Names <- paste(colnames(Adjacency_Layer),i,sep="_")
    Layer_Row_Names <- paste(rownames(Adjacency_Layer),i,sep="_")
    Col_Node_Names <- c(Col_Node_Names,Layer_Col_Names)
    Row_Node_Names <- c(Row_Node_Names,Layer_Row_Names)
    
    ## We fill the diagonal blocks with the adjacencies matrix of each layer.
    Position_ini_row <- 1 + (i-1)*N
    Position_end_row <- N + (i-1)*N
    SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_row:Position_end_row)] <- (1-delta)*(Adjacency_Layer)
    
    ## We fill the off-diagonal blocks with the transition probability among layers.
    for (j in 1:L){
      Position_ini_col <- 1 + (j-1)*N
      Position_end_col <- N + (j-1)*N
      if (j != i){
        SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (delta/(L-1))*Idem_Matrix
      }
    }
  }
  
  rownames(SupraAdjacencyMatrix) <- Row_Node_Names
  colnames(SupraAdjacencyMatrix) <- Col_Node_Names
  
  return(SupraAdjacencyMatrix)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 8.0- We build the supra adjacency matrix (not normalised yet) without edge weights. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

get.supra.adj.multiplexOneLayer <- function(Layers,delta,N){
  
  ## IDEM_MATRIX.
  Idem_Matrix <- Diagonal(N, x = 1)
  L <- length(Layers)
  
  SupraAdjacencyMatrix <- Matrix(0,ncol=N*L,nrow=N*L,sparse = TRUE)
  
  Col_Node_Names <- character()
  Row_Node_Names <- character()
  
  for (i in 1:L){
    Adjacency_Layer <-  as_adjacency_matrix(Layers[[1]],sparse = TRUE)
    
    ## We order the matrix by the node name. This way all the matrix will have the same. Additionally we include a label with the layer number for each node name.
    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),order(colnames(Adjacency_Layer))]
    Layer_Col_Names <- paste(colnames(Adjacency_Layer),i,sep="_")
    Layer_Row_Names <- paste(rownames(Adjacency_Layer),i,sep="_")
    Col_Node_Names <- c(Col_Node_Names,Layer_Col_Names)
    Row_Node_Names <- c(Row_Node_Names,Layer_Row_Names)
    
    ## We fill the diagonal blocks with the adjacencies matrix of each layer.
    Position_ini_row <- 1 + (i-1)*N
    Position_end_row <- N + (i-1)*N
    SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_row:Position_end_row)] <- (1-delta)*(Adjacency_Layer)
    
    ## We fill the off-diagonal blocks with the transition probability among layers.
    for (j in 1:L){
      Position_ini_col <- 1 + (j-1)*N
      Position_end_col <- N + (j-1)*N
      if (j != i){
        SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (delta/(L-1))*Idem_Matrix
      }
    }
  }
  
  rownames(SupraAdjacencyMatrix) <- Row_Node_Names
  colnames(SupraAdjacencyMatrix) <- Col_Node_Names
  
  return(SupraAdjacencyMatrix)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 8.5- We build the supra adjacency matrix (not normalised yet) with edge weights.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

get.supra.adj.multiplexWeighted <- function(Layers,delta,N){
  
  ## IDEM_MATRIX.
  Idem_Matrix <- Diagonal(N, x = 1)
  L <- length(Layers)
  SupraAdjacencyMatrix <- Matrix(0,ncol=N*L,nrow=N*L,sparse = TRUE)
  
  Col_Node_Names <- character()
  Row_Node_Names <- character()
  
  for (i in 1:L){
    Adjacency_Layer <-  as_adjacency_matrix(Layers[[i]], attr = "weight",sparse = TRUE)
    
    ## We order the matrix by the node name. This way all the matrix will have the same. Additionally we include a label with the layer number for each node name.
    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),order(colnames(Adjacency_Layer))]
    Layer_Col_Names <- paste(colnames(Adjacency_Layer),i,sep="_")
    Layer_Row_Names <- paste(rownames(Adjacency_Layer),i,sep="_")
    Col_Node_Names <- c(Col_Node_Names,Layer_Col_Names)
    Row_Node_Names <- c(Row_Node_Names,Layer_Row_Names)
    
    ## We fill the diagonal blocks with the adjacencies matrix of each layer.
    Position_ini_row <- 1 + (i-1)*N
    Position_end_row <- N + (i-1)*N
    SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_row:Position_end_row)] <- (1-delta)*(Adjacency_Layer)
    
    ## We fill the off-diagonal blocks with the transition probability among layers.
    for (j in 1:L){
      Position_ini_col <- 1 + (j-1)*N
      Position_end_col <- N + (j-1)*N
      if (j != i){
        SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (delta/(L-1))*Idem_Matrix
      }
    }
  }
  
  rownames(SupraAdjacencyMatrix) <- Row_Node_Names
  colnames(SupraAdjacencyMatrix) <- Col_Node_Names
  
  return(SupraAdjacencyMatrix)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 8.6- We build the supra adjacency matrix (not normalised yet) with edge weights and top K neighbor filtering.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

get.supra.adj.multiplexWeightedTopK <- function(Layers,delta,N,K){
  
  ## IDEM_MATRIX.
  Idem_Matrix <- Diagonal(N, x = 1)
  L <- length(Layers)

  SupraAdjacencyMatrix <- Matrix(0,ncol=N*L,nrow=N*L,sparse = TRUE)
  
  Col_Node_Names <- character()
  Row_Node_Names <- character()
  
  for (i in 1:L){
    Adjacency_Layer <-  as_adjacency_matrix(Layers[[i]], attr = "weight",sparse = TRUE)
    
    ## We order the matrix by the node name. This way all the matrix will have the same. Additionally we include a label with the layer number for each node name.
    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),order(colnames(Adjacency_Layer))]
    # --- Top-K filtering to keep only K relevant nodes---
    # Because the top K neighbors are not always the same in each layer, we filter the top K neighbors for each layer.

    # print("Filtering top K neighbors for layer:", i)
    # print("total edges before")
    # print(sum(Adjacency_Layer != 0))

    #found that running TopK in the gene layer is better than all the layers because the others are really sparce hence the i == 1
    if (!is.null(K) && K > 0 && i == 1){ 
      degree_vector <- rowSums(Adjacency_Layer != 0)
      # This is to remove the zero degree nodes that skew the average degree. You can tweak it to remove outliers.
      degree_nonzero <- degree_vector[degree_vector > 0]
      #calculating the avg degree of non-zero nodes 
      avg_nonzero <- mean(degree_nonzero) 
      TopK = ceiling(avg_nonzero*K)
      # needs to be an dense matrix for the filtering function
      Adjacency_Layer <- filter_topK_neighbors(as.matrix(Adjacency_Layer),TopK)
      # convert back to sparse matrix
      Adjacency_Layer <- Matrix(Adjacency_Layer, sparse = TRUE)

    }

    # print("total edges after")
    # print(sum(Adjacency_Layer != 0))

    # -----------------------
    Layer_Col_Names <- paste(colnames(Adjacency_Layer),i,sep="_")
    Layer_Row_Names <- paste(rownames(Adjacency_Layer),i,sep="_")
    Col_Node_Names <- c(Col_Node_Names,Layer_Col_Names)
    Row_Node_Names <- c(Row_Node_Names,Layer_Row_Names)
    
    ## We fill the diagonal blocks with the adjacencies matrix of each layer.
    Position_ini_row <- 1 + (i-1)*N
    Position_end_row <- N + (i-1)*N
    SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_row:Position_end_row)] <- (1-delta)*(Adjacency_Layer)
    
    ## We fill the off-diagonal blocks with the transition probability among layers.
    for (j in 1:L){
      Position_ini_col <- 1 + (j-1)*N
      Position_end_col <- N + (j-1)*N
      if (j != i){
        SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (delta/(L-1))*Idem_Matrix
      }
    }
  }
  
  rownames(SupraAdjacencyMatrix) <- Row_Node_Names
  colnames(SupraAdjacencyMatrix) <- Col_Node_Names
  
  return(SupraAdjacencyMatrix)
}
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 8.7- Top K neighbor filtering for the supra adjacency matrix. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
filter_topK_neighbors <- function(adj_matrix, topK) {
  n <- nrow(adj_matrix)
  #   We create a boolean matrix to store the top K neighbors
  topK_rows <- matrix(FALSE, nrow = n, ncol = n)
  topK_cols <- matrix(FALSE, nrow = n, ncol = n)

  if (!isSymmetric(adj_matrix)) {
    stop("This topK assumes undirected graph and the adjacency matrix must be symmetric.")
  }
  for (i in 1:n) {
    # For each row, we find the indices of the non-zero elements
    row_vals <- adj_matrix[i, ]
    nz_idx <- which(row_vals != 0)
    # If there are more than K non-zero elements, we keep only the top K
    # Else we keep all of them
    if (length(nz_idx) > topK) {
      top_idx <- nz_idx[order(row_vals[nz_idx], decreasing = TRUE)][1:topK]
    } else {
      top_idx <- nz_idx
    }
    # We set the top idx in the boolean matrix to true
    topK_rows[i, top_idx] <- TRUE
  }

  for (j in 1:n) {
    col_vals <- adj_matrix[, j]
    nz_idx <- which(col_vals != 0)
    # Same logic as above for columns because the adjacency matrix is symmetric
    if (length(nz_idx) > topK) {
      top_idx <- nz_idx[order(col_vals[nz_idx], decreasing = TRUE)][1:topK]
    } else {
      top_idx <- nz_idx
    }
    topK_cols[top_idx, j] <- TRUE
  }
  # We multiply the original adjacency matrix with the boolean matrices to keep only the top K neighbors
  # This will zero out all the non-top K neighbors
  result <- adj_matrix * (topK_rows & topK_cols)
  return(result)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 8.8- Elasticnet filtering for the supra adjacency matrix. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
get.supra.adj.multiplexWeightedEN <- function(Layers, delta, N, lambda, alpha) {
  
  Idem_Matrix <- Diagonal(N, x = 1)
  L <- length(Layers)
  SupraAdjacencyMatrix <- Matrix(0, ncol=N*L, nrow=N*L, sparse=TRUE)
  
  Col_Node_Names <- character()
  Row_Node_Names <- character()

  # adj_genes <- as_adjacency_matrix(Layers[[1]], attr="weight", sparse=TRUE)
  # summary(adj_genes@x)  # distribution of non-zero edge weights
  # adj_genes <- as_adjacency_matrix(Layers[[1]], attr="weight", sparse=TRUE)
  # cat("Edges > 0.1:", sum(adj_genes@x > 0.1), "\n")
  # cat("Edges > 1.0:", sum(adj_genes@x > 1.0), "\n")
  # cat("Edges > 0.005:", sum(adj_genes@x > 0.005), "\n")
  # so here the implementation is elasticnet but for this dataset ridge is redundant so we set a = 1
  for (i in 1:L) {
    Adjacency_Layer <- as_adjacency_matrix(Layers[[i]], attr="weight", sparse=TRUE)
    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),order(colnames(Adjacency_Layer))]
    
    if (i == 1) {
      lambda1 <- lambda * alpha
      lambda2 <- lambda * (1 - alpha) / 2

      # shrinks weak edges to zero and scales strong edges by L2 which for a = 1 it does not scale them
      Adjacency_Layer <- sign(Adjacency_Layer) * pmax(abs(Adjacency_Layer) - lambda1, 0) / (1 + 2 * lambda2)
      Adjacency_Layer <- drop0(Adjacency_Layer)
  
    }
    
    Layer_Col_Names <- paste(colnames(Adjacency_Layer), i, sep="_")
    Layer_Row_Names <- paste(rownames(Adjacency_Layer), i, sep="_")
    Col_Node_Names  <- c(Col_Node_Names, Layer_Col_Names)
    Row_Node_Names  <- c(Row_Node_Names, Layer_Row_Names)
    
    Position_ini_row <- 1 + (i-1)*N
    Position_end_row <- N + (i-1)*N
    SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_row:Position_end_row)] <- (1-delta)*Adjacency_Layer
    
    for (j in 1:L) {
      Position_ini_col <- 1 + (j-1)*N
      Position_end_col <- N + (j-1)*N
      if (j != i) {
        SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (delta/(L-1))*Idem_Matrix
      }
    }
  }
  
  rownames(SupraAdjacencyMatrix) <- Row_Node_Names
  colnames(SupraAdjacencyMatrix) <- Col_Node_Names
  
  return(SupraAdjacencyMatrix)
}
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 9.- GET THE SCORES OF THE DIFFERENT SEEDS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

get.seed.scores <- function(Genes,Diseases,eta,Number_Layers,tau) {
  
  
  if ((length(Genes) != 0 && length(Diseases)!= 0)){
    Seed_Genes_Layer_Labeled <- character(length = length(Genes)*Number_Layers)
    Seeds_Genes_Scores <- numeric(length = length(Genes)*Number_Layers)
    
    Current_Gene_Labeled <- character()
    
    for (k in 1:L){
      Current_Gene_Labeled <- c(Current_Gene_Labeled,paste(Genes[k],k,sep="_",collapse = "") )
      for (j in 1:length(Genes)){
        Seed_Genes_Layer_Labeled[((k-1)*length(Genes))+ j] <- paste(Genes[j],k,sep="_",collapse = "") 
        Seeds_Genes_Scores[((k-1)*length(Genes))+ j] <- ((1-eta) * tau[k])/length(Genes)
      }
    }  
    
    Disease_Seeds_Score <- eta/length(Diseases) 
    
  } else {
    eta <- 1
    if (length(Genes) == 0){
      Seed_Genes_Layer_Labeled <- character()
      Seeds_Genes_Scores <- numeric()
      Disease_Seeds_Score <- eta/length(Diseases) 
    } else {
      
      Seed_Genes_Layer_Labeled <- character(length = length(Genes)*Number_Layers)
      Seeds_Genes_Scores <- numeric(length = length(Genes)*Number_Layers)
      
      Current_Gene_Labeled <- character()
      
      for (k in 1:L){
        Current_Gene_Labeled <- c(Current_Gene_Labeled,paste(Genes[k],k,sep="_",collapse = "") )
        for (j in 1:length(Genes)){
          Seed_Genes_Layer_Labeled[((k-1)*length(Genes))+ j] <- paste(Genes[j],k,sep="_",collapse = "") 
          Seeds_Genes_Scores[((k-1)*length(Genes))+ j] <- tau[k]/length(Genes)
          Disease_Seeds_Score <- numeric()
        }
      } 
    }
  }
 ### We prepare a data frame with the seeds.
  Seeds_Score <- data.frame(Seeds_ID = c(Seed_Genes_Layer_Labeled,Diseases),
                            Score = c(Seeds_Genes_Scores,rep(Disease_Seeds_Score,length(Diseases)))  ,stringsAsFactors = FALSE)
  return(Seeds_Score)
}
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 10.- RANDOM WALK WITH RESTART
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

Random_Walk_Restart <- function(Network_Matrix, r,SeedGenes ){
  
  ### We define the threshold and the number maximum of iterations for the randon walker.
  Threeshold <- 1e-10
  NetworkSize <- ncol(Network_Matrix)
  
  ### We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1
  
  #### We define the prox_vector(The vector we will move after the first RW iteration. We start from The seed. We have to take in account
  #### that the walker with restart in some of the Seed genes, depending on the score we gave in that file).
  prox_vector <- matrix(0,nrow = NetworkSize,ncol=1)
  
  prox_vector[which(colnames(Network_Matrix) %in% SeedGenes[,1])] <- (SeedGenes[,2])
  
  prox_vector  <- prox_vector/sum(prox_vector)
  restart_vector <-  prox_vector
  
  while(residue >= Threeshold){
    
    old_prox_vector <- prox_vector
    prox_vector <- (1-r)*(Network_Matrix %*% prox_vector) + r*restart_vector
    residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
    iter <- iter + 1; 
  }
  return(prox_vector) 
} 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 14.- RANKING OF PROTEINS AFTER RANDOM WALK. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
rank_proteins <- function(Number_Proteins, Number_Layers,Results,Seeds){
  
  ## We sort the score to obtain the ranking of Proteins and Diseases.
  proteins_rank <- data.frame(GeneNames = character(length = Number_Proteins), Score = 0)
  proteins_rank$GeneNames <- gsub("_1", "", row.names(Results)[1:Number_Proteins])
  
  ## We calculate the Geometric Mean among the proteins in the different layers.
  proteins_rank$Score <- Geometric_Mean(as.vector(Results[,1]),Number_Layers,Number_Proteins)
  
  proteins_rank_sort <- proteins_rank[with(proteins_rank, order(-Score, GeneNames)), ]
  
  ### We remove the seed genes from the Ranking
  proteins_rank_sort_NoSeeds <- proteins_rank_sort[which(!proteins_rank_sort$GeneNames %in% Seeds),]
  
  return(proteins_rank_sort_NoSeeds)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 15.- RANKING OF DISEASES AFTER RANDOM WALK. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
rank_diseases <- function(Number_Proteins,Number_Layers,Number_Diseases,Results,Seeds){
  
  ## rank_diseases
  diseases_rank <- data.frame(DiseaseID = character(length = Number_Diseases), Score = 0)
  diseases_rank$DiseaseID <- row.names(Results)[((Number_Proteins*Number_Layers)+1):nrow(Results)]
  diseases_rank$Score <- Results[((Number_Proteins*Number_Layers)+1):nrow(Results),1]
  
  diseases_rank_sort <- diseases_rank[with(diseases_rank, order(-Score, DiseaseID)), ]
  diseases_rank_sort_NoSeeds <- diseases_rank_sort[which(!diseases_rank_sort$DiseaseID %in% Seeds),]
  
  return(diseases_rank_sort_NoSeeds)
}

check.seeds <- function(Seeds, All_proteins,All_Diseases){
  
  Genes_Seeds_Ok <- Seeds[which(Seeds %in% All_proteins)]
  Disease_Seeds_Ok <- Seeds[which(Seeds %in% All_Diseases)]
  All_seeds_ok <- c(Genes_Seeds_Ok,Disease_Seeds_Ok)
  All_seeds_ko <- Seeds[which(!Seeds %in% All_seeds_ok)]

  list_Seeds_Ok <- list(Genes_Seeds_Ok,Disease_Seeds_Ok)
  
  print("Seeds OK: ")
  print(paste(All_seeds_ok, sep=" "))
  print("Seeds KO: ")
  print(paste(All_seeds_ko, sep=" "))
  
  if ((length(Genes_Seeds_Ok) == 0) &&  (length(Disease_Seeds_Ok) ==0)){
    stop("Seeds not found in our network")
  } else {
    return(list_Seeds_Ok)
  }
  
}

get.disease.gene.relations <- function(proteins){
  
  Gene_Phenotype_relation <- read.delim("Input_files/Gene_Phenotype_relation.txt",
                                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
  
  ## We remove Gene-phenotype relations that are not in our proteins file. 
  #Gene_Phenotype_relation <- Gene_Phenotype_relation[which(Gene_Phenotype_relation$HMDB %in% OMIM) , ]
  
  return(Gene_Phenotype_relation)
}

get.bipartite.graph <- function(proteins_sorted, disease_names, Gene_Phenotype_relation,Number_Proteins,Number_Diseases){
  Bipartite_matrix <- Matrix(data=0, nrow=Number_Proteins, ncol=Number_Diseases)
  rownames(Bipartite_matrix) <- proteins_sorted
  colnames(Bipartite_matrix) <- disease_names
  log <- character()
  
  for (i in 1:Number_Proteins){
    current_gene <- proteins_sorted[i]
    current_mim <- Gene_Phenotype_relation$mim_morbid[which(Gene_Phenotype_relation$hgnc_symbol == current_gene)]
    
    # We have to check if the gene is the Gene_Phenotype_relation downloaded from Biomart. (Some weird HGNC symbols, not retrieved.)
    if (length(current_mim) > 0){
      for (j in 1:length(current_mim)){
        if (!is.na(current_mim[j])){
          # We need to identify the phenotypes position on the matrix.
          index_disease <- which(colnames(Bipartite_matrix) %in%  current_mim[j])
          # We have to check if that index is present in the matrix.
          if (length(index_disease) == 1){ 
            Bipartite_matrix[i,index_disease] <- 1
          } else {
            error_message <- paste("MIM_CODE", current_mim[j], length(index_disease), "No phenotype found",sep=";", collapse = NULL)
            log <- c(log,error_message)
          }           
        }
      }  
    } else {
      error_message <- paste("HGNC_Symbol", current_gene, "No HGNC found in Biomart",sep=";", collapse = NULL)
      log <- c(log,error_message)
    }
  }
  Bipartite_and_errorlog <- list(Bipartite_matrix,log)
  return(Bipartite_and_errorlog)
}

expand.bipartite.graph <- function(Number_Proteins,Number_Layers,Number_Diseases,Bipartite_matrix){
  
  SupraBipartiteMatrix <- Matrix(0,nrow=Number_Proteins*Number_Layers,ncol=Number_Diseases,sparse = TRUE)
  Row_Node_Names <- character()
  
  for (i in 1:Number_Layers){
    Layer_Row_Names <- paste(rownames(Bipartite_matrix),i,sep="_")
    Row_Node_Names <- c(Row_Node_Names,Layer_Row_Names)
    Position_ini_row <- 1 + (i-1)*Number_Proteins
    Position_end_row <- Number_Proteins + (i-1)*Number_Proteins
    SupraBipartiteMatrix[(Position_ini_row:Position_end_row),] <- Bipartite_matrix
  }  
  
  rownames(SupraBipartiteMatrix) <- Row_Node_Names
  colnames(SupraBipartiteMatrix) <- colnames(Bipartite_matrix)
  return(SupraBipartiteMatrix)
}

get.transition.protein.disease <- function(Number_Proteins,Number_Layers,Number_Diseases,SupraBipartiteMatrix,lambda){
  
  Transition_Protein_Disease <- Matrix(0,nrow=Number_Proteins*Number_Layers,ncol=Number_Diseases,sparse = TRUE)
  colnames(Transition_Protein_Disease) <- colnames(SupraBipartiteMatrix)
  rownames(Transition_Protein_Disease) <- rownames(SupraBipartiteMatrix)
  
  Col_Sum_Bipartite <- colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
  
  for (j in 1:Number_Diseases){
    if (Col_Sum_Bipartite[j] != 0){
      Transition_Protein_Disease[,j] <- (lambda*SupraBipartiteMatrix[,j]) /Col_Sum_Bipartite[j]
    }
  }
  return(Transition_Protein_Disease)
}

# 11.2.-Disease-Protein Transition Matrix.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
get.transition.disease.protein <- function(Number_Proteins,Number_Layers,Number_Diseases,SupraBipartiteMatrix,lambda){
  
  Transition_Disease_Protein <- Matrix(0,nrow=Number_Diseases,ncol=Number_Proteins*Number_Layers,sparse = TRUE)
  
  colnames(Transition_Disease_Protein) <- rownames(SupraBipartiteMatrix)
  rownames(Transition_Disease_Protein) <- colnames(SupraBipartiteMatrix)
  
  Row_Sum_Bipartite <- rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
  
  for (i in 1:(Number_Proteins*Number_Layers)){
    if (Row_Sum_Bipartite[i] != 0){
      Transition_Disease_Protein[,i] <- (lambda*SupraBipartiteMatrix[i,])/Row_Sum_Bipartite[i]
    }
  }
  return(Transition_Disease_Protein)
}

# 11.3.-Multiplex intra-transition Matrix. It's very slow... Transform to c++
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
get.transition.multiplex <- function(Number_Proteins,Number_Layers,lambda,SupraAdjacencyMatrix,SupraBipartiteMatrix){
  
  Transition_Multiplex_Network <- Matrix(0,nrow=Number_Proteins*Number_Layers,ncol=Number_Proteins*Number_Layers,sparse = TRUE)
  
  rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
  colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)
  
  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrix,na.rm=FALSE,dims=1, sparseResult=FALSE)
  Row_Sum_Bipartite <- rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE) 
  
  for (j in 1:(Number_Proteins*Number_Layers)){
    if(Row_Sum_Bipartite[j] != 0){
      Transition_Multiplex_Network[,j] <- ((1-lambda)*SupraAdjacencyMatrix[,j]) /Col_Sum_Multiplex[j]
    } else {
      Transition_Multiplex_Network[,j] <- SupraAdjacencyMatrix[,j] /Col_Sum_Multiplex[j]
    }
  }
  return(Transition_Multiplex_Network)
}

# 11.3.5-Multiplex intra-transition Matrix but using matrices to make it really quick
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
get.transition.multiplexFAST <- function(Number_Proteins, Number_Layers, lambda, SupraAdjacencyMatrix, SupraBipartiteMatrix) {
  
  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrix, na.rm=FALSE, sparseResult=FALSE)
  Row_Sum_Bipartite <- rowSums(SupraBipartiteMatrix, na.rm=FALSE, sparseResult=FALSE)
  
  # Columns where bipartite row sum != 0 get (1-lambda) scaling
  scale_factors <- ifelse(Row_Sum_Bipartite != 0, (1 - lambda), 1) / Col_Sum_Multiplex
  
  Transition_Multiplex_Network <- SupraAdjacencyMatrix %*% Diagonal(x = scale_factors)
  
  rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
  colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)
  
  return(Transition_Multiplex_Network)
}
# 11.4.-DiseaseSimilarity intra-transition Matrix.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
get.transition.disease <- function(Number_Diseases,lambda,AdjMatrix,SupraBipartiteMatrix){
  
  Transition_Disease_Network <- Matrix(0,nrow=Number_Diseases,ncol=Number_Diseases,sparse = TRUE)
  
  rownames(Transition_Disease_Network) <- rownames(AdjMatrix)
  colnames(Transition_Disease_Network) <- colnames(AdjMatrix)
  
  Col_Sum_Disease <- colSums (AdjMatrix,na.rm=FALSE,dims=1, sparseResult=FALSE)
  Col_Sum_Bipartite <- colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
  
  for (j in 1:Number_Diseases){
    if(Col_Sum_Bipartite[j] != 0){
      Transition_Disease_Network[,j] <- ((1-lambda)*AdjMatrix[,j]) /Col_Sum_Disease[j]
    } else {
      Transition_Disease_Network[,j] <- AdjMatrix[,j] /Col_Sum_Disease[j]
    }
  }
  return(Transition_Disease_Network)
}