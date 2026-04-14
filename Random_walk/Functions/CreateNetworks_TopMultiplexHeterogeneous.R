#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
###                       NETWORK GENERATION                               ####  
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
############# 1.- We define the function which is going to read the different networks.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

CreateNetworks_TopMultiplexHeterogeneous <- function(Network_List,All_Seeds,Top_Results_Genes,Top_Results_Diseases ){
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  ####### 1.1.- MULTIPLEX NETWORK 
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  
  Number_Networks <- length(Network_List)
  
  #### We do some sanity checks. 
  try(if(Number_Networks > 4) stop("Too many networks"))
  #try(if(Number_Networks < 1) stop("Introduce at least 1 network"))
  try(if(anyDuplicated(Network_List)) stop("Duplicated Networks Names"))
  
  ###### We read all the networks and we assign to their interactions an specific weight
  PPI_table <- data.frame()
  Pathway_table <- data.frame()
  Coexpresion_table <- data.frame()
  error <- FALSE
  for (i in 1:Number_Networks){
    if (Network_List[i] == "GENES"){
      PPI_table <- read.table("Networks/gene_network.tsv",sep="\t")
      PPI_table <- PPI_table[, -3]
      PPI_table$type = "G"
      PPI_table$color = "blue"
      
    } else {
      if (Network_List[i] == "METABOL"){
        Pathway_table <- read.delim("Networks/metabolomics_network.tsv")
        Pathway_table$type = "M"
        Pathway_table$color ="Yellow"
        
      } else {
        if (Network_List[i] == "COEX"){
          Coexpresion_table <- read.table("Networks/Co-Expression_2016-11-23.gr", sep= " ")
          Coexpresion_table$type = "Coexpresion"
          Coexpresion_table$color ="orange"
          
        } else {
          stop(paste("Not valid Network Introduced" ,Network_List[i], ". Please introduce one of those: PPI, PATH, COMP, COEX"))
          error <- TRUE
        }
      }
    }
  }
  
  # Multiplex_Network_df <-  rbind(PPI_table,Pathway_table,Complexes_table,Coexpresion_table)
  # Multiplex_Network <- graph.data.frame(Multiplex_Network_df,directed=FALSE)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  ####### 1.2.- DISEASE-DISEASE SIMILARITY NETWORK.
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  Disease_table <- read.table("Networks/DiseaseSimilarity_2016-12-06.gr", sep= " ")
  Disease_table$V1 <- as.character(Disease_table$V1)
  Disease_table$V2 <- as.character(Disease_table$V2)
  Disease_table$type = "Disease"
  Disease_table$color ="black"
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  ####### 1.3.- BIPARTITE GRAPH.
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  Gene_Phenotype_relation <- read.table("Input_files/Gene_Phenotype_relation.txt", 
                                        sep="", header=FALSE,stringsAsFactors = FALSE,skip = 1)
  Gene_Phenotype_relation$type <- "Bipartite"
  Gene_Phenotype_relation$color <- "grey30"
  
  Gene_Phenotype_relation <- Gene_Phenotype_relation[which(!is.na(Gene_Phenotype_relation$V2 )),]
  
  Gene_Phenotype_relation$V2 <- as.character(Gene_Phenotype_relation$V2)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  ####### 1.4.- WE MERGE THE NETWORKS
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  print("Merging all the networks...")
  Multiplex_Heterogeneous_Network_df <- rbind(PPI_table,
                                              Disease_table,Gene_Phenotype_relation)
  Multiplex_Heterogeneous_Network <- graph.data.frame(Multiplex_Heterogeneous_Network_df,directed=FALSE)
  
  
  Query_Genes <- c(All_Seeds,Top_Results_Genes,Top_Results_Diseases)
  
  
  print("Inducing Network with top results...")
  Induced_Network <- dNetInduce(g=Multiplex_Heterogeneous_Network, nodes_query=Query_Genes, knn=0, remove.loops=F, largest.comp=F)
  
  return(Induced_Network)
  
}





