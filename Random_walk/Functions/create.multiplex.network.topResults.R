### Function to generate a network with the TOP candidates of the random walk. 

create.multiplex.network.topResults <- function(Network_List,All_Seeds,Top_Results_Genes){
  
  print("Creating Network with Top results...")
  
  Number_Networks <- length(Network_List)
  
  ###### We read all the monolpex network and we assign them some features.
  PPI_table <- data.frame()
  Pathway_table <- data.frame()
  Coexpresion_table <- data.frame()
  error <- FALSE
  
  for (i in 1:Number_Networks){
    if (Network_List[i] == "GENES"){
      PPI_table <- read.table("Networks/gene_network.tsv",sep="\t")
      PPI_table$type = "PPI"
      PPI_table$color = "blue"
      
    } else {
      if (Network_List[i] == "METABOL"){
        Pathway_table <- read.delim("Networks/metabolomics_network.tsv")
        Pathway_table$type = "Pathway"
        Pathway_table$color ="orange"
        
      } else {
        if (Network_List[i] == "COEX"){
          Coexpresion_table <- read.table("Networks/Co-Expression_2016-11-23.gr", sep= " ")
          Coexpresion_table$type = "Coexpresion"
          Coexpresion_table$color ="red"
          
        } else {
          stop(paste("Not valid Network Introduced" ,Network_List[i], ". Please introduce one of those: PPI, PATH, COMP, COEX"))
          error <- TRUE
        }
      }
    }
  }
  

  #Multiplex_Network_df <- rbind(PPI_table,Pathway_table,Coexpresion_table)
  Multiplex_Network_df <- rbind(PPI_table,Pathway_table)
  Multiplex_Network <- graph_from_data_frame(Multiplex_Network_df,directed=FALSE)
  
  ## We include in the network the seeds.
  Query_Genes <- c(All_Seeds,Top_Results_Genes)
   
  Induced_Network <- dNetInduce(g=Multiplex_Network, nodes_query=Query_Genes, knn=0, remove.loops=F, largest.comp=F)
  return(Induced_Network)
}