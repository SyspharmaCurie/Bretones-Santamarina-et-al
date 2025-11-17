##################### THIRD GENERATION TOPOLOGY-BASED ENRICHMENT WITH RONTOTOOLS
#install packages, if required

if (!require("tidyverse", character.only = TRUE)) {
  install.packages("tidyverse")
}

if (!require("readr", character.only = TRUE)) {
  install.packages("readr")
}
if (!require("GSA", character.only = TRUE)) {
  install.packages("GSA")}
  
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("ROntoTools", character.only = TRUE))
  BiocManager::install("ROntoTools")

if (!require("graphite", character.only = TRUE))
  BiocManager::install("graphite")

if (!require("graph", character.only = TRUE))
  BiocManager::install("graph")

if (!require("org.Hs.eg.db", character.only = TRUE))
  BiocManager::install("org.Hs.eg.db")

#load packages

library(ROntoTools)
library(graphite)
library(graph)
library(tidyverse)
library(GSA)
library(org.Hs.eg.db)

#Set current directory to script directory in R Studio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Get current dir
pdir <- getwd()

#################### PROCESS THE REACTOME GRAPH FORMAT TO USE IT AS INPUT FOR RONTOTOOLS

#The ROntoTools documentation uses the KEGG database as pathway inputs. In this script, we will retrieve the Reactome
#graph database and adapt it to be used for ROntoTools


############################# UPLOAD GMT FILES

path_reactome_table <- "Inputs/Reactome_table_twice_pruned.rds"

reactome_table <- readRDS(path_reactome_table)

################################# 1. PROCESS THE REACTOME GRAPH DATABASE

#We load the GMT files and use them to select only the graph objects corresponding to interesting pathways.
#By default, graphite will retrieve the entire reactome database
min_size <- 10
max_size <- 500
reactome_pathway_ids <- unique(reactome_table$pathway_id)
reactome_pathway_names <- unique(reactome_table$pathway_name)

#First we retrieve all Reactome pathways for Human in the graph database using graphite's pathways function
humanReactome <- pathways("hsapiens", "reactome")
humanReactome

#Selecting elements of the list by pathway name is tricky, because some of the characters (e.g., greek letters) are wronly encoded in
#our table or in theirs, so we have no match. Instead we will subtract the pathway ids and compare to our list to
#match the indices and select the good pathways
reactome_graph_pathway_ids <- sapply(X = humanReactome,FUN = function(x){pathwayId(x)})
#Use match
indices_matching <- match(reactome_pathway_ids,reactome_graph_pathway_ids)
#Check if there are NA matches, which would mean that graphite has not uploaded a pathway in our table
sum(is.na(indices_matching))

blind_pathways_reactome <- sort(reactome_pathway_names[is.na(indices_matching)])

#Neutrophil degranulation does not have a topological version, so we eliminate it
indices_matching <- indices_matching[!is.na(indices_matching)]
human_Reactome_filtered <- humanReactome[indices_matching]
human_Reactome_filtered
reactome_graph_filtered_pathway_ids <- sapply(X = human_Reactome_filtered,FUN = function(x){pathwayId(x)})

#We need to change pathway_names by pathway_ids and run enrichment on IDS. This is the only way to map to our reactome_table
#This is because there is probably a shift of versions between the graph package and the Reactome database, so pathway_names change
names(human_Reactome_filtered@entries) <- reactome_graph_filtered_pathway_ids

#The identifiers for the nodes and the edges are in UNIPROT codes. However, we need to convert to gene SYMBOL
#to perform the gene enrichment analysis
reactome_symbols <- convertIdentifiers(human_Reactome_filtered, "SYMBOL")

#Convert the list of pathways to a graphNEL object
human_Reactome_filtered_graph <- lapply(X = reactome_symbols,FUN = function(x){pathwayGraph(x)})

#Identify those pathways that have 0 nodes (e.g., "Reversible hydration of carbon dioxide")
#If those pathways only contain metabolites and no proteins, then we need to eliminate them
pathways_zero_nodes <- sapply(X = human_Reactome_filtered_graph,FUN = function(x){ifelse(length(x@nodes) == 0,T,F)})
human_Reactome_filtered_graph <- human_Reactome_filtered_graph[!pathways_zero_nodes]
length(human_Reactome_filtered_graph)

#Out of the 14 edge types only 4 are used: activation = 1, inhibition = -1, expression = 1, repression = -1
#The rest of the edges are given weight = 0
#Retrieve all the edge types in reactome
#We will modify the edge information to substitute ";" by "," and to recode the edgeTypes into activation, inhibition
edge_types <- unique(unlist(sapply(human_Reactome_filtered_graph,FUN = function(pathway){
   unique(sapply(pathway@edgeData@data,FUN = function(pairs){
      pairs$edgeType
   }))
})))

#Separate by ; and get the unique building blocks
edge_types_sep <- unique(unlist(str_split(edge_types,pattern = ";")))
edge_types_sep

###################### RUN IN PARALLEL

#The overall data structure is a list, which allows as to use future_map to split the task between different workers
#However, inside each pathway, the data structure is a Bioconductor data type that we need to respect. 
#Using lapply or map returns lists so we can only use a normal for loop to modify the corresponding fields
human_Reactome_filtered_graph_edges <- human_Reactome_filtered_graph
library(mgsub)
mod_path <- function(curr_path){
   #For each edge, we modify it
   for (j in 1:length(curr_path@edgeData@data)){
      #cat(paste0("j = ",j,"\n"))
      #Retrieve the edgeType
      temp <- curr_path@edgeData@data[[j]]$edgeType
      #Substitute ";" by ","
      temp <- stringr::str_replace_all(temp,pattern = ";",replacement = ",")
      #Recode the edgeTypes to be recognised by setEdgeWeights()
      temp <-  mgsub(temp, pattern = c(
         "Control\\(Out\\: INHIBITION of BiochemicalReaction\\)",
         "Control\\(In\\: INHIBITION of BiochemicalReaction\\)",
         "Control\\(Out\\: ACTIVATION of BiochemicalReaction\\)",
         "Control\\(In\\: ACTIVATION of BiochemicalReaction\\)",
         "Control\\(Out\\: ACTIVATION of TemplateReaction\\)",
         "Control\\(In\\: ACTIVATION of Degradation\\)"), replacement = c(rep("inhibition",2),
                                                                          rep("activation",4)))
      #Substitute the value
      curr_path@edgeData@data[[j]]$edgeType <- temp
   }
   return(curr_path)
}

library(furrr)
plan(multisession(workers = 19))
human_Reactome_filtered_graph_edges = future_map(human_Reactome_filtered_graph,.f = mod_path)

#Check some edges
head(edgeData(human_Reactome_filtered_graph_edges$`R-HSA-418592`, attr = "edgeType"))
head(edgeData(human_Reactome_filtered_graph_edges$`R-HSA-429958`, attr = "weight"))

#Get again all unique edges to see everything was changed correctly
edge_types <- unique(unlist(sapply(human_Reactome_filtered_graph_edges,FUN = function(pathway){
   unique(sapply(pathway@edgeData@data,FUN = function(pairs){
      pairs$edgeType
   }))
})))
edge_types_sep <- unique(unlist(str_split(edge_types,pattern = ",")))
edge_types_sep

#Save graph database in RDS format to upload in the topology enrichment script
saveRDS(object = human_Reactome_filtered_graph_edges,paste0(pdir,"Outputs/Reactome_pruned_in_graph.rds"))

#################### UPLOAD DIFFERENTIAL EXPRESSION ANALYSIS OUTPUT

set_123 <- readRDS("Inputs/proteomics_DE_results_sets123.rds") 
set_4 <- readRDS("Inputs/proteomics_DE_results_set4.rds")

#Put into a list
proteomics_list <- c(set_123,set_4)

######################################### RUN ANALYSIS WITH REACTOME DATABASE ENCODED AS A GRAPH OBJECT

############################# UPLOAD GMT FILES

#Decide if we want the original or the pruned gmt files

path_reactome_gmt <- "Inputs/reactome_gmt_symbol_twice_pruned.gmt"
path_reactome_table <- "Inputs/Reactome_table_twice_pruned.rds"

reactome_table <- readRDS(path_reactome_table)
reactome_pathway_ids <- unique(reactome_table$pathway_id)
reactome_pathway_names <- unique(reactome_table$pathway_name)

###################################### DIRECTLY UPLOAD THE PROCESSED GRAPH DATABASE

reactome_graphs <- readRDS("Outputs/Reactome_pruned_in_graph.rds")

###### CREATE FUNCTION TO RUN THE ENRICHMENT

#ROntoTools can do an overrpresentation-like analysis where differentially expressed genes are previously selected based on 
#an adjusted p-value and/or log2FC cutoff. However, in order to avoid an arbitrary selection of genes, we will run the algorithm 
#using all genes, which will be weighted according to the differential expression metrics. 
#This makes the analysis comparable to GSEA, which also uses the whole ranked gene list.

run_topology_enrichment <- function(mutant_list,database,alpha){
   #Select differentially expressed genes with an adjusted p-value < alpha and save their log2FC in an array
   temp <- mutant_list %>%
      dplyr::filter(qval <= alpha)
   #Get Log2FoldChange of significant genes
   fc <- temp$logFC
   #Get names of genes with "SYMBOL:" in front of it, to match the graph format
   gene_names <- paste0("SYMBOL:",temp$gene_symbol_human)
   names(fc) <- gene_names
   #Get adjusted p-values
   padj <- temp$qval
   names(padj) <- gene_names
   #Get the reference
   ref <- paste0("SYMBOL:",mutant_list$gene_symbol_human)
   #Once we have modified the edgeTypes, we can use edgeWeights() to assign weights based on the edgeTypes
   current_datab_edges <- setEdgeWeights(database,edgeTypeAttr = "edgeType",
                                         edgeWeightByType = list(activation = 1, inhibition = -1,
                                                                 expression = 1, repression = -1),defaultWeight = 0)
   #Once the edgeWeights have been assigned, we add the node weights based on the adjusted p-values of the genes
   current_datab_nodes <- setNodeWeights(current_datab_edges, weights = alphaMLG(padj),defaultWeight = 1)
   #Perform the pathway analysis. To obtain accurate results the number of boostraps, nboot, should be increase to a number like 2000
   #according to the documentation. We will use 4000. Increasing it more takes a very long computational time.
   peRes <- pe(fc, graphs = current_datab_nodes, ref = ref,nboot = 4000, verbose = TRUE,seed = 123)
   #Extract results. We cannot save everything because it weights a lot in the memory, so we just call summary to obtain a dataframe
   summary(peRes,p.adjust.methods = "BH") %>%
      rownames_to_column(var = "pathway_id") %>%
      #We order by increasing PERTURBATION FDR and if there are similar values, by descending total Perturbation norm
      arrange(pPert.fdr,desc(totalPertNorm)) %>%
      as_tibble()
}

#We will select alpha = 1 to use all genes in the differential expression analysis. 
my_alpha <- 1
my_alpha_text <- str_replace(my_alpha,"\\.","")

res <- lapply(X = proteomics_list,FUN = function(mutant){
   #For each mutant, we run topological enrichment
   run_topology_enrichment(mutant,reactome_graphs,my_alpha)
})

###################### MAP REACTOME PATHWAY IDS TO PATHWAY NAMES

res_data_frame <- map(.x = res,.f = function(x){
   x <- x %>%
      #Join to Reactome table to retrieve pathway names
      left_join(reactome_table[,c("pathway_name","pathway_id")],by = "pathway_id") %>%
      select(pathway_name,-pathway_id,everything())
})

########################## SAVE RESULTS

saveRDS(res_data_frame,file = paste0("Ouputs/ROntoTools_proteomics_results.rds"))
