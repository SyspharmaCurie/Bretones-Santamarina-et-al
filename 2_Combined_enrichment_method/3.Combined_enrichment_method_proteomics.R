############################### SCRIPT TO MERGE GSEA AND RONTOTOOLS PROTEOMICS RESULTS USING THE GEOMETRIC MEAN OF RANKS
#install packages, if required
if (!require("tidyverse", character.only = TRUE)) {
  install.packages("tidyverse")
}
if (!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
}
if (!require("rlist", character.only = TRUE)) {
  install.packages("rlist")
}

if (!require("rrapply", character.only = TRUE))
  install.packages("rrapply")

if (!require("psych", character.only = TRUE))
  install.packages("psych")

if (!require("GSA", character.only = TRUE)) {
  install.packages("GSA")}


#Load required packages
library(tidyverse)
library(rlist)
library(data.table)
library(rrapply)
library(psych)
library(GSA)

#Set current directory to script directory in R Studio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

options(java.parameters = "- Xmx8000m")

#Get current dir
pdir <- getwd()


################## LOAD THE DATA

database_array <- list(reactome = "reactome")

#Load the data for proteomics

##GSEA
paths_proteomics_gsea <- "Outputs/gsea_proteomics_results.rds"

#ROntoTools
paths_proteomics_topology <- "Outputs/ROntoTools_proteomics_results.rds"

#Load the GMT files that were used for the enrichment
path_reactome_gmt <- "Inputs/reactome_gmt_symbol_twice_pruned.gmt"
path_reactome_table <- "Inputs/Reactome_table_twice_pruned_genes.rds"

############################ DEFINE FUNCTIONS TO PROCESS THE ENRICHMENT RESULTS OF EACH METHOD INTO A COMMON FORMAT

#GSEA (gene set enrichment analysis) results
#Order by padj and keep leading edge genes (this gives a list of the most important genes to compute Enrichment Score (ES))
process_gsea <- function(x,analysis){
   temp <- x %>%
      as.data.frame() %>%
      #arrange(!!as.name(p_measure)) %>%
      arrange(padj) %>%
      mutate(abs_NES = abs(NES))
   #Add ranks based on padj and abs(NES)
   temp$rank <- frank(temp,padj,-abs_NES) 
   temp %>%
      #Change name of rank
      rename(!!paste0("rank_gsea_",analysis) := rank) %>%
      #We call "p" either the p-value or the adjusted p-value
      #rename(!!paste0("p_gsea_",analysis) := !!as.name(p_measure)) %>%
      rename(!!paste0("p_gsea_",analysis) := padj) %>%
      rename(!!paste0("NES_gsea_",analysis) := NES) %>%
      rename(!!paste0("leadingEdge_gsea_",analysis):= leadingEdge) %>%
      distinct(pathway_name,.keep_all = T) %>%
      as_tibble()
}

#ROnToTools (TOPOLOGY)

#VERY IMPORTANT! There are pathways where all columns are NA. Those are going to get automatically the worst rank, which is bad when merging
#with GSEA. Therefore we will make the rank = NA, instead of the worst one
#This function can also analyse transcriptomic data
process_topology <- function(x,analysis){
   #Create a sublist called reactome where we put the results and another empty called bioplanet
   temp <- x %>%
      #We turn into a dataframe because if not rleid does not work
      as.data.frame() %>%
      #Use pPert.fdr as FDR
      arrange(`pPert.fdr`,desc(totalPertNorm)) 
      #Add ranks based on pPert.fdr and totalPertNorm
   temp$rank <- frank(temp,pPert.fdr,-totalPertNorm)
   temp %>%
      mutate(rank = ifelse(is.na(totalPertNorm),NA,rank)) %>%
      rename(!!paste0("rank_ROntoTools_",analysis) := rank) %>%
      #mutate(!!paste0("rank_camera-",algorithm,"_",analysis) := rleid(FDR)) %>%
      #We call "p" either the p-value or the adjusted p-value
      rename(!!paste0("p_ROntoTools_",analysis) := `pPert.fdr`) %>%
      distinct(pathway_name,.keep_all = T) %>%
      #Transform back to tibble
      as_tibble()
}

############################ READ PROTEOMICS DATA FOR EACH ENRICHMENT METHOD

#GSEA (gene set enrichmentanalysis) results
proteomics_gsea <- readRDS(paths_proteomics_gsea)
proteomics_gsea <-map_depth(proteomics_gsea,.depth = 2,.f = process_gsea,analysis = "proteome")

#ROntoTools (topology)
proteomics_topology <- readRDS(paths_proteomics_topology)
#IMPORTANT! We use map and not map_depth because for topology enrichment we only tested REACTOME.
#Bioplanet is not available as a graph object that we can test with ROntoTools
proteomics_topology <- map_depth(proteomics_topology,.depth = 2,.f = process_topology,analysis = "proteome")

########## GET THOSE PATHWAYS WITH NA RESULTS IN RONTOTOOLS

#These pathways cannot be analysed with ROntoTools and systematically give NA
#The reasons is that the edges do not include activations, inhibitions, expressions or repressions

map(.x = proteomics_topology,.f = ~.x$reactome %>% filter(is.na(rank_ROntoTools_proteome))) %>%
   bind_rows() %>%
   distinct(pathway_name,pathway_id) %>%
   print(n = 100)

########################### LOAD TABLES WITH GENES TO GET THE PATHWAYS THAT CONTAIN THE GENES IN OUR KO MUTANTS

reactome_table <- readRDS(path_reactome_table)

#Get names of all mutants from proteomics data
mutant_names <- map(.x = names(proteomics_gsea),.f = ~.x) %>% setNames(names(proteomics_gsea))

#Search for those names in reactome and get pathways containing them
pathways_with_KO_genes <- map(.x = mutant_names,.f = function(x){
   reactome_table %>%
      group_by(pathway_name) %>%
      mutate(contains_KO_gene = ifelse(x %in% gene_name,T,F)) %>%
      distinct(pathway_name,contains_KO_gene,length) %>%
      mutate(line = x,database = "reactome") 
}) 

#Collapse it
pathways_with_KO_genes_collapsed <- map(.x = pathways_with_KO_genes,.f = ~bind_rows(.x)) %>% bind_rows()

################## CALCULATION OF RANKINGS, GEOMETRIC MEAN OF PADJ AND GEOMETRIC MEAN OF RANKS FOR EACH DATABASE

#IMPORTANT! The ranks themselves are useless, as they are derived from the adjusted p-value of each method and alone they
#do not give hints about the significance (E.g.: a pathway may be ranked 1 and still have a very high p-value). 
#Hence, the solution is to average the adjusted p-values and then use a significance threshold (5%) to decide which pathways are overall
#significant in our enrichment analysis.
#The simple mean is a very bad way of merging p-values, as the putative range of p-values is very wide and it ignores large
#differences between numbers. The fisher method says that p-values must be independent, which are not. The harmonic mean and the
#geometric mean are better choices when there is a high degree of dependency between p-values. Between both, the geometric
#mean preserves best the differences between extreme numbers, so it will be our choice.

#library("psych")

#First, we take the list of pathways in the gmt files and we sort them alphabetically
#library(GSA)
path_gmt_list <- list(reactome = path_reactome_gmt)
gmt_list <- map(.x = path_gmt_list,.f = ~sort(GSA.read.gmt(.x)$geneset.names))
gmt_list <- gmt_list["reactome"]

#Create a list with all enrichment results
all_list <- list(proteomics_gsea = proteomics_gsea,proteomics_topology = proteomics_topology)

#Get only results of the databases in database_array
all_list <- map_depth(.x = all_list,.depth = 2,.f = ~.x[names(database_array)])

#Eliminate NULL things out of the nested list structure

all_list <- rrapply(all_list, condition = Negate(is.null), how = "prune")

#Save 
all_list_only_padj <- map_depth(.x = all_list,.depth = 3,.f = ~.x %>% select(pathway_name,contains(c("p_","NES_","rank","lead"))) %>% distinct())

########################### TRANSFORM ALL_LIST TO MERGE ALL 8 ENRICHMENT METHODS IN MUTANT AND FOR EACH DATABASE

#Create 2 lists: one with mutant and the other with database names
database_names <- database_array

########### MERGE RESULTS AND COMPUTE THE GEOMETRIC MEAN OF RANKS FOR PROTEOMICS

compute_res <- function(my_list){
   map(.x = mutant_names,.f = function(y){
      #For each database x, we go through each mutant 
      imap(.x = database_names,.f = function(z,datab_name){
         #For each mutant, we take the corresponding results and operate with them
         temp_data <- map(.x = my_list,.f = ~.x[[y]][[z]])
         #Eliminate any NULL entry (for example, for the BAP1 mutant, in RNA seq enrichments we will get NULL)
         temp_data <- temp_data[!sapply(temp_data,is.null)]
         #Append data.frame with gene_list from corresponding database
         #Like that, we can join all 4 methods starting with the list of all pathways tested for enrichment. This way,
         #we can use a left_join and no information will be lost. If a pathway is not present for a certain method, there will
         #be a NA value in the corresponding columns (padj,p_value...)
         current_gmt <- gmt_list[[z]]
         temp_data<- c(list(datab = data.frame(pathway_name = current_gmt,stringsAsFactors = F)),temp_data)
         names(temp_data)[1] <- paste0(z,"_gmt")
         #x is the temp list with 5 sublists
         joined_data <- temp_data %>% 
            reduce(left_join, by = "pathway_name") 
         #If we do not eliminate ORA, then the geometric mean of rank will include it
         joined_data$geometric_mean_rank_proteome <- apply(select(joined_data,intersect(starts_with("r"),ends_with("proteome"))),1,geometric.mean,na.rm = T)
         joined_data %>%
            mutate(database = datab_name) %>%
            as_tibble()
      })
   })
}

#Include all columns of results or only the padj + NES + direction columns (interesting information)

include_all <- F
if(include_all){
   include_all_text <- "all"
   res <- compute_res(all_list)
}else{
   include_all_text <- "only_padj"
   res <- compute_res(all_list_only_padj)
}

################################ TURN JOINED DATA INTO TIDY STRUCTURE TO MANIPULATE IT MORE EASILY

#Make it a tidy data structure

#Create function to turn data into tidy structure
make_tidy <- function(x){
   x %>%
      #We need to GATHER NES, direction and padj into a SINGLE VARIABLE
      pivot_longer(cols = starts_with(c("p_","NES","rank")),
                   names_to = c("variable","enrichment","omics"),
                   names_pattern = "^(.*)\\_{1}(.*)\\_(proteome)",
                   values_to = "values") %>%
      #Once we have gathered we can SPREAD the variable column to retrieve padj, NES and direction as separate cols
      #IMPORTANT! This 2 step process is necessary to have direction and NES columns with NA values for the rest of 
      #the enrichments that are not camera and GSEA (respectively). If we only do a GATHER with padj only, then
      #NES and direction will appear in all columns, no matter the type of enrichment
      pivot_wider(names_from = variable,values_from = "values") %>%
      #Pivot also geometric_mean_ columns
      pivot_longer(cols = starts_with("geometric_mean"),
                   names_to = "omics2",
                   names_pattern = "^geometric_mean_rank_(.*)$",
                   values_to = "geometric_mean_rank") %>%
      #Now we have duplicated rows, so we will only select those rows where omics == omics2
      filter(omics == omics2) %>%
      #Eliminate omics2 column, as the information is already contained in omics
      select(-omics2) %>%
      rename(padj = p) %>%
      mutate(enrichment = ifelse(enrichment == "gsea","GSEA",enrichment)) %>%
      #Reorder levels of factors
      mutate(enrichment = factor(enrichment,levels = c("GSEA","ROntoTools"))) %>%
      mutate(omics = factor(omics,levels = "proteome")) %>%
      #Turn pathway_name into a factor and preserve order from left to right by geometric mean
      mutate(pathway_name = factor(pathway_name,levels = unique(pathway_name))) %>%
      #Order by increasing geometric mean of rank
      arrange(geometric_mean_rank)
}

############# TURN ALL RESULTS INTO TIDY FORMAT AND SAVE

#The output is a list containing the combined enrichment results for each mutant as a sublist
#Inside each tibble we have every pathway duplicated in tidy format, as we can isolate GSEA and ROntoTools statistics
#The geometric mean of ranks is displayed in the geometric_mean_rank column
#The total perturbation from ROntoTools was used to rank pathways but it is not interesting, so it was discarded
#Instead, the NES score from GSEA indicates pathway up or downregulation, so it was kept

joined_data_tidy <- map_depth(.x = res,.depth = 2,.f = make_tidy)

#Save joined_data tidy for other scripts
saveRDS(joined_data_tidy,file = paste0("Outputs/Combined_enrichment_results.rds"))

############ SELECT SIGNIFICANT PATHWAYS BY AT LEAST ONE ENRICHMENT METHOD (OR RULE)

joined_data_significant_tidy <- map_depth(.x = joined_data_tidy,.depth = 2,.f = function(x){
   x %>%
      group_by(pathway_name) %>%
      filter(any(padj < 0.05)) %>%
      ungroup()
})

