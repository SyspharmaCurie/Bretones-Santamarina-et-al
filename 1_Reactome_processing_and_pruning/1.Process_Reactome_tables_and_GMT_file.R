########################################   PROCESS AND TIDY REACTOME GMT FILE AND CREATE TABLE OF PATHWAYS, GENES AND CATEGORIES

#Install packages, if required
if (!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
}
if (!require("tidyverse", character.only = TRUE)) {
  install.packages("tidyverse")
}

if (!require("readr", character.only = TRUE)) {
  install.packages("readr")
}
if (!require("GSA", character.only = TRUE)) {
  install.packages("GSA")
}
if (!require("rlist", character.only = TRUE)) {
  install.packages("rlist")
}
#Load packages
library(data.table)
library(tidyverse)
library(readr)
library(GSA)
library(rlist)

#Set current directory to script directory in R Studio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Get current dir
pdir <- getwd()


#Load mapping of all pathways and their Reactome identifier
name_identifier <- read_tsv(file = paste0(pdir,"/Inputs/ReactomePathways.txt"),col_names = c("pathway_id","pathway_name","species"))

#Load mapping of each parent pathway to the child in the pathway hierarchy
mapping_tree <- read_tsv(file = paste0(pdir,"/Inputs/ReactomePathwaysRelation.txt"),col_names = c("parent_id","child_id"))
head(mapping_tree)

#The idea is that the HIGHEST ORDER HIERARCHY is represented in Reactome as a pathway with a certain identifier. 
#Hence, we can use that information to climb up the hierarchy tree and get the highest hierarchy to which each pathway belongs.
#But first, let's eliminate information about other species which are not Homo sapiens:
name_identifier <- name_identifier %>%
    filter(species == "Homo sapiens") %>%
    as.data.frame()

#To check duplicated pathways with 1 name but more than 1 ID
name_identifier %>% group_by(pathway_name) %>% filter(n()>1) %>% distinct(pathway_name) %>% pull()

#Use regexp to get only those rows containing HSA in the name:
mapping_tree <- mapping_tree %>%
    filter(grepl(pattern = "-HSA-",x = parent_id)) %>%
    as.data.frame()

#Now, we first get those parent ids that ARE NOT childs. Those are the highest level pathways we need:
unique_parent_ids <- unique(mapping_tree$parent_id)
#We go through each pathway and check they are not children of any other pathway. Then, they are highest order pathways that define a category
highest_order_parents <- map_chr(.x = unique_parent_ids,.f = function(i){
    index_child <- which(mapping_tree$child_id == i)
    #If there is no match, we want it
    ifelse(any(index_child) == 0,i,NA)
}) 
#Select those that are not NA and map to retrieve the pathway_names corresponding to those pathways
highest_order_parents <- data.frame(pathway_id = highest_order_parents[!is.na(highest_order_parents)]) %>%
    left_join(name_identifier,by = "pathway_id")

#Now we can just follow the tree down and store any pathway that belongs to that category:
#First, we create a store list and initialise each category with the pathway_id that defines it:
store_list <- lapply(highest_order_parents$pathway_id,function(x) x) %>% setNames(highest_order_parents$pathway_name)
#Map starting from each of them and retrieve the children one by one
for(i in 1:length(highest_order_parents$pathway_id)){
    #If there is still one child, it means we did not get to the bottom
    id_child = 1
    input <- highest_order_parents$pathway_id[i]
    #While there is at least 1 child, we keep going down the tree
    while (length(id_child) > 0){
        #Get index of children of that parent
        index_parent <- which(mapping_tree$parent_id %in% input)
        id_child <- mapping_tree$child_id[index_parent]
        #If there is at least 1 child, store it
        if(length(id_child) > 0){
            #Store the child
            store_list[[i]] <- c(store_list[[i]],id_child)
            #The new input is the children we just stored
            input <- id_child
        }
    }
    #Once we got to the bottom of 1 category, we go for the next in the for loop
}
str(store_list)

#Check the sum of all elements
sum(sapply(store_list,length)) 

#Translate those identifiers into pathway names and add category names, retrieving a data.frame inside each category list
reactome_table <- imap(.x = store_list,.f = function(pathway_ids,category_name){
    tibble(pathway_id = pathway_ids) %>%
        left_join(name_identifier,by = "pathway_id") %>%
        mutate(category = category_name)
}) %>% bind_rows() %>% mutate(category = factor(category))

############################# ELIMINATE REPLICATES INSIDE THE SAME CATEGORY

#There are pathways, inside the same category, that are replicated. 
#This means that they appear as 2 different nodes of the hierarchy tree inside the same category.
#Example: "polymerase switching"
reactome_table_no_duplicates <- reactome_table[!duplicated(reactome_table),]

########################### ELIMINATE PATHWAYS WITH MULTIPLE PATHWAY IDs

#We need as well to eliminate those pathways that are repeated. 
#There are pathways like "Maturation of protein E" that appear under the same name, with different pathway ids. 
#This does not make any sense, so we will just pick one of both when they are duplicated.
#This pipeline does not employ pathway IDs and focuses on name, so we do not care what ID is chosen for these pathways.

#Print all examples
pathways_duplicated_ID <- reactome_table_no_duplicates %>% 
   distinct(pathway_id,pathway_name) %>% 
   group_by(pathway_name) %>% 
   filter(n()>1) %>% 
   arrange(pathway_name) 

#Get second pathway name for elimination
path_dup_ID_remove <- pathways_duplicated_ID %>%
   filter(row_number()==2)

#Remove duplicates
reactome_table_no_duplicates <- reactome_table_no_duplicates %>%
   anti_join(path_dup_ID_remove)

############################ UPLOAD REACTOME GMT FILE AND ELIMINATE GENES WITH SPACES IN THEIR NAMES

temp_gmt <- GSA.read.gmt(paste0(pdir,"/Inputs/ReactomePathways.gmt"))

#Obtain descriptions and make the link with pathway names
descriptions_and_names <- data.frame(Pathway = temp_gmt$geneset.names,description = temp_gmt$geneset.descriptions,stringsAsFactors = F)
#Now we change the format in temp_gmt to a list of lists where each sublist is a pathway and inside of it we got the genes
reactome_gmt <- temp_gmt$genesets
names(reactome_gmt) <- temp_gmt$geneset.names

#We eliminate from our gene sets those genes with spaces, which are not allowed
reactome_gmt <- sapply(reactome_gmt,function(x){
   #We eliminate all entries containing spaces
   indexes <- grep("[[:space:]]",x)
   if(length(indexes) > 0){
      out <- x[-indexes]
   }else{
      out <- x
   }
   invisible(out)
})

reactome_gmt_final <- list()
reactome_gmt_final$genesets <- reactome_gmt

#We add the descriptions
reactome_gmt_final <- mapply(c,temp_gmt$geneset.descriptions, reactome_gmt_final, SIMPLIFY=FALSE)

#Add correct list names
names(reactome_gmt_final) <- temp_gmt$geneset.names

############################ ADD THE GENES THAT DEFINE EACH PATHWAY (REPLICATE ROWS) FROM GMT FILE

#We build a table with pathway ids and genes
temp_genesets <- lapply(reactome_gmt,function(x) data.frame(gene_name = x,stringsAsFactors = F))
names(temp_genesets) <- unlist(temp_gmt$geneset.descriptions)
reactome_gmt_table <- list.stack(temp_genesets,idcol = "pathway_id",use.names = F)

#Check number of pathways in reactome_gmt_table
reactome_gmt_table %>% distinct(pathway_id) %>% as_tibble()

########################### DISCARD PATHWAYS IN MAPPING FILES BUT NOT IN THE GMT FILE

#We need to discard them, as only pathways in the GMT file have gene information to run pathway enrichment
reactome_table_no_duplicates %>% anti_join(reactome_gmt_table,by = "pathway_id") %>% as.data.frame()

#Eliminate those pathways that are not in the GMT file
reactome_table_no_duplicates <- reactome_table_no_duplicates %>%
    semi_join(reactome_gmt_table,by = "pathway_id")

#Join reactome_table_no_duplicates to the GMT file to retrieve the information about a pathway (name, category) 
reactome_table_genes <- reactome_gmt_table %>%
    left_join(reactome_table_no_duplicates,by = "pathway_id",relationship = "many-to-many") 

########################### DISCARD PATHWAYS IN GMT FILE THAT CANNOT BE MAPPED TO REACTOME CATEGORIES

#These pathways are new pathways that cannot be automatically mapped to a category, so we eliminate them
#Some of them are related, for example, to a new category called Drug ADME

reactome_table_genes %>%
    filter(is.na(category)) %>%
    distinct(pathway_id)

#We eliminate them
reactome_table_genes <- reactome_table_genes %>%
    filter(!is.na(category)) 

##################### ADD PATHWAY LENGTH EXPRESSING THE NUMBER OF GENES

reactome_table_genes <- reactome_table_genes %>%
   group_by(pathway_name,category) %>%
   mutate(length = n_distinct(gene_name)) %>%
   ungroup()

#Compute reactome_table_final by just eliminating the genes from reactome_table_genes

reactome_table_final <- reactome_table_genes %>%
   distinct(pathway_id,pathway_name,length,species,category)

#Finally, we update the GMT file to eliminate those pathways that do not appear in the pathway relationships (tree)

################### STORE FINAL GMT FILE WITH GENE SYMBOLS

#We add the descriptions
#We need to split by an ordered factor. If not, split will reorder the IDs in alphabetic order and then the order of pathways
#will be changed with respect to the original gmt file
#We need to eliminate the categories because if not, we will get repeated gene_names in certain pathways
new_reactome_gmt <- reactome_table_genes %>% arrange(pathway_name) %>% distinct(pathway_id,gene_name) %>% split(f = factor(.$pathway_id, levels = unique(.$pathway_id))) %>% imap(.f = ~c(.y,.x$gene_name))

#Add correct list names
names(new_reactome_gmt) <- data.frame(pathway_id = names(new_reactome_gmt)) %>% left_join(distinct(reactome_table_genes,pathway_id,pathway_name),by = "pathway_id") %>% pull(pathway_name)

#Store the final gmt file
#Get named array with rows of the gmt file
my_names_symbol <- sapply(names(new_reactome_gmt),function(x) paste(x,paste(new_reactome_gmt[[x]],collapse="\t"),sep = "\t"))
write(file = paste0(pdir,"/Outputs/reactome_gmt_symbol_without_pathways_not_in_tree.gmt"),x = my_names_symbol, ncolumns=10000,sep = "\n")
    
################### STORE FILES

#Store reactome table final without genes
saveRDS(reactome_table_final,file = paste0(pdir,"/Outputs/Reactome_table_all_categories.rds"))
#Store reactome table with genes
saveRDS(reactome_table_genes,file = paste0(pdir,"/Outputs/Reactome_table_all_categories_with_genes.rds"))

