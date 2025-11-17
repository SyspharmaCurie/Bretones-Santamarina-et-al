############################## REACTOME: PRUNE PATHWAY TREE TO ELIMINATE REDUNDANT PATHWAYS MONITORING THE NUMBER OF GENES

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
if (!require("rlist", character.only = TRUE)) {
  install.packages("rlist")
}
if (!require("janitor", character.only = TRUE)) {
  install.packages("janitor")
}
if (!require("openxlsx", character.only = TRUE)) {
  install.packages("openxlsx")
}

#Load packages
library(data.table)
library(tidyverse)
library(readr)
library(rlist)
library(janitor)
library(openxlsx)

#Set current directory to script directory in R Studio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Get current dir
pdir <- getwd()
getwd()


#Load mapping of all pathways and their Reactome identifier
name_identifier <- read_tsv(file = paste0(pdir,"/Inputs/ReactomePathways.txt"),col_names = c("pathway_id","pathway_name","species"))

#Load mapping of each parent pathway to the child in the pathway hierarchy
mapping_tree <- read_tsv(file = paste0(pdir,"/Inputs/ReactomePathwaysRelation.txt"),col_names = c("parent_id","child_id"))

#Read reactome table with genes. We take it from Outputs because it is an output of the first script
table_with_genes <- readRDS(paste0(pdir,"/Outputs/Reactome_table_all_categories_with_genes.rds"))

#Compute pathway length
table_with_genes_no_categ <- table_with_genes %>%
   distinct(pathway_name,gene_name) %>%
   group_by(pathway_name) %>%
   mutate(length = n_distinct(gene_name))

#Add length to table_with_genes
table_with_genes <- table_with_genes %>%
   left_join(table_with_genes_no_categ)

#Create table with genes without categories
table_with_genes_no_categories <- table_with_genes %>%
   distinct(pathway_name,pathway_id,gene_name,length)

#Check pathway lengths
table_with_genes_lengths <- table_with_genes %>%
   distinct(pathway_name,pathway_id,length,gene_name) %>%
   group_by(pathway_name) %>%
   chop(gene_name) 

#Compute number of pathways per category
num_pathways_per_category <- table_with_genes %>%
   distinct(pathway_name,category) %>%
   group_by(category) %>%
   #Get unique genes
   distinct(pathway_name) %>%
   summarise(number_pathways = n())

#The highest order hierarchy is represented in Reactome as a pathway with a certain identifier. 
#Hence, we can use that information to climb up #the hierarchy tree and get the highest hierarchy to which each pathway belongs.

#But first, let's eliminate information about other species which are not Homo sapiens:
name_identifier <- name_identifier %>%
   filter(species == "Homo sapiens") %>%
   as.data.frame()
#Use regexp to get only those rows containing HSA in the name:
mapping_tree <- mapping_tree %>%
   filter(grepl(pattern = "-HSA-",x = parent_id)) %>%
   as.data.frame()

#We first get those parent ids that are not childs. Those are the highest level pathways we need:
unique_parent_ids <- unique(mapping_tree$parent_id)
highest_order_parents <- c()
#Now we go through each of them and check they are not childs:
for(i in unique_parent_ids){
   index_child <- which(mapping_tree$child_id == i)
   #If there is no match, we want it
   if(any(index_child) == 0){
      highest_order_parents <- c(highest_order_parents,i)
   }
}
#We map the pahtway_ids to pathway names:
highest_order_parents_tab <- data.frame(pathway_id = highest_order_parents) %>%
   left_join(select(name_identifier,-species),by = "pathway_id")

index <- match(highest_order_parents,name_identifier$pathway_id)
highest_order_parent_names <- name_identifier$pathway_name[index]
highest_order_parent_names

#The variable highest_order_parent_names contains pathway names representing full categories

#Function to map a pathway id to its genes. This function is made to SUM all genes coming from the children OR get the parent genes for each parent SEPARATELY
map_id_genes <- function(x){
   if(length(x) < 1){
      return(character(0))
   }
   x <- data.frame(pathway_id = x,stringsAsFactors = F,id = 1:length(x))
   res <- table_with_genes_no_categories[,c("pathway_id","gene_name")] %>%
      inner_join(x,by = "pathway_id") %>%
      distinct(gene_name) %>% pull()
}

#Use this fucntion to retrieve all information about the children. We will add this as a list in the results tibble. 
#The goal is to be able to verify that the algorithm is working correctly by accessing the children information at each step
get_children_info <- function(x,minimal_length,maximal_length){
   if(length(x) < 1){
      return(character(0))
   }
   x <- data.frame(pathway_id = x,stringsAsFactors = F,id = 1:length(x))
   res <- table_with_genes_no_categories %>%
      inner_join(x,by = "pathway_id") %>%
      group_by(pathway_name) %>%
      chop(gene_name) %>%
      ungroup() %>%
      rename(child_id = pathway_id,child_name = pathway_name) %>%
      mutate(child_length = map_int(.x = gene_name,.f = ~length(.x))) %>%
      #Add if every child is or not of correct size
      mutate(is_children_correct_size = ifelse(child_length >= minimal_length & child_length <= maximal_length,T,F))
}

#Create similar function to map_id_genes taking into account a minimal size threshold
#If some of the children is under this threshold, we will flag it for elimination, as we do not want to test it for enrichment
#This is a necessary function because in case a child has incorrect size and needs to be eliminated, there exists the possibility
#that when we eliminate it, the rest of its siblings still account for all the genes in the parent. 
map_id_genes_between_threshold <- function(x,min_threshold,max_threshold){
   #In this case, x is an array with the children's IDs and we do not want to paste their genes together
   #Instead, we need to check if any of those children is under the threshold for min size
   if(length(x) < 1){
      return(character(0))
   }
   x <- data.frame(pathway_id = x,stringsAsFactors = F,id = 1:length(x))
   res <- table_with_genes_no_categories[,c("pathway_id","gene_name")] %>%
      inner_join(x,by = "pathway_id") %>%
      #We add the number of genes of each pathway in res
      group_by(pathway_id) %>%
      mutate(size = n())
   #Get only the genes of those pathways that are above the min and below the max threshold
   res_between_threshold <- res %>%
      filter(size >= min_threshold & size <= max_threshold) %>%
      distinct(gene_name) %>% pull()
   #If the output is character(0) it means that all children were under the threshold, so the parent cannot be eliminated
}

#Function to map a pathway id to pathway name
map_id_name <- function(x){
   if(length(x) < 1){
      return(character(0))
   }
   x <- data.frame(pathway_id = x,stringsAsFactors = F,id = 1:length(x))
   res <- name_identifier[,c("pathway_id","pathway_name")] %>%
      inner_join(x,by = "pathway_id")  %>%
      distinct(pathway_name)
   }

#Check children size
check_children_correct_size <- function(x,minimal_length,maximal_length){
   #If the analysed pathway has NO CHILDREN, then we return F
   #This is because if the pathway is a tree leaf, as long as it has correct size, it will be kept
   if(length(x) < 1){
      return(F)
   }
   x <- data.frame(pathway_id = x,stringsAsFactors = F,id = 1:length(x))
   res <- x %>%
      left_join(select(table_with_genes_no_categories,pathway_id,gene_name),by = "pathway_id")  %>%
      #Eliminate possible duplicated genes
      distinct(pathway_id,id,gene_name) %>%
      #We add the number of genes of each pathway in res
      group_by(pathway_id) %>%
      summarise(size = n()) %>%
      #Create column saying if the pathway is WRONG SIZE
      mutate(is_correct_size = ifelse(size >= minimal_length & size <= maximal_length,T,F)) 
   #If all children have correct size, we output T. If ANY of them is of incorrect size, we output F
   ifelse(all(res$is_correct_size) == T,T,F)
}

#Function to compare parent gene list to the union of children genes (after eliminating children of wrong size)
check_parent_is_contained <- function(parent_list,children_list){
   o <- setdiff(parent_list,children_list)
   #If o is greater than 0, it means that the parent list has more genes. 
   #If parent_list is empty, then it is not contained in the children because the parent has no genes
   if(length(o) > 0 | is_empty(parent_list)){
      #Then it means that the parent_list has at least 1 gene that is not present in the sum of all children's genes
      #In this case we have to output that the parent IS NOT contained in the children
      F
   }else if(length(o) == 0){
      #In this case we make sure the difference is 0 genes but the lists are not empty 
      T
   }
}

#Now we can just follow the tree down and store any pathway that belongs to that category:
#Create NOT %in% operator
'%ni%' <- Negate('%in%')

#Special function to create result lists
create_list <- function(x){
   data.frame(parent_id = character(),parent_name = character(),parent_length = integer(),parent_genes = list(),
              children_length = integer(),children_information = list(),is_parent_contained_in_children = logical(),is_parent_correct_size = logical(),
              all_children_correct_size = logical(),eliminate_flag = logical())
}

#Stablish cutoff for minimal pathway length
minimal_length <- 10
maximal_length <- 500

############################## CREATE FUNCTION TO ELIMINATE ALL PATHWAYS DUPLICATED AFTER PRUNING

#1. Some pathways are just duplicated because they appear with the same pathway ID in different parts of a category branch.

#2. Some pathways just appear in the Reactome tree but they are not included in the GMT file. Hence, we need to eliminate them, they have 0 genes

rm_duplicates <- function(sublist){
   #x should be a table. We can apply this function inside lapply to multiple data.frames (reactome_all,reactome_stored,reactome_eliminated...)
   #1. First, we find duplicated rows, which correspond to pathways entirely replicated in the category tree (if id, name, length and category are the same, we rule it out)
   temp <- sublist
   simple_duplicates <- temp %>%
      group_by(parent_id,parent_name,parent_length,category) %>%
      filter(n() > 1)
   #Eliminate them from store
   temp <- temp %>%
      group_by(parent_id,parent_name,parent_length,category) %>%
      filter(!duplicated(parent_id,parent_name,parent_length,category)) %>%
      ungroup()
   #2. Second, we find and discard replicated pathways with 0 counts (they were not in the gmt file, so their length is 0 genes)
   zero_duplicates <- temp %>%
      filter(parent_length == 0)
   temp <- temp %>%
      filter(parent_length != 0)
   #3. Third, those pathways that have the same name but different ID are in Disease. We will not eliminate them manually, as we will 
   #rule them out when the entire Disease category goes away
   id_replicates <- temp %>%
      group_by(parent_name,category) %>%
      filter(n() > 1)
   list(simple_duplicates = simple_duplicates,zero_duplicates = zero_duplicates,id_replicates = id_replicates,keep = temp)
}

#Create a function to start from a highest order parent(s) and climb down the tree
top_bottom_pruning <- function(big_categories,big_category_names,cutoff){
   names(big_categories)<- big_category_names
   #Initialise 3 lists: eliminated pathways, stored pathways and all pathways
   store_list <- lapply(big_categories,create_list)
   all_list <- lapply(big_categories,create_list)
   eliminated_list <- lapply(big_categories,create_list)
   #For each of the parents we got (it can be first order sons of a big category)
   #IMPORTANT! We use for and not a vectorised function like lapply or map because we want to SAVE inside the while() loop the current eliminated data.
   #This is much more difficult to do inside lapply or map, that have their own environment! Instead, the for loop allows us to do this easily, by enlarging the lists we created above
   for(i in 1:length(big_categories)){
      #Get current Reactome category to prune down
      input <- big_categories[[i]]
      current_category_name <- big_category_names[[i]]
      #If there is still one child, it means we did not get to the bottom
      id_child = 1
      parent_names <- sapply(input,map_id_name)
      cat(paste0("Pruning category: ",current_category_name))
      cat(" \n")
      #While there is at least 1 child, we keep going down the tree
      while (length(id_child) > 0){
         #Get index of children for each parent separately
         index_parent <- lapply(input,function(x) which(mapping_tree$parent_id %in% x))
         parent_names <- lapply(input,map_id_name) 
         names(parent_names) <- unlist(input)
         #We get the gene sets for each parent separately
         parent_genes <- lapply(input,map_id_genes)
         names(parent_genes) <- parent_names
         #We get the id of the children for each parent separately
         id_child <- lapply(index_parent,function(x)mapping_tree$child_id[x])
         name_child <- lapply(id_child,map_id_name)
         #Compare genes of the parent to those of the children. We take the genes of all the children that have correct size
         genes_sons <- lapply(id_child,function(x) unique(map_id_genes_between_threshold(x,minimal_length,maximal_length)))
         #Check if any of the children has wrong size to have that additional information
         all_children_correct_size <- sapply(id_child,function(x) check_children_correct_size(x,minimal_length,maximal_length))
         #Obtain information about children: pathway_name and list of genes
         children_tab <- lapply(id_child,function(x) get_children_info(x,minimal_length,maximal_length))
         #If there is more than one parent, we compare each parent to its children
         if(length(parent_names) > 1){
            genes_sons_lengths <- sapply(genes_sons,function(x) unique(length(x)))
            parent_genes_length <- sapply(parent_genes,length)
         }else if(length(parent_names) == 1){
            #If there is only 1 parent, we collapse all genes from children
            genes_sons_lengths <- length(unique(unlist(genes_sons,length)))
            parent_genes_length <- length(unlist(parent_genes))
         }
         names(parent_genes_length) <- parent_names
         #Check if we need to eliminate something
         is_parent_contained_in_children <- map2_lgl(.x = parent_genes,.y = genes_sons,.f = check_parent_is_contained)
         #If the pathway is very small or very big, we flag it for elimination
         is_parent_correct_size <- parent_genes_length >= minimal_length & parent_genes_length <= maximal_length
         ############ Eliminate parent if one of the following conditions is met
         eliminate_flag <- #1. Parent contained in children and correct size
            (is_parent_contained_in_children  & is_parent_correct_size)|
            #2. Parent contained in children and incorrect size 
            (is_parent_contained_in_children & !is_parent_correct_size)|
            #3. Parent not contained in children and incorrect size
            (!is_parent_contained_in_children & !is_parent_correct_size)
         ########### Keep parent if the following condition is met
         keep_flag <- #4. Parent not contained in children and correct size
            (!is_parent_contained_in_children & is_parent_correct_size)
         #Get index to eliminate
         index_eliminate <- which(eliminate_flag)
         index_keep <- which(keep_flag)
         #Create table with all the information about the parent and its children and the result of the logical rules and elimination decision
         tab <- tibble(parent_name = unlist(parent_names),parent_id = unlist(input),parent_length = parent_genes_length,
                       parent_genes = map(.x = parent_genes,.f = ~.x),
                       children_length = genes_sons_lengths, children_information = map(.x = children_tab,.f = ~.x),
                       is_parent_contained_in_children = is_parent_contained_in_children,
                       is_parent_correct_size = is_parent_correct_size,
                       all_children_correct_size = all_children_correct_size,
                       eliminate_flag = eliminate_flag) 
         #Table of pathways to eliminate
         tab_eliminate <- tab[index_eliminate, ]
         #Table of pathways to keep
         tab_keep <- tab[index_keep, ]
         #Add the parent to eliminated_list (if there is nothing to add, index eliminate will be integer(0) and nothing will be added)
         eliminated_list[[current_category_name]] <- rbind(eliminated_list[[current_category_name]],tab_eliminate)
         #Add the kept parents to stored_list
         store_list[[current_category_name]] <- rbind(store_list[[current_category_name]],tab_keep)
         #If there is at least 1 child, we keep climbing down the tree
         if(length(id_child) > 0){
            #Before adding the id_child to all_list, we eliminate those children that were 0 (their parent was a leaf)
            if(is.list(id_child)){
               #If it is a list, it means that in sapply up there, some parents were a leaf, so the list could not be collapsed to character array
               id_child <- unlist(id_child)
            }else{
               id_child <- id_child[length(id_child) != 0]
            }
            #Their lengths will be added in the next iteration of the algorithm, when they become parents
            #The new input is the children we just stored. We eliminate children that are a leaf and that are not useful as new inputs
            input <- id_child
         }
         #We store tab as all_list with both eliminated and stored results
         all_list[[current_category_name]] <-rbind(all_list[[current_category_name]],tab)
      }
   }
   
   #Put all 3 result lists inside an output list
   out <- list(store_list = store_list,all_list = all_list,eliminated_list = eliminated_list)
   
   ########## Add category names to result lits and turn categories into factors
   #Create function to add highest hierarchy names
   reactome_original <- map(.x = out,.f = function(sublist){
      imap(.x = sublist,.f = ~.x %>% mutate(category = factor(.y))) %>%
         bind_rows()
   }) 
   
   #Use categories in all_list to add all possible levels to eliminated and store list. This is needed because in eliminated_list
   #there are some categories with 0 pathways eliminated, like chromatin remodelling. If we do not have all levels, when collapsing the list
   #into a data.frame, chromatin remodelling will disappear and the 3 lists will have different length if splitted again by category
   reactome_original[c("eliminated_list","store_list")] <- map(.x = reactome_original[c("eliminated_list","store_list")],.f = ~.x %>%
                                                                  mutate(category = factor(category,levels = unique(reactome_original$all_list$category))))
   
   ########### Eliminate replicated pathways
   out_reps <- map(.x = reactome_original,.f = rm_duplicates)
}

################################### FUNCTIONS FOR EXPLORATORY DATA ANALYSIS

subtract <- function(x,y){x - y}
divide <- function(x,y){round(x/y*100,2)}

#Compute function to compare all and store list with GENES or PATHWAYS to get number and percentages of GENES/PATHWAYS lost per category
compare_all_and_store_per_category<- function(input_list,input_name){
   #Use map2 to compare categories of both results lists
   my_cols <- c(paste0(c("num_","total_",""),input_name,c("_lost","_in_category","_lost")))
   map2(.x = input_list$all_list,.y = input_list$store_list,.f = function(all,store){
      diff <- setdiff(all,store)
      #Get length of the difference
      diff_length <- length(diff)
      #Get length of total number of genes in that category
      total_genes_in_category <- length(all)
      tibble(diff_len = diff_length,total = total_genes_in_category) %>%
         mutate(diff_list = list(diff)) %>%
         setNames(my_cols) 
   }) %>%
      imap(.f = ~.x %>% mutate(category = .y)) %>%
      bind_rows() %>%
      mutate(percentage_lost = divide(!!as.name(my_cols[1]),!!as.name(my_cols[2]))) %>%
      select(category,contains("num"),contains("total"),everything()) 
}

#Compare all and store in total
compare_all_and_store_total <- function(input_list,input_name){
   my_cols <- paste0(c("num_","total_",""),input_name,c("_lost","_across_all_categories","_lost"))
   input_list <- map(.x = input_list,.f = function(sublist){
      unique_in_list <- unique(unlist(sublist))
   })
   #Compute diff
   diff <- setdiff(input_list$all_list,input_list$store_list)
   diff_length <- length(setdiff(input_list$all_list,input_list$store_list))
   #Total in list
   tibble(diff_len = diff_length,total = length(input_list$all_list)) %>%
      mutate(diff_list = list(diff)) %>%
      setNames(my_cols) %>%
      mutate(percentage_lost = divide(!!as.name(my_cols[1]),!!as.name(my_cols[2]))) 
}

#Function to compute number of genes lost during pruning (inside categories and between them)
#This function needs a data structure with 3 sublists (3 results), each of them splitted by category.
#This way it is much easier to compare lists between them by category and not collapsed data.frames
check_genes_pathways_lost <- function(temp){
   #Retrieve all genes from each category for every sublist
   list_genes <- map_depth(.x = temp,.depth = 2,.f = function(x){
      unique(map_id_genes(unique(x$parent_id)))
   })
   #Retrieve all pathways from each category for every sublist
   list_pathways <- map_depth(.x = temp,.depth = 2,.f = function(x){
      unique(x$parent_id)
   })
   #Get comparison of store and all for GENES
   comparison_genes_per_categ <- compare_all_and_store_per_category(list_genes,"genes")
   comparison_genes_per_categ_sum <- comparison_genes_per_categ %>% summarise(mean_genes_lost = mean(num_genes_lost),
                                                                              sd_genes_lost = sd(num_genes_lost),
                                                                              mean_perc_lost = mean(percentage_lost),
                                                                              sd_perc_lost = sd(percentage_lost))
   #Get comparison of store and all for PATHWAYS
   comparison_pathways_per_categ <- compare_all_and_store_per_category(list_pathways,"pathways")
   comparison_pathways_per_categ_sum <- comparison_pathways_per_categ %>% summarise(mean_pathways_lost = mean(num_pathways_lost),
                                                                                    sd_pathways_lost = sd(num_pathways_lost),
                                                                                    mean_perc_lost = mean(percentage_lost),
                                                                                    sd_perc_lost = sd(percentage_lost))
   #Among all categories
   #GENES
   comparison_genes_tot <- compare_all_and_store_total(list_genes,"genes")
   #Pathways
   comparison_pathways_tot <- compare_all_and_store_total(list_pathways,"pathways")
   #IMPORTANT! Compare genes_lost in total with genes_lost per category to see which category is reponsible for the loss
   genes_lost_tot <- comparison_genes_tot$genes_lost[[1]]
   genes_lost_categ <- comparison_genes_per_categ %>% select(category,genes_lost) %>% unnest(genes_lost) %>%
      #Get those genes in genes_lost_tot
      filter(genes_lost %in% genes_lost_tot) %>%
      #Count per category
      group_by(category) %>%
      summarise(num_genes = n()) %>%
      ungroup() %>%
      arrange(desc(num_genes))
   #Output
   list(genes_lost_per_category = comparison_genes_per_categ,
        summary_genes_lost_per_category = comparison_genes_per_categ_sum,
        genes_lost_total = comparison_genes_tot,
        list_genes_lost_per_category = genes_lost_categ,
        pathways_lost_per_category = comparison_pathways_per_categ,
        summary_pathways_lost_per_category = comparison_pathways_per_categ_sum,
        pathways_lost_total = comparison_pathways_tot)
   #list(number_genes_lost_per_category = difference_between_length,percentage_genes_lost_per_category = percentage_between,number_genes_lost_total = difference_all,percentage_genes_lost_total = percentage_all)
}

#Function to check total number of pathways in sublists
sublist_num <- function(x,sublist){
   list_all_pathways <- sapply(x[[sublist]],function(y) y$parent_id)
   #Get the number of unique pathways across all categories of a certain list (all,eliminated or stored)
   length(unique(unlist(list_all_pathways)))
}

#Exploratory data function
exploratory_function <- function(results){
   #Split each sublist into categories again (easier to compare and manipulate)
   temp <- map(.x = results,.f = ~.x %>% split(f = .$category,drop = F))
   #Check that for every category, the sum of store + eliminated = all
   check_sum_lists <- pmap_lgl(.l = temp,.f = function(store_list,all_list,eliminated_list){
      ifelse(nrow(store_list) + nrow(eliminated_list) == nrow(all_list),T,F)
   })
   if(all(check_sum_lists)){
      cat("\n\nCORRECT! The sum of eliminated + stored = all pathways for every category!\n\n")
   }else{
      cat("PROBLEM! The sum of eliminated + stored does NOT equal all pathways!")
   }
   #Check if all pathways are in all_list (no duplicates)
   cat(paste0("The total number of pathways before pruning is ",sublist_num(temp,"all_list")))
   cat("\n\n")
   #Check number of pathways in eliminated_list
   cat(paste0("The total number of pruned (eliminated) pathways is ",sublist_num(temp,"eliminated_list")))
   cat("\n\n")
   #Check unique remaining pathways
   cat(paste0("The total number of remaining pathways is ",sublist_num(temp,"store_list")))
   cat("\n\n")
   if(sum(sublist_num(temp,"store_list") + sublist_num(temp,"eliminated_list") == sublist_num(temp,"all_list"))){
      cat("CORRECT! The total number of stored + eliminated = all pathways!\n\n")
   }else{
      cat("PROBLEM! The total number of stored + eliminated does NOT equal all pathways!\n\n")
   }
   #Check number of genes lost in pruning (we map unique pathways in all_pathways to their genes and compare to mapped genes in store_list)
   print(check_genes_pathways_lost(temp))
}

################################### PRUNING FOR ALL CATEGORIES 

#Apply the function to each category independently
res_list <- top_bottom_pruning(highest_order_parents,highest_order_parent_names,cutoff = 0)

#Save excell file with eliminated pathways per category
eliminated_pathways <- res_list$eliminated_list$keep %>%
   select(parent_name,parent_id,parent_length,category,parent_genes,children_information,children_length,is_parent_contained_in_children,is_parent_correct_size,all_children_correct_size) %>%
   #Select only pathway_names from children_information
   mutate(children_information = map(.x = children_information,.f = function(x){
      if(is.character(x)){
         ""
      }else{
         paste0(x$child_name,collapse = "; ")
      }
   })) %>%
   unnest(children_information) %>%
   #Collapse pathway_genes into a single value with commas
   unnest(parent_genes) %>%
   group_by(parent_name) %>%
   mutate(parent_genes = paste0(parent_genes,collapse = ", ")) %>%
   ungroup() %>%
   distinct(parent_name,parent_genes,.keep_all = T) %>%
   #Replace parent with pathway for final tabular version
   rename_with(.fn = ~gsub("parent","pathway",.x),.cols = contains("parent")) %>%
   #Replace long category names that cannot be excel sheet names
   mutate(category = case_when(
      category == "Extracellular matrix organization" ~ "ECM organization",
      category == "Organelle biogenesis and maintenance" ~ "Organelle biogenesis",
      TRUE ~ category
   )) %>%
   split(f = .$category)

## Create new workbooks for transcriptomics and another for proteomics
wb <- createWorkbook() 
imap(.x = eliminated_pathways,.f = function(x,category){
   x = as.data.frame(x)
   ## Create the worksheets
   addWorksheet(wb, sheetName = category)
   ## Write the data
   writeData(wb,category,x)
})

#Once we have filled in all sheets of the workbook, we save it
saveWorkbook(wb, file = paste0(pdir,"/Outputs/Top_bottom_pruning_eliminated_pathways.xlsx"), overwrite = TRUE)

#Obtain the duplicates on one side and the results without duplicates on another structure
res_duplicates <- map(.x = res_list,.f = ~.x[!grepl("keep",names(.x))])

#Results without duplicates
reactome_pruned <- map(.x = res_list,.f = ~.x$keep)

reactome_pruned_save <- reactome_pruned$store_list %>%
   #Select only the important rows we need 
   select(parent_name,parent_id,parent_length,category) %>%
   setNames(c("pathway_name","pathway_id","length","category")) %>%
   left_join(table_with_genes[,c("pathway_id","gene_name","category")],by = c("pathway_id","category"))

########################## COUNT NUMBER OF PATHWAYS AND GENES LOST IN PRUNING DIRECTLY

#Number of pathways lost
reactome_pruned_pathways_lost <- reactome_pruned$eliminated_list %>% distinct(parent_name) %>% nrow()
#Number of pathways kept
reactome_pruned_pathways_kept <- reactome_pruned$store_list %>% distinct(parent_name) %>% nrow()
#Exactly 1214 are discarded with top-down pruning
#Number of genes
reactome_pruned_genes_kept <- reactome_pruned$store_list %>% unnest(parent_genes) %>% distinct(parent_genes)
genes_in_eliminated_top_down <- reactome_pruned$eliminated_list %>% unnest(parent_genes) %>% distinct(parent_genes)
reactome_pruned_genes_lost <-  genes_in_eliminated_top_down %>% anti_join(reactome_pruned_genes_kept)
#Exactly 361 genes are lost with top-down pruning
#Check number of pathways per category according to the algorithm
computed_num_pathways_per_category <- reactome_pruned$all_list %>% group_by(category) %>% summarise(algorithm_number_pathways = n())
computed_num_pathways_per_category <- computed_num_pathways_per_category[order(names(computed_num_pathways_per_category))]
#Add to num_pathways_per_category
num_pathways_per_category <- num_pathways_per_category %>% left_join(computed_num_pathways_per_category,by = "category") %>% as.data.frame()
num_pathways_per_category

################################## EXPLORATORY ANALYSIS TO CHECK NUMBER OF GENES AND PATHWAYS PRUNED

#Exploratory analysis
exploratory_data <- exploratory_function(reactome_pruned)

################################## GET HOW MANY GENES ARE LOST UNIQUELY TO ONE CATEGORY

#IMPORTANT! There are only 361 genes lost ACROSS ALL CATEGORIES, but if we sum the genes lost PER CATEGORY
#we get 446. This means that there are genes lost in one category but that are kept in another, so they are not really lost
exploratory_genes <- exploratory_data$genes_lost_per_category %>%
   distinct(category,genes_lost) %>%
   unnest(genes_lost) %>%
   #Compute number of categories per gene
   group_by(genes_lost) %>%
   mutate(num_categ = n()) %>%
   arrange(desc(num_categ))

################################## PLOT THE NUMBER OF LOST GENES/PATHWAYS LOST PER CATEGORY

exploratory_genes <- exploratory_data$genes_lost_per_category %>% mutate(data = "Genes")
#IMPORTANT! The problem with exploratory_genes is that there are repeated genes, because some are saved by another category
names(exploratory_genes) <- str_replace(names(exploratory_genes),"genes_","") 
exploratory_genes <- exploratory_genes %>%
   select(-num_lost,-lost,-percentage_lost) %>%
   left_join(exploratory_data$list_genes_lost_per_category,by = "category") %>%
   rename(num_lost = num_genes) %>%
   mutate(num_lost = ifelse(is.na(num_lost),0,num_lost),
          percentage_lost = round(num_lost/total_in_category*100,2))
exploratory_path <- exploratory_data$pathways_lost_per_category %>% mutate(data = "Pathways")
names(exploratory_path) <- str_replace(names(exploratory_path),"pathways_","")
exploratory_path <- exploratory_path %>%
   select(-lost)
#Put both datasets together
merged_things_lost_per_categ <- rbind(exploratory_genes,exploratory_path) %>%
   mutate(num_remaining = total_in_category - num_lost) 

#Save in excel file to put in a plot
merged_things_lost_per_categ_excel <- merged_things_lost_per_categ %>%
   pivot_wider(names_from = "data",values_from = c(num_lost,total_in_category,percentage_lost,num_remaining)) %>%
   select(category,num_lost_Genes,percentage_lost_Genes,total_in_category_Genes,num_lost_Pathways,percentage_lost_Pathways,total_in_category_Pathways) %>%
   rename("Lost genes" = num_lost_Genes,"Lost genes (%)" = percentage_lost_Genes,"Total genes" = total_in_category_Genes,
          "Lost pathways" = num_lost_Pathways,"Lost pathways (%)" = percentage_lost_Pathways,"Total pathways" = total_in_category_Pathways) %>%
   arrange(category)

write.xlsx(merged_things_lost_per_categ_excel,paste0(pdir,"/Outputs/Top_down_pruning_summary_table.xlsx"))

#Do a barplot showing categories and the number of genes or pathways (facetted) in the whole and pruned datasets
colors <- c( "Remaining" = "coral1","1.Top-down\npruning" = "dodgerblue")
plot_lost_genes_path_top_down <- ggplot(merged_things_lost_per_categ,aes(x = category,y = total_in_category)) + 
   geom_col(aes(fill = "1.Top-down\npruning"),position = "identity") + 
   geom_col(aes(y = num_remaining,fill = "Remaining"),position = "identity") + 
   facet_wrap(scales = "free_x",nrow = 1,ncol = 2,facets = "data") + 
   theme_bw() +
   coord_flip() +
   scale_x_discrete(limits = rev) + 
   labs(x = "",y = "Number of elements in Reactome") + 
   scale_fill_manual(name = "",breaks = c("Remaining","1.Top-down\npruning"),values = colors) + 
   theme(axis.text.x = element_text(size = 15),
         axis.text.y = element_text(size = 15),
         axis.title.x = element_text(size = 16),
         strip.text = element_text(size = 15),
         legend.text = element_text(size = 15),
         legend.position = "top")

#Save in pdf format
# pdf(file = "Plots/Genes_pathways_lost_top_down_Reactome.pdf",width = 8,height = 7)
# print(plot_lost_genes_path_top_down)
# dev.off()

############################################# PLOT THE NUMBER OF GENES/PATHWAYS LOST IN TOTAL

genes_lost_total <- exploratory_data$genes_lost_total %>% mutate(data = "Genes")
names(genes_lost_total) <- str_replace(names(genes_lost_total),"genes_","")
pathways_lost_total <- exploratory_data$pathways_lost_total %>% mutate(data = "Pathways")
names(pathways_lost_total) <- str_replace(names(pathways_lost_total),"pathways_","")
#Put both datasets together
merged_things_lost <- rbind(genes_lost_total,pathways_lost_total) %>%
   mutate(num_remaining = total_across_all_categories - num_lost) %>%
   select(num_lost,num_remaining,data) %>%
   pivot_longer(cols = -c("data"),names_to = "metric",values_to = "value") %>%
   mutate(metric = case_when(
      metric == "num_lost" ~ "1.Top-down\npruning",
      metric == "num_remaining" ~ "Remaining"
   ),
   metric = factor(metric,levels = c("Remaining","1.Top-down\npruning"))) %>%
   group_by(data) %>%
   mutate(perc = paste0(value," (",round(value/sum(value)*100,2),"%)")) %>%
   ungroup() %>%
   mutate()
#Do a barplot showing categories and the number of genes or pathways (facetted) in the whole and pruned datasets
my_plot <- ggplot(merged_things_lost,aes(x = "",y = value,fill = metric)) + 
   geom_col() +
   geom_text(aes(label = perc),position = position_stack(vjust = 0.5),size = 5) + 
   facet_wrap(scales = "free_y",nrow = 2,ncol = 1,facets = "data") + 
   theme_bw() + 
   labs(x = "",y = "Number of elements in Reactome",fill = "",title = "Top-down pruning") + 
   scale_fill_manual(values = c("coral1","dodgerblue")) + 
   theme(axis.text.x = element_text(angle = 90,size = 13),
         axis.text.y = element_text(size = 13),
         axis.title.y = element_text(size = 15),
         axis.ticks.x = element_blank(),
         strip.text = element_text(size = 13),
         legend.position = "top",
         legend.text = element_text(size = 13),
         plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

#Plot
plot(my_plot)
dev.off()

############################# GET LIST OF PATHWAYS WE WANT TO KEEP IN OUR DATA FROM THIS LIST

############# 1. There are two pathways in unwanted categories that we want to keep in the gmt file

pathways_to_keep_manual <- c("POU5F1 (OCT4), SOX2, NANOG activate genes related to proliferation",
                             "POU5F1 (OCT4), SOX2, NANOG repress genes related to differentiation")

#Check categories to which this pathways belong
tab_pathways_categories <- table_with_genes %>% distinct(pathway_name,category)
pathways_to_keep_manual_tab <- data.frame(pathway_name = pathways_to_keep_manual) %>%
   left_join(tab_pathways_categories,by = "pathway_name")
pathways_to_keep_manual_tab
pathways_to_keep <- pathways_to_keep_manual

############################# ELIMINATE ALL PATHWAYS IN UNWANTED CATEGORIES EXCEPT THOSE WE WANT TO KEEP

categories_to_eliminate <- c("Hemostasis","Digestion and absorption","Reproduction","Sensory Perception","Neuronal System","Developmental Biology","Muscle contraction","Disease")
#Get pathways in those categories that are not in pathways_to_keep
pathways_from_categories_to_eliminate <- reactome_pruned$store_list %>%
   filter(category %in% categories_to_eliminate) %>%
   filter(parent_name %ni% pathways_to_keep)

#Eliminate them
reactome_pruned_double_filtered <- reactome_pruned$store_list %>%
   anti_join(pathways_from_categories_to_eliminate,by = c("parent_name","category")) 
reactome_pruned_double_filtered

#Count remaining number of pathways
reactome_pruned_double_filtered %>% distinct(parent_name) %>% nrow()

################################ ELIMINATE PATHWAYS MANUALLY 

#There are a few pathways related to the metabolism of Abacavir and aflatoxins that we want to eliminate:
#Load reactome table with all pathways to get the full names
reactome_table <- readRDS(paste0(pdir,"/Outputs/Reactome_table_all_categories.rds"))
manual_to_eliminate <- c("Abacavir","Aflatoxin")
#Get pathway names
manual_pathways_to_eliminate <- reactome_table %>%
   filter(grepl(paste0(manual_to_eliminate,collapse = "|"),x = pathway_name))
#Eliminate them from reactome_pruned_doubled_filtered
reactome_pruned_double_filtered <- reactome_pruned_double_filtered %>%
   filter(!grepl(paste0(manual_to_eliminate,collapse = "|"),x = parent_name))

################################ STORE REACTOME TABLE PRUNED WITH UNWANTED CATEGORIES FILTERED

saveRDS(reactome_pruned_double_filtered,file = paste0(pdir,"/Outputs/Reactome_table_topdown_pruned_filtered.rds"))

############################## ADD THE GENES THAT DEFINE EACH PATHWAY

reactome_pruned_double_filtered_genes <- reactome_pruned_double_filtered %>%
   #Select only the important rows we need 
   select(parent_name,parent_id,parent_length,category) %>%
   setNames(c("pathway_name","pathway_id","length","category")) %>%
   left_join(table_with_genes[,c("pathway_id","gene_name","category")],by = c("pathway_id","category"))

saveRDS(reactome_pruned_double_filtered_genes,file = paste0(pdir,"/Outputs/Reactome_table_topdown_pruned_filtered_genes.rds"))

############################## STORE REACTOME TABLE PRUNED WITH ALL CATEGORIES

saveRDS(reactome_pruned$store_list,file = paste0(pdir,"/Outputs/Reactome_table_topdown_pruned_all_categories.rds"))

