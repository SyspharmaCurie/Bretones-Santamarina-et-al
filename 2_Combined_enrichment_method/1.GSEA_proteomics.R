########################### GENE SET ENRICHMENT ANALYSIS (GSEA) WITH FGSEA PACKAGE FOR ALL MUTANTS
#install packages, if required
if (!require("tidyverse", character.only = TRUE)) {
  install.packages("tidyverse")
}
if (!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
}
if (!require("ggplot2", character.only = TRUE)) {
  install.packages("ggplot2")
}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("fgsea", character.only = TRUE))
BiocManager::install("fgsea")


#Load packages
library(tidyverse)
library(fgsea)
library(data.table)

#Set current directory to script directory in R Studio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load ranked protein lists from differential expression analysis with limma for all mutants
set_123 <- readRDS("Inputs/proteomics_DE_results_sets123.rds")
set_4 <- readRDS("Inputs/proteomics_DE_results_set4.rds") 

#Put into a single list
ranked_lists <- c(sets123 = set_123,set4 = set_4) %>% setNames(str_replace_all(names(.),"(.*\\.)","")) %>%
   map(.f = ~.x %>% rename(stat = t,padj = qval,log2FoldChange = logFC) %>% select(Gene = gene_symbol_human,log2FoldChange,stat,padj,msd,line))

############################### A) ENRICHMENT WITH PRUNED GMT FILE FROM REACTOME

#In fgsea we need to input GMT files including pathway definitions
path_reactome_gmt <- "Inputs/reactome_gmt_symbol_twice_pruned.gmt"
path_reactome_table <- "Inputs/Reactome_table_twice_pruned.rds"

#Put all paths into a list and load at once
reactome_gmt_file <- gmtPathways(path_reactome_gmt)

#Put all paths into a list and load at once
reactome_table <- readRDS(path_reactome_table)

#Choose ranking metric
ranking_metric <- "stat"

library(rlist)
#IMPORTANT!! We must choose a seed to ensure the analysis is reproducible. This is because GSEA uses a 
#permutation test to span the null distribution of the statistic. Hence, each time we run it, results are 
#expected to change slightly.
my_seed <- 179
set.seed(my_seed)

############################## RUN FGSEA USING PARALLEL COMPUTING

#Use future_map for improved performance through parallel computing
#We will run 10 million permutations in order to obtain small p-values (see fgsea documentation)

library(furrr)
plan(multisession, workers = 18)

gsea_results_final <- future_imap(.x = ranked_lists,function(x,mutant_name){
   #Get gene list from ranking metric
   gene_list <- x %>% select(all_of(ranking_metric)) %>% pull()
   #Name the vector
   names(gene_list) <- x$Gene
   #Sort the list in decreasing order
   gene_list = sort(gene_list, decreasing = F)
   #GSEA analysis for the current cell line
   cat(paste0("Processing cell line ",mutant_name,"\n\n"))
   #Run 10 million permutations in order to estimate p-values as precisely as possible
   database_results <- fgseaMultilevel(pathways = reactome_gmt_file,gene_list,nPermSimple = 1e7,eps = 0)
},.options = furrr_options(seed = my_seed))

#Save raw results
saveRDS(gsea_results_final,file = paste0(pdir,"Outputs/GSEA_results_proteomics.rds"))
gsea_results_final <- readRDS(paste0("Final_GSEA_proteomics_gmt_",gmt_version,"_dataset_",dataset,"_metric_",ranking_metric,".rds"))

#Add Reactome pathway categories collapsed with commas (untidy format)
gsea_results_collapsed <- lapply(gsea_results_final,function(x){
   x <- x %>%
      rename(pathway_name = pathway) %>%
      #Get pathway ids and categories from the reactome table
      left_join(select(table_list$reactome,pathway_name,pathway_id,category),by = "pathway_name") %>%
      group_by(pathway_name) %>%
      mutate(category = paste0(category,collapse = ", ")) %>%
      distinct()
})

#Add Reactome pathway categories in tidy format (one row per category)
gsea_results_categories <- lapply(gsea_results_final,function(x){
   x <- x %>%
      rename(pathway_name = pathway) %>%
      #Get pathway ids and categories from the reactome table
      left_join(select(table_list$reactome,pathway_name,pathway_id,category),by = "pathway_name") 
})

############ STORE GSEA RESULTS

#Collapsed 
saveRDS(gsea_results_collapsed,file = paste0("GSEA_proteomics_results.rds"))

