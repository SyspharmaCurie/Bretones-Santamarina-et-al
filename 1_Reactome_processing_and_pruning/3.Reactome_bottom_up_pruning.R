############################# SCRIPT TO DO BOTTOM UP PRUNING FOR THE REACTOME DATABASE
#Install packages, if required

if (!require("tidyverse", character.only = TRUE)) {
  install.packages("tidyverse")
}
#Load packages
library(tidyverse)

#Set current directory to script directory in R Studio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Get current dir
pdir <- getwd()

#We will compute, for each pathway, 3 metrics and we will store them in a list for each comparison with other pathway:

#The intersecting genes with the rest of the pathways.
#The percentage of shared genes out of the total in the pathway
#Forward or reverse depending on whether pathway A is contained in B or the other way around.


#Load Reactome table after top-down pruning and unwanted category elimination
reactome_table <- readRDS(paste0(pdir,"/Outputs/Reactome_table_topdown_pruned_filtered_genes.rds"))

#Get the unique genes, we no longer care about pathways repeated in different categories
reactome_table_temp <- reactome_table %>% distinct(pathway_name,length,gene_name)

results <- reactome_table_temp %>%
   group_by(pathway_name,.drop = F) %>%
   #Create for each pathway a list with all genes inside
   summarize(genes = list(gene_name)) %>%
   rename(pathway = pathway_name) %>%
   #Use expand to add new columns to the data, expressing the pairwise operation between pathways
   expand(nesting(pathway,genes),nesting(pathway2 =pathway,genes2 = genes)) %>%
   #Eliminate rows corresponding to comparisons of each pathway with itself (we are not interested in this)
   filter(pathway != pathway2) %>%
   #Add pathway length
   mutate(length = map_int(genes,~length(.x)),
          length2 = map_int(genes2,~length(.x))) %>%
   #We need to eliminate replicated combinations. The nesting function in tidyr gives all possible combinations of pathways, so we have
   #a vs b and b vs a. We need to eliminate half of the entries, to get only the information we are looking for
   #To do that, we will only keep comparisons where pathway is before alphabetically compared to pathway2 (eliminate half of the results)
   filter(pathway > pathway2) %>%
   #Add new column computing the intersection of the genes
   mutate(intersection = map2(.x = genes,.y = genes2,.f = function(x,y){
      intersect(x,y)
   })) %>%
   #And the number of intersecting genes
   mutate(number_intersecting_genes = map_int(.x = intersection,.f = ~length(.x))) %>%
   #Add new column computing the percentage of shared genes
   mutate(percentage_shared = pmap_dbl(.l = list(number_intersecting_genes,length,length2),.f = function(inter,len1,len2){
      case_when(
         #If pathway A is <= B, then we compute it with respect to pathway A
         len1 <= len2 ~ inter/len1*100,
         #If pathway A is > B, then we compute it with respect to pathway B
         len1 > len2 ~ inter/len2*100
      )
   })) %>%
   mutate(longer_1 = ifelse(length > length2,1,0))

#Define list of unique pathways in our data
unique_pathways <- unique(reactome_table_temp$pathway_name)
`%ni%` = Negate(`%in%`)
#Initialise table to store necessary pathways (those that are not contained in others)
necessary_pathways <- data.frame()
#Initialise table to store unnecessary pathways (those contained in others)
unnecessary_pathways <- data.frame()
#Create temporal data, a table with all pairwise comparisons that we will start to prune
temp <- results
#Define threshold (from 1 to 100%). 100% means that we do not accept a single gene discarded when eliminating a pathway.
#If threshold < 100 we are saying that if a % of genes or more are included elsewhere, we rule out the pathway
threshold <- 100

######### Eliminate those pathway interactions where % shared < threshold

temp <- filter(temp,percentage_shared >= threshold)
#We just go through each pathway and see if it is not contained in other pathways at the threshold. If it is not, we keep it
for(i in 1:length(unique_pathways)){
   cat(paste0("\nProcessing pathway ",i,"\n"))
   #Get the current pathway under study
   current_pathway <- unique_pathways[i]
   #Check if the pathway is still in the temp data or was already eliminated
   if(current_pathway %ni% temp$pathway & current_pathway %ni% temp$pathway2){
      #If it is NOT in temp data, this pathway was already elimnated and we can go to the next
      next
   }else{
      #If the pathway is still in temp, we retrieve all the interactions for this pathway
      dat <- temp %>%
         filter(pathway == current_pathway | pathway2 == current_pathway) 
      #We check if the pathway contains other pathways
      contains_index <- as.logical(ifelse(dat$percentage_shared >= threshold & dat$pathway == current_pathway & dat$longer_1 == 1 |
                                             dat$percentage_shared >= threshold & dat$pathway2 == current_pathway & dat$longer_1 == 0,1,0))
      #If it contains at least 1 pathway, contains = T
      contains <- ifelse(any(contains_index),1,0)
      #If no 1 is present, this means the pathway does not contain any other pathway. Then we go to the next one
      #If it contains at least 1 pathway, we eliminate those that are contained and go to the next one
      if(contains == 1){
         #If the pathway contains other pathways, we eliminate those pathways it contains
         pathways_to_discard <- dat %>%
            filter(as.logical(contains_index)) 
         pathways_to_discard_1 <- pathways_to_discard %>% filter(pathway != current_pathway) %>% rename(eliminated_pathway = pathway,bigger_pathway = pathway2,
                                                                                                        eliminated_genes = genes,bigger_genes = genes2,eliminated_length = length,
                                                                                                        bigger_length = length2) 
         pathways_to_discard_2 <- pathways_to_discard %>% filter(pathway2 != current_pathway) %>% rename(eliminated_pathway = pathway2,bigger_pathway = pathway,
                                                                                                         eliminated_genes = genes2,bigger_genes = genes,eliminated_length = length2,
                                                                                                         bigger_length = length)
         #Put them together by joining both tables
         pathways_to_discard_final <- rbind(pathways_to_discard_1,pathways_to_discard_2)
         unnecessary_pathways <- rbind(unnecessary_pathways,pathways_to_discard_final)
         #Eliminate them from temp
         temp <- temp %>%
            filter(pathway %ni% pathways_to_discard) %>%
            filter(pathway2 %ni% pathways_to_discard_2)
      }
   }
}

eliminated_pathways <- unique(unnecessary_pathways$eliminated_pathway)
necessary_pathways <- setdiff(unique_pathways,eliminated_pathways)
final_necessary <- necessary_pathways

#Check that all bigger_pathway from unncessary_pathways are found in final_necessary or are also eliminated because they were also contained elsewhere
unnecessary_pathways <- unnecessary_pathways %>%
   mutate(is_bigger_final_pathway = map_lgl(.x = bigger_pathway,.f = ~.x %in% final_necessary & .x %ni% eliminated_pathway)) %>%
   select(-longer_1,-intersection,-eliminated_genes,-bigger_genes) %>%
   #Get only those rows corresponding to big pathways that are final pathways
   filter(is_bigger_final_pathway)

#Save in excell format
library(openxlsx)
write.xlsx(unnecessary_pathways,paste0(pdir,"/Outputs/Bottom_up_pruning_eliminated_pathways.xlsx"),overwrite = T)

#Get all information about necessary pathways
necessary_pathways_genes <- reactome_table_temp %>%
   filter(pathway_name %in% final_necessary) 

#Get number of necessary genes (unique genes in necessary pathways)
necessary_num_genes <- necessary_pathways_genes %>%
   distinct(gene_name) %>%
   pull() %>%
   length()
#Get number ot total genes (unique genes in the original data of this script)
total_genes_reactome <- reactome_table_temp %>%
   distinct(gene_name) %>%
   pull() %>%
   length()

cbind(necessary_genes = necessary_num_genes,total_genes_reactome = total_genes_reactome)
num_lost_genes <- total_genes_reactome - necessary_num_genes
num_lost_pathways <- length(eliminated_pathways)

############## PLOT NUMBER OF GENES/PATHWAYS LOST DURING BOTTOM-UP PRUNING

genes_lost_total <- tibble(num_lost = num_lost_genes,total_genes = total_genes_reactome,data = "Genes")
pathways_lost_total <- tibble(num_lost = num_lost_pathways,total_pathways = length(unique_pathways),data = "Pathways")
names(genes_lost_total) <- str_replace(names(genes_lost_total),"_genes","")
names(pathways_lost_total) <- str_replace(names(pathways_lost_total),"_pathways","")
#Put both datasets together
merged_things_lost <- rbind(genes_lost_total,pathways_lost_total) %>%
   mutate(num_remaining = total - num_lost) %>%
   select(num_lost,num_remaining,data) %>%
   pivot_longer(cols = -c("data"),names_to = "metric",values_to = "value") %>%
   mutate(metric = case_when(
      metric == "num_lost" ~ "Discarded during\nbottom-up pruning",
      metric == "num_remaining" ~ "Remaining after\nbottom-up pruning"
   ),
   metric = factor(metric,levels = c("Discarded during\nbottom-up pruning","Remaining after\nbottom-up pruning"))) %>%
   group_by(data) %>%
   mutate(perc = paste0(value," (",round(value/sum(value)*100,2),"%)")) %>%
   ungroup() %>%
   mutate()
#Do a barplot showing categories and the number of genes or pathways (facetted) in the whole and pruned datasets
ggplot(merged_things_lost,aes(x = "",y = value,fill = metric)) + 
   geom_col() +
   geom_text(aes(label = perc),position = position_stack(vjust = 0.5),size = 5) + 
   facet_wrap(scales = "free_y",nrow = 2,ncol = 1,facets = "data") + 
   theme_bw() + 
   labs(x = "",y = "Number of elements in Reactome top-down pruned",
        fill = "",title = "Bottom-up pruning") + 
   scale_fill_manual(values = c("coral1","dodgerblue")) + 
   theme(axis.text.x = element_text(angle = 90,size = 13),
         axis.text.y = element_text(size = 13),
         axis.title.y = element_text(size = 15),
         axis.ticks.x = element_blank(),
         strip.text = element_text(size = 13),
         legend.position = "top",
         legend.text = element_text(size = 13),
         plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

########## ELIMINATE THE UNNECESSARY PATHWAYS FROM REACTOME TABLE

#Filter the necessary pathways from reactome_table
reactome_table_pruned <- reactome_table %>%
   filter(pathway_name %in% final_necessary)

#Final number of pathways
unique_pathways <- unique(reactome_table_pruned$pathway_name)
unique_pathways

######## Save Reactome table after pruning with genes
saveRDS(reactome_table_pruned,paste0(pdir,"/Outputs/Reactome_table_twice_pruned_genes.rds"))
   
####### Save filtered table without genes
reactome_table_pruned_no_genes <- reactome_table_pruned %>%
   distinct(pathway_name,pathway_id,length,category)

saveRDS(reactome_table_pruned_no_genes,paste0(pdir,"/Outputs/Reactome_table_twice_pruned.rds"))

############################################ CREATE FINAL GMT FILE AFTER PRUNING 

#Upload the Reactome table with all categories to know the category to which pathways belong in the gmt file
reactome_table_original <- readRDS(paste0(pdir,"/Outputs/Reactome_table_all_categories.rds"))
reactome_table_original_genes <- readRDS(paste0(pdir,"/Outputs/Reactome_table_all_categories_with_genes.rds"))

#We upload the gmt file with the full Reactome database after processing
library(GSA)
temp_gmt <- GSA.read.gmt(paste0(pdir,"/Outputs/reactome_gmt_symbol_without_pathways_not_in_tree.gmt"))

#Obtain descriptions and make the link with pathway names
descriptions_and_names <- data.frame(Pathway = temp_gmt$geneset.names,description = temp_gmt$geneset.descriptions,stringsAsFactors = F)
#Now we change the format in temp_gmt to a list of lists where each sublist is a pathway and inside of it we got the genes
reactome_gmt <- temp_gmt$genesets
names(reactome_gmt) <- temp_gmt$geneset.names

###################### CREATE FUNCTION TO ASSEMBLE AND STORE A GMT FILE

create_save_GMT_file <- function(loaded_gmt,file_name){
   #Now we just assemble back the gmt file
   loaded_gmt_store <- loaded_gmt$genesets
   #We add the descriptions
   loaded_gmt_store <- mapply(c,loaded_gmt$geneset.descriptions,loaded_gmt_store, SIMPLIFY=FALSE)
   #Add correct list names
   names(loaded_gmt_store) <- loaded_gmt$geneset.names
   #Store the final gmt file
   my_names <- sapply(names(loaded_gmt_store),function(x) paste(x,paste(loaded_gmt_store[[x]],collapse="\t"),sep = "\t"))
   write(file = file_name,x = my_names, ncolumns=10000,sep = "\n")
}

########################### ELIMINATE FROM THE GMT FILE ALL PATHWAYS NOT INCLUDED IN REACTOME_TABLE_PRUNED

#Get pathways that were eliminated from reactome_table_filtered. They include:

######## 1. PATHWAYS FROM UNWANTED CATEGORIES

####### 2. PATHWAYS ELIMINATED DURING THE TOP-DOWN PRUNING (EITHER <10 OR >500 GENES) OR CONTAINED IN THE SUM OF THE CHILDREN

###### 3. PATHWAYS ELIMINATED DURING THE BOTTOM-UP PRUNING, WHICH ARE PATHWAYS CONTAINED IN A BIGGER PATHWAY

pathways_eliminated_from_reactome <- reactome_table_original %>%
   anti_join(reactome_table_pruned,by = "pathway_name")

#We take those pathways and we eliminate them from the gmt file
reactome_gmt_no_unwanted <- lapply(temp_gmt,function(x) x[temp_gmt$geneset.names %ni% pathways_eliminated_from_reactome$pathway_name])

########################### STORE PRUNED GMT FILE

create_save_GMT_file(reactome_gmt_no_unwanted,paste0(pdir,"/Outputs/reactome_gmt_symbol_twice_pruned.gmt"))
