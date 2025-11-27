########################### SCRIPT OF SLI ALGORITHM

#####Dataset download ########
# 1. Download transcriptomics of selected DepMap cell lines at "https://drive.google.com/file/d/1nBJ9qc5dgELHOls8PNUdUvpUh1Dur3wx/view?usp=sharing"

# 2. Download CRISPR screen of selected DepMap cell lines, at "https://drive.google.com/file/d/1Pm-G2kiorDezAu6mY-sI1KpJEWl1ORkx/view?usp=sharing"
  
# 3. Place both files in the Inputs folder

########################################################################


#Set current directory to script directory in R Studio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Obtain current dir
pdir <- getwd()


#install packages, if required

if (!require("tidyverse", character.only = TRUE)) {
  install.packages("tidyverse")
}

if (!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
}
if (!require("GSA", character.only = TRUE)) {
  install.packages("GSA")}

if (!require("broom", character.only = TRUE)) {
  install.packages("broom")}

if (!require("Hmisc", character.only = TRUE)) {
  install.packages("Hmisc")}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("fgsea", character.only = TRUE))
  BiocManager::install("fgsea")

if (!require("furrr", character.only = TRUE)) {
  install.packages("furrr")}

if (!require("broom", character.only = TRUE)) {
  install.packages("broom")}

if (!require("psych", character.only = TRUE)) {
  install.packages("psych")}

if (!require("openxlsx", character.only = TRUE)) {
  install.packages("openxlsx")}

if (!require("ggrepel", character.only = TRUE)) {
  install.packages("ggrepel")}

library(tidyverse)
library(data.table)
library(broom)
library(Hmisc)
library(ggpubr)
library(fgsea)
library(furrr)
library(broom)
library(psych)
library(openxlsx)
library(ggrepel)




################### UPLOAD AND PROCESS DEPMAP TABLE WITH GENE EXPRESSION

ccle_expression <- fread("Inputs/CCLE_expression.csv",sep = ",",header = T) %>% as_tibble() %>% rename(line = V1) %>%
   pivot_longer(cols = -line,names_to = "gene_name",values_to = "expression_value") %>%
   mutate(gene_name = str_replace_all(gene_name," \\(.*",""))



################## UPLOAD AND PROCESS DEPMAP TABLE WITH CRISPR SCORES

gene_effect <- fread("Inputs/CRISPR_gene_effect.csv",sep = ",",header = T) %>% as_tibble()

#Create function to tidy the data
tidy_data <- function(data,new_name){
   #Eliminate the parenthesis with the numbers from the column names containing gene_names
   colnames(data) <- str_replace_all(colnames(data)," \\(.*","")
   data %>%
      pivot_longer(cols = -DepMap_ID,names_to = "gene_name",values_to = new_name) %>%
      #REVERSE SIGN of gene effect to make the biggest effect positive
      mutate(gene_effect_CRISPR = -gene_effect_CRISPR) %>%
      rename(line = DepMap_ID)
}

#Make it tidy
gene_effect_tidy <- tidy_data(gene_effect,"gene_effect_CRISPR")

#Join the table with cell line expression to the table with the CRISPR gene effects
final_table <- ccle_expression %>%
   full_join(gene_effect_tidy,by = c("line","gene_name")) 

#IMPORTANT! By doing so, we are adding additional rows to ccle_expression, so now there are rows where expression = NA

#Basically, if we eliminate NA CRISPR_metrics from final_table, we lose 416 cell lines that do not have CRISPR scores.
#We cannot afford to do so: we want to keep ALL cell lines for which EXPRESSION WAS MEASURED! Even if there is no CRISPR score
#associated to them.

################## SELECT THE CRISPR VARIABLE WE WANT TO USE (EFFECT OR PROBABILITY)

crispr_metric <- "gene_effect_CRISPR"
metric_name <- "CRISPR_effect"

#We add a new column to the table specifying the crispr_metric we want to use
final_table <- final_table %>%
   mutate(CRISPR_metric = !!as.name(crispr_metric))

#Check all genes that have a non NA CRISPR effect for at least 1 cell line
crispr_effect_gene_list <- final_table %>%
   filter(!is.na(gene_effect_CRISPR)) %>%
   distinct(gene_name)

############### RUN WILCOX TEST TO COMPARE THE CRISPR SCORE BETWEEN BOTH POPULATIONS

#Create a function to perform the one_sided_test, testing that the average gene_effect is smaller in the control population
one_sided_test <- function(gene,data,mutant){
   if(nrow(data) != 0){
      #If there is data left on the table to perform the test
      data <- as.data.table(unique(data, by = c("line", "gene_name","CRISPR_metric","population")))
      #We take the unique levels of populations, to see if there are still 2 populations after all the filtering
      unique_populations <- unique(data$population)
      #If there is more than 1 population present in the data, then we can compare both populations
      if(length(unique_populations) > 1){
         #Count the number of observations in each of the populations we are comparing
         population_numbers <- data[, .N, by = population]
         population_numbers <- data.table(high = population_numbers$N[which(population_numbers$population == "high")],
                                          low = population_numbers$N[which(population_numbers$population == "low")])
         data$population <- factor(data$population,levels = unique_populations)
         #Detect if there are ties in both populations
         # ties <- data %>%
         #    group_by(population) %>%
         #    filter(duplicated(CRISPR_metric))
         #We use conf.int = T to be able to retrieve the sample estimate for the difference in location between both distributions
         #IMPORTANT! Wilcox is NOT testig a difference in the median between both populations, just a shift between both distributions!
         out <- wilcox.test(CRISPR_metric ~ population,data = data,alternative = "less",conf.int = T) %>%
            broom::tidy() %>% as.data.table()
         out[,c("gene","mutant","high","low") := .(gene,mutant,population_numbers$high,population_numbers$low)]
      }else{
         print("There are no 2 unique populations, so we cannot compare them")
      }
   }else{
      print("The data has 0 rows, so no hypothesis test was conducted")
   }
}

#Create function to plot boxplot of gene_effects_CRISPR of both populations (low and high expressing)
#Create a boxplot showing the CRISPR_metric for both populations, expressing and not expressing our SWI/SNF gene of interest
plot_boxplot_quantile <- function(data,gene,mutant,quantile){
   data <- data %>% mutate(population = ifelse(population == "high","high expression","low expression"))
   temp <- ggplot(data,aes(x = population,y = CRISPR_metric)) + 
      geom_boxplot() + 
      theme_classic() + 
      labs(x = "Cell line population",
           y = paste0("-",str_replace_all(metric_name,"_"," ")),
           title = paste0("Effect of KO ",gene," on high and low expressing ",mutant," lines"),
           subtitle = paste0("Quantile ",quantile)) + 
      theme(plot.title = element_text(hjust = 0.5,color = "red",face = "bold"),
            plot.subtitle = element_text(hjust = 0.5,color = "blue"))
   plot(temp)
}

#Create a function to plot the CRISPR displacement distributions corresponding to each expression group
plot_CRISPR_effect_distributions <- function(data,gene,mutant,quantile){
   #Get gene name
   gene_name <- unique(data$gene_name)
   data <- data %>%
      mutate(population = case_when(
         # population == "high" ~ paste0("High ",mutant," expression"),
         # population == "low" ~ paste0("Low ",mutant," expression")
         population == "low" ~ "Low expression",
         population == "high" ~ "High expression")
      ) %>%
      mutate(population = factor(population,levels = c("Low expression","High expression"))) %>%
      #Flip the sign of the CRISPR metric
      mutate(CRISPR_metric = -CRISPR_metric)
   #Define nbins for histogram
   nbins <- 16
   max_x <- round(max(data$CRISPR_metric),2) + 0.2
   min_x <- round(min(data$CRISPR_metric),2) - 0.2
   #Get median per expression group
   exp_median <- data %>%
      group_by(population) %>%
      dplyr::summarize(median_pop = median(CRISPR_metric))
   #Get max density values reached in the Y axis
   density <- data %>%
      group_by(population) %>%
      reframe(density = list(density(CRISPR_metric)))
   #Get max density in total
   max_y <- max(map_dbl(.x = density$density,.f = ~max(.x$y) *.x$n)) * diff(range(data$CRISPR_metric))/nbins
   #Get number of observations
   num_obs <- data %>% drop_na() %>% nrow()
   temp <- ggplot(data,aes(x = CRISPR_metric,color = population,fill = population)) + 
      #geom_histogram(position = "identity",alpha = 0.4,bins = 20) + 
      #geom_density(aes(y = stage(nbins, after_stat = count * diff(range(x))/nbins)),alpha = 0.4) +
      geom_density(aes(y = ..count.. * diff(range(x))/nbins),alpha = 0.4) 
   #Round to the closest higher number divisible by 5
   max_y <- ceiling(max_y/5) * 5 + 5
   temp_f <- temp +
      geom_vline(data = exp_median,mapping = aes(xintercept = median_pop,color = population),
                 linetype = "dashed",linewidth = 1.5) + 
      theme_classic() + 
      labs(x = "CRISPR gene KO effect",
           y = "Number of cell lines",
           title = paste0("Effect of KO ",gene," on high and low expressing ",mutant," lines"),
           subtitle = paste0("Quantile ",quantile),fill = "") + 
      scale_color_manual(values = c("blue", "red")) +
      scale_fill_manual(values = c("blue", "red"))+ 
      guides(color = "none") + 
      xlim(min_x,max_x) + 
      theme(plot.title = element_text(hjust = 0.5,color = "red",face = "bold"),
            plot.subtitle = element_text(hjust = 0.5,color = "blue"),
            axis.text.x = element_text(size = 29),
            axis.text.y = element_text(size = 29),
            axis.title.x = element_text(size = 31),
            axis.title.y = element_text(size = 31),
            legend.position = "top",
            legend.text = element_text(size = 31))
   
   plot(temp_f)
}

####################### CREATE A FUNCTION TO COMPUTE EXPRESSION QUANTILES AFTER FILTERING THE DATA FOR THE CRISPR GENE OF INTEREST

#IMPORTANT! We need first to select ONLY THE CELL LINES where our CRISPR gene was mutated and THEN we can split into LOW and HIGH expression populations
#This is the only way to wind up with BALANCED populations according to the number of cell lines that integrate them

compute_expression_quantiles <- function(mutant_tab,quantile_lim = NULL){
   a <- mutant_tab$expression_value
   a <- split(x = a,f = cut2(a,g = 20,digits = 5)) 
   #Get number of observations in each of the 20 groups
   num_obs_per_group <- map_dbl(.x = a,.f = ~length(.x))
   #Cummulative sum from 0-10 and 20 to 11
   cum_sum_1_10 <- cumsum(num_obs_per_group[1:10])
   cum_sum_20_11 <- rev(cumsum(num_obs_per_group[20:11]))
   #The names of each sublist is an interval like this [a,b)
   #We want to separate the LOWER and UPPER bound into 2 columns
   #IMPORTANT! If quantile_lim != NULL, then we just use a single quantile to define lower and upper populations
   #E.g., quantile_lim = 3 defines the 3rd and 18th quantiles to select lower and upper populations
   if(!is.null(quantile_lim)){
      correct_positions <- c(quantile_lim,length(names(a))+1-quantile_lim)
      int_names <- names(a)[correct_positions]
      cum_sum_range <- c(cum_sum_1_10[quantile_lim],cum_sum_20_11[10+1-quantile_lim])
   }else{
      correct_positions <- 1:20
      int_names <- names(a)
      cum_sum_range <- c(cum_sum_1_10,cum_sum_20_11)
   }
   separated_bounds <- rbindlist(map(int_names,.f = function(y){
      temp <- str_extract_all(string = y,pattern = "(\\d+\\.\\d{3})",simplify = T) %>%
         as.double() %>%
         t() %>%
         as.data.table()
      names(temp) <- c("lower_bound","upper_bound")
      temp
   }))
   separated_bounds[,c("partition_number","number_observations","cum_sum_num_observations"):=.(correct_positions,num_obs_per_group[correct_positions],cum_sum_range)] 
}

############### CREATE FUNCTION TO PLOT MUTANT GENE DISTRIBUTION WITH 2 POPULATIONS

plot_expression_distribution <- function(mutant_tab,z){
   #Get gene name
   gene_name <- unique(mutant_tab$gene_name)
   mutant_tab <- mutant_tab %>%
      mutate(population = as.factor(case_when(
         population == "high" ~ "High expression",
         population == "low" ~ "Low expression",
         is.na(population) ~ "middle expression",
         TRUE ~ as.character(population)
      ))) %>%
      mutate(population = factor(population,levels = c("Low expression","High expression")))
   temp <- ggplot(data = mutant_tab,aes(x = expression_value)) +
      geom_histogram(aes(fill = population,color = population),position = "identity",bins = 40,alpha = 0.4) + 
      labs(y = "Number of cell lines",
           x =expression("Gene expression value (log"[2]*"(TPM + 1))"),
           title = gene_name,
           fill = "")+
      theme_classic() + 
      guides(color = "none") + 
      scale_fill_manual(breaks = c("Low expression","High expression"),values = c("blue", "red","grey80")) +
      scale_color_manual(breaks = c("Low expression","High expression"),values = c("blue", "red","grey80")) +
      theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black",size = 33),
            plot.subtitle = element_text(hjust = 0.5,color = "blue"),
            axis.text.x = element_text(size = 29),
            axis.text.y = element_text(size = 29),
            axis.title.x = element_text(size = 31),
            axis.title.y = element_text(size = 31),
            legend.position = "top",
            legend.text = element_text(size = 31)
            )
   plot(temp)
}

#################### CREATE A FUNCTION TO SPLIT USING EVERY QUANTILE INTO TWO EXPRESSION POPULATIONS

#Define function to traverse quantiles
traverse_quantiles <- function(number_sequence,expression_tab,mutant_tab,plot_expression = F){
   temp <- expression_tab
   #Add population factor to mutant data
   mutant_tab$population <- fcase(
      mutant_tab$expression_value < temp[1,2]$upper_bound,"low",
      mutant_tab$expression_value >= temp[2,1]$lower_bound,"high"
   )
   ###### Optionally, plot the distribution of gene expression for the mutant gene, coloring the high and low populations
   if(plot_expression){
      plot_expression_distribution(mutant_tab,z)
   }
   #Once we have our lower and upper percentile bounds, we can divide our cell population into 2
   o <- mutant_tab[mutant_tab$expression_value < temp[1,2]$upper_bound | mutant_tab$expression_value >= temp[2,1]$lower_bound]
   drop_na(o)
}

split_populations <- function(mutant_tab,expression_tab,unique_percentile = F,check_size = F,plot_expression = F){
   number_sequence <- list(unique_percentile) %>% setNames(unique_percentile)
   #Traverse every quantile independently
   traverse_quantiles(number_sequence,expression_tab,mutant_tab,plot_expression = plot_expression)
}

#Create subtract function
subtract <- function(array){
   max(array) - min(array)
}

####### Create function to plot the distribution of expression of a mutant gene across cell lines

#Create function to build 2 plots in 1 figure:
#A) adjusted p-values vs quantiles
#B) diff median AND diff mean vs quantiles

plot_percentile_test <- function(data,ylimit = c(0,1)){
   #Padj vs quantiles
   subplot_a <- ggplot(data = data,aes(x = percentile, y = p_adj,color = test,shape = test)) + 
      geom_point() + 
      geom_line(aes(group = test)) + 
      lims(y = ylimit) + 
      labs(x = "Quantile",y = "Adjusted p-value (BH)",
           subtitle = paste0("Wilcox test padj vs quantile with ",metric_name)) + 
      geom_hline(aes(yintercept = 0.05),color = "red") + 
      scale_shape_manual(values=1:nlevels(data$test)) +
      scale_x_continuous(breaks = seq(0.05,0.5,by = 0.05)) +
      theme_classic() + 
      theme(plot.subtitle = element_text(hjust = 0.5,face = "bold",color = "blue"))
   #Diff mean and diff median vs quantiles
   #Prepare the data
   temp <- data %>%
      ungroup() %>%
      select(percentile,test,diff_mean,diff_median) %>% 
      pivot_longer(cols = -c(percentile,test),names_to = "statistic",values_to = "value")
   subplot_b <- ggplot(temp,aes(x = percentile,y = value)) + 
      geom_col(aes(fill = statistic),position = "dodge") + 
      facet_grid(~test) +
      labs(x = "Quantile",y = "Difference in CRISPR gene effect between populations") + 
      scale_x_continuous(breaks = seq(0.05,0.5,by = 0.05)) +
      scale_y_continuous(breaks = seq(0,0.6,by = 0.05)) + 
      theme(axis.text.x = element_text(angle = 90)) 
   #Arrange both plots into a single frame
   #ggarrange(subplot_a,subplot_b)
   print(subplot_a)
   print(subplot_b)
   invisible()
}

######################################## TEST ALL PATHWAYS IN THE REACTOME DATABASE FOR SWI/SNF MUTANTS

#Load our 3 databases gmt files and obtain the list of genes of every pathway


path_reactome_gmt <- paste0(pdir,"/Inputs/reactome_gmt_symbol.gmt")
path_reactome_table <- paste0(pdir,"/Inputs/Reactome_table_all_categories_with_genes.rds")

#Put all paths into a list and load at once
gmt_list <- gmtPathways(path_reactome_gmt)

#Load tables with genes
table_list <- readRDS(path_reactome_table)

#Input list of mutants for which we want to test the database (SWI/SNF in this case)
input_mutants <- c("ARID1A","ARID1B","ARID2","BAP1","CREBBP","EED","KMT2C","KMT2D","PBRM1","SETD2","SMARCA2","SMARCA4","SMARCB1") 
input_mutants <- map(.x = input_mutants,.f = ~.x) %>% setNames(input_mutants)



#Reactome genes
#reactome_genes <- data.frame(gene_name = reduce(gmt_list$reactome,union)) #Modified by Annabelle
reactome_genes <- data.frame(gene_name = reduce(gmt_list,union))

#Create function to take all genes from all pathways in a database and eliminate the mutant gene
get_genes_eliminate_mutant_database <- function(database,gene){
   #Collapse all pathways in database, eliminating duplicates
   all_genes <- reduce(database,union)
   #Eliminate current mutant gene from the list
   setdiff(all_genes,gene)
   #Add extra genes not present in the database that we would like to test
}

############################## CREATE FUNCTION TO RUN TEST FOR A UNIQUE PERCENTILE

#First, we will transform final_table into a data.table and we will set keys
#We do this because the table is very big and normal filtering with dplyr takes a long time
#However, if we set keys, we can do better than linear search and we can speed up the code
final_table <- as.data.table(final_table)
#Set key on the important columns we will use for search
setindex(final_table,gene_name)
#Select percentile to split populations: 2 = 10% extremes, 3 = 15 % and so on
my_percentile <- 2
run_hypothesis_test_unique_percentile <- function(mutant,gene,unique_percentile = my_percentile,create_boxplot = F,plot_expression = F,create_CRISPR_effect_dist = F){
   #Print the mutant we are processing
   #print(paste0("Processing mutant ",mutant," and gene ",gene))
   #1. We filter the information in final_table corresponding to the CRISPR gene of interest 
   #If there is NO CRISPR_metric for it, we can avoid looping
   #IMPORTANT! This step is KEY to gain computational time. Final_table is a big table, so we first filter the info
   #corresponding to the CRISPR gene of interest BEFORE left_joining. If not, the left_join operation takes longer
   crispr_data <- final_table[gene, on ="gene_name"] %>%
      #Only take the non NA values in CRISPR_metric
      filter(!is.na(CRISPR_metric))
   #Set temporary index on the line column
   setindex(crispr_data,line)
   ##IMPORTANT! At this point, if crispr_data is EMPTY, it means that the target gene of interest was not KO by CRISPR in any cell line
   #In this case, we abort the function and go to the next gene
   is_gene_KO_CRISPR <- any(!is.na(crispr_data$CRISPR_metric))
   if(is_gene_KO_CRISPR == F){
      #We create an output similar to the one given by one_sided_test but with NA and go to the next iteration
      data.frame(p.value = NA,method = NA,alternative = NA,gene = gene,mutant = mutant,high = NA,low = NA, percentile = unique_percentile)
   }else{
      #2. We filter the information in final_table corresponding to the mutant gene
      #From mutant_data we ELIMINATE those cell lines for which the CRISPR gene was NOT KO. We don't want them in our gene expression distribution
      #We can do a semi_join to take just the common lines with crispr_data, where the NA crispr_effects were eliminated
      mutant_data <- fsetdiff(final_table[mutant, on ="gene_name"],final_table[mutant, on ="gene_name"][!crispr_data, on="line"], all=T)
      #Eliminate NA values
      mutant_data <- mutant_data[!is.na(expression_value),]
      #3. We have the correct cell lines in mutant data and we can SPLIT into quantiles, to build LOW and HIGH expressing populations
      expression_quantiles <- compute_expression_quantiles(mutant_data,quantile_lim = unique_percentile)
      #4. Create two populations (LOW and HIGH) for which we want to test the difference in mean CRISPR_metric
      expression_populations <- split_populations(mutant_data,expression_quantiles,unique_percentile,plot_expression = plot_expression)
      #5. Once the two populations have been built we come back to crispr_data and ADD THE POPULATION INFORMATION for every cell line
      crispr_population_data <- crispr_data %>%
         right_join(select(expression_populations,line,population),by = "line")
      #6. Compute difference between MEDIAN of both groups and MEAN for each quantile independently
      mean_median_crispr <- crispr_population_data[, .(mean = mean(CRISPR_metric),median = median(CRISPR_metric)), by = "population"]
      mean_median_crispr <- data.table(diff_mean = subtract(mean_median_crispr$mean),diff_median = subtract(mean_median_crispr$median))
      #7. Finally, with the correct crispr_data and the population information, we can go through each quantile and compute the test
      if(create_CRISPR_effect_dist == T){
         plot_CRISPR_effect_distributions(crispr_population_data,gene,mutant,unique_percentile)
      }
      #Run the Wilcox test and get the output
      o <- one_sided_test(gene,crispr_population_data,mutant)
      #Check if the output is a successfull test or not
      if(is.character(o)){
         #If the test is not successfull, then we do not execute the rest of the code
         break
      }else{
         o <- as.data.table(o)
         o[,"percentile"] <- unique_percentile
         o[,"statistic":=NULL]
         o[,"diff_mean"] <- mean_median_crispr$diff_mean
         o[,"diff_median"] <- mean_median_crispr$diff_median
         o
      }
   }
}

# Set a "plan" for how the code should run.
plan(multisession, workers = 8)
options(future.globals.maxSize= Inf)

#Run for Reactome pruned + few additional genes of interest not present in the GMT file
genes_to_add <- c("PDCD10","CEP55","CDAN1","ESF1","TTF2","CDIN1","ZMYND8","COQ4","MARK2","RBM18","UFM1","EFR3A","UFC1","TTC27","COA5")

#Modified gmt Reactome file
gmt_list_mod <- gmt_list
gmt_list_mod$additional <- genes_to_add

run_CRISPR_database <- function(database_name,input_mutants,my_percentile){
   #Extract current database GMT
   input_database <- gmt_list_mod[[database_name]]
   #Loop through every mutant and run the entire database for each 
   future_map(.x = input_mutants,.f = function(mutant){
      print(paste0("Processing mutant ",mutant))
      #For each mutant, we eliminate the mutant gene from every pathway in the database
      temp_database <- get_genes_eliminate_mutant_database(input_database,mutant)
      #Find synthetic lethal targets for each gene independently
      curr_mutant_res <- map(.x = temp_database,.f = function(curr_gene,mutant){
         run_hypothesis_test_unique_percentile(mutant = mutant,gene = curr_gene,unique_percentile = my_percentile)
      },mutant = mutant)
      #Bind all tests of all genes together and apply a multiple testing correction with BH for all genes
      curr_mutant_res_proc <- bind_rows(curr_mutant_res) %>%
         #IMPORTANT! There are some genes which have NA in every column. These are genes that were not KO with CRISPR, so we can discard them
         filter(!is.na(estimate)) %>%
         #Correct for multiple testing for all genes in a mutant using BH method
         mutate(p_adj = p.adjust(p.value,method = "BH")) 
   })
   #.progress = T,#.options = furrr_options(packages = "broom.mixed"))
}

my_percentile <- 2
#Run algorithm
crispr_results <- run_CRISPR_database("reactome",input_mutants,my_percentile)
#Store results
saveRDS(crispr_results,paste0(pdir,"/Outputs/CRISPR_raw_output.rds"))

####################### SAVE FINAL RESULTS IN EXCEL

#Reorder mutants before saving
crispr_results_for_excell <- crispr_results %>%
   bind_rows() %>%
   mutate(mutant = factor(mutant,levels = c("ARID1A","ARID1B","ARID2","PBRM1","SMARCA2","SMARCA4","SMARCB1","BAP1","CREBBP","EED","KMT2C","KMT2D","SETD2"))) %>%
   split(.$mutant)

#Create function to save results
save_excell_res <- function(my_list){
   wb <- createWorkbook() 
   map(.x = names(my_list),.f = function(workbook_nam){
      #Select only the corresponding colums
      temp <- my_list[[workbook_nam]] %>% as_tibble() %>% dplyr::select(mutant,gene,high,low,estimate,p_adj) %>%
         dplyr::rename(num_lines_high = high,num_lines_low = low) %>%
         arrange(p_adj)
      ## Create the worksheets
      addWorksheet(wb, sheetName = workbook_nam)
      ## Write the data
      writeData(wb,workbook_nam,temp,keepNA = T)
   })
   #Once we have filled in all sheets of the workbook, we save it
   ## Save workbook to working directory 
   saveWorkbook(wb, file = paste0(pdir,"/Outputs/CRISPR_DepMap_final_analysis.xlsx"), overwrite = TRUE)
}

save_excell_res(crispr_results_for_excell)

#################### CREATE VOLCANO PLOT WITH THE DIFFERENCE IN LOCATION ESTIMATES (X AXIS) AND -LOG10(PADJ) (Y AXIS)


#Get volcano data
volcano_data <- map(.x = crispr_results,.f = function(res){
   res %>%
      #Compute -log10(p_adj) and significance scale
      mutate(`-log10_p_adj` = -log10(p_adj),
             significance_scale = case_when(
                p_adj >= 0.001 ~ "adjusted p-value > 1e-3",
                p_adj < 1e-3 & p_adj >= 1e-4 ~ "1e-3 > adjusted p-value > 1e-4",
                p_adj < 1e-4 ~ "adjusted p-value < 1e-4"
              )) %>%
      mutate(significance_scale = factor(significance_scale,levels = c("adjusted p-value > 1e-3",
                                                                       "1e-3 > adjusted p-value > 1e-4",
                                                                       "adjusted p-value < 1e-4"
                                                                       ))) %>%
      #Filter only displacements smaller than 0, as the rest are not interesting
      filter(estimate < 0) %>%
      #Also flip the sign of the Wilcox estimate
      mutate(estimate = -estimate) %>%
      #Label pathways with p_adj < 1e-3 AND estimate > 0.1
      mutate(add_label = ifelse(`-log10_p_adj` > 3 & estimate > 0.15,gene,NA)) 
})

#Choose color code for significant points
color_code <- c("grey","deepskyblue","darkorchid")

#Same but no colors and no legend
volcano_plots <- imap(.x = volcano_data,.f = function(a,mutant_name){
   #Get the max value of the -log10 of padj
   max_y_axis <- ceiling(max(a$`-log10_p_adj`))
   max_x_axis <- ceiling(max(a$estimate)*10)/10
   ggplot(a,aes(x = estimate,y = `-log10_p_adj`),color = "grey") + 
      #Add significance line at -log10_padj = 3
      #Add estimate line at estimate = 0.2
      geom_vline(aes(xintercept = 0.15),linetype = "dashed",color = "purple",size = 1.7) + 
      geom_hline(aes(yintercept = 3),linetype = "dashed",color = "purple",size = 1.7) + 
      geom_point(alpha = 0.7,shape = 16,size = 3) + 
      geom_label_repel(aes(label = add_label),show.legend = F,max.overlaps = 100,size = 6,min.segment.length = unit(0, 'lines')) + 
      labs(x = "Wilcox test displacement estimate", 
           y = bquote(-Log[10]("adjusted p-value")), 
           color = "",
           title = mutant_name) + 
      guides(color = "none") + 
      theme_classic()+ 
      scale_color_manual(values = color_code) + 
      scale_y_continuous(breaks = seq(0,max_y_axis,by = 1),limits = c(0,max_y_axis)) + 
      scale_x_continuous(breaks = seq(0,max_x_axis,by = 0.1),limits = c(0,max_x_axis)) + 
      theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "red"),
            axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 22),
            legend.text = element_text(size = 21),
            axis.text.x = element_text(size = 21),
            axis.text.y = element_text(size = 21),
            legend.position="bottom", legend.box = "horizontal")
})


############## EXTRACT LIST OF GENE HITS AND SUPER HITS FOR EVERY MUTANT

#HITS
#1. Adjusted p-value < 1e-3
#2. Displacement estimate > 0.15

#SUPER HITS
#1. Adjusted p-value < 1e-4
#2. Displacement estimate > 0.15

hits_list <- map(.x = volcano_data,.f = ~.x %>% filter(p_adj < 1E-3 & estimate > 0.15))
super_hits_list <- map(.x = volcano_data,.f = ~.x %>% filter(p_adj < 1e-4 & estimate > 0.15) %>% arrange(significance_scale))

