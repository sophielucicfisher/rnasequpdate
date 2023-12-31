library(data.table)
library(dplyr)
library(tidyr)
library(limma)
pacman::p_load(dplyr,ggplot2,data.table,tidyr,limma)


setwd("E:/Sophie/Proteomics/Spectronaut Report")
nopivot <- fread("20231212_112836_HeartTissueMaternal_MLReport.tsv")

newdf <- data.frame(nopivot$R.FileName, nopivot$PG.Genes, nopivot$PG.Quantity)



newdf1 <- newdf %>% pivot_wider(names_from = nopivot.R.FileName, 
                                values_from = nopivot.PG.Quantity)
# Delete all isoforms [START]
first_value <- function(input){
  return(sub(";.*", "", input))
}

first_gene <- function(input){
  return(sub(";.*", "", input))
}

testdf_clean <- data.frame()

testdf_clean <- newdf %>% 
  mutate(
    quantity = sapply(nopivot.PG.Quantity, first_value))
#gene = sapply(nopivot.PG.Genes, first_gene))

# [END]


clean1 <- testdf_clean %>% select(1,2,4) # working

testdf_clean_pivot <- reshape(data = clean1,idvar = "nopivot.PG.Genes", v.names = "quantity", timevar = "nopivot.R.FileName", direction = "wide")

colnames(testdf_clean_pivot)[colnames(testdf_clean_pivot) == "nopivot.PG.Genes"] = "SampleID"

group = read.csv("samples.csv")
exp = testdf_clean_pivot %>% as.matrix()

#This is to make the protein name as row name
rownames(exp) = exp[,1] 
exp = exp[,-1]

#Average the replicated proteins
exp = avereps(exp)

exp =  exp[which(rowMeans(!is.na(exp)) > 0.75), ] #delete proteins with more than 75% missing values
exp = data.frame(exp)
str(exp)


for (i in 1:ncol(exp)) {
  exp[,i]  = as.numeric(exp[,i]) #change to numeric
}

expr = log2(exp+1) # log transfomation
#expr = normalizeBetweenArrays(expr)  #normalize the data - not necessary (if needed can take off the hashtag)

dietary_group <- factor(group$group)

design <- model.matrix(~0 + dietary_group + dietary_group:dietary_group)

#design <- model.matrix(~0+dietary_group)

colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr)



#contrast.matrix <- makeContrasts(Disease-Control,levels = design) #set group



# Create for loop for every diet interaction 
interaction_diet <- c("Chow - Sucrose", "Chow - CS", "Chow - SC", "Sucrose - SC", "Sucrose - CS", "CS - SC")
for (diet in interaction_diet) {
  fit <- lmFit(expr, design)
  contrast_matrix <- makeContrasts(diet, levels = colnames(design))
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  DEG <- topTable(fit2, coef = diet, n = Inf, sort.by = "logFC")
  DEG <- na.omit(DEG)
  DEG$Change <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 0.5, "up",
                              ifelse(DEG$logFC < -0.5, "down", "unchanged")))
  
  print(table(DEG$Change))
  results_list[[diet]] <- DEG
}

results_list$"Chow - Sucrose"
results_list$"Chow - CS"
