# WORKING BUT INEFFICIENT
# Michael's Coding Purgatory (thanks Sophie!)

#---- set the environment ---- 
setwd("E:/Sophie/Proteomics/Spectronaut Report")


#---- call libraries ----

library(data.table)
library(dplyr)
library(tidyr)
library(limma)
library(openxlsx)
library(EnhancedVolcano)
pacman::p_load(dplyr,ggplot2,data.table,tidyr,limma)

#---- in house functions ----
# Delete all isoformRT
# function only keeps the first value
first_value <- function(input){
  return(sub(";.*", "", input))
}

# function only keeps the first gene
first_gene <- function(input){
  return(sub(";.*", "", input))
}

# function generates a volcano plot, requirement for topTable
gvolcanop <- function(result_table, interaction){
  volcanop <- EnhancedVolcano(result_table,
                              lab = rownames(fit2df),
                              x = "logFC",
                              y = "adj.P.Val",
                              title = dietary_interactions[a],
                              pCutoff = 0.05,
                              FCcutoff = 1,
                              drawConnectors = TRUE,
                              widthConnectors = 0.75,
                              max.overlaps = 100,
                              labSize = 6.0
  )
  return(volcanop)
}

#---- main code ----
nopivot <- fread("20231212_112836_HeartTissueMaternal_MLReport.tsv")

# pull only the required information to reduce size
newdf <- data.frame(nopivot$R.FileName, nopivot$PG.Genes, nopivot$PG.Quantity) 

# pivot the df to only take file name (sample name) and quantity
newdf1 <- newdf %>% pivot_wider(
  names_from = nopivot.R.FileName, 
  values_from = nopivot.PG.Quantity
)

# only pull the first result and ignore the rest
testdf_clean <- newdf %>% 
  mutate(
    quantity = sapply(nopivot.PG.Quantity, first_value)
  )

# working, only select sample name, gene, quantity
clean1 <- testdf_clean %>% select(1,2,4) 

# loading the various groups
group = read.csv("samples.csv") 

# pivot the data to a wide format
clean_pivot <- reshape(data = clean1,idvar = "nopivot.PG.Genes", v.names = "quantity", timevar = "nopivot.R.FileName", direction = "wide")

# pull the sample id and make it as a row name, once completed make it a null
colnames(clean_pivot)[colnames(clean_pivot) == "nopivot.PG.Genes"] <- "SampleID"
rownames(clean_pivot) <- clean_pivot$SampleID
clean_pivot$SampleID <- NULL

# ensure that clean_pivot df is a matrix
exp <- clean_pivot %>% as.matrix() 

#Average the replicated proteins
exp <- avereps(exp)

exp = data.frame(exp)

# create combined experimental factor with totalGroup
combineFactor <- factor(paste(group$group,group$Sex,sep="."))
totalGroup <- cbind(combineFactor,group=group)

#change exp to numeric
exp <- apply(
  X = exp, 
  MARGIN = 1, 
  FUN = as.numeric
)

logexp <- log2(exp) # log transfomation

# possible factor levels of dietary group
dietary_group <- factor(totalGroup$group.group, levels = c("Chow", "Sucrose", "CS", "SC"))
table(dietary_group)

# create a model matrix using the dietary_group factor
design <- model.matrix(~0+dietary_group)
colnames(design) <- levels(factor(dietary_group))
#rownames(design) <- colnames(logexp)

# specify the dietary interactions as a list
dietary_interactions <- as.factor(c("ChowvSucrose","ChowvCS","ChowvSC","SucrosevCS","SucrosevSC","CSvsSC"))

# transpose the results after log for limma
t_logexp <- t(logexp)

# need to create an empty list first
results_list <- list()

# dictate all the possible contrasts/dietary interactions
contrasts_fix <- makeContrasts(
  ChowvSucrose = Chow-Sucrose,
  ChowvCS = Chow-CS,
  ChowvSC = Chow-SC,
  SucrosevCS = Sucrose-CS,
  SucrosevSC = Sucrose-SC,
  CSvsSC = CS-SC,
  levels = design
)

# conduct all the necessary processing
fit <- lmFit(t_logexp, design)
fit2 <- contrasts.fit(fit, contrasts_fix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,coef=1, n = Inf, adjust.method = "BH") # only chowvsucrose results are taken due to coef=1
DEG <- na.omit(DEG)
#DEG$Change <- ifelse(DEG$P.Value > 0.05, "unchanged",
#ifelse(DEG$logFC > 0.5, "up",
#ifelse(DEG$logFC < -0.5, "down", "unchanged")))

# creating for loops to generate volcano plot for every dietary interactions (working!)
plot_list <- list()
for (a in 1:nlevels(dietary_interactions)){
  interaction <- as.character(dietary_interactions[a])
  result_table <- as.data.frame(topTable(fit2, coef = a, number = 5000, adjust.method = "BH"))
  volcanop <- gvolcanop(result_table, interaction)
  pdf(paste0("plot_",gsub("/","-", interaction),".pdf")) # ensure the results are not stored in a corrupted manner
  print(volcanop)
  dev.off()
  plot_list[[interaction]] <- volcanop
}

print(table(DEG$Change))
results_list[[diet]] <- DEG

# CSM <- results_list$"ChowvsSucrose"
# results_list$"CS - SC"

##---- export dataframe as excel 2 (working!, ignore)----

interaction_diet <- c("Chow - Sucrose", "Chow - CS", "Chow - SC", "Sucrose - SC", "Sucrose - CS", "CS - SC")
results_list <- list() #need to create an empty list first
for (diet in interaction_diet) {
  fit <- lmFit(expr, design = design)
  write.xlsx(fit, file = paste0(diet, "_prelmfit.xlsx"), rowNames = TRUE)
  contrast_matrix <- makeContrasts(diet, levels = colnames(design))
  fit2 <- contrasts.fit(fit, contrast_matrix)
  write.xlsx(fit2, file = paste0(diet, "_contrastfit.xlsx"), rowNames = TRUE)
  fit2 <- eBayes(fit2)
  DEG <- topTable(fit2, coef = diet, n = Inf, sort.by = "logFC", adjust.method = "BH")
  DEG <- na.omit(DEG)
  DEG$Change <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 0.5, "up",
                              ifelse(DEG$logFC < -0.5, "down", "unchanged")))
  
  print(table(DEG$Change))
  write.xlsx(DEG, file = paste0(diet, "_DEG.xlsx"), rowNames = TRUE)
  results_list[[diet]] <- DEG
}







#---- fix required for dietary interaction (fail the pub test)----
combineFactor <- factor(paste(group$group,group$Sex,sep=".")) #create combined experimental factor
totalGroup <- cbind(combineFactor,group=group) # create a total group df

#use biobase expression set?
clean1_naomit <- clean1 #copy df
clean1_naomit$quantity <- as.numeric(clean1_naomit$quantity) #changing quant to numeric
clean1_naomit <- clean1_naomit[complete.cases(clean1_naomit),] #only keep complete cases as DGE doesn't accept NA
clean1_pivot <- reshape(data = clean1_naomit,idvar = "nopivot.PG.Genes", v.names = "quantity", timevar = "nopivot.R.FileName", direction = "wide") #pivot
#designNew <- model.matrix(~0+totalGroup$combineFactor)


designNew <- model.matrix(~totalGroup$group.Sex*totalGroup$group.group) #notsure of this, need to check for diet and sex interactions


rownames(designNew) <- totalGroup$group.SampleID #add the rownames to the designNew as the sampleID
rownames(clean1_pivot) <- clean1_pivot$nopivot.PG.Genes #rownames of clean1_pivot, changed to genes
clean1_pivot <- clean1_pivot[-c(1)] #remove gene names
groups <- as.factor(totalGroup$combineFactor) #change both sex and diet as factors

y <- DGEList(counts=clean1_pivot,group=groups) #possibly working


#--

#---- not sure of the fit testing (discard) ----

y <- normLibSizes(y)
y <- estimateDisp(y,designNew)


fit <- glmQLFit(y,designNew)
qlf <- glmQLFTest(fit,coef=c(1,7)) #comparison chow
topTags(qlf)


#---- unused ----

nopivot <- fread("20231212_112836_HeartTissueMaternal_MLReport.tsv")

newdf <- data.frame(nopivot$R.FileName, nopivot$PG.Genes, nopivot$PG.Quantity)

newdf1 <- newdf %>% pivot_wider(
  names_from = nopivot.R.FileName, 
  values_from = nopivot.PG.Quantity
)



#yo try stringr (very technical :( ), rebus (if give-up) library
#only pickup the first value (function, remove after ";"), ngl too evasive here
testdf_clean <- newdf %>% 
  mutate(
    quantity = sapply(nopivot.PG.Quantity, first_value)
  )

clean1 <- testdf_clean %>% select(1,2,4) # working, only select sample name, gene, quantity


#gene = sapply(nopivot.PG.Genes, first_gene))

# [END]
group = read.csv("samples.csv") #groups loading



#---- SEX FILTER (REPLACE SEX in line 41) ----
# i think this is wrong as sex should be a separate group as it may interact with the results (confounders)
# DONT SPLIT IT UP
#filtered_sexgroup <- group[group$Sex == "F",] 
#sampleid_add <- sub("^quantity\\.","",filtered_sexgroup$SampleID)
#sampleid_add_df <- data.frame(SampleID = sampleid_add)
#clean1_filtered <- clean1[clean1$nopivot.R.FileName %in% sampleid_add,]
# END SEX FILTER
# Note: filtered_sexgroup is group and clean1_filtered is the value list
#---- old results, possible discard ----
testdf_clean_pivot <- reshape(data = clean1_filtered,idvar = "nopivot.PG.Genes", v.names = "quantity", timevar = "nopivot.R.FileName", direction = "wide")

colnames(testdf_clean_pivot)[colnames(testdf_clean_pivot) == "nopivot.PG.Genes"] <- "SampleID" #rename sampleid

# discard, possible big mistake, sex is ignored

exp <- testdf_clean_pivot %>% as.matrix() #ensure the exp is from testdfcleanpivot and matrix

#This is to make the protein name as row name
rownames(exp) = exp[,1] 
exp <- exp[,-1]

#Average the replicated proteins
exp <- avereps(exp)

#exp <-  exp[which(rowMeans(!is.na(exp)) > 0.80), ] #delete proteins with more than 80% missing values, modify to filter by groups

exp = data.frame(exp) #ensure dataframe iscreated


exp <- apply(
  X = exp, 
  MARGIN = 1, 
  FUN = as.numeric
)

for (i in 1:ncol(exp)) {
  exp[,i]  = as.numeric(exp[,i]) #change to numeric
}

expr <- log2(exp) # log transfomation
#expr = normalizeBetweenArrays(expr)  #normalize the data - not necessary (if needed can take off the hashtag)

dietary_group <- factor(filtered_sexgroup$group, levels = c("Sucrose", "Chow", "CS", "SC"))

table(dietary_group)

design <- model.matrix(~0 + dietary_group + dietary_group:dietary_group)

#design <- model.matrix(~0+dietary_group)

colnames(design) <- levels(factor(filtered_sexgroup$group))
rownames(design) <- colnames(expr)



#contrast.matrix <- makeContrasts(Disease-Control,levels = design) #set group



# Create for loop for every diet interaction 
interaction_diet <- c("Chow - Sucrose", "Chow - CS", "Chow - SC", "Sucrose - SC", "Sucrose - CS", "CS - SC")
results_list <- list() #need to create an empty list first
for (diet in interaction_diet) {
  contrast_matrix <- makeContrasts(diet, levels = colnames(design))
  fit <- lmFit(expr, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  DEG <- topTable(fit2, coef = diet, n = Inf, adjust.method = "BH")
  DEG <- na.omit(DEG)
  #DEG$Change <- ifelse(DEG$P.Value > 0.05, "unchanged",
  #ifelse(DEG$logFC > 0.5, "up",
  #ifelse(DEG$logFC < -0.5, "down", "unchanged")))
  
  print(table(DEG$Change))
  results_list[[diet]] <- DEG
}

CSM <- results_list$"Chow - Sucrose"
results_list$"CS - SC"

# end discard

### volcano plot
library(EnhancedVolcano)
library(tidyr)


colnames(CSM)
row.names(CSM)

head(CSM)

EnhancedVolcano(CSM,
                lab = rownames(CSM),
                x = "logFC",
                y = "adj.P.Val",
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                max.overlaps = 100,
                labSize = 6.0)


fix_fit <- glmQLFit(y,design=designNew)

contrasts_fix <- makeContrasts(
  ChowvSucrose = Chow-Sucrose,
  ChowvsCS = Chow-CS,
  ChowvsSC = Chow-SC,
  SucrosevsSC = Sucrose-CS,
  CSvsSC = CS-SC,
  levels = design
)

design = model.matrix(~0+Group)
fix_fit <- lmFit(expr, design = design)
fix_fit2 <- contrasts.fit(fit, contrasts_fix)
fix_fit3 <- eBayes(fix_fit2)
fix_deg <- topTable(fix_fit3, coef = diet, n = Inf, sort.by = "logFC", adjust.method = "BH")
fix_deg <- na.omit(fix_deg)
fix_deg$Change <- ifelse(DEG$P.Value > 0.05, "unchanged",
                     ifelse(DEG$logFC > 0.5, "up",
                            ifelse(DEG$logFC < -0.5, "down", "unchanged")))


y <- DGEList(genes = clean1_pivot$nopivot.PG.Genes)

y <- DGEList(counts = clean1_naomit$quantity, samples = clean1_naomit$nopivot.R.FileName, genes = clean1_naomit$nopivot.PG.Genes)



clean1_naomit$quantity <- na.omit(clean1_naomit)
clean1_naomit$quantity <- as.numeric(clean1_naomit$quantity)
summary(clean1_naomit$quantity)
y <- DGEList(counts=as.numeric(clean1_naomit$quantity), samples=clean1_pivot)

clean1_pivot <- clean1_pivot[complete.cases(clean1_pivot),] 
clean1_pivot <- clean1_pivot %>% mutate(across(-c(nopivot.PG.Genes),as.numeric))
y <- DGEList(counts=clean1_pivot)
