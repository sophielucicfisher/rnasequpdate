library(EnhancedVolcano)
library(tidyr)

# setwd("Z:/PRJ-MMM/Proteomics/MSDAP Result/2023-12-08_23-25-11") 
setwd("E:/Sophie/Proteomics/MSDAP/2023-12-08_23-25-11")

dea <- read.csv("differential_abundance_analysis.csv", check.names = FALSE, encoding = "UCS-2LE")

colnames(dea)

head(dea)


dea.df<- as.data.frame(dea)

dea.df <-dea.df %>%
  filter("foldchange.log2_deqms_contrast: Chow vs Sucrose" !="NA") %>% 
  filter("pvalue_deqms_contrast: Chow vs Sucrose" != "NA") %>%
  filter("pvalue_deqms_contrast: Chow vs Sucrose" != "0") 



dea.df<- as.data.frame(filtered)


dea.df$"foldchange.log2_deqms_contrast: Chow vs Sucrose" <- as.numeric(dea.df$"foldchange.log2_deqms_contrast: Chow vs Sucrose")

dea.df$"pvalue_deqms_contrast: Chow vs Sucrose" <- as.numeric(dea.df$"pvalue_deqms_contrast: Chow vs Sucrose")


EnhancedVolcano(dea.df,   x = "foldchange.log2_deqms_contrast: Chow vs Sucrose",
                y = -log("pvalue_deqms_contrast: Chow vs Sucrose"),
                lab = dea.df$gene_symbols_or_id)


EnhancedVolcano(dea.df,   x = "foldchange.log2_msempire_contrast: Chow vs Sucrose",
                y = -log("pvalue_deqms_contrast: Chow vs Sucrose"),
                         lab = dea.df$gene_symbols_or_id)

EnhancedVolcano(dea.df,
                lab = rownames(dea),
                x = 'foldchange.log2_msempire_contrast: Chow vs Sucrose',
                y = -log('pvalue_deqms_contrast: Chow vs Sucrose'))



## MICHAEL'S UNKNOWN PROBLEMS
dea$test1 <- as.numeric(dea$`foldchange.log2_msempire_contrast: Chow vs Sucrose`, na.omit = TRUE)

dea$test_contrast <- as.numeric(dea$`pvalue_deqms_contrast: Chow vs Sucrose`, na.omit = TRUE)



dea_naomit <- dea[complete.cases(dea[,c("foldchange.log2_msempire_contrast: Chow vs Sucrose","pvalue_deqms_contrast: Chow vs Sucrose")]),]
#dea_naomit <- na.omit(dea[, c("foldchange.log2_msempire_contrast: Chow vs Sucrose","pvalue_deqms_contrast: Chow vs Sucrose")])
dea_naomit_hard <- na.omit(dea)


EnhancedVolcano(dea,
                x = dea$test1,
                y = -log(dea$test_contrast))

EnhancedVolcano(dea_naomit,
                lab = rownames(dea_naomit),
                x = "foldchange.log2_msempire_contrast: Chow vs Sucrose",
                y = -log("pvalue_deqms_contrast: Chow vs Sucrose"))

dea_naomit$pvalue <- as.numeric(as.character(dea_naomit$`pvalue_deqms_contrast: Chow vs Sucrose`))
dea_naomit$foldchange <- as.numeric(as.character(dea_naomit$`foldchange.log2_msempire_contrast: Chow vs Sucrose`))
dea_naomit$row <- as.numeric(rownames(dea_naomit))

dea_naomit$gene_symbols_or_id <- as.character(dea_naomit_hard$gene_symbols_or_id)

EnhancedVolcano(dea_naomit,
                lab = dea_naomit$row,
                x = dea_naomit$foldchange,
                y = dea_naomit$pvalue)

### HARDCORE NA DELETION
dea_naomit_hard <- na.omit(dea)
dea_naomit_hard$pvalue <- as.numeric(as.character(dea_naomit_hard$`pvalue_msempire_contrast: Chow vs Sucrose`))
dea_naomit_hard$foldchange <- as.numeric(as.character(dea_naomit_hard$`foldchange.log2_msempire_contrast: Chow vs Sucrose`))
dea_naomit_hard$row <- as.numeric(rownames(dea_naomit_hard))
EnhancedVolcano(dea_naomit_hard,
                lab = dea_naomit_hard$row,
                x = dea_naomit_hard$foldchange,
                y = dea_naomit_hard$pvalue)

### WORKING PROTOTYPE DONT KILL IT
EnhancedVolcano(dea_naomit_hard,
                lab = "row",
                x = "foldchange",
                y = "pvalue")

## MICHAEL'S DEVIL TRIAL (partial na.omission)

##msempire successful
dea_naomit$pvalue <- as.numeric(as.character(dea_naomit$`pvalue_msempire_contrast: Chow vs Sucrose`))
dea_naomit$foldchange <- as.numeric(as.character(dea_naomit$`foldchange.log2_msempire_contrast: Chow vs Sucrose`))
dea_naomit$gene <- as.character(dea_naomit$gene_symbols_or_id)
dea_naomit$row <- as.numeric(rownames(dea_naomit))

EnhancedVolcano(dea_naomit,
                lab = dea_naomit$gene,
                x = "foldchange",
                y = "pvalue",
                FCcutoff = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                boxedLabels = TRUE,
                labSize = 6.0,
                title = "MSEmpire")

##deqms successful
dea_naomit$pvalue_deqms <- as.numeric(as.character(dea_naomit$`pvalue_deqms_contrast: Chow vs Sucrose`))
dea_naomit$foldchange_deqms <- as.numeric(as.character(dea_naomit$`foldchange.log2_deqms_contrast: Chow vs Sucrose`))
dea_naomit$gene_deqms <- as.character(dea_naomit$gene_symbols_or_id)
dea_naomit$row <- as.numeric(rownames(dea_naomit))

EnhancedVolcano(dea_naomit,
                lab = dea_naomit$gene,
                x = "foldchange_deqms",
                y = "pvalue_deqms",
                FCcutoff = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                boxedLabels = TRUE,
                labSize = 6.0,
                title = "DEQMS")

##msqrob
dea_naomit$pvalue_msqrob <- as.numeric(as.character(dea_naomit$`pvalue_msqrob_contrast: Chow vs Sucrose`))
dea_naomit$foldchange_msqrob <- as.numeric(as.character(dea_naomit$`foldchange.log2_msqrob_contrast: Chow vs Sucrose`))
dea_naomit$gene_msqrob <- as.character(dea_naomit$gene_symbols_or_id)
dea_naomit$row <- as.numeric(rownames(dea_naomit))

EnhancedVolcano(dea_naomit,
                lab = dea_naomit$gene,
                x = "foldchange_msqrob",
                y = "pvalue_msqrob",
                FCcutoff = 0.221,
                pCutoff = 0.05,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                boxedLabels = TRUE,
                labSize = 6.0,
                title = "MSQROB")
