library(data.table)
library(ggplot2)
setwd("Z:/PRJ-MMM/Proteomics/MSDAP Result/2023-12-08_23-25-11") 

dea <- read.csv("differential_abundance_analysis.csv", check.names = FALSE, encoding = "UCS-2LE")

colnames(dea)

head(dea)

volcano_plot <-dea %>% filter(foldchange.log2_msempire_contrast: Chow vs Sucrose!="NA")%>% 
  ggplot(aes(x = 'foldchange.log2_msempire_contrast: Chow vs Sucrose',
                y = -log('pvalue_deqms_contrast: Chow vs Sucrose'))) +
  geom_point(color="grey") +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw() 

## MICHAEL'S FAILURE

dea$test1 <- as.numeric(dea$`foldchange.log2_deqms_contrast: Chow vs Sucrose`, na.rm = TRUE)

dea$test_contrast <- as.numeric(dea$`pvalue_deqms_contrast: Chow vs Sucrose`, na.rm = TRUE)

volcano_plot <- ggplot(data = dea, aes(x = test1,
             y = -log(test_contrast))) +
  geom_point(color="grey") +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw() 

volcano_plot()

## END

