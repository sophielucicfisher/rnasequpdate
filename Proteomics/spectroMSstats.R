library(MSstats)
library(data.table)


setwd("E:/Sophie/Proteomics") 


# raw <- read.csv("20231207_181253_HeartTissueMaternal_CompleteReport.tsv", sep="\t")
raw <- fread("20231208_084424_HeartTissueMaternal_Schema2Report.tsv", sep="\t")


annot <- read.csv("SLF_Spectronaut_annotation.csv", header = TRUE)

head(raw)

quant <- SpectronauttoMSstatsFormat(raw,
                                    annotation = annot,
                                    filter_with_Qvalue = TRUE, ## same as default
                                    qvalue_cutoff = 0.01, ## same as default
                                    removeProtein_with1Feature = TRUE)

spectronaut.proposed <- dataProcess(raw,
                                    normalization='equalizeMedians',
                                    summaryMethod="TMP",
                                    censoredInt="NA",
                                    MBimpute=TRUE,
                                    maxQuantileforCensored=0.999)

