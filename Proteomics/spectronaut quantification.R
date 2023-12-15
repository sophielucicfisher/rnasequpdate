library(data.table)
setwd("E:/Sophie/Proteomics/Spectronaut Report")
ibaqdf <- fread("20231212_112836_HeartTissueMaternal_MLPivotReport.tsv")
nopivot <- fread("20231212_112836_HeartTissueMaternal_MLReport.tsv")

`# Clean acession names (thanks Jun!), not working
has_repeating_strings <- function(row){
  parts <- strsplit(row, ";")[[1]]
  if (length(parts) > 1){
    first_part <- sub("-.*", "", parts[1])
    all(grepl(paste0("^", first_part, "-\\d+$"), parts[-1]))
  }else{
    FALSE
  }
}

rownames_matrix <- matrix(rownames(ibaqdf))
rownames(rownames_matrix)<- rownames(ibaqdf)
selected_rows <- ibaqdf[apply(ibaqdf, 1, has_repeating_strings),]
uncertain_rows <- grepl(";", row.names(ibaqdf$PG.ProteinGroups))
ibaqdf <- ibaqdf[!uncertain_rows,]
ibaqdf <- rbind(ibaqdf, selected_rows)
write.csv(ibaqdf, "MLPivotReport_Cleaned.csv")


# Michael's Hallucination Playground (Mid-Coding Crisis, destroy)

## Exclude all columns that are outside of the quantification columns (after column 101)
ibaqdf_quantification <- ibaqdf[,1:101]

## Delete protein group if after ; is different
### Bad concept but if the first three characters between the group (after the semicolon) are different, unlikely to be isomers so delete them
ibaqdf_clean <- ibaqdf_quantification[!grep(";",ibaqdf_quantification$PG.ProteinGroups)] # deletes if semicolon is different (too evasive?)

### Possibly cleaner concept but less efficient (seem to not work? same results as the bad concept)
repeating_proteins <- function(x){
  split_proteins <- strsplit(x, ";")[[1]]
  character_check <- substr(split_proteins[1], 1, 3) # check if first three characters are the same
  ### Using a for loop with if function to check the first three characters (very, very inefficient)
  for (part in split_proteins){
    if (substr(part, 1, 3) != character_check){
      return(TRUE)
    }
  }
  return(FALSE) # when they have the same 3 characters
  ### Row exclusion if meet above criteria for protein group
  ibaqdf_clean <- ibaqdf_quantification[!sapply(ibaqdf_quantification$PG.ProteinGroups, repeating_proteins),]
}

# Delete row if protein group is the same and 70% of the quantification results are the same
duplicate_value <- function(row){
  length(unique(row)) / length(row) <= 0.3
}
ibaqdf_clean_duplicate <- ibaqdf_clean[!duplicated(ibaqdf_clean$PG.ProteinGroups) & !apply(ibaqdf_clean[, 6:101], 1, duplicate_value),]

# Reshape df to long format
ibaqdf_reshape <- melt(ibaqdf_clean_duplicate, variable.name = "sample", value.name = "quantification")




# Michael's Hallucination Playground 2 (DATA PIVOTED)
qdf <- fread("20231212_112836_HeartTissueMaternal_MLReport.tsv")
metadata <- readxl::read_xlsx("sample_metadata.xlsx")

# Delete protein groups (replicated from failure, dammit NOT WORKING)
  repeating_proteins <- function(x){
    split_proteins <- strsplit(x, ";")[[1]] # split strings with ;
    character_check <- substr(split_proteins[1], 1, 3) # check if first three characters are the same
    ### Using a for loop with if function to check the first three characters (very, very inefficient)
    for (part in split_proteins){
      if (substr(part, 1, 3) != character_check){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }
  }
  ### Row exclusion if meet above criteria for protein group
  qdf_exclude <- qdf[!sapply(qdf$PG.ProteinAccessions, repeating_proteins),]

# Delete protein groups when they have semicolon (very destructive, to be resolved in the future)
qdf_exclude <- qdf[!grep(";",qdf$PG.ProteinAccessions)]

# Join metadata with sample
qdf_exclude_metadata <- merge(x = qdf_exclude, y = metadata,
                              by.y = c("sample_id"),
                              by.x = c("R.FileName"),
                              all.x = TRUE)

