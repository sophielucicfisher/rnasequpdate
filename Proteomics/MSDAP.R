install.packages(c("devtools", "tidyverse", "tinytex", "BiocManager"))
tinytex::install_tinytex()
# On Windows; say 'no' to optionally compile packages and during TinyTex installation you may see 2 popups; these can be dismissed
BiocManager::install(c('ProtGenerics', 'MSnbase', 'limma'), update=T, ask=F)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
devtools::install_github("ftwkoopmans/msdap", upgrade = "never") # don't update dependencies if not needed

library(msdap)
setwd("E:/Sophie/Proteomics")

# Starting MSDAP
library(msdap)
library(data.table)
dataset <- msdap::import_dataset_spectronaut("20231208_151017_HeartTissueMaternal_MSDAP2Report.tsv")
  protein_qval_threshold = 0.05, # default is used
  collapse_peptide_by = "sequence_modified" # default is used

dataset = import_fasta(
  dataset,
  files = c("uniprotkb_proteome_UP000000589_2023_11_22.fasta"))
  
###write_template_for_sample_metadata(dataset, "sample_metadata.xlsx")  

dataset = import_sample_metadata(dataset, "sample_metadata.xlsx")

dataset = setup_contrasts(dataset, 
                          contrast_list = list(c("Chow", "Sucrose"), c("Chow", "CS"), c("Chow", "SC"), c("Sucrose", "SC") , c("Sucrose", "CS"), c("CS", "SC"))
)

dataset = setup_contrasts(dataset, 
                          contrast_list = list(c("M","F"), c("Chow", "Sucrose"), c("Chow", "CS"), c("Chow", "SC"), c("Sucrose", "SC") , c("Sucrose", "CS"), c("CS", "SC"))
                          )
                          
dataset = analysis_quickstart(
  dataset,
  filter_min_detect = 3,
  filter_min_quant = 3,
  filter_fraction_detect = 0.75,
  filter_fraction_quant = 0.75,
  filter_topn_peptides = 0, 
  filter_min_peptide_per_prot = 1,
  filter_by_contrast = TRUE,
  norm_algorithm = c("vsn", "modebetween_protein"),
  rollup_algorithm = "maxlfq",N
  dea_algorithm = c("deqms", "msempire", "msqrob"),
  dea_qvalue_threshold = 0.05,
  dea_log2foldchange_threshold = NA,
  diffdetect_min_samples_observed = 3,
  output_qc_report = TRUE,
  output_abundance_tables = TRUE,
  output_dir = "E:/Sophie/Proteomics/MSDAPRerun",
  output_within_timestamped_subdirectory = TRUE,
  dump_all_data = FALSE
)


print_dataset_summary(dataset)

library(xlsx)

x = summarise_stats(
  dataset,
  return_dea = TRUE,
  return_diffdetect = TRUE,
  dea_logfc_as_effectsize = FALSE,
  diffdetect_zscore_threshold = 5,
  remove_ambiguous_proteingroups = TRUE
)
