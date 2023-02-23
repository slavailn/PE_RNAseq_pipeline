library(regionReport)
library(ggplot2)
setwd("<path/to/dir>")


list.files("expression_data/")
## Epididymis_Control_insulated_vs_Control_non-insulated
reportDir <- "Epididymis_Control_insulated_vs_Control_non-insulated"
DESeqObj <- "expression_data/Epididymis_Control_insulated_vs_Control_non-insulated_diff.RData"
projTitle <- "Epididymis_Control_insulated_vs_Control_non-insulated"

createRegionReport <- function(reportDir, DESeqObj, projTitle) {
  ## The output will be saved in the DESeq2 report directory
  dir.create(reportDir, showWarnings = FALSE, recursive = TRUE)
  
  load(DESeqObj)
  
  ## Generate the HTML report
  report <- DESeq2Report(dds_diff, project = projTitle, c("Group"),
                         outdir = reportDir)  
}

list.files("expression_data/")

# 1) "Epididymis_Control_insulated_vs_Control_non-insulated_diff.RData"

createRegionReport(reportDir="Epididymis_Control_insulated_vs_Control_non-insulated", 
                   DESeqObj="expression_data/Epididymis_Control_insulated_vs_Control_non-insulated_diff.RData",
                   projTitle="Epididymis_Control_insulated_vs_Control_non-insulated")


