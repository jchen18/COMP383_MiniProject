library(sleuth)
library(dplyr)
#setwd("~/MiniProject")
stab <- read.table("hcmv_table.txt", header = TRUE, stringsAsFactors = FALSE)
#initialize sleuth object
so <- sleuth_prep(stab)
#fitting model comparing the 2
so <- sleuth_fit(so, ~timepoint, 'full')
#fitting reduced model to compare in likelihood ratio test
so <- sleuth_fit(so, ~1, 'reduced')
#performing likelihood ratio test for differential expression between conditions
so <- sleuth_lrt(so,'reduced','full')

#extracting results for the most significant
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)
#filter most significant results and sort by pval
sleuth_significant <- dplyr::filter(sleuth_table, qval<=0.05) %>% dplyr::arrange(pval)
#write FDR < 0.05 transcripts to file
write.table(sleuth_significant[, c(1, 4, 2, 3)], file="miniProject_Jessie_Chen/fdr05_results.txt", quote = FALSE, row.names = FALSE,sep = '\t')

