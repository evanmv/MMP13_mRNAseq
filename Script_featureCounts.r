##Packages -----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

##featureCounts
fC <- featureCounts(files = c("alignment_output/B1_AlignedReads.sam","alignment_output/B2_AlignedReads.sam","alignment_output/B3_AlignedReads.sam", "alignment_output/C1_AlignedReads.sam","alignment_output/C2_AlignedReads.sam","alignment_output/C3_AlignedReads.sam"),
              annot.inbuilt = "hg38",
              useMetaFeatures = TRUE,
              isPairedEnd = TRUE,
              nthreads = 8)
View(fC$annotation)
View(fC$stat)
