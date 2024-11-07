##Packages -----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread) #featureCounts

##featureCounts -----
fc <- featureCounts(files = c("alignment_output/B1_AlignedReads.sam","alignment_output/B2_AlignedReads.sam","alignment_output/B3_AlignedReads.sam", "alignment_output/C1_AlignedReads.sam","alignment_output/C2_AlignedReads.sam","alignment_output/C3_AlignedReads.sam"),
              annot.inbuilt = "hg38",
              useMetaFeatures = TRUE,
              isPairedEnd = TRUE,
              nthreads = 8)
View(fc$annotation)
View(fc$stat)
View(fc$counts)

fc$counts["4322",] #MMP13 counts
fc$counts["22806",]
