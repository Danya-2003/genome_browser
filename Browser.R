### RP11-398K22.12_03_both

data1 <- "~/Downloads/Genome_browser/parsing/RP11_03_both_peak_minus_logfc_minus.csv"
data2 <- "~/Downloads/diffbind/RP11-398K22.12_03_no_adj.diffbind_sicer"
out_peaks_fdr <- "~/Downloads/Genome_browser/out/RP11_03_peaks_fdr.bed"
out_peaks <- "~/Downloads/Genome_browser/out/RP11_03_peaks.bed"
out_genes <- "~/Downloads/Genome_browser/out/RP11_03_genes.bed"

####################################

diffbind <- read.table(data2, header =TRUE, sep = "\t")

diffbind_logfc_minus_fdr <- diffbind[(diffbind$Fold < 0) & (diffbind$FDR < 0.05), ]
diffbind_logfc_minus <- diffbind[(diffbind$Fold < 0), ]

peaks_logfc_minus_fdr <- as.data.frame(diffbind_logfc_minus_fdr$seqnames)
peaks_logfc_minus_fdr$chromStart <- diffbind_logfc_minus_fdr$start
peaks_logfc_minus_fdr$chromEnd <- diffbind_logfc_minus_fdr$end
peaks_logfc_minus_fdr$name <- "peak_name"
peaks_logfc_minus_fdr$score <-  "100"
peaks_logfc_minus_fdr$strand <- "."
peaks_logfc_minus_fdr$thickStart <- diffbind_logfc_minus_fdr$start
peaks_logfc_minus_fdr$thickEnd <- diffbind_logfc_minus_fdr$end
peaks_logfc_minus_fdr$itemRgb <- "166,0,0"

write.table(data.frame(peaks_logfc_minus_fdr), out_peaks_fdr, sep = '\t',  quote = F,  row.names = F, col.names = F)

peaks_logfc_minus <- as.data.frame(diffbind_logfc_minus$seqnames)
peaks_logfc_minus$chromStart <- diffbind_logfc_minus$start
peaks_logfc_minus$chromEnd <- diffbind_logfc_minus$end
peaks_logfc_minus$name <- "peak_name"
peaks_logfc_minus$score <-  "100"
peaks_logfc_minus$strand <- "."
peaks_logfc_minus$thickStart <- diffbind_logfc_minus$start
peaks_logfc_minus$thickEnd <- diffbind_logfc_minus$end
peaks_logfc_minus$itemRgb <- "166,0,0"

write.table(data.frame(peaks_logfc_minus), out_peaks, sep = '\t',  quote = F,  row.names = F, col.names = F)

chippeakann <- read.table(data1, header =TRUE, sep = ",")

genes <- as.data.frame(chippeakann$seqnames)
colnames(genes)[colnames(genes) == 'chippeakann$seqnames'] <- 'num'
genes$chr <- 'chr'
genes$chrom <-paste0(genes$chr, genes$num)
genes$chromStart <- chippeakann$start
genes$chromEnd <- chippeakann$end
genes$name <- chippeakann$feature
genes$score <-  chippeakann$logfc
genes$strand <- chippeakann$feature_strand
genes$thickStart <- chippeakann$chromStart
genes$thickEnd <- chippeakann$chromEnd
genes$itemRgb <- "166,0,0"  
genes <- genes[,c(3,4,5,6,7,8,9)]

write.table(data.frame(genes), out_genes, sep = '\t',  quote = F,  row.names = F, col.names = F)
rm(list = ls())
################

### RP11-398K22.12_05_both

data1 <- "~/Downloads/Genome_browser/parsing/RP11_05_both_peak_minus_logfc_minus.csv"
data2 <- "~/Downloads/diffbind/RP11-398K22.12_05_no_adj.diffbind_sicer"
out_peaks_fdr <- "~/Downloads/Genome_browser/out/RP11_05_peaks_fdr.bed"
out_peaks <- "~/Downloads/Genome_browser/out/RP11_05_peaks.bed"
out_genes <- "~/Downloads/Genome_browser/out/RP11_05_genes.bed"

####################################

diffbind <- read.table(data2, header =TRUE, sep = "\t")

diffbind_logfc_minus_fdr <- diffbind[(diffbind$Fold < 0) & (diffbind$FDR < 0.05), ]
diffbind_logfc_minus <- diffbind[(diffbind$Fold < 0), ]

peaks_logfc_minus_fdr <- as.data.frame(diffbind_logfc_minus_fdr$seqnames)
peaks_logfc_minus_fdr$chromStart <- diffbind_logfc_minus_fdr$start
peaks_logfc_minus_fdr$chromEnd <- diffbind_logfc_minus_fdr$end
peaks_logfc_minus_fdr$name <- "peak_name"
peaks_logfc_minus_fdr$score <-  "100"
peaks_logfc_minus_fdr$strand <- "."
peaks_logfc_minus_fdr$thickStart <- diffbind_logfc_minus_fdr$start
peaks_logfc_minus_fdr$thickEnd <- diffbind_logfc_minus_fdr$end
peaks_logfc_minus_fdr$itemRgb <- "166,0,0"

write.table(data.frame(peaks_logfc_minus_fdr), out_peaks_fdr, sep = '\t',  quote = F,  row.names = F, col.names = F)

peaks_logfc_minus <- as.data.frame(diffbind_logfc_minus$seqnames)
peaks_logfc_minus$chromStart <- diffbind_logfc_minus$start
peaks_logfc_minus$chromEnd <- diffbind_logfc_minus$end
peaks_logfc_minus$name <- "peak_name"
peaks_logfc_minus$score <-  "100"
peaks_logfc_minus$strand <- "."
peaks_logfc_minus$thickStart <- diffbind_logfc_minus$start
peaks_logfc_minus$thickEnd <- diffbind_logfc_minus$end
peaks_logfc_minus$itemRgb <- "166,0,0"

write.table(data.frame(peaks_logfc_minus), out_peaks, sep = '\t',  quote = F,  row.names = F, col.names = F)

chippeakann <- read.table(data1, header =TRUE, sep = ",")

genes <- as.data.frame(chippeakann$seqnames)
colnames(genes)[colnames(genes) == 'chippeakann$seqnames'] <- 'num'
genes$chr <- 'chr'
genes$chrom <-paste0(genes$chr, genes$num)
genes$chromStart <- chippeakann$start
genes$chromEnd <- chippeakann$end
genes$name <- chippeakann$feature
genes$score <-  chippeakann$logfc
genes$strand <- chippeakann$feature_strand
genes$thickStart <- chippeakann$chromStart
genes$thickEnd <- chippeakann$chromEnd
genes$itemRgb <- "166,0,0"  
genes <- genes[,c(3,4,5,6,7,8,9)]

write.table(data.frame(genes), out_genes, sep = '\t',  quote = F,  row.names = F, col.names = F)
rm(list = ls())
################

### CTD_01_both

data1 <- "~/Downloads/Genome_browser/parsing/CTD_01_both_peak_minus_logfc_minus.csv"
data2 <- "~/Downloads/diffbind/CTD-2587H24.5_01_no_adj.diffbind_sicer"
out_peaks_fdr <- "~/Downloads/Genome_browser/out/CTD_01_peaks_fdr.bed"
out_peaks <- "~/Downloads/Genome_browser/out/CTD_01_peaks.bed"
out_genes <- "~/Downloads/Genome_browser/out/CTD_01_genes.bed"

####################################

diffbind <- read.table(data2, header =TRUE, sep = "\t")

diffbind_logfc_minus_fdr <- diffbind[(diffbind$Fold < 0) & (diffbind$FDR < 0.05), ]
diffbind_logfc_minus <- diffbind[(diffbind$Fold < 0), ]

peaks_logfc_minus_fdr <- as.data.frame(diffbind_logfc_minus_fdr$seqnames)
peaks_logfc_minus_fdr$chromStart <- diffbind_logfc_minus_fdr$start
peaks_logfc_minus_fdr$chromEnd <- diffbind_logfc_minus_fdr$end
peaks_logfc_minus_fdr$name <- "peak_name"
peaks_logfc_minus_fdr$score <-  "100"
peaks_logfc_minus_fdr$strand <- "."
peaks_logfc_minus_fdr$thickStart <- diffbind_logfc_minus_fdr$start
peaks_logfc_minus_fdr$thickEnd <- diffbind_logfc_minus_fdr$end
peaks_logfc_minus_fdr$itemRgb <- "166,0,0"

write.table(data.frame(peaks_logfc_minus_fdr), out_peaks_fdr, sep = '\t',  quote = F,  row.names = F, col.names = F)

peaks_logfc_minus <- as.data.frame(diffbind_logfc_minus$seqnames)
peaks_logfc_minus$chromStart <- diffbind_logfc_minus$start
peaks_logfc_minus$chromEnd <- diffbind_logfc_minus$end
peaks_logfc_minus$name <- "peak_name"
peaks_logfc_minus$score <-  "100"
peaks_logfc_minus$strand <- "."
peaks_logfc_minus$thickStart <- diffbind_logfc_minus$start
peaks_logfc_minus$thickEnd <- diffbind_logfc_minus$end
peaks_logfc_minus$itemRgb <- "166,0,0"

write.table(data.frame(peaks_logfc_minus), out_peaks, sep = '\t',  quote = F,  row.names = F, col.names = F)

chippeakann <- read.table(data1, header =TRUE, sep = ",")

genes <- as.data.frame(chippeakann$seqnames)
colnames(genes)[colnames(genes) == 'chippeakann$seqnames'] <- 'num'
genes$chr <- 'chr'
genes$chrom <-paste0(genes$chr, genes$num)
genes$chromStart <- chippeakann$start
genes$chromEnd <- chippeakann$end
genes$name <- chippeakann$feature
genes$score <-  chippeakann$logfc
genes$strand <- chippeakann$feature_strand
genes$thickStart <- chippeakann$chromStart
genes$thickEnd <- chippeakann$chromEnd
genes$itemRgb <- "166,0,0"  
genes <- genes[,c(3,4,5,6,7,8,9)]

write.table(data.frame(genes), out_genes, sep = '\t',  quote = F,  row.names = F, col.names = F)
rm(list = ls())
################

### CTD_03_both

data1 <- "~/Downloads/Genome_browser/parsing/CTD_03_both_peak_minus_logfc_minus.csv"
data2 <- "~/Downloads/diffbind/CTD-2587H24.5_03_no_adj.diffbind_sicer"
out_peaks_fdr <- "~/Downloads/Genome_browser/out/CTD_03_peaks_fdr.bed"
out_peaks <- "~/Downloads/Genome_browser/out/CTD_03_peaks.bed"
out_genes <- "~/Downloads/Genome_browser/out/CTD_03_genes.bed"

####################################

diffbind <- read.table(data2, header =TRUE, sep = "\t")

diffbind_logfc_minus_fdr <- diffbind[(diffbind$Fold < 0) & (diffbind$FDR < 0.05), ]
diffbind_logfc_minus <- diffbind[(diffbind$Fold < 0), ]

peaks_logfc_minus_fdr <- as.data.frame(diffbind_logfc_minus_fdr$seqnames)
peaks_logfc_minus_fdr$chromStart <- diffbind_logfc_minus_fdr$start
peaks_logfc_minus_fdr$chromEnd <- diffbind_logfc_minus_fdr$end
peaks_logfc_minus_fdr$name <- "peak_name"
peaks_logfc_minus_fdr$score <-  "100"
peaks_logfc_minus_fdr$strand <- "."
peaks_logfc_minus_fdr$thickStart <- diffbind_logfc_minus_fdr$start
peaks_logfc_minus_fdr$thickEnd <- diffbind_logfc_minus_fdr$end
peaks_logfc_minus_fdr$itemRgb <- "166,0,0"

write.table(data.frame(peaks_logfc_minus_fdr), out_peaks_fdr, sep = '\t',  quote = F,  row.names = F, col.names = F)

peaks_logfc_minus <- as.data.frame(diffbind_logfc_minus$seqnames)
peaks_logfc_minus$chromStart <- diffbind_logfc_minus$start
peaks_logfc_minus$chromEnd <- diffbind_logfc_minus$end
peaks_logfc_minus$name <- "peak_name"
peaks_logfc_minus$score <-  "100"
peaks_logfc_minus$strand <- "."
peaks_logfc_minus$thickStart <- diffbind_logfc_minus$start
peaks_logfc_minus$thickEnd <- diffbind_logfc_minus$end
peaks_logfc_minus$itemRgb <- "166,0,0"

write.table(data.frame(peaks_logfc_minus), out_peaks, sep = '\t',  quote = F,  row.names = F, col.names = F)

chippeakann <- read.table(data1, header =TRUE, sep = ",")

genes <- as.data.frame(chippeakann$seqnames)
colnames(genes)[colnames(genes) == 'chippeakann$seqnames'] <- 'num'
genes$chr <- 'chr'
genes$chrom <-paste0(genes$chr, genes$num)
genes$chromStart <- chippeakann$start
genes$chromEnd <- chippeakann$end
genes$name <- chippeakann$feature
genes$score <-  chippeakann$logfc
genes$strand <- chippeakann$feature_strand
genes$thickStart <- chippeakann$chromStart
genes$thickEnd <- chippeakann$chromEnd
genes$itemRgb <- "166,0,0"  
genes <- genes[,c(3,4,5,6,7,8,9)]

write.table(data.frame(genes), out_genes, sep = '\t',  quote = F,  row.names = F, col.names = F)
rm(list = ls())
################

data <- "~/Downloads/Genome_browser/Genome_browser/RP11_398K22.himorna_peaks.bed"
output1 <- "~/Downloads/Genome_browser/Genome_browser/RP11_himorna_minus.bed"
output2 <- "~/Downloads/Genome_browser/Genome_browser/RP11_himorna_plus.bed"

file <- read.table(data, header = FALSE)
peak_minus <- subset(file, file$V5 < 0)
peak_plus <- subset(file, file$V5 > 0)
write.table(data.frame(peak_minus), output1, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(peak_plus), output2, sep = '\t',  quote = F,  row.names = F, col.names = F)

data <- "~/Downloads/Genome_browser/Genome_browser/FGD5_AS1.himorna_peaks.bed"
output1 <- "~/Downloads/Genome_browser/Genome_browser/FGD5_himorna_minus.bed"
output2 <- "~/Downloads/Genome_browser/Genome_browser/FGD5_himorna_plus.bed"

file <- read.table(data, header = FALSE)
peak_minus <- subset(file, file$V5 < 0)
peak_plus <- subset(file, file$V5 > 0)
write.table(data.frame(peak_minus), output1, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(peak_plus), output2, sep = '\t',  quote = F,  row.names = F, col.names = F)












