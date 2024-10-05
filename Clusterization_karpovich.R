library (data.table)
library(readr)
library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(clusterProfiler)
library(org.Dm.eg.db)


coli_chip_path <- 'Dmel_immunity/ChIP/Ecoli_3_vs_Control_1_peak_compare_table.tsv'
mlut_chip_path <- 'Dmel_immunity/ChIP/Mlut_5_vs_Control_1_peak_compare_table.tsv'
ecoli_rna_path <- 'Dmel_immunity/RNA/Ecoli_XIvsControl_X/Ecoli_XIvsControl_X.mRNA.gene_level.deg.xls'
mlut_rna_path <- 'Dmel_immunity/RNA/Mluteus_XIIvsControl_X/Mluteus_XIIvsControl_X.mRNA.gene_level.deg.xls'
ecoli_faire_path <-'Dmel_immunity/FAIRE/Ecoli_8_narrow_peaks_table.tsv'
mlut_faire_path <-'Dmel_immunity/FAIRE/MLuteus_9_narrow_peaks_table.tsv'
control_faire_path <-'Dmel_immunity/FAIRE/Control_7_narrow_peaks_table.tsv'

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene


# load RNA DEseq data
ecoli_rna <- fread(ecoli_rna_path)
mlut_rna <- fread(mlut_rna_path)

# load diff-peak info
ecoli_chip <- data.table(read_tsv(ecoli_chip_path))
mlut_chip <- data.table(read_tsv(mlut_chip_path))

#load faire-peak info
ecoli_faire <- data.table(read_tsv(ecoli_faire_path))
mlut_faire <- data.table(read_tsv(mlut_faire_path))
control_faire <- data.table(read_tsv(control_faire_path))

ecoli_chip$Chromosome <- paste0('chr',ecoli_chip$Chromosome)
ecoli_chip$foldchange <- sapply(ecoli_chip$log2foldchange, function(x) 2^x)

ecoli_chip_granges <- makeGRangesFromDataFrame(df = ecoli_chip,
                      ignore.strand = TRUE,
                      keep.extra.columns = TRUE,
                      seqnames.field = "Chromosome",
                      start.field = "PeakStart",
                      end.field = "PeakEnd")


covplot(ecoli_chip_granges, weightCol="foldchange")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(ecoli_chip_granges, windows=promoter)

heatmap <- peak_Profile_Heatmap(peak = ecoli_chip_granges,
                     upstream = 2000,
                     downstream = 2000,
                     by = "gene",
                     type = "start_site",
                     TxDb = txdb,
                     nbin = 800)

plotAvgProf2(ecoli_chip_granges, TxDb=txdb, upstream=2000, downstream=2000,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",conf = 0.95)

peakAnno <- annotatePeak(ecoli_chip_granges, tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Dm.eg.db")

plotAnnoPie(peakAnno)

plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

class(peakAnno)
ecoli_chip_anno <- as.data.frame(peakAnno@anno)

ecoli_chip_anno["is_promoter"] <- grepl('<=1kb', ecoli_chip_anno$annotation)

ecoli_promoters <- ecoli_chip_anno[ecoli_chip_anno$is_promoter == 'TRUE',]




