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
library(EnhancedVolcano)


ecoli_chip_path <- 'Dmel_immunity/ChIP/Ecoli_3_vs_Control_1_peak_compare_table.tsv'
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

length(unique(ecoli_chip$ClosestTSS_ID)) # посмотрели на количество аннотированных генов в обработке до меня
nrow(ecoli_chip) # получается 7124 индивидуальных гена, но при этом некоторые пики аннотированны дважды

ecoli_chip$Chromosome <- paste0('chr',ecoli_chip$Chromosome) # поменяли немного названия
mlut_chip$Chromosome <- paste0('chr',mlut_chip$Chromosome)

ecoli_chip$foldchange <- sapply(ecoli_chip$log2foldchange, function(x) 2^x) # превратили log2FC в FC
mlut_chip$foldchange <- sapply(mlut_chip$log2foldchange, function(x) 2^x)

length(unique(ecoli_chip$PeakStart)) # посмотрели на кол-во уникальных пиков
ecoli_chip <- ecoli_chip[!duplicated(ecoli_chip$PeakStart), ] # убрали дупликаты среди пиков
mlut_chip <- mlut_chip[!duplicated(mlut_chip$PeakStart), ]


ecoli_chip_granges <- makeGRangesFromDataFrame(df = ecoli_chip, # превратили в объект grange для работы с chipseeker
                      ignore.strand = TRUE,
                      keep.extra.columns = TRUE,
                      seqnames.field = "Chromosome",
                      start.field = "PeakStart",
                      end.field = "PeakEnd")

mlut_chip_granges <- makeGRangesFromDataFrame(df = mlut_chip,
                                               ignore.strand = TRUE,
                                               keep.extra.columns = TRUE,
                                               seqnames.field = "Chromosome",
                                               start.field = "PeakStart",
                                               end.field = "PeakEnd")


#covplot(ecoli_chip_granges, weightCol="foldchange") # график покрытия хромосом пиками chip-seq

promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000) # получили список промотеров 
tagMatrixEcoli <- getTagMatrix(ecoli_chip_granges, windows=promoter) # что такео TagMatrix я пока не разбирался
tagMatrixMlut <- getTagMatrix(mlut_chip_granges, windows=promoter)


#peak_Profile_Heatmap(peak = ecoli_chip_granges, # график расстояний до TSS
 #                     upstream = 2000,
  #                     downstream = 2000,
   #                    by = "gene",
    #                   type = "start_site",
     #                  TxDb = txdb,
      #                 nbin = 800)

#plotAvgProf2(ecoli_chip_granges, TxDb=txdb, upstream=2000, downstream=2000, 
#             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",conf = 0.95)

peakAnnoEcoli <- annotatePeak(ecoli_chip_granges, tssRegion=c(-2000, 2000), # аннотировали наши пики
                         TxDb=txdb, annoDb="org.Dm.eg.db")
peakAnnoMlut <- annotatePeak(mlut_chip_granges, tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Dm.eg.db")

#plotAnnoPie(peakAnno)

#plotDistToTSS(peakAnno, # графики для понимания происходящего
 #             title="Distribution of transcription factor-binding loci\nrelative to TSS")

ecoli_chip_anno <- as.data.frame(peakAnnoEcoli@anno)  # берем только часть с аннотациями
mlut_chip_anno <- as.data.frame(peakAnnoMlut@anno)

length(unique(ecoli_chip_anno$geneId)) # 6974
nrow(ecoli_chip_anno) # 10772 радуемся тому что кол-во генов меньше чем кол-ва пиков


# получаем только пики лежащие в областях близких к промотерам
ecoli_chip_anno["is_promoter"] <- grepl('<=1kb', ecoli_chip_anno$annotation) 
mlut_chip_anno["is_promoter"] <- grepl('<=1kb', mlut_chip_anno$annotation)


ecoli_promoters <- ecoli_chip_anno[ecoli_chip_anno$is_promoter == 'TRUE',]
mlut_promoters <- mlut_chip_anno[mlut_chip_anno$is_promoter == 'TRUE',]

nrow(ecoli_promoters)



# получаем пики которые выросли в ~1.4 (корень из 2) раза по сравнению с контролем
ecoli_promoters_increase <- ecoli_promoters[ecoli_promoters$log2foldchange > 0.5,]
mlut_promoters_increase <- mlut_promoters[mlut_promoters$log2foldchange > 0.5,]

# смотрим на их количество
nrow(ecoli_promoters_increase)
nrow(mlut_promoters_increase)

# получаем набор генов,рядом с которыми находятся наиболее возросшие пики
ecoli_genes_FC05 <- unique(ecoli_promoters_increase$geneId)
mlut_genes_FC05 <- unique(mlut_promoters_increase$geneId)

length(unique(ecoli_genes_FC05))
length(intersect(mlut_genes_FC05,ecoli_genes_FC05)) # пересечение в 406 генах для Ecoli и Mluteus

candidates <- intersect(mlut_genes_FC05,ecoli_genes_FC05) # получим список генов 
# рядом с которыми растут пики chip-seq у Ecoli и Mluteus

# top_rna <- ecoli_rna[['gene_id']] %in% candidates

# посмотрим что происходит с экспрессией в этих генах
top_rna <- ecoli_rna[ecoli_rna$gene_id %in% candidates,] # получаем те гены, про которые есть информация в rna-seq у Ecoli 
dim(top_rna) # в рнк есть информация про 309 из 406 интересующих нас генов
length(unique(candidates)) 
length(unique(ecoli_rna$gene_id)) # всего у нас в рнк секе генов


## Simple function for plotting a Volcano plot, returns a ggplot object
deseq.volcano <- function(res, datasetName) {
  return(EnhancedVolcano(res, x = 'log2FoldChange', y = 'padj',
                         lab=res$gene_name,
                         title = paste(datasetName, "RNA genes near high FC chip-seq peaks"),
                         subtitle = bquote(italic('FDR <= 0.05 and absolute FC >= 2')),
                         # Change text and icon sizes
                         labSize = 3, pointSize = 1.5, axisLabSize=10, titleLabSize=12,
                         subtitleLabSize=8, captionLabSize=10,
                         # Disable legend
                         legendPosition = "none",
                         # Set cutoffs
                         pCutoff = 0.05, FCcutoff = 2))
}

## Note: input data is the corrected DESeq2 output using the 'lfcShrink' function (see chapter 4)
deseq.volcano(res = top_rna, datasetName = "Drosophila")
