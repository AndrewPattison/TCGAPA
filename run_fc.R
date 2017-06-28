#!/usr/bin/Rscript
suppressMessages(library(limma))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(edgeR))
suppressMessages(library(Rsubread))

number <- commandArgs(trailingOnly = T)[2]

setwd(paste0("~/bioinformatics/Full_TCGA_BRCA_analysis/downloaded_bams/", number))

apa_db <-read.delim("~/bioinformatics/Full_TCGA_BRCA_analysis/poly_a_site_db/HG_38_full_gene_overlap_removed.bed", header=F, comment.char="",stringsAsFactors=F)
# 
# apa_db$V7 <- gsub("\\..*", "", apa_db$V4)


# gff <- data.frame(apa_db$V1, "Rprog", "exon", as.integer(apa_db$V2), as.integer(apa_db$V3), ".", apa_db$V6,".", paste0("gene_id=",apa_db$V4))

# write.table(x = gff,sep = "\t", quote = F,file = "~/bioinformatics/Full_TCGA_BRCA_analysis/poly_a_site_db/HG_38_full_gene_overlap_removed.gff",row.names = F, col.names = F)

uniq <- unique(apa_db[,4])

apadb_unique <- apa_db [!duplicated(apa_db$V4),]

annotatum <- data.frame(GeneID = apadb_unique[,4],
                        Chr = apadb_unique[,1],
                        Start = apadb_unique[,2],
                        End = apadb_unique[,3],
                        Strand = apadb_unique[,6],
                        stringsAsFactors=FALSE)

# write.table(x = annotatum,sep = "\t", quote = F,file = "~/bioinformatics/Full_TCGA_BRCA_analysis/poly_a_site_db/processed.gtf",row.names = F, col.names = F)
file_run <- list.files(paste0("/data/home/apattison/bioinformatics/Full_TCGA_BRCA_analysis/new_bams/", number), pattern = "*.bam$",full.names = T)[1]

file2 <- gsub("/data/home/apattison/bioinformatics/Full_TCGA_BRCA_analysis/new_bams*.*/", "",file_run)

all_APA <- featureCounts(minFragLength = 40,largestOverlap = T, annot.ext = annotatum , strandSpecific = 0,useMetaFeatures = F, nthreads = 1, files = file_run, isPairedEnd = T, autosort = T,countMultiMappingReads = F)

saveRDS(object = all_APA, paste0("/data/home/apattison/bioinformatics/Full_TCGA_BRCA_analysis/output_RDS/", file2, "_counts.RDS"))

cat("Feature counts has counted", file2)
