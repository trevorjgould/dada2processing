library(dada2)
args = commandArgs(trailingOnly=TRUE)
path <- (".")
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# ITSs
quality=args[1]
# good quality
if (quality == "good"){
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# poor quality
}
# bad quality
if (quality == "bad"){
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 4), truncQ = 2, minLen = 100, truncLen=c(240,175), rm.phix = TRUE, compress = TRUE, multithread = TRUE)
}
head(out)
#dereplicate reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

# error models
errF <- learnErrors(derep_forward, multithread=8, randomize=TRUE)
errR <- learnErrors(derep_reverse, multithread=8, randomize=TRUE)

dadaFs <- dada(derep_forward, err=errF, multithread=8, pool="pseudo")
dadaRs <- dada(derep_reverse, err=errR, multithread=8, pool="pseudo")

merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, trimOverhang=TRUE, minOverlap=20)

seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, "../dada2output/seqtab.rds")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=8, verbose=TRUE)
dim(seqtab.nochim)
lname <- nchar(colnames(seqtab.nochim))
seqtab.nochim <- seqtab.nochim[,(lname > 280)]
saveRDS(seqtab.nochim, "../dada2output/seqtab_nochim.rds")

uniquesToFasta(seqtab.nochim, fout = "../dada2output/sequences.fasta")

  # set a little function
getN <- function(x) sum(getUniques(x))

  # making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
               filtered=out[,2], dada_f=sapply(dadaFs, getN),
               dada_r=sapply(dadaRs, getN), merged=sapply(merged_amplicons, getN),nonchim=rowSums(seqtab.nochim))
write.table(summary_tab, file = "../dada2output/sequence_process_summary.txt", sep = "\t", quote=FALSE)

ref <- "/home/umii/public/dada2_taxonomy_references/sh_general_release_dynamic_all_04.04.2024.fasta"
taxa <- assignTaxonomy(seqtab.nochim, ref, multithread = TRUE, outputBootstraps = TRUE)

taxout <- taxa$tax
bootout <- taxa$boot
saveRDS(taxout, file = "../dada2output/taxID.rds")
saveRDS(bootout, file = "../dada2output/taxID_bootstrap.rds")

#saveRDS(taxa$tax, file = "taxa.rds")
#saveRDS(taxa$boot, file = "taxa_bootstrap.rds")
#
both1 <- cbind(t(seqtab.nochim),taxa$tax, taxa$boot)
write.table(both1, file = "../dada2output/ITS_combined_sequences_taxa_bootstrap.txt", sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxa$tax)
write.table(both2, file = "../dada2output/ITS_combined_sequences_taxa.txt", sep = "\t", quote = FALSE, col.names=NA)