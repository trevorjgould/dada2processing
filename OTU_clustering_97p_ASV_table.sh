# install vsearch
module load conda
module load R/4.4.0-openblas-rocky8
source /common/software/install/migrated/anaconda/python3-2023.03-libmamba/etc/profile.d/conda.sh
conda create --name vsearch -c conda-forge -c bioconda bioconda::vsearch
conda activate vsearch

R libraries required
tidyr
dada2
dplyr
data.table

    YucSoil16S 
    YucSoil18S 
    CLIMUSHITS > vsearch done
    fabsoil18S > done
    fabsoilITS > vsearch done
    fabsoil16S > done





getSeqs.R
########################################
# gets sequences as a fasta from .rds file
# arg[1] = seqtab_nochim.R
args = commandArgs(trailingOnly=TRUE)
seqtab.nochim <- readRDS(args[1])
dada2::uniquesToFasta(seqtab.nochim, fout = "sequences.fasta")
########################################

otuSeqsID.R
########################################
# get OTU ID and sequence matching ID
# arg[1] = all.otutab.txt
args = commandArgs(trailingOnly=TRUE)
all.otutab <- read.delim(args[1])
library(tidyr)
library(dplyr)
all.otutab[all.otutab == 0] <- NA
melted <- all.otutab %>% pivot_longer(cols = !X.OTU.ID, names_to = c("variable"), values_to = c("value"), values_drop_na = TRUE)
melted <- as.data.frame(melted)
melted$variable <- gsub("sq","",melted$variable)
melted <- melted[order(as.numeric(melted$variable)),]
write.table(melted, file = "OTUID_SEQID.txt", sep = "\t", quote = FALSE)
########################################

OTUmerge.R
########################################
# arg1 = "seqtab_nochim.rds"
# arg2 = "OTUID_SEQID.txt"
args = commandArgs(trailingOnly=TRUE)
seqtab <- readRDS(args[1])
seqtab = t(seqtab)
OTUID <- read.delim(args[2])
seqtab.OTU <- data.frame(OTUID$X.OTU.ID,seqtab)
colnames(seqtab.OTU)[1] <- "OTU.ID"
out <- setDT(seqtab.OTU)[, lapply(.SD, sum), keyby = OTU.ID]
write.table(out, file = "OTU97percentID.txt", sep = "\t", quote = FALSE, row.names = FALSE)
########################################

taxaFromFasta.R
########################################
#' get taxonomy prediction of a fasta file
#' arg[1] is fasta file name
#' arg[2] is one of [16S|ITS|18S]
#' 
#' @export

args = commandArgs(trailingOnly=TRUE)

seqs <- Biostrings::readDNAStringSet(args[1])
#16S
ref16S <- "/home/umii/public/dada2_taxonomy_references/silva_nr99_v138.1_train_set.fa"
#18S
ref18S <- "/home/umii/public/dada2_taxonomy_references/maarjam_dada2.txt"
#ITS
ITSref <- "/home/umii/public/dada2_taxonomy_references/sh_general_release_dynamic_all_04.04.2024.fasta"

if (args[2]=="16S") {
taxa <- dada2::assignTaxonomy(seqs, ref16S, multithread = TRUE, outputBootstraps = TRUE)
} else if(args[2]=="18S") {
taxa <- dada2::assignTaxonomy(seqs, ref18S, multithread = TRUE, outputBootstraps = TRUE)
} else if(args[2]=="ITS"){
taxa <- dada2::assignTaxonomy(seqs, ITSref, multithread = TRUE, outputBootstraps = TRUE)
}
taxout <- taxa$tax
bootout <- taxa$boot
write.table(taxout, file = "taxaID.txt", sep = "\t", quote = FALSE)
write.table(bootout, file = "taxaIDbootstrap.txt", sep = "\t", quote = FALSE)
########################################

Running pipeline
########################################
module load conda
module load R/4.4.0-openblas-rocky8
source /common/software/install/migrated/anaconda/python3-2023.03-libmamba/etc/profile.d/conda.sh
conda activate vsearch

Rscript /home/umii/goul0109/scripts/getSeqs.R "seqtab_nochim.rds"
vsearch --cluster_size sequences.fasta --threads 8 --id 0.97 --strand plus --sizein --sizeout --fasta_width 0 --uc all.clustered.uc --relabel OTU_ --centroids all.otus.fasta --otutabout all.otutab.txt
Rscript /home/umii/goul0109/scripts/otuSeqsID.R "all.otutab.txt"
Rscript /home/umii/goul0109/scripts/OTUmerge.R "seqtab_nochim.rds" "OTUID_SEQID.txt"
# here you need to specify the [ITS/18S/16S] reference type
Rscript /home/umii/goul0109/scripts/taxaFromFasta.R "all.otus.fasta" "18S"