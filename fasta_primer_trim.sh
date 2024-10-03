#!/bin/bash
# single end fasta trim intended for taxonomy database trimming
# requires bbtools in your PATH
# usage
# To call the script:
# bash fasta_primer_trim.sh -p forwardprimer -q reverseprimer
# primer_trim.sh -f file.fasta -p GTGCCAGCMGCCGCGGTAA -q GGACTACHVGGGTWTCTAAT
while getopts "f:p:q:" OPTION;
do
        case $OPTION in
        f)
        fasta=$OPTARG
        ;;
        p)
        forward=$OPTARG
        ;;
        q)
        reverse=$OPTARG
        ;;
  esac
done
#primer trimming
echo "Primer 1: ${forward}"
echo "Primer 2: ${reverse}"

# remove adapters
module load cutadapt
module load parallel

# remove primers
cutadapt --cores 8 --minimum-length 100 --discard-untrimmed -g ${forward} -G ${reverse} -o exact.fasta ${fasta} > cutadapt.primer.log.txt

chmod +x run_cutadapt2.cmd
parallel < run_cutadapt2.cmd

mkdir primer_logs
mv *log.txt primer_logs
grep "passing" primer_logs/* > summary_primer_trimming.txt
echo "primers trimmed from fasta"

# dedupe

dedupe.sh in=exact.fasta out=dedupe.fasta mergenames=t exact=t
#/home/umii/goul0109/bbmap/

echo "sequences deduplicated and fasta headers merged"

# headers will now look like this with > separating taxonomy of merged sequences
# >Bacteria;Actinobacteria;Actinobacteria;Pseudonocardiales;Pseudonocardiaceae;Saccharopolyspora;>Bacteria;Actinobacteria;Actinobacteria;Pseudonocardiales;Pseudonocardiaceae;Pseudonocardia;
# to
# >Bacteria;Actinobacteria;Actinobacteria;Pseudonocardiales;Pseudonocardiaceae;Saccharopolyspora;>;Pseudonocardia;

awk '{delete seen; c=0; for (i=1;i<=NF;i++) if (!seen[$i]++) printf "%s%s", (++c>1?OFS:""), $i; print ""}' dedupe.fasta > exactDB_${fasta}

echo "output complete: exactDB_${fasta}"