#!/bin/bash
Help(){
# Display Help
	echo "This script checks what 16S, ITS, or 18S primers you have"
	echo "primercheck.sh sample_R1_L001.fastq.gz sample_R2_L001.fastq.gz"
}	
while getopts ":h" OPTION;
do
    case $OPTION in
        h | * )
        Help >&2
        exit 1
        ;;
    esac
done

        while IFS=' ' read -r PRIME LOC || [ -n "${PRIME}" ]; do
        check=$(zgrep --max-count=100 -c ${PRIME} $1);
        if [ $check == 100 ]; then
          echo "$PRIME"
          R1=$(echo $PRIME | sed 's/\.\*/N/g')
          break 2
        exit
fi
done < /home/umii/goul0109/scripts/list_to_searchR1.txt

        while IFS=' ' read -r PRIME LOC || [ -n "${PRIME}" ]; do
        check=$(zgrep --max-count=100 -c ${PRIME} $2);
        if [ $check == 100 ]; then
            R2=$(echo $PRIME | sed 's/\.\*/N/g')
            echo "$PRIME"
            break 2
         exit
    fi
    done < /home/umii/goul0109/scripts/list_to_searchR2.txt

echo "./primer_trim.sh -p ${R1} -q ${R2}"

# usage primercheck3.sh forwardread reverseread
# list_to_searchR1.txt
AGAGTTTGATC.*TGGCTCAG #V1F_27F
CCTACGGG.*GGC.*GCAG #V3F_341F
GTGCCAGC.*GCCGCGGTAA #V4F_515F
.*GGATTAGATACCC #V5F_784F
TCGATGAAGAACGCAGCG #ITS1F
GCTGCGTTCTTCATCGATGC #ITS1F
TCCTCCGCTTATTGATATGC #ITS4_Nextera_(ITS2)
CAGCCGCGGTAATTCCAGCT #18SF WANDA1

# list_to_searchR2.txt note reverse order of forward file
CGAC**CCATGCA.*CACCT #V6R_1064R
GGACTAC.*.*GGGT.*TCTAAT #V4R_806R
ATTACCGCGGCTGCTGG #V3R_518R
TCGATGAAGAACGCAGCG #5.8SR_Nextera_(ITS2)
GAACCCAAACACTTTGGTTTCC #18SR AML2
