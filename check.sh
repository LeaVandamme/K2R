#!/bin/bash

# A PARTIR DU FICHIER DE SORTIE
# VERIFIER À LA MAIN SI SEQ OU KMER DE SEQ PRESENTS OU PAS

# ON PREND LES SEQ À QUERY
awk 'BEGIN {RS=">"; ORS=""} NR > 1 {print ">"$0}' $2 | grep -v "^>" > target_sequences.txt

# NUMEROTATION READS
awk '/^>/ {print ">" i++; next} {print}' $1 > numeroted.fasta

# RECHERCHE DES READS

# Extraire les séquences correspondant aux identifiants avec leurs positions
line=$(head -n 1 target_sequences.txt)
seqkit grep -s -p $line numeroted.fasta -o selected.fasta
seqtk seq -l 0 selected.fasta > selected_sequences.fasta


# TRI DES DEUX FICHIERS

awk '/^>/ {if (seq) print seq; print ""; print $0; seq=""; next} {seq=seq $0} END {if (seq) print seq}' $3 | awk 'BEGIN {RS=""; FS="\n"} {header=$1; seq=""; for (i=2; i<=NF; i++) seq=seq $i; print header, seq}' | sort -k1,1 | awk '{header=$1; seq=$2; print header; print seq}' > queryresults_sorted.fa
awk '/^>/ {if (seq) print seq; print ""; print $0; seq=""; next} {seq=seq $0} END {if (seq) print seq}' selected_sequences.fasta | awk 'BEGIN {RS=""; FS="\n"} {header=$1; seq=""; for (i=2; i<=NF; i++) seq=seq $i; print header, seq}' | sort -k1,1 | awk '{header=$1; seq=$2; print header; print seq}' > selectedsequences_sorted.fa

# DIFF

diff -s -B  queryresults_sorted.fa selectedsequences_sorted.fa > diff.txt 

num_lines=$(wc -l < "diff.txt")

if [ "$num_lines" -gt 1 ]; then
    echo "Test FAIL"
else
    echo "Test OK"
fi