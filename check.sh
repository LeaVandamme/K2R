#!/bin/bash


##################################################################################################
if [ "$1" = "seq" ]; then

    echo "TEST QUERY 1 FASTA"
    rm check_query/*

    # ON PREND LES SEQ À QUERY
    awk 'BEGIN {RS=">"; ORS=""} NR > 1 {print ">"$0}' $3 | grep -v "^>" > check_query/target_sequences.txt

    # NUMEROTATION READS
    awk '/^>/ {print ">" i++; next} {print}' $2 > check_query/numeroted.fasta

    # RECHERCHE DES READS

    # Extraire les séquences correspondant aux identifiants avec leurs positions
    line=$(head -n 1 check_query/target_sequences.txt)
    seqkit grep -P -s -p $line check_query/numeroted.fasta -o check_query/selected.fasta
    seqtk seq -l 0 check_query/selected.fasta > check_query/selected_sequences.fasta


    # TRI DES DEUX FICHIERS

    awk '/^>/ {if (seq) print seq; print ""; print $0; seq=""; next} {seq=seq $0} END {if (seq) print seq}' $4 | awk 'BEGIN {RS=""; FS="\n"} {header=$1; seq=""; for (i=2; i<=NF; i++) seq=seq $i; print header, seq}' | sort -k1,1 | awk '{header=$1; seq=$2; print header; print seq}' > check_query/queryresults_sorted.fa
    awk '/^>/ {if (seq) print seq; print ""; print $0; seq=""; next} {seq=seq $0} END {if (seq) print seq}' check_query/selected_sequences.fasta | awk 'BEGIN {RS=""; FS="\n"} {header=$1; seq=""; for (i=2; i<=NF; i++) seq=seq $i; print header, seq}' | sort -k1,1 | awk '{header=$1; seq=$2; print header; print seq}' > check_query/selectedsequences_sorted.fa

    # DIFF

    diff -s -B  check_query/queryresults_sorted.fa check_query/selectedsequences_sorted.fa > check_query/diff_1seq.txt 

    num_lines=$(wc -l < "check_query/diff_1seq.txt")

    if [ "$num_lines" -gt 1 ]; then
        echo "Test FAIL"
    else
        echo "Test OK"
    fi


##################################################################################################
else
    echo "TEST QUERY FILE OF FILE"
    rm check_query_fof/*.txt check_query_fof/*.fasta check_query_fof/*.fa
    
    # ON PREND LES SEQ À QUERY
    while IFS= read line; do
        filename="${line##*/}"
        basename="${filename%%.*}"
        awk 'BEGIN {RS=">"; ORS=""} NR > 1 {print ">"$0}' $line | grep -v "^>" > check_query_fof/target_sequences_$basename.txt
    done < $3
    

    # NUMEROTATION READS
    awk '/^>/ {print ">" i++; next} {print}' $2 > check_query_fof/numeroted.fasta

    # RECHERCHE DES READS

    # Extraire les séquences correspondant aux identifiants avec leurs positions
    while IFS= read -r line; do
        filename="${line##*/}"
        basename="${filename%%.*}"
        line_read=$(head -n 1 check_query_fof/target_sequences_$basename.txt)
        seqkit grep -P -s -p $line_read check_query_fof/numeroted.fasta -o check_query_fof/selected_$basename.fasta
        seqtk seq -l 0 check_query_fof/selected_$basename.fasta > check_query_fof/selected_sequences_$basename.fasta
    done < $3


    # TRI DES FICHIERS

    while IFS= read -r line; do
        filename="${line##*/}"
        basename="${filename%%.*}"
        awk '/^>/ {if (seq) print seq; print ""; print $0; seq=""; next} {seq=seq $0} END {if (seq) print seq}' $4_$basename.fasta | awk 'BEGIN {RS=""; FS="\n"} {header=$1; seq=""; for (i=2; i<=NF; i++) seq=seq $i; print header, seq}' | sort -k1,1 | awk '{header=$1; seq=$2; print header; print seq}' > check_query_fof/queryresults_sorted_$basename.fa
        awk '/^>/ {if (seq) print seq; print ""; print $0; seq=""; next} {seq=seq $0} END {if (seq) print seq}' check_query_fof/selected_sequences_$basename.fasta | awk 'BEGIN {RS=""; FS="\n"} {header=$1; seq=""; for (i=2; i<=NF; i++) seq=seq $i; print header, seq}' | sort -k1,1 | awk '{header=$1; seq=$2; print header; print seq}' > check_query_fof/selectedsequences_sorted_$basename.fa
    done < $3

    # DIFF
    
    while IFS= read -r line; do
        filename="${line##*/}"
        basename="${filename%%.*}"
        diff -s -B  check_query_fof/queryresults_sorted_$basename.fa check_query_fof/selectedsequences_sorted_$basename.fa > check_query_fof/diff_$basename.txt
        num_lines=$(wc -l < "check_query_fof/diff_$basename.txt")
        if [ "$num_lines" -gt 1 ]; then
            echo "Test FAIL : $basename"
        else
            echo "Test OK : $basename"
        fi
    done < $3
fi