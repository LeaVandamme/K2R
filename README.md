# K2R : index k-mer to reads

K2R is a tool capable of indexing a set of read data (FASTA) by associating each k-mer with the reads in which they appear.

## Compilation

```
git clone --recurse-submodules https://github.com/LeaVandamme/K2R.git
cd K2R/ ;
make -j
```

## Build an index

K2R can index a set of reads from a FASTA file.

```
./k2r_index -r example/reads/reads.fa -b example/output/output_binary
```

Or as a complete example :
```
./k2r_index -r example/reads/reads.fa -b example/output/output_binary -t 5 -k 31 -m 15 -s 26 --min-ab 2 --max-ab 2000
```

- reads.fa is a simulated dataset created from the E.Coli genome having a 10X coverage, an error rate of 1% and a read length of 10.000
- 5 threads are used
- k-mer length : 31
- minimizer length : 15
- counting blooom filter size : 2^26 (sufficient for the small dataset being used, this size must be superior to the number of minimizers to index, to avoid a full structure)
- each minimizer have to be seen at least 2 times and at most 2000 times to be kept


Here is a description of the options to adapt the indexation :

```
Arguments :

-r  :  Build index from file (FASTA)
-b  :  Write index in binary files (default : binary_index)

Options :

-t           :  Number of threads used (default: 1, maximum: 5)
-k           :  k-mer size (default : 31)
-m           :  Minimizer size (default : 15, advanced parameter)
-s           :  Counting bloom filter size (log(2), maximum: 32, default: 32). The size must be superior to the number of minimizers to index, to avoid a full structure.
-h           :  Homocompression of reads (delete repeted nucleotides in reads, ex : AAAGGTAAAACC becomes AGTAC, advanced parameter)
--min-ab     :  Minimizers minimum abundance (default: 2)
--max-ab     :  Minimizers maximum abundance (default: 1000)
```


## Launch queries

K2R can launch several queries from a file of file (fof), containing paths to fasta files. In our example, fof.txt contains 5 sequences to query :

```
example/sequences/sequence1.fasta
example/sequences/sequence2.fasta
example/sequences/sequence3.fasta
example/sequences/sequence4.fasta
example/sequences/sequence5.fasta
```

and the minimal command line can be :

```
./k2r_query -f example/sequences/fof.txt -o example/output/query_fof_output -b example/output/output_binary --match
```

But you can also launch a single request (1 sequence) from a FASTA file :

```
./k2r_query -s example/sequences/sequence1.fasta -o example/output/query_seq_output -b example/output/output_binary --match
```

As complete examples :
```
./k2r_query -f example/sequences/fof.txt -o example/output/query_fof_output -b example/output/output_binary --match -r 0.6 -t 1
./k2r_query -s example/sequences/sequence1.fasta -o example/output/query_seq_output -b example/output/output_binary --match -r 0.6 -t 1
```

- the same parameters are used for the two types of queries
- it only returns the number of matches for each sequence
- rate of minimizers found in the read to be kept : 0.6 
- 1 thread is used



Here is a description of the options to adapt the queries :

```
Arguments :

-s OR -f  :     Sequence file (FASTA) if a unique sequence is queried (-s), file of file if several sequences are queried (-f)
--reads OR --match : Output format (default: match, --reads create a fasta file for each query sequence, containing the matching reads (can be slower for a high number of sequences) ; --match only report the number of matches for each sequence, in a single file)

Option :

-o        :     Write output reads in fasta file (default : query_output)
-b        :     Index binary files prefix (default: binary_index)
-r        :     Rate of minimizer found in the read to keep it in results (between 0 and 1, default : 0.4)
-t        :     Number of threads used (default: 1)
```

## Complete test

You can test the different complete indexations and queries directly by using this line : 

```
make test
```

All the results will be available in example/output/.

## Output files

K2R creates 2 binary files for each index :

- [binary_prefix]_mmer.bin, which contains the association between each minimizer and its color identifier.
- [binary_prefix]_color.bin, which contains the association between each color identifier and its color.


Two options are available for queries : 

- K2R creates 1 file for each query (for example if the queries are launched on a file of file containing 100 paths, 100 files will be created containing the reads in which the sequences appear)
- K2R only creates a unique file, containing the number of matches (number of reads in which the sequences are seen). This file should look like this (the file can be unsorted due to multithreading): 

```
example/sequences/sequence1.fasta : 5
example/sequences/sequence2.fasta : 9
example/sequences/sequence3.fasta : 0
example/sequences/sequence4.fasta : 9
example/sequences/sequence5.fasta : 4
```

## Citation

To find out more about the structure, a preprint is available here :

> Lea Vandamme, Bastien Cazaux, Antoine Limasset : Tinted de Bruijn Graphs for efficient read extraction from sequencing datasets (<https://doi.org/10.1101/2024.02.15.580442>)

Contact : lea.vandamme@univ-lille.fr
