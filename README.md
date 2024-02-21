# K2R : index k-mer to reads

K2R is a tool capable of indexing a set of read data (FASTA) by associating each k-mer with the reads in which they appear.

## Compilation

```
git clone --recurse-submodules https://github.com/LeaVandamme/K2R.git
cd K2R/ ;
cd TurboPFor-Integer-Compression/ ; make ;
cd .. ; make -j
```

## Build an index

K2R can index a set of reads from a FASTA file.
Several options and parameters are available to adapt the creation of the index to your needs :

```
./k2r_index -r read_file.fasta -b path_to_binary/binary_prefix -t nb_threads [OPTIONS]
```

```
Arguments : 

-r  :  Build index from file (FASTA allowed)
-b  :  Write index in binary files (default : binary_index)
-t  :  Number of threads used (default: 1)

Options : 

-k  :  k-mer size (default : 31)
-m  :  Minimizer size (default : 15)
-s  :  Counting bloom filter size (log(2), default 32)
-a  :  Minimum occurence of minimizers (default : 2)
-h  :  Homocompression of reads
```


## Launch queries

K2R can launch several queries from a file of file (fof), containing paths to fasta files.

```
./k2r_query -f /path_to_fof/fof.txt -o output_prefix -b path_to_index/index_prefix -t nb_threads
```

But can also launch a single request from a FASTA file : 

```
./k2r_query -s /path_to_fof/seq.fasta -o output_prefix -b path_to_index/index_prefix -t nb_threads
```

Several options and parameters are available to adapt the queries :

```
Arguments : 

-s OR -f  :     Sequence file (FASTA) if a unique sequence is queried (-s), file of file if several sequences are queried (-f)
-o        :     Write output reads in fasta file (default : query_output)
-b        :     Index binary files prefix (default: binary_index)
-t        :     Number of threads used (default: 1)

Option : 

-r        :     Rate of minimizer found in the read to keep it in results (between 0 and 1, default : 0.4)
```

### Output files

K2R creates 2 binary files for each index : 

- [binary_prefix]_mmer.bin, which contains the association between each minimizer and its color identifier.
- [binary_prefix]_color.bin, which contains the association between each color identifier and its color.

K2R creates 1 file for each query. 

For example if the queries are launched on a file of file containing 100 paths, 100 files will be created containing the reads in which the sequences appear.

## Example

Here is an example of how to use K2R, using the data provided in the example folder. The dataset contains reads (HiFi) from the E.Coli genome, with a 10X coverage.

We choose to reduce the filter size to 2^26, as the dataset is sufficiently small.


```
./k2r_index --read-file example/reads/reads.fasta --binary-prefix example/output/output_binary -b 26

```

Once the index has been created, the sequences can be queried using the following command : 
  
```
./k2r_query --query-file example/sequences/fof.txt --output-prefix example/output/query_output --binary-prefix example/output/output_binary --th 1 -t 0.2 -f

```

## Citation

To find out more about the structure, a preprint is available here :

> Lea Vandamme, Bastien Cazaux, Antoine Limasset : Tinted de Bruijn Graphs for efficient read extraction from sequencing datasets (<https://doi.org/10.1101/2024.02.15.580442>)


