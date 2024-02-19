# K2R : index k-mer to reads

K2R is a tool capable of indexing a set of read data (FASTA) by associating each k-mer with the reads in which they appear.

## Compilation

git clone --recurse-submodules https://github.com/LeaVandamme/K2R_.git

cd index_kmer_to_reads/ ;

cd TurboPFor-Integer-Compression/ ; make ;

cd .. ; make -j

## Build an index

K2R can index a set of reads from a FASTA file.
Several options and parameters are available to adapt the creation of the index to your needs :

```
./k2r_index --read-file read_file.fasta --binary-prefix path_to_binary/binary_prefix --th nb_threads [OPTIONS]

--keep-unique : Keep all k-mers included unique ones (default : False)
--homocompression : Remove repetitions from sequences (example : AGGTTATG -> AGTATG; default : False)
-k : K-mer length (default : 31)
-m : Minimizer length (default : 15)
-b : Counting bloom filter size (log(2), default : 32)
-o : Minimum of occurence of each minimizer (default : 2)
```

### Output files

K2R creates 2 binary files for each index : 

- [binary_prefix]_mmer.bin, which contains the association between each minimizer and its color identifier.
- [binary_prefix]_color.bin, which contains the association between each color identifier and its color.


## Launch queries

K2R can launch several queries from a file of file (fof), containing path to fasta files.
Several options and parameters are available to adapt the queries :

```
./k2r_query --query-file /path_to_fof/fof.txt --output-prefix prefix --binary-prefix path_to_index/index_prefix --thread nb_threads -f

--query_file : File of file containing path to sequence Fasta file
--output-prefix : prefix of output read file (default : query_output)
--binary-prefix : Path and prefix of the index (default : binary_index)
--thread : Number of threads (default : 1)
-t : Rate of minimizer found in the read to keep it in results (between 0 and 1, default: 1)
-s : query a unique sequence
-f : query several sequences from a file of file
```

### Output files

K2R creates 1 file for each query. For example if the queries are launched on a file of file containing 100 paths, 100 files will be created containing the reads in which the sequences appear.

## Example

Here is an example of how to use K2R, using the data provided in the example folder. The dataset contains reads (HiFi) from the E.Coli genome, with a 10X coverage.

We choose to reduce the filter size to 2**26, as the dataset is sufficiently small.


```
./k2r_index --read-file example/reads/reads.fasta --binary-prefix example/output/output_binary -b 26

```

Once the index has been created, the sequences can be queried using the following command : 
  
```
./k2r_query --query-file example/sequences/fof.txt --output-prefix example/output/query_output --binary-prefix example/output/output_binary --th 1 -t 0.2 -f

``` 
