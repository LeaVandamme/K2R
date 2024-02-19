# Indexing K-mers to long reads

## Compilation

git clone --recurse-submodules https://github.com/LeaVandamme/index_kmer_to_reads.git

cd index_kmer_to_reads/TurboPFor-Integer-Compression/ ; make

make -j

## Launch 

#### Arguments (Creation) :

  Example : 
<code>
./k2r_index --read-file log/reads/phage10.fa --binary-prefix test_binary --th 8 -k 31 -m 15 -o 1
</code>


#### Arguments (Query) :

  Example : 
<code>
./k2r_query --query-file log/sequences/testPhage.fa --output-prefix output --binary-prefix test_binary -s
</code>


<br/> 
