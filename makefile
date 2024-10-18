CC = g++

CFLAGS= -O3 -g -std=c++17 -Wall  -w -fopenmp -flto
CFLAGS_DELTA= -O3 -std=c++17 -Wall -w -flto
LDFLAGS = -LTurboPFor-Integer-Compression -l:libic.a -pthread  -lpthread -lz -fopenmp -flto

PREXEC = turboPFor_compile
EXEC = k2r_index k2r_query test_colors
OBJ = utils.o index_color.o Decycling.o color.o test_colors.o
OBJ_INDEX = k2r_index.o utils.o index_color.o color.o Decycling.o TurboPFor-Integer-Compression/libic.a
OBJ_QUERY = k2r_query.o utils.o index_color.o color.o Decycling.o TurboPFor-Integer-Compression/libic.a
OBJ_TEST_COLOR = color.o test_colors.o

all :  $(EXEC)

TurboPFor-Integer-Compression/libic.a:
	make -C TurboPFor-Integer-Compression/

k2r_index : $(OBJ_INDEX)
	$(CC) -o k2r_index $^ $(LDFLAGS)

k2r_query : $(OBJ_QUERY)
	$(CC) -o k2r_query $^ $(LDFLAGS)

k2r_index.o: sources/k2r_index.cpp headers/MinimizerLister.h
	$(CC) -o $@ -c $< $(CFLAGS)

k2r_query.o: sources/k2r_query.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: sources/utils.cpp headers/utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

index_color.o: sources/index_color.cpp headers/index_color.h
	$(CC) -o $@ -c $< $(CFLAGS)

Decycling.o: sources/Decycling.cpp headers/Decycling.h
	$(CC) -o $@ -c $< $(CFLAGS)

color.o: sources/color.cpp headers/color.h
	$(CC) -o $@ -c $< $(CFLAGS)

test_colors : $(OBJ_TEST_COLOR)
	$(CC) -o test_colors $^ $(LDFLAGS)

test_colors.o: tests/test_colors.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf $(EXEC) $(OBJ) $(OBJ_INDEX) $(OBJ_QUERY)*.dat .snakemake .vscode log/memory_res/ log/reads_res/

rebuild: clean $(EXEC)

SMALL_FILE_TO_INDEX = example/reads/reads1000.fasta #log/celill50k_noN.fa log/bench_coverage_ecoli/ecoli01p_10000_0.01_50.fa log/bench_coverage_cel/cel_10000_0.01_10.fa
FILE_TO_INDEX = log/bench_coverage_ecoli/ecoli01p_10000_0.01_200.fa #log/bench_coverage_cel/cel_10000_0.01_30.fa log/celill500k_noN.fa 
LARGE_TO_INDEX = log/celill50M_noN.fa #log/bench_coverage_ecoli/ecoli01p_10000_0.01_500.fa log/bench_coverage_cel/cel_10000_0.01_70.fa

test_index_mono:
	./k2r_index -r $(FILE_TO_INDEX) -b example/output/output_binary -s 32 --min-ab 1 --max-ab 50000 -t 1

test_index_multi:
	./k2r_index -r $(FILE_TO_INDEX) -b example/output/output_binary -s 32 --min-ab 1 --max-ab 50000 -t 8

test_index_multiVSmono:
	bash test_multi.sh 

test_all_index_mono:
	./k2r_index -r $(SMALL_FILE_TO_INDEX) -b example/output/output_binary -s 32 --min-ab 1 --max-ab 50000 -t 1; \
	head -n 2 $(SMALL_FILE_TO_INDEX) | sed -n '1p; 2s/\(.\{40\}\).*/\1/p' > example/sequences/sequence1.fasta; \
	head -n 4 $(SMALL_FILE_TO_INDEX) | tail -n 2 | sed -n '1p; 2s/\(.\{40\}\).*/\1/p' > example/sequences/sequence2.fasta; \
	./k2r_query -f example/sequences/fof.txt -o out_test -b example/output/output_binary -r 1 -t 1 --match; \
	bash check.sh fof $(SMALL_FILE_TO_INDEX) example/sequences/fof.txt out_test; \

test_all_index_mono_loop:
	for dataset in $(SMALL_FILE_TO_INDEX); do \
		./k2r_index -r $$dataset -b example/output/output_binary -s 32 --min-ab 1 --max-ab 50000 -t 1; \
		head -n 2 $$dataset | sed -n '1p; 2s/\(.\{40\}\).*/\1/p' > example/sequences/sequence1.fasta; \
		head -n 4 $$dataset | tail -n 2 | sed -n '1p; 2s/\(.\{40\}\).*/\1/p' > example/sequences/sequence2.fasta; \
		./k2r_query -f example/sequences/fof.txt -o out_test -b example/output/output_binary -r 1 -t 1 --reads; \
		bash check.sh fof $$dataset example/sequences/fof.txt out_test; \
	done

test_all_index_multi:
	./k2r_index -r $(LARGE_TO_INDEX) -b example/output/output_binary -s 32 --min-ab 1 --max-ab 50000 -t 8
	head -n 2 $(LARGE_TO_INDEX) | sed -n '1p; 2s/\(.\{40\}\).*/\1/p' > example/sequences/sequence1.fasta
	head -n 4 $(LARGE_TO_INDEX) | tail -n 2 | sed -n '1p; 2s/\(.\{40\}\).*/\1/p' > example/sequences/sequence2.fasta
	./k2r_query -f example/sequences/fof.txt -o out_test -b example/output/output_binary -r 1 -t 1 --reads
	bash check.sh fof $(LARGE_TO_INDEX) example/sequences/fof.txt out_test

test_all_index_multi_loop:
	for dataset in $(LARGE_TO_INDEX); do \
		./k2r_index -r $$dataset -b example/output/output_binary -s 32 --min-ab 1 --max-ab 50000 -t 8; \
		head -n 2 $$dataset | sed -n '1p; 2s/\(.\{40\}\).*/\1/p' > example/sequences/sequence1.fasta; \
		head -n 4 $$dataset | tail -n 2 | sed -n '1p; 2s/\(.\{40\}\).*/\1/p' > example/sequences/sequence2.fasta; \
		./k2r_query -f example/sequences/fof.txt -o out_test -b example/output/output_binary -r 1 -t 1 --reads; \
		bash check.sh fof $$dataset example/sequences/fof.txt out_test; \
	done


