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

test:
	./k2r_index -r example/reads/reads.fa -b example/output/output_binary -t 5 -k 31 -m 15 -s 26 --min-ab 2 --max-ab 2000;
	./k2r_query -f example/sequences/fof.txt -o example/output/query_fof_output -b example/output/output_binary --match -r 0.6 -t 1;
	./k2r_query -s example/sequences/sequence1.fasta -o example/output/query_seq_output -b example/output/output_binary --match -r 0.6 -t 1;
	