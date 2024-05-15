CC = g++

CFLAGS= -Ofast -std=c++17 -Wall  -w -march=native -fopenmp -flto
CFLAGS_DELTA= -Ofast -std=c++17 -Wall -w -march=native
LDFLAGS = -LTurboPFor-Integer-Compression -l:libic.a -pthread -flto -lpthread -lz -fopenmp

PREXEC = turboPFor_compile
EXEC = k2r_index k2r_query k2r_stat
OBJ = utils.o index_color.o Decycling.o
OBJ_INDEX = k2r_index.o utils.o index_color.o Decycling.o TurboPFor-Integer-Compression/libic.a
OBJ_QUERY = k2r_query.o utils.o index_color.o Decycling.o TurboPFor-Integer-Compression/libic.a
OBJ_STAT = k2r_stat.o utils.o index_color.o Decycling.o  TurboPFor-Integer-Compression/libic.a

all :  $(EXEC)

TurboPFor-Integer-Compression/libic.a:
	make -C TurboPFor-Integer-Compression/

k2r_index : $(OBJ_INDEX)
	$(CC) -o k2r_index $^ $(LDFLAGS)

k2r_query : $(OBJ_QUERY)
	$(CC) -o k2r_query $^ $(LDFLAGS)

k2r_stat : $(OBJ_STAT)
	$(CC) -o k2r_stat $^ $(LDFLAGS)

k2r_index.o: sources/k2r_index.cpp headers/MinimizerLister.h
	$(CC) -o $@ -c $< $(CFLAGS)

k2r_query.o: sources/k2r_query.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

k2r_stat.o: sources/k2r_stat.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: sources/utils.cpp headers/utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

index_color.o: sources/index_color.cpp headers/index_color.h
	$(CC) -o $@ -c $< $(CFLAGS)

Decycling.o: sources/Decycling.cpp headers/Decycling.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf $(EXEC) $(OBJ) $(OBJ_INDEX) $(OBJ_QUERY)*.dat .snakemake .vscode log/memory_res/ log/reads_res/

rebuild: clean $(EXEC)

test:
	./k2r_index -r example/reads/reads.fasta -b example/output/output_binary -s 26
	./k2r_query -f example/sequences/fof.txt -o example/output/query_output -b example/output/output_binary -t 1 -r 0.2
