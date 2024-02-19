CC = g++

CFLAGS= -Ofast -std=c++17 -Wall -g -w -march=native -fopenmp -flto
CFLAGS_DELTA= -Ofast -std=c++17 -Wall -g -w -march=native
LDFLAGS = -LTurboPFor-Integer-Compression -l:libic.a -pthread -flto -lpthread -lz -fopenmp

EXEC = k2r_index k2r_query k2r_stat
OBJ = utils.o fastDelta.o index_color.o Decycling.o
OBJ_INDEX = k2r_index.o utils.o fastDelta.o index_color.o Decycling.o
OBJ_QUERY = k2r_query.o utils.o fastDelta.o index_color.o Decycling.o
OBJ_STAT = k2r_stat.o utils.o fastDelta.o index_color.o Decycling.o

all : $(EXEC)

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

fastDelta.o: include/fastDelta.c
	$(CC) -o $@ -c $< $(CFLAGS_DELTA)

clean:
	rm -rf $(EXEC) $(OBJ) $(OBJ_INDEX) $(OBJ_QUERY)*.dat .snakemake .vscode log/memory_res/ log/reads_res/

rebuild: clean $(EXEC)
