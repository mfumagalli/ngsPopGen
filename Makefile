
CC=g++
CFLAGS=-lm -lz -O3 -Wall

all: ngsCovar ngs2dSFS ngsFST ngsStat

ngsCovar: ngsCovar.cpp
	$(CC) $(CFLAGS) ngsCovar.cpp -o ngsCovar

ngs2dSFS: ngs2dSFS.cpp
	$(CC) $(CFLAGS) ngs2dSFS.cpp -o ngs2dSFS

ngsFST: ngsFST.cpp
	$(CC) ngsFST.cpp -o ngsFST $(CFLAGS)

ngsStat: ngsStat.cpp
	$(CC) ngsStat.cpp -o ngsStat $(CFLAGS)

clean:
	rm -rf *o ngsFST ngsCovar ngs2dSFS ngsStat
