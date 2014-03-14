
TOOLS = ngs2dSFS ngsCovar ngsFST ngsStat

CC = g++
CFLAGS = -lm -lz -O3 -Wall

all: $(TOOLS)

$(TOOLS):
	$(CC) $(CFLAGS) $@.cpp -o bin/$@

test:
	@cd examples/; sh ./test.sh 2> /dev/null; cd ../

clean:
	@rm -rf bin/ngsFST bin/ngsCovar bin/ngs2dSFS bin/ngsStat *.o examples/testA*
