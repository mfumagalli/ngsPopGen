
TOOLS = ngs2dSFS ngsCovar ngsFST ngsStat

CC = g++
CFLAGS = -lm -lz -O3 -Wall

all: $(TOOLS)

$(TOOLS):
	$(CC) $(CFLAGS) $@.cpp -o $@

test:
	@cd examples/; bash test.sh 2> test.log; cd ../

clean:
	@rm -rf ngsFST ngsCovar ngs2dSFS ngsStat *.o examples/testA*
