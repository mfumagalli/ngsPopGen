CXX ?= g++

TOOLS = ngs2dSFS ngsCovar ngsFST ngsStat

CXXFLAGS := -O3 -Wall $(CXXFLAGS)
LDLIBS = -lm -lz

all: $(TOOLS)

$(TOOLS): %: %.cpp %.hpp shared.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $@.cpp $(LDFLAGS) $(LDLIBS) -o $@

test:
	@cd examples/; bash test.sh 2> test.log; cd ../

clean:
	@rm -rf $(TOOLS) *.o examples/testA*

.PHONY: all clean test
