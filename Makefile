CXX := g++
#CXXFLAGS := -O0 -g -ggdb -fsanitize=address
CXXFLAGS := -O3 -march=native -g
PROGRAMS := mcsplit-si mcsplit-si-ll mcsplit-si-adjmat

all: $(PROGRAMS)

mcsplit-si: mcsplit-si.c sparse_graph.c sparse_graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o $@ $< sparse_graph.c -pthread

mcsplit-si-ll: mcsplit-si-ll.c sparse_graph.c sparse_graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o $@ $< sparse_graph.c -pthread

mcsplit-si-adjmat: mcsplit-si-adjmat.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o $@ $< graph.c -pthread

clean:
	rm $(PROGRAMS)
