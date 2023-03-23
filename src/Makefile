
# add -g for use with gdb
CPPFLAGS=-std=c++17 -O3 -march=native -Wall -Wno-sign-compare

LIBS=-lboost_timer
#LIBS= -lboost_timer -pthread

# uncomment next line to output timing information 
# (by default only number of queue operations is output)
# DFLAGS=-DOUTPUT_TIMING

# specify which data structure is used for priority queue 
#HEAP=-DBINOMIAL
HEAP=-DFIBONACCI


all: random testset fortunate


graph.o: graph.cpp bucket_list.hpp mcpu_timer.hpp
	g++ $(CPPFLAGS) $(HEAP) -c graph.cpp

nn.o: nn.hpp nn.cpp
	g++ $(CPPFLAGS) -c nn.cpp

random_graphs.o: random_graphs.cpp graph.hpp bucket_list.hpp mcpu_timer.hpp
	g++ $(CPPFLAGS) $(DFLAGS) -c random_graphs.cpp

test_set_graphs.o: test_set_graphs.cpp graph.hpp bucket_list.hpp mcpu_timer.hpp
	g++ $(CPPFLAGS) $(DFLAGS) -c test_set_graphs.cpp

fortunate_graphs.o: fortunate_graphs.cpp graph.hpp bucket_list.hpp mcpu_timer.hpp
	g++ $(CPPFLAGS) $(DFLAGS) -c fortunate_graphs.cpp 


random: random_graphs.o graph.o nn.o
	g++ -o random $(CPPFLAGS) $^ $(LIBS)

testset: test_set_graphs.o graph.o nn.o
	g++ -o testset $(CPPFLAGS) $^ $(LIBS)

fortunate: fortunate_graphs.o graph.o nn.o
	g++ -o fortunate $(CPPFLAGS) $^ $(LIBS)


clean:
	rm *.o 

clean-all:
	rm *.o random testset fortunate


