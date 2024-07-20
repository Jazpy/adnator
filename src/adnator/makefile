CXX = g++
read_sim.so: read_simulation.o
	$(CXX) -shared $^ -o $@ -fopenmp

read_simulation.o : CXXFLAGS+=-fPIC -std=c++17 -Wall -Wextra -Werror -fopenmp -O3
