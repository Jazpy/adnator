CXX = clang++
read_sim.so: read_simulation.o
	$(CXX) -shared $^ -o $@ -ltbb -fopenmp

read_simulation.o : CXXFLAGS+=-fPIC -std=c++17 -Wpedantic -fopenmp -O3
