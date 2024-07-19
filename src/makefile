read_sim.so: read_simulation.o
	$(CXX) -shared $^ -o $@ -Wl,--no-whole-archive

read_simulation.o : CXXFLAGS+=-fPIC
