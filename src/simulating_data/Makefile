program: simulate_data

CXXFLAGS = -Wall -pedantic -O2 -g -std=c++11

simulate_data: simulate_data.cc functionobjects.o functionobjects.hh felsenstein_pruning.o ../MCMC/felsenstein_pruning.hh mutationclass.o mutationclass.hh back_tree.o back_tree.hh parameters.o parameters.hh
	g++ -o simulate_Bock_data simulate_Bock_data.cc functionobjects.o felsenstein_pruning.o mutationclass.o back_tree.o parameters.o $(CXXFLAGS)

functionobjects.o: functionobjects.cc functionobjects.hh
	g++ -c functionobjects.cc $(CXXFLAGS)

felsenstein_pruning.o: ../MCMC/felsenstein_pruning.cc ../MCMC/felsenstein_pruning.hh
	g++ -c ../MCMC/felsenstein_pruning.cc $(CXXFLAGS)

back_tree.o: back_tree.cc back_tree.hh
	g++ -c back_tree.cc $(CXXFLAGS)

mutationclass.o: mutationclass.cc mutationclass.hh
	g++ -c mutationclass.cc $(CXXFLAGS)

parameters.o: parameters.cc parameters.hh
	g++ -c parameters.cc $(CXXFLAGS)

clean:
	rm ~.o methylation_simulator 
