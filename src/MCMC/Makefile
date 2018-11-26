program: IWEinference

CXXFLAGS = -Wall -pedantic -O2 -g -std=c++11

IWEinference: IWEinference.cc functionobjects.o functionobjects.hh felsenstein_pruning.o felsenstein_pruning.hh MCMC_step.o MCMC_step.hh mutationclass.o mutationclass.hh back_tree.o back_tree.hh parameters.o parameters.hh
	g++ -o IWEinference IWEinference.cc functionobjects.o felsenstein_pruning.o MCMC_step.o mutationclass.o back_tree.o parameters.o $(CXXFLAGS)

functionobjects.o: functionobjects.cc functionobjects.hh
	g++ -c functionobjects.cc $(CXXFLAGS)

MCMC_step.o: MCMC_step.cc MCMC_step.hh
	g++ -c MCMC_step.cc $(CXXFLAGS)

felsenstein_pruning.o: felsenstein_pruning.cc felsenstein_pruning.hh
	g++ -c felsenstein_pruning.cc $(CXXFLAGS)

back_tree.o: back_tree.cc back_tree.hh
	g++ -c back_tree.cc $(CXXFLAGS)

mutationclass.o: mutationclass.cc mutationclass.hh
	g++ -c mutationclass.cc $(CXXFLAGS)

parameters.o: parameters.cc parameters.hh
	g++ -c parameters.cc $(CXXFLAGS)

clean:
	rm ~.o IWEinference 
