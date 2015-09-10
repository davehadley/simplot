
all: simplot prob3pp

simplot:
	python setup.py build_ext --inplace

prob3pp:
	cd simplot/rootprob3pp/Prob3++.20121225 && make all

clean:
	-rm simplot/rootplot/cmerge.cpp
	-rm simplot/rootplot/cmerge.so
	-rm simplot/mc/likelihood/clikelihood.cpp
	-rm simplot/mc/likelihood/clikelihood.so
	-rm simplot/sparsehist/*.cpp
	-rm simplot/sparsehist/*.so
