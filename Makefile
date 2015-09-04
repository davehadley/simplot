
all:
	python setup.py build_ext --inplace

clean:
	-rm simplot/rootplot/cmerge.cpp
	-rm simplot/rootplot/cmerge.so
	-rm simplot/mc/likelihood/clikelihood.cpp
	-rm simplot/mc/likelihood/clikelihood.so
	-rm simplot/sparsehist/*.cpp
	-rm simplot/sparsehist/*.so
