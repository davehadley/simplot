
.PHONY: all
all: simplot prob3pp doc

.PHONY: simplot
simplot: prob3pp
	python setup.py build_ext --inplace

.PHONY: prob3pp
prob3pp:
	cd simplot/rootprob3pp/Prob3++.20121225 && make all

.PHONY: clean
clean:
	-rm simplot/rootplot/cmerge.cpp
	-rm simplot/rootplot/cmerge.so
	-rm simplot/mc/likelihood/clikelihood.cpp
	-rm simplot/mc/likelihood/clikelihood.so
	-rm simplot/sparsehist/*.cpp
	-rm simplot/sparsehist/*.so
	-rm simplot/binnedmodel/*.cpp
	-rm simplot/binnedmodel/*.so

.PHONY: doc
doc: simplot
	mkdir -p doc && cd doc && epydoc simplot && echo "simplot documentation: doc/html/index.html"
