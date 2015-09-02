
all:
	python setup.py build_ext --inplace

clean:
	-rm *.cpp
	-rm *.so
	-rm *.pyc
