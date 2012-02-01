venicepy

This is a python module for the Venice code by Jean Coupon.  This
wrapper is based on a C++ class with swig to generate the python
interface.  The C++ routines and swig interface file were written by
Ben Granett (granett@gmail.com).

The convention used when testing whether a point falls within the mask
is: the value True (1) indicates outside the mask and False (0) means
inside the mask.

The files are:

mask.cpp  ~ C++ class that wraps Venice functions
mask.i    ~ a swig interface file


*** How to compile:

make swig
this should run something like

	g++ $(CFLAGS) -c mask.cpp -o mask.o
	swig -python -c++ mask.i
	g++ $(CFLAGS) -c mask_wrap.cxx -I/usr/include/python2.7 
	g++ $(CFLAGS) -shared  mask_wrap.o mask.o main.o -o _venicepy.so  $(LDFLAGS)

It is important that main.o is compiled with the -fPIC option.

*** How to use it

See the test script testvenice.py.
Here is an example python session:

>>> import pylab,numpy as N
>>> import venicepy
>>> 
>>> ra = N.random.uniform(35,36,1e4)
>>> dec = N.random.uniform(-5,-4,1e4)
>>> 
>>> M = venicepy.Mask("CFHTLS_W1_i_T0005_Jean.reg")
npoly: 16098
x0: 30.799420, -7.843890
>>> 
>>> inside = M.check_point(ra,dec)
>>> 
>>> print len(ra),len(ra[inside])
10000 9065
