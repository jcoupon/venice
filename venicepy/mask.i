/* 
 * mask.i
 * The swig interface file to build a python module for venice.
 *
 * Jan 2012 Ben Granett granett@gmail.com
 *
 */


%module (docstring="A wrapper for Venice") venicepy

%feature("autodoc","0");
%feature("autodoc","get_weight(self, ra, dec, n) -> int") get_weight_array;


%rename (_check_point) _check_point_array(double*, int, double*, int, int*, int);
%ignore _check_point(int, double*, double*, int*);

%rename (_random_cat) _random_cat_array(double bounds[4], int inout,  double*, int, double*, int);
%ignore _random_cat(int n, double bounds[4], int inout, double *out_ra, double *out_dec);


%{
#define SWIG_FILE_WITH_INIT
#include "mask.h"
%}

%include "numpy/numpy.i"

%init %{
 	import_array();
%}



%apply (double* INPLACE_ARRAY1, int DIM1) {(double* ra, int nra),
                                           (double* dec, int ndec)};
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* out, int nout)};


%apply (double IN_ARRAY1[ANY]){(double bounds[4])};

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* out_ra, int nra),(double* out_dec, int ndec)};


%exception _check_point_array {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}



%extend Mask {
  void _check_point_array(double* ra, int nra, double* dec, int ndec, int* out, int nout) {
    if (nra != ndec) {
      PyErr_Format(PyExc_ValueError,
		   "Arrays of lengths (%d,%d) given",
		   nra, ndec);
    }
    $self->Mask::_check_point(nra, ra, dec, out);
  }
 
%pythoncode %{

def check_point(self, ra, dec):
    """ ra,dec -> bool 
The convention used in venice is False means INSIDE the mask, True OUTSIDE"""
    return self._check_point(ra,dec,len(ra)) > 0  # return a boolean type

%}


  int _random_cat_array(double bounds[4], int inout, double * out_ra, int nra, double * out_dec, int ndec) {
    return $self->Mask::_random_cat(nra, bounds, inout, out_ra, out_dec);
  }

%pythoncode %{
def random_cat(self, bounds, nrandom, inout='outside'):
    assert(len(bounds)==4)
    assert(bounds[0]<bounds[2])
    assert(bounds[1]<bounds[3])
    n = int(nrandom)

    print type(inout),type('ciao')
    if type(inout)==type('ciao'):
        print inout
        if inout[0]=='i':
            inout = 0
        elif inout[0]=='o':
            inout = 1
        else:
            print "Warning! random_cat: inout parameter should start with 'i' or 'o' or be integer 0 or 1. Using default outside."
            inout = 1
    
    inout = int(inout)
    assert(inout==0 or inout==1)

    n,ra,dec = self._random_cat(bounds,inout,n,n)
    return ra,dec
%}


 }

%include "mask.h"
