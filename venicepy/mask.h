extern "C" {
#include "main.h"
}

class Mask{

 private:
  double x0[2], xmin[2], xmax[2];

  Node *polyTree;
  Polygon *polys;
  int Npolys;

 public:
  Mask(char *mask_file);
  ~Mask();
  void _check_point(int, double*, double*, int*);
  int _random_cat(int n, double bounds[4], int inout, double *out_ra, double *out_dec);

};
