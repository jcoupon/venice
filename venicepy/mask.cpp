/* 
 * mask.cpp
 * a C++ wrapper for the venice code by Jean Coupon.  
 *
 * Jan 2012 Ben Granett granett@gmail.com
 *
 */

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "mask.h"


extern "C" {
#include "main.h"
}

using namespace std;


Mask::Mask(char *reg_file){

  FILE *fileRegIn = fopenAndCheck(reg_file,"r");

  //Reads the mask file (.reg DS9 file)
  //and construct the array of polygons.
  polyTree = readPolygonFileTree(fileRegIn,xmin,xmax);
  polys = (Polygon *)polyTree->polysAll;
  Npolys = polyTree->Npolys;

  //Reference point. It must be outside the mask;
  x0[0] = xmin[0] - 1.0; x0[1] = xmin[1] - 1.0;

  fprintf(stderr,"npoly: %d\n",Npolys);
  fprintf(stderr,"x0: %f, %f\n",x0[0], x0[1]);

}

Mask::~Mask(){
}

void Mask::_check_point(int n, double *ra, double *dec, int *out_flag){
  int i, flag, poly_id;
  double x[2];

  for (i=0; i<n; i++){

    x[0] = ra[i];
    x[1] = dec[i];

    flag = 0;
    //if (insidePolygonTree(polyTree,x0,x,&poly_id)) flag=1;
    if(flag=0, !insidePolygonTree(polyTree,x0,x,&poly_id)) flag = 1;

    out_flag[i] = flag;
  }
}

int Mask::_random_cat(int n, double bounds[4], int inout, double *out_ra, double *out_dec){
  int Npolys,poly_id,flag;
  int i, cati=0;
  double  x[2], xmin[2], xmax[2];
  gsl_rng *r = randomInitialize(1);

  xmin[0] = bounds[0];
  xmin[1] = bounds[1];
  xmax[0] = bounds[2];
  xmax[1] = bounds[3];

  while (cati < n){
    for (i=0;i<n;i++){
      if (cati >= n) break;
      x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
      x[1] = gsl_ran_flat(r,xmin[1],xmax[1]);

      x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
      x[1] = gsl_ran_flat(r,sin(xmin[1]*PI/180.0),sin(xmax[1]*PI/180.0));
      x[1] = asin(x[1])*180.0/PI;
    
      if(flag=0,!insidePolygonTree(polyTree,x0,x,&poly_id)) flag = 1;
    
      if (flag==inout){
	out_ra[cati] = x[0];
	out_dec[cati] = x[1];
	cati++;
      }
    }
  }


  return 0;
}
