/*
 *    fits.h
 *    venice
 *    Jean Coupon (2012-2017)
 */

#include "utils.h"
#include "fitsio.h"

#ifndef FITS_H
#define FITS_H


/*
 *    FITS
 */

 void *readFits(const Config *para, int *bitpix, int *typecode, char *tform, int *status, long naxes[2], double (**toDouble)(void *,long ), int *size);

double toDoubleCHAR(void *table, long i);
double toDoubleSHORT(void *table, long i);
double toDoubleLONG(void *table, long i);
double toDoubleFLOAT(void *table, long i);
double toDoubleDOUBLE(void *table, long i);

void readColFits(fitsfile *fileIn, int id_num, long N, double *x);

#endif
