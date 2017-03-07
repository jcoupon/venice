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

void *readFits(const Config *para, int *bitpix, int *status, long naxes[2], double (**convert)(void *,long ));


double convertCHAR(void *table, long i);
double convertSHORT(void *table, long i);
double convertLONG(void *table, long i);
double convertFLOAT(void *table, long i);
double convertDOUBLE(void *table, long i);

char *getFormatFromImageType_string(int bitpix);
int getFormatFromImageType_number(int bitpix);


#endif
