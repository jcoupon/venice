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

#endif
