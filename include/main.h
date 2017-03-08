/*
 *    main.h
 *    venice
 *    Jean Coupon (2012-2017)
 */

#ifndef MAIN_H
#define MAIN_H

#include "utils.h"
#include "init.h"
#include "fits.h"

void testPython();

int mask2d(const Config *para);
int flagCatFits(const Config *para);
int flagCat(const Config *para);
int randomCat(const Config *para);


#endif
