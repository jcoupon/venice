/*
 *    utils.c
 *    venice
 *    Jean Coupon (2012-2017)
 */

#include "fits.h"

void *readFits(const Config *para, int *bitpix, int *status, long naxes[2], double (**convert)(void *,long )){
  /* Reads fits file and return info. "convert" is a pointer to a pointer of a function. */

  if(para->coordType != CART){
    fprintf(stderr,"%s: fits file detected. coord should be set to cart for image coordinates. Exiting...\n",MYNAME);
    exit(EXIT_FAILURE);
  }
  fitsfile *fptr;

  fits_open_file(&fptr, para->fileRegInName, READONLY, status);
  if(*status == FILE_NOT_OPENED){
    fprintf(stderr,"%s: %s not found. Exiting...\n",MYNAME,para->fileRegInName);
    exit(EXIT_FAILURE);
  }
  fits_get_img_type(fptr, bitpix, status);
  fits_get_img_size(fptr, 2, naxes, status);

  long fpixel[2];
  fpixel[0] = fpixel[1] = 1;

  void   *result;
  char   *table_CHAR;
  short  *table_SHORT;
  long   *table_LONG;
  float  *table_FLOAT;
  double *table_DOUBLE;

  fprintf(stderr,"Reading fits file (format detected: ");
  switch (*bitpix){
  case BYTE_IMG:
    fprintf(stderr,"BYTE)...\n");
    table_CHAR = (char *)malloc(naxes[0]*naxes[1]*sizeof(char));
    fits_read_pix(fptr,TBYTE,fpixel,naxes[0]*naxes[1],NULL,table_CHAR,NULL,status);
    result = (void *)table_CHAR;
    *convert = convertCHAR;
    break;
  case SHORT_IMG:
    fprintf(stderr,"SHORT)...\n");
    table_SHORT = (short *)malloc(naxes[0]*naxes[1]*sizeof(short));
    fits_read_pix(fptr,TSHORT,fpixel,naxes[0]*naxes[1],NULL,table_SHORT,NULL,status);
    result = (void *)table_SHORT;
    *convert = convertSHORT;
    break;
  case LONG_IMG:
    fprintf(stderr,"LONG)...\n");
    table_LONG = (long *)malloc(naxes[0]*naxes[1]*sizeof(long));
    fits_read_pix(fptr,TLONG,fpixel,naxes[0]*naxes[1],NULL,table_LONG,NULL,status);
    result = (void *)table_LONG;
    *convert = convertLONG;
    break;
  case FLOAT_IMG:
    fprintf(stderr,"FLOAT)...\n");
    table_FLOAT = (float *)malloc(naxes[0]*naxes[1]*sizeof(float));
    fits_read_pix(fptr,TFLOAT,fpixel,naxes[0]*naxes[1],NULL,table_FLOAT,NULL,status);
    result = (void *)table_FLOAT;
    *convert = convertFLOAT;
    break;
  case DOUBLE_IMG:
    fprintf(stderr,"DOUBLE)...\n");
    table_DOUBLE = (double *)malloc(naxes[0]*naxes[1]*sizeof(double));
    fits_read_pix(fptr,TDOUBLE,fpixel,naxes[0]*naxes[1],NULL,table_DOUBLE,NULL,status);
    result = (void *)table_DOUBLE;
    *convert = convertDOUBLE;
    break;
  default:
    fprintf(stderr,"NULL) \n%s: fits format not recognized. Exiting...\n",MYNAME);
    exit(EXIT_FAILURE);
  }

  fits_close_file(fptr, status);

  if (*status)
    fits_report_error(stderr, *status);

  return result;
}



double convertCHAR(void *table, long i){
	char *result = (char *)table;
	return (double)(result[i]);
}
double convertSHORT(void *table, long i){
	short *result = (short *)table;
	return (double)(result[i]);
}
double convertLONG(void *table, long i){
	long *result = (long *)table;
	return (double)(result[i]);
}
double convertFLOAT(void *table, long i){
	float *result = (float *)table;
	return (double)(result[i]);
}
double convertDOUBLE(void *table, long i){
	double *result = (double *)table;
	return result[i];
}



char *getFormatFromImageType_string(int bitpix){
	/* 	Takes bitpix form image and convert into TFORM value */

	switch (bitpix){
		case BYTE_IMG:
			return "1B";
			break;
		case SHORT_IMG:
			return "1I";
			break;
		case LONG_IMG:
			return "1K";
			break;
		case FLOAT_IMG:
			return "1E";
			break;
		case DOUBLE_IMG:
			return "1D";
			break;
		default:
		fprintf(stderr,"NULL) \n%s: fits format not recognized. Exiting...\n",MYNAME);
		exit(EXIT_FAILURE);
	}

	return "\0";
}



int getFormatFromImageType_number(int bitpix){
	/* 	Takes bitpix form image and convert into TFORM value */

	switch (bitpix){
		case BYTE_IMG:
			return TBYTE;
			break;
		case SHORT_IMG:
			return TSHORT;
			break;
		case LONG_IMG:
			return TLONGLONG;
			break;
		case FLOAT_IMG:
			return TFLOAT;
			break;
		case DOUBLE_IMG:
			return TDOUBLE;
			break;
		default:
		fprintf(stderr,"NULL) \n%s: fits format not recognized. Exiting...\n",MYNAME);
		exit(EXIT_FAILURE);
	}

	return 0;
}
