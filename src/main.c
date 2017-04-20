/*
 *
 *    main.c
 *    venice
 *    Jean Coupon (2012-2017)
 *
 *    Program that reads a mask file (DS9 type or fits mask) and a catalogue
 *    of objects and computes one of the following tasks:
 *    1. Creates a pixelized mask.
 *    venice -m mask.reg [OPTIONS]
 *    2. Finds objects inside/outside a mask.
 *    venice -m mask.reg -cat file.cat [OPTIONS]
 *    3. Generates a random catalogue of objects inside/outside a mask.
 *    venice -m mask.reg -r [OPTIONS]
 *
 *    TODO:
 *    - adapt to be used with python
 *    - use RA DEC as input even with fits masks (use wcs functions)
 *    - correct "Definied" typos
 *    - allow different input and output files
 */

#include "main.h"

void testPython(){
   /*    test for python */
   fprintf(stderr,"Hello world\n");
}

int main(int argc, char **argv)
{
   /*    initialization */
   srand((unsigned int)time(NULL));
   EPS   = determineMachineEpsilon();
   IDERR = determineSize_tError();
   Config para;

   /*    tasks */
   switch (readParameters(argc,argv,&para)){
      case 1:
         mask2d(&para);     /* binary mask for visualization */
         break;
      case 2:
         if(para.oFileType == FITS){
            flagCatFits(&para);    /* objects in/out of mask */
         }else{
            flagCat(&para);    /* objects in/out of mask */
         }
         break;
      case 3:
         randomCat(&para);  /* random catalogue */
         break;
   }
   return EXIT_SUCCESS;
}

/*
 *    Main routines
 */

int  mask2d(const Config *para){
   /*     Returns the mask in gsl histogram format and writes the mask in fileOut.
    *     The limits are the extrema of the extreme polygons in fileRegIn.
    *     The pixel is set to 0 when inside the mask and 1 otherwise.
    *     For fits format, it writes the pixel value.
    */

   int Npolys,poly_id,flag;
   size_t i, j, count, total;
   double x[2], x0[2], xmin[2], xmax[2];

   FILE *fileOut = fopenAndCheck(para->fileOutName,"w");

   gsl_histogram2d *mask = gsl_histogram2d_alloc(para->nx,para->ny);
   total = para->nx*para->ny;

   if(checkFileExt(para->fileRegInName,".fits")){           /*     fits file */
      long fpixel[2], naxes[2];
      int bitpix, status = 0;
      double (*toDouble)(void *,long ) = NULL;

      /*   read fits file and put in table */
      void *table = readFits(para,&bitpix,NULL,NULL,&status,naxes,&toDouble,NULL);

      /* define limits */
      xmin[0] = xmin[1] = 1.0;
      xmax[0] = naxes[0];
      xmax[1] = naxes[1];
      if(para->minDefinied[0]) xmin[0] = para->min[0];
      if(para->maxDefinied[0]) xmax[0] = para->max[0];
      if(para->minDefinied[1]) xmin[1] = para->min[1];
      if(para->maxDefinied[1]) xmax[1] = para->max[1];
      /*    print out limits */
      fprintf(stderr,"limits:\n");
      fprintf(stderr,"-xmin %g -xmax %g -ymin %g -ymax %g\n",xmin[0],xmax[0],xmin[1],xmax[1]);

      gsl_histogram2d_set_ranges_uniform(mask,xmin[0],xmax[0],xmin[1],xmax[1]);

      /*    ATTENTION "toDouble" converts everything into double */
      fprintf(stderr,"\nProgress =     ");
      for(i=0; i<mask->nx; i++){
         for(j=0; j<mask->ny; j++){
            count = i*para->ny+j;
            printCount(&count,&total,1000);
            x[0] = (mask->xrange[i]+mask->xrange[i+1])/2.0;  /* center of the pixel */
            x[1] = (mask->yrange[j]+mask->yrange[j+1])/2.0;
            fpixel[0] = roundToNi(x[0]) - 1;
            fpixel[1] = roundToNi(x[1]) - 1;
            mask->bin[i*mask->ny+j] = toDouble(table,fpixel[1]*naxes[0]+fpixel[0]);
         }
      }
      fprintf(stderr,"\b\b\b\b100%%\n");

      free(table);

   }else if(checkFileExt(para->fileRegInName,".reg")){          /* ds9 file */

      FILE *fileRegIn = fopenAndCheck(para->fileRegInName,"r");
      Node *polyTree  = readPolygonFileTree(fileRegIn,xmin,xmax);
      Polygon *polys  = (Polygon *)polyTree->polysAll;
      Npolys          = polyTree->Npolys;
      fclose(fileRegIn);

      /*    define limits */
      if(para->minDefinied[0]) xmin[0] = para->min[0];
      if(para->maxDefinied[0]) xmax[0] = para->max[0];
      if(para->minDefinied[1]) xmin[1] = para->min[1];
      if(para->maxDefinied[1]) xmax[1] = para->max[1];
      /*    print out limits */
      fprintf(stderr,"limits:\n");
      fprintf(stderr,"-xmin %g -xmax %g -ymin %g -ymax %g\n",xmin[0],xmax[0],xmin[1],xmax[1]);

      /*    reference point. It must be outside the mask */
      x0[0] = xmin[0] - 1.0; x0[1] = xmin[1] - 1.0;

      gsl_histogram2d_set_ranges_uniform(mask,xmin[0],xmax[0],xmin[1],xmax[1]);

      total = para->nx*para->ny;

      fprintf(stderr,"Progress =     ");
      for(i=0; i<mask->nx; i++){
         for(j=0; j<mask->ny; j++){
            count = i*para->ny+j;
            printCount(&count,&total,1000);
            x[0] = (mask->xrange[i]+mask->xrange[i+1])/2.0; /* center of the pixel */
            x[1] = (mask->yrange[j]+mask->yrange[j+1])/2.0;
            /* 1 = outside the mask, 0 = inside the mask */
            if(!insidePolygonTree(polyTree,x0,x,&poly_id)) mask->bin[i*mask->ny+j] = 1.0;
         }
      }
      fflush(stdout);
      fprintf(stderr,"\b\b\b\b100%%\n");
   }else{
      fprintf(stderr,"%s: mask file format not recognized. Please provide .reg, .fits or no mask but input limits. Exiting...\n",MYNAME);
      exit(EXIT_FAILURE);
   }

   /*    write the output file */
   for(j=0; j<mask->ny; j++){
      for(i=0; i<mask->nx; i++){
         fprintf(fileOut,"%g ",mask->bin[i*mask->ny+j]);
      }
      fprintf(fileOut,"\n");
   }
   fclose(fileOut);

   return EXIT_SUCCESS;
}

int flagCatFits(const Config *para){
   /*
    *    Input fits catalogue version
    *
    *    Reads fileCatIn and add a flag at the end of the line. 1 is outside
    *    the mask and 0 is inside the mask. xcol and ycol are the column ids
    *    of resp. x coordinate and y coordinate.
    *    For mask fits format, it writes the pixel value.
    */


   int Npolys, poly_id, flag, verbose = 1, size, firstelem=1, firstrow=1;
   double x[2], x0[2], xmin[2], xmax[2];
   size_t i, Ncol;
   long N;

	fitsfile *fileCatIn;
	int status = 0, datatype, id_num[2];

	fits_open_table(&fileCatIn, para->fileCatInName, READONLY, &status);
	if (status) {
      fits_report_error(stderr, status);
      exit(EXIT_FAILURE);
   }
   fits_get_num_rows(fileCatIn, &N, &status);
	if (status) {
      fits_report_error(stderr, status);
      exit(EXIT_FAILURE);
   }
   fprintf(stderr,"Nobjects = %zd\n", N);

   int xcol = atoi(para->xcol);
   int ycol = atoi(para->ycol);
	if(xcol == 0){ /* 	if input column name is a string it will return "0" */
			fits_get_colnum(fileCatIn, CASEINSEN, para->xcol, &(xcol), &status);
			if (status){
            fits_report_error(stderr, status);
            exit(EXIT_FAILURE);
         }
   }
	if(ycol == 0){ /* 	if input column name is a string it will return "0" */
			fits_get_colnum(fileCatIn, CASEINSEN, para->ycol, &(ycol), &status);
			if (status) {
            fits_report_error(stderr, status);
            exit(EXIT_FAILURE);
         }
	}

   fitsfile *fileOutFits;     /* pointer to the FITS file, defined in fitsio.h */

   fits_create_file(&fileOutFits, para->fileOutName, &status);
	if (status){
      fits_report_error(stderr, status);
      if (status){
         fprintf(stderr, "Add \"!\" in front of the file name to overwrite: -o \"!FILEOUT\"\n");
      }
      exit(EXIT_FAILURE);
   }
   fits_copy_file(fileCatIn, fileOutFits, 1, 1, 1, &status);
	if (status){
      fits_report_error(stderr, status);
      exit(EXIT_FAILURE);
   }
   int ncols;

   fits_get_num_cols(fileOutFits, &ncols, &status);
	if (status){
      fits_report_error(stderr, status);
      exit(EXIT_FAILURE);
   }

   if(checkFileExt(para->fileRegInName,".fits")){
      if(para->coordType != CART){
         fprintf(stderr,"%s: fits file detected. coord should be set to cart for image coordinates. Exiting...\n",MYNAME);
         exit(EXIT_FAILURE);
      }

      long fpixel[2], naxes[2];
      int bitpix, typecode;
      char tform_char;
      double null = -99.0;
      double (*toDouble)(void *,long ) = NULL;

      /*    read fits file and put in table */
      void *table = readFits(para,&bitpix,&typecode,&tform_char,&status,naxes,&toDouble,&size);

      /*    define limits */
      xmin[0] = xmin[1] = 0.5;
      xmax[0] = naxes[0]+0.5;
      xmax[1] = naxes[1]+0.5;
      if(para->minDefinied[0]) xmin[0] = para->min[0];
      if(para->maxDefinied[0]) xmax[0] = para->max[0];
      if(para->minDefinied[1]) xmin[1] = para->min[1];
      if(para->maxDefinied[1]) xmax[1] = para->max[1];
      /*    print out limits */
      fprintf(stderr,"Mask limits:\n");
      fprintf(stderr,"-xmin %g -xmax %g -ymin %g -ymax %g\n",xmin[0],xmax[0],xmin[1],xmax[1]);

      fits_insert_col(fileOutFits, ncols+1, para->flagName,  concat("1", &tform_char), &status);
	   if (status){
         fits_report_error(stderr, status);
         exit(EXIT_FAILURE);
      }

      double *xx = (double *)malloc(N*sizeof(double));
      double *yy = (double *)malloc(N*sizeof(double));

      readColFits(fileCatIn, xcol, N, xx);
      readColFits(fileCatIn, ycol, N, yy);

      size_t N_size_t = (size_t)N;

      if(verbose) fprintf(stderr,"Progress =     ");
      for (i=0; i<N;i++){
         if(verbose) printCount(&i,&N_size_t,1000);

         fpixel[0] = roundToNi(xx[i]) - 1;
         fpixel[1] = roundToNi(yy[i]) - 1;

         if(xmin[0] < xx[i] && xx[i] < xmax[0] && xmin[1] < yy[i] && yy[i] < xmax[1]){
            fits_write_col(fileOutFits, typecode, ncols+1, firstrow+i, firstelem, 1, table + (fpixel[1]*naxes[0]+fpixel[0])*size, &status);
         }else{
            fits_write_col(fileOutFits, typecode, ncols+1, firstrow+i, firstelem, 1, &null, &status);
         }

      }
      fprintf(stderr,"\b\b\b\b100%%\n");
      free(xx);
      free(yy);

   }else if(checkFileExt(para->fileRegInName,".reg")){
      FILE *fileRegIn = fopenAndCheck(para->fileRegInName,"r");
      Node *polyTree  = readPolygonFileTree(fileRegIn,xmin,xmax);
      Polygon *polys = (Polygon *)polyTree->polysAll;
      Npolys          = polyTree->Npolys;
      fclose(fileRegIn);

      /*    or if the limits are defined by the user */
      if(para->minDefinied[0]) xmin[0] = para->min[0];
      if(para->maxDefinied[0]) xmax[0] = para->max[0];
      if(para->minDefinied[1]) xmin[1] = para->min[1];
      if(para->maxDefinied[1]) xmax[1] = para->max[1];

      /* print out limits */
      fprintf(stderr,"Mask limits:\n");
      fprintf(stderr,"-xmin %g -xmax %g -ymin %g -ymax %g\n",xmin[0],xmax[0],xmin[1],xmax[1]);

      /*    reference point. It must be outside the mask */
      x0[0] = xmin[0] - 1.0; x0[1] = xmin[1] - 1.0;

      fits_insert_col(fileOutFits, ncols+1, para->flagName,  "1I", &status);
	   if (status) {
         fits_report_error(stderr, status);
         exit(EXIT_FAILURE);
      }


      double *xx = (double *)malloc(N*sizeof(double));
      double *yy = (double *)malloc(N*sizeof(double));
      long *rowlist = (long *)malloc(N*sizeof(long));


      readColFits(fileCatIn, xcol, N, xx);
      readColFits(fileCatIn, ycol, N, yy);

      size_t N_size_t = (size_t)N;


      long count = 0;
      if(verbose) fprintf(stderr,"Progress =     ");
      for (i=0; i<N;i++){

         if(verbose) printCount(&i,&N_size_t,1000);

         x[0] = xx[i];
         x[1] = yy[i];

         if(flag=0,!insidePolygonTree(polyTree,x0,x,&poly_id)) flag = 1;
         fits_write_col(fileOutFits, TSHORT, ncols+1, firstrow+i, firstelem, 1, &flag, &status);
         if (status){
            fits_report_error(stderr, status);
            exit(EXIT_FAILURE);
         }

         if(para->format == 1 && !flag){
            rowlist[count] = i+1;
            count++;
         }
         if(para->format == 2 && flag){
            rowlist[count] = i+1;
            count++;
         }
      }

      fprintf(stderr,"\b\b\b\b100%%\n");
      free(xx);
      free(yy);

      if(para->format == 1 || para->format == 2){
         fits_delete_rowlist(fileOutFits, rowlist, count, &status);
         fits_delete_col(fileOutFits, ncols+1, &status);
         if (status){
            fits_report_error(stderr, status);
            exit(EXIT_FAILURE);
         }

      }

   }

	fits_close_file(fileOutFits, &status);
	if (status) {
      fits_report_error(stderr, status);
      exit(EXIT_FAILURE);
   }

	fits_close_file(fileCatIn, &status);
	if (status) {
      fits_report_error(stderr, status);
      exit(EXIT_FAILURE);
   }

   return EXIT_SUCCESS;

}


int flagCat(const Config *para){
   /*
    *    Reads fileCatIn and add a flag at the end of the line. 1 is outside
    *    the mask and 0 is inside the mask. xcol and ycol are the column ids
    *    of resp. x coordinate and y coordinate.
    *    For mask fits format, it writes the pixel value.
    */

   int Npolys, poly_id, flag, verbose = 1;
   size_t i,N, Ncol;
   double x[2], x0[2], xmin[2], xmax[2];
   char line[NFIELD*NCHAR], item[NFIELD*NCHAR],*str_end;

   FILE *fileOut   = fopenAndCheck(para->fileOutName,"w");
   FILE *fileCatIn = fopenAndCheck(para->fileCatInName,"r");

   N = 0;
   if(fileCatIn != stdin){
      while(fgets(line,NFIELD*NCHAR,fileCatIn) != NULL)
      if(getStrings(line,item," ",&Ncol))  N++;
      rewind(fileCatIn);
      fprintf(stderr,"Nobjects = %zd\n", N);
   }else{
      verbose = 0;
   }

   int xcol = atoi(para->xcol);
   int ycol = atoi(para->ycol);

   if(checkFileExt(para->fileRegInName,".fits")){

      if(para->coordType != CART){
         fprintf(stderr,"%s: fits file detected. coord should be set to cart for image coordinates. Exiting...\n",MYNAME);
         exit(EXIT_FAILURE);
      }


      long fpixel[2], naxes[2];
      int bitpix, status = 0;
      double (*toDouble)(void *,long ) = NULL;

      /*    read fits file and put in table */
      void *table = readFits(para,&bitpix,NULL,NULL,&status,naxes,&toDouble,NULL);

      /*    define limits */
      xmin[0] = xmin[1] = 0.5;
      xmax[0] = naxes[0]+0.5;
      xmax[1] = naxes[1]+0.5;
      if(para->minDefinied[0]) xmin[0] = para->min[0];
      if(para->maxDefinied[0]) xmax[0] = para->max[0];
      if(para->minDefinied[1]) xmin[1] = para->min[1];
      if(para->maxDefinied[1]) xmax[1] = para->max[1];
      /*    print out limits */
      fprintf(stderr,"Mask limits:\n");
      fprintf(stderr,"-xmin %g -xmax %g -ymin %g -ymax %g\n",xmin[0],xmax[0],xmin[1],xmax[1]);

      /*    ATTENTION "toDouble" converts everything into double */
      i = 0;
      if(verbose) fprintf(stderr,"Progress =     ");
      while(fgets(line,NFIELD*NCHAR,fileCatIn) != NULL){

         /*    keep commented lines */
         if (line[0] == '#') fprintf(fileOut,"%s",line);

         if(getStrings(line,item," ",&Ncol)){
            i++;
            if(verbose) printCount(&i,&N,1000);
            x[0] = getDoubleValue(item,xcol);
            x[1] = getDoubleValue(item,ycol);
            fpixel[0] = roundToNi(x[0]) - 1;
            fpixel[1] = roundToNi(x[1]) - 1;
            str_end = strstr(line,"\n");/* cariage return to the end of the line */
            strcpy(str_end,"\0");       /* "end" symbol to the line              */
            if(xmin[0] < x[0] && x[0] < xmax[0] && xmin[1] < x[1] && x[1] < xmax[1]){
               fprintf(fileOut,"%s %g\n",line,toDouble(table,fpixel[1]*naxes[0]+fpixel[0]));
            }else{
               fprintf(fileOut,"%s %d\n",line,-99);
            }
         }
      }
      if(verbose) fprintf(stderr,"\b\b\b\b100%%\n");

      free(table);

   }else if(checkFileExt(para->fileRegInName,".reg")){
      FILE *fileRegIn = fopenAndCheck(para->fileRegInName,"r");
      Node *polyTree  = readPolygonFileTree(fileRegIn,xmin,xmax);
      Polygon *polys = (Polygon *)polyTree->polysAll;
      Npolys          = polyTree->Npolys;
      fclose(fileRegIn);

      /*    or if the limits are defined by the user */
      if(para->minDefinied[0]) xmin[0] = para->min[0];
      if(para->maxDefinied[0]) xmax[0] = para->max[0];
      if(para->minDefinied[1]) xmin[1] = para->min[1];
      if(para->maxDefinied[1]) xmax[1] = para->max[1];
      /* print out limits */
      fprintf(stderr,"Mask limits:\n");
      fprintf(stderr,"-xmin %g -xmax %g -ymin %g -ymax %g\n",xmin[0],xmax[0],xmin[1],xmax[1]);

      /*    reference point. It must be outside the mask */
      x0[0] = xmin[0] - 1.0; x0[1] = xmin[1] - 1.0;

      i = 0;
      if(verbose) fprintf(stderr,"Progress =     ");
      while(fgets(line,NFIELD*NCHAR,fileCatIn) != NULL){

         /*    keep commented lines */
         if (line[0] == '#') fprintf(fileOut,"%s",line);

         if(getStrings(line,item," ",&Ncol)){
            i++;
            if(verbose) printCount(&i,&N,1000);
            x[0] = getDoubleValue(item,xcol);
            x[1] = getDoubleValue(item,ycol);

            if(flag=0, !insidePolygonTree(polyTree,x0,x,&poly_id)) flag = 1;

            str_end = strstr(line,"\n");/*   cariage return to the end of the line */
            strcpy(str_end,"\0");       /*   "end" symbol to the line */
            switch (para->format){
               case 1: /*  only objects outside the mask and inside the user's definied limits */
               if(flag) fprintf(fileOut,"%s\n",line);
               break;
               case 2: /*  only objects inside the mask or outside the user's definied limits */
               if(!flag) fprintf(fileOut,"%s\n",line);
               break;
               case 3: /*  all objects with the flag */
               fprintf(fileOut,"%s %d\n",line,flag);
            }
         }
      }
      fflush(stdout);
      if(verbose) fprintf(stderr,"\b\b\b\b100%%\n");
   }else{
      fprintf(stderr,"%s: mask file format not recognized. Please provide .reg, .fits or no mask with input limits. Exiting...\n",MYNAME);
      exit(EXIT_FAILURE);
   }

   fclose(fileOut);
   fclose(fileCatIn);

   return(EXIT_SUCCESS);
}



int randomCat(const Config *para){
   /*    Generates a random catalogue inside the mask (uniform PDF).
    *    If "all", it puts all objects and add a flag such as:
    *    outside the mask:1, inside the mask:0. Otherwise it puts only objects
    *    outside the mask.
    */

   int Npolys,poly_id,flag, verbose = 1, status = 0, size;
   long firstrow =1, firstelem = 1;

   size_t i, npart;
   double x[2], x0[2], xmin[2], xmax[2], z, area;
   gsl_rng *r = randomInitialize(para->seed);

   FILE *fileOut;             /* pointer to the ascii file */
   fitsfile *fileOutFits;     /* pointer to the FITS file, defined in fitsio.h */


   /*    redshift distribution
    *    GSL convention: bin[i] corresponds to range[i] <= x < range[i+1]
    */

   FILE *fileNofZ;
   size_t Nbins, Ncol;
   char line[NFIELD*NCHAR], item[NFIELD*NCHAR];

   gsl_histogram *nz;
   gsl_histogram_pdf *nz_PDF;

   if(para->nz){
      fileNofZ = fopenAndCheck(para->fileNofZName,"r");
      Nbins= 0;
      while(fgets(line,NFIELD*NCHAR, fileNofZ) != NULL)
      if(getStrings(line,item," ",&Ncol))  Nbins++;
      rewind(fileNofZ);
      fprintf(stderr,"Nbins = %zd\n",Nbins);

      nz     = gsl_histogram_alloc(Nbins);
      nz_PDF = gsl_histogram_pdf_alloc (Nbins);

      i = 0;
      while(fgets(line,NFIELD*NCHAR, fileNofZ) != NULL){
         if(getStrings(line,item," ",&Ncol)){
            nz->range[i] = getDoubleValue(item, 1);
            nz->bin[i]   = getDoubleValue(item, 2);
            i++;
         }
      }
      /*    gsl structure requires an upper limit for the last bin */
      nz->range[i]   = 100.0;
      gsl_histogram_pdf_init (nz_PDF, nz);
   }

   if(para->zrange){
      Nbins= 1000;
      double delta = (para->zmax - para->zmin)/(double)Nbins;

      nz     = gsl_histogram_alloc(Nbins);
      nz_PDF = gsl_histogram_pdf_alloc (Nbins);

      for(i=0;i<Nbins;i++){
         nz->range[i] = para->zmin + delta*(double)i;
         nz->bin[i]   = dvdz(nz->range[i] + delta/2.0, para->a);
      }
      /*    gsl structure requires an upper limit for the last bin*/
      nz->range[Nbins]   = para->zmin + delta*(double)Nbins;
      gsl_histogram_pdf_init (nz_PDF, nz);
   }

   /*
    * input mask is a fits image
    */
   if(checkFileExt(para->fileRegInName,".fits")){
      if(para->coordType != CART){
         fprintf(stderr,"%s: fits file detected. coord should be set to cart for image coordinates. Exiting...\n",MYNAME);
         exit(EXIT_FAILURE);
      }

      long fpixel[2], naxes[2];
      int bitpix, typecode;
      char tform_char;
      double (*toDouble)(void *,long ) = NULL;

      /*    read fits file and put in table */
      void *table = readFits(para, &bitpix, &typecode, &tform_char, &status, naxes, &toDouble, &size);

      /*    define limits */
      xmin[0] = xmin[1] = 0.5;
      xmax[0] = naxes[0]+0.5;
      xmax[1] = naxes[1]+0.5;
      if(para->minDefinied[0]) xmin[0] = para->min[0];
      if(para->maxDefinied[0]) xmax[0] = para->max[0];
      if(para->minDefinied[1]) xmin[1] = para->min[1];
      if(para->maxDefinied[1]) xmax[1] = para->max[1];

      /*    print out limits */
      fprintf(stderr,"Mask limits:\n");
      fprintf(stderr,"-xmin %g -xmax %g -ymin %g -ymax %g\n",xmin[0],xmax[0],xmin[1],xmax[1]);

      if(para->oFileType == FITS){
         fprintf(stderr, "Outpout file or stdout format: fits\n");

         fits_create_file(&fileOutFits, para->fileOutName, &status);
   		if (status) {
            fits_report_error(stderr, status);
            if (status){
               fprintf(stderr, "Add \"!\" in front of the file name to overwrite: -o \"!FILEOUT\"\n");
            }
            exit(EXIT_FAILURE);
         }
         /*    define the name, datatype, and physical units for the columns */
         int tfields   = 4;   /* table will have 3 columns */
         char *ttype[] = { "x", "y", para->flagName, "z" };
         char *tform[] = { "1D", "1D", concat("1", &tform_char), "1D"};
         char *tunit[] = { "pix", "pix", "\0", "\0" };
         if(para->nz || para->zrange){
            fits_create_tbl(fileOutFits, BINARY_TBL, 0, tfields, ttype, tform, tunit, "DATA", &status);
         }else{
            fits_create_tbl(fileOutFits, BINARY_TBL, 0, tfields-1, ttype, tform, tunit, "DATA", &status);
         }

         area = (xmax[0] - xmin[0])*(xmax[1] - xmin[1]);
         fprintf(stderr, "Area = %f (pix^2)\n", area);

         if(para->constDen){
            npart =  (size_t)round((double)para->npart*area / 1.e9);
         }else{
            npart = para->npart;
         }
         fprintf(stderr,"Creates a random catalogue with N = %zd objects. Format = %d\n",npart,para->format);


         fprintf(stderr,"Progress =     ");
         for(i=0;i<npart;i++){
            printCount(&i,&para->npart, 1000);
            x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
            x[1] = gsl_ran_flat(r,xmin[1],xmax[1]);
            fpixel[0] = roundToNi(x[0]) - 1;
            fpixel[1] = roundToNi(x[1]) - 1;
            fits_write_col(fileOutFits, TDOUBLE,  1, firstrow+i, firstelem, 1, &(x[0]), &status);
            fits_write_col(fileOutFits, TDOUBLE,  2, firstrow+i, firstelem, 1, &(x[1]), &status);
            fits_write_col(fileOutFits, typecode, 3, firstrow+i, firstelem, 1, table + (fpixel[1]*naxes[0]+fpixel[0])*size, &status);
            if(para->nz || para->zrange){
               z = gsl_histogram_pdf_sample (nz_PDF, gsl_ran_flat(r, 0.0, 1.0));
               fits_write_col(fileOutFits, TDOUBLE, 4, firstrow+i, firstelem, 1, &z, &status);
            }
            if (status) {
               fits_report_error(stderr, status);
               exit(EXIT_FAILURE);
            }
         }

      }else{
         fprintf(stderr, "Outpout file or stdout format: ascii\n");

         fileOut = fopenAndCheck(para->fileOutName,"w");
         /*    ATTENTION "toDouble" converts everything into double */
         fprintf(stderr,"Progress =     ");
         for(i=0;i<para->npart;i++){
            printCount(&i,&para->npart,1000);
            x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
            x[1] = gsl_ran_flat(r,xmin[1],xmax[1]);
            fpixel[0] = roundToNi(x[0]) - 1;
            fpixel[1] = roundToNi(x[1]) - 1;
            if(para->nz || para->zrange){
               z = gsl_histogram_pdf_sample (nz_PDF, gsl_ran_flat(r, 0.0, 1.0));
               fprintf(fileOut,"%f %f %g %f\n", x[0], x[1], toDouble(table,fpixel[1]*naxes[0]+fpixel[0]), z);
            }else{
               fprintf(fileOut,"%f %f %g\n", x[0], x[1], toDouble(table,fpixel[1]*naxes[0]+fpixel[0]));
            }
         }
      }

      fprintf(stderr,"\b\b\b\b100%%\n");
      free(table);

   }else if(checkFileExt(para->fileRegInName,".reg")){

      long count;

      FILE *fileRegIn = fopenAndCheck(para->fileRegInName,"r");
      Node *polyTree  = readPolygonFileTree(fileRegIn,xmin,xmax);
      Polygon *polys  = (Polygon *)polyTree->polysAll;
      Npolys          = polyTree->Npolys;
      fclose(fileRegIn);

      /*    or if the limits are defined by the user */
      if(para->minDefinied[0]) xmin[0] = para->min[0];
      if(para->maxDefinied[0]) xmax[0] = para->max[0];
      if(para->minDefinied[1]) xmin[1] = para->min[1];
      if(para->maxDefinied[1]) xmax[1] = para->max[1];
      /* print out limits */
      fprintf(stderr,"Mask limits:\n");
      fprintf(stderr,"-xmin %g -xmax %g -ymin %g -ymax %g\n",xmin[0],xmax[0],xmin[1],xmax[1]);

      /*    reference point. It must be outside the mask */
      x0[0] = xmin[0] - 1.0; x0[1] = xmin[1] - 1.0;

      //fprintf(stderr,"xmin = %f \nxmax = %f \nymin = %f \nymax = %f\n",xmin[0],xmax[0],xmin[1],xmax[1]);
      if(para->coordType == RADEC){
         area = (xmax[0] - xmin[0])*(sin(xmax[1]*PI/180.0) - sin(xmin[1]*PI/180.0))*180.0/PI;
      }else{
         area = (xmax[0] - xmin[0])*(xmax[1] - xmin[1]);
      }
      fprintf(stderr, "Area = %f\n", area);

      if(para->constDen){
         npart =  (size_t)round((double)para->npart*area);
      }else{
         npart = para->npart;
      }
      fprintf(stderr,"Creates a random catalogue with N = %zd objects. Format = %d\n",npart,para->format);

      if(para->oFileType == FITS){
         fprintf(stderr, "Outpout file or stdout format: fits\n");

         fits_create_file(&fileOutFits, para->fileOutName, &status);
   		if (status){
            fits_report_error(stderr, status);
            if (status){
               fprintf(stderr, "Add \"!\" in front of the file name to overwrite: -o \"!FILEOUT\"\n");
            }
            exit(EXIT_FAILURE);
         }

         /*    define the name, datatype, and physical units for the columns */
         if(para->nz || para->zrange){
            int tfields   = 4;   /* table will have 3 or 4 columns */
            char *ttype[] = { "ra", "dec", "z", para->flagName};
            char *tform[] = { "1D", "1D", "1D", "1I"};
            char *tunit[] = { "deg", "deg", "\0", "\0" };
            if (para->format == 1 || para->format == 2) {
               fits_create_tbl(fileOutFits, BINARY_TBL, 0, tfields-1, ttype, tform, tunit, "DATA", &status);
            }else{
               fits_create_tbl(fileOutFits, BINARY_TBL, 0, tfields, ttype, tform, tunit, "DATA", &status);
            }
         }else{
            int tfields   = 3;   /* table will have 2 or 3 columns */
            char *ttype[] = { "ra", "dec", para->flagName };
            char *tform[] = { "1D", "1D", "1I"};
            char *tunit[] = { "deg", "deg", "\0" };
            if (para->format == 1 || para->format == 2) {
               fits_create_tbl(fileOutFits, BINARY_TBL, 0, tfields-1, ttype, tform, tunit, "DATA", &status);
            }else{
               fits_create_tbl(fileOutFits, BINARY_TBL, 0, tfields, ttype, tform, tunit, "DATA", &status);
            }
         }

         if (status) {
            fits_report_error(stderr, status);
            exit(EXIT_FAILURE);
         }

         fprintf(stderr,"Progress =     ");
         count = 0;
         for(i=0;i<npart;i++){
            printCount(&i,&npart,1000);
            if(para->coordType == CART){
               x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
               x[1] = gsl_ran_flat(r,xmin[1],xmax[1]);
            }else{
               x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
               x[1] = gsl_ran_flat(r,sin(xmin[1]*PI/180.0),sin(xmax[1]*PI/180.0));
               x[1] = asin(x[1])*180.0/PI;
            }
            /*    1 = outside the mask, 0 = inside the mask */
            if(flag=0,!insidePolygonTree(polyTree,x0,x,&poly_id)) flag = 1;


            switch (para->format){
               case 1: /* only objects outside the mask */
                  if(flag){
                     fits_write_col(fileOutFits, TDOUBLE, 1, firstrow+count, firstelem, 1, &(x[0]), &status);
                     fits_write_col(fileOutFits, TDOUBLE, 2, firstrow+count, firstelem, 1, &(x[1]), &status);
                     if(para->nz || para->zrange){
                        z = gsl_histogram_pdf_sample (nz_PDF, gsl_ran_flat(r, 0.0, 1.0));
                        fits_write_col(fileOutFits, TDOUBLE, 3, firstrow+count, firstelem, 1, &z, &status);
                     }
                     if (status){
                        fits_report_error(stderr, status);
                        exit(EXIT_FAILURE);
                     }
                     count++;
                  }
                  break;
               case 2: /* only objects inside the mask */
                  if(!flag){
                     fits_write_col(fileOutFits, TDOUBLE, 1, firstrow+count, firstelem, 1, &(x[0]), &status);
                     fits_write_col(fileOutFits, TDOUBLE, 2, firstrow+count, firstelem, 1, &(x[1]), &status);
                     if(para->nz || para->zrange){
                        z = gsl_histogram_pdf_sample (nz_PDF, gsl_ran_flat(r, 0.0, 1.0));
                        fits_write_col(fileOutFits, TDOUBLE, 3, firstrow+count, firstelem, 1, &z, &status);
                     }
                     if (status) {
                        fits_report_error(stderr, status);
                        exit(EXIT_FAILURE);
                     }
                     count++;
                  }
                  break;
               case 3: /* all objects with the flag */
                  fits_write_col(fileOutFits, TDOUBLE, 1, firstrow+i, firstelem, 1, &(x[0]), &status);
                  fits_write_col(fileOutFits, TDOUBLE, 2, firstrow+i, firstelem, 1, &(x[1]), &status);
                  if(para->nz || para->zrange){
                     z = gsl_histogram_pdf_sample (nz_PDF, gsl_ran_flat(r, 0.0, 1.0));
                     fits_write_col(fileOutFits, TDOUBLE, 3, firstrow+count, firstelem, 1, &z, &status);
                     fits_write_col(fileOutFits, TSHORT, 4, firstrow+count, firstelem, 1, &flag, &status);
                  }else{
                     fits_write_col(fileOutFits, TSHORT, 3, firstrow+count, firstelem, 1, &flag, &status);
                  }
                  if (status){
                     fits_report_error(stderr, status);
                     exit(EXIT_FAILURE);
                  }
                  count++;
                  break;
            }

         }
         fprintf(stderr,"\b\b\b\b100%%\n");


      }else{
         fprintf(stderr, "Outpout file or stdout format: ascii\n");

         fileOut = fopenAndCheck(para->fileOutName,"w");

         fprintf(fileOut,"# %f\n", area);

         fprintf(stderr,"Progress =     ");
         for(i=0;i<npart;i++){
            printCount(&i,&npart,1000);
            if(para->coordType == CART){
               x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
               x[1] = gsl_ran_flat(r,xmin[1],xmax[1]);
            }else{
               x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
               x[1] = gsl_ran_flat(r,sin(xmin[1]*PI/180.0),sin(xmax[1]*PI/180.0));
               x[1] = asin(x[1])*180.0/PI;
            }
            /*    1 = outside the mask, 0 = inside the mask */
            if(flag=0,!insidePolygonTree(polyTree,x0,x,&poly_id)) flag = 1;
            if(para->nz || para->zrange){
               z = gsl_histogram_pdf_sample (nz_PDF, gsl_ran_flat(r, 0.0, 1.0));
               switch (para->format){
                  case 1: /* only objects outside the mask */
                     if(flag) fprintf(fileOut,"%f %f %f\n",x[0],x[1],z);
                     break;
                  case 2: /* only objects inside the mask */
                     if(!flag) fprintf(fileOut,"%f %f %f\n",x[0],x[1],z);
                     break;
                  case 3: /* all objects with the flag */
                     fprintf(fileOut,"%f %f %d %f\n",x[0],x[1],flag,z);
                     break;
               }
            }else{
               switch (para->format){
                  case 1: /* only objects outside the mask */
                     if(flag) fprintf(fileOut,"%f %f\n",x[0],x[1]);
                     break;
                  case 2: /* only objects inside the mask */
                     if(!flag) fprintf(fileOut,"%f %f\n",x[0],x[1]);
                     break;
                  case 3: /* all objects with the flag */
                     fprintf(fileOut,"%f %f %d\n",x[0],x[1],flag);
                     break;
               }
            }
         }
         fprintf(stderr,"\b\b\b\b100%%\n");

      }


   }else if(!strcmp(para->fileRegInName,"\0")){

      fprintf(stderr,"Generating catalogue with no mask...\n");
      xmin[0] = para->min[0];
      xmax[0] = para->max[0];
      xmin[1] = para->min[1];
      xmax[1] = para->max[1];

      /*    print out limits */
      fprintf(stderr,"limits:\n");
      fprintf(stderr,"-xmin %g -xmax %g -ymin %g -ymax %g\n",xmin[0],xmax[0],xmin[1],xmax[1]);
      if(para->coordType == RADEC){
         area = (xmax[0] - xmin[0])*(sin(xmax[1]*PI/180.0) - sin(xmin[1]*PI/180.0))*180.0/PI;
      }else{
         area = (xmax[0] - xmin[0])*(xmax[1] - xmin[1]);
      }
      fprintf(stderr, "Area = %f\n", area);

      if(para->constDen){
         npart = (size_t)round((double)para->npart*area);
      }else{
         npart = para->npart;
      }

      if(para->oFileType == FITS){
         fprintf(stderr, "Outpout file or stdout format: fits\n");

         fits_create_file(&fileOutFits, para->fileOutName, &status);
   		if (status){
            fits_report_error(stderr, status);
            exit(EXIT_FAILURE);
         }

         /*    define the name, datatype, and physical units for the columns */
         int tfields   = 3;   /* table will have 3 columns */
         char *ttype[] = { "ra", "dec", "z" };
         char *tform[] = { "1D", "1D", "1D"};
         char *tunit[] = { "deg", "deg", "\0" };
         if(para->coordType == CART){
            ttype[0] = "x";
            ttype[1] = "y";
            tunit[0] = "pix";
            tunit[1] = "pix";
         }
         if(para->nz || para->zrange){
            fits_create_tbl(fileOutFits, BINARY_TBL, 0, tfields, ttype, tform, tunit, "DATA", &status);
         }else{
            fits_create_tbl(fileOutFits, BINARY_TBL, 0, tfields-1, ttype, tform, tunit, "DATA", &status);
         }

         fprintf(stderr,"Creates a random catalogue with N = %zd objects.\n",npart);
         fprintf(stderr,"Progress =     ");
         for(i=0;i<npart;i++){
            printCount(&i,&npart,1000);
            if(para->coordType == CART){
               x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
               x[1] = gsl_ran_flat(r,xmin[1],xmax[1]);
            }else{
               x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
               x[1] = gsl_ran_flat(r,sin(xmin[1]*PI/180.0),sin(xmax[1]*PI/180.0));
               x[1] = asin(x[1])*180.0/PI;
            }
            fits_write_col(fileOutFits, TDOUBLE,  1, firstrow+i, firstelem, 1, &(x[0]), &status);
            fits_write_col(fileOutFits, TDOUBLE,  2, firstrow+i, firstelem, 1, &(x[1]), &status);
            if(para->nz || para->zrange){
               z = gsl_histogram_pdf_sample (nz_PDF, gsl_ran_flat(r, 0.0, 1.0));
               fits_write_col(fileOutFits, TDOUBLE, 3, firstrow+i, firstelem, 1, &z, &status);
            }
            if (status) {
               fits_report_error(stderr, status);
               exit(EXIT_FAILURE);
            }
         }

         fprintf(stderr,"\b\b\b\b100%%\n");

      }else{

         fprintf(stderr, "Outpout file or stdout format: ascii\n");

         fileOut = fopenAndCheck(para->fileOutName,"w");

         fprintf(fileOut, "# %f\n", area);

         fprintf(stderr,"Creates a random catalogue with N = %zd objects.\n",npart);
         fprintf(stderr,"Progress =     ");
         for(i=0;i<npart;i++){
            printCount(&i,&npart,1000);
            if(para->coordType == CART){
               x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
               x[1] = gsl_ran_flat(r,xmin[1],xmax[1]);
            }else{
               x[0] = gsl_ran_flat(r,xmin[0],xmax[0]);
               x[1] = gsl_ran_flat(r,sin(xmin[1]*PI/180.0),sin(xmax[1]*PI/180.0));
               x[1] = asin(x[1])*180.0/PI;
            }
            if(para->nz || para->zrange){
               z = gsl_histogram_pdf_sample (nz_PDF, gsl_ran_flat(r, 0.0, 1.0));
               fprintf(fileOut,"%f %f %f\n", x[0], x[1], z);
            }else{
               fprintf(fileOut,"%f %f\n", x[0], x[1]);
            }
         }
         fflush(stdout);
         fprintf(stderr,"\b\b\b\b100%%\n");

      }

   }else{

      fprintf(stderr,"%s: mask file format not recognized. Please provide .reg, .fits or no mask with input limits. Exiting...\n",MYNAME);
      exit(EXIT_FAILURE);

   }

   if(para->nz || para->zrange){
      gsl_histogram_pdf_free (nz_PDF);
      gsl_histogram_free (nz);
   }

   if(para->oFileType == FITS){
   	fits_close_file(fileOutFits, &status);
   	if (status){
         fits_report_error(stderr, status);
         exit(EXIT_FAILURE);
      }
   }else{
      fclose(fileOut);
   }

   return(EXIT_SUCCESS);
}
