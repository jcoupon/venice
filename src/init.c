/*
 *    init.c
 *    venice
 *    Jean Coupon (2012-2017)
 */

#include "init.h"

/*
 *		Initialization
 */

int readParameters(int argc, char **argv, Config *para){
	int i,task,nomask;
	char list[NFIELD*NCHAR];
	size_t Ncol;

	/* 	Default parameters */
	nomask          = 1;
	task            = 1;
	para->nx        = 512;
	para->ny        = 512;
	para->xcol = malloc((72+1) * sizeof(char));
	para->ycol = malloc((72+1) * sizeof(char));
	strcpy(para->xcol,"1");
	strcpy(para->ycol,"2");
	para->npart     = 1000000;
	para->format    = 1;
	para->coordType = CART;
	para->seed      = 20091982;
	para->constDen  = 0;
	para->nz        = 0;
	para->zrange    = 0;
	para->catFileType = FITS;
	para->oFileType = FITS;


	/* 	default cosmology <=> WMAP5 */
	/*		TODO: add it as an option */
	para->a[0] = H0;
	para->a[1] = Omega_M;
	para->a[2] = Omega_L;
	para->a[3] = c;

	for(i=0;i<2;i++){
		para->minDefinied[i] = 0;
		para->maxDefinied[i] = 0;
		para->min[i] = 0.0;
		para->max[i] = 0.0;
	}


	strcpy(MYNAME,"venice");
	strcpy(para->fileOutName,"");
	strcpy(para->fileCatInName,"\0");
	strcpy(para->fileRegInName,"\0");
	strcpy(para->fileNofZName,"\0");
	strcpy(para->flagName,"flag");



	for(i=0;i<argc;i++){
		/*		help */
		if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") || argc == 1){
			fprintf(stderr,"\n\n                   V E N I C E\n\n");
			fprintf(stderr,"           mask utility program version 4.0.3 \n\n");
			fprintf(stderr,"Usage: %s -m mask.[reg,fits]               [OPTIONS] -> binary mask for visualization\n",argv[0]);
			fprintf(stderr,"    or %s -m mask.[reg,fits] -cat file.cat [OPTIONS] -> objects in/out of mask\n",argv[0]);
			fprintf(stderr,"    or %s -m mask.[reg,fits] -cat -        [OPTIONS] -> objects in/out of mask (from stdin)\n",argv[0]);
			fprintf(stderr,"    or %s -m mask.[reg,fits] -r            [OPTIONS] -> random catalogue\n\n",argv[0]);
			fprintf(stderr,"Options:\n");
			fprintf(stderr,"    -r                       create random catalogue\n");
			fprintf(stderr,"    -cat FILE                input catalogue file name, default:stdin\n");
			fprintf(stderr,"    -o FILE                  output file name, default:stdout\n");
			fprintf(stderr,"    -catfmt [ascii,fits]     input catalogue format, default:fits if stdin\n");
			fprintf(stderr,"    -ofmt [ascii,fits]       output file format, default:fits if stdout\n");
			fprintf(stderr,"    -f [outside,inside,all]  output format, default:outside\n");
			fprintf(stderr,"    -[x,y]col N              column id for x and y (starts at 1)\n");
			fprintf(stderr,"                             One may use column names for input fits catalogue\n");
			fprintf(stderr,"    -coord [cart,spher]      coordinate type, default:cart\n");
			fprintf(stderr,"    -[x,y]min value          lower limit for x and y\n");
			fprintf(stderr,"    -[x,y]max value          upper limit for x and y\n");
			fprintf(stderr,"    -nz file_nz.in           redshift distribution for random objects\n");
			fprintf(stderr,"    -z zmin,zmax             redshift range for random objects (if volume limited)\n");
			fprintf(stderr,"    -seed  N                 random seed\n");
			fprintf(stderr,"    -npart N                 number of random objects\n");
			fprintf(stderr,"    -cd                      multiply npart by the mask area (for constant density)\n");
			fprintf(stderr,"    -flagName NAME           name of the flag colum for fits files. default: flag\n");
			fprintf(stderr,"    -h, --help               this message\n\n");
			fprintf(stderr,"For .reg files, the region must be \"polygon\", \"circle\", \"ellipse\" or \"box\".\n");
			fprintf(stderr,"Notice: 0 means inside the mask, 1 outside; for .fits files,\n");
			fprintf(stderr,"the pixel value is added at the end of the line\n");
			exit(EXIT_FAILURE);
		}
		/*		polygon file in */
		if(!strcmp(argv[i],"-m")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			strcpy(para->fileRegInName,argv[i+1]);
			nomask = 0;
		}
		/*		input catalogue (if -cat set) */
		if(!strcmp(argv[i],"-cat")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			strcpy(para->fileCatInName,argv[i+1]);
			task = 2;
		}
		/*		random catalogue */
		if(!strcmp(argv[i],"-r")){
			task = 3;
		}
		/*		constant density */
		if(!strcmp(argv[i],"-cd")){
			para->constDen = 1;
		}
		/*		output file */
		if(!strcmp(argv[i],"-o")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			strcpy(para->fileOutName,argv[i+1]);
		}
		/*		dimensions of the pixel mask */
		if(!strcmp(argv[i],"-nx")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			para->nx = atoi(argv[i+1]);
		}
		if(!strcmp(argv[i],"-ny")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			para->ny = atoi(argv[i+1]);
		}
		/*		column ids for the coordinates in FILE_CAT_IN */
		if(!strcmp(argv[i],"-xcol")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			strcpy(para->xcol,argv[i+1]);
			// para->xcol = atoi(argv[i+1]);
		}
		if(!strcmp(argv[i],"-ycol")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			strcpy(para->ycol,argv[i+1]);
			// para->ycol = atoi(argv[i+1]);
		}
		/*		NPART for the random catalogue */
		if(!strcmp(argv[i],"-npart")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			para->npart = atoi(argv[i+1]);
		}
		/*		option for the random catalogue */
		if(!strcmp(argv[i],"-f")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			if(!strcmp(argv[i+1],"outside")) para->format = 1;
			if(!strcmp(argv[i+1],"inside"))  para->format = 2;
			if(!strcmp(argv[i+1],"all"))     para->format = 3;
		}
		/*		limits for the random catalogue */
		if(!strcmp(argv[i],"-xmin")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			para->minDefinied[0] = 1;
			para->min[0] = atof(argv[i+1]);
		}
		if(!strcmp(argv[i],"-xmax")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			para->maxDefinied[0] = 1;
			para->max[0] = atof(argv[i+1]);
		}
		if(!strcmp(argv[i],"-ymin")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			para->minDefinied[1] = 1;
			para->min[1] = atof(argv[i+1]);
		}
		if(!strcmp(argv[i],"-ymax")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			para->maxDefinied[1] = 1;
			para->max[1] = atof(argv[i+1]);
		}
		/*		input catalogue (if -cat set) */
		if(!strcmp(argv[i],"-nz")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			strcpy(para->fileNofZName,argv[i+1]);
			if(para->zrange == 1){
				fprintf(stderr,"Please choose between nz (n(z) file) and z (redshift range, volume limited case)\n");
				exit(-1);
			}
			para->nz = 1;
		}
		if(!strcmp(argv[i],"-z")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			getStrings(argv[i+1],list,",",&Ncol);
			para->zmin = getDoubleValue(list,1);
			para->zmax = getDoubleValue(list,2);
			if(para->zmax < para->zmin){
				fprintf(stderr,"zmin must be lower than zmax...Exiting...\n");
				exit(-1);
			}
			if(para->nz == 1){
				fprintf(stderr,"Please choose between nz (n(z) file) and z (redshift range, volume limited case)\n");
				exit(-1);
			}
			para->zrange = 1;
		}
		/* 	coordinates type */
		if(!strcmp(argv[i],"-coord")) {
			if(!strcmp(argv[i+1],"spher")){
				para->coordType = RADEC;
			}else if(!strcmp(argv[i+1],"cart")) {
				para->coordType = CART;
			}
		}
		/*		random seed */
		if(!strcmp(argv[i],"-seed")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			if(atoi(argv[i+1]) > 0) para->seed = atoi(argv[i+1]);
		}

		/* 	input file type */
		if(!strcmp(argv[i],"-catfmt")) {
			if(!strcmp(argv[i+1],"ascii")){
				para->catFileType = ASCII;
			}else if(!strcmp(argv[i+1],"fits")) {
				para->catFileType = FITS;
			}
		}

		/* 	output file type */
		if(!strcmp(argv[i],"-ofmt")) {
			if(!strcmp(argv[i+1],"ascii")){
				para->oFileType = ASCII;
			}else if(!strcmp(argv[i+1],"fits")) {
				para->oFileType = FITS;
			}
		}

		/*		output file */
		if(!strcmp(argv[i],"-flagName")){
			if(argv[i+1] == NULL){
				fprintf(stderr,"Missing argument after %s\nExiting...\n",argv[i]);
				exit(-1);
			}
			strcpy(para->flagName,argv[i+1]);
		}



	}




	/* 	if no mask file is provided */
	if (nomask){
		/*		check if all the limits are definied in this case; */
		if(task == 3 && (para->minDefinied[0] + para->minDefinied[1] + para->maxDefinied[0] + para->maxDefinied[1] < 4)) {
			fprintf(stderr,"If you want to generate a random catalogue with no mask,\n");
			fprintf(stderr,"please provide all the coordinate limits:\n");
			fprintf(stderr,"%s -r -xmin value -xmax value -ymin value -ymax value [OPTIONS]\n",argv[0]);
			exit(-1);
		}
		/*		Checks if the limits are realistics (min < max) */
		if(task == 3 && ( para->min[0] > para->max[0] || para->min[1] > para->max[1])){
			fprintf(stderr,"Please put realistic limits (xmin < xmax and ymin < ymax).\n");
			exit(EXIT_FAILURE);
		}
		if(task != 3){
			fprintf(stderr,"please provide a mask file.\n");
			fprintf(stderr,"Usage: %s -m mask.reg               [OPTIONS]\n",argv[0]);
			fprintf(stderr,"    or %s -m mask.reg -cat file.cat [OPTIONS]\n",argv[0]);
			fprintf(stderr,"    or %s -m mask.reg -r            [OPTIONS]\n",argv[0]);
			exit(EXIT_FAILURE);
		}
	}


	/*		if input file ends with .ascii or .cat set ascii format */
   if(checkFileExt(para->fileCatInName,".ascii") || checkFileExt(para->fileCatInName,".cat")){
		para->catFileType = ASCII;
   }
	/*		if output file ends with .ascii or .cat set ascii format  */
   if(checkFileExt(para->fileOutName,".ascii") || checkFileExt(para->fileOutName,".cat")){
		para->oFileType = ASCII;
   }

	/* 	stdout if no output file given */
	if( !strcmp(para->fileOutName, "\0")){
		strcpy(para->fileOutName,"-");
	}
	if(task == 2){
		/* 	input and output format must be the same */
		if( (para->catFileType == ASCII && para->oFileType == FITS) || (para->catFileType == FITS && para->oFileType == ASCII)){
			fprintf(stderr,"The input and output catalogues must have the same format. Exiting...\n");
			exit(EXIT_FAILURE);
		}
	}


	return task;
}
