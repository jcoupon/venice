/*
 *    utils.c
 *    venice
 *    Jean Coupon (2012-2017)
 */

#include "utils.h"

/*
 *		Utils - geometric
 */

int insidePolygonTree(Node *polyTree, double x0[2], double x[2], int *poly_id){
	/*		Returns 1 if the point (x,y) is inside one of the polygons in
	 *		polys. Returns 0 if the object is oustide of any polygon or outside the
	 *		mask limits. See insidePolygon() for the algorithm explanations.
	 */

	int i,j,k,Ncross,result;
	double s,t,D;

	if(polyTree->Npolys == 0){
		*poly_id = -1;
		return 0;
	}

	if(polyTree->type != LEAF){
		if(x[polyTree->SplitDim] < polyTree->SplitValue){
			result = insidePolygonTree(polyTree->Left,  x0, x, poly_id);
		}else{
			result = insidePolygonTree(polyTree->Right, x0, x, poly_id);
		}
	}else{
		Polygon *polys = (Polygon *)polyTree->polysAll;
		for(k=0;k<polyTree->Npolys;k++){
			i = polyTree->poly_id[k];
			if(polys[i].xmin[0] < x[0] && x[0] < polys[i].xmax[0] && polys[i].xmin[1] < x[1] && x[1] < polys[i].xmax[1]){
				/* the object is inside the square around the polygon */
				Ncross=0;
				for(j=0;j<polys[i].N;j++){
					if(j<polys[i].N-1){
						D = (polys[i].x[j+1]-polys[i].x[j])*(x[1]-x0[1])-(polys[i].y[j+1]-polys[i].y[j])*(x[0]-x0[0]);
						s = ((x[0]-x0[0])*(polys[i].y[j]-x[1])-(x[1]-x0[1])*(polys[i].x[j]-x[0]))/D;
						t = ((polys[i].x[j]-x[0])*(polys[i].y[j+1]-polys[i].y[j])-(polys[i].y[j]-x[1])*(polys[i].x[j+1]-polys[i].x[j]))/D;
					}else{
						D = (polys[i].x[0]-polys[i].x[j])*(x[1]-x0[1])-(polys[i].y[0]-polys[i].y[j])*(x[0]-x0[0]);
						s = ((x[0]-x0[0])*(polys[i].y[j]-x[1])-(x[1]-x0[1])*(polys[i].x[j]-x[0]))/D;
						t = ((polys[i].x[j]-x[0])*(polys[i].y[0]-polys[i].y[j])-(polys[i].y[j]-x[1])*(polys[i].x[0]-polys[i].x[j]))/D;
					}
					if(0.0 < s && s < 1.0 + EPS && 0.0 < t && t < 1.0 + EPS) Ncross++;
				}
				if(GSL_IS_ODD(Ncross)){
					*poly_id = i;
					return 1;
				}
			}
		}
		*poly_id = -1;
		return 0;
	}

	return result;
}

int insidePolygon(Polygon *polys, int Npolys, double x0, double y0, double x, double y, int *poly_id){

	/**************************** NO LONGER USED ****************************** /

	Returns 1 if the point (x,y) is inside one of the polygons in
	polys. The id of the first polygon in which the point is
	found is also returned in poly_id. If the point is found to be outside
	all polygons, the function returns 0 and poly_id is set to -1.
	The function first checks if the point is inside the square drawn
	by the extrema of each polygon ("computationaly" more effecient).
	If yes, it counts how many times the segment {(x0,y0),(x,y)} crosses
	the sides of the polygon (Ncross). (x,y) inside the polygon
	implies Ncross is an odd number if the point (x0,y0) is outside the polygon
	(then the point (x0,y0) must be chosen outside any polygon).*/

	int i,j,Ncross;
	double s,t,D;


	for(i=0;i<Npolys;i++){

		if(polys[i].xmin[0] < x && x < polys[i].xmax[0] && polys[i].xmin[1] < y && y < polys[i].xmax[1]){
			/*		the object is inside the square around the polygon */
			Ncross=0;
			for(j=0;j<polys[i].N;j++){
				if(j<polys[i].N-1){
					D = (polys[i].x[j+1]-polys[i].x[j])*(y-y0)-(polys[i].y[j+1]-polys[i].y[j])*(x-x0);
					s = ((x-x0)*(polys[i].y[j]-y)-(y-y0)*(polys[i].x[j]-x))/D;
					t = ((polys[i].x[j]-x)*(polys[i].y[j+1]-polys[i].y[j])-(polys[i].y[j]-y)*(polys[i].x[j+1]-polys[i].x[j]))/D;
				}else{
					D = (polys[i].x[0]-polys[i].x[j])*(y-y0)-(polys[i].y[0]-polys[i].y[j])*(x-x0);
					s = ((x-x0)*(polys[i].y[j]-y)-(y-y0)*(polys[i].x[j]-x))/D;
					t = ((polys[i].x[j]-x)*(polys[i].y[0]-polys[i].y[j])-(polys[i].y[j]-y)*(polys[i].x[0]-polys[i].x[j]))/D;
				}
				if(0.0 < s && s < 1.0 + EPS && 0.0 < t && t < 1.0 + EPS) Ncross++;
			}
			if(GSL_IS_ODD(Ncross)){
				*poly_id = i;
				return 1;
			}
		}
	}
	*poly_id = -1;
	return 0;
}

Polygon *readPolygonFile(FILE *fileIn, int *Npolys, Node *polyTree){
	Polygon *result = malloc(sizeof(Polygon));
	return result;
}

Node *readPolygonFileTree(FILE *fileIn, double xmin[2], double xmax[2]){
	/* 	Reads the file file_in and returns the polygons tree.
	 * 	See http://hea-www.harvard.edu/RD/ds9/ref/region.html
	 */

	 fprintf(stderr,"Reading region mask file...");


	char line[NFIELD*NCHAR], item[NFIELD*NCHAR],*str_begin,*str_end;
	int i,j, spherical;
	size_t N, NpolysAll;
	double x0, y0, x, y, r, rx, ry, alpha, angle;

	NpolysAll = 0;
	/*		Read the entire file and count the total number of polygons, NpolysAll. */
	while(fgets(line,NFIELD*NCHAR,fileIn) != NULL)
	if(strstr(line,"polygon") != NULL || strstr(line,"circle") != NULL || strstr(line,"ellipse") != NULL || strstr(line,"box") != NULL) NpolysAll += 1;
	rewind(fileIn);
	Polygon *polysAll = (Polygon *)malloc(NpolysAll*sizeof(Polygon));
	// Polygon *polysAll;

	/*
	for(i=0; i<NpolysAll; i++){
		polysAll[i].x    = (double *)malloc(NVERTICES*sizeof(double));
		polysAll[i].y    = (double *)malloc(NVERTICES*sizeof(double));
		polysAll[i].xmin = (double *)malloc(2*sizeof(double));
		polysAll[i].xmax = (double *)malloc(2*sizeof(double));
	}
	*/
	i=0;
	/*		Read the file and fill the array with polygons. */
	while(fgets(line,NFIELD*NCHAR,fileIn) != NULL){

		if(strstr(line,"polygon") != NULL){

			str_begin = strstr(line,"(")+sizeof(char);
			str_end = strstr(line,")");
			strcpy(str_end,"\n\0");
			//strcpy(line,str_begin);
			//getStrings(line,item,",",&N);
			//DEBUGGING
			getStrings(str_begin, item, ",", &N);
			/*
			 * 	get all coordinates separated by comas.
			 */
			polysAll[i].N = N/2;
			if(N/2 > NVERTICES){
				fprintf(stderr,"%s: %zd = too many points for polygon %d (%d maxi). Exiting...\n",MYNAME,N/2,i,NVERTICES);
				exit(EXIT_FAILURE);
			}

			polysAll[i].x    = (double *)malloc(polysAll[i].N*sizeof(double));
			polysAll[i].y    = (double *)malloc(polysAll[i].N*sizeof(double));
			polysAll[i].xmin = (double *)malloc(2*sizeof(double));
			polysAll[i].xmax = (double *)malloc(2*sizeof(double));

			polysAll[i].id      = i;
			polysAll[i].xmin[0] = atof(item);
			polysAll[i].xmax[0] = atof(item);
			polysAll[i].xmin[1] = atof(item+NCHAR);
			polysAll[i].xmax[1] = atof(item+NCHAR);
			for(j=0;j<N/2;j++){
				polysAll[i].x[j]    = atof(item+NCHAR*2*j);
				polysAll[i].y[j]    = atof(item+NCHAR*(2*j+1));
				polysAll[i].xmin[0] = MIN(polysAll[i].xmin[0], polysAll[i].x[j]);
				polysAll[i].xmax[0] = MAX(polysAll[i].xmax[0], polysAll[i].x[j]);
				polysAll[i].xmin[1] = MIN(polysAll[i].xmin[1], polysAll[i].y[j]);
				polysAll[i].xmax[1] = MAX(polysAll[i].xmax[1], polysAll[i].y[j]);
			}

			//      fprintf(stderr,"%f %f %f %f\n", polysAll[i].xmin[0], polysAll[i].xmax[0],polysAll[i].xmin[1], polysAll[i].xmax[1]);

			i++;
		}else if(strstr(line,"circle") != NULL){

			str_begin = strstr(line,"(")+sizeof(char);
			str_end = strstr(line,")");
			strcpy(str_end,"\n\0");
			//strcpy(line,str_begin);
			//getStrings(line,item,",",&N);
			//DEBUGGING
			getStrings(str_begin, item, ",", &N);

			polysAll[i].N = 40;
			x0 = atof(item+NCHAR*0);
			y0 = atof(item+NCHAR*1);

			if(strstr(item+NCHAR*2,"\"") != NULL){
				spherical = 1;
				r  =  atof(item+NCHAR*2)/3600.0;
			}else if(strstr(item+NCHAR*2,"\'") != NULL){
				spherical = 1;
				r  =  atof(item+NCHAR*2)/60.0;
			}else if(strstr(item+NCHAR*2,"d") != NULL){
				spherical = 1;
				r  =  atof(item+NCHAR*2);
			}else{
				spherical = 0;
				r  =  atof(item+NCHAR*2);
			}

			polysAll[i].x    = (double *)malloc(polysAll[i].N*sizeof(double));
			polysAll[i].y    = (double *)malloc(polysAll[i].N*sizeof(double));
			polysAll[i].xmin = (double *)malloc(2*sizeof(double));
			polysAll[i].xmax = (double *)malloc(2*sizeof(double));

			polysAll[i].id      = i;
			polysAll[i].xmin[0] = x0;
			polysAll[i].xmax[0] = x0;
			polysAll[i].xmin[1] = y0;
			polysAll[i].xmax[1] = y0;
			for(j=0;j<polysAll[i].N;j++){
				alpha            = TWOPI*(double)j/((double)polysAll[i].N);
				polysAll[i].x[j] = r*cos(alpha) + x0;
				polysAll[i].y[j] = r*sin(alpha) + y0;

				if(spherical){
					polysAll[i].x[j] = 2.0*asin(sin((polysAll[i].x[j]-x0)*PI/180/2.0)/cos(polysAll[i].y[j]*PI/180))*180.0/PI+x0;
				}

				/* DEBUGGING
					fprintf(stderr,"%f %f\n",r, distAngSpher(x0, y0, polysAll[i].x[j], polysAll[i].y[j]));
				*/

				polysAll[i].xmin[0] = MIN(polysAll[i].xmin[0], polysAll[i].x[j]);
				polysAll[i].xmax[0] = MAX(polysAll[i].xmax[0], polysAll[i].x[j]);
				polysAll[i].xmin[1] = MIN(polysAll[i].xmin[1], polysAll[i].y[j]);
				polysAll[i].xmax[1] = MAX(polysAll[i].xmax[1], polysAll[i].y[j]);
			}
			i++;

		}else if(strstr(line,"ellipse") != NULL){

			str_begin = strstr(line,"(")+sizeof(char);
			str_end = strstr(line,")");
			strcpy(str_end,"\n\0");
			//strcpy(line,str_begin);
			//getStrings(line,item,",",&N);
			//DEBUGGING
			getStrings(str_begin, item, ",", &N);

			polysAll[i].N = 40;
			x0 = atof(item+NCHAR*0);
			y0 = atof(item+NCHAR*1);
			if(strstr(item+NCHAR*2,"\"") != NULL){
				spherical = 1;
				rx  =  atof(item+NCHAR*2)/3600.0;
				ry  =  atof(item+NCHAR*3)/3600.0;
			}else if(strstr(item+NCHAR*2,"\'") != NULL){
				spherical = 1;
				rx  =  atof(item+NCHAR*2)/60.0;
				ry  =  atof(item+NCHAR*3)/60.0;
			}else if(strstr(item+NCHAR*2,"d") != NULL){
				spherical = 1;
				rx  =  atof(item+NCHAR*2);
				ry  =  atof(item+NCHAR*3);
			}else{
				spherical = 0;
				rx  =  atof(item+NCHAR*2);
				ry  =  atof(item+NCHAR*3);
			}
			angle = atof(item+NCHAR*4)*PI/180.0;

			polysAll[i].x    = (double *)malloc(polysAll[i].N*sizeof(double));
			polysAll[i].y    = (double *)malloc(polysAll[i].N*sizeof(double));
			polysAll[i].xmin = (double *)malloc(2*sizeof(double));
			polysAll[i].xmax = (double *)malloc(2*sizeof(double));

			polysAll[i].id      = i;
			polysAll[i].xmin[0] = x0;
			polysAll[i].xmax[0] = x0;
			polysAll[i].xmin[1] = y0;
			polysAll[i].xmax[1] = y0;
			for(j=0;j<polysAll[i].N;j++){
				alpha            = TWOPI*(double)j/((double)polysAll[i].N);
				x = rx*cos(alpha) + x0;
				y = ry*sin(alpha) + y0;

				if(spherical){
					x = 2.0*asin(sin((x-x0)*PI/180/2.0)/cos(y*PI/180))*180.0/PI+x0;
				}

				rotate(x0, y0, x, y, &(polysAll[i].x[j]), &(polysAll[i].y[j]), angle, spherical);

				polysAll[i].xmin[0] = MIN(polysAll[i].xmin[0], polysAll[i].x[j]);
				polysAll[i].xmax[0] = MAX(polysAll[i].xmax[0], polysAll[i].x[j]);
				polysAll[i].xmin[1] = MIN(polysAll[i].xmin[1], polysAll[i].y[j]);
				polysAll[i].xmax[1] = MAX(polysAll[i].xmax[1], polysAll[i].y[j]);
			}
			i++;

		}

		/*
		else if(strstr(line,"newbox") != NULL){

			str_begin = strstr(line,"(")+sizeof(char);
			str_end = strstr(line,")");
			strcpy(str_end,"\n\0");

			getStrings(str_begin, item, ",", &N);

			polysAll[i].N = 40;
			x0 = atof(item+NCHAR*0);
			y0 = atof(item+NCHAR*1);
			if(strstr(item+NCHAR*2,"\"") != NULL){
				spherical = 1;
				rx  =  atof(item+NCHAR*2)/3600.0;
				ry  =  atof(item+NCHAR*3)/3600.0;
				rx  = 2.0*asin(sin(rx*PI/180.0/2.0)/cos(y0*PI/180.0))*180.0/PI;
			}else if(strstr(item+NCHAR*2,"\'") != NULL){
				spherical = 1;
				rx  =  atof(item+NCHAR*2)/60.0;
				ry  =  atof(item+NCHAR*3)/60.0;
				rx  = 2.0*asin(sin(rx*PI/180.0/2.0)/cos(y0*PI/180.0))*180.0/PI;
			}else if(strstr(item+NCHAR*2,"d") != NULL){
				spherical = 1;
				rx  =  atof(item+NCHAR*2);
				ry  =  atof(item+NCHAR*3);
				rx  = 2.0*asin(sin(rx*PI/180.0/2.0)/cos(y0*PI/180.0))*180.0/PI;
			}else{
				spherical = 0;
				rx  =  atof(item+NCHAR*2);
				ry  =  atof(item+NCHAR*3);
			}
			angle = atof(item+NCHAR*4)*PI/180.0;

			polysAll[i].x    = (double *)malloc(polysAll[i].N*sizeof(double));
			polysAll[i].y    = (double *)malloc(polysAll[i].N*sizeof(double));
			polysAll[i].xmin = (double *)malloc(2*sizeof(double));
			polysAll[i].xmax = (double *)malloc(2*sizeof(double));

			for(j=0;j<polysAll[i].N;j++){
				alpha            = TWOPI*(double)j/((double)polysAll[i].N);
				x = rx*cos(alpha) + x0;
				y = ry*sin(alpha) + y0;

				if(spherical){
					x = 2.0*asin(sin((x-x0)*PI/180/2.0)/cos(y*PI/180))*180.0/PI+x0;
				}
				rotate(x0, y0, x, y, &(polysAll[i].x[j]), &(polysAll[i].y[j]), angle, spherical);
			}

			polysAll[i].id      = i;
			polysAll[i].xmin[0] = x0;
			polysAll[i].xmax[0] = x0;
			polysAll[i].xmin[1] = y0;
			polysAll[i].xmax[1] = y0;
			for(j=0;j<polysAll[i].N;j++){
				polysAll[i].xmin[0] = MIN(polysAll[i].xmin[0], polysAll[i].x[j]);
				polysAll[i].xmax[0] = MAX(polysAll[i].xmax[0], polysAll[i].x[j]);
				polysAll[i].xmin[1] = MIN(polysAll[i].xmin[1], polysAll[i].y[j]);
				polysAll[i].xmax[1] = MAX(polysAll[i].xmax[1], polysAll[i].y[j]);
			}

			//fprintf(stderr,"\nxmin = %f \nxmax = %f \nymin = %f \nymax = %f\n",polysAll[i].xmin[0],polysAll[i].xmax[0],polysAll[i].xmin[1],polysAll[i].xmax[1]);

			i++;
		}
*/

		else if(strstr(line,"box") != NULL){

			str_begin = strstr(line,"(")+sizeof(char);
			str_end = strstr(line,")");
			strcpy(str_end,"\n\0");

			//strcpy(line,str_begin);
			//getStrings(line,item,",",&N);
			//DEBUGGING
			getStrings(str_begin, item, ",", &N);

			polysAll[i].N = 4;
			x0 = atof(item+NCHAR*0);
			y0 = atof(item+NCHAR*1);
			if(strstr(item+NCHAR*2,"\"") != NULL){
				spherical = 1;
				rx  =  atof(item+NCHAR*2)/3600;
				ry  =  atof(item+NCHAR*3)/3600;
				rx  = 2.0*asin(sin(rx*PI/180.0/2.0)/cos(y0*PI/180.0))*180.0/PI;
			}else if(strstr(item+NCHAR*2,"\'") != NULL){
				spherical = 1;
				rx  =  atof(item+NCHAR*2)/60.0;
				ry  =  atof(item+NCHAR*3)/60.0;
				rx  = 2.0*asin(sin(rx*PI/180.0/2.0)/cos(y0*PI/180.0))*180.0/PI;
			}else if(strstr(item+NCHAR*2,"d") != NULL){
				spherical = 1;
				rx  =  atof(item+NCHAR*2);
				ry  =  atof(item+NCHAR*3);
				rx  = 2.0*asin(sin(rx*PI/180.0/2.0)/cos(y0*PI/180.0))*180.0/PI;
			}else{
				spherical = 0;
				rx  =  atof(item+NCHAR*2);
				ry  =  atof(item+NCHAR*3);
			}
			angle = atof(item+NCHAR*4)*PI/180.0;

			polysAll[i].x    = (double *)malloc(polysAll[i].N*sizeof(double));
			polysAll[i].y    = (double *)malloc(polysAll[i].N*sizeof(double));
			polysAll[i].xmin = (double *)malloc(2*sizeof(double));
			polysAll[i].xmax = (double *)malloc(2*sizeof(double));

			rotate(x0, y0, +rx/2.0+x0, +ry/2.0+y0, &(polysAll[i].x[0]), &(polysAll[i].y[0]), angle, spherical);
			rotate(x0, y0, -rx/2.0+x0, +ry/2.0+y0, &(polysAll[i].x[1]), &(polysAll[i].y[1]), angle, spherical);
			rotate(x0, y0, -rx/2.0+x0, -ry/2.0+y0, &(polysAll[i].x[2]), &(polysAll[i].y[2]), angle, spherical);
			rotate(x0, y0, +rx/2.0+x0, -ry/2.0+y0, &(polysAll[i].x[3]), &(polysAll[i].y[3]), angle, spherical);

			polysAll[i].id      = i;
			polysAll[i].xmin[0] = x0;
			polysAll[i].xmax[0] = x0;
			polysAll[i].xmin[1] = y0;
			polysAll[i].xmax[1] = y0;
			for(j=0;j<polysAll[i].N;j++){
				polysAll[i].xmin[0] = MIN(polysAll[i].xmin[0], polysAll[i].x[j]);
				polysAll[i].xmax[0] = MAX(polysAll[i].xmax[0], polysAll[i].x[j]);
				polysAll[i].xmin[1] = MIN(polysAll[i].xmin[1], polysAll[i].y[j]);
				polysAll[i].xmax[1] = MAX(polysAll[i].xmax[1], polysAll[i].y[j]);
			}

			//fprintf(stderr,"\nxmin = %f \nxmax = %f \nymin = %f \nymax = %f\n",polysAll[i].xmin[0],polysAll[i].xmax[0],polysAll[i].xmin[1],polysAll[i].xmax[1]);

			i++;
		}



	}
	if(i==0){
		fprintf(stderr,"%s: 0 polygon found, check input file. Exiting...\n",MYNAME);
		exit(EXIT_FAILURE);
	}else{
		fprintf(stderr,"%d polygon(s) found\n",i);
	}

	double minArea;
	xmin[0] = polysAll[0].xmin[0];
	xmax[0] = polysAll[0].xmax[0];
	xmin[1] = polysAll[0].xmin[1];
	xmax[1] = polysAll[0].xmax[1];

	//fprintf(stderr,"xmin = %f \nxmax = %f \nymin = %f \nymax = %f\n",xmin[0],xmax[0],xmin[1],xmax[1]);
	if(0){
		minArea = (polysAll[0].xmax[0] - polysAll[0].xmin[0])*(sin(polysAll[0].xmax[1]*PI/180.0) - sin(polysAll[0].xmin[1]*PI/180.0))*180.0/PI;
	}else{
		minArea = (polysAll[0].xmax[0] - polysAll[0].xmin[0])*(polysAll[0].xmax[1] - polysAll[0].xmin[1]);
	}

	for(i=1;i<NpolysAll;i++){
		xmin[0] = MIN(xmin[0],polysAll[i].xmin[0]);
		xmax[0] = MAX(xmax[0],polysAll[i].xmax[0]);
		xmin[1] = MIN(xmin[1],polysAll[i].xmin[1]);
		xmax[1] = MAX(xmax[1],polysAll[i].xmax[1]);
		if(0){
			minArea += (polysAll[i].xmax[0] - polysAll[i].xmin[0])*(sin(polysAll[i].xmax[1]*PI/180.0) - sin(polysAll[i].xmin[1]*PI/180.0))*180.0/PI;
		}else{
			minArea += (polysAll[i].xmax[0] - polysAll[i].xmin[0])*(polysAll[i].xmax[1] - polysAll[i].xmin[1]);
		}


	}
	minArea /= 1.0*(double)NpolysAll;

	int SplitDim = 0, firstCall = 1;
	return createNode(polysAll,NpolysAll, minArea, SplitDim, xmin, xmax, firstCall);
}

Node *createNode(Polygon *polys, size_t Npolys, double minArea, int SplitDim, double xmin[2], double xmax[2], int firstCall){
	size_t i,j,SplitDim_new;

	/* 	Allocate memory for THIS node */
	Node *result = (Node *)malloc(sizeof(Node));
	static size_t countNodes, NpolysAll;
	static void   *root, *polysAll;

	Polygon *polysChild;

	/* 	Root & node id */
	if(firstCall){
		fprintf(stderr,"Building region mask tree...");

		root         = result;
		countNodes   = 0;
		polysAll     = polys;
		NpolysAll    = Npolys;



	}
	result->root      = root;
	result->id        = countNodes;
	result->Npolys    = Npolys;
	result->NpolysAll = NpolysAll;
	countNodes++;

	/*		Copy address of the complete polygon sample and
	 *		save ids of polygons inside the node
	 */
	result->polysAll     = polysAll;
	result->poly_id      = (int *)malloc(Npolys*sizeof(int));
	for(i=0;i<Npolys;i++){
		result->poly_id[i] = polys[i].id;
	}

	double area;
	if(0){
		area = (xmax[0] - xmin[0])*(sin(xmax[1]*PI/180.0) - sin(xmin[1]*PI/180.0))*180.0/PI;
	}else{
		area = (xmax[0] - xmin[0])*(xmax[1] - xmin[1]);
	}


	//double area = (xmax[0]-xmin[0])*(xmax[1]-xmin[1]);
	/* 	Leaf: either no polygon or cell smaller than minArea */
	if(result->Npolys == 0 || area < minArea) {
		result->type     = LEAF;
		result->Left     = NULL;
		result->Right    = NULL;
	}else{
		result->type       = NODE;
		result->SplitDim   = SplitDim;
		result->SplitValue = (xmax[result->SplitDim] + xmin[result->SplitDim])/2.0;

		/* 	Temporary data */
		polysChild = (Polygon *)malloc(Npolys*sizeof(Polygon));
		//for(i=0; i<Npolys; i++){
		//	polysChild[i].x    = (double *)malloc(NVERTICES*sizeof(double));
		//	polysChild[i].y    = (double *)malloc(NVERTICES*sizeof(double));
		//	polysChild[i].xmin = (double *)malloc(2*sizeof(double));
		//	polysChild[i].xmax = (double *)malloc(2*sizeof(double));
		//}

		double xminChild[2], xmaxChild[2];
		for(i=0;i<2;i++){
			xminChild[i] = xmin[i];
			xmaxChild[i] = xmax[i];
		}

		/* 	New splitDim for children */
		SplitDim++;
		if(SplitDim > 1)  SplitDim = 0;

		/*		"left" */

		/*		Set new right limit */
		xmaxChild[result->SplitDim] = result->SplitValue;
		j=0;
		for(i=0;i<Npolys;i++){
			if(polys[i].xmin[result->SplitDim] < result->SplitValue){
				cpyPolygonAddress(&polysChild[j],&polys[i]);
				j++;
			}
		}
		result->Left = createNode(polysChild,j,minArea,SplitDim,xminChild,xmaxChild,0);
		free(polysChild);

		/* 	Temporary data */
		polysChild = (Polygon *)malloc(Npolys*sizeof(Polygon));



		/*		"right" */

		/*		restore right limit and set new left limit */
		xmaxChild[result->SplitDim] = xmax[result->SplitDim];
		xminChild[result->SplitDim] = result->SplitValue;
		j=0;
		for(i=0;i<Npolys;i++){
			if(polys[i].xmax[result->SplitDim] > result->SplitValue){
				cpyPolygonAddress(&polysChild[j],&polys[i]);
				j++;
			}
		}
		result->Right = createNode(polysChild,j,minArea,SplitDim,xminChild,xmaxChild,0);

		free(polysChild);

		//free_Polygon(polysChild,Npolys);
	}

	result->Nnodes=countNodes;

	if (firstCall){

		fprintf(stderr,"Done.\n");

	}


	return result;
}

void free_Polygon(Polygon *polygon, size_t N){
	size_t i;
	for(i=0;i<N;i++){
		free(polygon[i].x);
		free(polygon[i].y);
		free(polygon[i].xmin);
		free(polygon[i].xmax);
	}
	free(polygon);
}

void free_Node(Node *node){
	if(node->type != LEAF){
		free_Node(node->Left);
		free_Node(node->Right);
	}
	free(node->poly_id);
	free(node);
	return;
}


void cpyPolygonAddress(Polygon *a, Polygon *b){
	/*		Copies b into a	*/
	size_t i;

	a->N  = b->N;
	a->id = b->id;

	a->x = b->x;
	a->y = b->y;
	a->xmin = b->xmin;
	a->xmax = b->xmax;
}


void cpyPolygon(Polygon *a, Polygon *b){
	/*		Copies b into a	*/
	size_t i;

	a->N  = b->N;
	a->id = b->id;
	for(i=0;i<NVERTICES;i++){
		a->x[i] = b->x[i];
		a->y[i] = b->y[i];
	}
	for(i=0;i<2;i++){
		a->xmin[i] = b->xmin[i];
		a->xmax[i] = b->xmax[i];
	}
}

/*
 *		Utils - numeric
 */

double distComo(double z, const double a[4]){
	/* Return the comoving distance of z, given the cosmological
	 * parameters a.
	 */

	int n = 1000;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc (n);

	gsl_function F;
	F.function = &drdz;
	F.params   = (void *)a;

	double result, error;
	gsl_integration_qag(&F, 0.0, z, 0, 1e-7, n, 6, w, &result, &error);
	gsl_integration_workspace_free (w);

	return result;
}

double drdz(double x, void * params){
	/* 	Inverse of E(z). With following cosmological parameters:
	 *	 	a[0] = H0;
	 * 	a[1] = Omega_M;
	 * 	a[2] = Omega_L;
	 * 	a[3] = c;
	 */

	double *a = (double *) params;
	return a[3]/(a[0]*sqrt(a[1]*pow(1+x,3.0)+a[2]));
}

double dvdz(double z, const double a[4]){
	/* 	Differential comoving volume per unit solid angle */
	return SQUARE(distComo(z, a))*a[3]/(a[0]*sqrt(a[1]*pow(1+z,3.0)+a[2]));
}


double distAngSpher(const double RA1, double DEC1, double RA2, double DEC2){
	/*		Returns the angular distance between points with RA,DEC (in degree)
	 */

	double sin2_ra  = 0.5*(1.0 - cos(RA1*PI/180.0)*cos(RA2*PI/180.0)  - sin(RA1*PI/180.0)*sin(RA2*PI/180.0));
	double sin2_dec = 0.5*(1.0 - cos(DEC1*PI/180.0)*cos(DEC2*PI/180.0)- sin(DEC1*PI/180.0)*sin(DEC2*PI/180.0));

	return 2.0*asin(sqrt(MAX(EPS/100.0, sin2_dec + cos(DEC1*PI/180.0)*cos(DEC2*PI/180.0)*sin2_ra)))*180.0/PI;

}

void rotate(double x0, double y0, double x, double y, double *xrot, double *yrot, double angle, int spherical){
	/* 	Takes positions (in degree if spherical = 1) and returns the rotated coordinates. */

	if(spherical){
		double sign  = x - x0 < 0.0 ? -1.0 : 1.0;
		*yrot = -sign*(distAngSpher(x0, y0, x, y0))*sin(angle) + (y-y0)*cos(angle)+y0;
		*xrot =  sign*(distAngSpher(x0, y0, x, y0))*cos(angle) + (y-y0)*sin(angle)+x0;
		*xrot =  2.0*asin(sin((*xrot-x0)*PI/180.0/2.0)/cos(*yrot*PI/180.0))*180.0/PI+x0;

	}else{
		*xrot = +(x-x0)*cos(angle) - (y-y0)*sin(angle)+x0;
		*yrot =  (x-x0)*sin(angle) + (y-y0)*cos(angle)+y0;
	}

}

gsl_rng *randomInitialize(size_t seed){
	/*		Random number generator.
	 *		Define here which type of random generator you want
	 */

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set (r,seed);

	return r;
}

double determineMachineEpsilon()
{
	double u, den;

	u = 1.0;
	do {
		u /= 2.0;
		den = 1.0 + u;
	} while(den>1.0);

	return(100.0 * u);
}

size_t determineSize_tError(){
	/*		Equals size_t_max in fact. */
	size_t count=0;
	count--;
	//count = (size_t)(-1);
	return count;
}

FILE *fopenAndCheck(const char *fileName,char *mode){
	/*		Checks if fileName exists and opens it. Exits otherwise.;	*/

	if (  !(strcmp(mode,"w")) && (!(strcmp(fileName,"")) || !(strcmp(fileName,"-")))){
		return stdout;
	}
	if (  !(strcmp(mode,"r")) && (!(strcmp(fileName,"")) || !(strcmp(fileName,"-")))){
		return stdin;
	}

	FILE *fileTmp = fopen(fileName,mode);

	if (fileTmp == NULL){
		fprintf(stderr,"%s: %s not found. Exiting...\n",MYNAME,fileName);
		exit(EXIT_FAILURE);
	}
	return fileTmp;
}

int getStrings(char *line, char *strings, char *delimit, size_t *N){
	/*		Extract each word/number in line separated by delimit and returns
	 *		the array of items in strings.
	 */
	int i,j,begin,length;

	if(line == NULL || line[0] == '\n' || line[0] == '#' || line[0] == '\0') return 0;

	i = 0;
	j = 0;
	while(line[i] != '\0' && line[i] != '#' && line[i] != '\n'){
		begin = i;
		while((line[i] == *delimit || line[i] == '\t') && (line[i] != '\0' || line[i] != '#' || line[i] != '\n')) i++;
		begin = i;
		while(line[i] != *delimit && line[i] != '\t' && line[i] != '\0' && line[i] != '#' && line[i] != '\n') i++;
		length = i - begin;
		if(length > 0){
			strncpy(strings+NCHAR*j,&line[begin],length);
			strcpy(strings+NCHAR*j+length,"\0");
			j++;
		}
	}

	(*N) = j;

	if(*N > 0){
		return 1;
	}else{
		return 0;
	}
}


void printCount(const size_t *count, const size_t *total, const size_t step){
	if((*count)%step == 0){
		fflush(stdout);
		fprintf(stderr,"\b\b\b\b%3.0f%%",100.0*(double)(*count)/(double)(*total));
	}
	return;
}

int checkFileExt(const char *s1, const char *s2){
	/*		Checks if s2 matches s1 extension */


	char *filterPtr = NULL;
	char s1Tmp[FILENAMESIZE];

	strcpy(s1Tmp, s1);

	/* 	look for filter string and remove it */
	filterPtr = strstr(s1Tmp, "[");
	if (filterPtr != NULL){
		*filterPtr = '\0';
	}

	int ext = strlen(s1Tmp) - strlen(s2);
	if(strcmp(s1Tmp+ext, s2) == 0){
		return 1;
	}else{
		return 0;
	}
}


int roundToNi(double a){
	/*		round a to the next integer */
	return (int)round(a);
}

int compareDoubles(const void *a,const void *b){
	/* 	Compares two double precision numbers */
	if (*(double *)a > *(double *)b)
	return 1;
	else if (*(double *)a < *(double *)b)
	return -1;
	else
	return 0;
}


char* concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the zero-terminato
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}
