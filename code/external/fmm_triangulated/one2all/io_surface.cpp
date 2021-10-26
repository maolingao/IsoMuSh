/***************************************************************************/
/* init_stencil build the triangulation for the geodesic eikonal 	   */
/* 	Fast Marching method This program was written by Ronny Kimmel      */
/***************************************************************************/
#include "y_stand.h" /* declaration file */
#include "y_stand2.h" /* declaration file */
#include "mex.h"
#include <stdlib.h>

/***************************************************************************/
/* Write the triangles to a file so they can be viewed by other programs.  */
/***************************************************************************/
void
WriteTriangles(Triangle *T, Vertex * V, int Tnum, char *filename)
// 	struct Triangle *T;	/* triangles, numbered 1..Tnum */
// 	struct Vertex  *V;	/* vertices coordinates 1..VNum */
// 	int             Tnum;
// 	char           *filename;
{
	int             i, j;
	FILE           *out_file = fopen(filename, "w");

	if (out_file == NULL) {
		fprintf(stderr, "Cannot open file %s for write - exiting\n",
		       filename);
		return;
	}
	for (i = 0; i < Tnum; i++){			/* For every triangle. */
		for (j = 0; j < 3; j++){		/* For every vertex in the triangle */
			fprintf(out_file, "%g %g %g ",
				V[T[i].Vind[j]].x,
				V[T[i].Vind[j]].y,
				V[T[i].Vind[j]].z);
		}
		fprintf(out_file, "\n");
	}
	fclose(out_file);

}
/***************************************************************************/
/* Write the edges for viewing by other programs.						   */
/* Thew edges are written to the file "filename".						   */
/***************************************************************************/
void
WriteGraph(Vertex  *V, int Vnum, char *filename)
// 	struct Vertex  *V;	/* vertices coordinates 1..VNum */
// 	int             Vnum;
// 	char           *filename;
{
	int             i, j;
	FILE           *out_file = fopen(filename, "w");

	if (out_file == NULL) {	/* If the file cannot be opended. */
		fprintf(stderr, "Cannot open file %s for write - exiting \n",
		       filename);
		return;
	}
	for (i = 0; i < Vnum; i++) {	/* For every vertex */
		for (j = 0; j < V[i].vn; j++){ /* For every neighbour of the vertex */
									/* Write the coordinates of the vertex  */
			fprintf(out_file, "%g %g %g ",V[i].x,V[i].y,V[i].z);
									/* and of its neighbour.				*/
			fprintf(out_file, "%g %g %g ",
			V[V[i].VN[j]].x,
			V[V[i].VN[j]].y,
			V[V[i].VN[j]].z);
			fprintf(out_file, "\n");
		}
	}
	fclose(out_file);

}
/***************************************************************************/
/* Write the 9 coordinates of each triangle, but replace the Z coordinate  */
/* the U value, for viewing by other programs.							   */
/***************************************************************************/
void
WriteUTriangles(Triangle *T, Vertex * V, int Tnum, char *filename)
// 	struct Triangle *T;	/* triangles, numbered 1..Tnum	*/
// 	struct Vertex  *V;	/* vertices coordinates 1..VNum */
// 	int             Tnum;
// 	char           *filename;
{
	int             i, j;
	FILE           *out_file = fopen(filename, "w");

	if (out_file == NULL) {		/* If the file cannot be opened. */
		fprintf(stderr, "Cannot open file %s for write - exiting \n",
		       filename);
		return;
	}

	for (i = 0; i < Tnum; i++) { /* For every triangle.			 */
		for (j = 0; j < 3; j++){ /* For every vertex in the triangle. */
			fprintf(out_file, "%g %g %g ", /* Write coordinates of the vertex. */
				V[T[i].Vind[j]].x, V[T[i].Vind[j]].y,
				V[T[i].Vind[j]].U);
		}
		fprintf(out_file, "\n");
	}
	fclose(out_file);

}
/***************************************************************************/
/* Write the avarage U value for each traingle to be viewed by other	   */
/* programs.															   */
/***************************************************************************/
void
WriteU(Triangle *T, Vertex * V, int Tnum, char *filename)
// 	struct Triangle *T;	/* triangles, numbered 1..Tnum */
// 	struct Vertex  *V;	/* vertices coordinates 1..VNum */
// 	int             Tnum;
// 	char           *filename;
{
	int             i, j;
	FILE           *out_file = fopen(filename, "w");
	double		u;

	if (out_file == 0) {
		fprintf(stderr, "Cannot open file %s for write - exiting \n",
		       filename);
		return;
	}
	
	for (i = 0; i < Tnum; i++) {		/* For every vertex */
		u = 0;							/* Sum up the U values / 3 */
		for (j = 0; j < 3; j++)			/* For every vertex in the triangle. */
			u += V[T[i].Vind[j]].U/3.0;	/* Add the next vertex to the avarage. */
		fprintf(out_file, "%g ", u);
	}
	fclose(out_file);

}
/***************************************************************************/
/*Add2Path append a segment in 3D to the path for mathematica usage   */
/***************************************************************************/
void
Add2Path(double xi,double yi,double zi,double xo,double yo,double zo, FILE *out_file)
// 	double	xi,yi,zi,xo,yo,zo;	
// 	FILE           *out_file;
{
	double		epsi = 0.03/((double)(XSIZE*XSIZE));

	/* shift the segment in 4 dirrections, so that it is shown on the 
		3D object */
	fprintf(out_file, "%g %g %g %g %g %g  ", xi+epsi,yi,zi, xo+epsi,yo,zo);
	fprintf(out_file, "%g %g %g %g %g %g  ", xi-epsi,yi,zi, xo-epsi,yo,zo);
	fprintf(out_file, "%g %g %g %g %g %g  ", xi,yi+epsi,zi, xo,yo+epsi,zo);
	fprintf(out_file, "%g %g %g %g %g %g  ", xi,yi-epsi,zi, xo,yo-epsi,zo);
	fprintf(out_file, "%g %g %g %g %g %g  ", xi,yi,zi+epsi, xo,yo,zo+epsi);
	fprintf(out_file, "%g %g %g %g %g %g \n", xi,yi,zi-epsi, xo,yo,zo-epsi);
}

/***************************************************************************/
/*OpenPath for Viewing by other programs.								   */
/***************************************************************************/
FILE
*OpenPath(char *filename)
// 	char           *filename;
{
	FILE           *out_file = fopen(filename, "w");
	if (out_file == NULL) {	/* If the opening fail. */
		fprintf(stderr, "Cannot open file %s for write - exiting \n",
		       filename);
		exit(1);
	}

	fprintf(out_file, " ");
	return (out_file);
}

/***************************************************************************/
/*ClosePath colse pathes opended by OpenPath.							   */
/***************************************************************************/
void
ClosePath(FILE *out_file)
{
	fclose(out_file);
}

/***************************************************************************/
/*AddAPoint append a vertex `sead' point in 3D for mathematica usage       */
/***************************************************************************/
void
AddPoint(double xi,double yi,double zi, FILE *out_file)
// 	double	xi,yi,zi;	
// 	FILE           *out_file;
{
	double		epsi = 0.03/((double)XSIZE * XSIZE);

	/* shift the segment in 4 dirrections, so that it is shown on the 
		3D object */
	fprintf(out_file, "%g %g %g  ", xi+epsi,yi,zi);
	fprintf(out_file, "%g %g %g  ", xi,yi+epsi,zi);
	fprintf(out_file, "%g %g %g  ", xi,yi,zi+epsi);
	fprintf(out_file, "%g %g %g  ", xi-epsi,yi,zi);
	fprintf(out_file, "%g %g %g  ", xi,yi-epsi,zi);
	fprintf(out_file, "%g %g %g  ", xi,yi,zi-epsi);
}

/************************************************************************************************/
/* WriteMatrix: This function prints the matrix it gets into a file whose name is also given	*/
/* to the fucntion. The function also gets the size of each dimension of the matrix.			*/
/************************************************************************************************/

void WriteMatrix(double **matrix, const char *filename, int height, int width){

	FILE *output_file;								/* File pointer to the file where the matrix
													   will be printed.							*/
	int i, j;										/* Indexes for the matrix.					*/

	output_file = fopen(filename, "w");
	if(output_file == NULL){
		fprintf(stderr, "Cannot open file %s for writing\n", filename);
		exit(-1);
	}

	for(i = 0; i < height; i++){					/* Print the matrix.						*/
		for(j = 0; j < width; j++)
			fprintf(output_file, "%-15lg", MIN(matrix[i][j], matrix[j][i]));
		fprintf(output_file, "\n\n");
	}
	fclose(output_file);
}

/***************************************************************************/
/* WriteManifold  x,y,z and U in a row scan (by Alon Spira)				   */
/***************************************************************************/
void
WriteManifold(struct Vertex *V, int Vnum, char *filename)
// 	struct Vertex  *V;	/* vertixes coordinates 1..VNum */
// 	int				Vnum; /* The number of vertixes in V. */
// 	char           *filename;
{
	int            ind;			/* Index of a vertex. */
	double		   *b_ind;		/* The current place in the buffer for the buffer. */
	FILE           *out_file = fopen(filename, "wb");
	double		   *buf;		/* Here all the picture data is saved. */

	if (out_file == NULL) {
		fprintf(stderr, "Cannot open file %s for write - exiting \n",
		       filename);
		exit(1);
	}

	b_ind=buf=(double *)malloc(4*Vnum*sizeof(double)); /* Allocating the buffer */
	if(buf == NULL){
		fprintf(stderr, "ERROR: Not enoguh memory for the output buffer.\n");
		exit(1);
	}

	/* Collecting the data about all the vertexes into the buffer. */
	for (ind = 0; ind < Vnum; ind++, b_ind++){
		*b_ind = V[ind].U;
	}
	for (ind = 0; ind < Vnum; ind++, b_ind++){
		*b_ind = V[ind].x;
	}
	for (ind = 0; ind < Vnum; ind++, b_ind++){
		*b_ind = V[ind].y;
	}
	for (ind = 0; ind < Vnum; ind++, b_ind++){
		*b_ind = V[ind].z;
	}
	/* Write the result. */
	fwrite((void *)buf, sizeof(double), 4*Vnum, out_file); /* Writing the data into the file */
	fclose(out_file);
	free(buf);
}

/************************************************************************************************/
/* This function gets the matrix with the splitting data and prints it into a file with the		*/
/* name it gets. This data should be printed if and only if the FULL_OUTPUT flag is on and		*/
/* therefore the function exists only in this case.												*/
/************************************************************************************************/

#ifdef FULL_OUTPUT

void print_split_data(split_data **Array, char filename[]){

	FILE *output_file;								/* File pointer to the file where the array
													   will be printed.							*/
	int *buffer;									/* The buffer where the output will be
													   writen.									*/
	int *buffer_ind;								/* Point to the next place in the buffer.	*/
	int i, j;										/* Indexes for the array.					*/

	output_file = fopen(filename, "wb");
	if(output_file == NULL){						/* If the file could not be opened.			*/
		fprintf(stderr, "Cannot open file %s for writing\n", filename);
		exit(-1);
	}

	buffer = malloc( 6 * (SIZE) * sizeof(int));
													/* Allocate memory for the buffer for
													   printing the data.						*/
	if(buffer == NULL){
		fprintf(stderr, "Not enough memory for output buffer\n");
		exit(-1);
	}

	buffer_ind = buffer;
	for(i = 0; i < XSIZE; i++)						/* Print the array into the buffer.			*/
		for(j = 0; j < YSIZE; j++){
			*(buffer_ind++) = Array[i][j].virtual_neighbour1_X;
			*(buffer_ind++) = Array[i][j].virtual_neighbour1_Y;
			*(buffer_ind++) = Array[i][j].virtual_neighbour2_X;
			*(buffer_ind++) = Array[i][j].virtual_neighbour2_Y;
			*(buffer_ind++) = Array[i][j].virtual_neighbour3_X;
			*(buffer_ind++) = Array[i][j].virtual_neighbour3_Y;
		}
													/* Print the buffer to the file.			*/
	fwrite(buffer, sizeof(int), 6 * (SIZE - 2 * XSIZE - 2 * YSIZE + 4) , output_file);
	fclose(output_file);


}

#endif

/************************************************************************************************/
/* This function gets the array with the update data and prints it into a file with the			*/
/* name it gets. This data should be printed if and only if the FULL_OUTPUT flag is on and		*/
/* therefore the function exists only in this case.												*/
/* ArrayEnd is the next free record in the array, the array is printed up to this place.		*/
/************************************************************************************************/

#ifdef FULL_OUTPUT

void print_update_data(update_data *Array, update_data *ArrayEnd, char filename[]){

	FILE *output_file;								/* File pointer to the file where the array
													   will be printed.							*/
	double *buffer;									/* The buffer where the output will be
													   writen.									*/
	double *buffer_ind;								/* Point to the next place in the buffer.	*/
	update_data *ArrayInd;							/* Point to a place in Array.				*/

	output_file = fopen(filename, "wb");
	if(output_file == NULL){						/* If the file could not be opened.			*/
		fprintf(stderr, "Cannot open file %s for writing\n", filename);
		exit(-1);
	}

	buffer = malloc( 9 * (ArrayEnd - Array) * sizeof(double));
													/* Allocate memory for the buffer for
													   printing the data.						*/
	if(buffer == NULL){
		fprintf(stderr, "Not enough memory for output buffer\n");
		exit(-1);
	}

	buffer_ind = buffer;							/* Print the array into the buffer.			*/
	for(ArrayInd = Array; ArrayInd != ArrayEnd; ArrayInd++){
		*(buffer_ind++) = ArrayInd->NewAliveX;
		*(buffer_ind++) = ArrayInd->NewAliveY;
		*(buffer_ind++) = ArrayInd->UpdatedX;
		*(buffer_ind++) = ArrayInd->UpdatedY;
		*(buffer_ind++) = ArrayInd->ThirdX;
		*(buffer_ind++) = ArrayInd->ThirdY;
		*(buffer_ind++) = ArrayInd->NewU;
#ifndef NO_OLD_U_DATA
		*(buffer_ind++) = ArrayInd->NewAliveU;
		*(buffer_ind++) = ArrayInd->ThirdU;
#endif

	}
													/* Print the buffer to the file.			*/
#ifndef NO_OLD_U_DATA
	fwrite(buffer, sizeof(double), 9 * (ArrayEnd - Array) , output_file);
#else
	fwrite(buffer, sizeof(double), 7 * (ArrayEnd - Array) , output_file);
#endif
	fclose(output_file);


}

#endif

/********************************************************************************************/
/* This function read the traingles and vertixes from the input matrix and accordingly fill */
/* in T, V, Tnum and Vnum (T and V are also allocated accordingly) by this function.		*/
/********************************************************************************************/

// void ReadSurface(struct Vertex *(V[]), struct Triangle *(T[]), int *Vnum, int *NonSplitTnum, int *Tnum, const double *CoordinatesData){
void ReadSurface(struct Vertex **V, struct Triangle **T, int *Vnum, int *NonSplitTnum, int *Tnum, const double *CoordinatesData){

	int i;											/* Index for the loop.					*/

	*Vnum = (int)*(CoordinatesData++);				/* Read Tnum and Vnum.					*/
	*NonSplitTnum = (int)*(CoordinatesData++);
	*Tnum = *NonSplitTnum * 3;
													/* Allocate memory for both triangles and
												       vertixes.							*/
    *T = (struct Triangle *) malloc(sizeof(struct Triangle) * *Tnum);
    if(*T == NULL){
		fprintf(stderr, "Out of memory for triangles - exiting.\n");
		exit(-1);
	}
	*V = (struct Vertex *)   malloc(sizeof(struct Vertex) * *Vnum);
    if(*V == NULL){
		free(T);
		fprintf(stderr, "Out of memory for vertiexes - exiting.\n");
		exit(-1);
	}

													/* Move the data to V and T.			*/
	for(i = 0; i < *Vnum; i++){
		(*V)[i].x = *(CoordinatesData++);
		(*V)[i].y = *(CoordinatesData++);
		(*V)[i].z = *(CoordinatesData++);
	}
	for(i = 0; i < *NonSplitTnum; i++){
		(*T)[i].Vind[0] = (int)*(CoordinatesData++);
		(*T)[i].Vind[1] = (int)*(CoordinatesData++);
		(*T)[i].Vind[2] = (int)*(CoordinatesData++);
	}

	for(i = 0; i < *Vnum; i++){						/* Add every triangle to its vertixes.	*/
		(*V)[i].ti = 0;
	}
	for(i = 0; i < *NonSplitTnum; i++){
		if((*T)[i].Vind[0] != *NonSplitTnum){
			(*V)[(*T)[i].Vind[0]].Tind[(*V)[(*T)[i].Vind[0]].ti++] = i;
			(*V)[(*T)[i].Vind[1]].Tind[(*V)[(*T)[i].Vind[1]].ti++] = i;
			(*V)[(*T)[i].Vind[2]].Tind[(*V)[(*T)[i].Vind[2]].ti++] = i;
		}
	}
	
	return;

}

/********************************************************************************************/
/* This function is the equivalent of ReadSurface for the case when the program gets		*/
/* coordinates vectors and not a surface array. This function allocate the memory for the	*/
/* vertexes and the triangles, fill the coordinates of the vertexes and fill the non split	*/
/* triangles.																				*/
/********************************************************************************************/

void CombineVectors(struct Vertex *(V[]), struct Triangle *(T[]), int *Vnum, int *NonSplitTnum, int *Tnum, const mxArray *plhs[]){

	int i, j;										/* Indexes for the loops.				*/
	const double *CoordinatesData[3];				/* Pointer to the data in the three 
													   vector arrays.						*/
	const int *TrianglesData;						/* Pointer to the triangles data.		*/


	*Vnum = mxGetNumberOfElements(plhs[2]);			/* Read Tnum and Vnum.					*/
	*NonSplitTnum = mxGetNumberOfElements(plhs[0]) / 3;
	*Tnum = *NonSplitTnum * 3;

													/* Allocate memory for both triangles and
												       vertixes.							*/
    *T = (struct Triangle *) malloc(sizeof(struct Triangle) * *Tnum);
    if(*T == NULL){
		fprintf(stderr, "Out of memory for triangles - exiting.\n");
		exit(-1);
	}
	*V = (struct Vertex *)   malloc(sizeof(struct Vertex) * *Vnum);
    if(*V == NULL){
		free(T);
		fprintf(stderr, "Out of memory for vertiexes - exiting.\n");
		exit(-1);
	}

	for(i = 0; i < 3; i++)							/* Get pointers to the data in the vectors arrays.	*/
		CoordinatesData[i] = (const double *)mxGetData(plhs[2 + i]);
	TrianglesData = (const int *)mxGetData(plhs[0]);				/* Get the data pointer of the triangles array.		*/
	
	for(i = 0; i < *Vnum; i++){						/* Move the data to V and T.			*/
		(*V)[i].x = *(CoordinatesData[0]++);
		(*V)[i].y = *(CoordinatesData[1]++);
		(*V)[i].z = *(CoordinatesData[2]++);
	}
	for(i = 0; i < 3; i++)
		for(j = 0; j < *NonSplitTnum; j++){
			(*T)[j].Vind[i] = *(TrianglesData++);
		}
		
	for(i = 0; i < *Vnum; i++){						/* Add every triangle to its vertixes.	*/
		(*V)[i].ti = 0;
	}
	/* Can be greatly improved! */
	for(i = 0; i < *NonSplitTnum; i++){
		if((*T)[i].Vind[0] != *NonSplitTnum){
			(*V)[(*T)[i].Vind[0]].Tind[(*V)[(*T)[i].Vind[0]].ti++] = i;
			(*V)[(*T)[i].Vind[1]].Tind[(*V)[(*T)[i].Vind[1]].ti++] = i;
			(*V)[(*T)[i].Vind[2]].Tind[(*V)[(*T)[i].Vind[2]].ti++] = i;
		}
	}
	
	return;

}