/***************************************************************************/
/* Triangle.c : test fast marching Eikonal on acute (non obtuse)           */ 
/* triangulated grids	   											  	   */
/* This program was written by Ronny Kimmel    						       */
/***************************************************************************/

// #include "windows.h"


#include "y_stand.h"									/* declaration file */
#include "mex.h"
#include <stdlib.h>
#include "string.h"

#define LastUpdate "25-Nov-2005"						/* The last update
														   of the program.	*/
extern void WriteManifold(struct Vertex V[], int Vnum, char file_name[]);
extern void WriteMatrix(double **matrix, const char *filename, int height, int width);
extern void ReadSurface(struct Vertex **V, struct Triangle **T, int *Vnum, int *NonSplitTnum, int *Tnum, const double CoordinatesData[]);
extern void CombineVectors(struct Vertex *(V[]), struct Triangle *(T[]), int *Vnum, int *NonSplitTnum, int *Tnum, const mxArray *plhs[]);

extern void gaussj(double a[7][7], double b[7]);

extern void		init_heap(int Vnum);
extern void		free_heap(void);
extern void		SetSourceSet(int xCoordinates[], int yCoordinates[], int points_no, struct Vertex VertexesArray[], int Vnum);
extern void   	_fastcall MarchAlg(struct Vertex  *V, int Vnum, int start_from); 
extern double 	z_function(double x,double y); 
extern double 	f_potential(double x,double y); 

extern void LevelSets(int Ci, struct Triangle *T, struct Vertex  *V,
	int Tnum, int Vnum, char FileName[]);

extern void WriteTriangles(struct Triangle *T,struct Vertex  *V, 
	int Tnum, char *filename);

extern void WriteGraph(struct Vertex  *V, int Vnum, char *filename);


extern void WriteUTriangles(struct Triangle *T,struct Vertex  *V, 
	int Tnum, char *filename);

extern void WriteU(struct Triangle *T,struct Vertex  *V, 
	int Tnum, char *filename);

extern void WriteTV(struct Triangle *T,struct Vertex  *V, 
	int Tnum, int Vnum);

extern FILE * OpenPath(char   *filename);
extern void ClosePath(FILE *file_out);
extern void write_vertices_2_pgm( char *filename, struct Vertex  *V,
							int xs,int ys);

#if defined(LINEAR_GEODESIC) || defined(SMOOTH_GEODESIC)
extern void ConicGrad(struct Triangle *T,struct Vertex  *V, 
	int Tnum, int Vnum, int IllegalVertex);
extern void FindSmoothGeodesic(int vi, struct Triangle *T, struct Vertex  *V,
	int Tnum, int Vnum, FILE* out_file);
extern void FindGeodesic(int vi, struct Triangle *T, struct Vertex  *V,
	int Tnum, int Vnum, FILE *out_file);
#endif

#ifdef FULL_OUTPUT

extern void print_split_data(split_data **Array, char filename[]);
extern void print_update_data(update_data *Array, update_data *ArrayEnd, char filename[]);

#endif



void switch_priority(boolean TOP){
	
	static int process_priority;	/* Storge for the normal priority	*/
	static int thread_priority;
	
	return;


	if(TOP){
		/* Save the normal state.										*/
		process_priority = GetPriorityClass(GetCurrentProcess());
		thread_priority = GetThreadPriority(GetCurrentThread());

		/* Increase the priority of the application and the thread.		*/
//		SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
//		SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL);
	
	}
	else{
		/* Restore the normal priority.									*/
		SetPriorityClass(GetCurrentProcess(), process_priority);
		SetThreadPriority(GetCurrentThread(), thread_priority);
		
	}
}


/***************************************************************************/
/* L2  error the diff between the distance and U*/
/***************************************************************************/
double
L2(V, Vnum)
	struct Vertex  *V;	/* vertices coordinates 1..VNum */
	int             Vnum;
{
	int             i;
	double          d, x, y, z;
	double          sum = 0, dxdy = 1.0/((double)Vnum);
	/* find the zero point */
	for (i = 0; i < Vnum; i++)
		if (V[i].U < 0.00000001) {
			x = V[i].x;
			y = V[i].y;
			z = V[i].z;
//			mexPrintf ("Source at: (%g,%g,%g)\n",x,y,z);
		}
	for (i = 0; i < Vnum; i++) {
		d = sqrt(Sqr(x - V[i].x) + Sqr(y - V[i].y) + Sqr(z - V[i].z));
		sum += Sqr(V[i].U - d)*dxdy; 
	}
	return (sqrt(sum));
}
/***************************************************************************/
/* L1  error norm the diff between the distance and U*/
/***************************************************************************/
double
L1(V, Vnum)
	struct Vertex  *V;	/* vertices coordinates 1..VNum */
	int             Vnum;
{
	int             i;
	double          d, x, y, z;
	double          sum = 0.0, dxdy = 1.0/((double)Vnum);
	/* find the zero point */
	for (i = 0; i < Vnum; i++)
		if (V[i].U < 0.00000001) {
			x = V[i].x;
			y = V[i].y;
			z = V[i].z;
		}
	for (i = 0; i < Vnum; i++) {
		d = sqrt(Sqr(x - V[i].x) + Sqr(y - V[i].y) + Sqr(z - V[i].z));
		sum += ABS(V[i].U - d)*dxdy;
	}
	return (sum);
}

/***************************************************************************/
/* CosAngle the cos of the angle at the vertex v0, between v1 and v2	   */
/***************************************************************************/
double
CosAngle(int v0,int v1,int v2,struct Vertex *V)
{
	double x1,x2,y1,y2,z1,z2,res;
	if(v0 != -1 and v1 != -1 and v2 != -1){
		x1 = V[v1].x - V[v0].x;
		x2 = V[v2].x - V[v0].x;
		y1 = V[v1].y - V[v0].y;
		y2 = V[v2].y - V[v0].y;
		z1 = V[v1].z - V[v0].z;
		z2 = V[v2].z - V[v0].z;
		res = x1*x2+y1*y2+z1*z2;		/* dot product */
		res /= sqrt(x1*x1+y1*y1+z1*z1); /* normalize */
		res /= sqrt(x2*x2+y2*y2+z2*z2);
		return(res);
	}
	else
		return 0;
}
/***************************************************************************/
/* Length between the vertex v0 and v1									   */
/***************************************************************************/
double
Length(int v0,int v1,struct Vertex *V)
{
	double x1,y1,z1,res;
	if(v0 != -1 and v1 != -1){
		x1 = V[v1].x - V[v0].x;
		y1 = V[v1].y - V[v0].y;
		z1 = V[v1].z - V[v0].z;
		res = sqrt(x1*x1+y1*y1+z1*z1); /* distance */
		return(res);
	}
	else
		return MAX_DISTANCE;
}
/***************************************************************************/
/* nextT next triangle to be unfolded. find the triangle #, and the vertex */
/* Returns true is the next triangle was found.							   */
/* v1 and v2 indicate the edge that is common to the original triangle	   */
/* and the triangle to be unfolded.										   */
/* ti is the original triangle and v3 is the other vertex of the triangle  */
/* to be unfolded.														   */
/* vn is the index of the triangle to be unfolded.						   */
/***************************************************************************/
boolean
nextT(int ti, int v1, int v2, int *v3, int *tn,
      struct Triangle * T, struct Vertex * V, int Tnum)
{
	boolean				found = FALSE;	/* Indicates whether we found the next triangle. */
	int					i,				/* Index for the loop.							 */
						tj,				/* A candidate tp be the next triangle.			 */
						k;				/* Index for the inner loop.					 */
	/* scan every triangle of vi */
	for (i = 0; i < V[v1].ti and not found; i++) {
		tj = V[v1].Tind[i];
		if (tj < Tnum and tj != ti)
			/* search for tj in the list of v2 */
			for (k = 0; k < V[v2].ti and not found; k++) 
				if (V[v2].Tind[k] == tj && !T[tj].Split) {
					found = TRUE;
					*tn = tj;
				}
	}
	if (found){ /* find v3, the other vertex in the triangle to be unfolded.			 */
		if(T[*tn].Vind[0] == v1){
			if(T[*tn].Vind[1] == v2)
				*v3 = T[*tn].Vind[2];
			else
				*v3 = T[*tn].Vind[1];
		}
		else if(T[*tn].Vind[1] == v1){
			if(T[*tn].Vind[0] == v2)
				*v3 = T[*tn].Vind[2];
			else
				*v3 = T[*tn].Vind[0];
		}
		else{
			if(T[*tn].Vind[0] == v2)
				*v3 = T[*tn].Vind[1];
			else
				*v3 = T[*tn].Vind[0];
		}
	}
	return (found);
}

/***************************************************************************/
/* Split obtuse angles by unfolding splitting and connecting, return the   */
/* number of unfoldings that were nessesery.							   */
/* ti is the tirnalge to be splittined, V0 is the vertex with the obtuse   */
/* angle while V1 and V2 are the other vertexes of ti.					   */
/***************************************************************************/
int
Split(int ti, int V0, int V1, int V2, struct Triangle * T, struct Vertex * V,
      int NonSplitTnum, int Vnum)
{
	double          xv1,x1, y1, yv1,
					x2,				/* The distance between V0 and V2 */
					y2,
					x3,y3, xt2, xt3, yt3,
					e0,				/* The distance between v1 and v2 */
					e1,				/* The distance between v2 and v3.*/
					e2,				/* The distance between v0 and v1 */
					ta,				/* Tan of alpha.				  */
					cb, sb;			/* Cos and Sin of beta.			  */
	int             v1 = V1, v2 = V2, v3,
					tm=ti,			/* The current triangle we are
									   working on.					  */
					tn,				/* The triangle returned by NextT */
					count = 0;		/* The number of triangles unfolded
									   so far.						  */
									/* Becomes true when the split was done */
	boolean         splitter = FALSE;
	x2 = Length(V0, V2, V);
	y2 = 0;
	e0 = Length(V1, V2, V);
	e2 =  Length(V0, V1, V);
	xv1 = x1 = (x2 * x2 + e2 * e2 - e0 * e0) / (2.0 * x2);/* translation */
	yv1 = y1 = sqrt(e2 * e2 - x1 * x1);
	ta = -x1 / y1;		/* tan (alpha) in Fig. 1 */
	/* if there is a next triangle and not splited */
	while (nextT(tm, v1, v2, &v3, &tn, T, V, NonSplitTnum) and (not splitter) ) {
		count++;
		tm = tn;		/* Update the wording triangle. */
		cb = (x2 - x1) / sqrt(Sqr(x2 - x1) + Sqr(y2 - y1));	/* cos beta */
		sb = sqrt(1 - cb * cb);								/* sin beta */
		if (y2 < y1)	/* Adjast the sign of SIN(beta).				*/
			sb *= -1.0;
		xt2 = Length(v1, v2, V);
		e1 = Length(v2, v3, V);
		e2 = Length(v1, v3, V);
		xt3 = (xt2 * xt2 + e2 * e2 - e1 * e1) / (2.0 * xt2);
		yt3 = sqrt(e2 * e2 - xt3 * xt3);
		x3 = cb * xt3 - sb * yt3 + x1;
		y3 = sb * xt3 + cb * yt3 + y1;

		if (x3 > 0 and y3/x3 > ta) {		/* if we found a splitter */
			splitter = TRUE;
											/* Add the stencils involving the
											   splitting edge.			*/
			V[V0].ST[V[V0].si].Ctheta = (x3*xv1+y3*yv1)/
				sqrt((xv1*xv1+yv1*yv1)*(x3*x3+y3*y3));
			V[V0].ST[V[V0].si].v1 = V1;
			V[V0].ST[V[V0].si].v2 = v3;
			V[V0].ST[V[V0].si].l1 = Length(V0, V1, V);
			if(V[V0].si == MaxTri - 1)
			{ // mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n"); 
			}
			else
				V[V0].ST[V[V0].si++].l2 = sqrt(x3*x3+y3*y3);


			V[V0].ST[V[V0].si].Ctheta = x3/sqrt(x3*x3+y3*y3);
			V[V0].ST[V[V0].si].v1 = v3;
			V[V0].ST[V[V0].si].v2 = V2;
			V[V0].ST[V[V0].si].l1 =sqrt(x3*x3+y3*y3);
			if(V[V0].si == MaxTri - 1)
			{ // mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V0].ST[V[V0].si++].l2 =Length(V0, V2, V);

			if(V[v3].vn == 3 * MaxTri - 1)
			{
				//mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[v3].VN[V[v3].vn++] = V0; /* add dirrectional edge
											  to v3					*/

			T[NonSplitTnum + (ti * 2)].Vind[0] = V1;	/* Add the triangles of the splitting. */
			T[NonSplitTnum + (ti * 2)].Vind[1] = V0;
			T[NonSplitTnum + (ti * 2)].Vind[2] = v3;
			T[NonSplitTnum + (ti * 2)].Split = TRUE;
			if(V[V1].ti == MaxTri - 1)
			{
				//mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V1].Tind[V[V1].ti++] = NonSplitTnum + (ti * 2);
			if(V[V0].ti == MaxTri - 1)
			{
				//mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V0].Tind[V[V0].ti++] = NonSplitTnum + (ti * 2);
			if(V[v3].ti == MaxTri - 1)
			{
				//mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[v3].Tind[V[v3].ti++] = NonSplitTnum + (ti * 2);
			T[NonSplitTnum + (ti * 2) + 1].Vind[0] = v3;
			T[NonSplitTnum + (ti * 2) + 1].Vind[1] = V0;
			T[NonSplitTnum + (ti * 2) + 1].Vind[2] = V2;
			T[NonSplitTnum + (ti * 2) + 1].Split = TRUE;
			if(V[v3].ti == MaxTri - 1)
			{
				//mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[v3].Tind[V[v3].ti++] = NonSplitTnum + (ti * 2) + 1;
			if(V[V0].ti == MaxTri - 1)
			{
				//mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V0].Tind[V[V0].ti++] = NonSplitTnum + (ti * 2) + 1;
			if(V[V2].ti == MaxTri - 1)
			{
				//mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V2].Tind[V[V2].ti++] = NonSplitTnum + (ti * 2) + 1;


#ifdef FULL_OUTPUT						/* Add the split data to the split file. */
			if(SplitDataArray[V[V0].u1][V[V0].u2].virtual_neighbour1_X == NoData){
				SplitDataArray[V[V0].u1][V[V0].u2].virtual_neighbour1_X = V[v3].u1;
				SplitDataArray[V[V0].u1][V[V0].u2].virtual_neighbour1_Y = V[v3].u2;
			}
			else if(SplitDataArray[V[V0].u1][V[V0].u2].virtual_neighbour2_X == NoData){
				SplitDataArray[V[V0].u1][V[V0].u2].virtual_neighbour2_X = V[v3].u1;
				SplitDataArray[V[V0].u1][V[V0].u2].virtual_neighbour2_Y = V[v3].u2;
			}
			else{
				SplitDataArray[V[V0].u1][V[V0].u2].virtual_neighbour3_X = V[v3].u1;
				SplitDataArray[V[V0].u1][V[V0].u2].virtual_neighbour3_Y = V[v3].u2;
			}
#endif
		}
		else {						   /* we have not found a splitter,
										  continue unfolding			*/
			if (x3 < 0){
				v1 = v3; x1 = x3; y1 = y3;
			} else {
				v2 = v3; x2 = x3; y2 = y3;
			}
		}	
	}
	return(count);				/* Return the number of triangles that were
								   unfolded.							   */
}

/***************************************************************************/
/* InitGraph init the numerical graph V, unfold in order to split obtuse   */
/* triangles.															   */
/***************************************************************************/
void
InitGraph(struct Triangle * T, struct Vertex * V, int Vnum, int NonSplitTnum, int Tnum)
{
	int             i, k,			/* Indexes for the loops.				*/
#if defined(SMOOTH_GEODESIC) || defined(LINEAR_GEODESIC)
					j,
#endif
					ti,				/* Index for the current triangle.		*/
					v1,v2,			/* Indexes for the neighbours in the
									   triangle of the the current vertex.	*/
					count= 0,		/* The number of unfolding in one split.*/
					mcount=0;		/* The maximum value count got.			*/
	boolean         found;			/* Used for adding the neighbours to
									   every triangle.						*/
	double          ca;				/* The cosin value that determine
									   whether this triangle is obtuse.		*/
	struct Stencil *p_st;			/* Pointer to a stencil.				*/
	int ind, si_count;				/* Used for precalculting values for the
									   stencil.								*/


	/* Initialize the vertixes. */
	for (i = 0; i < Vnum; i++) {	/* zero counters of all vertices */
		V[i].si = 0;				/* # of connections to other triangles */
		V[i].vn = 0;				/* # of connections to other vertices */
	}

	/* Set the split field of the triangles that exist now, before the splitting to false.	*/
	for(i = 0; i < NonSplitTnum; i++)
		T[i].Split = FALSE;
	
	for (i = 0; i < Vnum; i++) {		/* scan all vertices */
		for (ti = 0; ti < V[i].ti; ti++){/* scan connected triangles */
			if (V[i].Tind[ti] < NonSplitTnum) {	/* if valid triangle */
												/* Make v1 and v2 the neighbours.			*/
				if (T[V[i].Tind[ti]].Vind[0] == i){
					v1 = T[V[i].Tind[ti]].Vind[1];
					v2 = T[V[i].Tind[ti]].Vind[2];
				}
				else if (T[V[i].Tind[ti]].Vind[1] == i){
					v1 = T[V[i].Tind[ti]].Vind[2];
					v2 = T[V[i].Tind[ti]].Vind[0];
				}
				else if (T[V[i].Tind[ti]].Vind[2] == i){
					v1 = T[V[i].Tind[ti]].Vind[0];
					v2 = T[V[i].Tind[ti]].Vind[1];
				}

				found = FALSE;					/* Add v1 as a neighbour if it is not already
												   a neighbour.								*/
				for (k = 0; k < V[i].vn; k++)
					if (v1 == V[i].VN[k])
						found = TRUE;
				if (not found)
					V[i].VN[V[i].vn++] = v1;

				found = FALSE;					/* Add v2 as a neigbour if it is not already
												   a neighbour.								*/
				for (k = 0; k < V[i].vn; k++)
					if (v2 == V[i].VN[k])
						found = TRUE;
				if (not found)
					V[i].VN[V[i].vn++] = v2;
			
				ca = CosAngle(i,v1,v2,V);
				if (ca < 0){					/* If this triangle is an obtuse angle		*/
					count = Split(V[i].Tind[ti],i,v1,v2,
						T,V,NonSplitTnum,Vnum);
					if (count > mcount)			/* Update m count.							*/
						mcount = count;
				} 
				else {							/* If no splitting was nessesery create
												   the stencil for this vertex and triangle.*/
					V[i].ST[V[i].si].Ctheta = ca;
					V[i].ST[V[i].si].v1 = v1;
					V[i].ST[V[i].si].l1 = Length(i,v1,V);
					V[i].ST[V[i].si].v2 = v2;
					V[i].ST[V[i].si].l2 = Length(i,v2,V);
					V[i].si++;
				}
			}
		}
	}

	for(ind = 0; ind < Vnum; ind++)					/* Calculate the data for each stencil.	*/
		for(p_st = V[ind].ST, si_count = V[ind].si - 1; si_count >= 0; p_st++, si_count--){
			p_st->Stheta =				1 - Sqr(p_st->Ctheta);
			p_st->Ctheta_mul_l1_div_l2_minus_one_mul_2 =p_st->Ctheta * p_st->l1 * 2 / p_st->l2 - 2;
			p_st->Ctheta_mul_l2_div_l1_minus_one_mul_2 =p_st->Ctheta * p_st->l2 * 2 / p_st->l1 - 2;
			p_st->sqr_l1 = p_st->l1 * p_st->l1;
			p_st->sqr_l2 = p_st->l2 * p_st->l2;
			p_st->shortcut1_1 = 1 - p_st->sqr_l1 * (p_st->Ctheta_mul_l2_div_l1_minus_one_mul_2 + 1) / p_st->sqr_l2;
			p_st->shortcut1_2 = 1 - p_st->sqr_l2 * (p_st->Ctheta_mul_l1_div_l2_minus_one_mul_2 + 1) / p_st->sqr_l1;
			p_st->shortcut2_1 = - p_st->Stheta * p_st->sqr_l1;
			p_st->shortcut2_2 = - p_st->Stheta * p_st->sqr_l2;
		}

	//mexPrintf("\nNumber of unfoldings = %d\n",mcount);
}

/************************************************************************************************/
/* This function check the arguments of maxFunction and returns the data that can be extracted	*/
/* from these arguments.																		*/
/* When this function return FALSE maxFunction should terminate immidiatly.						*/
/************************************************************************************************/

boolean phrase_args(int nlhs,			/* Number of expected output arrays */
			int nrhs,					/* Number of input arrays.			*/
			const mxArray *prhs[],		/* array of pointers to input arguments */
			int *number_of_sources,		/* Returned number of sources.		*/
			int **xCoordinates,			/* Returned array of sources X coordinates	*/
			int **yCoordinates,			/* Returned array of sources Y cooridnates	*/
			int *XSIZE, int *YSIZE,		/* Returnes sizes of the surfaces.	*/
			const double **CoordinatesData	/* The coordinates of the surface as a one-D array.	*/
			){

	const int *DimensionsArray;			/* Pointer to arrays of the size of each dimension. */

										/* Check input arguemnts.			*/
	if(nlhs != 1 && nlhs != 2){
		mexPrintf("Function returns exactly one value, an array of distances\n");
		return FALSE;
	}
	else if(nrhs != 2  && nrhs != 3 && nrhs != 5){
		mexPrintf("Function should get 2, 3 or 5 arrays, see documantation for the types of arrays expected.\n");
		return FALSE;
	}
	else if(nrhs != 5 && mxGetClassID(prhs[0]) != mxDOUBLE_CLASS){
		mexPrintf("The surface array should contain double data.\n");
		return FALSE;
	}
	else if(nrhs == 5 && mxGetClassID(prhs[0]) != mxINT32_CLASS){
		mexPrintf("The triangles array should contain 'int32' data.\n");
		return FALSE;
	}
	else if(mxGetClassID(prhs[1]) != mxINT32_CLASS){
		mexPrintf("The sources array should contain integer data.\n");
		return FALSE;
	}
	else if(nrhs == 3 && mxGetClassID(prhs[2]) != mxINT32_CLASS){
		mexPrintf("The area dimensions array should contain integer data.\n");
		return FALSE;
	}
	else if(nrhs != 5 && (mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetDimensions(prhs[0])[1] != 1) && mxGetNumberOfDimensions(prhs[0]) != 1){
		mexPrintf("The surface array should be a one dimensional array.\n");
		return FALSE;
	}
	else if(nrhs == 3 && (mxGetNumberOfDimensions(prhs[1]) != 2 || mxGetDimensions(prhs[1])[1] != 2)){
		mexPrintf("The sources array should be a two dimensional array with two coloums when the dimensions are supplied\n");
		return FALSE;
	}
	else if((nrhs == 2 || nrhs == 5) && mxGetNumberOfDimensions(prhs[1]) != 1 && (mxGetNumberOfDimensions(prhs[1]) != 2 || mxGetDimensions(prhs[1])[1] != 1)){
		mexPrintf("The sources array should be a one dimensional array when the dimensions are not supplied\n");
		return FALSE;
	}
	else if(nrhs == 3 && mxGetNumberOfElements(prhs[2]) != 2){
		mexPrintf("The dimensions array must contian exactly two elements.\n");
		return FALSE;
	}
	else if(nrhs == 5 && (mxGetNumberOfElements(prhs[2]) != mxGetNumberOfElements(prhs[3]) ||
		mxGetNumberOfElements(prhs[3]) != mxGetNumberOfElements(prhs[4]))){
		mexPrintf("All the coordinates vector must agree with each other in the number of elements.\n");
		return FALSE;
	}
	else if(nrhs == 5 && (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS || mxGetClassID(prhs[3]) != mxDOUBLE_CLASS || mxGetClassID(prhs[4]) != mxDOUBLE_CLASS)){
		mexPrintf("All the coordinates vector must be of type 'int32'.\n");
		return FALSE;
	}
	else if(nrhs == 5 && (mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetDimensions(prhs[0])[1] != 3)){
		mexPrintf("The triangles array should be a two dimensional array with three coloums.\n");
		return FALSE;
	}

	DimensionsArray = mxGetDimensions(prhs[1]);
	*number_of_sources = DimensionsArray[0];			/* The number of sources is the number of rows of the sources array.	*/
	*xCoordinates = mxGetData(prhs[1]);					/* Set the x and y coordinates arrays to the place in the data of		*/
														/* of the sources array where their colom start. If there are only 2,
														   arguments, the yCooridnates array should be an array of zeros.		*/
/*	if(nrhs == 3)
		*yCoordinates = *xCoordinates + *number_of_sources;	
	else{
		*yCoordinates = calloc(*number_of_sources, sizeof(int));
		if(yCoordinates == NULL){
			fprintf(stderr, "Out of memory\n");
			exit(-1);
		}
	}*/
	*yCoordinates = NULL;

	if(nrhs == 3 || nrhs == 2){						/* If the first array is a surface array.								*/
		*CoordinatesData = mxGetData(prhs[0]);			/* The coordinates are the data of the surface array.					*/

		if(nrhs == 3){
			DimensionsArray = mxGetData(prhs[2]);		/* Get the sizes of the surface											*/
			*XSIZE = DimensionsArray[1];
			*YSIZE = DimensionsArray[0];
			if(**CoordinatesData != *XSIZE * *YSIZE){	/* Check the coralation of the surface and dimensions arrays */
				mexPrintf("The surface array contains different number of vertexes then expected by dimensions array\n");
				return FALSE;			
			}
		}
		else{
			*XSIZE = (int)**CoordinatesData;			/* The number of vertexes.												*/
			*YSIZE = 1;
		}

														/* Check that the size of the surface array agree with the numbers in
														   it. */
		if(3 * (**CoordinatesData + *(*CoordinatesData + 1)) != mxGetNumberOfElements(prhs[0]) - 2){
			mexPrintf("Illegal surface array, size does not much vertexes and triangles numbers.\n");
			return FALSE;
		}
	}
	else{												/* If the coordinates vectors method is used.							*/
		*XSIZE = mxGetNumberOfElements(prhs[2]);
		*YSIZE = 1;
	}

	return TRUE;

}





/***************************************************************************/
/* This function adds a u1 and u2 value to every point, so like on		   */
/* parametric sufface, it is possible to name a point or an area.		   */
/***************************************************************************/
void ParametricSurface(struct Vertex V[], int XSIZE, int YSIZE){
#ifdef FULL_OUTPUT	
	int i, j;
	for(i = 0; i < XSIZE; i++)
		for(j = 0; j < YSIZE; j++){
			V[i * YSIZE + j ].u1 = i;
			V[i * YSIZE + j ].u2 = j;
		}
#endif
}
/***************************************************************************/
/* This function is the entry point of the dll, matlab shuld call this	   */
/* function with a surface matrix and a sources matrix, the dll return the */
/* distances matrix and can print various types of output, according to	   */
/* the compilation flags that are defined.								   */
/***************************************************************************/
void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[]         /* array of pointers to input arguments */
)
 {

	struct Triangle *T; 	/* triangles, numbered 1..Tnum */
							/* with the 3 indices for vertices */
	struct Vertex   *V; 	/*  vertices coordinates 1..VNum */

	int	Tnum, Vnum,k=257,
#if defined(SMOOTH_GEODESIC) || defined(LINEAR_GEODESIC)
		l,
#endif
	m=0;
	int NonSplitTnum;		/* The number of triangles with no splitting. */
							/* Names of files that are used for writing the results. */

	int *xStart, *yStart;	/* The x and y coordinates of the sources */


	int number_of_sources;	/* The number of sources in the file	*/

	int i, j;				/* Indexes for the loops.				*/

	const double *CoordinatesData;	/* A pointer to the one dimensional array of input data from the surface file.	*/

	double *DistancesData;			/* A pointer to the the one dimensional data of the result matrix.				*/

	double tm;

//	switch_priority(TRUE);

	//mexPrintf("%s\n", LastUpdate);
							/* Parse the command line arguments.*/
	if(!phrase_args(nlhs, nrhs, prhs, &number_of_sources, &xStart, &yStart, &XSIZE, &YSIZE, &CoordinatesData))
		return;

	PICSIZE = YSIZE*XSIZE;		/* Calculate the number of vertexes
							   nessesery.							*/

	//mexPrintf("Read Surface\n");
	if(nrhs != 5)			/* If there is a surface array present. */
		ReadSurface(&V, &T, &Vnum, &NonSplitTnum, &Tnum, CoordinatesData);
	else
		CombineVectors(&V, &T, &Vnum, &NonSplitTnum, &Tnum, prhs);

	//mexPrintf("InitGraph\n");
	InitGraph(T,V, Vnum,NonSplitTnum,Tnum); /* init V to isolated vertices */



	//mexPrintf("March\n");
	init_heap(Vnum);			/* Initialize the heap */
	SetSourceSet(xStart, yStart, number_of_sources, V, Vnum);
	

	for(i = number_of_sources - 1; i >= 0; i--){
		MarchAlg(V,Vnum, i);	/* Perform the march */			
	}

//	printf("\nELAPSED TIME: %f ms\n", tm);

														/* Write the distances matrix to
														   a file.							*/
									/* Allocate the result matrix.									*/
	plhs[0] = mxCreateDoubleMatrix(number_of_sources, number_of_sources, mxREAL);
	DistancesData = mxGetData(plhs[0]);
	for(i = 0; i < number_of_sources; i++){	/* Fill up the result matrix with data.					*/
		for(j = 0; j < number_of_sources; j++){
			*(DistancesData++) = MIN(sources_dis[i][j], sources_dis[j][i]);;
		}
	}

/*	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	*((double *)mxGetPr(plhs[1])) = tm;
*/


	/* Free the resources used by the program. */
	free_heap();
	free(T);
	free(V);
	matrix_free(sources_dis, number_of_sources);

//	switch_priority(FALSE);

	return;
}
